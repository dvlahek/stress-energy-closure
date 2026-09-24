#!/usr/bin/env python3
"""Pre-odd eBOSS weighted RR pilot with sparse exact KD-tree pair geometry.

This is deliberately NOT the final survey window. It uses two predeclared
independent seeds, nested random samples (12k, 40k) and the two *simplified*
published flat-LambdaCDM distance backgrounds. It reports fixed-scale RR
subsampling error and a separate distance-convention sensitivity test. All
inputs are random FITS catalogues; no observed galaxy pair or odd vector is
read, and no physical wake template is fitted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

from inspect_eboss_dr16_selection import WEIGHTS, download
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import (
    NUMERICAL_ZERO_WEIGHT_TOL, fetch_with_retry,
)
from check_eboss_dr16_fiducials import (
    C_KM_S, MODELS, distance_mpc_over_h_quad,
)
from check_eboss_rr_scalar import independent_scalar_rr

ROOT = Path(__file__).resolve().parents[1]
REAL_REF = ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json"
MOCK_REF = (
    ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
)
MOCK_ID = 1
ZLO, ZHI = 0.8, 0.9
SEPARATION_EDGES = np.asarray((20., 40., 60., 80., 100., 120., 140.))
N_MU = 48
THETA_MIN_DEG = 0.05
SMALL, LARGE = 12000, 40000
SEEDS = (4071, 8159)
COSMO_PRIMARY = "published_configuration_space_LRG_ELG_multitracer"
COSMO_SENSITIVITY = "published_DR16_ELG_BAO"


def distance_grid(model_name: str) -> tuple[np.ndarray, np.ndarray, dict]:
    params = MODELS[model_name]
    om = params["Omega_m"]
    zgrid = np.linspace(0.0, 1.0, 20001)
    integrand = 1.0 / np.sqrt(om * (1.0 + zgrid)**3 + 1.0 - om)
    # c/H0 in Mpc times h = c/100 in Mpc/h.
    chi = (C_KM_S / 100.0) * cumulative_trapezoid(
        integrand, zgrid, initial=0.0)
    errors = []
    for z in (0.6, 0.8, 0.85, 0.9, 1.0):
        exact = distance_mpc_over_h_quad(z, params["h"], om)
        estimate = float(np.interp(z, zgrid, chi))
        errors.append(abs(estimate - exact) / exact)
    maxerr = max(errors)
    if maxerr > 1e-8:
        raise RuntimeError(f"Distance interpolation does not match quadrature: {maxerr}")
    return zgrid, chi, {
        "model": model_name,
        "h": params["h"], "Omega_m": om,
        "max_relative_interpolation_error_on_fixed_grid": maxerr,
        "radiation_or_massive_neutrino_terms_included": False,
    }


def cartesian(cat: tuple[np.ndarray, ...], distance) -> tuple[np.ndarray, ...]:
    ra, dec, z, w = cat
    lon, lat = np.deg2rad(ra), np.deg2rad(dec)
    c = np.cos(lat)
    unit = np.column_stack((
        c * np.cos(lon), c * np.sin(lon), np.sin(lat)))
    r = np.asarray(distance(z), dtype="f8")
    return r[:, None] * unit, unit, np.asarray(w, dtype="f8")


def kdtree_rr(
    cat1: tuple[np.ndarray, ...], cat2: tuple[np.ndarray, ...],
    sedges: np.ndarray, muedges: np.ndarray,
    distance, theta_min_deg: float, chunk: int = 128,
    max_neighbour_pairs: int = 150_000_000,
) -> tuple[np.ndarray, dict]:
    """Exactly count all candidate pairs in these *sampled* random catalogues."""
    x1, u1, w1 = cartesian(cat1, distance)
    x2, u2, w2 = cartesian(cat2, distance)
    tree = cKDTree(x2, leafsize=32)
    radius = float(sedges[-1])
    pairs_possible = tree.query_ball_point(
        x1, radius, return_length=True, workers=2)
    candidates = int(np.sum(pairs_possible, dtype=np.int64))
    if candidates > max_neighbour_pairs:
        raise RuntimeError(
            f"RR candidate-pair budget exceeded: {candidates} > {max_neighbour_pairs}")
    nbins = (len(sedges) - 1, len(muedges) - 1)
    out = np.zeros(nbins, dtype="f8")
    cos_min_theta = float(np.cos(np.deg2rad(theta_min_deg)))
    accepted = 0
    for first in range(0, len(x1), chunk):
        last = min(first + chunk, len(x1))
        neighbour_lists = tree.query_ball_point(
            x1[first:last], radius, workers=2)
        count = np.fromiter((len(v) for v in neighbour_lists),
                            dtype=np.int64, count=last - first)
        if not np.any(count):
            continue
        ii = np.repeat(np.arange(first, last, dtype=np.int64), count)
        jj = np.concatenate(
            [np.asarray(v, dtype=np.int64)
             for v in neighbour_lists if len(v)])
        delta = x2[jj] - x1[ii]
        s2 = np.einsum("ij,ij->i", delta, delta)
        s = np.sqrt(s2)
        midpoint2 = x2[jj] + x1[ii]
        midpoint_norm = np.linalg.norm(midpoint2, axis=1)
        denom = s * midpoint_norm
        good = denom > 0
        mu = np.zeros(len(s), dtype="f8")
        dot = np.einsum("ij,ij->i", delta, midpoint2)
        np.divide(dot, denom, out=mu, where=good)
        cosangle = np.einsum("ij,ij->i", u1[ii], u2[jj])
        si = np.searchsorted(sedges, s, side="right") - 1
        mi = np.searchsorted(muedges, mu, side="right") - 1
        ok = (
            good & (si >= 0) & (si < nbins[0])
            & (mi >= 0) & (mi < nbins[1])
            & (cosangle <= cos_min_theta) & (cosangle > -1.0))
        if not np.any(ok):
            continue
        accepted += int(np.count_nonzero(ok))
        flattened = si[ok] * nbins[1] + mi[ok]
        weighted = w1[ii[ok]] * w2[jj[ok]]
        out += np.bincount(
            flattened, weights=weighted,
            minlength=nbins[0] * nbins[1]).reshape(nbins)
    normalization = float(
        np.sum(w1, dtype="f8") * np.sum(w2, dtype="f8"))
    if not np.isfinite(normalization) or normalization <= 0:
        raise ValueError("Invalid random-weight normalization")
    return out, {
        "first_randoms": len(w1), "second_randoms": len(w2),
        "candidate_neighbour_pairs": candidates,
        "accepted_weighted_pair_count": accepted,
        "sum_first_weights": float(np.sum(w1, dtype="f8")),
        "sum_second_weights": float(np.sum(w2, dtype="f8")),
        "pair_normalization": normalization,
    }


def histogram_distance(a: np.ndarray, b: np.ndarray) -> dict:
    total = float(np.sum(np.abs(b)))
    if total <= 0:
        raise ValueError("No RR weight in benchmark histogram")
    return {
        "weighted_normalized_RR_L1_relative": float(
            np.sum(np.abs(a - b)) / total),
        "per_sep_bin_L1_relative": [
            (float(np.sum(np.abs(x - y)) / np.sum(np.abs(y)))
             if np.sum(np.abs(y)) > 0 else None)
            for x, y in zip(a, b)
        ],
    }


def select_catalogue_draws(path: Path, maxrows: int, seeds: tuple[int, ...],
                           seed_offset: int) -> tuple[dict, dict]:
    selections = {}
    with fits.open(path, memmap=path.suffix != ".gz") as hdus:
        hdus.verify("exception")
        hdu = hdus[1]
        table = hdu.data
        if not {"RA", "DEC", "Z", *WEIGHTS}.issubset(table.columns.names):
            raise ValueError("Missing required FITS columns")
        ra = np.asarray(table["RA"], dtype="f8")
        dec = np.asarray(table["DEC"], dtype="f8")
        z = np.asarray(table["Z"], dtype="f8")
        candidate = (
            np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
            & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90)
            & (z >= ZLO) & (z < ZHI))
        eligible = candidate.copy()
        diagnostic = {}
        for column in WEIGHTS:
            w = np.asarray(table[column], dtype="f8")
            finite = np.isfinite(w)
            zero = (candidate & finite &
                    (np.abs(w) <= NUMERICAL_ZERO_WEIGHT_TOL))
            tiny = candidate & finite & (np.abs(w) <= 1e-20)
            if column == "WEIGHT_SYSTOT" and np.count_nonzero(zero) != np.count_nonzero(tiny):
                raise ValueError("The predeclared systematic-weight numerical gap failed")
            if np.any(candidate & (~finite | (w < -NUMERICAL_ZERO_WEIGHT_TOL))):
                raise ValueError(f"Nonfinite or significant negative {column}")
            floor = NUMERICAL_ZERO_WEIGHT_TOL if column == "WEIGHT_SYSTOT" else 0.0
            eligible &= finite & (w > floor)
            diagnostic[column] = {
                "numerical_zero_candidate_rows": int(np.count_nonzero(zero)),
                "candidate_nonpositive_rows": int(np.count_nonzero(
                    candidate & finite & (w <= 0))),
            }
        population = np.flatnonzero(eligible)
        if len(population) < maxrows:
            raise RuntimeError(f"Insufficient eligible random rows: {len(population)} < {maxrows}")
        for seed in seeds:
            rng = np.random.default_rng(seed + seed_offset)
            draw = rng.choice(population, size=maxrows, replace=False)
            data = (
                np.asarray(ra[draw], dtype="f8").copy(),
                np.asarray(dec[draw], dtype="f8").copy(),
                np.asarray(z[draw], dtype="f8").copy(),
            )
            weights = np.ones(maxrows, dtype="f8")
            for column in WEIGHTS:
                weights *= np.asarray(table[column][draw], dtype="f8")
            if not (np.all(np.isfinite(weights)) and np.all(weights > 0)):
                raise ValueError("Nonfinite or nonpositive sampled random weight product")
            # Permute the independently fixed maximum draw once: nested
            # prefix of 12k is contained in the corresponding 40k draw.
            ordering = np.random.default_rng(
                seed + seed_offset + 100000).permutation(maxrows)
            sample = tuple(
                np.asarray(arr[ordering], dtype="f8")
                for arr in (*data, weights))
            digest = hashlib.sha256()
            for arr in sample:
                digest.update(np.ascontiguousarray(arr).tobytes())
            selections[seed] = {
                "catalogue": sample, "sample_sha256": digest.hexdigest(),
            }
    return selections, {
        "FITS_header_rows": int(hdu.header["NAXIS2"]),
        "fixed_diagnostic_z_window": [ZLO, ZHI],
        "eligible_random_rows_in_window": int(len(population)),
        "candidate_weight_column_diagnostics_in_window": diagnostic,
    }


def self_test() -> None:
    rng = np.random.default_rng(5821)
    first = (
        rng.uniform(120, 125, 31), rng.uniform(10, 15, 31),
        rng.uniform(0.82, 0.88, 31), rng.uniform(0.7, 1.8, 31))
    second = (
        rng.uniform(121, 126, 29), rng.uniform(11, 16, 29),
        rng.uniform(0.82, 0.88, 29), rng.uniform(0.6, 1.9, 29))
    sep = np.asarray((0., 20., 40., 60., 90., 140., 200.))
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, 49)
    test_distance = lambda z: np.asarray(z, dtype="f8") * 2800.0
    fast, norm = kdtree_rr(
        first, second, sep, mu, test_distance, 0.05, chunk=9)
    scalar, count, weight_norm = independent_scalar_rr(
        first, second, sep, mu, 0.05)
    np.testing.assert_allclose(fast, scalar, rtol=1e-11, atol=5e-10)
    assert norm["accepted_weighted_pair_count"] == count
    assert math.isclose(norm["pair_normalization"], weight_norm, rel_tol=1e-13)
    for model in (COSMO_PRIMARY, COSMO_SENSITIVITY):
        distance_grid(model)
    print("Sparse KD-tree and independent scalar RR closure passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/rr_kdtree_pilot.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/rr_kdtree_catalogues")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--chunk", type=int, default=128)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.chunk < 8 or args.timeout <= 0:
        parser.error("Invalid pair chunk size or HTTP timeout")
    real_ref = json.loads(REAL_REF.read_text(encoding="utf-8"))
    mock_ref = json.loads(MOCK_REF.read_text(encoding="utf-8"))
    if (real_ref["status"] != "four_random_catalogues_checked" or
            mock_ref["status"] != "mock_random_sample_selection_compatible"):
        raise RuntimeError("Upstream random catalogues have not passed input audits")
    real_sha = {(r["tracer"], r["cap"]): r["sha256"]
                for r in real_ref["randoms"]}
    mock_sha = {
        (int(r["id"]), r["tracer"], r["cap"]): r["compressed_file_sha256"]
        for r in mock_ref["mock_random_catalogues"]
    }
    distance_models = {}
    distances = {}
    for model in (COSMO_PRIMARY, COSMO_SENSITIVITY):
        zgrid, chi, details = distance_grid(model)
        distance_models[model] = details
        distances[model] = lambda z, grid=zgrid, vals=chi: np.interp(z, grid, vals)
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, N_MU + 1)
    cache = Path(args.cache_dir)
    report_cases, errors = [], []
    for cap_idx, cap in enumerate(("NGC", "SGC")):
        for kind_idx, kind in enumerate(("observed", "realistic_mock")):
            input_files, draws = {}, {}
            for tracer_idx, tracer in enumerate(("LRG", "ELG")):
                real = kind == "observed"
                if real:
                    filename, expected_rows = RANDOMS[(tracer, cap)]
                    local = cache / kind / filename
                else:
                    relative = mock_path(
                        "eBOSS_" + tracer, cap, "ran", MOCK_ID)
                    local = cache / kind / Path(relative).name
                try:
                    if real:
                        local, digest, nbytes = download(
                            filename, local.parent, 1024 * 1024 * 1024, args.timeout)
                        if digest != real_sha[(tracer, cap)]:
                            raise RuntimeError("Observed random SHA256 changed")
                    else:
                        digest, nbytes = fetch_with_retry(
                            MOCK_BASE + relative, local, 512 * 1024 * 1024,
                            args.timeout)
                        if digest != mock_sha[(MOCK_ID, tracer, cap)]:
                            raise RuntimeError("Realistic-mock random SHA256 changed")
                    seed_offset = (
                        100 * cap_idx + 10 * kind_idx + tracer_idx)
                    selected, diagnostic = select_catalogue_draws(
                        local, LARGE, SEEDS, seed_offset)
                    if real and diagnostic["FITS_header_rows"] != expected_rows:
                        raise RuntimeError("Public FITS NAXIS2 changed")
                    draws[tracer] = selected
                    input_files[tracer] = {
                        "filename": local.name, "file_SHA256": digest,
                        "file_bytes": nbytes, **diagnostic,
                        "per_seed_sample_SHA256": {
                            str(seed): selected[seed]["sample_sha256"]
                            for seed in SEEDS},
                    }
                    print("KDTREE_RR_INPUT_OK", kind, cap, tracer,
                          "eligible", diagnostic["eligible_random_rows_in_window"],
                          flush=True)
                except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                    errors.append(f"{kind}/{cap}/{tracer}: {exc}")
                    print("KDTREE_RR_INPUT_ERROR", errors[-1], flush=True)
                finally:
                    if local.exists():
                        local.unlink()
            if set(draws) != {"LRG", "ELG"}:
                continue
            seeds_output, high_seed_histograms = [], {}
            for seed in SEEDS:
                per_size = {}
                for n in (SMALL, LARGE):
                    a = tuple(x[:n] for x in draws["LRG"][seed]["catalogue"])
                    b = tuple(x[:n] for x in draws["ELG"][seed]["catalogue"])
                    hist, norm = kdtree_rr(
                        a, b, SEPARATION_EDGES, muedges,
                        distances[COSMO_PRIMARY], THETA_MIN_DEG,
                        chunk=args.chunk)
                    normalized = hist / norm["pair_normalization"]
                    per_size[n] = {"normalized": normalized, "norm": norm}
                    print("KDTREE_RR", kind, cap, "seed", seed,
                          "sample", n, "accepted",
                          norm["accepted_weighted_pair_count"], flush=True)
                high_seed_histograms[seed] = per_size[LARGE]["normalized"]
                distance_sensitivity_hist, sensitivity_norm = kdtree_rr(
                    tuple(x[:SMALL] for x in draws["LRG"][seed]["catalogue"]),
                    tuple(x[:SMALL] for x in draws["ELG"][seed]["catalogue"]),
                    SEPARATION_EDGES, muedges,
                    distances[COSMO_SENSITIVITY], THETA_MIN_DEG,
                    chunk=args.chunk)
                distances_compared = histogram_distance(
                    distance_sensitivity_hist / sensitivity_norm["pair_normalization"],
                    per_size[SMALL]["normalized"])
                seeds_output.append({
                    "seed": seed,
                    "small_random_objects_per_tracer": SMALL,
                    "large_random_objects_per_tracer": LARGE,
                    "small_pair_count": per_size[SMALL]["norm"][
                        "accepted_weighted_pair_count"],
                    "large_pair_count": per_size[LARGE]["norm"][
                        "accepted_weighted_pair_count"],
                    "small_to_large_normalized_RR":
                        histogram_distance(
                            per_size[SMALL]["normalized"],
                            per_size[LARGE]["normalized"]),
                    "fixed_small_sample_distance_convention_sensitivity":
                        distances_compared,
                    "large_normalized_RR_s_mu":
                        per_size[LARGE]["normalized"].tolist(),
                    "large_weighted_RR_sum_fraction_in_test_bins": float(
                        np.sum(per_size[LARGE]["normalized"])),
                })
            between = histogram_distance(
                high_seed_histograms[SEEDS[0]], high_seed_histograms[SEEDS[1]])
            report_cases.append({
                "catalogue": kind, "cap": cap, "input_files": input_files,
                "two_seed_large_sample_comparison": between,
                "seeds": seeds_output,
            })
            print("KDTREE_RR_SEED_COMPARISON", kind, cap,
                  "n", LARGE,
                  "relative_l1", between[
                      "weighted_normalized_RR_L1_relative"], flush=True)
    complete = len(report_cases) == 4 and not errors
    report = {
        "study": "eBOSS pre-odd high-density sparse KD-tree random-pair window pilot",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": ("high_density_RR_pilot_complete" if complete else "partial"),
        "fixed_upstream_real_reference": str(REAL_REF.relative_to(ROOT)),
        "fixed_upstream_mock_reference": str(MOCK_REF.relative_to(ROOT)),
        "mock_realization_id": MOCK_ID,
        "candidate_interval_frozen_for_inference": False,
        "pilot_z_interval": [ZLO, ZHI],
        "pilot_separation_edges_Mpc_over_h": SEPARATION_EDGES.tolist(),
        "pilot_mu_bins": N_MU, "theta_min_deg": THETA_MIN_DEG,
        "orientation": "first LRG then ELG", "line_of_sight": "midpoint",
        "sample_size_small_and_large_per_tracer": [SMALL, LARGE],
        "fixed_seeds": list(SEEDS),
        "weight_model": "provisional product of four released weight columns",
        "cosmologies": distance_models,
        "primary_cosmology": COSMO_PRIMARY,
        "predeclared_distance_sensitivity_cosmology": COSMO_SENSITIVITY,
        "cosmology_radiation_or_neutrino_terms_finalized": False,
        "results": report_cases, "errors": errors,
        "observed_galaxy_catalogues_read": False,
        "observed_odd_data_vector_read": False,
        "physical_template_tuned": False,
        "full_random_density_window_computed": False,
        "full_mock_ensemble_validated": False,
        "mock_cross_covariance_computed": False,
        "note": (
            "All RR counts are exact for each fixed sampled random catalogue "
            "within the stated geometry and simplified flat-LCDM mapping. "
            "They are not full-density RR counts or an inferentially frozen "
            "survey window. The two fixed seeds and size levels quantify "
            "subsampling error. The alternative published fiducial is an "
            "input-only distance sensitivity, never a template fit."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("KDTREE_RR_PILOT", report["status"], out, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
