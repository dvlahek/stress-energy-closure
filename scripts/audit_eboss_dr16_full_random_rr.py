#!/usr/bin/env python3
"""Full-density, random-only eBOSS LRG->ELG RR(s,mu) candidate-bin audit.

Only the four SHA-pinned public eBOSS DR16 random FITS files enter. We count
all eligible LRG x ELG random pairs for each cap and four predeclared
candidate redshift bins (0.6-0.7, 0.7-0.8, 0.8-0.9, 0.9-1.0) using
midpoint LOS and the validated pure-eBOSS multiplicative random weights.

Two disjoint fixed random half-samples plus a nested quarter-sample quantify
sampling sensitivity on the same (s,mu) grid. The predeclared secondary
distance mapping is checked on the complete high-z bin only. These RR
counts are a *candidate window input*, not a frozen survey-window
convolution, mock covariance, or observed galaxy odd-sector measurement.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from eboss_dr16_fiducial import (
    DISTANCE_CONVENTIONS, PRIMARY_GEOMETRY, SECONDARY_GEOMETRY,
    WEIGHT_COLUMNS, comoving_mpc_over_h, validated_weight_product,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import Z_EDGES, download
from audit_eboss_dr16_rr_density import histogram_distance

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "source_data/eboss_dr16_weighted_normalization_2026-09-24.json"
SEPARATION_EDGES = (20., 40., 60., 80., 100., 120., 140.)
MU_BINS = 240
THETA_MIN_DEG = 0.05
RNG_SEED = 4071
BINS = tuple((float(Z_EDGES[i]), float(Z_EDGES[i + 1])) for i in range(4))


def read_random_bins(
    path: Path, tracer: str, cap: str, expected_rows: int,
) -> tuple[list[tuple[np.ndarray, ...]], dict]:
    result = []
    with fits.open(path, memmap=True) as h:
        h.verify("exception")
        tables = [x for x in h if isinstance(x, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected one observed random FITS BINTABLE")
        table = tables[0]
        if int(table.header["NAXIS2"]) != expected_rows:
            raise ValueError("Full observed random row count changed")
        missing = {"RA", "DEC", "Z", *WEIGHT_COLUMNS} - set(
            table.columns.names)
        if missing:
            raise ValueError(f"Missing public random fields: {sorted(missing)}")
        z = np.asarray(table.data["Z"], dtype="f8")
        ra = np.asarray(table.data["RA"], dtype="f8")
        dec = np.asarray(table.data["DEC"], dtype="f8")
        if (not np.isfinite(z).all() or not np.isfinite(ra).all()
                or not np.isfinite(dec).all()
                or np.any(ra < 0) or np.any(ra >= 360)
                or np.any(dec < -90) or np.any(dec > 90)):
            raise ValueError("Public random coordinate or redshift invalid")
        for lo, hi in BINS:
            selected = (z >= lo) & (z < hi)
            weight_fields = {
                name: np.asarray(table.data[name][selected], dtype="f8")
                for name in WEIGHT_COLUMNS
            }
            weight, retained = validated_weight_product(weight_fields)
            x = (
                np.asarray(ra[selected][retained], dtype="f8"),
                np.asarray(dec[selected][retained], dtype="f8"),
                np.asarray(z[selected][retained], dtype="f8"),
                np.asarray(weight[retained], dtype="f8"),
            )
            if not len(x[0]) or np.any(x[3] <= 0):
                raise ValueError(
                    f"Empty/invalid observed random slice {tracer} {cap} {lo}-{hi}")
            result.append(x)
            print(f"FULL_RR_RANDOM_BIN_OK {cap} {tracer} "
                  f"z={lo:.1f}-{hi:.1f} raw={np.count_nonzero(selected)} "
                  f"kept={len(x[0])}", flush=True)
    summary = {
        "cap": cap, "tracer": tracer,
        "retained_per_candidate_bin": [int(len(x[0])) for x in result],
        "weighted_sum_per_candidate_bin": [
            float(np.sum(x[3], dtype="f8")) for x in result
        ],
    }
    return result, summary


def coordinates(cat, profile: str):
    return (
        cat[0], cat[1], comoving_mpc_over_h(cat[2], profile), cat[3]
    )


def compile_rr(
    first, second, sedges: np.ndarray, muedges: np.ndarray, nthreads: int,
) -> tuple[np.ndarray, dict]:
    from pycorr import TwoPointCounter

    result = TwoPointCounter(
        mode="smu", edges=(sedges, muedges),
        positions1=[first[0], first[1], first[2]],
        positions2=[second[0], second[1], second[2]],
        weights1=first[3], weights2=second[3],
        position_type="rdd", los="midpoint",
        selection_attrs={"theta": (THETA_MIN_DEG, 180.0)},
        engine="corrfunc", nthreads=nthreads, compute_sepsavg=False,
    )
    raw = np.asarray(result.wcounts, dtype="f8")
    if raw.shape != (len(sedges) - 1, len(muedges) - 1):
        raise ValueError(f"Unexpected Corrfunc RR grid: {raw.shape}")
    if not np.isfinite(raw).all() or np.any(raw < 0):
        raise ValueError("Compiled RR is nonfinite or negative")
    norm = float(result.wnorm)
    direct_norm = float(
        np.sum(first[3], dtype="f8") * np.sum(second[3], dtype="f8"))
    err = abs(norm - direct_norm) / direct_norm
    if err > 1e-9:
        raise ValueError("RR wnorm disagrees with direct catalogue weight sums")
    if not np.any(raw):
        raise ValueError("No RR support in declared candidate separation range")
    return raw, {
        "sum_weights_first": float(np.sum(first[3], dtype="f8")),
        "sum_weights_second": float(np.sum(second[3], dtype="f8")),
        "pair_weight_normalization": norm,
        "norm_relative_difference": err,
        "raw_rr_weight_sum_in_test_bins": float(np.sum(raw, dtype="f8")),
        "normalized_rr_weight_fraction_in_test_bins":
            float(np.sum(raw, dtype="f8") / norm),
        "positive_s_mu_cells": int(np.count_nonzero(raw)),
    }


def diagnostic_odd_rr(raw: np.ndarray, muedges: np.ndarray) -> list[dict]:
    """Raw odd window geometry only, never galaxy odd multipoles."""
    widths = np.diff(muedges)
    p1 = ((muedges[1:] ** 2 - muedges[:-1] ** 2) / 2) / widths
    p3 = (
        5 / 8 * (muedges[1:] ** 4 - muedges[:-1] ** 4)
        - 3 / 4 * (muedges[1:] ** 2 - muedges[:-1] ** 2)
    ) / widths
    out = []
    for arr in raw:
        total = float(np.sum(arr, dtype="f8"))
        out.append({
            "raw_weighted_rr_total": total,
            "raw_rr_mu_mirror_l1_fraction": (
                float(np.sum(np.abs(arr - arr[::-1])) / total)
                if total else None),
            "raw_rr_mean_P1": float(np.dot(arr, p1) / total)
            if total else None,
            "raw_rr_mean_P3": float(np.dot(arr, p3) / total)
            if total else None,
        })
    return out


def deterministic_subsets(cat, cap_index: int, tracer_index: int,
                          bin_index: int):
    n = len(cat[0])
    if n < 8:
        raise ValueError("Too few randoms for predeclared half-sample checks")
    seed = RNG_SEED + 1000 * cap_index + 100 * bin_index + tracer_index
    ix = np.random.default_rng(seed).permutation(n)
    half = n // 2
    quarter = n // 4
    return {
        "full": cat,
        "half_A": tuple(np.asarray(v[ix[:half]]) for v in cat),
        "half_B": tuple(np.asarray(v[ix[half:]]) for v in cat),
        "quarter_A": tuple(np.asarray(v[ix[:quarter]]) for v in cat),
    }, seed


def self_test() -> None:
    from audit_eboss_dr16_rr_pair_closure import (
        rr_histogram, odd_rr_moments,
    )
    rng = np.random.default_rng(7367)
    def cat(n):
        return (
            rng.uniform(120, 124, n),
            rng.uniform(10, 15, n),
            rng.uniform(0.85, 0.89, n),
            rng.uniform(0.8, 1.7, n),
        )
    a, b = cat(107), cat(113)
    s = np.asarray(SEPARATION_EDGES)
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    ac, bc = coordinates(a, PRIMARY_GEOMETRY), coordinates(b, PRIMARY_GEOMETRY)
    compiled, info = compile_rr(ac, bc, s, mu, 2)
    reference, norm = rr_histogram(
        a, b, s, mu, THETA_MIN_DEG,
        lambda z: comoving_mpc_over_h(z, PRIMARY_GEOMETRY),
        block=37,
    )
    np.testing.assert_allclose(compiled, reference, atol=1e-7, rtol=2e-8)
    assert abs(info["pair_weight_normalization"] /
               norm["pair_normalization"] - 1) < 1e-10
    subs, seed = deterministic_subsets(a, 0, 0, 1)
    assert seed == RNG_SEED + 100
    assert len(subs["quarter_A"][0]) == 107 // 4
    assert len(subs["half_A"][0]) == 107 // 2
    assert sum(len(subs[key][0]) for key in ("half_A", "half_B")) == 107
    result = diagnostic_odd_rr(compiled, mu)
    assert len(result) == len(s) - 1
    print("EBOSS_FULL_DENSITY_RR_INDEPENDENT_SELF_TEST_OK", flush=True)


def write_artifacts(output: Path, arrays: dict, report: dict) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(
        output.with_suffix(".npz"),
        **arrays,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/full_random_rr_audit.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/full_rr_randoms")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout <= 0 or args.threads < 1:
        parser.error("Invalid timeout or Corrfunc thread count")
    reference = json.loads(REFERENCE.read_text())
    if reference["status"] != "sample_weighted_normalization_audited":
        raise ValueError("The full observed data/random normalization gate is missing")
    expected = {
        (r["cap"], r["tracer"]): r
        for r in reference["per_file_records"]
        if r["survey"] == "observed" and r["role"] == "random"
    }
    if len(expected) != 4:
        raise ValueError("Four observed random SHA256 identities are required")
    sedges = np.asarray(SEPARATION_EDGES, dtype="f8")
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    arrays = {
        "sedges_mpc_over_h": sedges, "muedges": muedges,
        "zlo": np.array([v[0] for v in BINS], dtype="f8"),
        "zhi": np.array([v[1] for v in BINS], dtype="f8"),
    }
    reports, errors, inputs = [], [], []
    report = {
        "study": "Full observed eBOSS LRG-ELG random-only candidate-bin RR",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "partial",
        "input_normalization_source": str(REFERENCE.relative_to(ROOT)),
        "caps": ["NGC", "SGC"],
        "candidate_redshift_bins": [list(v) for v in BINS],
        "candidate_redshift_and_separation_selection_frozen": False,
        "separation_edges_mpc_over_h": sedges.tolist(),
        "mu_bin_count": MU_BINS,
        "theta_min_deg": THETA_MIN_DEG,
        "orientation": "LRG->ELG",
        "los": "midpoint",
        "weights": "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP",
        "numerical_zero_selection": "abs(WEIGHT_SYSTOT)<=1e-12",
        "primary_distance": PRIMARY_GEOMETRY,
        "secondary_distance": SECONDARY_GEOMETRY,
        "distance_conventions": DISTANCE_CONVENTIONS,
        "distance_approximations": (
            "flat homogeneous matter+Lambda; Tcmb0=0, no radiation/neutrino "
            "background correction"),
        "fixed_halving_seed_root": RNG_SEED,
        "pairs_computed": "random_R1R2_only",
        "records": reports, "input_random_file_summaries": inputs,
        "errors": errors,
        "observed_galaxy_positions_read": False,
        "observed_odd_vector_read": False,
        "observed_DD_or_DR_computed": False,
        "mock_covariance_computed": False,
        "even_to_odd_forward_model_convolved": False,
        "full_realistic_mock_RR_computed": False,
        "note": (
            "The full observed-random RR histograms cover fixed candidate "
            "bins and are window-input diagnostics only. Their presence "
            "does not freeze the physical model, analysis cuts or survey "
            "selection; observed galaxy odd statistics stay unopened."
        ),
    }
    output = Path(args.out)
    root = Path(args.cache_dir)
    for ci, cap in enumerate(("NGC", "SGC")):
        by_tracer = {}
        for ti, tracer in enumerate(("LRG", "ELG")):
            filename, nrows = RANDOMS[(tracer, cap)]
            local = root / filename
            try:
                local, sha, size = download(
                    filename, root, 1024 * 1024 * 1024, args.timeout)
                if sha != expected[(cap, tracer)]["sha256"]:
                    raise ValueError(f"SHA256 of {filename} changed")
                catalogues, summary = read_random_bins(
                    local, tracer, cap, nrows)
                summary.update({
                    "filename": filename, "full_fits_file_sha256": sha,
                    "downloaded_bytes": size,
                })
                rounded_sum = expected[(cap, tracer)]["weighted_sum_rounded_for_log"]
                actual_sum = sum(summary["weighted_sum_per_candidate_bin"])
                if abs(actual_sum - rounded_sum) / rounded_sum > 5e-6:
                    raise ValueError(
                        "Weighted full-candidate sum changed from retained audit")
                by_tracer[tracer] = catalogues
                inputs.append(summary)
            except (OSError, ValueError, KeyError, MemoryError) as exc:
                errors.append(f"{cap}/{tracer}: {exc}")
                print("FULL_RR_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if local.exists():
                    local.unlink()
        if set(by_tracer) != {"LRG", "ELG"}:
            report["status"] = "partial"
            write_artifacts(output, arrays, report)
            continue
        for iz, (lo, hi) in enumerate(BINS):
            try:
                lrg, elg = by_tracer["LRG"][iz], by_tracer["ELG"][iz]
                lm, lseed = deterministic_subsets(lrg, ci, 0, iz)
                em, eseed = deterministic_subsets(elg, ci, 1, iz)
                results = {}
                # One distance conversion per tracer/level; the cached full
                # comoving coordinates also determine all sub-sample arrays.
                lc = coordinates(lrg, PRIMARY_GEOMETRY)
                ec = coordinates(elg, PRIMARY_GEOMETRY)
                subsets = {}
                # Use the same deterministic permutation as above to retain
                # matching RDD and weights in the compiled counter.
                for tracer_id, cat, profile in (
                    (0, lc, "LRG"), (1, ec, "ELG")
                ):
                    _, seed = deterministic_subsets(
                        lrg if tracer_id == 0 else elg, ci, tracer_id, iz)
                    ix = np.random.default_rng(seed).permutation(len(cat[0]))
                    half = len(ix) // 2
                    quarter = len(ix) // 4
                    subsets[profile] = {
                        "full": cat,
                        "half_A": tuple(v[ix[:half]] for v in cat),
                        "half_B": tuple(v[ix[half:]] for v in cat),
                        "quarter_A": tuple(v[ix[:quarter]] for v in cat),
                    }
                for label in ("full", "half_A", "half_B", "quarter_A"):
                    raw, info = compile_rr(
                        subsets["LRG"][label],
                        subsets["ELG"][label],
                        sedges, muedges, args.threads,
                    )
                    arrays[f"rr_{cap}_z{iz}_{label}"] = raw
                    arrays[f"rrnorm_{cap}_z{iz}_{label}"] = np.array([
                        info["pair_weight_normalization"]
                    ], dtype="f8")
                    results[label] = {
                        "meta": info,
                        "random_rows_lrg": len(subsets["LRG"][label][0]),
                        "random_rows_elg": len(subsets["ELG"][label][0]),
                        "raw_window_odd_moment_diagnostic": (
                            diagnostic_odd_rr(raw, muedges)
                            if label == "full" else None),
                        "npz_key": f"rr_{cap}_z{iz}_{label}",
                    }
                    print(f"FULL_RANDOM_RR_OK cap={cap} z={lo:.1f}-{hi:.1f} "
                          f"sample={label} "
                          f"rr_fraction={info['normalized_rr_weight_fraction_in_test_bins']:.8g}",
                          flush=True)
                full = arrays[f"rr_{cap}_z{iz}_full"]
                norm_full = results["full"]["meta"]["pair_weight_normalization"]
                differences = {}
                for label in ("half_A", "half_B", "quarter_A"):
                    sample = arrays[f"rr_{cap}_z{iz}_{label}"]
                    sample_norm = results[label]["meta"]["pair_weight_normalization"]
                    differences[label] = histogram_distance(
                        sample / sample_norm, full / norm_full)
                differences["half_A_vs_half_B"] = histogram_distance(
                    arrays[f"rr_{cap}_z{iz}_half_A"] /
                    results["half_A"]["meta"]["pair_weight_normalization"],
                    arrays[f"rr_{cap}_z{iz}_half_B"] /
                    results["half_B"]["meta"]["pair_weight_normalization"],
                )
                alternate = None
                if iz == len(BINS) - 1:
                    alt_lrg = coordinates(lrg, SECONDARY_GEOMETRY)
                    alt_elg = coordinates(elg, SECONDARY_GEOMETRY)
                    other, oinfo = compile_rr(
                        alt_lrg, alt_elg, sedges, muedges, args.threads)
                    arrays[f"rr_{cap}_z{iz}_secondary_fiducial_full"] = other
                    arrays[f"rrnorm_{cap}_z{iz}_secondary_fiducial_full"] = (
                        np.array([oinfo["pair_weight_normalization"]], dtype="f8"))
                    alternate = {
                        "meta": oinfo,
                        "secondary_primary_l1": histogram_distance(
                            other / oinfo["pair_weight_normalization"],
                            full / norm_full),
                        "npz_key": f"rr_{cap}_z{iz}_secondary_fiducial_full",
                    }
                record = {
                    "cap": cap, "zlo": lo, "zhi": hi, "bin_index": iz,
                    "candidate_only": True, "sampling_seed_lrg": lseed,
                    "sampling_seed_elg": eseed,
                    "lrg_input_random_rows": len(lrg[0]),
                    "elg_input_random_rows": len(elg[0]),
                    "full_and_thinned": results,
                    "sampling_differences": differences,
                    "secondary_distance_high_z_only": alternate,
                }
                reports.append(record)
                print("FULL_RANDOM_RR_CONVERGENCE", cap, lo, hi,
                      "halfA_l1",
                      differences["half_A"][
                          "normalized_rr_l1_over_larger_sample_l1"],
                      "halfB_l1",
                      differences["half_B"][
                          "normalized_rr_l1_over_larger_sample_l1"],
                      flush=True)
            except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                errors.append(f"{cap}/z{iz}: {exc}")
                print("FULL_RANDOM_RR_ERROR", errors[-1], flush=True)
        report["status"] = "partial"
        write_artifacts(output, arrays, report)
        del by_tracer
        gc.collect()
    complete = len(reports) == 2 * len(BINS) and not errors
    report["status"] = (
        "full_observed_random_candidate_bin_RR_computed"
        if complete else "partial")
    report["full_observed_random_candidate_bin_RR_computed"] = complete
    write_artifacts(output, arrays, report)
    print("EBOSS_FULL_RANDOM_RR_AUDIT", report["status"], output, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
