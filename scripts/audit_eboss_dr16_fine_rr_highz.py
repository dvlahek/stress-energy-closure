#!/usr/bin/env python3
"""One-Mpc/h fine eBOSS random-only RR pilot in the fixed high-z diagnostic bin.

This uses only the four SHA-pinned observed DR16 random catalogues. It
independently counts the same full LRG->ELG random pairs on coarse
20-Mpc/h and fine 1-Mpc/h separation grids and checks that rebinning the
fine histogram reproduces the coarse histogram bin-by-bin. Fixed disjoint
random halves diagnose the fine-grid sampling error. The predeclared second
fiducial is checked on the full fine grid. No observed galaxy data, odd
data vector, mock covariance or physical wake template is read.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre

from eboss_dr16_fiducial import (
    DISTANCE_CONVENTIONS, PRIMARY_GEOMETRY, SECONDARY_GEOMETRY,
)
from audit_eboss_dr16_full_random_rr import (
    ROOT, REFERENCE, BINS, SEPARATION_EDGES, MU_BINS, THETA_MIN_DEG,
    RNG_SEED, read_random_bins, coordinates, compile_rr,
    deterministic_subsets,
)
from audit_eboss_dr16_rr_density import histogram_distance
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import download

FULL_REF = ROOT / "source_data/eboss_dr16_full_observed_rr_2026-09-24.json"
HIGH_Z_INDEX = 3
FINE_STEP = 1.0
THINNING = ("full", "half_A", "half_B")


def fine_edges(step=FINE_STEP) -> np.ndarray:
    lo, hi = float(SEPARATION_EDGES[0]), float(SEPARATION_EDGES[-1])
    n = int(round((hi - lo) / step))
    if n < 1 or abs(lo + n * step - hi) > 1e-12:
        raise ValueError("Fine step does not divide the declared separation interval")
    return np.linspace(lo, hi, n + 1, dtype="f8")


def rebin_to_coarse(raw: np.ndarray, fine: np.ndarray,
                    coarse: np.ndarray) -> np.ndarray:
    if raw.shape[0] != len(fine) - 1:
        raise ValueError("Fine RR array dimensions do not match fine edges")
    result = np.zeros((len(coarse) - 1, raw.shape[1]), dtype="f8")
    for k, (lo, hi) in enumerate(zip(coarse[:-1], coarse[1:])):
        keep = (fine[:-1] >= lo - 1e-12) & (fine[1:] <= hi + 1e-12)
        if not np.any(keep) or abs(fine[keep][0] - lo) > 1e-12:
            raise ValueError("Fine edges do not cover an output separation bin")
        if abs(fine[np.flatnonzero(keep)[-1] + 1] - hi) > 1e-12:
            raise ValueError("Fine edges end before an output separation bin")
        result[k] = np.sum(raw[keep], axis=0, dtype="f8")
    return result


def odd_rr_geometry(raw: np.ndarray, muedges: np.ndarray) -> dict:
    """Raw RR angular moments only: no galaxy multipole or leakage inference."""
    moment = {}
    count = np.sum(raw, axis=1, dtype="f8")
    for ell in (1, 3):
        antideriv = Legendre.basis(ell).integ()
        averaged = (
            np.diff(antideriv(muedges)) / np.diff(muedges)
        )
        numerator = raw @ averaged
        value = np.full(len(count), np.nan)
        np.divide(numerator, count, out=value, where=(count > 0))
        moment[str(ell)] = value
    return moment


def unit_test() -> None:
    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    synthetic = np.ones((len(fine) - 1, 4), dtype="f8")
    combined = rebin_to_coarse(synthetic, fine, coarse)
    assert combined.shape == (6, 4)
    assert np.all(combined == 20)
    m = np.linspace(-1, 1, 25)
    symmetrical = np.ones((len(fine) - 1, 24))
    odd = odd_rr_geometry(symmetrical, m)
    assert all(float(np.max(np.abs(a))) < 1e-14 for a in odd.values())
    rng = np.random.default_rng(49963)
    def cat(n):
        return (
            rng.uniform(140, 144, n),
            rng.uniform(10, 14, n),
            rng.uniform(0.91, 0.99, n),
            rng.uniform(0.7, 1.5, n),
        )
    lrg, elg = cat(91), cat(103)
    ac, bc = coordinates(lrg, PRIMARY_GEOMETRY), coordinates(
        elg, PRIMARY_GEOMETRY)
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    fine_result, fine_info = compile_rr(ac, bc, fine, mu, 2)
    coarse_result, coarse_info = compile_rr(ac, bc, coarse, mu, 2)
    aggregated = rebin_to_coarse(fine_result, fine, coarse)
    np.testing.assert_allclose(
        aggregated, coarse_result, rtol=2e-11, atol=1e-7)
    assert abs(fine_info["pair_weight_normalization"] /
               coarse_info["pair_weight_normalization"] - 1) < 1e-12
    mirrored, _ = compile_rr(bc, ac, fine, mu, 2)
    assert np.sum(np.abs(mirrored[:, ::-1] - fine_result)) / np.sum(
        fine_result) < 1e-11
    print("EBOSS_FINE_RR_SYNTHETIC_REBIN_AND_ORIENTATION_OK", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/fine_rr_highz_audit.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/fine_rr_randoms")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        unit_test()
        return 0
    if args.threads < 1 or args.timeout <= 0:
        parser.error("Invalid thread count or download timeout")
    ref = json.loads(REFERENCE.read_text())
    full_ref = json.loads(FULL_REF.read_text())
    if (ref["status"] != "sample_weighted_normalization_audited"
            or full_ref["status"] !=
            "full_observed_random_candidate_bin_RR_computed"):
        raise ValueError("Full observed input and coarse RR evidence is required")
    expect = {
        (v["cap"], v["tracer"]): v
        for v in ref["per_file_records"]
        if v["survey"] == "observed" and v["role"] == "random"
    }
    assert len(expect) == 4
    high_z = BINS[HIGH_Z_INDEX]
    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    arrays = {
        "fine_sedges_mpc_over_h": fine,
        "coarse_sedges_mpc_over_h": coarse,
        "muedges": mu,
    }
    cases, errors = [], []
    root = Path(args.cache_dir)
    for ci, cap in enumerate(("NGC", "SGC")):
        inputs, cats = {}, {}
        for ti, tracer in enumerate(("LRG", "ELG")):
            filename, expected_rows = RANDOMS[(tracer, cap)]
            path = root / filename
            try:
                path, sha, size = download(
                    filename, root, 1024 * 1024 * 1024, args.timeout)
                if sha != expect[(cap, tracer)]["sha256"]:
                    raise ValueError("Published random FITS SHA changed")
                bins, meta = read_random_bins(path, tracer, cap, expected_rows)
                target = bins[HIGH_Z_INDEX]
                retained = int(len(target[0]))
                expected = next(
                    x for x in full_ref["input_random_redshift_counts"]
                    if x["cap"] == cap and x["tracer"] == tracer
                    and x["zlo"] == high_z[0])
                if retained != expected["retained_rows"]:
                    raise ValueError("High-z retained random count changed")
                cats[tracer] = target
                inputs[tracer] = {
                    "filename": filename, "sha256": sha,
                    "downloaded_bytes": size,
                    "high_z_retained_rows": retained,
                    "high_z_weight_sum": float(np.sum(target[3], dtype="f8")),
                }
                print("FINE_RR_INPUT_OK", cap, tracer, retained, sha, flush=True)
            except (OSError, ValueError, KeyError, MemoryError) as exc:
                errors.append(f"{cap}/{tracer}: {exc}")
                print("FINE_RR_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()
        if set(cats) != {"LRG", "ELG"}:
            continue
        try:
            lrg, elg = cats["LRG"], cats["ELG"]
            lrgsub, lseed = deterministic_subsets(
                lrg, ci, 0, HIGH_Z_INDEX)
            elgsub, eseed = deterministic_subsets(
                elg, ci, 1, HIGH_Z_INDEX)
            fine_rr = {}
            metadata = {}
            for label in THINNING:
                first = coordinates(lrgsub[label], PRIMARY_GEOMETRY)
                second = coordinates(elgsub[label], PRIMARY_GEOMETRY)
                raw, information = compile_rr(
                    first, second, fine, mu, args.threads)
                arrays[f"fine_rr_{cap}_{label}"] = raw
                arrays[f"fine_norm_{cap}_{label}"] = np.array(
                    [information["pair_weight_normalization"]])
                fine_rr[label] = raw
                metadata[label] = information
                print(f"FINE_RR_OK {cap} {label} "
                      f"positive_cells={information['positive_s_mu_cells']}"
                      f"/{raw.size}", flush=True)
            full = fine_rr["full"]
            full_norm = metadata["full"]["pair_weight_normalization"]
            compiled_coarse, coarse_meta = compile_rr(
                coordinates(lrg, PRIMARY_GEOMETRY),
                coordinates(elg, PRIMARY_GEOMETRY),
                coarse, mu, args.threads)
            rebinned = rebin_to_coarse(full, fine, coarse)
            residual = histogram_distance(
                rebinned / full_norm,
                compiled_coarse / coarse_meta["pair_weight_normalization"])
            l1_rebin = residual["normalized_rr_l1_over_larger_sample_l1"]
            if l1_rebin > 1e-8:
                raise ValueError(
                    f"Fine-to-coarse pair-conservation test failed: {l1_rebin}")
            # Check the published coarse pilot's total RR sum using its
            # retained, log-rounded value. Never use odd-sector results.
            published = next(
                x for x in full_ref["full_vs_deterministic_disjoint_halves_relative_l1"]
                if x["cap"] == cap and x["zlo"] == high_z[0])
            if not published:
                raise ValueError("Missing pre-existing coarse-bin evidence")
            thinning = {
                label: histogram_distance(
                    fine_rr[label] / metadata[label]["pair_weight_normalization"],
                    full / full_norm)
                for label in ("half_A", "half_B")
            }
            angular_full = odd_rr_geometry(
                rebin_to_coarse(full, fine, coarse), mu)
            angular_half = {
                label: odd_rr_geometry(
                    rebin_to_coarse(fine_rr[label], fine, coarse), mu)
                for label in ("half_A", "half_B")
            }
            odd_diagnostics = {}
            for ell in ("1", "3"):
                base = angular_full[ell]
                deviations = [
                    np.max(np.abs(angular_half[label][ell] - base))
                    for label in ("half_A", "half_B")
                ]
                odd_diagnostics[f"raw_rr_P{ell}"] = {
                    "coarse_bin_values_full": base.tolist(),
                    "largest_abs_half_minus_full": float(max(deviations)),
                    "no_half_A_vs_half_B_sign_flips": bool(np.all(
                        np.sign(angular_half["half_A"][ell]) ==
                        np.sign(angular_half["half_B"][ell]))),
                }
            secondary, secondary_meta = compile_rr(
                coordinates(lrg, SECONDARY_GEOMETRY),
                coordinates(elg, SECONDARY_GEOMETRY),
                fine, mu, args.threads)
            arrays[f"fine_rr_{cap}_secondary_fiducial"] = secondary
            arrays[f"fine_norm_{cap}_secondary_fiducial"] = np.array(
                [secondary_meta["pair_weight_normalization"]])
            alt_dist = histogram_distance(
                secondary / secondary_meta["pair_weight_normalization"],
                full / full_norm)
            case = {
                "cap": cap, "zlo": high_z[0], "zhi": high_z[1],
                "inputs": inputs,
                "subsampling_seeds": {"LRG": lseed, "ELG": eseed},
                "fine_edges_mpc_over_h": fine.tolist(),
                "fine_positive_cells": int(np.count_nonzero(full)),
                "fine_total_cells": int(full.size),
                "coarse_positive_cells": int(np.count_nonzero(compiled_coarse)),
                "coarse_total_cells": int(compiled_coarse.size),
                "rebin_coarse_relative_l1": l1_rebin,
                "full_fine_rr_normalization": full_norm,
                "fine_half_vs_full_relative_l1": thinning,
                "raw_rr_angular_moment_checks": odd_diagnostics,
                "secondary_primary_fine_relative_l1": alt_dist,
                "full_fine_rr_npz_key": f"fine_rr_{cap}_full",
            }
            cases.append(case)
            print("FINE_RR_REBIN_CLOSURE", cap, "relative_l1",
                  l1_rebin, "fine_positive", np.count_nonzero(full),
                  "of", full.size, flush=True)
            print("FINE_RR_HIGHZ_SAMPLING", cap, "halfA",
                  thinning["half_A"]["normalized_rr_l1_over_larger_sample_l1"],
                  "halfB",
                  thinning["half_B"]["normalized_rr_l1_over_larger_sample_l1"],
                  "secondary",
                  alt_dist["normalized_rr_l1_over_larger_sample_l1"],
                  flush=True)
        except (OSError, ValueError, RuntimeError, KeyError,
                MemoryError) as exc:
            errors.append(f"{cap}: {exc}")
            print("FINE_RR_ERROR", errors[-1], flush=True)
    complete = len(cases) == 2 and not errors
    report = {
        "study": "Observed eBOSS high-z random-only fine-RR geometric closure",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "fine_RR_high_z_geometric_closure_passed" if complete else "partial",
        "prior_full_coarse_rr_reference": str(FULL_REF.relative_to(ROOT)),
        "input_weight_reference": str(REFERENCE.relative_to(ROOT)),
        "caps": ["NGC", "SGC"],
        "candidate_z_interval": list(high_z),
        "candidate_interval_frozen": False,
        "fine_s_step_mpc_over_h": FINE_STEP,
        "coarse_s_edges_mpc_over_h": coarse.tolist(),
        "mu_bins": MU_BINS,
        "theta_min_deg": THETA_MIN_DEG,
        "weights": "audited product of four pure-eBOSS random columns",
        "primary_distance": DISTANCE_CONVENTIONS[PRIMARY_GEOMETRY],
        "secondary_distance": DISTANCE_CONVENTIONS[SECONDARY_GEOMETRY],
        "primary_distance_is_production_exact": False,
        "orientation": "LRG->ELG",
        "los": "midpoint",
        "cases": cases, "errors": errors,
        "observed_galaxies_read": False,
        "observed_odd_data_vector_read": False,
        "physical_even_to_odd_window_response_built": False,
        "all_candidate_z_bins_fine_RR_computed": False,
        "mock_RR_window_computed": False,
        "note": (
            "Fine RR and rebin closure are input-only geometric checks. "
            "Raw RR odd moments measure survey-pair asymmetry, not an "
            "observed galaxy multipole. Full four-bin fine counts, theory "
            "convolution and mock-window closure remain separate gates."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if arrays:
        np.savez_compressed(out.with_suffix(".npz"), **arrays)
    print("EBOSS_FINE_RR_HIGHZ_AUDIT", report["status"], out, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
