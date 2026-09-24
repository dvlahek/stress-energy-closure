#!/usr/bin/env python3
"""Full four-bin one-Mpc/h eBOSS observed-random RR geometry audit.

This is the all-candidate-bin extension of the already validated high-z
fine-RR audit. It opens only the four SHA-pinned public eBOSS DR16 random
catalogues, never observed galaxy data. For each cap and each predeclared
candidate redshift bin 0.6-0.7, 0.7-0.8, 0.8-0.9, 0.9-1.0 it counts the
complete weighted LRG->ELG RR(s,mu) histogram on a 1 Mpc/h separation grid
with 240 signed-mu bins and the fixed midpoint LOS / 0.05 degree cut.

For every cap/z case, two deterministic disjoint half-samples are counted
with the same geometry and their own pair-weight normalization. The fine
histogram is rebinned to the coarse 20 Mpc/h grid and compared with a fresh
independent compiled coarse count on the identical objects. This certifies
finite-bin geometry and random sampling sensitivity only; it is not an
observed odd-sector measurement, a covariance, or a physical wake fit.
"""
from __future__ import annotations

import argparse
import gc
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_fine_rr_highz import (
    fine_edges, rebin_to_coarse, odd_rr_geometry,
)
from audit_eboss_dr16_full_random_rr import (
    ROOT, REFERENCE, BINS, SEPARATION_EDGES, MU_BINS,
    read_random_bins, coordinates, compile_rr, deterministic_subsets,
)
from audit_eboss_dr16_rr_density import histogram_distance
from eboss_dr16_fiducial import PRIMARY_GEOMETRY
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import download

FULL_REF = ROOT / "source_data/eboss_dr16_full_observed_rr_2026-09-24.json"
HIGHZ_REF = ROOT / "source_data/eboss_dr16_fine_rr_highz_2026-09-24.json"
THINNING = ("full", "half_A", "half_B")


def self_test() -> None:
    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    rng = np.random.default_rng(24683)

    def cat(n, lo, hi):
        return (
            rng.uniform(140, 145, n),
            rng.uniform(5, 12, n),
            rng.uniform(lo + 0.01, hi - 0.01, n),
            rng.uniform(0.6, 1.7, n),
        )

    for iz, (lo, hi) in enumerate(BINS):
        lrg, elg = cat(83 + iz, lo, hi), cat(91 + iz, lo, hi)
        lc, ec = coordinates(lrg, PRIMARY_GEOMETRY), coordinates(
            elg, PRIMARY_GEOMETRY)
        fine_rr, fine_meta = compile_rr(lc, ec, fine, mu, 2)
        coarse_rr, coarse_meta = compile_rr(lc, ec, coarse, mu, 2)
        np.testing.assert_allclose(
            rebin_to_coarse(fine_rr, fine, coarse),
            coarse_rr, atol=2e-7, rtol=2e-11)
        assert abs(
            fine_meta["pair_weight_normalization"]
            / coarse_meta["pair_weight_normalization"] - 1
        ) < 1e-12
    print("EBOSS_ALLBIN_FINE_RR_SYNTHETIC_REBIN_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out", default="eboss_workspace/fine_rr_allbins_audit.json")
    ap.add_argument(
        "--cache-dir", default="eboss_workspace/fine_rr_allbins_randoms")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
        return 0
    if args.threads < 1 or args.timeout <= 0:
        ap.error("Invalid thread count or timeout")

    input_ref = json.loads(REFERENCE.read_text())
    coarse_ref = json.loads(FULL_REF.read_text())
    highz_ref = json.loads(HIGHZ_REF.read_text())
    if (input_ref["status"] != "sample_weighted_normalization_audited"
            or coarse_ref["status"]
            != "full_observed_random_candidate_bin_RR_computed"
            or highz_ref["status"] != "fine_RR_high_z_geometric_closure_passed"):
        raise ValueError("Required blinded RR provenance gates are incomplete")

    known = {
        (x["cap"], x["tracer"]): x
        for x in input_ref["per_file_records"]
        if x["survey"] == "observed" and x["role"] == "random"
    }
    retained_ref = {
        (x["cap"], x["tracer"], float(x["zlo"]), float(x["zhi"])):
        int(x["retained_rows"])
        for x in coarse_ref["input_random_redshift_counts"]
    }
    if len(known) != 4 or len(retained_ref) != 16:
        raise ValueError("Incomplete retained observed-random provenance")

    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    arrays = {
        "fine_sedges_mpc_over_h": fine,
        "coarse_sedges_mpc_over_h": coarse,
        "muedges": mu,
        "zlo": np.asarray([x[0] for x in BINS], dtype="f8"),
        "zhi": np.asarray([x[1] for x in BINS], dtype="f8"),
    }
    cases, input_files, errors = [], [], []
    cache = Path(args.cache_dir)

    for ci, cap in enumerate(("NGC", "SGC")):
        by_tracer = {}
        for tracer in ("LRG", "ELG"):
            filename, expected_rows = RANDOMS[(tracer, cap)]
            path = cache / filename
            try:
                path, sha, size = download(
                    filename, cache, 1024 * 1024 * 1024, args.timeout)
                if sha != known[(cap, tracer)]["sha256"]:
                    raise ValueError("Observed random FITS SHA256 changed")
                all_bins, summary = read_random_bins(
                    path, tracer, cap, expected_rows)
                for iz, (lo, hi) in enumerate(BINS):
                    expected = retained_ref[(cap, tracer, lo, hi)]
                    if len(all_bins[iz][0]) != expected:
                        raise ValueError(
                            f"Retained row count changed in {cap}/{tracer}/{lo}-{hi}")
                by_tracer[tracer] = all_bins
                input_files.append({
                    "cap": cap, "tracer": tracer, "filename": filename,
                    "sha256": sha, "downloaded_bytes": size,
                    "retained_rows_per_candidate_bin":
                        summary["retained_per_candidate_bin"],
                    "weighted_sum_per_candidate_bin":
                        summary["weighted_sum_per_candidate_bin"],
                })
                print("ALLBIN_FINE_INPUT_OK", cap, tracer, sha, flush=True)
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError) as exc:
                errors.append(f"{cap}/{tracer}: {exc}")
                print("ALLBIN_FINE_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()

        if set(by_tracer) != {"LRG", "ELG"}:
            continue

        for iz, (lo, hi) in enumerate(BINS):
            try:
                lrg, elg = by_tracer["LRG"][iz], by_tracer["ELG"][iz]
                lsubs, lseed = deterministic_subsets(lrg, ci, 0, iz)
                esubs, eseed = deterministic_subsets(elg, ci, 1, iz)

                # Convert the complete catalogues once. Reproduce exact
                # deterministic subset indices used by deterministic_subsets
                # on the converted RDD arrays to avoid repeated cosmology work.
                lc, ec = coordinates(lrg, PRIMARY_GEOMETRY), coordinates(
                    elg, PRIMARY_GEOMETRY)
                converted = {}
                for tracer_id, fullcat, rawcat, label in (
                    (0, lc, lrg, "LRG"), (1, ec, elg, "ELG")
                ):
                    _, seed = deterministic_subsets(
                        rawcat, ci, tracer_id, iz)
                    order = np.random.default_rng(seed).permutation(
                        len(fullcat[0]))
                    half = len(order) // 2
                    converted[label] = {
                        "full": fullcat,
                        "half_A": tuple(v[order[:half]] for v in fullcat),
                        "half_B": tuple(v[order[half:]] for v in fullcat),
                    }

                result = {}
                normalized = {}
                for sample in THINNING:
                    rr, meta = compile_rr(
                        converted["LRG"][sample],
                        converted["ELG"][sample],
                        fine, mu, args.threads)
                    norm = meta["pair_weight_normalization"]
                    normalized[sample] = rr / norm
                    arrays[f"fine_rr_{cap}_z{iz}_{sample}"] = rr
                    arrays[f"fine_norm_{cap}_z{iz}_{sample}"] = np.asarray(
                        [norm], dtype="f8")
                    result[sample] = {
                        "lrg_rows": int(len(converted["LRG"][sample][0])),
                        "elg_rows": int(len(converted["ELG"][sample][0])),
                        "positive_cells": int(np.count_nonzero(rr)),
                        "total_cells": int(rr.size),
                        "pair_weight_normalization": float(norm),
                        "normalized_rr_fraction_in_candidate_s_bins":
                            meta["normalized_rr_weight_fraction_in_test_bins"],
                        "raw_rr_odd_geometry_coarse":
                            odd_rr_geometry(
                                rebin_to_coarse(rr, fine, coarse), mu),
                        "npz_rr_key": f"fine_rr_{cap}_z{iz}_{sample}",
                    }
                    print(
                        "ALLBIN_FINE_RR_OK", cap, f"{lo:.1f}-{hi:.1f}",
                        sample, "positive", np.count_nonzero(rr), "of",
                        rr.size, flush=True)

                direct_coarse, direct_meta = compile_rr(
                    lc, ec, coarse, mu, args.threads)
                fine_rebinned = rebin_to_coarse(
                    arrays[f"fine_rr_{cap}_z{iz}_full"], fine, coarse)
                fine_norm = result["full"]["pair_weight_normalization"]
                direct_norm = direct_meta["pair_weight_normalization"]
                closure = histogram_distance(
                    fine_rebinned / fine_norm,
                    direct_coarse / direct_norm)
                rebin_l1 = closure[
                    "normalized_rr_l1_over_larger_sample_l1"]
                if rebin_l1 > 1e-8:
                    raise ValueError(
                        f"Fine-to-coarse closure failed: {rebin_l1:g}")

                half_checks = {
                    sample: histogram_distance(
                        normalized[sample], normalized["full"])
                    for sample in ("half_A", "half_B")
                }
                half_ab = histogram_distance(
                    normalized["half_A"], normalized["half_B"])

                # Cross-check the already certified high-z result numerically,
                # but do not use it to choose settings or discard a case.
                previous_highz = None
                if iz == 3:
                    old = next(x for x in highz_ref["cap_results"]
                               if x["cap"] == cap)
                    previous_highz = {
                        "previous_rebin_l1":
                            old["rebin_fine_to_coarse_relative_l1"],
                        "new_minus_previous_rebin_l1": float(
                            rebin_l1
                            - old["rebin_fine_to_coarse_relative_l1"]),
                    }
                    if abs(
                        result["full"]["lrg_rows"]
                        - highz_ref[
                            "full_observed_random_catalogue_counts_from_prior_full_rr"
                        ][cap]["LRG"]
                    ) > 0:
                        raise ValueError("High-z LRG input count changed")

                cases.append({
                    "cap": cap, "bin_index": iz,
                    "zlo": lo, "zhi": hi,
                    "lrg_full_rows": len(lrg[0]),
                    "elg_full_rows": len(elg[0]),
                    "subsampling_seeds": {"LRG": lseed, "ELG": eseed},
                    "samples": result,
                    "fine_to_independent_coarse": closure,
                    "half_A_vs_full": half_checks["half_A"],
                    "half_B_vs_full": half_checks["half_B"],
                    "half_A_vs_half_B": half_ab,
                    "high_z_prior_crosscheck": previous_highz,
                })
                print(
                    "ALLBIN_FINE_CLOSURE", cap, f"{lo:.1f}-{hi:.1f}",
                    "rebin_l1", rebin_l1,
                    "halfA_l1",
                    half_checks["half_A"][
                        "normalized_rr_l1_over_larger_sample_l1"],
                    "halfB_l1",
                    half_checks["half_B"][
                        "normalized_rr_l1_over_larger_sample_l1"],
                    flush=True)
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError) as exc:
                errors.append(f"{cap}/z{iz}: {exc}")
                print("ALLBIN_FINE_RR_ERROR", errors[-1], flush=True)

        del by_tracer
        gc.collect()

    complete = len(cases) == 8 and not errors
    report = {
        "study": "Full four-bin observed eBOSS random-only fine RR geometry",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": (
            "full_four_bin_fine_RR_geometric_closure_passed"
            if complete else "partial"),
        "candidate_redshift_bins": [list(x) for x in BINS],
        "candidate_bins_frozen_for_inference": False,
        "caps": ["NGC", "SGC"],
        "fine_separation_step_mpc_over_h": 1.0,
        "fine_separation_edges_mpc_over_h": fine.tolist(),
        "coarse_separation_edges_mpc_over_h": coarse.tolist(),
        "mu_bins": MU_BINS,
        "orientation": "LRG->ELG",
        "line_of_sight": "midpoint",
        "distance_geometry": PRIMARY_GEOMETRY,
        "input_file_evidence": input_files,
        "cases": cases,
        "errors": errors,
        "uses_observed_galaxy_positions": False,
        "uses_observed_odd_data_vector": False,
        "uses_mock_galaxy_positions": False,
        "wake_template_fit_or_tuning": False,
        "physical_window_response_built": False,
        "statistical_covariance_built": False,
        "scope_note": (
            "All four redshift slices remain candidate bins. This workflow "
            "certifies fine observed-random RR support, independent fine-to-"
            "coarse pair conservation, and deterministic half-sample "
            "sensitivity. It does not freeze the analysis protocol or "
            "measure a galaxy odd multipole."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(out.with_suffix(".npz"), **arrays)
    print("EBOSS_ALLBIN_FINE_RR_AUDIT", report["status"], out, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
