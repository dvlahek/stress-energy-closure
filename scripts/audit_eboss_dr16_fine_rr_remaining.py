#!/usr/bin/env python3
"""Full observed-random fine RR for the three remaining eBOSS candidate z bins.

Input-only continuation of the completed z=0.9-1.0 fine-RR audit. Each run
takes a PREDECLARED z-bin index 0, 1 or 2 and counts both NGC and SGC
LRG->ELG random-only RR(s,mu) on the same 1 Mpc/h x 240 signed-mu grid.
The primary full and disjoint half counts are independently rebinned to
the published 20 Mpc/h coarse grid. We independently recount coarse
pairs and compare to the SHA-pinned completed *prior workflow* RR array.
The previously declared secondary fiducial is computed on the same full
random input, never chosen from observed galaxy odd statistics.

Neither real galaxies, mock galaxy positions, odd vectors, wake templates,
nor a physical signal/covariance is read by this script. The candidate
z cuts are not the frozen prospective analysis protocol.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_fine_rr_highz import (
    fine_edges, rebin_to_coarse, odd_rr_geometry, unit_test as highz_unit_test,
)
from audit_eboss_dr16_full_random_rr import (
    BINS, MU_BINS, REFERENCE, ROOT, SEPARATION_EDGES, THETA_MIN_DEG,
    compile_rr, coordinates, deterministic_subsets, read_random_bins,
)
from audit_eboss_dr16_rr_density import histogram_distance
from eboss_dr16_fiducial import (
    DISTANCE_CONVENTIONS, PRIMARY_GEOMETRY, SECONDARY_GEOMETRY,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import download

COARSE_COMPACT = (
    ROOT / "source_data/eboss_dr16_full_observed_rr_2026-09-24.json"
)
PINNED_COARSE_SHA = "eab30356c5cab01cba8131e7c460bf8e32a5a90c"
REMAINING_BIN_INDICES = (0, 1, 2)
THINNING = ("full", "half_A", "half_B")
COARSE_REBIN_TOL = 1e-8
PRIMARY_PAIR_NORM_TOL = 1e-9


def compare_normalized(
    left: np.ndarray, leftnorm: float, right: np.ndarray, rightnorm: float,
) -> float:
    if (left.shape != right.shape or leftnorm <= 0 or rightnorm <= 0
            or not np.isfinite(leftnorm) or not np.isfinite(rightnorm)):
        raise ValueError("Incompatible RR shapes or pair normalizations")
    return float(histogram_distance(
        left / leftnorm, right / rightnorm,
    )["normalized_rr_l1_over_larger_sample_l1"])


def validate_prior_coarse(audit: dict, compact: dict,
                          npz: np.lib.npyio.NpzFile) -> None:
    if (audit.get("status") != "full_observed_random_candidate_bin_RR_computed"
            or audit.get("revision_commit") != PINNED_COARSE_SHA
            or audit.get("errors") != []
            or audit.get("observed_odd_vector_read") is not False
            or compact.get("source_commit") != PINNED_COARSE_SHA
            or compact.get("status") !=
            "full_observed_random_candidate_bin_RR_computed"
            or compact.get("all_8_full_histograms_have_all_1440_smu_cells_nonzero")
            is not True):
        raise ValueError("Published coarse RR provenance/selection gate missing")
    expected = {(x["cap"], int(x["bin_index"])) for x in audit["records"]}
    if expected != {(cap, iz) for cap in ("NGC", "SGC") for iz in range(4)}:
        raise ValueError("Prior coarse RR audit does not cover all candidate bins")
    if not np.array_equal(np.asarray(npz["sedges_mpc_over_h"]),
                          np.asarray(SEPARATION_EDGES)):
        raise ValueError("Previous coarse separation edges changed")
    if np.asarray(npz["muedges"]).shape != (MU_BINS + 1,):
        raise ValueError("Previous coarse signed mu edges changed")


def independently_compare_prior(
    raw_fine: np.ndarray, fine: np.ndarray, coarse: np.ndarray,
    cap: str, bin_index: int, expected: np.lib.npyio.NpzFile,
    pairnorm: float,
) -> dict:
    rebinned = rebin_to_coarse(raw_fine, fine, coarse)
    old_raw = np.asarray(expected[f"rr_{cap}_z{bin_index}_full"], dtype="f8")
    old_norm = float(
        np.asarray(expected[f"rrnorm_{cap}_z{bin_index}_full"]).item())
    relative = compare_normalized(rebinned, pairnorm, old_raw, old_norm)
    if relative > COARSE_REBIN_TOL:
        raise ValueError(
            f"Fine RR disagrees with previous SHA-pinned coarse RR: {relative}")
    norm_delta = abs(pairnorm / old_norm - 1)
    if norm_delta > PRIMARY_PAIR_NORM_TOL:
        raise ValueError("Fine RR normalization differs from previous full RR")
    if not np.all(old_raw > 0) or not np.all(rebinned > 0):
        raise ValueError("Previous or rebinned coarse RR has missing (s,mu) cells")
    return {
        "fine_rebin_vs_prior_coarse_normalized_relative_l1": relative,
        "prior_coarse_pairnorm_relative_difference": norm_delta,
        "positive_prior_coarse_cells": int(np.count_nonzero(old_raw)),
        "positive_rebinned_coarse_cells": int(np.count_nonzero(rebinned)),
        "coarse_total_cells": int(old_raw.size),
    }


def self_test() -> None:
    highz_unit_test()
    a = np.array([[1., 2.], [3., 4.]])
    assert compare_normalized(a, 10., 3 * a, 30.) < 1e-15
    try:
        compare_normalized(a, -1, a, 1)
    except ValueError:
        pass
    else:
        raise AssertionError("A nonpositive pair normalization must fail")
    assert REMAINING_BIN_INDICES == (0, 1, 2)
    print("EBOSS_REMAINING_FINE_RR_INDEPENDENT_SYNTHETIC_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bin-index", type=int, choices=REMAINING_BIN_INDICES)
    ap.add_argument("--coarse-audit")
    ap.add_argument("--coarse-rr")
    ap.add_argument("--out")
    ap.add_argument("--cache-dir", default="eboss_workspace/fine_remaining_randoms")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (args.bin_index is None or not args.coarse_audit or
            not args.coarse_rr or not args.out
            or args.threads < 1 or args.timeout <= 0):
        ap.error("Specify --bin-index, --coarse-audit, --coarse-rr, --out, "
                 "and valid runtime options")
    zindex = int(args.bin_index)
    lo, hi = BINS[zindex]
    input_ref = json.loads(REFERENCE.read_text())
    compact = json.loads(COARSE_COMPACT.read_text())
    prior = json.loads(Path(args.coarse_audit).read_text())
    coarse_counts = np.load(args.coarse_rr, allow_pickle=False)
    validate_prior_coarse(prior, compact, coarse_counts)
    if input_ref.get("status") != "sample_weighted_normalization_audited":
        raise ValueError("Four observed random input weights are not audited")
    file_records = {
        (v["cap"], v["tracer"]): v
        for v in input_ref["per_file_records"]
        if v["survey"] == "observed" and v["role"] == "random"
    }
    if len(file_records) != 4:
        raise ValueError("Expected exactly four SHA-pinned observed random files")
    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    if not np.array_equal(mu, np.asarray(coarse_counts["muedges"])):
        raise ValueError("Fine and pinned coarse signed-mu binning differs")
    arrays = {
        "fine_sedges_mpc_over_h": fine,
        "coarse_sedges_mpc_over_h": coarse,
        "muedges": mu,
        "zlo": np.asarray([lo], dtype="f8"),
        "zhi": np.asarray([hi], dtype="f8"),
    }
    cases, errors = [], []
    root = Path(args.cache_dir) / f"z{zindex}"
    for ci, cap in enumerate(("NGC", "SGC")):
        catalogues, sources = {}, {}
        for tracer in ("LRG", "ELG"):
            filename, header_rows = RANDOMS[(tracer, cap)]
            path = root / filename
            try:
                path, sha, nbytes = download(
                    filename, root, 1024 * 1024 * 1024, args.timeout)
                if sha != file_records[(cap, tracer)]["sha256"]:
                    raise ValueError("Public observed random FITS SHA256 changed")
                all_bins, _ = read_random_bins(
                    path, tracer, cap, header_rows)
                selected = all_bins[zindex]
                expected = next(
                    x for x in compact["input_random_redshift_counts"]
                    if x["cap"] == cap and x["tracer"] == tracer
                    and x["zlo"] == lo and x["zhi"] == hi
                )
                if len(selected[0]) != expected["retained_rows"]:
                    raise ValueError("Candidate-bin random count changed")
                catalogues[tracer] = selected
                sources[tracer] = {
                    "filename": filename,
                    "full_sha256": sha, "full_bytes": nbytes,
                    "candidate_retained_rows": len(selected[0]),
                    "candidate_weight_sum": float(
                        np.sum(selected[3], dtype="f8")),
                }
                print("REMAINING_FINE_RR_INPUT_OK", zindex, cap, tracer,
                      len(selected[0]), sha, flush=True)
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError, StopIteration) as exc:
                errors.append(f"z{zindex}/{cap}/{tracer}: {exc}")
                print("REMAINING_FINE_RR_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()
        if set(catalogues) != {"LRG", "ELG"}:
            continue
        try:
            lrg = catalogues["LRG"]
            elg = catalogues["ELG"]
            lrg_thin, lrg_seed = deterministic_subsets(
                lrg, ci, 0, zindex)
            elg_thin, elg_seed = deterministic_subsets(
                elg, ci, 1, zindex)
            rr, infos = {}, {}
            for label in THINNING:
                first = coordinates(lrg_thin[label], PRIMARY_GEOMETRY)
                second = coordinates(elg_thin[label], PRIMARY_GEOMETRY)
                raw, info = compile_rr(
                    first, second, fine, mu, args.threads)
                arrays[f"fine_rr_{cap}_z{zindex}_{label}"] = raw
                arrays[f"fine_norm_{cap}_z{zindex}_{label}"] = np.asarray(
                    [info["pair_weight_normalization"]], dtype="f8")
                rr[label], infos[label] = raw, info
                print("REMAINING_FINE_RR_COUNT_OK", zindex, cap, label,
                      "positive_cells", info["positive_s_mu_cells"],
                      "of", raw.size, flush=True)
            full = rr["full"]
            full_norm = infos["full"]["pair_weight_normalization"]
            fresh_coarse, coarse_info = compile_rr(
                coordinates(lrg, PRIMARY_GEOMETRY),
                coordinates(elg, PRIMARY_GEOMETRY),
                coarse, mu, args.threads)
            binned = rebin_to_coarse(full, fine, coarse)
            independent = compare_normalized(
                binned, full_norm,
                fresh_coarse, coarse_info["pair_weight_normalization"])
            if independent > COARSE_REBIN_TOL:
                raise ValueError(
                    f"Direct independent coarse RR differs after rebin: {independent}")
            old_closure = independently_compare_prior(
                full, fine, coarse, cap, zindex, coarse_counts, full_norm)
            thinning = {
                name: histogram_distance(
                    rr[name] / infos[name]["pair_weight_normalization"],
                    full / full_norm,
                )
                for name in ("half_A", "half_B")
            }
            angular_full = odd_rr_geometry(binned, mu)
            angular_half = {
                name: odd_rr_geometry(
                    rebin_to_coarse(rr[name], fine, coarse), mu)
                for name in ("half_A", "half_B")
            }
            angular = {
                f"raw_rr_P{ell}": {
                    "coarse_bin_values_full": angular_full[ell].tolist(),
                    "max_abs_half_minus_full": max(
                        float(np.max(np.abs(
                            angular_half[name][ell] - angular_full[ell])))
                        for name in ("half_A", "half_B")),
                }
                for ell in ("1", "3")
            }
            other, other_info = compile_rr(
                coordinates(lrg, SECONDARY_GEOMETRY),
                coordinates(elg, SECONDARY_GEOMETRY),
                fine, mu, args.threads)
            arrays[f"fine_rr_{cap}_z{zindex}_secondary_fiducial"] = other
            arrays[f"fine_norm_{cap}_z{zindex}_secondary_fiducial"] = (
                np.asarray([other_info["pair_weight_normalization"]]))
            alternate = histogram_distance(
                other / other_info["pair_weight_normalization"],
                full / full_norm,
            )
            case = {
                "cap": cap, "bin_index": zindex,
                "zlo": lo, "zhi": hi,
                "input_random_file_evidence": sources,
                "subsampling_seeds": {
                    "LRG": lrg_seed, "ELG": elg_seed,
                },
                "full_fine_weighted_rr_metadata": infos["full"],
                "full_fine_positive_cells": int(np.count_nonzero(full)),
                "full_fine_total_cells": int(full.size),
                "full_fine_coarse_aggregated_positive_cells":
                    int(np.count_nonzero(binned)),
                "full_fine_coarse_aggregated_total_cells": int(binned.size),
                "fine_to_independent_fresh_coarse_relative_l1": independent,
                "fine_to_prior_pinned_coarse": old_closure,
                "fine_disjoint_halves_vs_full_normalized_relative_l1": thinning,
                "raw_rr_odd_geometry_only": angular,
                "predeclared_secondary_vs_primary_fine_rr": alternate,
                "full_fine_rr_npz_key": f"fine_rr_{cap}_z{zindex}_full",
            }
            cases.append(case)
            print("REMAINING_FINE_RR_CLOSURE", zindex, cap,
                  "direct_rebin_l1", independent,
                  "prior_coarse_l1",
                  old_closure["fine_rebin_vs_prior_coarse_normalized_relative_l1"],
                  "halfA_l1",
                  thinning["half_A"]["normalized_rr_l1_over_larger_sample_l1"],
                  "halfB_l1",
                  thinning["half_B"]["normalized_rr_l1_over_larger_sample_l1"],
                  "alt_geometry_l1",
                  alternate["normalized_rr_l1_over_larger_sample_l1"],
                  flush=True)
        except (OSError, ValueError, RuntimeError, KeyError,
                MemoryError, StopIteration) as exc:
            errors.append(f"z{zindex}/{cap}: {exc}")
            print("REMAINING_FINE_RR_CASE_ERROR", errors[-1], flush=True)
    completed = len(cases) == 2 and not errors
    audit = {
        "study": "Full observed eBOSS random-only fine RR for one remaining predeclared candidate z slice",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": (
            "remaining_candidate_bin_fine_RR_geometric_closure_passed"
            if completed else "partial"),
        "bin_index": zindex,
        "candidate_zlo": lo, "candidate_zhi": hi,
        "candidate_selection_frozen_for_inference": False,
        "pinned_prior_coarse_rr_workflow_sha": PINNED_COARSE_SHA,
        "pinned_prior_coarse_rr_status": prior["status"],
        "pinned_random_weight_reference":
            str(REFERENCE.relative_to(ROOT)),
        "caps": ["NGC", "SGC"],
        "fine_separation_step_mpc_over_h": 1.0,
        "fine_s_edges_mpc_over_h": fine.tolist(),
        "coarse_s_edges_mpc_over_h": coarse.tolist(),
        "signed_mu_edges": mu.tolist(),
        "theta_min_deg": THETA_MIN_DEG,
        "orientation": "LRG->ELG", "los": "midpoint",
        "four_factor_weight_product": "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP",
        "primary_fiducial": DISTANCE_CONVENTIONS[PRIMARY_GEOMETRY],
        "secondary_fiducial": DISTANCE_CONVENTIONS[SECONDARY_GEOMETRY],
        "production_fiducial_distance_mapping_certified": False,
        "cases": cases, "errors": errors,
        "observed_galaxy_positions_read": False,
        "observed_odd_data_vector_read": False,
        "observed_DD_or_DR_computed": False,
        "wake_template_fit_or_tuning": False,
        "statistical_covariance_estimated": False,
        "fully_validated_physical_survey_window_built": False,
        "note": (
            "The full and half-sample fine RR are independent random-only "
            "input geometry checks for candidate bins. Full physical "
            "window, common mask, production fiducial and mock ensemble "
            "remain outside this certified scope."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    if arrays:
        np.savez_compressed(out.with_suffix(".npz"), **arrays)
    print("EBOSS_REMAINING_FINE_RR_AUDIT", audit["status"],
          "z_index", zindex, "cases", len(cases), out, flush=True)
    return 0 if completed else 2


if __name__ == "__main__":
    raise SystemExit(main())
