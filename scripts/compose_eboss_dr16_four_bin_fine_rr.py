#!/usr/bin/env python3
"""Merge the four SHA-pinned, random-only eBOSS candidate fine RR z slices.

The three independent lower-z matrix jobs are joined only when every cap
and candidate interval passes the direct-coarse and previous-workflow
coarse RR comparison. The fourth slice is the already validated high-z
fine-RR source. Output preserves *separate* NGC and SGC weighted RR
counts and separate RR pair-weight norms; it does not silently average
or mix caps and does not read galaxy, odd-sector or wake data.

This is a complete random-only fine-pair grid over the *candidate* bins,
not a frozen physical eBOSS selection or fully calibrated survey window.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_fine_rr_remaining import (
    COARSE_COMPACT, PINNED_COARSE_SHA,
    compare_normalized, independently_compare_prior, validate_prior_coarse,
)
from audit_eboss_dr16_fine_rr_highz import rebin_to_coarse
from audit_eboss_dr16_full_random_rr import (
    BINS, MU_BINS, SEPARATION_EDGES,
)

CAPS = ("NGC", "SGC")
SAMPLES = ("full", "half_A", "half_B")
REMAINING_SHA = "d64fdfd5b2ab6aa00ad964b5fc0de961ff1acbb3"
HIGHZ_SHA = "5dbfca45c25c7af57bf05693022bfa98b499209d"


def input_paths(root: Path, index: int) -> tuple[Path, Path]:
    base = root / f"z{index}"
    return (
        base / f"fine_rr_z{index}_audit.json",
        base / f"fine_rr_z{index}_audit.npz",
    )


def check_case_geometry(
    rr: np.ndarray, fine: np.ndarray, coarse: np.ndarray,
    mu: np.ndarray,
) -> dict:
    if rr.shape != (120, MU_BINS):
        raise ValueError("Fine RR array must have 120 separation x 240 mu cells")
    if not np.isfinite(rr).all() or np.any(rr < 0):
        raise ValueError("Fine RR has invalid weighted pair counts")
    binned = rebin_to_coarse(rr, fine, coarse)
    if not np.all(binned > 0):
        raise ValueError(
            "A coarse (s,mu) cell has no random support for the LS estimator")
    return {
        "positive_fine_cells": int(np.count_nonzero(rr)),
        "total_fine_cells": int(rr.size),
        "positive_coarse_cells": int(np.count_nonzero(binned)),
        "total_coarse_cells": int(binned.size),
        "all_coarse_mu_cells_supported": True,
        "minimum_positive_coarse_cell_weight": float(np.min(binned)),
    }


def synthetic_self_test() -> None:
    fine = np.linspace(20, 140, 121)
    coarse = np.asarray(SEPARATION_EDGES)
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    rr = np.ones((120, MU_BINS), dtype="f8")
    counts = check_case_geometry(rr, fine, coarse, mu)
    assert counts["positive_fine_cells"] == 120 * MU_BINS
    assert counts["positive_coarse_cells"] == 6 * MU_BINS
    assert compare_normalized(rr, 2, 5 * rr, 10) < 1e-15
    rr[0, 0] = 0.0
    assert check_case_geometry(rr, fine, coarse, mu)["positive_coarse_cells"] == 1440
    rr[:20, 0] = 0.0
    try:
        check_case_geometry(rr, fine, coarse, mu)
    except ValueError:
        pass
    else:
        raise AssertionError("Missing coarse signed-mu support must fail")
    print("EBOSS_FOUR_BIN_FINE_RR_COMPOSITION_SELF_TEST_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--remaining-root")
    ap.add_argument("--highz-audit")
    ap.add_argument("--highz-rr")
    ap.add_argument("--coarse-audit")
    ap.add_argument("--coarse-rr")
    ap.add_argument("--out", default="eboss_workspace/four_bin_fine_rr_audit.json")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        synthetic_self_test()
        return 0
    if not all((args.remaining_root, args.highz_audit,
                args.highz_rr, args.coarse_audit, args.coarse_rr)):
        ap.error("Provide the three independent lower-z, high-z and prior coarse RR artifacts")
    prior = json.loads(Path(args.coarse_audit).read_text())
    compact = json.loads(COARSE_COMPACT.read_text())
    old = np.load(args.coarse_rr, allow_pickle=False)
    validate_prior_coarse(prior, compact, old)
    highz = json.loads(Path(args.highz_audit).read_text())
    if (highz.get("status") != "fine_RR_high_z_geometric_closure_passed"
            or highz.get("revision_commit") != HIGHZ_SHA
            or highz.get("observed_odd_data_vector_read") is not False
            or highz.get("errors") != []
            or tuple(highz["candidate_z_interval"]) != BINS[3]):
        raise ValueError("Previously validated high-z fine-RR source mismatch")
    highz_npz = np.load(args.highz_rr, allow_pickle=False)
    lower = {}
    for iz in range(3):
        jpath, npath = input_paths(Path(args.remaining_root), iz)
        audit = json.loads(jpath.read_text())
        if (audit.get("status") !=
                "remaining_candidate_bin_fine_RR_geometric_closure_passed"
                or audit.get("revision_commit") != REMAINING_SHA
                or audit.get("pinned_prior_coarse_rr_workflow_sha")
                != PINNED_COARSE_SHA
                or audit.get("errors") != []
                or audit.get("observed_odd_data_vector_read") is not False
                or audit.get("observed_galaxy_positions_read") is not False
                or audit.get("bin_index") != iz
                or (audit.get("candidate_zlo"), audit.get("candidate_zhi"))
                != BINS[iz]):
            raise ValueError(f"Input provenance for lower candidate z{iz} invalid")
        cases = [v["cap"] for v in audit["cases"]]
        if len(cases) != 2 or set(cases) != set(CAPS):
            raise ValueError(f"Lower z{iz} does not cover two complete caps")
        lower[iz] = (audit, np.load(npath, allow_pickle=False))
    fine = np.asarray(highz_npz["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(highz_npz["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(highz_npz["muedges"], dtype="f8")
    if (fine.shape != (121,) or coarse.shape != (7,)
            or mu.shape != (MU_BINS + 1,)
            or not np.array_equal(coarse, np.asarray(SEPARATION_EDGES))
            or not np.allclose(np.diff(fine), 1, rtol=0, atol=1e-12)
            or not np.array_equal(mu, np.asarray(old["muedges"]))):
        raise ValueError("Fine observed high-z geometry differs from the prior coarse grid")
    for iz, (_, source) in lower.items():
        if not all(np.array_equal(
                source[key], match) for key, match in (
                    ("fine_sedges_mpc_over_h", fine),
                    ("coarse_sedges_mpc_over_h", coarse),
                    ("muedges", mu))):
            raise ValueError(f"z{iz} grid differs from pinned high-z grid")

    arrays = {
        "zlo": np.asarray([x[0] for x in BINS], dtype="f8"),
        "zhi": np.asarray([x[1] for x in BINS], dtype="f8"),
        "fine_sedges_mpc_over_h": fine,
        "coarse_sedges_mpc_over_h": coarse,
        "muedges": mu,
    }
    summary = []
    for cap in CAPS:
        gathered = {label: [] for label in SAMPLES}
        norms = {label: [] for label in SAMPLES}
        alternates, alternate_norms = [], []
        for iz in range(4):
            is_high = iz == 3
            source = highz_npz if is_high else lower[iz][1]
            record = highz if is_high else lower[iz][0]
            matches = [x for x in record["cases"] if x["cap"] == cap]
            if len(matches) != 1:
                raise ValueError(f"Missing unique {cap} z{iz} audit case")
            case = matches[0]
            if (case.get("zlo"), case.get("zhi")) != BINS[iz]:
                raise ValueError("Incorrect audit redshift interval")
            for label in SAMPLES:
                key = (f"fine_rr_{cap}_{label}" if is_high
                       else f"fine_rr_{cap}_z{iz}_{label}")
                norm_key = (f"fine_norm_{cap}_{label}" if is_high
                            else f"fine_norm_{cap}_z{iz}_{label}")
                raw = np.asarray(source[key], dtype="f8")
                norm = float(np.asarray(source[norm_key]).item())
                if norm <= 0 or not np.isfinite(norm):
                    raise ValueError("Fine RR normalization invalid")
                geom = check_case_geometry(raw, fine, coarse, mu)
                gathered[label].append(raw)
                norms[label].append(norm)
                if label == "full":
                    closure = independently_compare_prior(
                        raw, fine, coarse, cap, iz, old, norm)
                    previous = (case[
                        "fine_to_prior_pinned_coarse"][
                            "fine_rebin_vs_prior_coarse_normalized_relative_l1"]
                        if not is_high else None)
                    if previous is not None and abs(
                        previous - closure[
                            "fine_rebin_vs_prior_coarse_normalized_relative_l1"]
                    ) > 1e-10:
                        raise ValueError(
                            "New combined fine RR differs from per-bin prior closure")
                    summary.append({
                        "cap": cap, "bin_index": iz,
                        "zlo": BINS[iz][0], "zhi": BINS[iz][1],
                        "full_fine": geom,
                        "comparison_to_prior_full_coarse": closure,
                        "disjoint_half_A_vs_full_relative_l1": None,
                        "disjoint_half_B_vs_full_relative_l1": None,
                        "secondary_vs_primary_relative_l1": None,
                    })
            full, full_norm = gathered["full"][-1], norms["full"][-1]
            for label, key in (
                ("half_A", "disjoint_half_A_vs_full_relative_l1"),
                ("half_B", "disjoint_half_B_vs_full_relative_l1"),
            ):
                relative = compare_normalized(
                    gathered[label][-1], norms[label][-1],
                    full, full_norm)
                summary[-1][key] = relative
            alternative_key = (
                f"fine_rr_{cap}_secondary_fiducial"
                if is_high else
                f"fine_rr_{cap}_z{iz}_secondary_fiducial")
            alternative_norm_key = (
                f"fine_norm_{cap}_secondary_fiducial"
                if is_high else
                f"fine_norm_{cap}_z{iz}_secondary_fiducial")
            alternative = np.asarray(source[alternative_key], dtype="f8")
            alternative_norm = float(
                np.asarray(source[alternative_norm_key]).item())
            check_case_geometry(alternative, fine, coarse, mu)
            summary[-1]["secondary_vs_primary_relative_l1"] = (
                compare_normalized(
                    alternative, alternative_norm, full, full_norm))
            alternates.append(alternative)
            alternate_norms.append(alternative_norm)
            print("EBOSS_COMPOSED_FINE_RR", cap, iz,
                  "prior_coarse_l1",
                  summary[-1]["comparison_to_prior_full_coarse"][
                      "fine_rebin_vs_prior_coarse_normalized_relative_l1"],
                  "half_A_l1",
                  summary[-1]["disjoint_half_A_vs_full_relative_l1"],
                  "half_B_l1",
                  summary[-1]["disjoint_half_B_vs_full_relative_l1"],
                  "secondary_l1",
                  summary[-1]["secondary_vs_primary_relative_l1"],
                  "positive_fine_cells", summary[-1]["full_fine"]["positive_fine_cells"],
                  flush=True)
        for label in SAMPLES:
            arrays[f"rr_wcounts_{cap}_{label}"] = np.asarray(
                gathered[label], dtype="f8")
            arrays[f"rr_wnorm_{cap}_{label}"] = np.asarray(
                norms[label], dtype="f8")
        arrays[f"rr_wcounts_{cap}_secondary"] = np.asarray(
            alternates, dtype="f8")
        arrays[f"rr_wnorm_{cap}_secondary"] = np.asarray(
            alternate_norms, dtype="f8")
    if len(summary) != 8:
        raise ValueError("Full fine RR composition requires eight cap/z cases")
    report = {
        "study": "All four candidate eBOSS LRG->ELG observed-random fine RR window inputs",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "full_four_candidate_bin_observed_random_fine_RR_geometric_closure_passed",
        "candidate_redshift_bins": [list(x) for x in BINS],
        "candidate_selection_frozen_for_inference": False,
        "coarse_source_workflow_sha": PINNED_COARSE_SHA,
        "remaining_lower_z_source_workflow_sha": REMAINING_SHA,
        "highz_source_workflow_sha": HIGHZ_SHA,
        "caps_stored_separately": list(CAPS),
        "caps_combined_for_estimator": False,
        "fine_separation_step_mpc_over_h": 1.0,
        "coarse_separation_edges_mpc_over_h": coarse.tolist(),
        "signed_mu_edges": mu.tolist(),
        "orientation": "LRG->ELG",
        "line_of_sight": "midpoint",
        "weight_normalization": (
            "Per cap and per redshift bin separately: "
            "sum(w_LRG) * sum(w_ELG); each disjoint half uses its own "
            "sum of audited four-factor weights"),
        "cases": summary,
        "observed_galaxy_data_read": False,
        "observed_galaxy_pairs_computed": False,
        "observed_odd_data_vector_read": False,
        "wake_template_evaluated_or_tuned": False,
        "physical_even_to_odd_forward_model_certified": False,
        "realistic_mock_ensemble_covariance_computed": False,
        "common_footprint_and_production_distance_frozen": False,
        "note": (
            "These are full random-only fine RR candidate-window inputs with "
            "two independent coarse closure checks for the lower bins and "
            "an additional pinned coarse check for the high-z bin. "
            "The per-cap RR arrays are intentionally not combined. "
            "No observed odd-sector data or statistical fit was opened."),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(target.with_suffix(".npz"), **arrays)
    print("EBOSS_FOUR_CANDIDATE_BIN_FINE_RR_AUDIT",
          report["status"], target, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
