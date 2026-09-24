#!/usr/bin/env python3
"""Aggregate fixed nine-mock eBOSS random-only window and joint-sky pilot.

Requires all 18 per-cap matched mock RANDOM artifacts from the predeclared
nine-ID matrix, plus the separately SHA-pinned successful four-bin observed
RANDOM fine RR artifact. Independently downloads only public observed
RANDOM FITS to construct an identically gridded tracer-intersection proxy.
Retains all mock IDs and caps, disjoint-half noise, and response scatter.
This does not establish the exact selection mask or a mock galaxy covariance.
No observed/mock galaxy positions, odd data vector, or wake fit enter.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_full_random_rr import (
    ROOT, BINS, REFERENCE, read_random_bins,
)
from build_eboss_dr16_conditional_window import (
    OUT_ELLS, IN_ELLS, build_blocks, closure_summary,
    CONSTANT_CLOSURE_TOL, probe_even_to_odd,
)
from inspect_eboss_dr16_joint_randoms import (
    RANDOMS, grid_indices, N_PIX, RA_STEP, SIN_DEC_STEP,
)
from inspect_eboss_dr16_selection import download

IDS = (1, 125, 250, 375, 500, 625, 750, 875, 1000)
CAPS = ("NGC", "SGC")
TRACERS = ("LRG", "ELG")
THRESHOLDS = (1, 5, 10, 15, 25)
PREDECLARATION_COMMIT = "5fe7a0017da6e8a0b0b9de9cf2114d7c0ca9d0bd"
PINNED_OBS_RUN = "35995805473"
PINNED_OBS_SHA = "c842da035f4522926466cc28af5e8956e8286270"
MANIFEST = ROOT / "source_data/eboss_dr16_predeclared_9mock_window_pilot_2026-09-24.json"


def measure_grid_overlap(obs_l, obs_e, mock_l, mock_e):
    result = {}
    for threshold in THRESHOLDS:
        a = (obs_l >= threshold) & (obs_e >= threshold)
        b = (mock_l >= threshold) & (mock_e >= threshold)
        overlap = int(np.count_nonzero(a & b))
        union = int(np.count_nonzero(a | b))
        na = int(np.count_nonzero(a))
        nb = int(np.count_nonzero(b))
        result[str(threshold)] = {
            "observed_random_joint_supported_pixels": na,
            "mock_random_joint_supported_pixels": nb,
            "joint_intersection_pixels": overlap,
            "joint_union_pixels": union,
            "jaccard_of_thresholded_joint_random_grids":
                overlap / union if union else None,
            "fraction_of_observed_joint_pixels_in_mock":
                overlap / na if na else None,
            "fraction_of_mock_joint_pixels_in_observed":
                overlap / nb if nb else None,
        }
    return result


def self_test():
    a = np.array([0, 0, 20, 20, 21, 1], dtype="i8")
    b = np.array([0, 0, 20, 21, 20, 2], dtype="i8")
    c = np.array([0, 0, 20, 0, 21, 0], dtype="i8")
    d = np.array([0, 0, 0, 20, 21, 0], dtype="i8")
    x = measure_grid_overlap(a, b, c, d)
    assert x["15"]["observed_random_joint_supported_pixels"] == 3
    assert x["15"]["mock_random_joint_supported_pixels"] == 1
    assert x["15"]["joint_intersection_pixels"] == 1
    assert x["15"]["joint_union_pixels"] == 3
    assert abs(x["15"]["jaccard_of_thresholded_joint_random_grids"]
               - 1 / 3) < 1e-14
    print("EBOSS_NINEMOCK_AGGREGATION_SELF_TEST_OK", flush=True)


def descriptive(values):
    x = np.asarray(values, dtype="f8")
    if x.shape != (9,) or not np.isfinite(x).all():
        raise ValueError("All nine finite mock windows are required")
    return {
        "n_predeclared_realizations": len(x),
        "minimum": float(np.min(x)),
        "median": float(np.median(x)),
        "maximum": float(np.max(x)),
        "mean": float(np.mean(x)),
        "sample_std_descriptive_only": float(np.std(x, ddof=1)),
        "by_predeclared_id": {
            f"{mid:04d}": float(value) for mid, value in zip(IDS, x)
        },
    }


def get_observed_sky_grid(cache_dir, timeout):
    """Only four SHA-pinned public observed random FITS, never data FITS."""
    ref = json.loads(REFERENCE.read_text())
    if ref["status"] != "sample_weighted_normalization_audited":
        raise ValueError("Observed random input SHA provenance unavailable")
    known = {
        (x["cap"], x["tracer"]): x["sha256"]
        for x in ref["per_file_records"]
        if x["survey"] == "observed" and x["role"] == "random"
    }
    if len(known) != 4:
        raise ValueError("Observed reference lacks four random SHA256s")
    result, evidence = {}, []
    for cap in CAPS:
        for tracer in TRACERS:
            filename, expected_rows = RANDOMS[(tracer, cap)]
            path = Path(cache_dir) / filename
            try:
                path, sha, nbytes = download(
                    filename, Path(cache_dir), 1024 * 1024 * 1024, timeout)
                if sha != known[(cap, tracer)]:
                    raise ValueError("Published observed random SHA256 changed")
                cats, summary = read_random_bins(
                    path, tracer, cap, expected_rows)
                grids = []
                for iz, cat in enumerate(cats):
                    pixel = np.bincount(
                        grid_indices(cat[0], cat[1]), minlength=N_PIX
                    ).astype("<i4")
                    grids.append(pixel)
                    print("NINEMOCK_OBS_SKY_GRID_OK", cap, tracer, iz,
                          "rows", len(cat[0]), flush=True)
                candidate = np.sum(grids, axis=0, dtype="i8")
                result[(cap, tracer)] = (candidate, grids)
                evidence.append({
                    "cap": cap, "tracer": tracer, "filename": filename,
                    "sha256": sha, "downloaded_bytes": nbytes,
                    "retained_per_candidate_bin":
                        summary["retained_per_candidate_bin"],
                })
            finally:
                if path.exists():
                    path.unlink()
    return result, evidence


def run(args, output):
    manifest = json.loads(MANIFEST.read_text())
    if (tuple(manifest["predeclared_mock_realization_ids"]) != IDS
            or tuple(manifest["diagnostic_sky_grid"][
                "minimum_randoms_per_tracer_per_pixel"]) != THRESHOLDS
            or manifest["observed_odd_data_vector_read"] is not False):
        raise ValueError("Fixed nine-mock random-only declaration changed")
    obs_json = json.loads(Path(args.observed_audit).read_text())
    if (obs_json["status"]
            != "full_four_bin_fine_RR_geometric_closure_passed"
            or obs_json["revision_commit"] != PINNED_OBS_SHA
            or obs_json.get("uses_observed_odd_data_vector") is not False
            or obs_json.get("uses_observed_galaxy_positions") is not False):
        raise ValueError("Incorrect full-four-bin blinded observed RR source")
    expected_mock_sha = os.environ.get("SOURCE_MOCK_SHA")
    fine_npz = np.load(args.observed_rr, allow_pickle=False)
    fine = np.asarray(fine_npz["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(fine_npz["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(fine_npz["muedges"], dtype="f8")
    if fine.shape != (121,) or coarse.shape != (7,) or mu.shape != (241,):
        raise ValueError("Observed fine bin geometry changed")
    observed, obs_evidence = get_observed_sky_grid(
        args.observed_cache_dir, args.timeout)
    observed_cases = {}
    matrices = {}
    arrays = {"fine_sedges_mpc_over_h": fine,
              "output_sedges_mpc_over_h": coarse, "muedges": mu}
    for cap in CAPS:
        full_rr = np.asarray(
            fine_npz[f"fine_rr_{cap}_z3_full"], dtype="f8")
        full_norm = float(np.asarray(
            fine_npz[f"fine_norm_{cap}_z3_full"]).item())
        if full_rr.shape != (120, 240) or full_norm <= 0:
            raise ValueError("Observed high-z RR shape or norm incorrect")
        observed_cases[cap] = {
            "normalized_full_fine_rr": full_rr / full_norm,
            "blocks": build_blocks(full_rr, fine, coarse, mu),
        }
        if (max(x["max_constant_closure_abs_error"]
                for x in closure_summary(
                    observed_cases[cap]["blocks"]).values())
                >= CONSTANT_CLOSURE_TOL):
            raise ValueError("Observed conditional window constant closure failed")
        for (o, i), mat in observed_cases[cap]["blocks"].items():
            arrays[f"observed_M_{cap}_out{o}_in{i}"] = mat

    cases, errors = [], []
    for mid in IDS:
        for cap in CAPS:
            folder = Path(args.shards_dir) / f"eboss-ninemock-window-id{mid}"
            path = folder / f"{cap}.json"
            rr_path = folder / f"{cap}.npz"
            if not path.is_file() or not rr_path.is_file():
                errors.append(f"{mid}/{cap}: shard artifact missing")
                continue
            case = json.loads(path.read_text())
            if (case.get("status") != "mock_random_rr_shard_complete"
                    or case.get("predeclaration_commit")
                    != PREDECLARATION_COMMIT
                    or case.get("mock_id") != mid or case.get("cap") != cap
                    or case.get("observed_odd_data_vector_read") is not False
                    or case.get("mock_galaxy_data_read") is not False
                    or (expected_mock_sha and
                        case.get("revision_commit") != expected_mock_sha)):
                errors.append(f"{mid}/{cap}: input provenance invalid")
                continue
            arr = np.load(rr_path, allow_pickle=False)
            if not all(np.array_equal(
                arr[k], ref) for k, ref in (
                    ("fine_sedges_mpc_over_h", fine),
                    ("coarse_sedges_mpc_over_h", coarse),
                    ("muedges", mu))):
                errors.append(f"{mid}/{cap}: fine geometry differs")
                continue
            mock_rr = np.asarray(arr["fine_rr_full"], dtype="f8")
            norm = float(np.asarray(arr["fine_norm_full"]).item())
            obs_rr = observed_cases[cap]["normalized_full_fine_rr"]
            if mock_rr.shape != (120, 240) or norm <= 0:
                errors.append(f"{mid}/{cap}: mock RR invalid")
                continue
            mock_rr = mock_rr / norm
            rr_l1 = float(np.sum(np.abs(mock_rr - obs_rr)) /
                          max(np.sum(np.abs(mock_rr)), np.sum(np.abs(obs_rr))))
            blocks = {}
            matrix_scatter = {}
            probe_scatter = {}
            for o in OUT_ELLS:
                ref_diagonal = float(np.linalg.norm(
                    observed_cases[cap]["blocks"][(o, o)]))
                for i in IN_ELLS:
                    key = f"{o}<-{i}"
                    obs_mat = observed_cases[cap]["blocks"][(o, i)]
                    mock_mat = np.asarray(
                        arr[f"M_full_out{o}_in{i}"], dtype="f8")
                    if mock_mat.shape != obs_mat.shape or ref_diagonal <= 0:
                        raise ValueError("Invalid saved mock RR window matrix")
                    blocks[(o, i)] = mock_mat
                    arrays[f"mock_M_{cap}_id{mid:04d}_out{o}_in{i}"] = mock_mat
                    matrix_scatter[key] = {
                        "normalized_frobenius_difference_over_observed_odd_diagonal":
                            float(np.linalg.norm(mock_mat - obs_mat)
                                  / ref_diagonal),
                        "max_absolute_cell_difference": float(np.max(
                            np.abs(mock_mat - obs_mat))),
                    }
                    if i in (0, 2, 4):
                        observed_probe = probe_even_to_odd(
                            observed_cases[cap]["blocks"], fine, coarse, o, i)
                        mock_probe = probe_even_to_odd(
                            blocks if key in blocks else
                            {(o, i): mock_mat}, fine, coarse, o, i)
                        # Build is only required for the requested (o,i) key.
                        actual_reported = np.asarray(
                            case["predeclared_unitless_radial_ramp_probes"][
                                key]["full"])
                        if not np.allclose(mock_probe, actual_reported,
                                           atol=2e-12, rtol=0):
                            raise ValueError(
                                "Mock saved response and input RR probe disagree")
                        probe_scatter[key] = {
                            "max_abs_mock_minus_observed": float(
                                np.max(np.abs(mock_probe - observed_probe))),
                            "mock_minus_observed_by_coarse_s":
                                (mock_probe - observed_probe).tolist(),
                            "max_abs_half_A_minus_full_mock":
                                case["half_to_full_max_abs_ramp_probe_change"][
                                    key]["half_A"],
                            "max_abs_half_B_minus_full_mock":
                                case["half_to_full_max_abs_ramp_probe_change"][
                                    key]["half_B"],
                        }
            mask = {}
            for label, iz in (("candidate", None), *(
                    (f"z{j}", j) for j in range(4))):
                mock_l = np.asarray(
                    arr["sky_LRG_candidate" if iz is None
                        else f"sky_LRG_z{iz}"], dtype="i8")
                mock_e = np.asarray(
                    arr["sky_ELG_candidate" if iz is None
                        else f"sky_ELG_z{iz}"], dtype="i8")
                obs_l = (observed[(cap, "LRG")][0] if iz is None else
                         observed[(cap, "LRG")][1][iz])
                obs_e = (observed[(cap, "ELG")][0] if iz is None else
                         observed[(cap, "ELG")][1][iz])
                if any(a.shape != (N_PIX,) for a in (
                        mock_l, mock_e, obs_l, obs_e)):
                    raise ValueError("Mock or observed sky grid shape changed")
                mask[label] = measure_grid_overlap(
                    obs_l, obs_e, mock_l, mock_e)
            cases.append({
                "mock_id": mid, "cap": cap,
                "full_mock_vs_observed_normalized_fine_rr_relative_l1": rr_l1,
                "input_file_sha256": {
                    t: case["inputs"][t]["compressed_file_sha256"]
                    for t in TRACERS},
                "full_mock_RR_supported_cells":
                    case["full_and_half_samples"]["full"][
                        "positive_fine_smu_cells"],
                "half_A_vs_full_mock_normalized_fine_rr_l1":
                    case["full_and_half_samples"]["half_A"][
                        "normalized_fine_l1_vs_full"],
                "half_B_vs_full_mock_normalized_fine_rr_l1":
                    case["full_and_half_samples"]["half_B"][
                        "normalized_fine_l1_vs_full"],
                "conditional_window_block_differences": matrix_scatter,
                "predeclared_unitless_ramp_probe_differences": probe_scatter,
                "diagnostic_joint_random_sky_grid_overlap": mask,
            })
            print("NINEMOCK_AGGREGATE_CASE", mid, cap,
                  "RR_l1", rr_l1, "highz_jaccard_t15",
                  mask["z3"]["15"][
                      "jaccard_of_thresholded_joint_random_grids"],
                  flush=True)

    if errors or len(cases) != 18:
        raise ValueError("Incomplete fixed cohort: " + "; ".join(errors))
    summaries = {}
    for cap in CAPS:
        own = [c for c in cases if c["cap"] == cap]
        if [c["mock_id"] for c in own] != list(IDS):
            raise ValueError("Ensemble order differs from preregistration")
        summaries[cap] = {
            "mock_vs_observed_fine_rr_relative_l1": descriptive([
                c["full_mock_vs_observed_normalized_fine_rr_relative_l1"]
                for c in own]),
            "highz_joint_grid_jaccard_at_threshold_15": descriptive([
                c["diagnostic_joint_random_sky_grid_overlap"]["z3"]["15"][
                    "jaccard_of_thresholded_joint_random_grids"]
                for c in own]),
            "predeclared_operator_probe_max_abs_difference": {
                key: descriptive([
                    c["predeclared_unitless_ramp_probe_differences"][key][
                        "max_abs_mock_minus_observed"]
                    for c in own])
                for key in ("1<-0", "1<-2", "1<-4",
                            "3<-0", "3<-2", "3<-4")
            },
            "all_ten_full_operator_blocks_frobenius_scatter": {
                f"{o}<-{i}": descriptive([
                    c["conditional_window_block_differences"][f"{o}<-{i}"][
                        "normalized_frobenius_difference_over_observed_odd_diagonal"]
                    for c in own])
                for o in OUT_ELLS for i in IN_ELLS
            },
        }
    report = {
        "study": "Predeclared nine-mock eBOSS random-only window and joint-grid pilot",
        "status": "nine_mock_random_only_window_pilot_complete",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "predeclaration_commit": PREDECLARATION_COMMIT,
        "source_mock_shard_workflow_run": os.environ.get("SOURCE_MOCK_RUN"),
        "source_mock_shard_commit": expected_mock_sha,
        "source_observed_fine_rr_workflow_run": PINNED_OBS_RUN,
        "source_observed_fine_rr_commit": PINNED_OBS_SHA,
        "predeclared_mock_ids": list(IDS),
        "n_individual_mock_realizations": 9,
        "n_cap_by_realization_shards": len(cases),
        "candidate_z_bins": [list(v) for v in BINS],
        "highz_window_diagnostic": [0.9, 1.0],
        "caps": list(CAPS),
        "diagnostic_sky_grid": {
            "ra_step_deg": RA_STEP,
            "sin_dec_step": SIN_DEC_STEP,
            "thresholds_randoms_per_tracer_per_pixel": list(THRESHOLDS),
            "declared_reporting_threshold": 15,
            "exact_common_mask": False,
        },
        "observed_random_input_file_evidence": obs_evidence,
        "cases": cases, "descriptive_ensemble_summaries": summaries,
        "errors": errors,
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "wake_template_tuned_or_fitted": False,
        "joint_mock_galaxy_covariance_estimated": False,
        "physical_window_convolution_validated": False,
        "analysis_protocol_frozen_for_inference": False,
        "scope_note": (
            "Nine predeclared mock realizations describe window-input scatter "
            "and density-dependent thresholded sky-grid overlap. Disjoint "
            "halves are correlated with their full mock; they are not "
            "independent mock realizations. Neither a sky grid nor its "
            "Jaccard index is an exact joint selection mask, and nine mock "
            "random windows cannot calibrate observational significance."),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_suffix(".tmp.npz")
    np.savez_compressed(temp, **arrays)
    temp.replace(output.with_suffix(".npz"))
    temp_json = output.with_suffix(".tmp.json")
    temp_json.write_text(json.dumps(report, indent=2) + "\n",
                         encoding="utf-8")
    temp_json.replace(output)
    print("EBOSS_NINEMOCK_RANDOM_ONLY_ENSEMBLE_COMPLETE",
          len(cases), output, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shards-dir", help="gh run download output directory")
    ap.add_argument("--observed-audit", help="Pinned full-four-bin observed RR JSON")
    ap.add_argument("--observed-rr", help="Matching full-four-bin observed RR NPZ")
    ap.add_argument("--observed-cache-dir",
                    default="eboss_workspace/ninemock_observed_randoms")
    ap.add_argument("--out",
                    default="eboss_workspace/ninemock_window_ensemble.json")
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.shards_dir or not args.observed_audit or not args.observed_rr:
        ap.error("Provide all 18 pinned mock shards and observed RR inputs")
    output = Path(args.out)
    try:
        run(args, output)
    except (OSError, KeyError, ValueError, RuntimeError,
            MemoryError) as exc:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({
            "status": "nine_mock_random_only_ensemble_incomplete",
            "errors": [str(exc)],
            "observed_galaxy_data_read": False,
            "observed_odd_data_vector_read": False,
            "mock_galaxy_data_read": False,
        }, indent=2) + "\n", encoding="utf-8")
        print("NINEMOCK_AGGREGATION_ERROR", exc, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
