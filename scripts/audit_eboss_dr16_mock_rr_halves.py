#!/usr/bin/env python3
"""Blinded mock-internal finite-random control for eBOSS high-z fine RR.

For the previously fixed realistic mock IDs 0001, 0500, 1000 and both caps,
load the *same SHA-pinned random-only* catalogues used by the preceding
observed-vs-mock fine RR comparison. Keep the full mock RR from that validated
artifact. Independently permute each tracer's mock random rows using four
predeclared seeds; count the LRG_A x ELG_A and LRG_B x ELG_B disjoint halves,
normalizing each by its own product of weighted sums. Compare each half with
the corresponding full mock RR in fine (s,mu), and repeat the predeclared
conditional even->odd radial-ramp probe.

The earlier exact-count observed thinning control is a separate reference.
Neither control supplies a covariance, a formal test of survey-mask
adequacy, or a physical odd-sector result. No observed or mock galaxy
positions, wake templates, or real odd data vector are read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_full_random_rr import ROOT, coordinates, compile_rr
from audit_eboss_dr16_rr_density import histogram_distance
from audit_eboss_dr16_density_matched_fine_rr import (
    MOCK_IDS, CAPS, TRACERS, SEEDS, HIGH_Z, THREE_MOCK_REF, MOCK_REF,
    extract_mock_cases, verify_mock_input_sha, ramp_probes,
    summarize_relative,
)
from compare_eboss_dr16_mock_fine_rr import extract_mock_high_z
from eboss_dr16_fiducial import PRIMARY_GEOMETRY
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry

PINNED_DENSITY_REF_SHA = "3959cf5173299a9c7957c20d4b184c7573b75ce2"
PINNED_MOCK_COMMITS = {
    1: "009ff8873ae149e622a9ef22a33ada73d59fcc86",
    500: "b4e8434f389b399ab7142af1894b14c82907d5b8",
    1000: "b4e8434f389b399ab7142af1894b14c82907d5b8",
}


def disjoint_half_draw(cat, seed):
    n = len(cat[0])
    if n < 200 or any(len(x) != n for x in cat):
        raise ValueError("Invalid mock half-sample inputs")
    order = np.random.default_rng(seed).permutation(n)
    na = n // 2
    pieces = {"half_A": order[:na], "half_B": order[na:]}
    if len(pieces["half_A"]) + len(pieces["half_B"]) != n:
        raise AssertionError("Half-samples do not partition the full catalogue")
    out = {}
    for label, idx in pieces.items():
        part = tuple(np.asarray(x[idx], dtype="f8") for x in cat)
        w = part[3]
        if np.any(~np.isfinite(w)) or np.any(w <= 0):
            raise ValueError("Nonfinite or nonpositive half-sample weights")
        digest = hashlib.sha256()
        digest.update(np.ascontiguousarray(idx.astype("<i8")).tobytes())
        out[label] = (part, {
            "raw_rows": len(idx),
            "row_indices_sha256": digest.hexdigest(),
            "sum_weights": float(np.sum(w, dtype="f8")),
            "kish_effective_weighted_count": float(
                np.sum(w, dtype="f8") ** 2 / np.dot(w, w)),
        })
    return out


def self_test():
    rng = np.random.default_rng(713)
    cat = (
        rng.uniform(120, 126, 401),
        rng.uniform(0, 10, 401),
        rng.uniform(0.91, 0.99, 401),
        rng.uniform(0.8, 1.6, 401),
    )
    a = disjoint_half_draw(cat, 4071)
    b = disjoint_half_draw(cat, 4071)
    c = disjoint_half_draw(cat, 8159)
    assert len(a["half_A"][0][0]) == 200
    assert len(a["half_B"][0][0]) == 201
    for label in ("half_A", "half_B"):
        for x, y in zip(a[label][0], b[label][0]):
            np.testing.assert_array_equal(x, y)
    assert a["half_A"][1]["row_indices_sha256"] != c[
        "half_A"][1]["row_indices_sha256"]
    original = np.sort(cat[0])
    reunified = np.sort(np.r_[
        a["half_A"][0][0], a["half_B"][0][0]
    ])
    np.testing.assert_array_equal(original, reunified)
    assert abs(
        a["half_A"][1]["sum_weights"] + a["half_B"][1]["sum_weights"]
        - np.sum(cat[3], dtype="f8")
    ) < 1e-10
    print("EBOSS_MOCK_RANDOM_DISJOINT_HALF_SELF_TEST_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mock-0001-rr", required=False)
    ap.add_argument("--mock-0001-audit", required=False)
    ap.add_argument("--mock-stratified-dir", required=False)
    ap.add_argument("--density-control", required=False)
    ap.add_argument("--out", default="eboss_workspace/mock_random_half_control.json")
    ap.add_argument("--cache-dir", default="eboss_workspace/mock_random_half_files")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (args.threads < 1 or args.timeout <= 0
            or not all((args.mock_0001_rr, args.mock_0001_audit,
                        args.mock_stratified_dir, args.density_control))):
        ap.error("Supply the three pinned mock and matched-density audit inputs")

    density = json.loads(Path(args.density_control).read_text())
    if (density.get("status") != "density_matched_random_only_control_complete"
            or tuple(density["source_predeclared_mock_ids"]) != MOCK_IDS
            or tuple(density["predeclared_seeds"]) != SEEDS
            or density.get("revision_commit") != PINNED_DENSITY_REF_SHA
            or density.get("uses_observed_odd_sector_data_vector") is not False):
        raise ValueError("Input density control is not the pinned, blinded result")
    match_cases = {
        (int(v["mock_realization_id"]), v["cap"]): v
        for v in density["cases"]
    }
    if set(match_cases) != {(m, c) for m in MOCK_IDS for c in CAPS}:
        raise ValueError("Incomplete exact-count observed random control")
    release = json.loads(MOCK_REF.read_text())
    triple = json.loads(THREE_MOCK_REF.read_text())
    if (release.get("status") != "mock_random_sample_selection_compatible"
            or triple.get("status") !=
            "three_preselected_mock_fine_rr_comparisons_complete"):
        raise ValueError("Previously verified mock selection is unavailable")
    known_mocks = {
        (int(v["id"]), v["tracer"], v["cap"]): v
        for v in release["mock_random_catalogues"]
    }
    mock_paths = {
        1: (Path(args.mock_0001_audit), Path(args.mock_0001_rr)),
        500: (
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_0500.json",
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_0500.npz"),
        1000: (
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_1000.json",
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_1000.npz"),
    }
    mock_audits = {
        mid: json.loads(paths[0].read_text()) for mid, paths in mock_paths.items()
    }
    mock_npz = {
        mid: np.load(paths[1], allow_pickle=False)
        for mid, paths in mock_paths.items()
    }
    # Every preexisting successful mock artifact shares these exact grids.
    rr0 = mock_npz[1]
    fine = np.asarray(rr0["fine_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(rr0["muedges"], dtype="f8")
    if fine.shape != (121,) or mu.shape != (241,):
        raise ValueError("Pinned mock fine RR grid changed")
    coarse = np.array([20., 40., 60., 80., 100., 120., 140.])
    cases, errors = [], []
    cache = Path(args.cache_dir)
    for mock_id in MOCK_IDS:
        for ci, cap in enumerate(CAPS):
            try:
                mock_full, source_case = extract_mock_cases(
                    mock_audits[mock_id], mock_npz[mock_id], mock_id, cap,
                    PINNED_MOCK_COMMITS[mock_id],
                )
                expected_targets = verify_mock_input_sha(
                    source_case, known_mocks, mock_id, cap)
                reference = match_cases[(mock_id, cap)]
                if (reference["mock_target_raw_retained_random_rows"]
                        != expected_targets):
                    raise ValueError("Target row counts differ from density control")
                full_norm = float(np.asarray(
                    mock_npz[mock_id][
                        f"mock_fine_norm_{cap}_{mock_id:04d}"]).item())
                full_raw = np.asarray(mock_npz[mock_id][
                    f"mock_fine_rr_{cap}_{mock_id:04d}"], dtype="f8")
                if not np.allclose(full_raw / full_norm, mock_full,
                                   atol=0, rtol=5e-15):
                    raise ValueError("Preexisting full mock RR changed")
                cats, identities = {}, {}
                for tracer in TRACERS:
                    relative = mock_path("eBOSS_" + tracer, cap, "ran", mock_id)
                    filename = Path(relative).name
                    local = cache / filename
                    expected = known_mocks[(mock_id, tracer, cap)]
                    try:
                        digest, size = fetch_with_retry(
                            MOCK_BASE + relative, local,
                            512 * 1024 * 1024, args.timeout)
                        if digest != expected["compressed_file_sha256"]:
                            raise ValueError("Mock random SHA256 differs from source audit")
                        cat, meta = extract_mock_high_z(local, expected)
                        if meta["high_z_retained_rows"] != expected_targets[tracer]:
                            raise ValueError("High-z mock random count changed")
                        if abs(
                            meta["high_z_weight_sum"] /
                            source_case["inputs"][tracer]["high_z_weight_sum"] - 1
                        ) > 2e-11:
                            raise ValueError("High-z total input weight changed")
                        cats[tracer] = coordinates(cat, PRIMARY_GEOMETRY)
                        identities[tracer] = {
                            **meta, "full_compressed_sha256": digest,
                            "downloaded_compressed_bytes": size,
                        }
                        print("MOCK_HALF_INPUT_OK", mock_id, cap, tracer,
                              meta["high_z_retained_rows"], digest, flush=True)
                    finally:
                        if local.exists():
                            local.unlink()
                direct = float(
                    np.sum(cats["LRG"][3], dtype="f8")
                    * np.sum(cats["ELG"][3], dtype="f8"))
                if abs(direct / full_norm - 1) > 1e-9:
                    raise ValueError("Source full RR pair normalization differs from mock weights")
                full_probe = ramp_probes(mock_full, fine, coarse, mu)
                split_runs = []
                for seed in SEEDS:
                    splits = {}
                    for ti, tracer in enumerate(TRACERS):
                        selected_seed = seed + 10000 * ci + 100 * mock_id + ti
                        splits[tracer] = disjoint_half_draw(
                            cats[tracer], selected_seed)
                    per_half = {}
                    results = {}
                    for label in ("half_A", "half_B"):
                        half_rr, meta = compile_rr(
                            splits["LRG"][label][0],
                            splits["ELG"][label][0],
                            fine, mu, args.threads,
                        )
                        normalized = half_rr / meta["pair_weight_normalization"]
                        response = ramp_probes(normalized, fine, coarse, mu)
                        per_half[label] = normalized
                        probe_noise = {}
                        for key in full_probe:
                            change = response[key] - full_probe[key]
                            probe_noise[key] = {
                                "max_abs_half_minus_full_mock": float(
                                    np.max(np.abs(change))),
                                "half_minus_full_by_coarse_separation": change.tolist(),
                            }
                        results[label] = {
                            "catalogue_half_metadata": {
                                t: splits[t][label][1] for t in TRACERS
                            },
                            "half_vs_full_mock_fine_rr": summarize_relative(
                                mock_full, normalized),
                            "conditional_radial_ramp_probe_noise": probe_noise,
                        }
                        print("MOCK_HALF_RR", mock_id, cap, seed, label,
                              "relative_l1",
                              results[label]["half_vs_full_mock_fine_rr"][
                                  "fine_rr_l1_relative_to_reference"],
                              flush=True)
                    results["half_A_vs_half_B"] = summarize_relative(
                        per_half["half_B"], per_half["half_A"])
                    split_runs.append({"base_seed": seed, "results": results})
                scatter = [
                    record["results"][label]["half_vs_full_mock_fine_rr"][
                        "fine_rr_l1_relative_to_reference"]
                    for record in split_runs for label in ("half_A", "half_B")
                ]
                probe_spread = {}
                for key in full_probe:
                    changes = [
                        record["results"][label][
                            "conditional_radial_ramp_probe_noise"][key][
                                "max_abs_half_minus_full_mock"]
                        for record in split_runs for label in ("half_A", "half_B")
                    ]
                    obs_probe = reference["predeclared_operator_probe_comparisons"][key]
                    probe_spread[key] = {
                        "mock_half_vs_full_max_abs_probe_shift_range": [
                            float(min(changes)), float(max(changes))],
                        "mock_full_vs_observed_full_max_abs_probe_shift":
                            obs_probe["max_abs_mock_minus_full"],
                        "matched_observed_thinning_max_abs_probe_shift_range":
                            obs_probe["matched_observed_thinning_max_abs_shift_range"],
                    }
                case = {
                    "cap": cap, "mock_realization_id": mock_id,
                    "full_mock_random_file_evidence": identities,
                    "full_mock_vs_full_observed_rr_relative_l1": reference[
                        "mock_vs_full_observed_fine_RR"][
                            "fine_rr_l1_relative_to_reference"],
                    "count_matched_observed_thinning_l1_range": reference[
                        "count_matched_observed_thinning_l1_range"],
                    "mock_internal_half_vs_full_l1_range": [
                        float(min(scatter)), float(max(scatter))],
                    "predeclared_response_probe_comparison": probe_spread,
                    "four_predeclared_disjoint_half_runs": split_runs,
                }
                cases.append(case)
                print("MOCK_HALF_CAP_ID", cap, mock_id,
                      "mock_vs_obs_full_l1",
                      case["full_mock_vs_full_observed_rr_relative_l1"],
                      "obs_countmatched_minmax",
                      case["count_matched_observed_thinning_l1_range"],
                      "mock_half_minmax",
                      case["mock_internal_half_vs_full_l1_range"],
                      flush=True)
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError) as exc:
                errors.append(f"{cap}/mock_{mock_id:04d}: {exc}")
                print("MOCK_HALF_ERROR", errors[-1], flush=True)
    ok = len(cases) == len(CAPS) * len(MOCK_IDS) and not errors
    report = {
        "study": "Three preselected mock-internal random-half density controls",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "three_mock_internal_random_half_controls_complete" if ok
        else "partial",
        "input_exact_count_matched_observed_control_revision":
            PINNED_DENSITY_REF_SHA,
        "input_exact_count_matched_observed_control_run":
            os.environ.get("PINNED_DENSITY_RUN"),
        "source_pinned_mock_artifact_runs": {
            str(mid): os.environ.get("PINNED_MOCK_0001_RUN")
            if mid == 1 else os.environ.get("PINNED_MOCK_STRATIFIED_RUN")
            for mid in MOCK_IDS
        },
        "mock_realization_ids_chosen_before_real_odd_vector": list(MOCK_IDS),
        "predeclared_halving_seeds": list(SEEDS),
        "candidate_redshift_slice": list(HIGH_Z),
        "redshift_or_analysis_protocol_frozen": False,
        "fine_separation_step_mpc_over_h": 1.0,
        "mu_bins": len(mu) - 1,
        "distance_mapping": PRIMARY_GEOMETRY,
        "oriented_pair": "LRG->ELG", "los": "midpoint",
        "cases": cases, "errors": errors,
        "observed_galaxies_read": False,
        "mock_galaxies_read": False,
        "observed_odd_data_vector_read": False,
        "wake_template_tuned": False,
        "statistical_covariance_estimated": False,
        "mock_ensemble_window_validated": False,
        "scope_note": (
            "Disjoint halves estimate finite-mock-random sensitivity on "
            "the tested fine grid. Multiple permutations of one mock are "
            "correlated, and mock half-vs-full and observed count-matched "
            "vs-full differences are not identically distributed. "
            "Never combine them as an inferential variance or equate a "
            "single realization with an ensemble selection-function test."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("EBOSS_MOCK_HALF_RR_AUDIT", report["status"], output, flush=True)
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
