#!/usr/bin/env python3
"""Density-matched random-only control for eBOSS high-z LRG x ELG fine RR.

The three realistic mock IDs were selected at input-audit time (0001, 0500,
1000), before any real eBOSS odd-sector data vector. For each mock and cap,
draw exactly the mock's *raw retained row counts* of LRG and ELG from the
verified observed DR16 random catalogues, with four predeclared independent
seeds and a fixed no-replacement rule. Use each draw's own sum(w_LRG) *
sum(w_ELG) to normalize compiled RR. Compare it to (a) the complete observed
RR and (b) the existing matched mock RR on the identical fine grid. Test
the small even->odd conditional-response blocks with a predefined
within-bin linear radial probe.

A count-matched draw need not have identical weighted effective counts,
redshift distribution or mask. This is a *finite-random control*, not a
significance test, covariance estimate, mock-adequacy verdict or physical
wake measurement. No observed or mock galaxy FITS is opened.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_full_random_rr import (
    ROOT, read_random_bins, coordinates, compile_rr,
)
from audit_eboss_dr16_rr_density import histogram_distance
from build_eboss_dr16_conditional_window import (
    OUT_ELLS, build_blocks, probe_even_to_odd,
)
from eboss_dr16_fiducial import PRIMARY_GEOMETRY
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import download

MOCK_IDS = (1, 500, 1000)
CAPS = ("NGC", "SGC")
TRACERS = ("LRG", "ELG")
SEEDS = (4071, 8159, 12239, 16381)
HIGH_Z_INDEX = 3
HIGH_Z = (0.9, 1.0)
INPUT_REF = ROOT / "source_data/eboss_dr16_weighted_normalization_2026-09-24.json"
MOCK_REF = ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
FULL_FINE_REF = ROOT / "source_data/eboss_dr16_fine_rr_highz_2026-09-24.json"
THREE_MOCK_REF = ROOT / "source_data/eboss_dr16_mock_fine_rr_3id_2026-09-24.json"
PROBE_BLOCKS = ((o, i) for o in (1, 3) for i in (0, 2, 4))


def fixed_draw(cat: tuple[np.ndarray, ...], n: int, seed: int):
    if not (100 <= n <= len(cat[0])):
        raise ValueError(f"Cannot draw {n} without replacement from {len(cat[0])}")
    ix = np.sort(np.random.default_rng(seed).choice(
        len(cat[0]), size=n, replace=False))
    arrays = tuple(np.asarray(v[ix], dtype="f8") for v in cat)
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(ix.astype("<i8")).tobytes())
    effective = float(
        np.sum(arrays[3], dtype="f8") ** 2 / np.dot(arrays[3], arrays[3])
    )
    return arrays, {
        "seed": seed, "selected_row_indices_sha256": digest.hexdigest(),
        "raw_selected_rows": n, "sum_weights": float(np.sum(arrays[3], dtype="f8")),
        "kish_effective_weighted_count": effective,
    }


def extract_mock_cases(mock_audit: dict, mock_counts: dict, mock_id: int,
                       cap: str, expected_commit: str) -> tuple[np.ndarray, dict]:
    status = f"mock_{mock_id:04d}_fine_rr_comparison_complete"
    if (mock_audit.get("status") != status
            or int(mock_audit["realization_id"]) != mock_id
            or mock_audit.get("revision_commit") != expected_commit
            or mock_audit.get("observed_odd_data_vector_read") is not False
            or mock_audit.get("mock_galaxy_data_read") is not False):
        raise ValueError("Untrusted or unmatched mock fine RR audit")
    matching = [c for c in mock_audit["cases"] if c["cap"] == cap]
    if len(matching) != 1:
        raise ValueError("Mock audit must have one case for each cap")
    case = matching[0]
    raw = np.asarray(
        mock_counts[f"mock_fine_rr_{cap}_{mock_id:04d}"], dtype="f8")
    norm = float(np.asarray(
        mock_counts[f"mock_fine_norm_{cap}_{mock_id:04d}"]).item())
    if raw.shape != (120, 240) or norm <= 0:
        raise ValueError("Mock high-z fine RR geometry or norm differs")
    return raw / norm, case


def verify_mock_input_sha(
    case: dict, known: dict, mock_id: int, cap: str,
) -> dict[str, int]:
    targets = {}
    for tracer in TRACERS:
        meta = case["inputs"][tracer]
        pinned = known[(mock_id, tracer, cap)]
        if meta["compressed_file_sha256"] != pinned["compressed_file_sha256"]:
            raise ValueError("Mock random input SHA changed")
        targets[tracer] = int(meta["high_z_retained_rows"])
        if targets[tracer] <= 0:
            raise ValueError("Matched mock target random count is nonpositive")
    return targets


def ramp_probes(raw: np.ndarray, fine: np.ndarray,
                coarse: np.ndarray, muedges: np.ndarray) -> dict:
    blocks = build_blocks(raw, fine, coarse, muedges)
    return {
        f"{ellout}<-{ellin}": probe_even_to_odd(
            blocks, fine, coarse, ellout, ellin)
        for ellout in OUT_ELLS for ellin in (0, 2, 4)
    }


def summarize_relative(ref: np.ndarray, item: np.ndarray) -> dict:
    same = histogram_distance(item, ref)
    return {
        "fine_rr_l1_relative_to_reference":
            same["normalized_rr_l1_over_larger_sample_l1"],
        "fine_rr_total_pair_fraction_shift":
            same["total_pair_fraction_difference"],
    }


def self_test():
    rng = np.random.default_rng(5312)
    cat = (
        rng.uniform(130, 145, 300),
        rng.uniform(-1, 12, 300),
        rng.uniform(0.9, 1.0, 300),
        rng.uniform(0.6, 1.5, 300),
    )
    a, ma = fixed_draw(cat, 180, 4071)
    b, mb = fixed_draw(cat, 180, 4071)
    c, mc = fixed_draw(cat, 180, 8159)
    for x, y in zip(a, b):
        np.testing.assert_array_equal(x, y)
    assert ma == mb and ma["selected_row_indices_sha256"] != mc[
        "selected_row_indices_sha256"]
    assert len(a[0]) == 180 and ma["kish_effective_weighted_count"] <= 180
    ref = np.array([[1., 2.], [3., 4.]])
    assert summarize_relative(ref, ref)["fine_rr_l1_relative_to_reference"] == 0
    try:
        fixed_draw(cat, 301, 4071)
    except ValueError:
        pass
    else:
        raise AssertionError("Sampling with replacement or too many rows forbidden")
    print("EBOSS_DENSITY_MATCHED_RANDOM_ONLY_SELF_TEST_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--observed-rr", help="SHA-pinned successful observed fine RR NPZ")
    ap.add_argument("--observed-audit", help="Matching successful fine RR JSON")
    ap.add_argument("--mock-0001-rr", help="Successful ID0001 mock RR NPZ")
    ap.add_argument("--mock-0001-audit", help="Successful ID0001 mock JSON")
    ap.add_argument("--mock-stratified-dir", help="Successful IDs 0500/1000 files")
    ap.add_argument("--out", default="eboss_workspace/density_matched_rr_audit.json")
    ap.add_argument("--cache-dir", default="eboss_workspace/density_matched_rr_input")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (args.threads < 1 or args.timeout <= 0 or
        not all([args.observed_rr, args.observed_audit, args.mock_0001_rr,
                 args.mock_0001_audit, args.mock_stratified_dir])):
        ap.error("Supply all three successful source artifacts and valid runtime options")

    observed_audit = json.loads(Path(args.observed_audit).read_text())
    source_obs_sha = os.environ.get(
        "PINNED_OBSERVED_SHA", "5dbfca45c25c7af57bf05693022bfa98b499209d")
    if (observed_audit.get("status") != "fine_RR_high_z_geometric_closure_passed"
            or observed_audit.get("revision_commit") != source_obs_sha
            or observed_audit.get("observed_odd_data_vector_read") is not False):
        raise ValueError("Fine observed RR provenance or data blindness mismatch")
    input_ref = json.loads(INPUT_REF.read_text())
    mock_ref = json.loads(MOCK_REF.read_text())
    three_ref = json.loads(THREE_MOCK_REF.read_text())
    full_ref = json.loads(FULL_FINE_REF.read_text())
    if (input_ref["status"] != "sample_weighted_normalization_audited"
            or mock_ref["status"] != "mock_random_sample_selection_compatible"
            or three_ref["status"] != "three_preselected_mock_fine_rr_comparisons_complete"
            or full_ref["status"] != "fine_RR_high_z_geometric_closure_passed"
            or tuple(three_ref["preselected_realistic_ids"]) != MOCK_IDS):
        raise ValueError("Predeclared input provenance is incomplete")
    known_mocks = {
        (int(x["id"]), x["tracer"], x["cap"]): x
        for x in mock_ref["mock_random_catalogues"]
    }
    known_obs = {
        (x["tracer"], x["cap"]): x
        for x in input_ref["per_file_records"]
        if x["survey"] == "observed" and x["role"] == "random"
    }
    observed_npz = np.load(args.observed_rr, allow_pickle=False)
    fine = np.asarray(observed_npz["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(observed_npz["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(observed_npz["muedges"], dtype="f8")
    if (fine.shape != (121,) or coarse.shape != (7,) or mu.shape != (241,)
            or not np.allclose(np.diff(fine), 1.0, rtol=0, atol=1e-12)):
        raise ValueError("Observed RR grid is not the pinned high-z fine grid")
    mock_files = {
        1: (Path(args.mock_0001_audit), Path(args.mock_0001_rr)),
        500: (
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_0500.json",
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_0500.npz"),
        1000: (
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_1000.json",
            Path(args.mock_stratified_dir) / "matched_mock_fine_rr_1000.npz"),
    }
    expected_mock_commit = {
        1: os.environ.get(
            "PINNED_MOCK_0001_SHA",
            "009ff8873ae149e622a9ef22a33ada73d59fcc86"),
        500: os.environ.get(
            "PINNED_MOCK_STRATIFIED_SHA",
            "b4e8434f389b399ab7142af1894b14c82907d5b8"),
        1000: os.environ.get(
            "PINNED_MOCK_STRATIFIED_SHA",
            "b4e8434f389b399ab7142af1894b14c82907d5b8"),
    }
    mock_audits = {
        k: json.loads(paths[0].read_text()) for k, paths in mock_files.items()
    }
    mock_arrays = {
        k: np.load(paths[1], allow_pickle=False)
        for k, paths in mock_files.items()
    }

    results, errors = [], []
    root = Path(args.cache_dir)
    for ci, cap in enumerate(CAPS):
        observed_full = np.asarray(
            observed_npz[f"fine_rr_{cap}_full"], dtype="f8")
        norm_full = float(np.asarray(
            observed_npz[f"fine_norm_{cap}_full"]).item())
        if observed_full.shape != (120, 240) or norm_full <= 0:
            raise ValueError("Observed reference RR shape or normalization changed")
        observed_full /= norm_full
        full_probe = ramp_probes(observed_full, fine, coarse, mu)
        candidates = {}
        for tracer in TRACERS:
            filename, nrows = RANDOMS[(tracer, cap)]
            path = root / filename
            try:
                local, checksum, size = download(
                    filename, root, 1024 * 1024 * 1024, args.timeout)
                if checksum != known_obs[(tracer, cap)]["sha256"]:
                    raise ValueError("Observed random FITS SHA256 differs from pinned input")
                selected, info = read_random_bins(
                    local, tracer, cap, nrows)
                high_z = selected[HIGH_Z_INDEX]
                exp = next(
                    x for x in full_ref["cap_results"]
                    if x["cap"] == cap
                )
                exp_count = int(
                    full_ref["full_observed_random_catalogue_counts_from_prior_full_rr"][
                        cap][tracer])
                if len(high_z[0]) != exp_count:
                    raise ValueError("Observed high-z retained random count changed")
                candidates[tracer] = high_z
                print("DENSITY_MATCH_INPUT_OK", cap, tracer,
                      len(high_z[0]), checksum, flush=True)
            except (OSError, ValueError, KeyError, RuntimeError,
                    MemoryError) as exc:
                errors.append(f"{cap}/{tracer}: {exc}")
                print("DENSITY_MATCH_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()
        if set(candidates) != set(TRACERS):
            continue
        # Transform positions once for each tracer, then thin exact same
        # array indices. This does not apply different distances to mocks.
        preconverted = {
            tracer: coordinates(candidates[tracer], PRIMARY_GEOMETRY)
            for tracer in TRACERS
        }
        for mock_id in MOCK_IDS:
            try:
                mock, mock_case = extract_mock_cases(
                    mock_audits[mock_id], mock_arrays[mock_id], mock_id,
                    cap, expected_mock_commit[mock_id],
                )
                targets = verify_mock_input_sha(
                    mock_case, known_mocks, mock_id, cap)
                if any(targets[t] > len(preconverted[t][0]) for t in TRACERS):
                    raise ValueError("Observed random supply below matched mock target")
                mock_probe = ramp_probes(mock, fine, coarse, mu)
                mock_vs_full = summarize_relative(observed_full, mock)
                previous = next(x for x in
                    three_ref["per_realization_cap"]
                    if int(x["id"]) == mock_id and x["cap"] == cap)
                if abs(
                    mock_vs_full["fine_rr_l1_relative_to_reference"]
                    - previous["normalized_rr_relative_l1"]
                ) > 1e-12:
                    raise ValueError("Mock RR numerical identity differs from pinned summary")
                controls = []
                for si, seed in enumerate(SEEDS):
                    samples, metadata = {}, {}
                    for ti, tracer in enumerate(TRACERS):
                        this_seed = seed + 10000 * ci + 100 * mock_id + ti
                        samples[tracer], metadata[tracer] = fixed_draw(
                            preconverted[tracer], targets[tracer], this_seed)
                    counts, stats = compile_rr(
                        samples["LRG"], samples["ELG"], fine, mu,
                        args.threads)
                    thinned = counts / stats["pair_weight_normalization"]
                    noise = summarize_relative(observed_full, thinned)
                    mock_to_matched = summarize_relative(thinned, mock)
                    thin_probe = ramp_probes(thinned, fine, coarse, mu)
                    probes = {}
                    for key in full_probe:
                        noise_probe = thin_probe[key] - full_probe[key]
                        mismatch_probe = mock_probe[key] - full_probe[key]
                        probes[key] = {
                            "matched_thin_minus_full_observed": noise_probe.tolist(),
                            "max_abs_matched_thin_minus_full": float(
                                np.max(np.abs(noise_probe))),
                            "max_abs_mock_minus_full": float(
                                np.max(np.abs(mismatch_probe))),
                        }
                    controls.append({
                        "base_seed": seed,
                        "catalogue_draws": metadata,
                        "normalized_RR_matched_observed_minus_full": noise,
                        "normalized_RR_mock_minus_matched_observed": mock_to_matched,
                        "predeclared_even_to_odd_ramp_probe": probes,
                        "normalized_RR_fraction_in_candidate_s_bins":
                            stats["normalized_rr_weight_fraction_in_test_bins"],
                    })
                    print(
                        "DENSITY_MATCHED_RR", cap, mock_id, seed,
                        "obs_thin_vs_full_l1",
                        noise["fine_rr_l1_relative_to_reference"],
                        "mock_vs_thin_l1",
                        mock_to_matched["fine_rr_l1_relative_to_reference"],
                        flush=True)
                baseline = np.array(
                    [x["normalized_RR_matched_observed_minus_full"][
                        "fine_rr_l1_relative_to_reference"] for x in controls],
                    dtype="f8",
                )
                cases_probes = {}
                for key in full_probe:
                    obs_noise = [
                        x["predeclared_even_to_odd_ramp_probe"][key][
                            "max_abs_matched_thin_minus_full"] for x in controls]
                    cases_probes[key] = {
                        "max_abs_mock_minus_full": float(np.max(np.abs(
                            mock_probe[key] - full_probe[key]))),
                        "matched_observed_thinning_max_abs_shift_range": [
                            float(min(obs_noise)), float(max(obs_noise))
                        ],
                    }
                results.append({
                    "cap": cap, "mock_realization_id": mock_id,
                    "mock_target_raw_retained_random_rows": targets,
                    "observed_raw_retained_random_rows": {
                        t: int(len(candidates[t][0])) for t in TRACERS
                    },
                    "mock_vs_full_observed_fine_RR": mock_vs_full,
                    "count_matched_observed_thinning_l1_range": [
                        float(np.min(baseline)), float(np.max(baseline))
                    ],
                    "predeclared_operator_probe_comparisons": cases_probes,
                    "four_predeclared_observed_draws": controls,
                    "count_matching_does_not_match_weighted_effective_density": True,
                })
                print(
                    "DENSITY_MATCHED_CAP_ID", cap, mock_id,
                    "mock_full_l1",
                    mock_vs_full["fine_rr_l1_relative_to_reference"],
                    "thin_control_min", float(np.min(baseline)),
                    "thin_control_max", float(np.max(baseline)), flush=True,
                )
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError) as exc:
                errors.append(f"{cap}/mock_{mock_id:04d}: {exc}")
                print("DENSITY_MATCHED_ERROR", errors[-1], flush=True)
    complete = len(results) == len(CAPS) * len(MOCK_IDS) and not errors
    out = {
        "study": "Predeclared density-matched observed-random control for eBOSS high-z fine RR",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "density_matched_random_only_control_complete" if complete else "partial",
        "source_run_observed_fine_RR": os.environ.get("PINNED_OBSERVED_RUN"),
        "source_run_mock_0001": os.environ.get("PINNED_MOCK_0001_RUN"),
        "source_run_mock_stratified": os.environ.get("PINNED_MOCK_STRATIFIED_RUN"),
        "source_predeclared_mock_ids": list(MOCK_IDS),
        "predeclared_seeds": list(SEEDS),
        "sampling": "Uniform without replacement, observed random row counts matched separately to each mock/cap/tracer; direct pair-weight normalization per draw.",
        "caps": list(CAPS),
        "candidate_high_z_interval": list(HIGH_Z),
        "candidate_interval_frozen": False,
        "geometry_primary": PRIMARY_GEOMETRY,
        "fine_s_edges_mpc_over_h": fine.tolist(),
        "coarse_s_edges_mpc_over_h": coarse.tolist(),
        "mu_edges": mu.tolist(),
        "orientation": "LRG->ELG", "los": "midpoint",
        "cases": results, "errors": errors,
        "uses_observed_galaxy_positions": False,
        "uses_mock_galaxy_positions": False,
        "uses_observed_odd_sector_data_vector": False,
        "physical_template_fit_or_tuning": False,
        "statistical_covariance_estimated": False,
        "full_mock_ensemble_validated": False,
        "conclusion_limit": (
            "Four count-matched observed thinning draws per fixed mock/cap are "
            "correlated because they reuse the same observed catalogue. "
            "The control is descriptive: matching raw row counts does not "
            "ensure equal sum of squared weights, effective redshift "
            "selection, mask, or identical sky support. No p-value, "
            "mock-adequacy decision or final physical window is reported."),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print("EBOSS_DENSITY_MATCHED_RR_AUDIT", out["status"], target, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
