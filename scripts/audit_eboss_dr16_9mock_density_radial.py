#!/usr/bin/env python3
"""Blinded exact-count nine-mock angular and coarse-radial selection control.

Uses only SHA-pinned observed random catalogues and the 18 SHA-audited
realistic EZmock random-only shard artefacts of the previously registered
nine fixed mock IDs. No galaxy positions, pair odd vector, template fit,
mock-galaxy covariance, or mask-dependent signal optimization enters.

Counts are compared on the identical declared (RA,sin(dec)) sky grid,
with two *parametric* exact-count multinomial replicates per catalogue
at the minimum full count of the observed+fixed mock cohort per
cap/tracer/z-bin. Within-catalogue replicate differences quantify the
noise of these resampled gridded inputs; they are not independent mock
realizations and do not validate a continuous exact common mask.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from aggregate_eboss_dr16_9mock_window import (
    IDS, CAPS, TRACERS, THRESHOLDS, N_PIX, RA_STEP, SIN_DEC_STEP,
    PREDECLARATION_COMMIT, get_observed_sky_grid,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_9mock_density_radial_protocol_2026-09-24.json"
SOURCE_MOCK_RUN = "36015246541"
SOURCE_MOCK_SHA = "9150992ac06c03cf1b16c5fdd4e3bb3856f89ff0"
OBS_SEED_ROOT = 7711
MOCK_SEED_ROOT = 12911
Z_LABELS = ("candidate", "z0", "z1", "z2", "z3")
N_REPLICATES = 2


def total_variation(a, b):
    a = np.asarray(a, dtype="f8")
    b = np.asarray(b, dtype="f8")
    if a.shape != b.shape or np.any(a < 0) or np.any(b < 0):
        raise ValueError("Invalid nonnegative matched sky histograms")
    an, bn = float(np.sum(a)), float(np.sum(b))
    if an <= 0 or bn <= 0:
        raise ValueError("Cannot normalize an empty sky grid")
    return float(0.5 * np.sum(np.abs(a / an - b / bn), dtype="f8"))


def jaccard_joint(a_l, a_e, b_l, b_e, threshold):
    a = (a_l >= threshold) & (a_e >= threshold)
    b = (b_l >= threshold) & (b_e >= threshold)
    union = int(np.count_nonzero(a | b))
    intersection = int(np.count_nonzero(a & b))
    return {
        "intersection": intersection, "union": union,
        "jaccard": intersection / union if union else None,
        "observed_joint_pixels": int(np.count_nonzero(a)),
        "mock_joint_pixels": int(np.count_nonzero(b)),
    }


def exact_count_replicates(grid, target, seeds):
    arr = np.asarray(grid, dtype="i8")
    if (arr.shape != (N_PIX,) or np.any(arr < 0)
            or target <= 0 or target > int(np.sum(arr))
            or len(seeds) != N_REPLICATES):
        raise ValueError("Invalid fixed exact-count resampling inputs")
    prob = arr.astype("f8") / float(np.sum(arr))
    draws = [
        np.random.default_rng(seed).multinomial(target, prob)
        for seed in seeds
    ]
    if any(int(np.sum(x)) != target for x in draws):
        raise AssertionError("Multinomial exact-count sampling failed")
    return draws


def self_test():
    a = np.zeros(N_PIX, dtype="i8")
    b = np.zeros(N_PIX, dtype="i8")
    a[:4] = [80, 40, 20, 60]
    b[:4] = [70, 50, 25, 55]
    assert total_variation(a, a) == 0
    assert 0 < total_variation(a, b) < 1
    assert jaccard_joint(a, a, a, a, 5)["jaccard"] == 1
    assert jaccard_joint(a, a, b, b, 85)["jaccard"] is None
    draws = exact_count_replicates(a, 100, [OBS_SEED_ROOT, OBS_SEED_ROOT + 1])
    assert all(int(np.sum(x)) == 100 for x in draws)
    assert not np.array_equal(draws[0], draws[1])
    for i in range(2):
        check = exact_count_replicates(a, 100, [OBS_SEED_ROOT,
                                                OBS_SEED_ROOT + 1])[i]
        np.testing.assert_array_equal(draws[i], check)
    print("EBOSS_NINEMOCK_MATCHED_COUNT_SKY_AND_RADIAL_SELF_TEST_OK",
          flush=True)


def label_key(tracer, label):
    return "sky_" + tracer + "_" + label


def summarize(nine):
    x = np.asarray(nine, dtype="f8")
    if x.shape != (len(IDS),) or not np.isfinite(x).all():
        raise ValueError("Descriptive summary requires all nine mocks")
    return {
        "min": float(np.min(x)),
        "median": float(np.median(x)),
        "max": float(np.max(x)),
        "mean": float(np.mean(x)),
        "by_predeclared_id": {
            str(mid): float(y) for mid, y in zip(IDS, x)
        },
    }


def run(args, output):
    protocol = json.loads(PROTOCOL.read_text())
    if (tuple(protocol["mock_ids"]) != IDS
            or tuple(protocol["caps"]) != CAPS
            or tuple(protocol["tracers"]) != TRACERS
            or tuple(protocol["fixed_angular_grid"][
                "thresholds_randoms_per_tracer_per_pixel"]) != THRESHOLDS
            or protocol["control_seeds"]["observed_seed_base"] != OBS_SEED_ROOT
            or protocol["control_seeds"]["mock_seed_base"] != MOCK_SEED_ROOT
            or protocol["observed_odd_data_vector_read"] is not False):
        raise ValueError("Registered blinded fixed-control protocol differs")
    ensemble = json.loads(Path(args.ensemble_audit).read_text())
    if (ensemble.get("status") != "nine_mock_random_only_window_pilot_complete"
            or ensemble.get("predeclaration_commit")
            != PREDECLARATION_COMMIT
            or ensemble.get("source_mock_shard_workflow_run")
            != SOURCE_MOCK_RUN
            or ensemble.get("source_mock_shard_commit")
            != SOURCE_MOCK_SHA
            or tuple(ensemble.get("predeclared_mock_ids", [])) != IDS
            or ensemble.get("observed_odd_data_vector_read") is not False
            or ensemble.get("observed_galaxy_data_read") is not False):
        raise ValueError("Nine-mock source ensemble provenance invalid")

    observed, observed_files = get_observed_sky_grid(
        args.observed_cache_dir, args.timeout)
    observed_ref = {
        (row["cap"], row["tracer"]): row
        for row in ensemble["observed_random_input_file_evidence"]
    }
    if len(observed_ref) != 4:
        raise ValueError("Expected all four observed-random file identities")
    for x in observed_files:
        ref = observed_ref[(x["cap"], x["tracer"])]
        if (x["sha256"] != ref["sha256"]
                or x["retained_per_candidate_bin"]
                != ref["retained_per_candidate_bin"]):
            raise ValueError("Observed random file/binned rows differ from ensemble")

    mock = {}
    for mid in IDS:
        for cap in CAPS:
            folder = Path(args.shards_dir) / f"eboss-ninemock-window-id{mid}"
            json_path, npz_path = folder / f"{cap}.json", folder / f"{cap}.npz"
            if not json_path.is_file() or not npz_path.is_file():
                raise ValueError(f"Missing fixed mock shard {mid}/{cap}")
            record = json.loads(json_path.read_text())
            if (record.get("status") != "mock_random_rr_shard_complete"
                    or record.get("mock_id") != mid
                    or record.get("cap") != cap
                    or record.get("revision_commit") != SOURCE_MOCK_SHA
                    or record.get("predeclaration_commit")
                    != PREDECLARATION_COMMIT
                    or record.get("observed_odd_data_vector_read") is not False
                    or record.get("mock_galaxy_data_read") is not False):
                raise ValueError(f"Invalid fixed mock input {mid}/{cap}")
            with np.load(npz_path, allow_pickle=False) as z:
                sky = {
                    (t, label): np.asarray(z[label_key(t, label)], dtype="i8")
                    for t in TRACERS for label in Z_LABELS
                }
            for tracer in TRACERS:
                counts = record["inputs"][tracer]["redshift_bin_retained_rows"]
                if len(counts) != 5:
                    raise ValueError("Expected five source redshift counts")
                if not np.array_equal(
                    sky[(tracer, "candidate")],
                    np.sum([sky[(tracer, f"z{i}")]
                            for i in range(4)], axis=0, dtype="i8")
                ):
                    raise ValueError(f"Mock sky candidate map not sum of z-bins: {mid}/{cap}/{tracer}")
                for iz in range(4):
                    if int(np.sum(sky[(tracer, f"z{iz}")])) != counts[iz]:
                        raise ValueError("Mock saved sky map disagrees with source FITS z count")
            mock[(mid, cap)] = (record, sky)

    z_selection, selection_summaries = [], {}
    for cap in CAPS:
        selection_summaries[cap] = {}
        for tracer in TRACERS:
            obs_z = np.asarray(
                observed_ref[(cap, tracer)]["retained_per_candidate_bin"],
                dtype="f8")
            if obs_z.shape != (4,) or np.any(obs_z <= 0):
                raise ValueError("Observed coarse z-bin count shape invalid")
            per_mock = []
            for mid in IDS:
                rec, _ = mock[(mid, cap)]
                m_z = np.asarray(
                    rec["inputs"][tracer]["redshift_bin_retained_rows"][:4],
                    dtype="f8")
                if np.any(m_z <= 0) or m_z.shape != obs_z.shape:
                    raise ValueError("Mock coarse radial support missing")
                per_mock.append({
                    "mock_id": mid,
                    "mock_bin_counts": m_z.astype("i8").tolist(),
                    "mock_bin_fractions": (m_z / sum(m_z)).tolist(),
                    "observed_bin_fractions": (obs_z / sum(obs_z)).tolist(),
                    "coarse_unweighted_z_fraction_total_variation":
                        total_variation(obs_z, m_z),
                })
            z_selection.extend(
                {"cap": cap, "tracer": tracer, **row}
                for row in per_mock)
            selection_summaries[cap][tracer] = summarize([
                v["coarse_unweighted_z_fraction_total_variation"]
                for v in per_mock])

    sky_cases, sky_summaries = [], {}
    for ci, cap in enumerate(CAPS):
        sky_summaries[cap] = {}
        for zi, label in enumerate(Z_LABELS):
            obs_grids = {
                t: (observed[(cap, t)][0] if label == "candidate" else
                    observed[(cap, t)][1][zi - 1])
                for t in TRACERS
            }
            for t in TRACERS:
                check = obs_grids[t]
                if check.shape != (N_PIX,):
                    raise ValueError("Observed sky pixel shape changed")
                if label != "candidate" and (
                    int(np.sum(check))
                    != observed_ref[(cap, t)]["retained_per_candidate_bin"][zi-1]
                ):
                    raise ValueError("Observed sky grid disagrees with source z count")
            # Matched count is fixed before examining any response metric:
            # smallest count of observed+all nine fixed mocks, per tracer.
            target = {
                t: min(
                    [int(np.sum(obs_grids[t]))] +
                    [int(np.sum(mock[(mid, cap)][1][(t, label)]))
                     for mid in IDS]
                ) for t in TRACERS
            }
            if any(x < 100 for x in target.values()):
                raise ValueError("Insufficient common input sample for matched-count diagnostic")
            obs_draw = {
                ti: exact_count_replicates(
                    obs_grids[ti], target[ti],
                    [OBS_SEED_ROOT + ci*100000 + zi*1000
                     + it*100 + k for k in range(N_REPLICATES)]
                )
                for it, ti in enumerate(TRACERS)
            }
            obs_within = {
                str(t): jaccard_joint(
                    obs_draw["LRG"][0], obs_draw["ELG"][0],
                    obs_draw["LRG"][1], obs_draw["ELG"][1], t)
                for t in THRESHOLDS
            }
            per_mock = []
            for mid in IDS:
                _, grid = mock[(mid, cap)]
                m = {t: grid[(t, label)] for t in TRACERS}
                draws = {
                    tr: exact_count_replicates(
                        m[tr], target[tr],
                        [MOCK_SEED_ROOT + ci*1000000 + zi*100000
                         + ti*10000 + mid*10 + k
                         for k in range(N_REPLICATES)]
                    )
                    for ti, tr in enumerate(TRACERS)
                }
                matched_tv = {
                    tr: {
                        "cross_matched_replicate_TV": [
                            total_variation(obs_draw[tr][k], draws[tr][k])
                            for k in range(N_REPLICATES)
                        ],
                        "observed_within_replicate_TV":
                            total_variation(obs_draw[tr][0], obs_draw[tr][1]),
                        "mock_within_replicate_TV":
                            total_variation(draws[tr][0], draws[tr][1]),
                    }
                    for tr in TRACERS
                }
                thresholds = {}
                for t in THRESHOLDS:
                    thresholds[str(t)] = {
                        "unthinned_joint": jaccard_joint(
                            obs_grids["LRG"], obs_grids["ELG"],
                            m["LRG"], m["ELG"], t),
                        "matched_cross": [
                            jaccard_joint(
                                obs_draw["LRG"][k], obs_draw["ELG"][k],
                                draws["LRG"][k], draws["ELG"][k], t)
                            for k in range(N_REPLICATES)
                        ],
                        "matched_observed_within": obs_within[str(t)],
                        "matched_mock_within": jaccard_joint(
                            draws["LRG"][0], draws["ELG"][0],
                            draws["LRG"][1], draws["ELG"][1], t),
                    }
                case = {
                    "mock_id": mid, "cap": cap, "z_label": label,
                    "matched_exact_target_count_per_tracer": target,
                    "full_normalized_tracer_sky_total_variation": {
                        tr: total_variation(obs_grids[tr], m[tr])
                        for tr in TRACERS},
                    "matched_count_tracer_sky_total_variation": matched_tv,
                    "fixed_joint_threshold_comparisons": thresholds,
                }
                sky_cases.append(case)
                per_mock.append(case)
                print("NINEMOCK_DENSITY_CASE", mid, cap, label,
                      "targets", target,
                      "matched_jaccard_t1",
                      thresholds["1"]["matched_cross"][0]["jaccard"],
                      "matched_jaccard_t15",
                      thresholds["15"]["matched_cross"][0]["jaccard"],
                      flush=True)
            sky_summaries[cap][label] = {
                "matched_exact_target_count_per_tracer": target,
                "tracer_full_sky_TV_by_mock": {
                    tr: summarize([
                        c["full_normalized_tracer_sky_total_variation"][tr]
                        for c in per_mock])
                    for tr in TRACERS},
                "threshold_1_matched_cross_jaccard_rep0": summarize([
                    c["fixed_joint_threshold_comparisons"]["1"][
                        "matched_cross"][0]["jaccard"]
                    for c in per_mock]),
                "threshold_15_matched_cross_jaccard_rep0": summarize([
                    c["fixed_joint_threshold_comparisons"]["15"][
                        "matched_cross"][0]["jaccard"]
                    for c in per_mock]),
                "threshold_15_mock_internal_within_jaccard": summarize([
                    c["fixed_joint_threshold_comparisons"]["15"][
                        "matched_mock_within"]["jaccard"]
                    for c in per_mock]),
                "threshold_15_observed_internal_within_jaccard":
                    obs_within["15"]["jaccard"],
            }
    if len(sky_cases) != 90 or len(z_selection) != 36:
        raise ValueError("Fixed nine-ID, two-cap, five-grid cohort incomplete")
    report = {
        "study": "Blinded eBOSS nine-mock matched-count grid and coarse radial controls",
        "status": "nine_mock_matched_count_angular_coarse_radial_control_complete",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "parent_predeclaration_commit": PREDECLARATION_COMMIT,
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "source_mock_shard_workflow_run": SOURCE_MOCK_RUN,
        "source_mock_shard_commit": SOURCE_MOCK_SHA,
        "source_window_ensemble_run": ensemble.get("revision_commit"),
        "predeclared_ids": list(IDS), "n_mock_ids": len(IDS),
        "n_cap_by_id_sky_slices": len(sky_cases),
        "n_cap_tracer_id_coarse_radial_cases": len(z_selection),
        "fixed_angular_grid": {"ra_step_deg": RA_STEP,
                               "sin_dec_step": SIN_DEC_STEP,
                               "thresholds": list(THRESHOLDS)},
        "parametric_resample_seed_bases": {
            "observed": OBS_SEED_ROOT, "mock": MOCK_SEED_ROOT
        },
        "matched_draws_per_catalogue": N_REPLICATES,
        "observed_sha_verified_random_files": observed_files,
        "coarse_unweighted_radial_selection": z_selection,
        "sky_matched_count_cases": sky_cases,
        "radial_selection_summaries": selection_summaries,
        "sky_selection_summaries": sky_summaries,
        "errors": [],
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "wake_template_fitted_or_tuned": False,
        "exact_joint_selection_mask_established": False,
        "fine_weighted_radial_nz_validated": False,
        "joint_mock_galaxy_covariance_estimated": False,
        "observational_significance_computed": False,
        "scope_note": (
            "Exact-count multinomial draws from audited unweighted coarse "
            "pixel maps compare observed/mock joint support at equal count. "
            "Replicates share empirical parent maps and are not independent "
            "mock realizations. Counts in four coarse z slices cannot "
            "certify fine or weighted radial selection. Grid support is "
            "not an exact continuous veto/selection mask. No p-value, "
            "odd measurement, or physical wake inference is supported."),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    tmp.replace(output)
    print("EBOSS_NINEMOCK_MATCHED_COUNT_RADIAL_CONTROL_COMPLETE",
          len(sky_cases), len(z_selection), output, flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ensemble-audit")
    ap.add_argument("--shards-dir")
    ap.add_argument("--observed-cache-dir",
                    default="eboss_workspace/ninemock_density_obs_randoms")
    ap.add_argument("--out",
                    default="eboss_workspace/ninemock_density_radial.json")
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.shards_dir or not args.ensemble_audit:
        ap.error("Require pinned nine-mock ensemble JSON and all 18 input shards")
    output = Path(args.out)
    try:
        run(args, output)
    except (OSError, KeyError, ValueError, RuntimeError,
            MemoryError, AssertionError) as exc:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({
            "status": "nine_mock_matched_count_angular_coarse_radial_incomplete",
            "errors": [str(exc)],
            "observed_galaxy_data_read": False,
            "observed_odd_data_vector_read": False,
            "mock_galaxy_data_read": False,
        }, indent=2) + "\n", encoding="utf-8")
        print("EBOSS_NINEMOCK_DENSITY_RADIAL_ERROR", exc, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
