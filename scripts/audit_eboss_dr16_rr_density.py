#!/usr/bin/env python3
"""Fixed-plan random-density study for blinded eBOSS weighted RR pair counting.

This samples 600, 1200 and 2400 nested random objects per tracer, with two
predeclared seeds, from the already audited observed and realistic-mock-0001
random FITS files. It measures RR subsampling noise only. No galaxy data,
odd data vector, mock covariance or production survey window is computed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_rr_pair_closure import (
    REAL_REF, MOCK_REF, ROOT, REAL_BASE, MOCK_BASE, RANDOMS,
    DIAGNOSTIC_Z, DIAGNOSTIC_SEP_EDGES, legacy_distance,
    rr_histogram, odd_rr_moments, select_fixed_random,
)
from inspect_eboss_dr16_mock_headers import mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry
from inspect_eboss_dr16_selection import download

LEVELS = (600, 1200, 2400)
SEEDS = (4071, 8159)
MOCK_ID = 1
THETA_MIN_DEG = 0.05
MU_BINS = 24


def normalized_rr(hist: np.ndarray, norm: dict) -> np.ndarray:
    value = norm["pair_normalization"]
    if value <= 0 or not np.isfinite(value):
        raise ValueError("Invalid random-pair normalization")
    return hist / value


def histogram_distance(sample: np.ndarray, reference: np.ndarray) -> dict:
    if sample.shape != reference.shape:
        raise ValueError("RR histogram geometry differs")
    denominator = float(np.sum(np.abs(reference)))
    if denominator <= 0:
        raise ValueError("The larger RR sample has no accepted pair weight")
    per_bin = []
    for a, b in zip(sample, reference):
        den = float(np.sum(np.abs(b)))
        per_bin.append(float(np.sum(np.abs(a - b)) / den) if den else None)
    return {
        "normalized_rr_l1_over_larger_sample_l1": float(
            np.sum(np.abs(sample - reference)) / denominator),
        "per_separation_bin_l1_relative": per_bin,
        "total_pair_fraction_difference": float(
            (np.sum(sample) - np.sum(reference)) / np.sum(reference)),
    }


def fixed_nested_sample(
    cat: tuple[np.ndarray, ...], seed: int
) -> tuple[dict[int, tuple[np.ndarray, ...]], str]:
    if not all(len(a) == LEVELS[-1] for a in cat):
        raise ValueError("Expected the predetermined maximum sample size")
    priority = np.random.default_rng(seed).permutation(LEVELS[-1])
    digest = hashlib.sha256(np.ascontiguousarray(
        priority.astype("<i8")).tobytes()).hexdigest()
    return {
        n: tuple(np.asarray(field[priority[:n]], dtype="f8")
                 for field in cat)
        for n in LEVELS
    }, digest


def self_test() -> None:
    rng = np.random.default_rng(113)
    cat = (
        rng.uniform(5, 15, LEVELS[-1]), rng.uniform(-5, 5, LEVELS[-1]),
        rng.uniform(0.8, 0.9, LEVELS[-1]), rng.uniform(0.5, 2, LEVELS[-1]),
    )
    nested, fingerprint = fixed_nested_sample(cat, 7821)
    assert len(fingerprint) == 64
    assert np.all(np.isin(nested[600][0], nested[1200][0]))
    assert np.all(np.isin(nested[1200][0], nested[2400][0]))
    h = np.array([[1., 2.], [0., 3.]])
    n = normalized_rr(h, {"pair_normalization": 2.0})
    assert histogram_distance(n, n)["normalized_rr_l1_over_larger_sample_l1"] == 0
    print("Fixed nested-sampling and normalized RR distance self-tests passed")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="eboss_workspace/rr_density_diagnostic.json")
    p.add_argument("--cache-dir", default="eboss_workspace/rr_density_catalogues")
    p.add_argument("--timeout", type=float, default=120)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout <= 0:
        p.error("timeout must be positive")
    real_ref = json.loads(REAL_REF.read_text(encoding="utf-8"))
    mock_ref = json.loads(MOCK_REF.read_text(encoding="utf-8"))
    if (real_ref["status"] != "four_random_catalogues_checked"
            or mock_ref["status"] != "mock_random_sample_selection_compatible"):
        raise RuntimeError("The pre-odd random input evidence is incomplete")
    real_sha = {(v["tracer"], v["cap"]): v["sha256"]
                for v in real_ref["randoms"]}
    mock_sha = {
        (int(v["id"]), v["tracer"], v["cap"]): v["compressed_file_sha256"]
        for v in mock_ref["mock_random_catalogues"]
    }
    sedges = np.array(DIAGNOSTIC_SEP_EDGES, dtype="f8")
    muedges = np.linspace(-1.0 - 1e-7, 1.0 + 1e-7, MU_BINS + 1)
    failures, reports = [], []
    root = Path(args.cache_dir)

    for cap_number, cap in enumerate(("NGC", "SGC")):
        for survey_number, survey in enumerate(("observed", "realistic_mock")):
            draws: dict[int, dict[str, dict]] = {seed: {} for seed in SEEDS}
            inputs = {}
            for tracer_number, tracer in enumerate(("LRG", "ELG")):
                real = survey == "observed"
                if real:
                    filename, expected_rows = RANDOMS[(tracer, cap)]
                    path = root / survey / filename
                else:
                    relative = mock_path(
                        "eBOSS_" + tracer, cap, "ran", MOCK_ID)
                    path = root / survey / Path(relative).name
                try:
                    if real:
                        path, sha, size = download(
                            filename, path.parent, 1024 * 1024 * 1024, args.timeout)
                        if sha != real_sha[(tracer, cap)]:
                            raise ValueError("Observed random SHA256 does not match retained input")
                    else:
                        sha, size = fetch_with_retry(
                            MOCK_BASE + relative, path,
                            512 * 1024 * 1024, args.timeout)
                        if sha != mock_sha[(MOCK_ID, tracer, cap)]:
                            raise ValueError("Realistic-mock random SHA256 does not match retained input")
                    inputs[tracer] = {
                        "file_sha256": sha, "compressed_or_fits_bytes": size,
                        "filename": path.name,
                    }
                    for seed in SEEDS:
                        base_seed = (
                            seed + 100 * cap_number
                            + 10 * survey_number + tracer_number
                        )
                        max_cat, audit = select_fixed_random(
                            path, LEVELS[-1], base_seed,
                            DIAGNOSTIC_Z[0], DIAGNOSTIC_Z[1])
                        nested, priority_sha = fixed_nested_sample(
                            max_cat, base_seed + 10000)
                        draws[seed][tracer] = {
                            "catalogues": nested, "maximum_sample_audit": audit,
                            "nested_priority_sha256": priority_sha,
                        }
                    print("RR_DENSITY_INPUT_OK", survey, cap, tracer, sha, flush=True)
                except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                    failures.append(f"{survey}/{cap}/{tracer}: {exc}")
                    print("RR_DENSITY_INPUT_ERROR", failures[-1], flush=True)
                finally:
                    if path.exists():
                        path.unlink()
            if not all(set(draws[s]) == {"LRG", "ELG"} for s in SEEDS):
                continue
            per_seed, max_results = [], {}
            for seed in SEEDS:
                hists, levels = {}, []
                for n in LEVELS:
                    first = draws[seed]["LRG"]["catalogues"][n]
                    second = draws[seed]["ELG"]["catalogues"][n]
                    raw, norm = rr_histogram(
                        first, second, sedges, muedges,
                        THETA_MIN_DEG, legacy_distance, block=128)
                    hists[n] = normalized_rr(raw, norm)
                    odd = odd_rr_moments(raw, muedges)
                    levels.append({
                        "random_objects_per_tracer": n,
                        "accepted_pair_count": norm["accepted_pairs"],
                        "weighted_rr_fraction_within_test_bins": float(np.sum(hists[n])),
                        "raw_rr_geometric_odd_moments": odd["per_separation_bin"],
                    })
                reference = hists[LEVELS[-1]]
                for entry in levels:
                    n = entry["random_objects_per_tracer"]
                    entry["distance_to_largest_sample"] = histogram_distance(
                        hists[n], reference)
                max_results[seed] = reference
                per_seed.append({
                    "seed": seed,
                    "lrg_max_sample_sha256": draws[seed]["LRG"][
                        "maximum_sample_audit"]["sample_arrays_sha256"],
                    "elg_max_sample_sha256": draws[seed]["ELG"][
                        "maximum_sample_audit"]["sample_arrays_sha256"],
                    "lrg_nested_priority_sha256": draws[seed]["LRG"][
                        "nested_priority_sha256"],
                    "elg_nested_priority_sha256": draws[seed]["ELG"][
                        "nested_priority_sha256"],
                    "levels": levels,
                })
                print(f"RR_DENSITY {survey} {cap} seed={seed} "
                      f"n={LEVELS[-1]} accepted={levels[-1]['accepted_pair_count']} "
                      f"small_to_large_l1="
                      f"{levels[0]['distance_to_largest_sample']['normalized_rr_l1_over_larger_sample_l1']:.6f}",
                      flush=True)
            between = histogram_distance(
                max_results[SEEDS[0]], max_results[SEEDS[1]])
            reports.append({
                "survey": survey, "cap": cap, "inputs": inputs,
                "per_seed": per_seed,
                "largest_sample_between_seed_comparison": between,
            })
            print(f"RR_DENSITY_SEEDS {survey} {cap} n={LEVELS[-1]} "
                  f"l1={between['normalized_rr_l1_over_larger_sample_l1']:.6f}",
                  flush=True)

    complete = len(reports) == 4 and not failures
    result = {
        "study": "eBOSS random-only fixed-plan RR subsampling density diagnostic",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "fixed_density_diagnostic_complete" if complete else "partial",
        "source_real_random_reference": str(REAL_REF.relative_to(ROOT)),
        "source_mock_random_reference": str(MOCK_REF.relative_to(ROOT)),
        "mock_realization_id": MOCK_ID,
        "sample_levels_per_tracer": list(LEVELS),
        "predeclared_sampling_seeds": list(SEEDS),
        "nested_sampling": (
            "The maximum sample is a uniform deterministic draw without "
            "replacement. One fixed permutation orders nested samples of "
            "600, 1200 and 2400 from that same maximum sample."
        ),
        "diagnostic_redshift_interval": list(DIAGNOSTIC_Z),
        "diagnostic_interval_frozen_for_inference": False,
        "separation_edges_Mpc_over_h": sedges.tolist(),
        "mu_edges": muedges.tolist(), "theta_min_deg": THETA_MIN_DEG,
        "los": "midpoint", "orientation": "LRG_to_ELG",
        "provisional_weights": (
            "WEIGHT_FKP*WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ"),
        "distance_mapping": (
            "legacy FlatLambdaCDM H0=67.4 Om0=0.315, diagnostic only"),
        "final_eboss_fiducial_and_weight_conventions_validated": False,
        "reports": reports, "errors": failures,
        "observed_odd_sector_vector_read": False,
        "galaxy_pairs_computed": False,
        "complete_RR_window_computed": False,
        "mock_covariance_computed": False,
        "scope_note": (
            "These fixed nested samples quantify random-pair sampling noise, "
            "not the physical eBOSS window or an inferential significance. "
            "No threshold is optimized and the two seeds are not a validated "
            "mock covariance ensemble."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("RR_DENSITY_AUDIT", result["status"], out, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
