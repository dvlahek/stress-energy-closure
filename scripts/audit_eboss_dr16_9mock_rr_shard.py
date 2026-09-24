#!/usr/bin/env python3
"""Predeclared nine-ID eBOSS realistic mock-random high-z RR shard.

One fixed realization and one cap per job. Reads ONLY matching public
LRG/ELG realistic mock RANDOM FITS (never galaxy/data FITS), validates
full compressed SHA256, selection, numerical-zero weights, and four-bin
diagnostic sky support before counting high-z 0.9-1.0 RR. Records full
and fixed disjoint-half 1-Mpc/h RR, independently counted coarse RR,
and the conditional finite-bin response. No odd data vector, wake
template, galaxy pairs, physical inference, or mock covariance enters.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_full_random_rr import (
    ROOT, BINS, SEPARATION_EDGES, MU_BINS, compile_rr, coordinates,
)
from audit_eboss_dr16_fine_rr_highz import fine_edges, rebin_to_coarse
from build_eboss_dr16_conditional_window import (
    build_blocks, closure_summary, CONSTANT_CLOSURE_TOL, IN_ELLS,
    OUT_ELLS, probe_even_to_odd,
)
from compare_eboss_dr16_mock_fine_rr import extract_mock_high_z, self_test as earlier_self_test
from eboss_dr16_fiducial import PRIMARY_GEOMETRY
from inspect_eboss_dr16_mock_headers import BASE, mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry
from inspect_eboss_dr16_mock_randoms import (
    inspect_mock_random, sample_is_healthy, joint_support,
)
from inspect_eboss_dr16_joint_randoms import N_PIX

MANIFEST = ROOT / "source_data/eboss_dr16_predeclared_9mock_window_pilot_2026-09-24.json"
PREDECLARATION_COMMIT = "5fe7a0017da6e8a0b0b9de9cf2114d7c0ca9d0bd"
PRIOR = ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
EXPECTED_IDS = (1, 125, 250, 375, 500, 625, 750, 875, 1000)
SEED_ROOT = 4071
CAPS = ("NGC", "SGC")


def test_fixed_scheme():
    earlier_self_test()
    rng = np.random.default_rng(5904)
    x = rng.normal(size=101)
    for mid in EXPECTED_IDS:
        for ci in range(2):
            order = np.random.default_rng(
                SEED_ROOT + 10000 * ci + 100 * mid
            ).permutation(len(x))
            assert len(np.unique(order)) == len(x)
            np.testing.assert_array_equal(np.sort(np.r_[order[:50], order[50:]]),
                                          np.arange(len(x)))
    print("EBOSS_NINEMOCK_SHARD_SELF_TEST_OK", flush=True)


def full_and_half(cat, seed):
    n = len(cat[0])
    if n < 200 or any(len(item) != n for item in cat):
        raise ValueError("Insufficient or inconsistent high-z mock randoms")
    order = np.random.default_rng(seed).permutation(n)
    middle = n // 2
    return {
        "full": cat,
        "half_A": tuple(v[order[:middle]] for v in cat),
        "half_B": tuple(v[order[middle:]] for v in cat),
    }


def run(args, out):
    manifest = json.loads(MANIFEST.read_text())
    if (tuple(manifest["predeclared_mock_realization_ids"]) != EXPECTED_IDS
            or manifest["highz_window_pilot"] != [0.9, 1.0]
            or manifest["observed_odd_data_vector_read"] is not False
            or manifest["weight_product"]
            != "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP"):
        raise ValueError("Nine-mock predeclaration changed or is not blinded")
    if args.mock_id not in EXPECTED_IDS or args.cap not in CAPS:
        raise ValueError("Mock shard outside fixed predeclared cohort")
    prior = json.loads(PRIOR.read_text())
    if prior["status"] != "mock_random_sample_selection_compatible":
        raise ValueError("Earlier three-mock random input gate unavailable")
    prior_files = {
        (int(x["id"]), x["tracer"], x["cap"]): x
        for x in prior["mock_random_catalogues"]
    }
    fine = fine_edges()
    coarse = np.asarray(SEPARATION_EDGES, dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    arrays = {"fine_sedges_mpc_over_h": fine,
              "coarse_sedges_mpc_over_h": coarse, "muedges": mu}
    records, cats, grids = {}, {}, {}
    ci = CAPS.index(args.cap)
    cache = Path(args.cache_dir)
    for ti, tracer in enumerate(("LRG", "ELG")):
        rel = mock_path("eBOSS_" + tracer, args.cap, "ran", args.mock_id)
        path = cache / Path(rel).name
        try:
            sha, nbytes = fetch_with_retry(
                BASE + rel, path, 512 * 1024 * 1024, args.timeout)
            report, sky, per_z = inspect_mock_random(
                path, tracer, args.cap, args.mock_id, sha, nbytes, 150000)
            if not sample_is_healthy([report]):
                raise ValueError("Mock random weight/selection health failed")
            old = prior_files.get((args.mock_id, tracer, args.cap))
            if old and (sha != old["compressed_file_sha256"]
                        or report["header_rows"] != old["header_rows"]
                        or report["candidate_random_rows"] != old["candidate_random_rows"]
                        or report["candidate_excluded_by_weight_convention"]
                        != old["candidate_zero_weight_exclusions"]):
                raise ValueError("Existing SHA-pinned mock random input changed")
            reference = {
                "header_rows": report["header_rows"],
                "candidate_random_rows": report["candidate_random_rows"],
                "candidate_zero_weight_exclusions":
                    report["candidate_excluded_by_weight_convention"],
            }
            cat, meta = extract_mock_high_z(path, reference)
            if meta["high_z_retained_rows"] != report["redshift_bin_retained_rows"][3]:
                raise ValueError("Two independent high-z mock selections disagree")
            cats[tracer] = cat
            records[tracer] = {**report, "high_z": meta,
                               "existing_three_mock_sha_pinned": old is not None}
            grids[tracer] = (sky, per_z)
            arrays["sky_" + tracer + "_candidate"] = sky.astype("<i4")
            for iz in range(4):
                arrays[f"sky_{tracer}_z{iz}"] = per_z[iz].astype("<i4")
            print("NINEMOCK_INPUT_OK", args.mock_id, args.cap, tracer,
                  sha, meta["high_z_retained_rows"], flush=True)
        finally:
            if path.exists():
                path.unlink()

    if set(cats) != {"LRG", "ELG"}:
        raise ValueError("Incomplete matched LRG/ELG mock random input")
    joint = {
        "candidate": joint_support(grids["LRG"][0], grids["ELG"][0]),
        "by_z": [joint_support(grids["LRG"][1][iz], grids["ELG"][1][iz])
                 for iz in range(4)]
    }
    converted = {}
    half_seeds = {}
    for ti, tracer in enumerate(("LRG", "ELG")):
        seed = SEED_ROOT + 10000 * ci + 100 * args.mock_id + ti
        half_seeds[tracer] = seed
        converted[tracer] = full_and_half(
            coordinates(cats[tracer], PRIMARY_GEOMETRY), seed)
    samples, matrices, probes = {}, {}, {}
    for sample in ("full", "half_A", "half_B"):
        rr, meta = compile_rr(converted["LRG"][sample],
                              converted["ELG"][sample], fine, mu, args.threads)
        if sample == "full" and np.count_nonzero(rr) != 28800:
            raise ValueError("Missing full fine-RR support")
        norm = float(meta["pair_weight_normalization"])
        arrays[f"fine_rr_{sample}"] = rr
        arrays[f"fine_norm_{sample}"] = np.array([norm], dtype="f8")
        samples[sample] = {
            "lrg_rows": len(converted["LRG"][sample][0]),
            "elg_rows": len(converted["ELG"][sample][0]),
            "positive_fine_smu_cells": int(np.count_nonzero(rr)),
            "pair_weight_normalization": norm,
            "normalized_pair_fraction": meta[
                "normalized_rr_weight_fraction_in_test_bins"],
        }
        blocks = build_blocks(rr, fine, coarse, mu)
        closure = max(val["max_constant_closure_abs_error"]
                      for val in closure_summary(blocks).values())
        if closure >= CONSTANT_CLOSURE_TOL:
            raise ValueError(f"Failed constant theory closure in {sample}")
        matrices[sample] = blocks
        samples[sample]["max_constant_theory_closure_error"] = closure
        for o in OUT_ELLS:
            for i in IN_ELLS:
                if sample == "full":
                    arrays[f"M_full_out{o}_in{i}"] = blocks[(o, i)]
                if i in (0, 2, 4):
                    probes.setdefault(f"{o}<-{i}", {})[sample] = (
                        probe_even_to_odd(blocks, fine, coarse, o, i).tolist())
        print("NINEMOCK_FINE_RR_OK", args.mock_id, args.cap, sample,
              "positive", np.count_nonzero(rr), "closure", closure, flush=True)

    direct, direct_info = compile_rr(
        converted["LRG"]["full"], converted["ELG"]["full"],
        coarse, mu, args.threads)
    rebin = rebin_to_coarse(arrays["fine_rr_full"], fine, coarse)
    norm = samples["full"]["pair_weight_normalization"]
    denom = max(float(np.sum(np.abs(rebin / norm))),
                float(np.sum(np.abs(
                    direct / direct_info["pair_weight_normalization"]))))
    closure_l1 = float(np.sum(np.abs(
        rebin / norm - direct / direct_info["pair_weight_normalization"])) / denom)
    if closure_l1 > 1e-8:
        raise ValueError("Independent full fine->coarse RR closure failed")
    mirror = {}
    for o in OUT_ELLS:
        for i in IN_ELLS:
            from build_eboss_dr16_conditional_window import conditional_window_block
            val = conditional_window_block(
                arrays["fine_rr_full"][:, ::-1], fine, coarse, mu, o, i)
            parity = -1 if (o + i) % 2 else 1
            mirror[f"{o}<-{i}"] = float(np.max(np.abs(
                val - parity * matrices["full"][(o, i)])))
    if max(mirror.values()) > 1e-10:
        raise ValueError("Conditional mirror-parity closure failed")
    half_noise = {}
    for key, val in probes.items():
        base = np.asarray(val["full"])
        half_noise[key] = {
            sample: float(np.max(np.abs(np.asarray(val[sample]) - base)))
            for sample in ("half_A", "half_B")
        }
    for half in ("half_A", "half_B"):
        rr = arrays[f"fine_rr_{half}"] / samples[half][
            "pair_weight_normalization"]
        full = arrays["fine_rr_full"] / norm
        samples[half]["normalized_fine_l1_vs_full"] = float(
            np.sum(np.abs(rr - full)) /
            max(np.sum(np.abs(rr)), np.sum(np.abs(full))))
    report = {
        "study": "Nine-predeclared-ID per-cap mock-random window and sky-grid shard",
        "status": "mock_random_rr_shard_complete",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "predeclaration_commit": PREDECLARATION_COMMIT,
        "mock_id": args.mock_id, "cap": args.cap,
        "candidate_bins_frozen_for_inference": False,
        "high_z_index": 3, "high_z": [0.9, 1.0],
        "inputs": records, "sky_grid_joint_support": joint,
        "sky_grid_thresholds_predeclared": [1, 5, 10, 15, 25],
        "full_and_half_samples": samples,
        "half_sample_seeds": half_seeds,
        "independent_coarse_rr_rebin_relative_l1": closure_l1,
        "max_mirror_parity_error": max(mirror.values()),
        "predeclared_unitless_radial_ramp_probes": probes,
        "half_to_full_max_abs_ramp_probe_change": half_noise,
        "errors": [],
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "wake_template_fitted_or_tuned": False,
        "statistical_mock_covariance_estimated": False,
        "exact_joint_mask_certified": False,
        "scope_note": (
            "Matched random-only per-realization input and finite-bin window "
            "diagnostics. Gridded joint support is not an exact angular mask. "
            "Nine realizations do not supply an observational covariance."),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    temp_npz = out.with_suffix(".tmp.npz")
    np.savez_compressed(temp_npz, **arrays)
    temp_npz.replace(out.with_suffix(".npz"))
    temp_json = out.with_suffix(".tmp.json")
    temp_json.write_text(json.dumps(report, indent=2) + "\n",
                         encoding="utf-8")
    temp_json.replace(out)
    print("NINEMOCK_SHARD_COMPLETE", args.mock_id, args.cap,
          "fine_rebin_l1", closure_l1,
          "joint_pixels_at_15",
          joint["by_z"][3]["15"]["joint_supported_pixels"],
          flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mock-id", type=int, choices=EXPECTED_IDS)
    ap.add_argument("--cap", choices=CAPS)
    ap.add_argument("--out", default="eboss_workspace/ninemock_shard.json")
    ap.add_argument("--cache-dir", default="eboss_workspace/ninemock_fits")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        test_fixed_scheme()
        return 0
    if args.mock_id is None or args.cap is None or args.threads < 1:
        ap.error("Supply one preregistered mock ID, one cap and positive threads")
    try:
        run(args, Path(args.out))
    except (OSError, ValueError, KeyError, RuntimeError,
            MemoryError, AssertionError) as exc:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({
            "status": "mock_random_rr_shard_failed",
            "mock_id": args.mock_id, "cap": args.cap,
            "predeclaration_commit": PREDECLARATION_COMMIT,
            "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "mock_galaxy_data_read": False,
        }, indent=2) + "\n", encoding="utf-8")
        print("NINEMOCK_SHARD_FAILED", args.mock_id, args.cap, exc, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
