#!/usr/bin/env python3
"""Corrfunc weighted eBOSS RR pilot on observed *randoms only*.

Predeclared pilot: NGC and SGC, diagnostic 0.9 <= z < 1.0, levels
6000/18000/54000 per tracer, two seeds 4071/8159. The chosen high-z
slice follows the independent random-support audit, not any odd statistic.
All samples are deterministic and nested within each seed. The same
catalogues are measured under the published multitracer distance
convention and (at maximum level only) the DR16 ELG fiducial alternative.
No observed galaxy catalogue, mock odd vector, wake template, D-R or D-D
pair is read. This is not a final eBOSS RR window or likelihood.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from eboss_dr16_fiducial import (
    DISTANCE_CONVENTIONS, PRIMARY_GEOMETRY, SECONDARY_GEOMETRY,
    comoving_mpc_over_h,
)
from audit_eboss_dr16_rr_pair_closure import (
    REAL_REF, ROOT, rr_histogram, select_fixed_random,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_selection import download
from audit_eboss_dr16_rr_density import histogram_distance

LEVELS = (6000, 18000, 54000)
SEEDS = (4071, 8159)
DIAGNOSTIC_Z = (0.9, 1.0)
SEPARATION_EDGES = (20., 40., 60., 80., 100., 120., 140.)
MU_BINS = 24
THETA_MIN_DEG = 0.05


def subset(cat, indices):
    return tuple(np.asarray(a[indices], dtype="f8") for a in cat)


def corrfunc_rr(cat1, cat2, sedges, muedges, distance_profile, nthreads):
    from pycorr import TwoPointCounter
    ra1, dec1, z1, w1 = cat1
    ra2, dec2, z2, w2 = cat2
    first = [
        ra1, dec1, comoving_mpc_over_h(z1, distance_profile)
    ]
    second = [
        ra2, dec2, comoving_mpc_over_h(z2, distance_profile)
    ]
    result = TwoPointCounter(
        mode="smu",
        edges=(sedges, muedges),
        positions1=first, positions2=second,
        weights1=w1, weights2=w2, position_type="rdd",
        los="midpoint",
        selection_attrs={"theta": (THETA_MIN_DEG, 180.0)},
        engine="corrfunc", nthreads=nthreads,
        compute_sepsavg=False,
    )
    counts = np.asarray(result.wcounts, dtype="f8")
    expected_shape = (len(sedges) - 1, len(muedges) - 1)
    if counts.shape != expected_shape:
        raise ValueError(f"Corrfunc RR shape {counts.shape} != {expected_shape}")
    if not np.isfinite(counts).all() or np.any(counts < 0):
        raise ValueError("Corrfunc RR contains nonfinite/negative weighted pairs")
    expected_norm = float(
        np.sum(w1, dtype="f8") * np.sum(w2, dtype="f8"))
    backend_norm = float(result.wnorm)
    norm_rel = abs(backend_norm - expected_norm) / expected_norm
    if norm_rel > 1e-9:
        raise ValueError(
            "Corrfunc pair normalization differs from direct weight sums: "
            f"relative={norm_rel:g}")
    if backend_norm <= 0:
        raise ValueError("RR pair normalization is nonpositive")
    normalized = counts / backend_norm
    return normalized, {
        "backend": "pycorr/orrfunc" .replace("orrfunc", "Corrfunc"),
        "weighted_rr_total": float(np.sum(counts, dtype="f8")),
        "weighted_pair_normalization": backend_norm,
        "pair_norm_relative_difference_vs_direct_sums": norm_rel,
        "normalized_rr_sum_in_test_bins": float(np.sum(normalized, dtype="f8")),
        "raw_rr_mu_asymmetry_l1_over_total": [
            float(np.sum(np.abs(row - row[::-1])) /
                  np.sum(row, dtype="f8")) if np.any(row) else None
            for row in counts
        ],
    }


def synthetic_self_test() -> None:
    rng = np.random.default_rng(31639)
    def make(n):
        return (
            rng.uniform(145, 149, n),
            rng.uniform(10, 14, n),
            rng.uniform(0.91, 0.99, n),
            rng.uniform(0.7, 1.5, n),
        )
    a, b = make(125), make(133)
    sedges = np.asarray(SEPARATION_EDGES, dtype="f8")
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    independent, norm = rr_histogram(
        a, b, sedges, muedges, THETA_MIN_DEG,
        lambda z: comoving_mpc_over_h(z, PRIMARY_GEOMETRY),
        block=23,
    )
    compiled, metadata = corrfunc_rr(
        a, b, sedges, muedges, PRIMARY_GEOMETRY, 2)
    reference = independent / norm["pair_normalization"]
    np.testing.assert_allclose(
        compiled, reference, atol=2e-9, rtol=2e-8,
        err_msg="Corrfunc RR differs from independent NumPy pair counter")
    assert metadata["normalized_rr_sum_in_test_bins"] > 0
    assert abs(np.sum(compiled) - np.sum(reference)) < 2e-9
    print("EBOSS_CORRFUNC_INDEPENDENT_RR_SELF_TEST_OK", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/eboss_corrfunc_rr_pilot.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/eboss_corrfunc_rr_input")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        synthetic_self_test()
        return 0
    if args.timeout <= 0 or args.threads < 1:
        parser.error("Invalid timeout or thread count")
    reference = json.loads(REAL_REF.read_text())
    if reference["status"] != "four_random_catalogues_checked":
        raise ValueError("Observed random provenance must be validated first")
    expected = {(r["tracer"], r["cap"]): r for r in reference["randoms"]}
    rng_names = ("NGC", "SGC")
    sedges = np.asarray(SEPARATION_EDGES, dtype="f8")
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    records, failures = [], []
    root = Path(args.cache_dir)
    for ci, cap in enumerate(rng_names):
        seed_samples = {seed: {} for seed in SEEDS}
        input_identity = {}
        for ti, tracer in enumerate(("LRG", "ELG")):
            filename, expected_nrows = RANDOMS[(tracer, cap)]
            path = root / filename
            try:
                path, checksum, nbytes = download(
                    filename, root, 1024 * 1024 * 1024, args.timeout)
                if checksum != expected[(tracer, cap)]["sha256"]:
                    raise ValueError("Random FITS SHA256 disagrees with retained audit")
                for seed in SEEDS:
                    full, aud = select_fixed_random(
                        path, LEVELS[-1],
                        seed + 100 * ci + ti,
                        DIAGNOSTIC_Z[0], DIAGNOSTIC_Z[1],
                    )
                    perm = np.random.default_rng(
                        seed + 100 * ci + ti + 10000).permutation(LEVELS[-1])
                    seed_samples[seed][tracer] = {
                        n: subset(full, perm[:n]) for n in LEVELS
                    }
                    input_identity[(tracer, seed)] = {
                        "full_file_sha256": checksum,
                        "full_file_nbytes": nbytes,
                        "sample_size": aud["sample_size"],
                        "sample_arrays_sha256": aud["sample_arrays_sha256"],
                        "candidate_eligible_rows": aud["candidate_eligible_rows"],
                    }
                print(f"EBOSS_RR_INPUT_OK cap={cap} tracer={tracer} "
                      f"eligible={aud['candidate_eligible_rows']} sha256={checksum}",
                      flush=True)
            except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                failures.append(f"{cap}/{tracer}: {exc}")
                print("EBOSS_RR_INPUT_ERROR", failures[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()
        if not all(set(seed_samples[s]) == {"LRG", "ELG"} for s in SEEDS):
            continue
        sample_reports, largest = [], {}
        for seed in SEEDS:
            per_level, raw = [], {}
            for n in LEVELS:
                counts, meta = corrfunc_rr(
                    seed_samples[seed]["LRG"][n],
                    seed_samples[seed]["ELG"][n],
                    sedges, muedges, PRIMARY_GEOMETRY, args.threads)
                raw[n] = counts
                per_level.append({"n_randoms_per_tracer": n, **meta})
                print(f"EBOSS_RR_PILOT {cap} seed={seed} n={n} "
                      f"fraction={meta['normalized_rr_sum_in_test_bins']:.8g}",
                      flush=True)
            maximum = raw[LEVELS[-1]]
            for item in per_level:
                item["distance_to_largest_same_seed"] = histogram_distance(
                    raw[item["n_randoms_per_tracer"]], maximum)
            sample_reports.append({"seed": seed, "levels": per_level})
            largest[seed] = maximum
        between = histogram_distance(largest[SEEDS[0]], largest[SEEDS[1]])
        secondary, meta_secondary = corrfunc_rr(
            seed_samples[SEEDS[0]]["LRG"][LEVELS[-1]],
            seed_samples[SEEDS[0]]["ELG"][LEVELS[-1]],
            sedges, muedges, SECONDARY_GEOMETRY, args.threads,
        )
        cosmology = histogram_distance(secondary, largest[SEEDS[0]])
        record = {
            "cap": cap, "file_and_sample_evidence": [
                {"tracer": t, "seed": s, **input_identity[(t, s)]}
                for t in ("LRG", "ELG") for s in SEEDS
            ],
            "samples": sample_reports,
            "between_seed_comparison_at_maximum": between,
            "secondary_distance_comparison_same_maximum_sample": {
                "primary": PRIMARY_GEOMETRY,
                "secondary": SECONDARY_GEOMETRY,
                "relative_histogram_difference": cosmology,
                "secondary_counter": meta_secondary,
            },
            "largest_seed_histograms_normalized": {
                str(seed): largest[seed].tolist() for seed in SEEDS
            },
        }
        records.append(record)
        print(f"EBOSS_RR_PILOT_CAP {cap} "
              f"between_seed_l1={between['normalized_rr_l1_over_larger_sample_l1']:.6g} "
              f"distance_variant_l1="
              f"{cosmology['normalized_rr_l1_over_larger_sample_l1']:.6g}",
              flush=True)
    complete = len(records) == len(rng_names) and not failures
    report = {
        "study": "eBOSS random-only compiled RR high-z density and geometry pilot",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "pilot_random_only_rr_checked" if complete else "partial",
        "tracers": ["eBOSS_LRG", "eBOSS_ELG"],
        "galactic_caps": list(rng_names),
        "diagnostic_z_interval": list(DIAGNOSTIC_Z),
        "redshift_and_pair_bins_frozen_for_real_data_inference": False,
        "nested_random_counts_per_tracer": list(LEVELS),
        "fixed_seeds": list(SEEDS),
        "separation_edges_mpc_over_h": sedges.tolist(),
        "mu_edges": muedges.tolist(),
        "theta_min_deg": THETA_MIN_DEG,
        "orientation": "LRG->ELG",
        "los": "midpoint",
        "weights": "product WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP",
        "input_gate": "candidate numerical-zero rule and public random SHA256",
        "fiducial_primary": DISTANCE_CONVENTIONS[PRIMARY_GEOMETRY],
        "fiducial_secondary": DISTANCE_CONVENTIONS[SECONDARY_GEOMETRY],
        "distance_background": "flat matter+Lambda, Tcmb0=0, no radiation or neutrino correction",
        "cases": records,
        "errors": failures,
        "observed_galaxy_data_read": False,
        "observed_odd_vector_read": False,
        "galaxy_pairs_computed": False,
        "full_density_RR_window_computed": False,
        "window_even_to_odd_leakage_closed": False,
        "mock_covariance_computed": False,
        "note": (
            "A compiled Corrfunc pilot at fixed, input-only redshift and "
            "sampling levels; full random-density convergence and final "
            "window convolution are separate gates. No sample size, z cut, "
            "cosmology or random weight is selected by an odd-sector statistic."
        ),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("EBOSS_CORRFUNC_RR_PILOT", target, report["status"], flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
