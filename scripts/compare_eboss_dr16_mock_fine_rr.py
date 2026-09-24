#!/usr/bin/env python3
"""Matched eBOSS realistic-mock 0001 fine-RR comparison, randoms only.

Read the four SHA-verified matched realistic EZmock random FITS catalogues,
apply the already audited input-only numerical-zero and multiplicative
weight rule, and compare their fine, weighted LRG->ELG RR(s,mu) histograms
with the public observed-random *RR* histograms from the successful high-z
fine-geometry workflow.

This is a one-realization selection/window-input diagnostic, not a
validated cross-covariance, template tuning, galaxy odd vector or a test
of the wake response. No observed or mock galaxy positions are read.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from audit_eboss_dr16_full_random_rr import ROOT, compile_rr, coordinates
from audit_eboss_dr16_rr_density import histogram_distance
from audit_eboss_dr16_fine_rr_highz import (
    fine_edges, rebin_to_coarse, odd_rr_geometry,
)
from build_eboss_dr16_conditional_window import (
    OUT_ELLS, IN_ELLS, build_blocks, closure_summary,
    CONSTANT_CLOSURE_TOL, probe_even_to_odd,
)
from eboss_dr16_fiducial import (
    PRIMARY_GEOMETRY, WEIGHT_COLUMNS, validated_weight_product,
)
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry

MOCK_REFERENCE = ROOT / (
    "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
)
MOCK_ID = 1
Z_CANDIDATE = (0.6, 1.0)
Z_HIGH = (0.9, 1.0)
MU_BINS = 240


def extract_mock_high_z(path: Path, reference: dict,
                        chunk_rows: int = 150000):
    values = [[], [], [], []]
    count_all = 0
    count_excluded = 0
    count_high = 0
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected one realistic mock random BINTABLE")
        tab = tables[0]
        if int(tab.header["NAXIS2"]) != reference["header_rows"]:
            raise ValueError("Mock random FITS row count changed")
        needed = {"RA", "DEC", "Z", *WEIGHT_COLUMNS}
        if not needed.issubset(set(tab.columns.names)):
            raise ValueError("Mock random missing required weight or position columns")
        for begin in range(0, int(tab.header["NAXIS2"]), chunk_rows):
            batch = tab.data[begin:begin + chunk_rows]
            ra = np.asarray(batch["RA"], dtype="f8")
            dec = np.asarray(batch["DEC"], dtype="f8")
            z = np.asarray(batch["Z"], dtype="f8")
            valid = (
                np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90)
            )
            if not np.all(valid):
                raise ValueError("Unhealthy coordinate/redshift in matched random")
            candidate = (z >= Z_CANDIDATE[0]) & (z < Z_CANDIDATE[1])
            count_all += int(np.count_nonzero(candidate))
            if not np.any(candidate):
                continue
            weights = {
                col: np.asarray(batch[col][candidate], dtype="f8")
                for col in WEIGHT_COLUMNS
            }
            total, keep = validated_weight_product(weights)
            count_excluded += int(np.count_nonzero(~keep))
            high = (
                (z[candidate] >= Z_HIGH[0])
                & (z[candidate] < Z_HIGH[1]) & keep
            )
            if np.any(high):
                for bucket, field in zip(
                    values,
                    (ra[candidate][high], dec[candidate][high],
                     z[candidate][high], total[high]),
                ):
                    bucket.append(np.asarray(field, dtype="f8"))
                count_high += int(np.count_nonzero(high))
    if (count_all != reference["candidate_random_rows"]
            or count_excluded != reference["candidate_zero_weight_exclusions"]):
        raise ValueError("Mock random candidate/exclusion counts changed")
    if not count_high:
        raise ValueError("No high-z random objects survived the fixed gate")
    result = tuple(np.concatenate(parts) for parts in values)
    return result, {
        "candidate_rows": count_all,
        "candidate_zero_exclusions": count_excluded,
        "high_z_retained_rows": count_high,
        "high_z_weight_sum": float(np.sum(result[3], dtype="f8")),
    }


def operator_comparison(
    observed: np.ndarray, mock: np.ndarray,
    fine: np.ndarray, coarse: np.ndarray, mu: np.ndarray,
) -> dict:
    obs_blocks = build_blocks(observed, fine, coarse, mu)
    mock_blocks = build_blocks(mock, fine, coarse, mu)
    for label, blocks in (("observed", obs_blocks), ("mock", mock_blocks)):
        worst = max(
            entry["max_constant_closure_abs_error"]
            for entry in closure_summary(blocks).values()
        )
        if worst > CONSTANT_CLOSURE_TOL:
            raise ValueError(f"{label} constant-theory window closure failed")
    result = {}
    for o in OUT_ELLS:
        diag = float(np.linalg.norm(obs_blocks[(o, o)]))
        if diag <= 0:
            raise ValueError("Degenerate reference diagonal response")
        for i in IN_ELLS:
            key = f"{o}<-{i}"
            a, b = obs_blocks[(o, i)], mock_blocks[(o, i)]
            row = {
                "mock_minus_observed_max_absolute_cell": float(
                    np.max(np.abs(b - a))),
                "frobenius_difference_over_observed_odd_diagonal": float(
                    np.linalg.norm(b - a) / diag),
            }
            if i in (0, 2, 4):
                observed_probe = probe_even_to_odd(
                    obs_blocks, fine, coarse, o, i)
                mock_probe = probe_even_to_odd(
                    mock_blocks, fine, coarse, o, i)
                row["predeclared_linear_radial_probe"] = {
                    "observed_input_rr_response": observed_probe.tolist(),
                    "matched_mock_input_rr_response": mock_probe.tolist(),
                    "max_abs_mock_minus_observed": float(
                        np.max(np.abs(mock_probe - observed_probe))),
                }
            result[key] = row
    return result


def self_test() -> None:
    rng = np.random.default_rng(2609)
    fine = fine_edges()
    coarse = np.array([20., 40., 60., 80., 100., 120., 140.])
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, MU_BINS + 1)
    base = rng.uniform(0.3, 2, (len(fine) - 1, MU_BINS))
    scaled = 1.9 * base
    compare = operator_comparison(base, scaled, fine, coarse, mu)
    for row in compare.values():
        assert row["mock_minus_observed_max_absolute_cell"] < 1e-12
        if "predeclared_linear_radial_probe" in row:
            assert row["predeclared_linear_radial_probe"][
                "max_abs_mock_minus_observed"] < 1e-12
    assert np.max(np.abs(
        rebin_to_coarse(base, fine, coarse) -
        rebin_to_coarse(base, fine, coarse))) == 0
    print("EBOSS_MOCK_FINE_RR_NORMALIZATION_AND_OPERATOR_SELF_TEST_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--observed-rr", help="Successful observed fine-RR NPZ")
    ap.add_argument("--observed-audit", help="Matching observed fine-RR JSON")
    ap.add_argument("--out", default="eboss_workspace/matched_mock_fine_rr.json")
    ap.add_argument("--cache-dir", default="eboss_workspace/matched_mock_fine_rr")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.observed_rr or not args.observed_audit:
        ap.error("Successful observed high-z fine-RR NPZ and JSON are required")
    if args.threads < 1 or args.timeout <= 0:
        ap.error("Invalid thread count or download timeout")
    provenance = json.loads(Path(args.observed_audit).read_text())
    if (provenance["status"] != "fine_RR_high_z_geometric_closure_passed"
            or provenance["observed_odd_data_vector_read"]):
        raise ValueError("Observed input must be validated random-only fine RR")
    expected_source = os.environ.get("SOURCE_WORKFLOW_SHA")
    if expected_source and provenance.get("revision_commit") != expected_source:
        raise ValueError("Fine-RR source artifact SHA does not match run event")
    release = json.loads(MOCK_REFERENCE.read_text())
    if release["status"] != "mock_random_sample_selection_compatible":
        raise ValueError("The matched mock random selection audit is incomplete")
    known = {
        (int(x["id"]), x["tracer"], x["cap"]): x
        for x in release["mock_random_catalogues"]
    }
    arr = np.load(args.observed_rr, allow_pickle=False)
    fine = np.asarray(arr["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(arr["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(arr["muedges"], dtype="f8")
    if len(fine) != 121 or len(mu) != 241 or len(coarse) != 7:
        raise ValueError("Observed fine-RR bin geometry differs from fixed pilot")
    arrays = {"fine_sedges_mpc_over_h": fine, "muedges": mu}
    cases, errors = [], []
    root = Path(args.cache_dir)
    for cap in ("NGC", "SGC"):
        cats, sources = {}, {}
        for tracer in ("LRG", "ELG"):
            item = known[(MOCK_ID, tracer, cap)]
            relative = mock_path(
                "eBOSS_" + tracer, cap, "ran", MOCK_ID)
            path = root / Path(relative).name
            try:
                sha, size = fetch_with_retry(
                    MOCK_BASE + relative, path,
                    512 * 1024 * 1024, args.timeout)
                if sha != item["compressed_file_sha256"]:
                    raise ValueError("Matched mock random SHA256 changed")
                cat, meta = extract_mock_high_z(path, item)
                cats[tracer] = cat
                sources[tracer] = {
                    **meta, "compressed_file_sha256": sha,
                    "compressed_bytes": size,
                }
                print("MOCK_FINE_RR_INPUT_OK", cap, tracer,
                      meta["high_z_retained_rows"], sha, flush=True)
            except (OSError, ValueError, RuntimeError, KeyError,
                    MemoryError) as exc:
                errors.append(f"{cap}/{tracer}: {exc}")
                print("MOCK_FINE_RR_INPUT_ERROR", errors[-1], flush=True)
            finally:
                if path.exists():
                    path.unlink()
        if set(cats) != {"LRG", "ELG"}:
            continue
        try:
            observed = np.asarray(arr[f"fine_rr_{cap}_full"], dtype="f8")
            observed_norm = float(
                np.asarray(arr[f"fine_norm_{cap}_full"]).item())
            raw, meta = compile_rr(
                coordinates(cats["LRG"], PRIMARY_GEOMETRY),
                coordinates(cats["ELG"], PRIMARY_GEOMETRY),
                fine, mu, args.threads,
            )
            arrays[f"mock_fine_rr_{cap}_0001"] = raw
            arrays[f"mock_fine_norm_{cap}_0001"] = np.array(
                [meta["pair_weight_normalization"]], dtype="f8")
            comparison = histogram_distance(
                raw / meta["pair_weight_normalization"],
                observed / observed_norm)
            count_comparison = {
                "observed_positive_cells": int(np.count_nonzero(observed)),
                "mock_positive_cells": int(np.count_nonzero(raw)),
                "total_cells": int(raw.size),
                "normalized_RR_histogram_difference": comparison,
            }
            operator = operator_comparison(
                observed, raw, fine, coarse, mu)
            odd = {}
            obs_p = odd_rr_geometry(
                rebin_to_coarse(observed, fine, coarse), mu)
            mock_p = odd_rr_geometry(
                rebin_to_coarse(raw, fine, coarse), mu)
            for ell in ("1", "3"):
                odd[ell] = {
                    "observed_raw_RR_angular_moment": obs_p[ell].tolist(),
                    "mock_raw_RR_angular_moment": mock_p[ell].tolist(),
                    "max_abs_mock_minus_observed": float(np.max(
                        np.abs(obs_p[ell] - mock_p[ell]))),
                }
            cases.append({
                "cap": cap, "mock_realization_id": MOCK_ID,
                "inputs": sources, "compiled_mock_counter": meta,
                "fine_rr_comparison": count_comparison,
                "conditional_response_difference": operator,
                "raw_rr_odd_geometry_only": odd,
            })
            print("MATCHED_MOCK_FINE_RR_COMPARISON", cap,
                  "relative_l1",
                  comparison["normalized_rr_l1_over_larger_sample_l1"],
                  "max_raw_RR_P1_shift",
                  odd["1"]["max_abs_mock_minus_observed"], flush=True)
        except (OSError, ValueError, RuntimeError, KeyError,
                MemoryError) as exc:
            errors.append(f"{cap}: {exc}")
            print("MOCK_FINE_RR_COMPARISON_ERROR", errors[-1], flush=True)
    complete = len(cases) == 2 and not errors
    output = {
        "study": "One-realization observed-vs-matched-mock fine RR selection audit",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "source_observed_fine_rr_run": os.environ.get("SOURCE_WORKFLOW_RUN"),
        "status": "mock_0001_fine_rr_comparison_complete" if complete else "partial",
        "mock_reference": str(MOCK_REFERENCE.relative_to(ROOT)),
        "realization_id": MOCK_ID, "tracers": ["LRG", "ELG"],
        "candidate_redshift_slice": list(Z_HIGH),
        "candidate_slice_frozen_for_inference": False,
        "separation_step_mpc_over_h": 1.0,
        "mu_bins": MU_BINS,
        "distance_profile": PRIMARY_GEOMETRY,
        "orientation": "LRG->ELG",
        "line_of_sight": "midpoint",
        "cases": cases, "errors": errors,
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "wake_template_fit_or_optimization": False,
        "ensemble_covariance_validated": False,
        "full_mock_window_calibrated": False,
        "note": (
            "One selected realistic realization is an input-only geometry "
            "comparison, not an ensemble uncertainty estimate. Any observed-"
            "vs-mock difference is reported without optimizing the tracer "
            "selection, fiducial cosmology, redshift bins or wake template."
        ),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    if arrays:
        np.savez_compressed(target.with_suffix(".npz"), **arrays)
    print("EBOSS_MOCK_FINE_RR_AUDIT", output["status"], target, flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
