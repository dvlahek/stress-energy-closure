#!/usr/bin/env python3
"""Compare observed and realistic-mock LRG/ELG random sky support.

Only public random catalogues and previously fixed mock IDs enter this
pre-unblinding diagnostic. Mock pixel counts are normalized by the ratio of
real to mock random totals, separately by tracer, cap and redshift interval.
No galaxy pair counting, odd multipole, wake fit or covariance is calculated.
The grid intersection is not an exact survey veto mask or pair window.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from inspect_eboss_dr16_joint_randoms import (
    RANDOMS, RA_STEP, SIN_DEC_STEP, examine as inspect_real,
)
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_randoms import (
    inspect_mock_random, SUPPORT_THRESHOLDS,
)
from inspect_eboss_dr16_mock_selection import fetch_with_retry
from inspect_eboss_dr16_selection import BASE as REAL_BASE, Z_EDGES, download

ROOT = Path(__file__).resolve().parents[1]
REAL_REFERENCE = (
    ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json"
)
MOCK_REFERENCE = (
    ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
)
THRESHOLDS = (1, 5, 10, 15, 25, 50)


def counts_overlap(real: np.ndarray, mock: np.ndarray,
                   target_scale: float, threshold: int) -> dict:
    if target_scale <= 0 or not np.isfinite(target_scale):
        raise ValueError("Normalization must be positive and finite")
    if real.shape != mock.shape:
        raise ValueError("Pixel grids do not share geometry")
    real_support = real >= threshold
    mock_support = (mock.astype(np.float64) * target_scale) >= threshold
    return masks_overlap(real_support, mock_support)


def masks_overlap(real_support: np.ndarray, mock_support: np.ndarray) -> dict:
    a, b = int(np.count_nonzero(real_support)), int(np.count_nonzero(mock_support))
    common = int(np.count_nonzero(real_support & mock_support))
    union = int(np.count_nonzero(real_support | mock_support))
    return {
        "observed_supported_pixels": a,
        "mock_normalized_supported_pixels": b,
        "intersecting_pixels": common,
        "union_pixels": union,
        "jaccard": float(common / union) if union else None,
        "observed_support_covered_by_mock": float(common / a) if a else None,
        "mock_support_inside_observed": float(common / b) if b else None,
        "observed_only_pixels": a - common,
        "mock_only_pixels": b - common,
    }


def tracer_comparison(real_counts: np.ndarray, mock_counts: np.ndarray,
                      real_rows: int, mock_rows: int) -> dict:
    if real_rows <= 0 or mock_rows <= 0:
        raise ValueError("Cannot compare empty random selection")
    factor = real_rows / mock_rows
    return {
        "observed_random_rows": real_rows,
        "mock_random_rows_after_input_rule": mock_rows,
        "random_density_scale_real_over_mock": factor,
        "thresholded_grid_overlap": {
            str(t): counts_overlap(real_counts, mock_counts, factor, t)
            for t in THRESHOLDS
        },
    }


def joint_comparison(real_lrg: np.ndarray, real_elg: np.ndarray,
                     mock_lrg: np.ndarray, mock_elg: np.ndarray,
                     real_rows: dict, mock_rows: dict) -> dict:
    if min(*real_rows.values(), *mock_rows.values()) <= 0:
        raise ValueError("A joint random normalization has zero rows")
    real_over_mock = {
        key: real_rows[key] / mock_rows[key] for key in ("LRG", "ELG")
    }
    tests = {}
    for threshold in THRESHOLDS:
        real_joint = ((real_lrg >= threshold) & (real_elg >= threshold))
        mock_joint = (
            (mock_lrg.astype(np.float64) * real_over_mock["LRG"] >= threshold)
            & (mock_elg.astype(np.float64) * real_over_mock["ELG"] >= threshold)
        )
        tests[str(threshold)] = masks_overlap(real_joint, mock_joint)
    return {
        "tracer_density_scale_real_over_mock": real_over_mock,
        "thresholded_joint_grid_overlap": tests,
    }


def self_test() -> None:
    real = np.array([0, 3, 4, 0], dtype=np.int64)
    mock = np.array([0, 2, 5, 1], dtype=np.int64)
    comparison = counts_overlap(real, mock, 1.5, 3)
    assert comparison["jaccard"] == 1.0
    assert comparison["observed_supported_pixels"] == 2
    bad = masks_overlap(real >= 3, mock >= 5)
    assert bad["intersecting_pixels"] == 1
    assert bad["jaccard"] == 0.5
    print("Normalized random support overlap self-tests passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/real_mock_random_geometry.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/real_mock_randoms")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--chunk-rows", type=int, default=200000)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout <= 0 or args.chunk_rows < 1000:
        parser.error("Timeout must be positive; chunk size must be >= 1000")
    real_ref = json.loads(REAL_REFERENCE.read_text(encoding="utf-8"))
    mock_ref = json.loads(MOCK_REFERENCE.read_text(encoding="utf-8"))
    if real_ref["status"] != "four_random_catalogues_checked":
        raise ValueError("The full observed-random reference was not validated")
    if mock_ref["status"] != "mock_random_sample_selection_compatible":
        raise ValueError("The realistic-mock random input reference is not validated")
    ids = mock_ref["sample_realization_ids"]
    real_expect = {(x["tracer"], x["cap"]): x for x in real_ref["randoms"]}
    mock_expect = {
        (int(x["id"]), x["tracer"], x["cap"]): x
        for x in mock_ref["mock_random_catalogues"]
    }
    cache = Path(args.cache_dir)
    real_grids, real_zgrids, real_records = {}, {}, []
    mock_grids, mock_zgrids, mock_records = {}, {}, []
    errors = []

    for (kind, cap), (filename, expected_rows) in RANDOMS.items():
        path = cache / "observed" / filename
        try:
            digest, size = None, None
            path, digest, size = download(filename, path.parent, 1024 * 1024 * 1024,
                                          args.timeout)
            ref = real_expect[(kind, cap)]
            if digest != ref["sha256"] or expected_rows != ref["rows"]:
                raise ValueError("Observed random file differs from the proven release SHA/rows")
            record, counts, zgrids = inspect_real(
                path, kind, cap, expected_rows, digest, size, args.chunk_rows)
            if record["candidate_random_rows"] != ref["candidate_rows"]:
                raise ValueError("Observed random candidate redshift count changed")
            if (record["invalid_coordinate_z_rows"] or any(
                    d["nonfinite"] or d["nonpositive"]
                    for d in record["weight_column_diagnostics"].values())):
                raise ValueError("Observed random input health changed")
            real_records.append({
                "tracer": kind, "cap": cap, "sha256": digest, "bytes": size,
                "candidate_random_rows": record["candidate_random_rows"],
                "zbin_counts": record["redshift_bin_counts"],
            })
            real_grids[(kind, cap)] = counts
            real_zgrids[(kind, cap)] = zgrids
            print(f"REAL_RANDOM_GEOMETRY_OK {kind} {cap} "
                  f"candidate={record['candidate_random_rows']}", flush=True)
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            errors.append(f"{filename}: {exc}")
            print("REAL_RANDOM_GEOMETRY_ERROR", errors[-1], flush=True)
        finally:
            if path.exists():
                path.unlink()

    for rid in ids:
        for cap in ("NGC", "SGC"):
            for kind in ("LRG", "ELG"):
                tracer = "eBOSS_" + kind
                relative = mock_path(tracer, cap, "ran", rid)
                path = cache / "realistic" / Path(relative).name
                try:
                    digest, size = fetch_with_retry(
                        MOCK_BASE + relative, path, 512 * 1024 * 1024,
                        args.timeout)
                    ref = mock_expect[(rid, kind, cap)]
                    if digest != ref["compressed_file_sha256"]:
                        raise ValueError("Mock random compressed SHA differs from proven release")
                    record, counts, zgrids = inspect_mock_random(
                        path, kind, cap, rid, digest, size, args.chunk_rows)
                    if (record["header_rows"] != ref["header_rows"]
                            or record["candidate_random_rows"] !=
                            ref["candidate_random_rows"]
                            or record["candidate_excluded_by_weight_convention"] !=
                            ref["candidate_zero_weight_exclusions"]):
                        raise ValueError("Mock random row counts differ from audited sample")
                    mock_records.append({
                        "id": rid, "tracer": kind, "cap": cap, "sha256": digest,
                        "candidate_random_rows": record["candidate_random_rows"],
                        "candidate_retained_random_rows":
                            record["candidate_retained_after_weight_convention"],
                        "candidate_excluded_random_rows":
                            record["candidate_excluded_by_weight_convention"],
                        "zbin_retained_rows": record["redshift_bin_retained_rows"],
                    })
                    mock_grids[(rid, kind, cap)] = counts
                    mock_zgrids[(rid, kind, cap)] = zgrids
                    print(f"MOCK_RANDOM_GEOMETRY_OK {rid:04d} {kind} {cap} "
                          f"retained={record['candidate_retained_after_weight_convention']}",
                          flush=True)
                except (OSError, ValueError, KeyError, RuntimeError,
                        MemoryError) as exc:
                    errors.append(f"{relative}: {exc}")
                    print("MOCK_RANDOM_GEOMETRY_ERROR", errors[-1], flush=True)
                finally:
                    if path.exists():
                        path.unlink()

    observed = {(x["tracer"], x["cap"]): x for x in real_records}
    mock = {(x["id"], x["tracer"], x["cap"]): x for x in mock_records}
    comparisons = []
    for rid in ids:
        for cap in ("NGC", "SGC"):
            keys = [(rid, kind, cap) for kind in ("LRG", "ELG")]
            if not all(key in mock_grids for key in keys):
                continue
            if not all((kind, cap) in real_grids for kind in ("LRG", "ELG")):
                continue
            l, e = ("LRG", cap), ("ELG", cap)
            ml, me = (rid, "LRG", cap), (rid, "ELG", cap)
            real_rows = {k: observed[(k, cap)]["candidate_random_rows"]
                         for k in ("LRG", "ELG")}
            mock_rows = {k: mock[(rid, k, cap)]["candidate_retained_random_rows"]
                         for k in ("LRG", "ELG")}
            by_tracer = {
                kind: tracer_comparison(
                    real_grids[(kind, cap)], mock_grids[(rid, kind, cap)],
                    real_rows[kind], mock_rows[kind])
                for kind in ("LRG", "ELG")
            }
            by_z = []
            for i in range(4):
                real_bin_rows = {
                    k: observed[(k, cap)]["zbin_counts"][i]
                    for k in ("LRG", "ELG")
                }
                mock_bin_rows = {
                    k: mock[(rid, k, cap)]["zbin_retained_rows"][i]
                    for k in ("LRG", "ELG")
                }
                by_z.append({
                    "zlo": Z_EDGES[i], "zhi": Z_EDGES[i + 1],
                    "joint": joint_comparison(
                        real_zgrids[l][i], real_zgrids[e][i],
                        mock_zgrids[ml][i], mock_zgrids[me][i],
                        real_bin_rows, mock_bin_rows),
                })
            joint = joint_comparison(
                real_grids[l], real_grids[e], mock_grids[ml], mock_grids[me],
                real_rows, mock_rows)
            comparisons.append({
                "id": rid, "cap": cap, "tracer_selections": by_tracer,
                "joint_candidate_selection": joint, "joint_by_redshift_bin": by_z,
            })
            j = joint["thresholded_joint_grid_overlap"]["15"]
            print(f"REAL_MOCK_JOINT_GEOMETRY id={rid:04d} {cap} "
                  f"threshold=15 jaccard={j['jaccard']:.6f} "
                  f"observed={j['observed_supported_pixels']} "
                  f"mock={j['mock_normalized_supported_pixels']}", flush=True)

    complete = (
        len(real_records) == 4 and len(mock_records) == len(ids) * 4
        and len(comparisons) == len(ids) * 2 and not errors
    )
    report = {
        "study": "eBOSS DR16 realistic EZmock versus observed random-only geometry",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "random_geometry_diagnostic_complete" if complete else "partial",
        "fixed_real_reference": str(REAL_REFERENCE.relative_to(ROOT)),
        "fixed_mock_reference": str(MOCK_REFERENCE.relative_to(ROOT)),
        "fixed_mock_ids": ids,
        "candidate_interval_frozen": False,
        "sky_grid": {"ra_step_deg": RA_STEP, "sin_dec_step": SIN_DEC_STEP},
        "thresholds_in_real_equivalent_random_counts": list(THRESHOLDS),
        "normalization": (
            "Each mock random cell count is multiplied by the ratio of real "
            "to selected mock random totals, separately for tracer, cap and "
            "candidate or redshift interval. This matches sampling density "
            "for diagnostic grid occupancy only, not survey completeness."
        ),
        "real_random_file_checks": real_records,
        "mock_random_file_checks": mock_records,
        "geometry_comparisons": comparisons,
        "errors": errors,
        "odd_sector_data_read": False,
        "galaxy_pair_counts_computed": False,
        "survey_window_computed": False,
        "wake_template_tuned": False,
        "mock_covariance_computed": False,
        "full_mock_ensemble_checked": False,
        "scope_note": (
            "The Jaccard indices quantify threshold- and resolution-dependent "
            "random-grid support only. They do not prove matching continuous "
            "survey masks or correct even-to-odd leakage. All matching input "
            "SHA256 digests are checked against previously retained evidence."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("REAL_MOCK_RANDOM_GEOMETRY", out, report["status"], flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
