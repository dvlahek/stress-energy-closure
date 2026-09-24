#!/usr/bin/env python3
"""Check matched realistic eBOSS LRG/ELG mock randoms before odd-sector analysis.

The sample is fixed at realization IDs 0001, 0500 and 1000. The script reads
the matching public random FITS files to inspect redshift/weight health and
diagnostic LRG-ELG joint random support. It never computes a pair-count,
multipole, physical template, likelihood, covariance or real odd data vector.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from inspect_eboss_dr16_joint_randoms import grid_indices, N_PIX, RA_STEP, SIN_DEC_STEP
from inspect_eboss_dr16_mock_headers import BASE, TRACERS, CAPS, mock_path
from inspect_eboss_dr16_mock_selection import (
    CANDIDATE_LO, CANDIDATE_HI, NUMERICAL_ZERO_WEIGHT_TOL,
    fetch_with_retry,
)
from inspect_eboss_dr16_selection import WEIGHTS, Z_EDGES

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "source_data/eboss_dr16_realistic_mock_selection_audit_2026-09-24.json"
SUPPORT_THRESHOLDS = (1, 5, 10, 15, 25)


def blank_weights() -> dict:
    return {
        name: {
            "all_nonfinite": 0, "all_nonpositive": 0,
            "candidate_nonfinite": 0, "candidate_nonpositive": 0,
            "candidate_abs_le_1e_20": 0, "candidate_abs_le_1e_12": 0,
            "candidate_significantly_negative_below_minus_1e_12": 0,
            "candidate_min_abs_above_1e_12": None,
        }
        for name in WEIGHTS
    }


def check_weights(batch, candidate: np.ndarray, diagnostic: dict
                  ) -> np.ndarray:
    retained = candidate.copy()
    for name in WEIGHTS:
        w = np.asarray(batch[name], dtype=np.float64)
        good = np.isfinite(w)
        d = diagnostic[name]
        d["all_nonfinite"] += int(np.count_nonzero(~good))
        d["all_nonpositive"] += int(np.count_nonzero(good & (w <= 0)))
        d["candidate_nonfinite"] += int(np.count_nonzero(candidate & ~good))
        d["candidate_nonpositive"] += int(np.count_nonzero(candidate & good & (w <= 0)))
        d["candidate_abs_le_1e_20"] += int(np.count_nonzero(
            candidate & good & (np.abs(w) <= 1e-20)))
        d["candidate_abs_le_1e_12"] += int(np.count_nonzero(
            candidate & good & (np.abs(w) <= NUMERICAL_ZERO_WEIGHT_TOL)))
        d["candidate_significantly_negative_below_minus_1e_12"] += int(
            np.count_nonzero(candidate & good & (w < -NUMERICAL_ZERO_WEIGHT_TOL)))
        other = candidate & good & (np.abs(w) > NUMERICAL_ZERO_WEIGHT_TOL)
        if np.any(other):
            m = float(np.min(np.abs(w[other])))
            d["candidate_min_abs_above_1e_12"] = (
                m if d["candidate_min_abs_above_1e_12"] is None
                else min(m, d["candidate_min_abs_above_1e_12"]))
        retained &= good & (w > (NUMERICAL_ZERO_WEIGHT_TOL
                                  if name == "WEIGHT_SYSTOT" else 0.0))
    return retained


def inspect_mock_random(path: Path, tracer: str, cap: str, rid: int,
                        checksum: str, nbytes: int, chunk_rows: int
                        ) -> tuple[dict, np.ndarray, list[np.ndarray]]:
    total_candidate = 0
    retained_candidate = 0
    invalid_coordinate_z = 0
    redshift_raw = np.zeros(len(Z_EDGES) - 1, dtype=np.int64)
    redshift_retained = np.zeros(len(Z_EDGES) - 1, dtype=np.int64)
    candidate_grid = np.zeros(N_PIX, dtype=np.int64)
    zgrids = [np.zeros(N_PIX, dtype=np.int64) for _ in range(4)]
    wdiag = blank_weights()

    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected one binary-table HDU in " + path.name)
        hdu = tables[0]
        missing = {"RA", "DEC", "Z", *WEIGHTS} - set(hdu.columns.names)
        if missing:
            raise ValueError(f"Missing mock random columns: {sorted(missing)}")
        declared = int(hdu.header["NAXIS2"])
        if declared <= 0:
            raise ValueError("The random FITS table has no rows")
        for first in range(0, declared, chunk_rows):
            last = min(declared, first + chunk_rows)
            batch = hdu.data[first:last]
            ra = np.asarray(batch["RA"], dtype=np.float64)
            dec = np.asarray(batch["DEC"], dtype=np.float64)
            z = np.asarray(batch["Z"], dtype=np.float64)
            valid = (np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                     & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
            invalid_coordinate_z += int(np.count_nonzero(~valid))
            candidate = valid & (z >= CANDIDATE_LO) & (z < CANDIDATE_HI)
            total_candidate += int(np.count_nonzero(candidate))
            retained = check_weights(batch, candidate, wdiag)
            retained_candidate += int(np.count_nonzero(retained))
            if np.any(retained):
                candidate_grid += np.bincount(
                    grid_indices(ra[retained], dec[retained]), minlength=N_PIX)
            for i, (lo, hi) in enumerate(zip(Z_EDGES[:-1], Z_EDGES[1:])):
                within = valid & (z >= lo) & (z < hi)
                redshift_raw[i] += int(np.count_nonzero(within))
                if i < 4:
                    accepted = retained & within
                    redshift_retained[i] += int(np.count_nonzero(accepted))
                    if np.any(accepted):
                        zgrids[i] += np.bincount(
                            grid_indices(ra[accepted], dec[accepted]),
                            minlength=N_PIX)

    report = {
        "realization_id": rid, "tracer": tracer, "cap": cap,
        "filename": path.name,
        "compressed_file_sha256": checksum,
        "compressed_bytes": nbytes, "header_rows": declared,
        "invalid_coordinate_z_rows": invalid_coordinate_z,
        "candidate_random_rows": total_candidate,
        "candidate_retained_after_weight_convention": retained_candidate,
        "candidate_excluded_by_weight_convention": (
            total_candidate - retained_candidate),
        "redshift_bin_edges": list(Z_EDGES),
        "redshift_bin_raw_rows": redshift_raw.tolist(),
        "redshift_bin_retained_rows": redshift_retained.tolist(),
        "candidate_sky_grid_support_pixels": {
            str(t): int(np.count_nonzero(candidate_grid >= t))
            for t in SUPPORT_THRESHOLDS
        },
        "weight_column_diagnostics": wdiag,
    }
    return report, candidate_grid, zgrids


def joint_support(a: np.ndarray, b: np.ndarray) -> dict:
    result = {}
    for threshold in SUPPORT_THRESHOLDS:
        la, eb = a >= threshold, b >= threshold
        inter = la & eb
        result[str(threshold)] = {
            "lrg_supported_pixels": int(np.count_nonzero(la)),
            "elg_supported_pixels": int(np.count_nonzero(eb)),
            "joint_supported_pixels": int(np.count_nonzero(inter)),
            "joint_fraction_of_elg_supported_pixels": (
                float(np.count_nonzero(inter) / np.count_nonzero(eb))
                if np.any(eb) else None),
        }
    return result


def sample_is_healthy(records: list[dict]) -> bool:
    for item in records:
        if item["invalid_coordinate_z_rows"] != 0:
            return False
        if item["candidate_random_rows"] <= 0:
            return False
        for name, d in item["weight_column_diagnostics"].items():
            if d["candidate_nonfinite"] != 0:
                return False
            if name == "WEIGHT_SYSTOT":
                if d["candidate_significantly_negative_below_minus_1e_12"]:
                    return False
                if d["candidate_abs_le_1e_20"] != d["candidate_abs_le_1e_12"]:
                    return False
            elif d["candidate_nonpositive"] != 0:
                return False
    return True


def self_test() -> None:
    fake = np.zeros(5, dtype=[(n, "f8") for n in WEIGHTS])
    for name in WEIGHTS:
        fake[name] = 1.0
    fake["WEIGHT_SYSTOT"] = [1.0, 1e-30, -1e-30, -1e-11, 1.4]
    d = blank_weights()
    kept = check_weights(fake, np.ones(5, dtype=bool), d)
    assert kept.tolist() == [True, False, False, False, True]
    assert d["WEIGHT_SYSTOT"]["candidate_abs_le_1e_12"] == 2
    assert d["WEIGHT_SYSTOT"][
        "candidate_significantly_negative_below_minus_1e_12"] == 1
    x = np.array([0, 5, 15, 25])
    y = np.array([25, 5, 14, 0])
    assert joint_support(x, y)["5"]["joint_supported_pixels"] == 2
    assert joint_support(x, y)["25"]["joint_supported_pixels"] == 0
    print("Matched mock random weight and gridded-support self-tests passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/mock_random_selection_audit.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/mock_random_catalogues")
    parser.add_argument("--ids", default="1,500,1000")
    parser.add_argument("--max-file-mib", type=int, default=512)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--chunk-rows", type=int, default=150000)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    ids = sorted(set(int(v.strip()) for v in args.ids.split(",")))
    if (not ids or any(i < 1 or i > 1000 for i in ids)
            or args.max_file_mib < 1 or args.timeout <= 0
            or args.chunk_rows < 1000):
        parser.error("Invalid IDs, limits or chunk size")
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    if ids != reference["mock_realization_ids"]:
        parser.error("Inspect the fixed mock data sample IDs 0001, 0500 and 1000")
    data_counts = reference["candidate_rows_raw_retained_excluded"]
    results, errors = [], []
    sky, by_z = {}, {}
    root = Path(args.cache_dir)

    for rid in ids:
        for cap in CAPS:
            for tracer in TRACERS:
                kind = "LRG" if tracer == "eBOSS_LRG" else "ELG"
                relative = mock_path(tracer, cap, "ran", rid)
                local = root / Path(relative).name
                try:
                    sha, nbytes = fetch_with_retry(
                        BASE + relative, local, args.max_file_mib * 1024 * 1024,
                        args.timeout)
                    record, counts, bins = inspect_mock_random(
                        local, kind, cap, rid, sha, nbytes, args.chunk_rows)
                    matched_data = data_counts[f"{rid:04d}"][f"{kind}_{cap}"]
                    record["matched_mock_data_candidate_rows"] = matched_data[0]
                    record["matched_mock_data_retained_rows"] = matched_data[1]
                    record["candidate_random_to_data_ratio"] = (
                        record["candidate_retained_after_weight_convention"] /
                        matched_data[1] if matched_data[1] else None)
                    sky[(rid, cap, kind)] = counts
                    by_z[(rid, cap, kind)] = bins
                    results.append(record)
                    print(f"MOCK_RANDOM_OK id={rid:04d} {kind} {cap} "
                          f"rows={record['header_rows']} "
                          f"candidate={record['candidate_random_rows']} "
                          f"excluded={record['candidate_excluded_by_weight_convention']} "
                          f"sha256={sha}", flush=True)
                except (OSError, ValueError, KeyError, RuntimeError,
                        MemoryError) as exc:
                    errors.append(f"{relative}: {exc}")
                    print("MOCK_RANDOM_ERROR", errors[-1], flush=True)
                finally:
                    if local.exists():
                        local.unlink()

    joint = []
    for rid in ids:
        for cap in CAPS:
            lrg, elg = (rid, cap, "LRG"), (rid, cap, "ELG")
            if lrg not in sky or elg not in sky:
                continue
            joint.append({
                "realization_id": rid,
                "cap": cap,
                "candidate_joint_random_support": joint_support(sky[lrg], sky[elg]),
                "redshift_bin_joint_random_support": [
                    {
                        "zlo": Z_EDGES[i], "zhi": Z_EDGES[i + 1],
                        "support": joint_support(by_z[lrg][i], by_z[elg][i]),
                    }
                    for i in range(4)
                ],
            })
    expected = len(ids) * len(CAPS) * len(TRACERS)
    complete = len(results) == expected and not errors and len(joint) == len(ids) * 2
    healthy = sample_is_healthy(results)
    status = (
        "mock_random_sample_selection_compatible"
        if complete and healthy else "mock_random_weight_or_selection_unresolved"
        if complete else "partial")
    report = {
        "study": "eBOSS DR16 realistic LRG-ELG matched random selection audit",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "release_url": BASE, "realization_ids": ids, "status": status,
        "candidate_z_interval": [CANDIDATE_LO, CANDIDATE_HI],
        "candidate_interval_frozen": False,
        "numeric_zero_systot_abs_tolerance": NUMERICAL_ZERO_WEIGHT_TOL,
        "sky_grid": {
            "ra_step_deg": RA_STEP, "sin_dec_step": SIN_DEC_STEP,
            "support_thresholds_randoms_per_tracer_per_cell":
                list(SUPPORT_THRESHOLDS),
        },
        "random_catalogues": results,
        "joint_support_by_id_and_cap": joint,
        "errors": errors,
        "mock_data_reference": str(REFERENCE.relative_to(ROOT)),
        "real_data_odd_vector_read": False, "pair_counts_computed": False,
        "pair_window_computed": False, "mock_covariance_computed": False,
        "mock_ensemble_all_1000_checked": False,
        "note": (
            "This is an input-only sample of matched random FITS catalogues. "
            "Gridded support is a resolution/threshold diagnostic, not an exact "
            "angular mask or survey-window matrix. The full ensemble's weights, "
            "cross-tracer covariance and estimator closure remain unvalidated."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("MOCK_RANDOM_AUDIT", out, status, flush=True)
    for record in joint:
        k = record["candidate_joint_random_support"]["15"]
        print(f"MOCK_RANDOM_JOINT {record['realization_id']:04d} {record['cap']} "
              f"threshold=15 joint={k['joint_supported_pixels']} "
              f"ELG={k['elg_supported_pixels']}", flush=True)
    return 0 if status == "mock_random_sample_selection_compatible" else 2


if __name__ == "__main__":
    raise SystemExit(main())
