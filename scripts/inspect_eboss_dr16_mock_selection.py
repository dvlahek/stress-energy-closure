#!/usr/bin/env python3
"""Validate a stratified eBOSS realistic EZmock data-only selection sample.

This is a pre-unblinding input audit. It reads only mock data FITS catalogues,
not mock randoms or real odd-sector measurements. It compares mock and real
redshift distributions and reports weight and coordinate health; it does not
fit a signal, select bins or estimate a covariance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits

from inspect_eboss_dr16_mock_headers import BASE, TRACERS, CAPS, mock_path
from inspect_eboss_dr16_selection import WEIGHTS, Z_EDGES

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "source_data/eboss_dr16_data_selection_audit_2026-09-24.json"
CANDIDATE_LO, CANDIDATE_HI = 0.6, 1.0
NUMERICAL_ZERO_WEIGHT_TOL = 1e-12


def fetch_data(url: str, destination: Path, byte_limit: int, timeout: float
               ) -> tuple[str, int]:
    allowed, target = urlsplit(BASE), urlsplit(url)
    if (target.scheme != "https" or target.netloc != allowed.netloc
            or not target.path.startswith(allowed.path) or target.query):
        raise ValueError("Mock data URL is outside the approved release")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    sha, count = hashlib.sha256(), 0
    try:
        request = Request(url, headers={
            "User-Agent": "eBOSS-DR16-realistic-mock-selection/1.0",
            "Accept-Encoding": "identity",
        })
        with urlopen(request, timeout=timeout) as response, partial.open("wb") as stream:
            final = urlsplit(response.geturl())
            if (final.scheme != "https" or final.netloc != allowed.netloc
                    or not final.path.startswith(allowed.path)):
                raise ValueError("Mock data redirected outside the allowed release")
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                count += len(block)
                if count > byte_limit:
                    raise ValueError("Mock data exceed the bounded byte limit")
                stream.write(block)
                sha.update(block)
        if count == 0:
            raise ValueError("Empty mock catalogue")
        partial.replace(destination)
        return sha.hexdigest(), count
    finally:
        if partial.exists():
            partial.unlink()


def fetch_with_retry(url: str, destination: Path, byte_limit: int, timeout: float
                     ) -> tuple[str, int]:
    for attempt in range(3):
        try:
            return fetch_data(url, destination, byte_limit, timeout)
        except (OSError, HTTPError, URLError, TimeoutError) as exc:
            if attempt == 2:
                raise
            print(f"RETRY_MOCK_DOWNLOAD {attempt + 1}/2 {destination.name}: {exc}",
                  flush=True)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("Unreachable mock download retry state")


def inspect(path: Path, tracer: str, cap: str, rid: int,
            sha: str, nbytes: int) -> tuple[dict, set[int]]:
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError(f"Expected one BINTABLE: {path.name}")
        table = tables[0]
        required = {"RA", "DEC", "Z", *WEIGHTS}
        actual = set(table.columns.names)
        if not required.issubset(actual):
            raise ValueError(f"Missing mock columns: {sorted(required - actual)}")
        ra = np.asarray(table.data["RA"], dtype=np.float64)
        dec = np.asarray(table.data["DEC"], dtype=np.float64)
        z = np.asarray(table.data["Z"], dtype=np.float64)
        valid = (np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                 & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
        candidate = valid & (z >= CANDIDATE_LO) & (z < CANDIDATE_HI)
        zcounts = [
            int(np.count_nonzero(valid & (z >= lo) & (z < hi)))
            for lo, hi in zip(Z_EDGES[:-1], Z_EDGES[1:])
        ]
        wdiag = {}
        retained_candidate = candidate.copy()
        for name in WEIGHTS:
            w = np.asarray(table.data[name], dtype=np.float64)
            finite = np.isfinite(w)
            lower = NUMERICAL_ZERO_WEIGHT_TOL if name == "WEIGHT_SYSTOT" else 0.0
            retained_candidate &= finite & (w > lower)
            wdiag[name] = {
                "nonfinite_rows": int(np.count_nonzero(~finite)),
                "nonpositive_rows": int(np.count_nonzero(finite & (w <= 0))),
                "nonfinite_candidate_rows": int(np.count_nonzero(candidate & ~finite)),
                "nonpositive_candidate_rows": int(np.count_nonzero(
                    candidate & finite & (w <= 0))),
                "significant_negative_candidate_rows_below_minus_1e_12": int(
                    np.count_nonzero(candidate & finite &
                                     (w < -NUMERICAL_ZERO_WEIGHT_TOL))),
                "near_zero_candidate_rows_abs_le_1e_20": int(np.count_nonzero(
                    candidate & finite & (np.abs(w) <= 1e-20))),
                "near_zero_candidate_rows_abs_le_1e_12": int(np.count_nonzero(
                    candidate & finite & (np.abs(w) <= 1e-12))),
                "min_abs_candidate_weight_above_1e_12": (
                    float(np.min(np.abs(w[candidate & finite & (np.abs(w) > 1e-12)])))
                    if np.any(candidate & finite & (np.abs(w) > 1e-12)) else None),
                "nonpositive_outside_candidate_rows": int(np.count_nonzero(
                    (~candidate) & finite & (w <= 0))),
                "min_finite": float(np.min(w[finite])) if finite.any() else None,
                "max_finite": float(np.max(w[finite])) if finite.any() else None,
            }
        ix = np.floor(ra[candidate] / 4.0).astype(np.int32)
        iy = np.floor((np.sin(np.deg2rad(dec[candidate])) + 1) / 0.05).astype(np.int32)
        pixels = set((ix * 41 + np.clip(iy, 0, 40)).tolist())
        result = {
            "realization_id": rid, "tracer": tracer, "cap": cap,
            "path": path.name, "compressed_file_sha256": sha,
            "downloaded_compressed_bytes": nbytes,
            "header_rows": int(table.header["NAXIS2"]),
            "invalid_coordinate_z_rows": int(np.count_nonzero(~valid)),
            "candidate_0p6_to_1p0_rows": int(np.count_nonzero(candidate)),
            "candidate_retained_after_numerical_zero_selection": int(
                np.count_nonzero(retained_candidate)),
            "candidate_excluded_by_numerical_zero_selection": int(
                np.count_nonzero(candidate & ~retained_candidate)),
            "raw_zbin_counts": zcounts,
            "weight_column_diagnostics": wdiag,
            "candidate_coarse_data_occupied_cells": len(pixels),
        }
    return result, pixels


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/mock_selection_audit.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/mock_selection_catalogues")
    parser.add_argument("--ids", default="1,500,1000")
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--max-file-mib", type=int, default=16)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        from inspect_eboss_dr16_selection import self_test
        self_test()
        assert mock_path("eBOSS_ELG", "SGC", "dat", 500).endswith(
            "EZmock_realistic_eBOSS_ELG_SGC_v7_0500.dat.fits.gz")
        print("Mock data selection and path self-tests passed")
        return 0
    ids = sorted(set(int(s.strip()) for s in args.ids.split(",")))
    if (not ids or any(not 1 <= x <= 1000 for x in ids)
            or args.max_file_mib < 1 or args.timeout <= 0):
        parser.error("Invalid mock IDs, file size or timeout")
    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    real = {(v["tracer"], v["cap"]): v for v in reference["catalogues"]}
    records, errors, pixels = [], [], {}
    cache = Path(args.cache_dir)
    for rid in ids:
        for cap in CAPS:
            for tracer in TRACERS:
                kind = "LRG" if tracer == "eBOSS_LRG" else "ELG"
                relative = mock_path(tracer, cap, "dat", rid)
                local = cache / Path(relative).name
                try:
                    sha, length = fetch_with_retry(
                        BASE + relative, local, args.max_file_mib * 1024 * 1024,
                        args.timeout)
                    info, occupied = inspect(local, kind, cap, rid, sha, length)
                    baseline = real[(kind, cap)]
                    realcounts = baseline["zbin_counts"]
                    denom = sum(realcounts[:4])
                    mockdenom = sum(info["raw_zbin_counts"][:4])
                    info["observed_data_reference_sha256"] = baseline["sha256"]
                    info["real_data_candidate_rows"] = baseline["candidate_rows"]
                    info["mock_to_real_candidate_row_ratio"] = (
                        info["candidate_0p6_to_1p0_rows"] / baseline["candidate_rows"])
                    info["fractional_zbin_comparison"] = [
                        {"zlo": lo, "zhi": hi,
                         "real_data_fraction_in_candidate": (realcounts[i] / denom if denom else None),
                         "mock_fraction_in_candidate": (info["raw_zbin_counts"][i] / mockdenom
                                                        if mockdenom else None)}
                        for i, (lo, hi) in enumerate(zip(Z_EDGES[:-2], Z_EDGES[1:-1]))
                    ]
                    records.append(info)
                    pixels[(rid, cap, kind)] = occupied
                    print("MOCK_SELECTION_OK", rid, kind, cap,
                          "candidate", info["candidate_0p6_to_1p0_rows"],
                          "ratio", round(info["mock_to_real_candidate_row_ratio"], 4),
                          flush=True)
                except (OSError, ValueError, KeyError, RuntimeError) as exc:
                    errors.append(f"{relative}: {exc}")
                    print("MOCK_SELECTION_ERROR", errors[-1], flush=True)
                finally:
                    if local.exists():
                        local.unlink()
    joint = []
    for rid in ids:
        for cap in CAPS:
            a, b = pixels.get((rid, cap, "LRG")), pixels.get((rid, cap, "ELG"))
            if a is not None and b is not None:
                joint.append({
                    "realization_id": rid, "cap": cap,
                    "lrg_coarse_data_cells": len(a), "elg_coarse_data_cells": len(b),
                    "joint_coarse_data_cells": len(a & b),
                    "scope": "data-object occupancy only; not a random-derived joint mask",
                })
    complete = len(records) == len(ids) * len(CAPS) * len(TRACERS) and not errors
    health = all(
        rec["invalid_coordinate_z_rows"] == 0 and
        all(
            x["nonfinite_candidate_rows"] == 0 and (
                (x["near_zero_candidate_rows_abs_le_1e_20"] ==
                 x["near_zero_candidate_rows_abs_le_1e_12"] and
                 x["significant_negative_candidate_rows_below_minus_1e_12"] == 0)
                if name == "WEIGHT_SYSTOT"
                else (x["nonpositive_candidate_rows"] == 0 and
                      x["near_zero_candidate_rows_abs_le_1e_12"] == 0))
            for name, x in rec["weight_column_diagnostics"].items())
        for rec in records
    )
    report = {
        "study": "eBOSS DR16 realistic EZmock data-only selection sample",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "release_url": BASE, "realization_ids_checked": ids,
        "status": ("sample_candidate_selection_validated_with_numerical_zero_convention"
                   if complete and health else "candidate_weight_policy_unresolved"
                   if complete else "partial"),
        "candidate_z_interval": [CANDIDATE_LO, CANDIDATE_HI],
        "candidate_z_interval_frozen": False,
        "numerical_zero_weight_convention": {
            "column": "WEIGHT_SYSTOT",
            "exclude_if_weight_le": NUMERICAL_ZERO_WEIGHT_TOL,
            "same_selection_required_on_real_and_mock_data_and_randoms": True,
            "applied_to_published_real_and_mock_randoms_here": False,
            "rule_is_pre_odd_sector_and_not_a_template_fit": True,
        },
        "real_data_selection_reference": str(REFERENCE.relative_to(ROOT)),
        "mock_data_catalogues": records,
        "coarse_data_occupancy_by_realization_cap": joint,
        "errors": errors, "mock_random_files_downloaded": False,
        "mock_joint_selection_verified": False,
        "mock_covariance_estimated": False,
        "real_data_odd_vector_inspected": False,
        "note": "The three checked realization IDs sample released mock data only. "
                "Weights with absolute magnitude <= 1e-12 are separated from the "
                "other candidate WEIGHT_SYSTOT values by the recorded numerical gap. "
                "This candidate-only rule is a proposed input-sanitization convention; "
                "it must also be checked against data and matched mock randoms before "
                "the protocol is frozen. Counts, n(z) and coarse occupancy do not "
                "certify the joint mock selection or covariance.",
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("MOCK_SELECTION_AUDIT", target, report["status"], flush=True)
    return 0 if complete and health else 2


if __name__ == "__main__":
    raise SystemExit(main())
