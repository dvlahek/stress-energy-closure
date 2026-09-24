#!/usr/bin/env python3
"""Pre-odd eBOSS LRG/ELG weighted data/random normalization audit.

The release and fixed-realization files are matched by audited SHA256 or
previously audited row counts. All catalogue rows are read solely to inspect
redshift, total input weight, weighted n(z), and ELG chunk normalization.
No data-data, data-random or random-random pairs, galaxy odd multipoles,
wake amplitudes or mock covariances are computed.

The candidate interval 0.6 <= z < 1.0 and bins of width 0.1 are diagnostics,
not a frozen selection. The numerical-zero systematic-weight convention
is identical for data and randoms, observed and mock.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from eboss_dr16_fiducial import (
    DISTANCE_CONVENTIONS, PRIMARY_GEOMETRY, SECONDARY_GEOMETRY,
    SYSTOT_NUMERICAL_ZERO_TOL, WEIGHT_COLUMNS, comoving_mpc_over_h,
    validated_weight_product,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry
from inspect_eboss_dr16_selection import FILES as REAL_DATA, Z_EDGES, download

ROOT = Path(__file__).resolve().parents[1]
REAL_DATA_REF = ROOT / "source_data/eboss_dr16_data_selection_audit_2026-09-24.json"
REAL_RANDOM_REF = (
    ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json"
)
MOCK_DATA_REF = (
    ROOT / "source_data/eboss_dr16_realistic_mock_selection_audit_2026-09-24.json"
)
MOCK_RANDOM_REF = (
    ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
)
BINS = tuple((float(lo), float(hi)) for lo, hi in
             zip(Z_EDGES[:4], Z_EDGES[1:5]))
MOCK_ID = 1


def empty_bucket() -> dict:
    return {
        "rows": 0, "retained_rows": 0, "numerical_zero_excluded": 0,
        "sum_weights": 0.0, "sum_squared_weights": 0.0,
        "sum_completeness_weight": 0.0, "sum_fkp_weight": 0.0,
    }


def accumulate_bucket(bucket: dict, total: np.ndarray,
                      completeness: np.ndarray, fkp: np.ndarray,
                      retained: np.ndarray) -> None:
    bucket["rows"] += int(len(total))
    bucket["retained_rows"] += int(np.count_nonzero(retained))
    bucket["numerical_zero_excluded"] += int(np.count_nonzero(~retained))
    bucket["sum_weights"] += float(np.sum(total, dtype="f8"))
    bucket["sum_squared_weights"] += float(np.dot(total, total))
    bucket["sum_completeness_weight"] += float(
        np.sum(completeness, dtype="f8"))
    bucket["sum_fkp_weight"] += float(np.sum(
        np.where(retained, fkp, 0.0), dtype="f8"))


def report_bucket(bucket: dict) -> dict:
    result = dict(bucket)
    sumw2 = bucket["sum_squared_weights"]
    result["kish_effective_random_count"] = (
        bucket["sum_weights"] ** 2 / sumw2 if sumw2 > 0 else None)
    return result


def inspect(path: Path, tracer: str, cap: str, role: str, survey: str,
            expected_full_rows: int, chunk_rows: int) -> dict:
    bins = [empty_bucket() for _ in BINS]
    chunks = defaultdict(empty_bucket)
    excluded_by_bin = [0 for _ in BINS]
    invalid_candidate_coords = 0
    finite_weight_error = 0
    declared = None
    chunk_column_available = False
    elg_random_cp_noz_not_one = 0
    sys_abs_le_1e_20 = 0
    sys_abs_le_1e_12 = 0
    sys_max_abs_below_tolerance = 0.0
    sys_min_above_tolerance = None
    with fits.open(path, memmap=(path.suffix != ".gz")) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected exactly one FITS binary table")
        table = tables[0]
        declared = int(table.header["NAXIS2"])
        if declared != expected_full_rows:
            raise ValueError(
                f"Published FITS row count changed {declared} != {expected_full_rows}")
        missing = {"RA", "DEC", "Z", *WEIGHT_COLUMNS} - set(table.columns.names)
        if missing:
            raise ValueError(f"FITS missing required columns: {sorted(missing)}")
        chunk_column_available = "chunk" in table.columns.names
        if tracer == "ELG" and not chunk_column_available:
            raise ValueError("ELG catalogue has no chunk labels")
        for first in range(0, declared, chunk_rows):
            data = table.data[first:min(first + chunk_rows, declared)]
            ra = np.asarray(data["RA"], dtype="f8")
            dec = np.asarray(data["DEC"], dtype="f8")
            z = np.asarray(data["Z"], dtype="f8")
            valid = (
                np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
            target_z = np.isfinite(z) & (z >= BINS[0][0]) & (z < BINS[-1][1])
            invalid_candidate_coords += int(np.count_nonzero(target_z & ~valid))
            for i, (lo, hi) in enumerate(BINS):
                mask = valid & (z >= lo) & (z < hi)
                if not np.any(mask):
                    continue
                w = {
                    name: np.asarray(data[name][mask], dtype="f8")
                    for name in WEIGHT_COLUMNS
                }
                # Audit numerical-zero separation independently of the gate.
                sy = np.abs(w["WEIGHT_SYSTOT"])
                sys_abs_le_1e_20 += int(np.count_nonzero(sy <= 1e-20))
                near_zero = sy <= SYSTOT_NUMERICAL_ZERO_TOL
                sys_abs_le_1e_12 += int(np.count_nonzero(near_zero))
                if np.any(near_zero):
                    sys_max_abs_below_tolerance = max(
                        sys_max_abs_below_tolerance,
                        float(np.max(sy[near_zero])))
                above = sy > SYSTOT_NUMERICAL_ZERO_TOL
                if np.any(above):
                    minimum = float(np.min(sy[above]))
                    sys_min_above_tolerance = (
                        minimum if sys_min_above_tolerance is None
                        else min(sys_min_above_tolerance, minimum))
                # The helper must fail rather than silently remove a
                # significant negative or nonfinite candidate weight.
                total, keep = validated_weight_product(w)
                completeness = np.where(
                    keep,
                    w["WEIGHT_SYSTOT"] * w["WEIGHT_CP"] * w["WEIGHT_NOZ"],
                    0.0,
                )
                accumulate_bucket(
                    bins[i], total, completeness, w["WEIGHT_FKP"], keep)
                excluded_by_bin[i] += int(np.count_nonzero(~keep))
                if tracer == "ELG":
                    chunk_labels = np.asarray(
                        data["chunk"][mask]).astype(str)
                    for chunk_name in np.unique(chunk_labels):
                        k = chunk_labels == chunk_name
                        accumulate_bucket(
                            chunks[str(chunk_name)], total[k], completeness[k],
                            w["WEIGHT_FKP"][k], keep[k])
                    if role == "random":
                        elg_random_cp_noz_not_one += int(np.count_nonzero(
                            (np.abs(w["WEIGHT_CP"] - 1.0) > 1e-10)
                            | (np.abs(w["WEIGHT_NOZ"] - 1.0) > 1e-10)
                        ))
    if invalid_candidate_coords:
        raise ValueError("Invalid RA/DEC within candidate redshift interval")
    if sys_abs_le_1e_20 != sys_abs_le_1e_12:
        raise ValueError(
            "Numerical-zero WEIGHT_SYSTOT entries have no verified magnitude gap")
    counts = [report_bucket(v) for v in bins]
    return {
        "survey": survey, "tracer": tracer, "cap": cap,
        "role": role, "filename": path.name, "full_fits_header_rows": declared,
        "redshift_bins": [
            {"zlo": lo, "zhi": hi, **counts[i]}
            for i, (lo, hi) in enumerate(BINS)
        ],
        "candidate_rows": sum(v["rows"] for v in counts),
        "candidate_retained_rows": sum(v["retained_rows"] for v in counts),
        "candidate_numerical_zero_exclusions": sum(
            v["numerical_zero_excluded"] for v in counts),
        "candidate_systot_abs_le_1e_20": sys_abs_le_1e_20,
        "candidate_systot_abs_le_1e_12": sys_abs_le_1e_12,
        "candidate_systot_max_abs_below_tolerance": (
            sys_max_abs_below_tolerance if sys_abs_le_1e_12 else None),
        "candidate_systot_min_abs_above_tolerance": sys_min_above_tolerance,
        "candidate_total_weight": sum(v["sum_weights"] for v in counts),
        "candidate_total_completeness_weight": sum(
            v["sum_completeness_weight"] for v in counts),
        "chunk_column_present": chunk_column_available,
        "elg_chunks": {
            key: report_bucket(value) for key, value in sorted(chunks.items())
        },
        "elg_random_cp_or_noz_nonunity_candidate_rows": (
            elg_random_cp_noz_not_one if tracer == "ELG" and role == "random"
            else None),
        "odd_vector_read": False,
    }


def alpha(data: dict, random: dict) -> dict:
    db, rb = data["redshift_bins"], random["redshift_bins"]
    comparisons = []
    for a, b in zip(db, rb):
        if (a["zlo"], a["zhi"]) != (b["zlo"], b["zhi"]):
            raise ValueError("Data and random redshift bins differ")
        if a["sum_weights"] <= 0 or b["sum_weights"] <= 0:
            raise ValueError("Data or random candidate bin has no weight")
        comparisons.append({
            "zlo": a["zlo"], "zhi": a["zhi"],
            "weighted_data_to_random_alpha": a["sum_weights"] / b["sum_weights"],
            "completeness_only_alpha": (
                a["sum_completeness_weight"] /
                b["sum_completeness_weight"]
                if b["sum_completeness_weight"] else None),
            "data_retained_rows": a["retained_rows"],
            "random_retained_rows": b["retained_rows"],
        })
    if not data["candidate_total_weight"] or not random["candidate_total_weight"]:
        raise ValueError("An aggregate candidate weight is zero")
    chunks = {}
    if data["tracer"] == "ELG":
        for k in sorted(set(data["elg_chunks"]) | set(random["elg_chunks"])):
            if k not in data["elg_chunks"] or k not in random["elg_chunks"]:
                raise ValueError(f"Unmatched ELG chunk in data/random: {k}")
            da, ra = data["elg_chunks"][k], random["elg_chunks"][k]
            if da["sum_weights"] <= 0 or ra["sum_weights"] <= 0:
                raise ValueError(f"Empty weighted ELG chunk: {k}")
            chunks[k] = {
                "data_rows": da["retained_rows"],
                "random_rows": ra["retained_rows"],
                "weighted_alpha": da["sum_weights"] / ra["sum_weights"],
                "completeness_only_alpha": (
                    da["sum_completeness_weight"] /
                    ra["sum_completeness_weight"]
                    if ra["sum_completeness_weight"] else None),
            }
    return {
        "tracer": data["tracer"], "cap": data["cap"], "survey": data["survey"],
        "candidate_weighted_alpha": (
            data["candidate_total_weight"] / random["candidate_total_weight"]),
        "per_bin": comparisons, "per_elg_chunk": chunks,
    }


def self_test() -> None:
    x = empty_bucket()
    total = np.array([2.0, 0.0, 3.0])
    keep = np.array([True, False, True])
    accumulate_bucket(
        x, total, np.array([1.0, 0.0, 1.5]),
        np.array([2.0, 1.0, 2.0]), keep)
    assert x["rows"] == 3 and x["retained_rows"] == 2
    assert x["numerical_zero_excluded"] == 1
    assert x["sum_weights"] == 5.0 and x["sum_squared_weights"] == 13.0
    assert abs(report_bucket(x)["kish_effective_random_count"] - 25 / 13) < 1e-12
    primary = comoving_mpc_over_h([0.6, 0.8, 1.0], PRIMARY_GEOMETRY)
    secondary = comoving_mpc_over_h([0.6, 0.8, 1.0], SECONDARY_GEOMETRY)
    assert np.all(np.diff(primary) > 0) and np.all(np.diff(secondary) > 0)
    print("EBOSS_WEIGHTED_NORMALIZATION_SELF_TEST_OK")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/eboss_weighted_normalization.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/eboss_normalization_fits")
    parser.add_argument("--chunk-rows", type=int, default=150000)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.chunk_rows < 1000 or args.timeout <= 0:
        parser.error("Invalid chunk size or timeout")
    refs = {
        "observed_data": json.loads(REAL_DATA_REF.read_text()),
        "observed_random": json.loads(REAL_RANDOM_REF.read_text()),
        "mock_data": json.loads(MOCK_DATA_REF.read_text()),
        "mock_random": json.loads(MOCK_RANDOM_REF.read_text()),
    }
    if (refs["observed_data"]["status"] != "four_data_catalogues_checked"
            or refs["observed_random"]["status"] != "four_random_catalogues_checked"
            or refs["mock_random"]["status"] != "mock_random_sample_selection_compatible"
            or refs["mock_data"]["status"] !=
            "sample_candidate_selection_validated_with_numerical_zero_convention"):
        raise ValueError("Source catalogue selection gates are incomplete")
    real_d = {(v["tracer"], v["cap"]): v
              for v in refs["observed_data"]["catalogues"]}
    real_r = {(v["tracer"], v["cap"]): v
              for v in refs["observed_random"]["randoms"]}
    mock_r = {
        (int(v["id"]), v["tracer"], v["cap"]): v
        for v in refs["mock_random"]["mock_random_catalogues"]
    }
    mock_d = refs["mock_data"]["candidate_rows_raw_retained_excluded"]["0001"]
    root = Path(args.cache_dir)
    records, errors = [], []
    for survey in ("observed", "realistic_mock"):
        for cap in ("NGC", "SGC"):
            for tracer in ("LRG", "ELG"):
                for role in ("data", "random"):
                    actual = survey == "observed"
                    if actual:
                        filename, declared = (REAL_DATA if role == "data"
                                              else RANDOMS)[(tracer, cap)]
                        path = root / survey / filename
                        evidence = (real_d if role == "data" else real_r)[(tracer, cap)]
                        expected_sha = evidence["sha256"]
                    else:
                        relative = mock_path(
                            "eBOSS_" + tracer, cap,
                            "dat" if role == "data" else "ran", MOCK_ID)
                        path = root / survey / Path(relative).name
                        if role == "random":
                            evidence = mock_r[(MOCK_ID, tracer, cap)]
                            declared = evidence["header_rows"]
                            expected_sha = evidence["compressed_file_sha256"]
                        else:
                            evidence = mock_d[f"{tracer}_{cap}"]
                            declared = None
                            expected_sha = None
                    try:
                        if actual:
                            path, checksum, nbytes = download(
                                filename, path.parent, 1024 * 1024 * 1024,
                                args.timeout)
                        else:
                            checksum, nbytes = fetch_with_retry(
                                MOCK_BASE + relative, path,
                                512 * 1024 * 1024, args.timeout)
                        if expected_sha and checksum != expected_sha:
                            raise ValueError("Catalogue SHA256 differs from verified snapshot")
                        if declared is None:
                            # The compact mock-data snapshot gives candidate
                            # counts; its original workflow artifact holds
                            # the original full-file SHA256 digests.
                            with fits.open(path, memmap=False) as hdr:
                                declared = int(hdr[1].header["NAXIS2"])
                        record = inspect(
                            path, tracer, cap, role, survey, declared,
                            args.chunk_rows)
                        if actual:
                            if record["candidate_rows"] != evidence["candidate_rows"]:
                                raise ValueError("Observed candidate count changed")
                        elif role == "random":
                            if (record["candidate_rows"] != evidence["candidate_random_rows"]
                                    or record["candidate_numerical_zero_exclusions"] !=
                                    evidence["candidate_zero_weight_exclusions"]):
                                raise ValueError("Mock random candidate or exclusions changed")
                        else:
                            if (record["candidate_rows"] != evidence[0]
                                    or record["candidate_retained_rows"] != evidence[1]
                                    or record["candidate_numerical_zero_exclusions"] != evidence[2]):
                                raise ValueError("Mock data candidate or exclusions changed")
                        record["full_file_sha256"] = checksum
                        record["downloaded_bytes"] = nbytes
                        records.append(record)
                        print("EBOSS_NORMALIZATION_INPUT_OK", survey, cap,
                              tracer, role,
                              "candidate", record["candidate_rows"],
                              "excluded", record["candidate_numerical_zero_exclusions"],
                              "sumw", f"{record['candidate_total_weight']:.7g}",
                              "sha256", checksum, flush=True)
                    except (OSError, ValueError, RuntimeError, KeyError,
                            MemoryError) as exc:
                        errors.append(f"{survey}/{cap}/{tracer}/{role}: {exc}")
                        print("EBOSS_NORMALIZATION_INPUT_ERROR", errors[-1],
                              flush=True)
                    finally:
                        if path.exists():
                            path.unlink()
    indexed = {(v["survey"], v["cap"], v["tracer"], v["role"]): v
               for v in records}
    comparisons = []
    for survey in ("observed", "realistic_mock"):
        for cap in ("NGC", "SGC"):
            for tracer in ("LRG", "ELG"):
                dk, rk = ((survey, cap, tracer, role)
                          for role in ("data", "random"))
                if dk in indexed and rk in indexed:
                    try:
                        comp = alpha(indexed[dk], indexed[rk])
                        comparisons.append(comp)
                        print("EBOSS_WEIGHTED_ALPHA", survey, cap, tracer,
                              f"{comp['candidate_weighted_alpha']:.10g}",
                              "chunks", len(comp["per_elg_chunk"]), flush=True)
                    except (OSError, ValueError, KeyError) as exc:
                        errors.append(f"{survey}/{cap}/{tracer}: {exc}")
    ok = len(records) == 16 and len(comparisons) == 8 and not errors
    ztest = np.linspace(0.6, 1.0, 9)
    a = comoving_mpc_over_h(ztest, PRIMARY_GEOMETRY)
    b = comoving_mpc_over_h(ztest, SECONDARY_GEOMETRY)
    report = {
        "study": "eBOSS data-and-random weighted input normalization pre-odd audit",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "sample_weighted_normalization_audited" if ok else "partial",
        "fixed_mock_realization_id": MOCK_ID,
        "candidate_z_edges": list(Z_EDGES[:5]),
        "candidate_selection_frozen": False,
        "weight_columns_multiplied": list(WEIGHT_COLUMNS),
        "numerical_zero_systot_abs_tolerance": SYSTOT_NUMERICAL_ZERO_TOL,
        "distance_conventions": DISTANCE_CONVENTIONS,
        "geometry_primary_for_random_only_tests": PRIMARY_GEOMETRY,
        "geometry_secondary_predeclared": SECONDARY_GEOMETRY,
        "distance_mapping_assumptions": (
            "flat matter+Lambda only; radiation and neutrino background "
            "neglected with Tcmb0=0; same distance choice on real and mocks"),
        "distance_test_z": ztest.tolist(),
        "distance_primary_mpc_over_h": a.tolist(),
        "distance_secondary_mpc_over_h": b.tolist(),
        "max_fractional_distance_shift_secondary_vs_primary":
            float(np.max(np.abs(b / a - 1))),
        "records": records,
        "normalizations": comparisons,
        "errors": errors,
        "mock_data_sha256_only_newly_recorded": True,
        "observed_galaxy_positions_used_only_for_input_weights": True,
        "observed_galaxy_pairs_computed": False,
        "observed_odd_data_vector_read": False,
        "rr_pair_counts_computed": False,
        "wake_template_tuned": False,
        "cross_covariance_computed": False,
        "full_mock_ensemble_checked": False,
        "note": (
            "Weighted alpha and n(z) are input-selection diagnostics, not "
            "a validated estimator window. The ELG chunk entries permit "
            "separate assessment of normalization. The exact analysis "
            "selection, weights, cosmology, window and covariance are not frozen."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("EBOSS_WEIGHTED_NORMALIZATION", report["status"], output, flush=True)
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
