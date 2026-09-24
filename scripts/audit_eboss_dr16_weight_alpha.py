#!/usr/bin/env python3
"""Pre-odd audit of DR16 LRG/ELG data and random weight normalization.

Reads *only* catalogue columns RA/DEC/Z and the four released weight columns.
Checks all four public data/random pairs and matched realistic mock 0001,
individually in NGC/SGC. Records selected weighted n(z), alpha=sum(Dw)/sum(Rw)
and Kish effective random counts. It does not count galaxy pairs, form an
odd data vector, estimate covariance or choose inferential binning.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from inspect_eboss_dr16_selection import (
    FILES as REAL_DATA, WEIGHTS, Z_EDGES, download,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS as REAL_RANDOM
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import (
    NUMERICAL_ZERO_WEIGHT_TOL, fetch_with_retry,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCES = {
    "observed_data": ROOT / "source_data/eboss_dr16_data_selection_audit_2026-09-24.json",
    "observed_random": ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json",
    "mock_data": ROOT / "source_data/eboss_dr16_realistic_mock_selection_audit_2026-09-24.json",
    "mock_random": ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json",
}
MOCK_ID = 1
CANDIDATE = (0.6, 1.0)
Z_WINDOWS = [CANDIDATE, *list(zip(Z_EDGES[:-2], Z_EDGES[1:-1]))]
LABELS = ["candidate_0p6_to_1p0", "0p6_to_0p7", "0p7_to_0p8",
          "0p8_to_0p9", "0p9_to_1p0"]


def summary_for_chunk(batch) -> tuple[np.ndarray, np.ndarray, dict]:
    ra = np.asarray(batch["RA"], dtype=np.float64)
    dec = np.asarray(batch["DEC"], dtype=np.float64)
    z = np.asarray(batch["Z"], dtype=np.float64)
    valid = (
        np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
        & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
    candidate = valid & (z >= CANDIDATE[0]) & (z < CANDIDATE[1])
    kept = candidate.copy()
    weight = np.ones(len(z), dtype="f8")
    diag = {}
    for name in WEIGHTS:
        w = np.asarray(batch[name], dtype="f8")
        finite = np.isfinite(w)
        nonfinite = int(np.count_nonzero(candidate & ~finite))
        sig_negative = int(np.count_nonzero(
            candidate & finite & (w < -NUMERICAL_ZERO_WEIGHT_TOL)))
        nonpositive = int(np.count_nonzero(candidate & finite & (w <= 0)))
        zero = int(np.count_nonzero(
            candidate & finite & (np.abs(w) <= NUMERICAL_ZERO_WEIGHT_TOL)))
        tiny = int(np.count_nonzero(
            candidate & finite & (np.abs(w) <= 1e-20)))
        unity = int(np.count_nonzero(candidate & finite & (w == 1)))
        diag[name] = {
            "candidate_nonfinite": nonfinite,
            "candidate_significant_negative": sig_negative,
            "candidate_nonpositive": nonpositive,
            "candidate_numerical_zero_le_1e_12": zero,
            "candidate_tiny_le_1e_20": tiny,
            "candidate_exact_unity": unity,
        }
        floor = NUMERICAL_ZERO_WEIGHT_TOL if name == "WEIGHT_SYSTOT" else 0.0
        kept &= finite & (w > floor)
        weight *= np.where(finite, w, 0.0)
    return z, kept, {
        "candidate_rows": int(np.count_nonzero(candidate)),
        "invalid_coordinate_z_rows": int(np.count_nonzero(~valid)),
        "weights": diag,
        "weighted_product": weight,
    }


def initial_accum() -> dict:
    return {
        "candidate_rows": 0, "candidate_retained_rows": 0,
        "invalid_coordinate_z_rows": 0,
        "weight_columns": {
            name: {
                "candidate_nonfinite": 0,
                "candidate_significant_negative": 0,
                "candidate_nonpositive": 0,
                "candidate_numerical_zero_le_1e_12": 0,
                "candidate_tiny_le_1e_20": 0,
                "candidate_exact_unity": 0,
            }
            for name in WEIGHTS
        },
        "windows": {
            label: {
                "rows": 0, "sum_weight": 0., "sum_weight_squared": 0.,
            }
            for label in LABELS
        },
    }


def inspect(path: Path, chunk_rows: int) -> dict:
    acc = initial_accum()
    with fits.open(path, memmap=path.suffix != ".gz") as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError(f"Expected one BINTABLE: {path.name}")
        table = tables[0]
        required = {"RA", "DEC", "Z", *WEIGHTS}
        if not required.issubset(set(table.columns.names)):
            raise ValueError(f"Missing required weight/coordinate columns in {path.name}")
        nrows = int(table.header["NAXIS2"])
        for first in range(0, nrows, chunk_rows):
            batch = table.data[first:min(nrows, first + chunk_rows)]
            z, kept, part = summary_for_chunk(batch)
            acc["candidate_rows"] += part["candidate_rows"]
            acc["candidate_retained_rows"] += int(np.count_nonzero(kept))
            acc["invalid_coordinate_z_rows"] += part["invalid_coordinate_z_rows"]
            for name, diagnostic in part["weights"].items():
                for metric, count in diagnostic.items():
                    acc["weight_columns"][name][metric] += count
            w = part["weighted_product"]
            for label, (lo, hi) in zip(LABELS, Z_WINDOWS):
                use = kept & (z >= lo) & (z < hi)
                if not use.any():
                    continue
                selected = w[use]
                item = acc["windows"][label]
                item["rows"] += int(np.count_nonzero(use))
                item["sum_weight"] += float(np.sum(selected, dtype="f8"))
                item["sum_weight_squared"] += float(np.dot(selected, selected))
    acc["header_rows"] = nrows
    for item in acc["windows"].values():
        item["Kish_effective_weighted_rows"] = (
            item["sum_weight"] ** 2 / item["sum_weight_squared"]
            if item["sum_weight_squared"] > 0 else None
        )
    return acc


def health(acc: dict, expected_candidate: int, expected_retained: int) -> list[str]:
    failures = []
    if acc["candidate_rows"] != expected_candidate:
        failures.append("Candidate row count disagrees with retained source audit")
    if acc["candidate_retained_rows"] != expected_retained:
        failures.append("Retained input row count disagrees with retained source audit")
    if acc["invalid_coordinate_z_rows"]:
        failures.append("Invalid coordinate/redshift")
    for name, info in acc["weight_columns"].items():
        if info["candidate_nonfinite"] or info["candidate_significant_negative"]:
            failures.append(f"Nonfinite or significant negative candidate {name}")
        if name == "WEIGHT_SYSTOT":
            if (info["candidate_numerical_zero_le_1e_12"]
                    != info["candidate_tiny_le_1e_20"]):
                failures.append("Numerical-zero gap fails at WEIGHT_SYSTOT")
        elif info["candidate_nonpositive"]:
            failures.append(f"Nonpositive {name}")
    if not math.isclose(
        acc["windows"][LABELS[0]]["sum_weight"],
        sum(acc["windows"][name]["sum_weight"] for name in LABELS[1:]),
        rel_tol=1e-11, abs_tol=1e-9,
    ):
        failures.append("Selected weight sum not conserved across redshift bins")
    if (sum(acc["windows"][name]["rows"] for name in LABELS[1:])
            != acc["candidate_retained_rows"]):
        failures.append("Retained row counts not conserved across redshift bins")
    if any(not np.isfinite(item["sum_weight"]) or item["sum_weight"] < 0
           or not np.isfinite(item["sum_weight_squared"])
           for item in acc["windows"].values()):
        failures.append("Invalid weighted sum")
    return failures


def combinations(data: dict, random: dict) -> dict:
    result = {}
    candidate_data_total = data["windows"][LABELS[0]]["sum_weight"]
    candidate_random_total = random["windows"][LABELS[0]]["sum_weight"]
    for name in LABELS:
        d, r = data["windows"][name], random["windows"][name]
        if min(d["sum_weight"], r["sum_weight"]) <= 0:
            raise ValueError(f"Nonpositive data or random weighted sum for {name}")
        result[name] = {
            "data_retained_rows": d["rows"],
            "random_retained_rows": r["rows"],
            "sum_data_weight": d["sum_weight"],
            "sum_random_weight": r["sum_weight"],
            "alpha_sum_data_weights_over_random_weights":
                d["sum_weight"] / r["sum_weight"],
            "data_Kish_effective_weighted_rows": d["Kish_effective_weighted_rows"],
            "random_Kish_effective_weighted_rows": r["Kish_effective_weighted_rows"],
            "fraction_of_candidate_data_weight": (
                d["sum_weight"] / candidate_data_total),
            "fraction_of_candidate_random_weight": (
                r["sum_weight"] / candidate_random_total),
        }
    return result


def self_test():
    array = np.zeros(5, dtype=[(n, "f8") for n in ("RA", "DEC", "Z", *WEIGHTS)])
    for w in WEIGHTS:
        array[w] = 1
    array["RA"][:] = 120.
    array["Z"][:] = [0.6, 0.7, 0.8, 0.9, 1.0]
    array["WEIGHT_SYSTOT"][:] = [1, 1e-30, 1.5, -1e-11, 1]
    z, good, partial = summary_for_chunk(array)
    assert partial["candidate_rows"] == 4
    assert good.tolist() == [True, False, True, False, False]
    assert partial["weights"]["WEIGHT_SYSTOT"]["candidate_significant_negative"] == 1
    assert partial["weights"]["WEIGHT_SYSTOT"]["candidate_numerical_zero_le_1e_12"] == 1
    assert np.array_equal(z, array["Z"])
    print("Full-catalogue weight and alpha audit self-tests passed")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="eboss_workspace/weight_alpha_audit.json")
    p.add_argument("--cache-dir", default="eboss_workspace/weight_alpha_catalogues")
    p.add_argument("--chunk-rows", type=int, default=150000)
    p.add_argument("--timeout", type=float, default=120)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.chunk_rows < 1000 or args.timeout <= 0:
        p.error("Invalid chunk size or timeout")
    refs = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in REFERENCES.items()
    }
    if (
        refs["observed_data"]["status"] != "four_data_catalogues_checked"
        or refs["observed_random"]["status"] != "four_random_catalogues_checked"
        or refs["mock_data"]["status"]
        != "sample_candidate_selection_validated_with_numerical_zero_convention"
        or refs["mock_random"]["status"] != "mock_random_sample_selection_compatible"
    ):
        raise RuntimeError("Input sample selection audits are not complete")
    real_data = {
        (a["tracer"], a["cap"]): a for a in refs["observed_data"]["catalogues"]
    }
    real_random = {
        (a["tracer"], a["cap"]): a for a in refs["observed_random"]["randoms"]
    }
    mock_random = {
        (int(a["id"]), a["tracer"], a["cap"]): a
        for a in refs["mock_random"]["mock_random_catalogues"]
    }
    records, combined, errors = [], [], []
    cache = Path(args.cache_dir)
    for catalogue_type in ("observed", "realistic_mock"):
        for cap in ("NGC", "SGC"):
            selected = {}
            for tracer in ("LRG", "ELG"):
                for role in ("data", "random"):
                    real = catalogue_type == "observed"
                    if real:
                        files = REAL_DATA if role == "data" else REAL_RANDOM
                        filename, expected_rows = files[(tracer, cap)]
                        expected_info = (
                            real_data[(tracer, cap)] if role == "data"
                            else real_random[(tracer, cap)]
                        )
                        expected_sha = expected_info["sha256"]
                        expected_candidate = (
                            expected_info["candidate_rows"]
                            if role == "data" else expected_info["candidate_rows"])
                        expected_retained = expected_candidate
                        local = cache / catalogue_type / filename
                    else:
                        kind = "dat" if role == "data" else "ran"
                        relative = mock_path("eBOSS_" + tracer, cap, kind, MOCK_ID)
                        local = cache / catalogue_type / Path(relative).name
                        key = f"{tracer}_{cap}"
                        if role == "data":
                            counts = refs["mock_data"][
                                "candidate_rows_raw_retained_excluded"][f"{MOCK_ID:04d}"][key]
                            expected_candidate, expected_retained = counts[:2]
                            expected_sha = None
                        else:
                            expected_info = mock_random[(MOCK_ID, tracer, cap)]
                            expected_sha = expected_info["compressed_file_sha256"]
                            expected_candidate = expected_info["candidate_random_rows"]
                            expected_retained = (
                                expected_candidate -
                                expected_info["candidate_zero_weight_exclusions"])
                        expected_rows = None
                    try:
                        if real:
                            path, sha, size = download(
                                filename, local.parent,
                                1024 * 1024 * 1024, args.timeout)
                        else:
                            sha, size = fetch_with_retry(
                                MOCK_BASE + relative, local,
                                512 * 1024 * 1024, args.timeout)
                            path = local
                        if expected_sha and sha != expected_sha:
                            raise ValueError("Catalog SHA256 differs from earlier full-file audit")
                        result = inspect(path, args.chunk_rows)
                        if expected_rows is not None and result["header_rows"] != expected_rows:
                            raise ValueError("FITS row count differs from published inventory")
                        bad = health(result, expected_candidate, expected_retained)
                        if bad:
                            print("WEIGHT_ALPHA_HEALTH_DETAILS",
                                  catalogue_type, cap, tracer, role,
                                  "candidate", result["candidate_rows"],
                                  "retained", result["candidate_retained_rows"],
                                  "expected", expected_retained,
                                  "SYSTOT", json.dumps(
                                      result["weight_columns"]["WEIGHT_SYSTOT"],
                                      sort_keys=True),
                                  "issues", bad, flush=True)
                            raise ValueError("; ".join(bad))
                        result.update({
                            "catalogue": catalogue_type, "tracer": tracer, "cap": cap,
                            "role": role, "filename": path.name,
                            "file_sha256": sha, "downloaded_bytes": size,
                        })
                        selected[(tracer, role)] = result
                        records.append(result)
                        print("WEIGHT_ALPHA_INPUT_OK", catalogue_type, cap, tracer,
                              role, "candidate", result["candidate_rows"],
                              "retained", result["candidate_retained_rows"],
                              "sha256", sha, flush=True)
                    except (OSError, ValueError, RuntimeError,
                            KeyError, MemoryError) as exc:
                        errors.append(f"{catalogue_type}/{cap}/{tracer}/{role}: {exc}")
                        print("WEIGHT_ALPHA_ERROR", errors[-1], flush=True)
                    finally:
                        if local.exists():
                            local.unlink()
            for tracer in ("LRG", "ELG"):
                if (tracer, "data") in selected and (tracer, "random") in selected:
                    try:
                        weight_alpha = combinations(
                            selected[(tracer, "data")], selected[(tracer, "random")])
                        combined.append({
                            "catalogue": catalogue_type, "cap": cap,
                            "tracer": tracer, "window_weight_normalizations": weight_alpha,
                        })
                        overall = weight_alpha[LABELS[0]]
                        print("WEIGHT_ALPHA", catalogue_type, cap, tracer,
                              "alpha", overall[
                                  "alpha_sum_data_weights_over_random_weights"],
                              "random_Neff", overall[
                                  "random_Kish_effective_weighted_rows"], flush=True)
                    except (ValueError, KeyError) as exc:
                        errors.append(f"{catalogue_type}/{cap}/{tracer}/alpha: {exc}")
    successful = len(records) == 16 and len(combined) == 8 and not errors
    report = {
        "study": "pre-odd full-catalogue DR16 data/random weight and alpha audit",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "observed_and_mock_sample_weight_alpha_checked"
            if successful else "partial",
        "source_audits": {
            k: str(v.relative_to(ROOT)) for k, v in REFERENCES.items()},
        "mock_id": MOCK_ID,
        "candidate_z_interval": list(CANDIDATE),
        "candidate_z_interval_frozen": False,
        "subwindows": [
            {"label": label, "zlo": lo, "zhi": hi}
            for label, (lo, hi) in zip(LABELS, Z_WINDOWS)
        ],
        "provisional_column_product": "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP",
        "numerical_zero_tolerance_on_mock_systot": NUMERICAL_ZERO_WEIGHT_TOL,
        "observed_galaxy_coordinates_read_for_input_metadata_only": True,
        "observed_galaxy_pair_counts_computed": False,
        "observed_odd_data_vector_read": False,
        "catalogue_inspections": records,
        "data_random_weight_normalizations": combined,
        "errors": errors,
        "physical_random_weight_normalization_frozen": False,
        "mock_ensemble_all_1000_validated": False,
        "RR_window_computed": False,
        "note": (
            "The alpha and effective random counts are input-only diagnostics. "
            "The per-tracer, per-cap, per-bin alpha does not by itself validate "
            "the joint angular mask, RR pair geometry or even-to-odd leakage. "
            "Production weights must be checked against eBOSS data models; "
            "the candidate redshift interval is not frozen."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("WEIGHT_ALPHA_AUDIT", report["status"], output, flush=True)
    return 0 if successful else 2


if __name__ == "__main__":
    raise SystemExit(main())
