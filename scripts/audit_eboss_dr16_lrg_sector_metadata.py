#!/usr/bin/env python3
"""SHA-gated blinded eBOSS DR16 LRG sector-metadata consistency audit.

Only four FITS row columns may be decoded: SECTOR, sector_TSR,
sector_SSR and COMP_BOSS. All catalogue redshifts, coordinates,
weights, target IDs, mock/random rows, masks, pairs and observed
odd data are inaccessible to this code path. Output aggregates only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import struct

import numpy as np

from inspect_eboss_dr16_fits_headers import hdu_data_length, parse_header

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_sector_metadata_protocol_2026-09-25.json"
PARENT = ROOT / "source_data/eboss_dr16_lrg_full_bytes_provenance_protocol_2026-09-25.json"
ARCHIVE = ROOT / "source_data/eboss_dr16_lrg_full_bytes_provenance_local_2026-09-25.json"
FITS_WIDTH = {
    "A": 1, "L": 1, "B": 1, "I": 2, "J": 4, "K": 8,
    "E": 4, "D": 8, "C": 8, "M": 16, "P": 8, "Q": 16,
}
TFORM = re.compile(r"^([1-9][0-9]*)?([A-Z])$")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def report_sha_and_json(path):
    data = path.read_bytes()
    return hashlib.sha256(data).hexdigest(), json.loads(data)


def frozen_gate(p):
    parent = read_json(PARENT)
    raw_sha, raw = report_sha_and_json(ROOT / p["parent_lrg_full_bytes_local_report"])
    header_sha, header = report_sha_and_json(
        ROOT / p["parent_lrg_header_local_report"]
    )
    archived = ARCHIVE.read_bytes()
    if (
        raw_sha != p["parent_lrg_full_bytes_report_sha256"]
        or hashlib.sha256(archived).hexdigest() != raw_sha
        or (ROOT / p["parent_lrg_full_bytes_local_report"]).read_bytes() != archived
        or header_sha != p["parent_lrg_header_report_sha256"]
        or p["first_seen_full_file_sha256_frozen_before_sector_rows"] !=
           "39b831801adec04fe6dc5d6ab76a4b303aa58bfd303548cb7fddfae7d7f9331d"
        or p["first_header_sha256"] !=
           "a48b0d1055b5f7638dc77a86eb5ff5c861b629d19ff75cd8ae8122caa7e9014b"
        or p["official_source_url"] != parent["official_source_url"]
        or p["local_source_file"] != parent["local_output"]
        or p["official_source_file_size_bytes"] != parent["expected_file_bytes"]
        or p["first_header_bytes"] != parent["expected_header_bytes"]
        or p["first_header_sha256"] != parent["expected_header_sha256"]
        or raw.get("status") != "OFFICIAL_LRG_FULL_BYTES_SHA_PINNED_ONLY"
        or raw.get("full_file_sha256") !=
           p["first_seen_full_file_sha256_frozen_before_sector_rows"]
        or raw.get("file_bytes") != p["official_source_file_size_bytes"]
        or raw.get("verified_first_header_sha256") != p["first_header_sha256"]
        or raw.get("verified_parent_header_report_sha256") != header_sha
        or raw.get("observed_galaxy_rows_parsed") is not False
        or raw.get("observed_sector_columns_parsed") is not False
        or raw.get("observed_odd_data_vector_read") is not False
        or header.get("status") !=
           "OFFICIAL_LRG_FULL_HEADER_PARSED_SECTOR_SEMANTICS_NOT_YET_CERTIFIED"
        or header.get("official_url") != p["official_source_url"]
        or header.get("header_range_sha256") != p["first_header_sha256"]
        or header.get("header_bytes_requested") != p["first_header_bytes"]
        or header.get("declared_bintable_rows_NOT_READ") != p["expected_bintable_rows"]
        or header.get("declared_bintable_row_bytes") != p["expected_bintable_row_bytes"]
        or header.get("column_count") != p["expected_bintable_columns"]
        or header.get("observed_galaxy_data_rows_read") is not False
        or header.get("observed_odd_data_vector_read") is not False
        or p["only_permitted_observed_row_columns"] !=
           {"SECTOR": "J", "sector_TSR": "D",
            "sector_SSR": "D", "COMP_BOSS": "D"}
        or p["fixed_numeric_abs_tolerance"] != 1e-10
        or p["observed_odd_data_vector_read"] is not False
        or p["new_selection_rule_applied"] is not False
    ):
        raise ValueError("Frozen LRG source/full-SHA/local-report/header contract mismatch")
    file = ROOT / p["local_source_file"]
    if raw.get("local_file_path") != str(file.resolve()):
        raise ValueError("Local official FITS path differs from byte-only report")
    if file.stat().st_size != p["official_source_file_size_bytes"]:
        raise ValueError("Official LRG file byte count changed")
    return file, header


def validate_full_sha_before_rows(file, p):
    digest, header_digest, total = hashlib.sha256(), hashlib.sha256(), 0
    with file.open("rb") as fd:
        while True:
            buf = fd.read(1024 * 1024)
            if not buf:
                break
            digest.update(buf)
            if total < p["first_header_bytes"]:
                header_digest.update(
                    buf[:min(len(buf), p["first_header_bytes"] - total)]
                )
            total += len(buf)
            if total > p["official_source_file_size_bytes"]:
                raise ValueError("Official LRG file exceeds exact size cap")
        if (
            total != p["official_source_file_size_bytes"]
            or digest.hexdigest() !=
               p["first_seen_full_file_sha256_frozen_before_sector_rows"]
            or header_digest.hexdigest() != p["first_header_sha256"]
        ):
            raise ValueError("Official LRG full/header SHA mismatch BEFORE any row decoding")
        fd.seek(0)
        header_bytes = fd.read(p["first_header_bytes"])
    return header_bytes, digest.hexdigest()


def fits_width(form):
    match = TFORM.fullmatch(form.strip().upper())
    if match is None:
        raise ValueError("Unsupported FITS table TFORM in SHA-pinned header")
    count = int(match.group(1) or "1")
    code = match.group(2)
    if code == "X":
        return (count + 7) // 8
    if code not in FITS_WIDTH:
        raise ValueError("Unsupported FITS table field type in SHA-pinned header")
    return count * FITS_WIDTH[code]


def parse_schema(header_bytes, p, prior):
    primary, after_primary = parse_header(header_bytes, 0)
    if primary.get("SIMPLE") is not True or hdu_data_length(primary) != 0:
        raise ValueError("Unexpected official primary FITS data")
    table, data_start = parse_header(header_bytes, after_primary)
    if (
        str(table.get("XTENSION", "")).upper() != "BINTABLE"
        or data_start != p["first_header_bytes"]
        or int(table["NAXIS1"]) != p["expected_bintable_row_bytes"]
        or int(table["NAXIS2"]) != p["expected_bintable_rows"]
        or int(table["TFIELDS"]) != p["expected_bintable_columns"]
    ):
        raise ValueError("Pinned LRG BINTABLE geometry changed")
    nfields = p["expected_bintable_columns"]
    names = [str(table.get(f"TTYPE{i}", "")).strip()
             for i in range(1, nfields + 1)]
    forms = [str(table.get(f"TFORM{i}", "")).strip().upper()
             for i in range(1, nfields + 1)]
    if names != prior["column_names"]:
        raise ValueError("LRG FITS field names differ from first header report")
    if len(forms) != len(names) or len(set(names)) != len(names):
        raise ValueError("Malformed LRG FITS column names/forms")
    mapping = {}
    offset = 0
    for name, form in zip(names, forms):
        mapping[name] = (offset, form)
        offset += fits_width(form)
    if offset != p["expected_bintable_row_bytes"]:
        raise ValueError("LRG FITS column widths not equal to NAXIS1")
    for name, expected in p["only_permitted_observed_row_columns"].items():
        if name not in mapping or mapping[name][1] != expected:
            raise ValueError("Required LRG sector-only field type missing: " + name)
    if (data_start + p["expected_bintable_rows"] *
            p["expected_bintable_row_bytes"]
            > p["official_source_file_size_bytes"]):
        raise ValueError("Table exceeds official full-file byte cap")
    return data_start, mapping


def read_four_fields(mm, data_start, p, offsets):
    # Explicit four-field whitelist: do not create whole-structured FITS rows.
    spec = [
        ("SECTOR", ">i4", np.int64),
        ("sector_TSR", ">f8", np.float64),
        ("sector_SSR", ">f8", np.float64),
        ("COMP_BOSS", ">f8", np.float64),
    ]
    if [x[0] for x in spec] != list(p["only_permitted_observed_row_columns"]):
        raise ValueError("Observed FITS field decoding whitelist changed")
    cols = {}
    for name, dtype, target_type in spec:
        pos, _ = offsets[name]
        view = np.ndarray(
            shape=(p["expected_bintable_rows"],),
            dtype=np.dtype(dtype),
            buffer=mm,
            offset=data_start + pos,
            strides=(p["expected_bintable_row_bytes"],),
        )
        cols[name] = view.astype(target_type, copy=True)
    return cols


def aggregates(cols, p):
    names = list(p["only_permitted_observed_row_columns"])
    if list(cols) != names:
        raise ValueError("Observed sector-only field names differ from whitelist")
    n = len(cols["SECTOR"])
    if n != p["expected_bintable_rows"] or any(
        x.ndim != 1 or len(x) != n for x in cols.values()
    ):
        raise ValueError("Invalid LRG sector-field shapes/number of rows")
    tol = p["fixed_numeric_abs_tolerance"]
    vals = [cols[name] for name in names[1:]]
    if any(not np.isfinite(x).all() for x in vals):
        raise ValueError("Nonfinite official sector numeric values; no imputation")
    sectors, inv, freq = np.unique(
        cols["SECTOR"], return_inverse=True, return_counts=True
    )
    if len(sectors) < 1 or int(np.sum(freq)) != n:
        raise ValueError("Sector grouping does not cover all official rows")
    stats = {}
    constancy = {}
    for name, x in zip(names[1:], vals):
        stats[name] = {
            "finite_count": int(np.count_nonzero(np.isfinite(x))),
            "minimum": float(np.min(x)),
            "median": float(np.median(x)),
            "maximum": float(np.max(x)),
            "exact_zero_count": int(np.count_nonzero(x == 0)),
            "less_or_equal_0_5_DIAGNOSTIC_NOT_A_CUT":
                int(np.count_nonzero(x <= 0.5)),
            "greater_than_1_count": int(np.count_nonzero(x > 1)),
            "less_than_0_count": int(np.count_nonzero(x < 0)),
        }
        small = np.full(len(sectors), np.inf)
        large = np.full(len(sectors), -np.inf)
        np.minimum.at(small, inv, x)
        np.maximum.at(large, inv, x)
        constancy[name] = {
            "sector_labels_with_span_gt_fixed_tol":
                int(np.count_nonzero(large - small > tol)),
            "largest_within_sector_span": float(np.max(large - small)),
        }
    tsr, ssr, comp = vals
    diffs = {
        "COMP_BOSS - sector_TSR": comp - tsr,
        "COMP_BOSS - sector_TSR * sector_SSR": comp - tsr * ssr,
        "sector_TSR - sector_SSR": tsr - ssr,
    }
    if list(diffs) != p["candidate_formula_diagnostics_NO_SEMANTIC_CONCLUSION"]:
        raise ValueError("Frozen numeric diagnostic expressions changed")
    relation = {
        label: {
            "abs_difference_le_fixed_tol_count":
                int(np.count_nonzero(np.abs(diff) <= tol)),
            "max_abs_difference": float(np.max(np.abs(diff))),
            "median_abs_difference": float(np.median(np.abs(diff))),
        }
        for label, diff in diffs.items()
    }
    return {
        "status": "OFFICIAL_LRG_SECTOR_METADATA_AGGREGATES_ONLY_NOT_SEMANTICS_CERTIFIED",
        "full_catalogue_rows": n,
        "distinct_integer_sector_labels": int(len(sectors)),
        "sector_label_less_or_equal_zero_count":
            int(np.count_nonzero(cols["SECTOR"] <= 0)),
        "rows_per_sector_minimum": int(np.min(freq)),
        "rows_per_sector_median": float(np.median(freq)),
        "rows_per_sector_maximum": int(np.max(freq)),
        "column_aggregate_diagnostics": stats,
        "within_exact_sector_label_constancy": constancy,
        "fixed_candidate_formula_diagnostics_NO_SEMANTIC_CONCLUSION": relation,
        "fixed_numeric_abs_tolerance": tol,
        "observed_row_columns_decoded": names,
        "individual_sector_ids_or_galaxy_rows_output": False,
        "observed_coordinates_or_redshifts_decoded": False,
        "observed_galaxy_weight_columns_decoded": False,
        "observed_random_or_mock_rows_decoded": False,
        "observed_odd_data_vector_read": False,
        "new_selection_rule_applied": False,
        "lrg_sector_completeness_semantics_certified": False,
        "official_lrg_veto_polygon_composition_certified": False,
        "physical_lrg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "note": (
            "Aggregate sector-field consistency of published full catalogue, "
            "not proof of C_eBOSS/C_z field equivalence, historical MANGLE "
            "production, random weights, published clustering cuts, physical "
            "continuous mask, pair window or odd-sector inference."
        ),
        "errors": [],
    }


def audit(p):
    path, header_report = frozen_gate(p)
    header, fullsha = validate_full_sha_before_rows(path, p)
    data_start, offsets = parse_schema(header, p, header_report)
    mm = np.memmap(path, dtype="u1", mode="r")
    try:
        cols = read_four_fields(mm, data_start, p, offsets)
        result = aggregates(cols, p)
    finally:
        del mm
    result.update({
        "full_file_sha256_reverified_before_rows": fullsha,
        "first_header_sha256_reverified_before_rows": p["first_header_sha256"],
        "frozen_previous_byte_only_report_sha256":
            p["parent_lrg_full_bytes_report_sha256"],
    })
    return result


def self_test():
    p = read_json(PROTOCOL)
    if (p["first_seen_full_file_sha256_frozen_before_sector_rows"] !=
            "39b831801adec04fe6dc5d6ab76a4b303aa58bfd303548cb7fddfae7d7f9331d"
            or p["only_permitted_observed_row_columns"] !=
               {"SECTOR": "J", "sector_TSR": "D",
                "sector_SSR": "D", "COMP_BOSS": "D"}):
        raise AssertionError("Frozen input SHA/field whitelist changed")
    # 630-byte synthetic rows: 4+8+8+8 sector fields plus 602 ignored bytes.
    fake = b"".join(
        struct.pack(">iddd", sector, tsr, ssr, comp) + b"x" * 602
        for sector, tsr, ssr, comp in [
            (3, 0.8, 0.7, 0.56), (3, 0.81, 0.7, 0.567),
            (4, 0.9, 0.9, 0.81), (4, 0.9, 0.9, 0.81),
        ]
    )
    synthetic = dict(p, expected_bintable_rows=4)
    offsets = {name: (i, form) for name, i, form in [
        ("SECTOR", 0, "J"),
        ("sector_TSR", 4, "D"),
        ("sector_SSR", 12, "D"),
        ("COMP_BOSS", 20, "D"),
    ]}
    cols = read_four_fields(fake, 0, synthetic, offsets)
    r = aggregates(cols, synthetic)
    assert r["full_catalogue_rows"] == 4
    assert r["distinct_integer_sector_labels"] == 2
    assert r["within_exact_sector_label_constancy"]["sector_TSR"][
        "sector_labels_with_span_gt_fixed_tol"
    ] == 1
    assert r["fixed_candidate_formula_diagnostics_NO_SEMANTIC_CONCLUSION"][
        "COMP_BOSS - sector_TSR * sector_SSR"
    ]["abs_difference_le_fixed_tol_count"] == 4
    assert r["observed_odd_data_vector_read"] is False
    try:
        broken = dict(cols, sector_TSR=np.array([0.8, np.nan, 0.9, 0.9]))
        aggregates(broken, synthetic)
    except ValueError:
        pass
    else:
        raise AssertionError("Nonfinite sector diagnostic must fail closed")
    print("EBOSS_LRG_SECTOR_FOUR_COLUMN_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    p = read_json(PROTOCOL)
    output = args.out or ROOT / p["report_output"]
    try:
        result = audit(p)
    except Exception as e:
        result = {
            "status": "OFFICIAL_LRG_SECTOR_METADATA_INCOMPLETE_STOP",
            "errors": [str(e)],
            "observed_odd_data_vector_read": False,
            "physical_lrg_mask_certified": False,
            "new_selection_rule_applied": False,
        }
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)
    print("EBOSS_LRG_SECTOR_METADATA", result["status"], flush=True)
    print("REPORT", output, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep="\n", flush=True)
        return 2
    print("ROWS", result["full_catalogue_rows"], flush=True)
    print("DISTINCT_SECTORS", result["distinct_integer_sector_labels"], flush=True)
    print("WITHIN_SECTOR_CONSTANCY",
          json.dumps(result["within_exact_sector_label_constancy"], sort_keys=True),
          flush=True)
    print("FIELD_DIAGNOSTICS",
          json.dumps(result["column_aggregate_diagnostics"], sort_keys=True),
          flush=True)
    print("FIXED_FORMULAE",
          json.dumps(result["fixed_candidate_formula_diagnostics_NO_SEMANTIC_CONCLUSION"],
                     sort_keys=True), flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
