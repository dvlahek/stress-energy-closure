#!/usr/bin/env python3
"""Source-informed Ross Table 2 LRG sector-count replay (aggregate-only).

This comparison was specified after the first-stage COMP_BOSS count was
seen. It is not an independent prospective test of that first stage.
Only the previously allowed four FITS columns are decoded after full SHA.
No released science selection, MANGLE geometry or odd vector is changed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_eboss_dr16_lrg_sector_metadata import (
    ROOT, frozen_gate, parse_schema, read_four_fields, read_json,
    validate_full_sha_before_rows,
)

PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_ross_table2_sector_replay_protocol_2026-09-25.json"
PARENT = ROOT / "source_data/eboss_dr16_lrg_sector_metadata_protocol_2026-09-25.json"
MANIFEST = ROOT / "source_data/eboss_dr16_lrg_sector_metadata_uploaded_manifest_2026-09-25.json"


def preflight(p, parent):
    manifest = read_json(MANIFEST)
    archive = (ROOT / p["uploaded_four_column_report"]).read_bytes()
    local = (ROOT / manifest["source_report_path"]).read_bytes()
    report_sha = hashlib.sha256(archive).hexdigest()
    if (
        archive != local
        or report_sha != p["exact_four_column_report_sha256"]
        or report_sha != manifest["uploaded_report_sha256"]
        or p["source_fits_sha256"] != parent["first_seen_full_file_sha256_frozen_before_sector_rows"]
        or p["source_header_sha256"] != parent["first_header_sha256"]
        or p["source_fits_path"] != parent["local_source_file"]
        or p["only_permitted_fields"] != list(parent["only_permitted_observed_row_columns"])
        or p["only_permitted_fields"] != ["SECTOR", "sector_TSR", "sector_SSR", "COMP_BOSS"]
        or p["first_stage_field"] != "COMP_BOSS"
        or p["second_stage_candidates"] != ["sector_SSR", "sector_TSR"]
        or p["threshold"] != 0.5
        or p["no_optimization_or_selection_changes"] is not True
        or p["observed_odd_data_vector_read"] is not False
    ):
        raise ValueError("Frozen source-informed follow-up provenance changed")
    r = json.loads(archive)
    if (
        r.get("status") != "OFFICIAL_LRG_SECTOR_METADATA_AGGREGATES_ONLY_NOT_SEMANTICS_CERTIFIED"
        or r.get("full_catalogue_rows") != p["published_reference"]["paper_full_postveto_target_count"]
        or r.get("observed_row_columns_decoded") != p["only_permitted_fields"]
        or r.get("full_file_sha256_reverified_before_rows") != p["source_fits_sha256"]
        or r.get("observed_odd_data_vector_read") is not False
        or r.get("new_selection_rule_applied") is not False
        or r.get("errors") != []
        or r.get("column_aggregate_diagnostics", {}).get("COMP_BOSS", {}).get(
            "less_or_equal_0_5_DIAGNOSTIC_NOT_A_CUT"
        ) != p["post_aggregate_observed_field_counts_prior_to_this_protocol"][
            "COMP_BOSS_less_or_equal_0_5_rows"
        ]
    ):
        raise ValueError("Archived first sector audit disagrees with frozen protocol")
    return report_sha


def counts(cols, p):
    ref = p["published_reference"]
    n = ref["paper_full_postveto_target_count"]
    if (
        list(cols) != p["only_permitted_fields"]
        or any(v.ndim != 1 or len(v) != n for v in cols.values())
        or any(
            not np.isfinite(cols[name]).all()
            or np.any(cols[name] < 0)
            or np.any(cols[name] > 1)
            for name in p["only_permitted_fields"][1:]
        )
    ):
        raise ValueError("Four-column whitelist, row count or finite [0,1] domain failed")
    first_keep = cols[p["first_stage_field"]] > p["threshold"]
    first_rejected = int(np.count_nonzero(~first_keep))
    first_retained = int(np.count_nonzero(first_keep))
    variants = {}
    for name in p["second_stage_candidates"]:
        rejected = int(np.count_nonzero(first_keep & (cols[name] <= p["threshold"])))
        retained = int(np.count_nonzero(first_keep & (cols[name] > p["threshold"])))
        if first_rejected + rejected + retained != n:
            raise ValueError("Sequential target tally is incomplete")
        variants[name] = {
            "rejected_after_first_stage": rejected,
            "retained_after_both_stages": retained,
            "rejected_difference_from_paper": rejected
                - ref["paper_subsequent_C_z_less_or_equal_0_5_removed_targets"],
            "retained_difference_from_paper": retained
                - ref["paper_after_both_target_count"],
        }
    return {
        "status": "LRG_ROSS_TABLE2_SOURCE_INFORMED_AGGREGATE_REPLAY_ONLY",
        "full_catalogue_rows": n,
        "first_stage_comp_boss_rejected": first_rejected,
        "first_stage_comp_boss_retained": first_retained,
        "paper_first_stage_rejected": ref[
            "paper_angular_C_eBOSS_less_or_equal_0_5_removed_targets"
        ],
        "first_stage_difference_from_paper": first_rejected
            - ref["paper_angular_C_eBOSS_less_or_equal_0_5_removed_targets"],
        "paper_subsequent_rejected": ref[
            "paper_subsequent_C_z_less_or_equal_0_5_removed_targets"
        ],
        "paper_after_both_retained": ref["paper_after_both_target_count"],
        "fixed_second_field_candidates": variants,
        "post_aggregate_source_informed": True,
        "historical_production_code_authenticated": False,
        "physical_lrg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "observed_row_columns_decoded": p["only_permitted_fields"],
        "observed_redshift_coordinates_weight_ID_read": False,
        "observed_random_or_mock_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_science_selection_applied": False,
        "errors": [],
    }


def self_test(p):
    synthetic = dict(p)
    synthetic["published_reference"] = dict(
        p["published_reference"],
        paper_full_postveto_target_count=4,
        paper_angular_C_eBOSS_less_or_equal_0_5_removed_targets=1,
        paper_subsequent_C_z_less_or_equal_0_5_removed_targets=1,
        paper_after_both_target_count=2,
    )
    fields = {
        "SECTOR": np.array([1, 2, 3, 4], dtype="i8"),
        "sector_TSR": np.array([0.9, 0.9, 0.9, 0.9]),
        "sector_SSR": np.array([0.1, 0.2, 0.9, 0.9]),
        "COMP_BOSS": np.array([0.5, 0.9, 0.9, 0.9]),
    }
    result = counts(fields, synthetic)
    assert result["first_stage_comp_boss_rejected"] == 1
    assert result["fixed_second_field_candidates"]["sector_SSR"][
        "rejected_after_first_stage"
    ] == 1
    assert result["fixed_second_field_candidates"]["sector_SSR"][
        "retained_after_both_stages"
    ] == 2
    assert result["fixed_second_field_candidates"]["sector_TSR"][
        "rejected_after_first_stage"
    ] == 0
    assert result["observed_odd_data_vector_read"] is False
    try:
        bad = dict(fields, sector_SSR=np.array([0.1, np.nan, 0.9, 0.9]))
        counts(bad, synthetic)
    except ValueError:
        pass
    else:
        raise AssertionError("Non-finite synthetic sector value was accepted")
    print("EBOSS_LRG_ROSS_TABLE2_SOURCE_INFORMED_SYNTHETIC_SELF_TEST_OK", flush=True)


def audit(p):
    parent = read_json(PARENT)
    report_sha = preflight(p, parent)
    path, header_report = frozen_gate(parent)
    raw_header, full_sha = validate_full_sha_before_rows(path, parent)
    start, offsets = parse_schema(raw_header, parent, header_report)
    mm = np.memmap(path, dtype="u1", mode="r")
    try:
        fields = read_four_fields(mm, start, parent, offsets)
        result = counts(fields, p)
    finally:
        del mm
    result["full_file_sha256_reverified_before_rows"] = full_sha
    result["first_sector_report_sha256_verified"] = report_sha
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    p = read_json(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    out = ROOT / p["report_output"]
    try:
        result = audit(p)
    except Exception as exc:
        result = {
            "status": "LRG_ROSS_TABLE2_SOURCE_INFORMED_REPLAY_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "new_science_selection_applied": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result, indent=2) + chr(10), encoding="utf-8")
    tmp.replace(out)
    print("EBOSS_LRG_ROSS_TABLE2_REPLAY", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep=chr(10), flush=True)
        return 2
    print("FIRST_STAGE_REJECTED", result["first_stage_comp_boss_rejected"], flush=True)
    print("SECOND_STAGE", json.dumps(result["fixed_second_field_candidates"], sort_keys=True), flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
