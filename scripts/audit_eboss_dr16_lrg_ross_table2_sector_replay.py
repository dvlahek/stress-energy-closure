#!/usr/bin/env python3
"""Fixed Ross Table 2 LRG sector-tally replay, not a new science selection.

This is a source-informed POST-AGGREGATE check: the first-stage count was
already observed before this follow-up was declared. Only four previously
allowed FITS columns are accessed; the observed odd sector remains closed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

from audit_eboss_dr16_lrg_sector_metadata import (
    ROOT, frozen_gate, parse_schema, read_four_fields,
    read_json, validate_full_sha_before_rows,
)

PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_ross_table2_sector_replay_protocol_2026-09-25.json"
PARENT = ROOT / "source_data/eboss_dr16_lrg_sector_metadata_protocol_2026-09-25.json"


def tally(cols, p):
    ref = p["published_reference"]
    n = ref["paper_full_postveto_target_count"]
    if (list(cols) != p["only_permitted_fields"]
            or len(cols["COMP_BOSS"]) != n
            or any(a.ndim != 1 or len(a) != n for a in cols.values())
            or any(not np.isfinite(cols[k]).all()
                   for k in p["only_permitted_fields"][1:])):
        raise ValueError("Only frozen four finite sector fields are permitted")
    keep = cols[p["first_stage_field"]] > p["threshold"]
    first_bad = int(np.count_nonzero(~keep))
    variants = {}
    for name in p["second_stage_candidates"]:
        second_bad = int(np.count_nonzero(keep & (cols[name] <= p["threshold"])))
        final = int(np.count_nonzero(keep & (cols[name] > p["threshold"])))
        if first_bad + second_bad + final != n:
            raise ValueError("Published cohort accounting failed")
        variants[name] = {
            "second_rejected_after_first": second_bad,
            "final_retained": final,
            "second_difference_from_paper":
                second_bad - ref["paper_subsequent_C_z_less_or_equal_0_5_removed_targets"],
            "final_difference_from_paper":
                final - ref["paper_after_both_target_count"],
        }
    return {
        "status": "LRG_PUBLISHED_TABLE2_SEQUENTIAL_AGGREGATE_CROSSCHECK_ONLY",
        "post_aggregate_source_informed": True,
        "full_catalogue_rows": n,
        "first_stage_rejected_COMP_BOSS_LE_0_5": first_bad,
        "first_stage_diff_from_paper":
            first_bad - ref["paper_angular_C_eBOSS_less_or_equal_0_5_removed_targets"],
        "first_stage_retained": n - first_bad,
        "fixed_second_candidate_results": variants,
        "production_code_authenticated": False,
        "new_science_cut_applied": False,
        "physical_mask_certified": False,
        "observed_odd_data_vector_read": False,
        "errors": [],
    }


def preflight(p):
    prior_file = ROOT / p["uploaded_four_column_report"]
    archived = prior_file.read_bytes()
    local = (ROOT / "eboss_workspace/official_mask_inventory/lrg_full_sector_metadata_audit.json").read_bytes()
    if (archived != local
            or hashlib.sha256(archived).hexdigest() != p["exact_four_column_report_sha256"]
            or hashlib.sha1(f"blob {len(archived)}\\0".encode() + archived).hexdigest()
                != p["four_column_report_git_blob_sha1"]):
        raise ValueError("Original uploaded sector report changed")
    prior = json.loads(archived)
    parent = read_json(PARENT)
    if (prior["full_catalogue_rows"] != p["published_reference"]["paper_full_postveto_target_count"]
            or prior["full_file_sha256_reverified_before_rows"] != p["source_fits_sha256"]
            or prior["observed_row_columns_decoded"] != p["only_permitted_fields"]
            or prior["observed_odd_data_vector_read"] is not False
            or parent["first_seen_full_file_sha256_frozen_before_sector_rows"] != p["source_fits_sha256"]
            or p["first_stage_field"] != "COMP_BOSS"
            or p["second_stage_candidates"] != ["sector_SSR", "sector_TSR"]
            or p["threshold"] != 0.5):
        raise ValueError("Fixed source-informed follow-up protocol mismatch")
    return parent


def self_test(p):
    x = dict(p)
    x["published_reference"] = dict(p["published_reference"],
        paper_full_postveto_target_count=4,
        paper_angular_C_eBOSS_less_or_equal_0_5_removed_targets=1,
        paper_subsequent_C_z_less_or_equal_0_5_removed_targets=1,
        paper_after_both_target_count=2)
    cols = {"SECTOR": np.array([1, 2, 3, 4]),
            "sector_TSR": np.array([.9, .9, .9, .9]),
            "sector_SSR": np.array([.1, .2, .9, .9]),
            "COMP_BOSS": np.array([.2, .9, .9, .9])}
    r = tally(cols, x)
    assert r["first_stage_rejected_COMP_BOSS_LE_0_5"] == 1
    assert r["fixed_second_candidate_results"]["sector_SSR"]["second_rejected_after_first"] == 1
    assert r["fixed_second_candidate_results"]["sector_TSR"]["second_rejected_after_first"] == 0
    print("EBOSS_LRG_ROSS_TABLE2_SYNTHETIC_SELF_TEST_OK", flush=True)


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
        parent = preflight(p)
        source, header_report = frozen_gate(parent)
        first_header, full_sha = validate_full_sha_before_rows(source, parent)
        start, offsets = parse_schema(first_header, parent, header_report)
        mm = np.memmap(source, dtype="u1", mode="r")
        try:
            cols = read_four_fields(mm, start, parent, offsets)
            report = tally(cols, p)
        finally:
            del mm
        report["full_file_sha256_reverified_before_rows"] = full_sha
        report["previous_sector_report_sha256"] = p["exact_four_column_report_sha256"]
    except Exception as exc:
        report = {"status": "LRG_TABLE2_REPLAY_INCOMPLETE_STOP",
                  "errors": [str(exc)], "observed_odd_data_vector_read": False,
                  "new_science_cut_applied": False}
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(report, indent=2) + "\\n", encoding="utf-8")
    tmp.replace(out)
    print("EBOSS_LRG_ROSS_TABLE2_REPLAY", report["status"], flush=True)
    print("REPORT", out, flush=True)
    if report["errors"]:
        print("ERRORS", *report["errors"], sep="\\n", flush=True)
        return 2
    print("FIRST_STAGE", report["first_stage_rejected_COMP_BOSS_LE_0_5"], flush=True)
    print("SECOND_STAGE", json.dumps(report["fixed_second_candidate_results"], sort_keys=True), flush=True)
    print("ODD_DATA_READ", report["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
