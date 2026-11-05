#!/usr/bin/env python3
"""Offline eBOSS DR16 paper-cited brickmask v1.0 source-version QA.

Only exact MIT-licensed vendored upstream C source at public tag v1.0
and later pinned public commit is inspected, plus synthetic bit values.
No original production invocation, FITS, image, galaxy, random,
observed odd-sector data or mask selection is accessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_brickmask_v1_0_source_comparison_protocol_2026-09-26.json"
PASS = "PAPER_CITED_BRICKMASK_V1_0_SOURCE_RECONCILED_ONLY"
STOP = "PAPER_CITED_BRICKMASK_V1_0_SOURCE_RECONCILIATION_INCOMPLETE_STOP"


def sha_git_blob(b):
    return hashlib.sha1(
        ("blob " + str(len(b))).encode("ascii") + bytes((0,)) + b
    ).hexdigest()


def strict_source_check(m, original=None, later=None):
    v = m["historical_tag"]
    n = m["later_public_checkout"]
    if (
        m["upstream_repository"] != "cheng-zhao/brickmask"
        or v["name"] != "v1.0"
        or v["annotated_tag_object_sha"] != "3a94b72f4ee39c215cab713cf464d54106899838"
        or v["tagged_commit"] != "4c0f940934ff8c4e6b0f0c709b3093383abff8d0"
        or v["tag_message"] != "Version for eBOSS DR16"
        or n["commit"] != "b9eb684a579b56ec3dbdb46549224be7e3fa2830"
        or n["required_compile_flag"] != "-DEBOSS"
        or m["new_selection_applied"] is not False
        or m["observed_odd_data_vector_read"] is not False
    ):
        raise ValueError("Paper-cited 2020 tag, newer commit or source-only contract changed")
    root_v = ROOT / m["archived_source_root"]
    root_n = ROOT / m["archived_later_source_root"]
    original = original or {}
    later = later or {}
    for rel, expected in v["source_files"].items():
        b = original.get(rel)
        if b is None:
            b = (root_v / rel).read_bytes()
        if sha_git_blob(b) != expected:
            raise ValueError("2020 v1.0 source Git blob changed: " + rel)
        original[rel] = b
    for rel, expected in n["source_files"].items():
        b = later.get(rel)
        if b is None:
            b = (root_n / rel).read_bytes()
        if sha_git_blob(b) != expected:
            raise ValueError("2023 public source Git blob changed: " + rel)
        later[rel] = b
    old_define = original["define.h"].decode("utf-8")
    old_main = original["brickmask.c"].decode("utf-8")
    old_wcs = original["read_data.c"].decode("utf-8")
    new_define = later["src/define.h"].decode("utf-8")
    new_assign = later["src/bit_code.c"].decode("utf-8")
    # The checks intentionally fail CLOSED upon source refactoring.
    required = [
        (old_define, r"#define\s+NUMSUB\s+4\b"),
        (old_define, r"#define\s+MASK_VALID\s+\(bit\s*&\s*1\)"),
        (old_define, r"#define\s+XYBUG\s+4\b"),
        (old_define, r"#define\s+XYBUG_VALID\s+\(bit\s*>>\s*2\s*&\s*1\)"),
        (old_define, r'#define\s+MSK_COL\s+"VETOMASK"'),
        (old_define, r'#define\s+CHK_COL\s+"MCHUNK"'),
        (old_main, r"j\s*=\s*\(int\)\s*round\(x\)\s*;"),
        (old_main, r"k\s*=\s*\(int\)\s*round\(y\)\s*;"),
        (old_main, r"if\s*\(MASK_VALID\)"),
        (old_main, r"data\[n\]\.mask\s*\+=\s*bit\s*-\s*XYBUG\s*;"),
        (old_main, r"j\s*=\s*\(int\)\s*x\s*;"),
        (old_main, r"k\s*=\s*\(int\)\s*y\s*;"),
        (old_main, r"if\s*\(XYBUG_VALID\)\s*data\[n\]\.mask\s*\+=\s*XYBUG\s*;"),
        (old_main, r"data\[n\]\.flag\s*\+=\s*1\s*<<\s*i\s*;"),
        (old_main, r"fits_read_tblbytes\(fptr,\s*1,\s*1,\s*n,\s*maskbit\.bit"),
        (old_wcs, r"only CTYPE1='RA---TAN' is supported"),
        (old_wcs, r"only CTYPE2='DEC--TAN' is supported"),
        (new_define, r"#define\s+EBOSS_MASK_VALID\(bit\)\s+\(\(bit\)\s*&\s*1\)"),
        (new_define, r"#define\s+EBOSS_XYBUG_BIT\s+4\b"),
        (new_define, r"#define\s+EBOSS_XYBUG_VALID\(bit\)\s+\(\(bit\)\s*&\s*\(EBOSS_XYBUG_BIT\)\)"),
        (new_assign, r"long\s+rx\s*=\s*round\(x\)\s*;"),
        (new_assign, r"long\s+ry\s*=\s*round\(y\)\s*;"),
        (new_assign, r"if\s*\(!\(EBOSS_MASK_VALID\(bit\)\)\)\s*continue\s*;"),
        (new_assign, r"data->mask\[i\]\s*\+=\s*bit\s*-\s*EBOSS_XYBUG_BIT\s*;"),
        (new_assign, r"long\s+ix\s*=\s*\(long\)\s*x\s*;"),
        (new_assign, r"long\s+iy\s*=\s*\(long\)\s*y\s*;"),
        (new_assign, r"if\s*\(EBOSS_XYBUG_VALID\(bit\)\)\s*data->mask\[i\]\s*\+=\s*EBOSS_XYBUG_BIT\s*;"),
        (new_assign, r"if\s*\(data->subid\)\s*data->subid\[i\]\s*=\s*subid\s*;"),
    ]
    for source, expr in required:
        if not re.search(expr, source):
            raise ValueError("Required source-level eBOSS version semantic anchor absent: " + expr)
    family = re.findall(r'"(mask-eboss[0-9]+-%s[.]fits[.]gz)"', old_define)
    if family != v["legacy_four_filenames"]:
        raise ValueError("2020 DR16 release tag eBOSS four-family definition changed")
    if v["legacy_mask_output_column"] != "VETOMASK" or v["legacy_chunk_output_column"] != "MCHUNK":
        raise ValueError("Legacy field names changed")
    return {
        "source_version_tag":"v1.0",
        "annotated_git_tag_object_sha":v["annotated_tag_object_sha"],
        "tagged_commit_sha":v["tagged_commit"],
        "tag_message":v["tag_message"],
        "later_public_commit_sha":n["commit"],
        "all_eight_vendored_MIT_source_blobs_byte_identical":True,
        "v1_0_exact_four_chunk_filename_order":family,
        "v1_0_field_code":"VETOMASK + MCHUNK 1<<i bitset",
        "later_public_field_code":"maskbit column configurable + optional last-valid SUBID",
        "both_reviewed_C_paths_round_validity_bit0_and_truncate_XYBUG_bit2":True,
        "complete_source_program_equivalence_claimed":False,
        "historical_DR16_production_invocation_authenticated":False,
    }


def compare_byte_code_paths():
    """Exhaustive uint8 pair algebra ONLY, conditioned on same two image bytes."""
    n=0
    for rounded in range(256):
        for truncated in range(256):
            old = None
            if rounded & 1:
                old = (rounded - 4 if ((rounded >> 2) & 1) else rounded)
                if (truncated >> 2) & 1:
                    old += 4
            new = None
            if rounded & 1:
                new = (rounded - 4 if rounded & 4 else rounded)
                if truncated & 4:
                    new += 4
            if old != new:
                raise AssertionError("Reviewed old/new uint8 bit0/bit2 paths disagree")
            n += 1
    # For one and the same set of valid chunk passes, mask codes sum
    # in both implementations. This does NOT prove WCS/pixel selection.
    overlap_legacy_flag=sum(1 << i for i in [0,1])
    later_optional_subid=1   # second valid family overwrites first.
    if overlap_legacy_flag != 3 or later_optional_subid != 1:
        raise AssertionError("Different multi-chunk output semantics were lost")
    return {
        "synthetic_bitbyte_pairs_exhaustively_compared":n,
        "rounded_and_truncated_byte_value_algebra_identical":True,
        "two_valid_chunk_example_legacy_MCHUNK":overlap_legacy_flag,
        "two_valid_chunk_example_later_SUBID":later_optional_subid,
        "overlap_outputs_not_interchangeable":True,
        "coordinate_WCS_equivalence_NOT_TESTED":True,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path)
    args=ap.parse_args()
    m=json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if args.self_test:
        res=strict_source_check(m)
        assert res["all_eight_vendored_MIT_source_blobs_byte_identical"]
        synth=compare_byte_code_paths()
        assert synth["synthetic_bitbyte_pairs_exhaustively_compared"]==65536
        old=(ROOT/m["archived_source_root"]/"define.h").read_bytes()
        try:
            strict_source_check(m,original={"define.h":old.replace(
                b'#define XYBUG 4',b'#define XYBUG 8'
            )})
        except ValueError:
            pass
        else:
            raise AssertionError("Tampered historical bit2 macro passed")
        print("EBOSS_ELG_BRICKMASK_V1_0_SOURCE_SYNTHETIC_SELF_TEST_OK",flush=True)
        return 0
    out=args.out or ROOT/"eboss_workspace/official_mask_inventory/brickmask_v1_0_source_reconciliation.json"
    try:
        source=strict_source_check(m)
        synth=compare_byte_code_paths()
        result={
            "status":PASS,
            "source":source,
            "synthetic_source_conditioned_uint8_pair_check":synth,
            "no_FITS_image_catalogue_random_or_mock_data_read":True,
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
            "continuous_ELG_mask_certified":False,
            "LRG_ELG_pair_window_certified":False,
            "errors":[],
        }
    except Exception as exc:
        result={
            "status":STOP,
            "errors":[str(exc)],
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
        }
    out.parent.mkdir(parents=True,exist_ok=True)
    tmp=out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result,indent=2)+chr(10),encoding="utf-8")
    tmp.replace(out)
    print("EBOSS_ELG_V1_0_SOURCE_RECONCILIATION",result["status"],flush=True)
    print("REPORT",out,flush=True)
    if result["errors"]:
        print("ERRORS",*result["errors"],sep=chr(10),flush=True)
        return 2
    print("SYNTHETIC_UINT8_PAIRS",
          result["synthetic_source_conditioned_uint8_pair_check"]["synthetic_bitbyte_pairs_exhaustively_compared"],
          flush=True)
    print("HISTORICAL_EXECUTION_AUTHENTICATED",
          result["source"]["historical_DR16_production_invocation_authenticated"],flush=True)
    print("ODD_DATA_READ",result["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
