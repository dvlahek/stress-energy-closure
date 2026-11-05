#!/usr/bin/env python3
"""E17D2b3b4: user-supplied checksum listing source-only audit. NO ASDF I/O."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
S = ROOT / "source_data"
PRE = S / "eboss_dr16_a03_e17d2b3b4_user_globus_checksum_intake_prereg_2026-09-28.json"
RAW = S / "eboss_dr16_a03_e17d2b3b4_user_supplied_checksums.crc32"
ARCHIVE = S / "eboss_dr16_a03_e17d2b3b4_archived_source_only_2026_09_28"
PRE_BLOB = "a6ad151ff116d8b20208bebe76d4423db589b76f"
RAW_BLOB = "0bca2fd0a28c9e163ad6198f67616da4355f4c81"
PAT = re.compile(rb"([0-9]+) ([0-9]+) (halo_info_[0-9]{3}[.]asdf)")
NAMES = tuple(f"halo_info_{i:03}.asdf" for i in range(34))

def need(condition, reason):
    if not condition:
        raise ValueError("E17D2B3B4_FAIL_CLOSED: " + reason)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def gitblob(data):
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def canonical(value):
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")

def original_regex(payload):
    out = {}
    for line in payload.splitlines():
        match = PAT.fullmatch(line)
        need(match is not None, "invalid original POSIX cksum syntax")
        crc_b, size_b, name_b = match.groups()
        name = name_b.decode("ascii")
        need(name not in out, "duplicate filename")
        crc, size = int(crc_b), int(size_b)
        need(0 <= crc <= 4294967295 and 0 < size < 8 * 1024**3,
             "invalid CRC32 or individual file size")
        out[name] = (crc, size)
    need(tuple(sorted(out)) == NAMES, "missing/extra/unordered filename set")
    return out

def independent_tokenizer(payload):
    # Independent tokenizer without the regex/parser and without original records.
    text = payload.decode("ascii")
    entries = []
    for line in text.split("\n"):
        if not line:
            continue
        tokens = line.split(" ")
        need(len(tokens) == 3 and all(tokens[:2]) and
             tokens[0].isascii() and tokens[0].isdecimal() and
             tokens[1].isascii() and tokens[1].isdecimal(),
             "independent invalid three-field cksum record")
        crc = int(tokens[0])
        size = int(tokens[1])
        filename = tokens[2]
        need(0 <= crc < 2**32 and 0 < size < 8 * 1024**3,
             "independent invalid CRC or bytes")
        need(filename in NAMES, "independent unknown filename")
        need(filename not in [row[0] for row in entries], "independent duplicate")
        entries.append((filename, crc, size))
    need(len(entries) == 34 and set(x[0] for x in entries) == set(NAMES),
         "independent missing/extra entry")
    return entries

def audit(payload, prereg, enforce_pin=True):
    if enforce_pin:
        need(sha(payload) == prereg["user_intake"]["supplied_manifest_sha256"],
             "supplied manifest SHA256 differs from prospectively pinned upload")
        need(gitblob(payload) == RAW_BLOB, "archived uploaded checksum bytes drift")
        need(len(payload) == prereg["user_intake"]["supplied_manifest_exact_size_bytes"],
             "archived upload size drift")
    o = original_regex(payload)
    i = independent_tokenizer(payload)
    need(sorted((name, *pair) for name, pair in o.items()) == sorted(i),
         "original and independent disagree on individual records")
    o_sum = sum(v[1] for v in o.values())
    i_sum = sum(row[2] for row in i)
    target = prereg["source_snapshot"]["published_table_expected_aggregate_bytes"]
    need(o_sum == i_sum and o_sum + len(payload) == target,
         "checksum 34-file-plus-manifest byte aggregate disagrees with source-pinned static portal")
    chosen = prereg["locked_selected_file"]
    selected = prereg["user_intake"]["intake_only_first_record"]
    need(chosen == selected["filename"] and
         o[chosen] == (selected["gnu_posix_cksum_crc"], selected["byte_count"]),
         "preregistered exact selected superslab checksum/size drift")
    need(len(o) + 1 == prereg["source_snapshot"]["published_table_expected_directory_entries"],
         "actual user listing count does not account for archived source directory entries")
    return o, i, o_sum

def negative_controls(raw, prereg):
    cases = {
        "ONE_BYTE_SHA_PIN_REJECT": raw.replace(b"3502514872", b"3502514873", 1),
        "DUPLICATE_ENTRY_REJECT": raw + raw.splitlines(keepends=True)[0],
        "MISSING_033_REJECT": b"".join(line for line in raw.splitlines(keepends=True)
                                     if b"halo_info_033.asdf" not in line),
        "CRC_OUT_OF_RANGE_REJECT": raw.replace(b"3502514872", b"4294967296", 1),
        "ALTERED_ONE_FILE_BYTECOUNT_REJECT": raw.replace(b"2223483833", b"2223483834", 1),
        "FAKE_35TH_ASDF_REJECT": raw + b"123 100 halo_info_034.asdf\n",
    }
    passed = []
    for label, payload in cases.items():
        need(payload != raw, "negative control did not mutate source")
        try:
            audit(payload, prereg, enforce_pin=(label == "ONE_BYTE_SHA_PIN_REJECT"))
        except (ValueError, UnicodeError):
            passed.append(label)
        else:
            raise AssertionError("negative control unexpectedly accepted: " + label)
    return passed

def reports():
    prereg_bytes = PRE.read_bytes()
    need(gitblob(prereg_bytes) == PRE_BLOB, "prospective prereg Git blob drift")
    p = json.loads(prereg_bytes)
    need(all(p["absolute_STOP"].values()) and
         p["user_intake"]["provider_origin_independently_attested"] is False and
         p["user_intake"]["real_halo_asdf_supplied_or_read"] is False,
         "no-real-halo/observed-odd/independent-provider-attestation gate")
    raw = RAW.read_bytes()
    o, i, asdf_bytes = audit(raw, p)
    passed = negative_controls(raw, p)
    base = {
        "stage": p["stage"],
        "source": "user-uploaded Globus checksum manifest, NOT independent provider attestation",
        "prospective_prereg_git_blob": PRE_BLOB,
        "archived_user_manifest_git_blob": RAW_BLOB,
        "archived_user_manifest_sha256": sha(raw),
        "archived_user_manifest_bytes": len(raw),
        "source_pinned_portal_commit": p["source_snapshot"]["portal_commit"],
        "source_pinned_portal_table_git_blob": p["source_snapshot"]["published_table_git_blob"],
        "simname": p["source_snapshot"]["simulation"],
        "nominal_directory": p["source_snapshot"]["nominal_directory"],
        "actual_halo_ASDF_read": False,
        "independent_provider_manifest_attestation": False,
        "observed_odd_SEALED": True,
        "physical_M200c_a_or_neutrino_wake_or_B": "NOT_EVALUATED",
    }
    original = {
        **base,
        "audit": "original_strict_regex",
        "confirmed_user_supplied_asdf_record_count": len(o),
        "additional_user_supplied_checksum_manifest_entry_count": 1,
        "full_record_ledger": [
            {"filename": n, "GNU_POSIX_cksum_crc32": o[n][0], "byte_count": o[n][1]}
            for n in NAMES
        ],
        "asdf_byte_sum": asdf_bytes,
        "checksum_manifest_bytes": len(raw),
        "total_directory_bytes_with_manifest": asdf_bytes + len(raw),
        "prior_published_portal_directory_bytes_exact_match": True,
        "chosen_file": {
            "filename": p["locked_selected_file"],
            "GNU_POSIX_cksum_crc32": o[p["locked_selected_file"]][0],
            "byte_count": o[p["locked_selected_file"]][1],
            "real_file_bytes_NOT_obtained": True,
        },
    }
    independent = {
        **base,
        "audit": "independent_nonregex_tokenizer",
        "confirmed_user_supplied_asdf_record_count": len(i),
        "independent_record_ledger_sha256": sha(canonical([
            {"filename": name, "GNU_POSIX_cksum_crc32": crc, "byte_count": size}
            for name, crc, size in sorted(i)
        ])),
        "original_record_ledger_sha256": sha(canonical(original["full_record_ledger"])),
        "all_34_individual_records_identical_to_original": True,
        "independent_asdf_byte_sum": sum(size for _, _, size in i),
        "independent_total_directory_bytes_with_manifest": sum(size for _, _, size in i) + len(raw),
        "prior_published_portal_directory_bytes_exact_match": True,
        "synthetic_manifest_only_negative_controls_passed": passed,
    }
    need(independent["independent_record_ledger_sha256"] ==
         independent["original_record_ledger_sha256"],
         "independent ledger SHA differs")
    return original, independent

def main():
    original, independent = reports()
    expected = {
        "e17d2b3b4_original_user_checksum_ledger.json": canonical(original),
        "e17d2b3b4_independent_user_checksum_ledger.json": canonical(independent),
    }
    manifest = {
        "stage": original["stage"],
        "archive_kind": "immutable source-only user-supplied checksum audit, NOT provider attestation",
        "prereg_git_blob": PRE_BLOB,
        "source_user_checksum_git_blob": RAW_BLOB,
        "source_user_checksum_sha256": original["archived_user_manifest_sha256"],
        "archived_reports_sha256": {
            name: sha(data) for name, data in expected.items()
        },
        "prior_portal_aggregate_exact_match": True,
        "no_real_ASDF_download_or_read": True,
        "observed_odd_SEALED": True,
    }
    expected["archive_manifest.json"] = canonical(manifest)
    if "--emit" in sys.argv:
        ARCHIVE.mkdir(exist_ok=True)
        for name, data in expected.items():
            (ARCHIVE / name).write_bytes(data)
    else:
        for name, data in expected.items():
            target = ARCHIVE / name
            need(target.is_file() and target.read_bytes() == data,
                 "archived report missing or differs from original+independent replay: " + name)
    print("E17D2B3B4_SOURCE_ONLY_ORIGINAL_INDEPENDENT_AUDIT_OK", flush=True)
    print("ASDF_RECORDS=34 MANIFEST_BYTES=1386 TOTAL_DIRECTORY_BYTES=75722740486", flush=True)
    print("SELECTED halo_info_000.asdf CRC=3502514872 BYTES=2223483833", flush=True)
    print("NEGATIVE_CONTROLS=6 PASS REAL_ASDF_READ=NO PROVIDER_ATTESTATION=NO", flush=True)
    print("REPORT_SHA256", manifest["archived_reports_sha256"], flush=True)

if __name__ == "__main__":
    main()
