#!/usr/bin/env python3
"""Acquire and SHA-pin official DR16 LRG full FITS raw bytes; NO row decoding.

The previous official LRG header-only report was uploaded and SHA-pinned.
Validate exact report bytes and header/official-index/polygon provenance
before streaming a quarantined 196-MB official FITS file. Reuse the
previously tested ELG raw-byte stream code (not its ELG-specific contract).
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

from audit_eboss_dr16_elg_full_bytes_provenance import (
    acquire, digest_stream, hash_existing,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_full_bytes_provenance_protocol_2026-09-25.json"
INDEX_PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json"
LRG_HEADER_PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_sector_header_protocol_2026-09-25.json"
POLYGON_MANIFEST = ROOT / "source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json"
INDEX_REPORT = ROOT / "eboss_workspace/official_mask_inventory/bit8_reference_index.json"


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def verify_previous(p, index, header, header_raw_sha):
    index_p = load(INDEX_PROTOCOL)
    header_p = load(LRG_HEADER_PROTOCOL)
    manifest = load(POLYGON_MANIFEST)
    if (
        p["official_source_url"] != header_p["official_url"]
        or p["official_source_url"] != index_p["official_root"] + header_p["file_name"]
        or p["expected_file_bytes"] != p["byte_cap"]
        or p["expected_file_bytes"] != 196485120
        or p["expected_header_bytes"] != 20160
        or p["expected_header_sha256"] !=
           "a48b0d1055b5f7638dc77a86eb5ff5c861b629d19ff75cd8ae8122caa7e9014b"
        or p["verified_user_uploaded_lrg_header_report"]["sha256"] !=
           "901c5486e817da57a162d9bec74171a0fbc09f4e0901d72b87d7f4495ff51308"
        or header_raw_sha != p["verified_user_uploaded_lrg_header_report"]["sha256"]
        or index.get("status") != "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED"
        or index.get("root_index_sha256") != index_p["root_listing_sha256"]
        or header.get("status") !=
           "OFFICIAL_LRG_FULL_HEADER_PARSED_SECTOR_SEMANTICS_NOT_YET_CERTIFIED"
        or header.get("official_url") != p["official_source_url"]
        or header.get("remote_bytes_reported") != p["expected_file_bytes"]
        or header.get("etag") != p["expected_etag"]
        or header.get("last_modified") != p["expected_last_modified"]
        or header.get("header_bytes_requested") != p["expected_header_bytes"]
        or header.get("header_range_sha256") != p["expected_header_sha256"]
        or header.get("declared_bintable_rows_NOT_READ") !=
           p["expected_bintable_header_rows_NOT_READ"]
        or header.get("declared_bintable_row_bytes") !=
           p["expected_bintable_row_bytes"]
        or header.get("column_count") != p["expected_column_count"]
        or header.get("observed_galaxy_data_rows_read") is not False
        or header.get("observed_odd_data_vector_read") is not False
        or p["observed_odd_data_vector_read"] is not False
        or p["observed_galaxy_rows_parsed"] is not False
        or p["full_source_sha256_NOT_YET_KNOWN"] is not True
        or manifest["status"] != "official_lrg_elg_polygon_bytes_pinned_only"
    ):
        raise ValueError("Prior uploaded header / official source / byte-only protocol mismatch")
    entries = [
        x for x in index["candidate_files"]
        if x.get("filename") == header_p["file_name"]
        and x.get("source_url") == p["official_source_url"]
        and x.get("is_directory") is False
        and x.get("classification") == "unverified_file_name_only"
    ]
    if len(entries) != 1:
        raise ValueError("Official LRG full catalogue absent from previously pinned index")
    cols = {
        x["name"]: x["format"]
        for x in header["interesting_sector_mask_columns"]
    }
    if any(cols.get(k) != v
           for k, v in p["expected_sector_header_columns"].items()):
        raise ValueError("LRG sector/completeness column types differ from uploaded header")
    expected_polygon_sha = {
        r["filename"]: r["sha256"]
        for r in manifest["products"]
        if r["role"] in ("lrg_published_polygon", "lrg_noveto_footprint_polygon")
    }
    if (
        len(expected_polygon_sha) != 8
        or header.get("pinned_official_lrg_polygon_sha256") != expected_polygon_sha
    ):
        raise ValueError("LRG official pinned polygon manifest does not match uploaded report")


def self_test(p):
    assert p["expected_file_bytes"] == p["byte_cap"] == 196485120
    assert p["expected_header_bytes"] == 20160
    assert p["observed_odd_data_vector_read"] is False
    raw = b"synthetic LRG FITS header" + b"synthetic data never parsed"
    n_header = len(b"synthetic LRG FITS header")
    out = io.BytesIO()
    digest, n = digest_stream(
        io.BytesIO(raw), out, expected_bytes=len(raw),
        header_bytes=n_header,
        header_sha=hashlib.sha256(raw[:n_header]).hexdigest(),
    )
    assert out.getvalue() == raw
    assert digest == hashlib.sha256(raw).hexdigest() and n == len(raw)
    for invalid in (raw[:-1], raw + b"!", b"X" + raw[1:]):
        try:
            digest_stream(
                io.BytesIO(invalid), io.BytesIO(), expected_bytes=len(raw),
                header_bytes=n_header,
                header_sha=hashlib.sha256(raw[:n_header]).hexdigest(),
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Bad synthetic byte stream was accepted")
    print("EBOSS_LRG_FULL_RAW_PROVENANCE_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--existing-file", type=Path,
                    help="Hash an already downloaded exact official LRG FITS file")
    ap.add_argument("--timeout", type=float, default=120)
    args = ap.parse_args()
    p = load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0

    out = ROOT / p["report_output"]
    try:
        if args.timeout <= 0:
            raise ValueError("Timeout must be positive")
        path_report = ROOT / p["verified_user_uploaded_lrg_header_report"]["local_path"]
        header_bytes = path_report.read_bytes()
        header_sha = hashlib.sha256(header_bytes).hexdigest()
        if len(header_bytes) != p["verified_user_uploaded_lrg_header_report"]["bytes"]:
            raise ValueError("Previously uploaded LRG header JSON byte count changed")
        verify_previous(p, load(INDEX_REPORT), json.loads(header_bytes), header_sha)
        path = args.existing_file or ROOT / p["local_output"]
        if args.existing_file:
            digest, nbytes = hash_existing(path, p)
            mode = "verified_existing_local_bytes"
        else:
            digest, nbytes, mode = acquire(p, path, args.timeout)
        if out.is_file():
            old = load(out)
            if (
                old.get("status") != "OFFICIAL_LRG_FULL_BYTES_SHA_PINNED_ONLY"
                or old.get("full_file_sha256") != digest
                or old.get("file_bytes") != nbytes
            ):
                raise ValueError("Previous full LRG checksum disagrees: no silent repin")
        result = {
            "status": "OFFICIAL_LRG_FULL_BYTES_SHA_PINNED_ONLY",
            "official_source_url": p["official_source_url"],
            "local_file_path": str(path.resolve()),
            "acquisition_mode": mode,
            "file_bytes": nbytes,
            "full_file_sha256": digest,
            "verified_parent_header_report_sha256": header_sha,
            "verified_first_header_sha256": p["expected_header_sha256"],
            "full_source_sha256_frozen_in_repository": False,
            "published_completeness_selection_applied": False,
            "official_lrg_sector_semantics_certified": False,
            "official_lrg_polygon_composition_certified": False,
            "observed_galaxy_rows_parsed": False,
            "observed_sector_columns_parsed": False,
            "observed_random_or_mock_rows_parsed": False,
            "observed_odd_data_vector_read": False,
            "physical_lrg_elg_pair_window_certified": False,
            "note": (
                "Complete official LRG FITS bytes SHA-hashed only. No FITS table "
                "column decoded. Freeze this first-seen SHA256 in a separate "
                "pre-row inspection protocol; it is not an SDSS-published checksum."
            ),
            "errors": [],
        }
    except Exception as exc:
        result = {
            "status": "OFFICIAL_LRG_FULL_BYTES_PROVENANCE_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "observed_galaxy_rows_parsed": False,
            "observed_odd_data_vector_read": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    tmp.replace(out)
    print("EBOSS_LRG_FULL_BYTES_PROVENANCE", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep="\n", flush=True)
        return 2
    print("BYTES", result["file_bytes"], flush=True)
    print("SHA256", result["full_file_sha256"], flush=True)
    print("OBSERVED_ODD_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
