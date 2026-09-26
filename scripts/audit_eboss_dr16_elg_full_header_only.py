#!/usr/bin/env python3
"""Inspect ONLY the FITS header of pinned DR16 ELG full_ALLdata candidate.

Every HTTP request asks for exactly ONE 2880-byte FITS header block.
Stop at the primary and first BINTABLE END cards; never request data rows.
Column names are structural metadata, NOT bit-8 production verification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from inspect_eboss_dr16_fits_headers import BLOCK, hdu_data_length, parse_header

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_full_header_only_protocol_2026-09-25.json"
INDEX_PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json"
DEFAULT_INDEX = Path("eboss_workspace/official_mask_inventory/bit8_reference_index.json")
DEFAULT_OUT = Path("eboss_workspace/official_mask_inventory/elg_full_header_only.json")
CONTENT_RANGE = re.compile(r"^bytes ([0-9]+)-([0-9]+)/([0-9]+)$")


def frozen(protocol, index_protocol, index):
    root = protocol["source_root"]
    file_name = protocol["source_file"]
    if (
        root != index_protocol["official_root"]
        or root != "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/"
        or file_name != "eBOSS_ELG_full_ALLdata-vDR16.fits"
        or protocol["max_total_header_blocks"] != 90
        or protocol["max_requests"] != 90
        or protocol["observed_odd_data_vector_read"] is not False
        or index["status"] != "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED"
        or index["root_index_sha256"] != index_protocol["root_listing_sha256"]
    ):
        raise ValueError("Pinned official index / header-only protocol mismatch")
    named = [
        item for item in index["candidate_files"]
        if item["filename"] == file_name
        and item["source_url"] == root + file_name
        and item["classification"] == "unverified_file_name_only"
        and item["is_directory"] is False
    ]
    if len(named) != 1:
        raise ValueError("Exact ELG full file absent from verified official index")
    return root + file_name


def read_header_only(fetch_block, offset, remaining_blocks):
    blocks = []
    for i in range(remaining_blocks):
        block = fetch_block(offset + i * BLOCK)
        if len(block) != BLOCK:
            raise ValueError("Incomplete requested FITS header block")
        if not all(b < 128 for b in block):
            raise ValueError("Non-ASCII bytes in FITS header block")
        blocks.append(block)
        end = any(
            block[j:j + 8].decode("ascii").strip() == "END"
            for j in range(0, BLOCK, 80)
        )
        if end:
            raw = b"".join(blocks)
            header, next_pos = parse_header(raw, 0)
            if next_pos != len(raw):
                raise ValueError("FITS END card did not terminate parsed header")
            return header, raw
    raise ValueError("FITS header exceeds frozen bounded block count")


def inspect(url, protocol, timeout, block_reader=None):
    requested_offsets = []
    headers = []
    remote_size = None
    first_headers = None
    if block_reader is None:
        def block_reader(offset):
            request = Request(
                url,
                headers={
                    "Range": f"bytes={offset}-{offset + BLOCK - 1}",
                    "Accept-Encoding": "identity",
                    "User-Agent": "blinded-eBOSS-ELG-header-only/1.0",
                },
            )
            with urlopen(request, timeout=timeout) as response:
                if response.status != 206:
                    raise ValueError("Server did not honor exact HTTP 206 range")
                if response.geturl() != url:
                    raise ValueError("Range request redirected away from exact official product")
                match = CONTENT_RANGE.fullmatch(
                    response.headers.get("Content-Range", "")
                )
                if not match:
                    raise ValueError("Missing or invalid Content-Range")
                begin, finish, total = map(int, match.groups())
                if (begin != offset or finish != offset + BLOCK - 1
                        or total < finish + 1):
                    raise ValueError("Content-Range differs from exact requested block")
                raw = response.read(BLOCK + 1)
                meta = {
                    "total": total, "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                }
                return raw, meta

    def get_block(offset):
        nonlocal remote_size, first_headers
        if len(requested_offsets) >= protocol["max_requests"]:
            raise ValueError("Frozen maximum HTTP range requests reached")
        raw, meta = block_reader(offset)
        if len(raw) != BLOCK:
            raise ValueError("Partial/overlong FITS range response")
        if remote_size is None:
            remote_size = meta["total"]
            first_headers = meta
        elif remote_size != meta["total"]:
            raise ValueError("Remote total size changed between header ranges")
        if offset + BLOCK > remote_size:
            raise ValueError("Attempt to read past reported official file size")
        requested_offsets.append(offset)
        return raw

    primary, raw_primary = read_header_only(
        get_block, 0, protocol["max_total_header_blocks"]
    )
    if primary.get("SIMPLE") is not True:
        raise ValueError("Missing FITS SIMPLE primary header")
    if hdu_data_length(primary) != 0:
        raise ValueError("Primary HDU contains data: cannot advance header-only")
    extension_offset = len(raw_primary)
    remaining = protocol["max_total_header_blocks"] - len(requested_offsets)
    if remaining < 1:
        raise ValueError("No remaining header blocks")
    extension, raw_extension = read_header_only(
        get_block, extension_offset, remaining
    )
    if str(extension.get("XTENSION", "")).upper() != protocol["expected_first_extension"]:
        raise ValueError("First extension is not expected BINTABLE")
    nfields = int(extension["TFIELDS"])
    if nfields <= 0 or nfields > 2048:
        raise ValueError("Invalid BINTABLE TFIELDS")
    names = [str(extension.get(f"TTYPE{i}", "")).strip()
             for i in range(1, nfields + 1)]
    forms = [str(extension.get(f"TFORM{i}", "")).strip()
             for i in range(1, nfields + 1)]
    if any(not x for x in names) or any(not x for x in forms):
        raise ValueError("Missing published column TTYPE/TFORM")
    if not set(protocol["must_include_columns"]).issubset(
        set(name.upper() for name in names)
    ):
        raise ValueError("Required RA/DEC columns missing from candidate")
    relevant = [
        {"name": name, "format": form}
        for name, form in zip(names, forms)
        if any(token in name.upper() for token in protocol["report_columns_matching_any"])
    ]
    header_end = extension_offset + len(raw_extension)
    if requested_offsets != list(range(0, header_end, BLOCK)):
        raise ValueError("Header-range offsets not contiguous from byte zero")
    return {
        "status": "ELG_FULL_HEADER_ONLY_PARSED_REFERENCE_NOT_VERIFIED",
        "source_url": url,
        "remote_bytes_reported": remote_size,
        "remote_etag": first_headers["etag"],
        "remote_last_modified": first_headers["last_modified"],
        "header_bytes_requested": header_end,
        "header_blocks_requested": len(requested_offsets),
        "header_range_sha256": hashlib.sha256(
            raw_primary + raw_extension
        ).hexdigest(),
        "full_source_file_sha256_verified": False,
        "primary_hdu_data_bytes": 0,
        "first_extension": "BINTABLE",
        "header_declared_rows_NOT_READ": int(extension["NAXIS2"]),
        "header_declared_row_bytes": int(extension["NAXIS1"]),
        "column_count": nfields,
        "column_names": names,
        "column_formats": forms,
        "candidate_mask_and_id_columns": relevant,
        "has_mskbit_column": "MSKBIT" in (x.upper() for x in names),
        "official_pre_veto_bit8_reference_verified": False,
        "official_bit8_production_convention_resolved": False,
        "official_bit8_applied": False,
        "physical_mask_certified": False,
        "observed_galaxy_data_rows_read": False,
        "mock_galaxy_data_rows_read": False,
        "observed_random_positions_read": False,
        "mock_random_positions_read": False,
        "observed_odd_data_vector_read": False,
        "note": (
            "The release header contains MSKBIT, but no table rows were read. "
            "Bit8-positive-row retention and the official production mapping "
            "are not inferable from header fields alone. No data-row byte was requested."
        ),
    }


def self_test():
    def card(name, value):
        return (f"{name:<8}= {value}".ljust(80)).encode("ascii")

    def head(cards):
        raw = b"".join(cards + [b"END".ljust(80, b" ")])
        return raw.ljust((len(raw) + BLOCK - 1) // BLOCK * BLOCK, b" ")

    primary = head([
        card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0"),
    ])
    table = head([
        card("XTENSION", "'BINTABLE'"), card("BITPIX", "8"),
        card("NAXIS", "2"), card("NAXIS1", "24"), card("NAXIS2", "7"),
        card("PCOUNT", "0"), card("GCOUNT", "1"), card("TFIELDS", "3"),
        card("TTYPE1", "'RA'"), card("TFORM1", "'1D'"),
        card("TTYPE2", "'DEC'"), card("TFORM2", "'1D'"),
        card("TTYPE3", "'MSKBIT'"), card("TFORM3", "'1J'"),
    ])
    fake_file = primary + table + b"R" * (7 * 24)
    offsets = []

    def fake_reader(offset):
        offsets.append(offset)
        if offset + BLOCK > len(primary + table):
            raise AssertionError("Test tried to download even one data row")
        return fake_file[offset:offset + BLOCK], {
            "total": len(fake_file), "etag": '"synthetic"',
            "last_modified": None,
        }

    p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    ip = json.loads(INDEX_PROTOCOL.read_text(encoding="utf-8"))
    index = {
        "status": "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED",
        "root_index_sha256": ip["root_listing_sha256"],
        "candidate_files": [{
            "filename": p["source_file"],
            "source_url": p["source_root"] + p["source_file"],
            "classification": "unverified_file_name_only",
            "is_directory": False,
        }],
    }
    url = frozen(p, ip, index)
    r = inspect(url, p, 1, block_reader=fake_reader)
    assert offsets == [0, BLOCK]
    assert r["has_mskbit_column"] and r["header_declared_rows_NOT_READ"] == 7
    assert r["observed_galaxy_data_rows_read"] is False
    print("EBOSS_ELG_FULL_HEADER_ONLY_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--index-report", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--timeout", type=float, default=30)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    try:
        if args.timeout <= 0:
            raise ValueError("Timeout must be positive")
        p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
        ip = json.loads(INDEX_PROTOCOL.read_text(encoding="utf-8"))
        idx = json.loads(args.index_report.read_text(encoding="utf-8"))
        source_url = frozen(p, ip, idx)
        report = inspect(source_url, p, args.timeout)
    except Exception as exc:
        report = {
            "status": "ELG_FULL_HEADER_ONLY_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "official_pre_veto_bit8_reference_verified": False,
            "observed_galaxy_data_rows_read": False,
            "observed_odd_data_vector_read": False,
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    tmp.replace(args.out)
    print("EBOSS_ELG_FULL_HEADER_ONLY", report["status"], flush=True)
    print("REPORT", args.out, flush=True)
    if report.get("errors"):
        print("ERRORS", *report["errors"], sep="\n", flush=True)
        return 2
    print("COLUMN_COUNT", report["column_count"], flush=True)
    print("HAS_MSKBIT", report["has_mskbit_column"], flush=True)
    print("CANDIDATE_COLUMNS", report["candidate_mask_and_id_columns"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
