#!/usr/bin/env python3
"""Official eBOSS DR16 LRG full-catalogue FITS HEADER ONLY sector inventory.

Reuse the already SHA-verified official DR16 directory index and pinned
11-polygon source manifest. Request only exact 2880-byte FITS header
blocks of one named full LRG catalogue; no FITS data rows, mask
membership, mock/data weights, randoms or odd-sector measurements.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

from audit_eboss_dr16_elg_full_header_only import read_header_only
from inspect_eboss_dr16_fits_headers import BLOCK, hdu_data_length

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_lrg_sector_header_protocol_2026-09-25.json"
OFFICIAL_INPUT = ROOT / "source_data/eboss_dr16_official_mask_source_protocol_2026-09-25.json"
INDEX_PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json"
PINNED_POLYGONS = ROOT / "source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json"
INDEX_REPORT = ROOT / "eboss_workspace/official_mask_inventory/bit8_reference_index.json"
RANGE = re.compile(r"^bytes ([0-9]+)-([0-9]+)/([0-9]+)$")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def validate_pinned_inputs(p, index):
    official, ip, polygon = (read_json(OFFICIAL_INPUT),
                             read_json(INDEX_PROTOCOL),
                             read_json(PINNED_POLYGONS))
    base = official["source_data_roots"]["nersc_mirror"]
    if (p["official_url"] != base + p["file_name"]
            or p["official_url"] !=
               "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/eBOSS_LRG_full_ALLdata-vDR16.fits"
            or p["source_index_sha256"] != ip["root_listing_sha256"]
            or p["source_index_sha256"] != official["source_index_sha256"]["dr16_root"]
            or index.get("status") != "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED"
            or index.get("root_index_sha256") != p["source_index_sha256"]
            or p["published_lrg_sector_thresholds_DESCRIPTIVE_NOT_NEW_CUT"]
               != {"C_eBOSS": "> 0.5", "C_z": "> 0.5"}
            or p["max_total_header_blocks"] != 90 or p["max_requests"] != 90
            or p["new_selection_rule_applied"] is not False
            or p["observed_odd_data_vector_read"] is not False):
        raise ValueError("Pinned LRG source/index/sector protocol mismatch")
    candidates = [
        v for v in index["candidate_files"]
        if v["filename"] == p["file_name"]
        and v["source_url"] == p["official_url"]
        and not v["is_directory"]
        and v["classification"] == "unverified_file_name_only"
    ]
    if len(candidates) != 1:
        raise ValueError("Exact LRG full file absent from already SHA-pinned release index")
    products = {x["filename"]: x for x in polygon["products"]}
    if len(products) != 11 or polygon["status"] != "official_lrg_elg_polygon_bytes_pinned_only":
        raise ValueError("Prior official 11-polygon manifest not complete")
    expected_lrg = official["lrg_published_polygons"]
    if (set(expected_lrg) !=
            set(p["exact_lrg_candidate_polygon_files_UNMAPPED"])
            | {p["qso_only_veto_DO_NOT_APPLY_TO_LRG"]}
            or official["lrg_root_geometry_polygon"] != p["official_footprint_file"]):
        raise ValueError("Published LRG/QSO polygon filenames are not as predeclared")
    for name in expected_lrg + [p["official_footprint_file"]]:
        record = products.get(name)
        if (not record or len(record["sha256"]) != 64
                or record["bytes"] <= 0
                or record["role"] != (
                    "lrg_noveto_footprint_polygon" if name == p["official_footprint_file"]
                    else "lrg_published_polygon"
                )):
            raise ValueError("Missing/malformed published polygon fingerprint: " + name)
    return {name: products[name]["sha256"]
            for name in expected_lrg + [p["official_footprint_file"]]}


def inspect(url, p, *, timeout, block_reader=None):
    n = 0
    meta_first = None
    if block_reader is None:
        def block_reader(offset):
            req = Request(url, headers={
                "Range": f"bytes={offset}-{offset + BLOCK - 1}",
                "Accept-Encoding": "identity",
                "User-Agent": "blinded-eBOSS-DR16-LRG-sector-header/1.0",
            })
            with urlopen(req, timeout=timeout) as response:
                if response.status != 206 or response.geturl() != url:
                    raise ValueError("Exact official LRG FITS HTTP 206 range not honored")
                if response.headers.get("Content-Encoding", "identity").lower() != "identity":
                    raise ValueError("Official LRG FITS response content-encoded")
                m = RANGE.fullmatch(response.headers.get("Content-Range", ""))
                if not m:
                    raise ValueError("Missing official LRG FITS Content-Range")
                first, last, length = map(int, m.groups())
                if first != offset or last != offset + BLOCK - 1 or length < last + 1:
                    raise ValueError("Official LRG FITS Content-Range does not match request")
                data = response.read(BLOCK + 1)
                return data, {
                    "total": length,
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                }

    def one_block(offset):
        nonlocal n, meta_first
        if n >= p["max_requests"]:
            raise ValueError("LRG header exceeds frozen HTTP range request bound")
        block, meta = block_reader(offset)
        if len(block) != BLOCK or not isinstance(meta.get("total"), int):
            raise ValueError("Non-exact LRG FITS header range")
        if meta_first is None:
            meta_first = meta
        elif meta != meta_first:
            raise ValueError("LRG FITS remote length/ETag/Last-Modified changed mid-header")
        if offset + BLOCK > meta_first["total"]:
            raise ValueError("LRG header request beyond official file size")
        n += 1
        return block

    primary, raw_first = read_header_only(one_block, 0, p["max_total_header_blocks"])
    if primary.get("SIMPLE") is not True or hdu_data_length(primary) != 0:
        raise ValueError("LRG primary FITS HDU unexpected or contains data")
    remaining = p["max_total_header_blocks"] - n
    if remaining < 1:
        raise ValueError("LRG BINTABLE header exceeds frozen block bound")
    table, raw_table = read_header_only(one_block, len(raw_first), remaining)
    if str(table.get("XTENSION", "")).upper() != "BINTABLE":
        raise ValueError("LRG first FITS extension is not BINTABLE")
    fields = int(table["TFIELDS"])
    if fields <= 0 or fields > 2048:
        raise ValueError("Malformed official LRG BINTABLE TFIELDS")
    names = [str(table.get(f"TTYPE{i}", "")).strip()
             for i in range(1, fields + 1)]
    formats = [str(table.get(f"TFORM{i}", "")).strip()
               for i in range(1, fields + 1)]
    if (any(not name or not form for name, form in zip(names, formats))
            or not set(p["required_header_columns"]).issubset(set(names))):
        raise ValueError("Required LRG RA/DEC or BINTABLE structure missing")
    interesting = [
        {"name": name, "format": form}
        for name, form in zip(names, formats)
        if any(token in name.upper() for token in p["interesting_header_tokens"])
    ]
    return {
        "status": "OFFICIAL_LRG_FULL_HEADER_PARSED_SECTOR_SEMANTICS_NOT_YET_CERTIFIED",
        "official_url": url,
        "remote_bytes_reported": meta_first["total"],
        "etag": meta_first["etag"],
        "last_modified": meta_first["last_modified"],
        "header_bytes_requested": len(raw_first + raw_table),
        "header_blocks_requested": n,
        "header_range_sha256": hashlib.sha256(raw_first + raw_table).hexdigest(),
        "full_official_fits_sha256_verified": False,
        "declared_bintable_rows_NOT_READ": int(table["NAXIS2"]),
        "declared_bintable_row_bytes": int(table["NAXIS1"]),
        "column_count": fields,
        "column_names": names,
        "interesting_sector_mask_columns": interesting,
        "published_lrg_sector_thresholds_DESCRIPTIVE_NOT_APPLIED":
            p["published_lrg_sector_thresholds_DESCRIPTIVE_NOT_NEW_CUT"],
        "published_lrg_veto_classes_NOT_YET_MAPPED_TO_FILENAMES":
            p["published_lrg_veto_classes_NOT_YET_MAPPED_TO_EXACT_FILES"],
        "qso_only_veto_NOT_APPLIED_TO_LRG": p["qso_only_veto_DO_NOT_APPLY_TO_LRG"],
        "physical_lrg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "observed_galaxy_data_rows_read": False,
        "observed_random_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_selection_rule_applied": False,
        "note": (
            "Full LRG FITS header and previously pinned MANGLE filenames only. "
            "Header column names do not identify exact sector completeness "
            "mapping, veto polygon composition or physical mask membership. "
            "No data-row byte was requested."
        ),
        "errors": [],
    }


def self_test(p):
    def card(key, val):
        return (f"{key:<8}= {val}".ljust(80)).encode("ascii")

    def header(cards):
        a = b"".join(cards + [b"END".ljust(80, b" ")])
        return a.ljust((len(a) + BLOCK - 1) // BLOCK * BLOCK, b" ")

    a = header([card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")])
    b = header([
        card("XTENSION", "'BINTABLE'"), card("BITPIX", "8"),
        card("NAXIS", "2"), card("NAXIS1", "36"), card("NAXIS2", "2"),
        card("PCOUNT", "0"), card("GCOUNT", "1"),
        card("TFIELDS", "4"), card("TTYPE1", "'RA'"),
        card("TFORM1", "'1D'"), card("TTYPE2", "'DEC'"),
        card("TFORM2", "'1D'"), card("TTYPE3", "'sector'"),
        card("TFORM3", "'1J'"), card("TTYPE4", "'COMP_BOSS'"),
        card("TFORM4", "'1D'"),
    ])
    data = a + b + b"X" * 72
    requested = []

    def fake_reader(offset):
        requested.append(offset)
        if offset + BLOCK > len(a + b):
            raise AssertionError("Attempt to request synthetic LRG table data row")
        return data[offset:offset + BLOCK], {
            "total": len(data), "etag": '"synthetic"', "last_modified": None,
        }

    report = inspect(p["official_url"], p, timeout=1, block_reader=fake_reader)
    assert requested == [0, BLOCK]
    assert report["column_count"] == 4
    assert report["observed_galaxy_data_rows_read"] is False
    assert any(x["name"] == "COMP_BOSS"
               for x in report["interesting_sector_mask_columns"])
    validate_pinned_inputs(p, {
        "status": "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED",
        "root_index_sha256": p["source_index_sha256"],
        "candidate_files": [{
            "filename": p["file_name"], "source_url": p["official_url"],
            "is_directory": False, "classification": "unverified_file_name_only",
        }],
    })
    print("EBOSS_LRG_FULL_HEADER_SECTOR_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--index-report", type=Path, default=INDEX_REPORT)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--timeout", type=float, default=40)
    args = ap.parse_args()
    p = read_json(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    out = args.out or ROOT / p["output"]
    try:
        if args.timeout <= 0:
            raise ValueError("Timeout must be positive")
        idx = read_json(args.index_report)
        veto_manifest = validate_pinned_inputs(p, idx)
        r = inspect(p["official_url"], p, timeout=args.timeout)
        r["pinned_official_lrg_polygon_sha256"] = veto_manifest
    except Exception as e:
        r = {
            "status": "OFFICIAL_LRG_FULL_HEADER_INCOMPLETE_STOP",
            "errors": [str(e)],
            "observed_galaxy_data_rows_read": False,
            "observed_odd_data_vector_read": False,
            "physical_lrg_mask_certified": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(".tmp.json")
    temp.write_text(json.dumps(r, indent=2) + "\n", encoding="utf-8")
    temp.replace(out)
    print("EBOSS_LRG_FULL_SECTOR_HEADER", r["status"], flush=True)
    print("REPORT", out, flush=True)
    if r.get("errors"):
        print("ERRORS", *r["errors"], sep="\n", flush=True)
        return 2
    print("COLUMN_COUNT", r["column_count"], flush=True)
    print("SECTOR_AND_MASK_COLUMNS",
          json.dumps(r["interesting_sector_mask_columns"]), flush=True)
    print("OBSERVED_ODD_READ", r["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
