#!/usr/bin/env python3
"""Bounded HTTP Range probe for two slow official eBOSS DR16 LRG polygons.

Records the Content-Range file length and a 64-KiB prefix SHA256;
does not download or mark either full polygon as SHA-certified.
No data catalogues, sky masks, odd measurements or inference are read.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

from audit_eboss_dr16_official_polygon_sha import (
    PROTOCOL, expected_products, validate_url, atomic_json,
)

LARGE_FILES = (
    "allsky_bright_star_mask_pix.ply",
    "badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply",
)
PROBE = 65536
RESPONSE = re.compile(r"bytes 0-(\d+)/(\d+)")


def check_range_header(status, content_range, payload):
    if status != 206:
        raise ValueError("Official source does not return HTTP 206 for bounded Range probe")
    match = RESPONSE.fullmatch(content_range or "")
    if match is None:
        raise ValueError("Official Content-Range is missing or malformed")
    high, total = map(int, match.groups())
    if total < 16 or high != min(total, PROBE) - 1:
        raise ValueError("Official Content-Range does not cover precisely the requested prefix")
    if len(payload) != high + 1:
        raise ValueError("Official source prefix is truncated")
    if not payload.startswith(b"#") and b"polygon" not in payload[:64].lower():
        raise ValueError("Published MANGLE polygon prefix is unexpected")
    return total


def probe(inventory, protocol, filename, timeout):
    matches = [
        p for p in expected_products(protocol,inventory)
        if p["filename"] == filename
    ]
    if len(matches) != 1 or filename not in LARGE_FILES:
        raise ValueError("Range probe not in fixed two-file cohort")
    product = matches[0]
    url = product["source_url"]
    validate_url(url, product["root"], filename)
    req = Request(url, headers={
        "Range": f"bytes=0-{PROBE-1}",
        "Accept-Encoding": "identity",
        "User-Agent":"eBOSS-DR16-official-LRG-range-size-probe/1.0"
    })
    with urlopen(req, timeout=timeout) as response:
        validate_url(response.geturl(), product["root"], filename)
        if response.status != 206:
            raise ValueError("Official server ignores bounded HTTP Range")
        data = response.read(PROBE + 1)
        length = check_range_header(
            response.status,response.headers.get("Content-Range"),data)
    return {
        "status":"official_LRG_range_prefix_only",
        "filename":filename,"source_url":url,
        "http_status":206,"prefix_bytes":len(data),
        "prefix_sha256":hashlib.sha256(data).hexdigest(),
        "advertised_total_bytes":length,
        "within_current_512MiB_limit":length <= 512*1024*1024,
        "full_file_downloaded":False,
        "full_sha256_certified":False,
        "observed_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
        "official_LRG_mask_applied":False,
    }


def self_test():
    test=b"# MANGLE polygons\n" + b"x"*(PROBE-len(b"# MANGLE polygons\n"))
    assert check_range_header(206,f"bytes 0-{PROBE-1}/{PROBE*2}",
                              test) == PROBE*2
    for code,header,blob in (
        (200,f"bytes 0-{PROBE-1}/{PROBE*2}",test),
        (206,None,test),
        (206,f"bytes 0-{PROBE-1}/{PROBE*2}",test[:-1]),
    ):
        try:
            check_range_header(code,header,blob)
        except ValueError:
            pass
        else:
            raise AssertionError("An invalid source Range response was accepted")
    print("EBOSS_OFFICIAL_LRG_RANGE_PROBE_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory-json",type=Path)
    ap.add_argument("--filename",choices=LARGE_FILES)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--timeout",type=float,default=100)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.inventory_json or not args.filename or not args.out:
        ap.error("--inventory-json, --filename and --out are required")
    if args.timeout <=0:
        ap.error("timeout must be positive")
    inv=json.loads(args.inventory_json.read_text())
    prot=json.loads(PROTOCOL.read_text())
    try:
        r=probe(inv,prot,args.filename,args.timeout)
    except Exception as exc:
        r={
            "status":"official_LRG_range_probe_failed",
            "filename":args.filename,"errors":[str(exc)],
            "full_file_downloaded":False,"full_sha256_certified":False,
            "observed_odd_data_vector_read":False,
        }
    atomic_json(args.out,r)
    print("EBOSS_OFFICIAL_LRG_RANGE_PROBE",r["status"],args.filename,
          r.get("advertised_total_bytes"),flush=True)
    return 2 if "errors" in r else 0


if __name__=="__main__":
    raise SystemExit(main())
