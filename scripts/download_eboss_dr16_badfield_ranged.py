#!/usr/bin/env python3
"""Resume-safe official DR16 large LRG MANGLE polygon SHA256 transfer.

Only the predeclared badfield polygon with SHA-pinned HTTP 206 prefix and
Content-Range total may be streamed. Chunks are atomically accepted only
after exact HTTP 206 range identity and byte-length checks. Incomplete
.part files are not mask products and must never be used in physics.
No galaxy FITS, odd data, mask membership or covariance enters.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from audit_eboss_dr16_official_polygon_sha import (
    PROTOCOL, expected_products, validate_url, digest_file, atomic_json,
)

ROOT=Path(__file__).resolve().parents[1]
PIN=ROOT/"source_data/eboss_dr16_badfield_range_probe_2026-09-25.json"
CHUNK=4*1024*1024
MAX_BYTES=512*1024*1024
PART_RETRIES=3
CONTENT_RANGE=re.compile(r"bytes (\d+)-(\d+)/(\d+)")


def validate_range(status, response_range, start, end, total, data, length_header=None):
    if status!=206:
        raise ValueError("Range download requires official HTTP 206")
    match=CONTENT_RANGE.fullmatch(response_range or "")
    if match is None or tuple(map(int,match.groups()))!=(start,end,total):
        raise ValueError("Official Content-Range differs from the pinned requested segment")
    if len(data)!=(end-start+1):
        raise ValueError("Incomplete official Range segment")
    if length_header is not None and int(length_header)!=len(data):
        raise ValueError("HTTP Content-Length differs from exact Range segment")
    return hashlib.sha256(data).hexdigest()


def transfer(inventory,protocol,pin,out_dir,timeout):
    if (pin.get("status")!="official_LRG_range_prefix_only"
        or pin.get("http_status")!=206
        or pin.get("filename")!="badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply"
        or pin.get("advertised_total_bytes")!=473778577
        or pin.get("prefix_bytes")!=65536
        or pin.get("range_chunk_bytes")!=CHUNK
        or pin.get("full_file_downloaded") is not False
        or pin.get("observed_odd_data_vector_read") is not False):
        raise ValueError("Large LRG source Range registration differs")
    products=[p for p in expected_products(protocol,inventory)
              if p["filename"]==pin["filename"]]
    if len(products)!=1 or products[0]["source_url"]!=pin["source_url"]:
        raise ValueError("Official source URL differs from immutable verified directory")
    src=products[0]
    validate_url(pin["source_url"],src["root"],pin["filename"])
    total=pin["advertised_total_bytes"]
    if total<65536 or total>MAX_BYTES:
        raise ValueError("Unexpected fixed total polygon file size")
    out_dir.mkdir(parents=True,exist_ok=True)
    path=out_dir/pin["filename"]
    partial=path.with_suffix(path.suffix+".part")
    report_path=out_dir/"official_badfield_ranged_sha_manifest.json"
    if path.exists():
        digest,nbytes=digest_file(path)
        if nbytes!=total:
            raise ValueError("Existing completed mask file size changed")
        previous=json.loads(report_path.read_text())
        if (previous.get("status")!="official_badfield_polygon_SHA256_pinned"
            or previous["sha256"]!=digest or previous["bytes"]!=nbytes):
            raise ValueError("Existing completed polygon lacks matching prior SHA provenance")
        print("EBOSS_LRG_BADFIELD_SHA_VERIFIED_CACHE",digest,total,flush=True)
        return previous
    current=partial.stat().st_size if partial.exists() else 0
    if current>total or current%CHUNK!=0:
        raise ValueError("Unverified local .part length (requires complete chunk boundary)")
    sha=hashlib.sha256()
    if current:
        with partial.open("rb") as f:
            first=f.read(65536)
            if hashlib.sha256(first).hexdigest()!=pin["prefix_sha256"]:
                raise ValueError("Existing partial file differs from pinned official prefix")
            sha.update(first)
            for block in iter(lambda:f.read(1024*1024),b""):
                sha.update(block)
        print("EBOSS_LRG_BADFIELD_RANGE_RESUME",current,total,flush=True)
    observed_prefix=(current>0)
    with partial.open("ab") as sink:
        for start in range(current,total,CHUNK):
            end=min(start+CHUNK,total)-1
            payload=None
            for attempt in range(1,PART_RETRIES+1):
                req=Request(pin["source_url"],headers={
                    "Range":f"bytes={start}-{end}",
                    "Accept-Encoding":"identity",
                    "User-Agent":"eBOSS-DR16-LRG-pinned-Range-SHA/1.0"})
                try:
                    with urlopen(req,timeout=timeout) as response:
                        validate_url(response.geturl(),src["root"],pin["filename"])
                        if response.status!=206:
                            raise ValueError("Official server ignored exact Range request")
                        response_range=response.headers.get("Content-Range")
                        response_len=response.headers.get("Content-Length")
                        payload=response.read(end-start+2)
                        validate_range(response.status,response_range,
                                       start,end,total,payload,response_len)
                    break
                except (OSError,HTTPError,URLError,TimeoutError) as exc:
                    if attempt==PART_RETRIES:
                        raise
                    print("EBOSS_LRG_BADFIELD_RANGE_RETRY",
                          start,end,attempt,type(exc).__name__,str(exc),
                          flush=True)
                    time.sleep(attempt*3)
            if payload is None:
                raise ValueError("Missing official Range payload")
            if not observed_prefix:
                if len(payload)<65536 or hashlib.sha256(payload[:65536]).hexdigest()!=pin["prefix_sha256"]:
                    raise ValueError("Official first 64KiB differs from pinned publication prefix")
                observed_prefix=True
            sink.write(payload)
            sink.flush()
            sha.update(payload)
            print("EBOSS_LRG_BADFIELD_RANGE_CHUNK_OK",end+1,total,flush=True)
    if partial.stat().st_size!=total or not observed_prefix:
        raise ValueError("Official polygon not downloaded completely")
    digest=sha.hexdigest()
    partial.replace(path)
    report={
        "status":"official_badfield_polygon_SHA256_pinned",
        "source_workflow_probe":pin["source_workflow"],
        "source_url":pin["source_url"],"filename":pin["filename"],
        "sha256":digest,"bytes":total,
        "prefix_sha256":pin["prefix_sha256"],
        "range_chunk_bytes":CHUNK,
        "verified_HTTP_206_Content_Range":True,
        "observed_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
        "actual_LRG_veto_mask_certified":False,
        "exact_cross_pair_window_validated":False,
    }
    atomic_json(report_path,report)
    print("EBOSS_LRG_BADFIELD_FULL_SHA256_PINNED",digest,total,flush=True)
    return report


def self_test():
    buf=b"a"*100
    assert validate_range(206,"bytes 0-99/200",0,99,200,buf,"100")
    for args in [
        (200,"bytes 0-99/200",0,99,200,buf,"100"),
        (206,"bytes 0-99/201",0,99,200,buf,"100"),
        (206,"bytes 0-99/200",0,99,200,buf[:99],"100"),
    ]:
        try:
            validate_range(*args)
        except ValueError:
            pass
        else:
            raise AssertionError("Malformed or truncated source segment accepted")
    pin=json.loads(PIN.read_text())
    assert pin["advertised_total_bytes"]==473778577
    assert pin["range_chunk_bytes"]==CHUNK
    assert pin["full_SHA256_certified"] is False
    print("EBOSS_LRG_BADFIELD_BOUNDED_RANGE_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory-json",type=Path)
    ap.add_argument("--out-dir",type=Path,
                    default=Path("eboss_workspace/large_lrg_badfield"))
    ap.add_argument("--timeout",type=float,default=110)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.inventory_json or args.timeout<=0:
        ap.error("--inventory-json and positive --timeout required")
    inventory=json.loads(args.inventory_json.read_text())
    protocol=json.loads(PROTOCOL.read_text())
    pin=json.loads(PIN.read_text())
    transfer(inventory,protocol,pin,args.out_dir,args.timeout)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
