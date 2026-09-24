#!/usr/bin/env python3
"""Inspect bounded gzip/FITS headers of matched public eBOSS realistic EZmocks.

A stratified, preselected sample of realization IDs is used. Only compressed
prefixes and FITS headers are inspected; galaxy and random rows, pair counts,
survey windows, and odd-sector statistics are not read or computed.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from urllib.request import Request, urlopen
import zlib

from inspect_eboss_dr16_fits_headers import metadata

BASE = "https://data.sdss.org/sas/dr17/eboss/lss/EZmocks/v1_0_0/realistic/"
TRACERS = ("eBOSS_LRG", "eBOSS_ELG")
CAPS = ("NGC", "SGC")
ROLES = ("dat", "ran")
DEFAULT_IDS = (1, 500, 1000)
MAX_COMPRESSED = 512 * 1024
MAX_DECOMPRESSED = 8 * 1024 * 1024


def mock_path(tracer: str, cap: str, role: str, realization: int) -> str:
    if tracer not in TRACERS or cap not in CAPS or role not in ROLES:
        raise ValueError("Unsupported tracer, cap or role")
    if not 1 <= realization <= 1000:
        raise ValueError("Mock ID must be in the published 0001-1000 range")
    filename = (
        f"EZmock_realistic_{tracer}_{cap}_v7_{realization:04d}."
        f"{role}.fits.gz"
    )
    return f"{tracer}/{role}/{filename}"


def inspect_prefix(url: str, timeout: float, max_compressed: int,
                   max_decompressed: int) -> dict:
    base, request_url = urlsplit(BASE), urlsplit(url)
    if (request_url.scheme != "https" or request_url.netloc != base.netloc
            or not request_url.path.startswith(base.path)
            or request_url.query or request_url.fragment):
        raise ValueError("Mock URL is outside the approved SDSS release")
    request = Request(url, headers={
        "Range": f"bytes=0-{max_compressed-1}",
        "User-Agent": "eBOSS-DR16-mock-header-inventory/1.0",
        "Accept-Encoding": "identity",
    })
    with urlopen(request, timeout=timeout) as response:
        final = urlsplit(response.geturl())
        if (final.scheme != "https" or final.netloc != base.netloc
                or not final.path.startswith(base.path)):
            raise ValueError("Mock URL redirected outside the approved release")
        if response.status != 206:
            raise ValueError(f"Expected a bounded HTTP 206 response; got {response.status}")
        content_range = response.headers.get("Content-Range", "")
        match = re.fullmatch(r"bytes 0-(\d+)/(\d+)", content_range)
        if not match:
            raise ValueError("Missing or malformed bounded Content-Range")
        expected_length = int(match.group(1)) + 1
        if expected_length > max_compressed:
            raise ValueError("HTTP Range response exceeds requested byte limit")
        payload = response.read(max_compressed + 1)
        if len(payload) != expected_length:
            raise ValueError("HTTP Range response length differs from header")
    if not payload.startswith(b"\x1f\x8b"):
        raise ValueError("Released mock does not start with a gzip header")
    decoder = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    output = decoder.decompress(payload, max_decompressed + 1)
    if len(output) > max_decompressed:
        raise ValueError("Uncompressed FITS prefix exceeds safe byte limit")
    result = metadata(output)
    if result["header_rows"] < 1 or result["row_bytes"] < 1:
        raise ValueError("Mock BINTABLE has no rows or no row storage")
    return {
        "http_status": 206,
        "remote_compressed_bytes_reported": int(match.group(2)),
        "compressed_prefix_bytes_read": len(payload),
        "fits_prefix_bytes_decoded": len(output),
        "fits_header": result,
    }


def compare_schemas(records: list[dict]) -> dict:
    groups = defaultdict(set)
    for record in records:
        if record.get("status") != "header_parsed":
            continue
        header = record["fits_header"]
        groups[(record["tracer"], record["role"])].add(
            (tuple(header["column_names"]), tuple(header["column_formats"])))
    return {
        tracer + "_" + role: {
            "schema_count_across_checked_caps_and_ids": len(groups[(tracer, role)]),
            "consistent_within_sample": len(groups[(tracer, role)]) == 1,
        }
        for tracer in TRACERS for role in ROLES
    }


def self_test() -> None:
    from io import BytesIO
    from astropy.io import fits
    import numpy as np
    table = fits.BinTableHDU.from_columns([
        fits.Column(name="RA", format="D", array=np.array([1.0, 2.0])),
        fits.Column(name="DEC", format="D", array=np.array([0.0, 1.0])),
        fits.Column(name="Z", format="D", array=np.array([0.7, 0.8])),
    ])
    b = BytesIO()
    fits.HDUList([fits.PrimaryHDU(), table]).writeto(b)
    pack = zlib.compressobj(wbits=16 + zlib.MAX_WBITS)
    compressed = pack.compress(b.getvalue()) + pack.flush()
    dec = zlib.decompressobj(wbits=16 + zlib.MAX_WBITS)
    h = metadata(dec.decompress(compressed[:len(compressed) // 2], MAX_DECOMPRESSED))
    assert h["header_rows"] == 2 and h["has_ra_dec_z"]
    assert mock_path("eBOSS_LRG", "NGC", "dat", 1).endswith(
        "EZmock_realistic_eBOSS_LRG_NGC_v7_0001.dat.fits.gz")
    print("Bounded mock gzip/FITS metadata self-tests passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/mock_header_audit.json")
    parser.add_argument("--ids", default="1,500,1000")
    parser.add_argument("--timeout", type=float, default=35)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    ids = sorted(set(int(x.strip()) for x in args.ids.split(",")))
    if not ids or any(not 1 <= i <= 1000 for i in ids) or args.timeout <= 0:
        parser.error("Choose mock IDs between 1 and 1000 and a positive timeout")
    records, errors = [], []
    for rid in ids:
        for cap in CAPS:
            for tracer in TRACERS:
                for role in ROLES:
                    relative = mock_path(tracer, cap, role, rid)
                    record = {
                        "realization_id": rid, "cap": cap,
                        "tracer": tracer, "role": role,
                        "path": "v1_0_0/realistic/" + relative,
                    }
                    try:
                        record.update(inspect_prefix(BASE + relative, args.timeout,
                                                     MAX_COMPRESSED,
                                                     MAX_DECOMPRESSED))
                        record["status"] = "header_parsed"
                        print("MOCK_HEADER_OK", relative, "rows",
                              record["fits_header"]["header_rows"], flush=True)
                    except (OSError, ValueError, KeyError, TimeoutError,
                            zlib.error) as exc:
                        record["status"] = "header_unavailable"
                        record["error"] = str(exc)
                        errors.append(f"{relative}: {exc}")
                        print("MOCK_HEADER_ERROR", errors[-1], flush=True)
                    records.append(record)
    schemas = compare_schemas(records)
    parsed = [r for r in records if r["status"] == "header_parsed"]
    structure_ok = (len(parsed) == len(ids) * len(CAPS) * len(TRACERS) * len(ROLES)
                    and all(r["fits_header"]["has_ra_dec_z"] for r in parsed)
                    and all(s["consistent_within_sample"] for s in schemas.values()))
    result = {
        "study": "eBOSS DR16 matched realistic EZmock bounded FITS header audit",
        "release_url": BASE, "realization_ids_checked": ids,
        "status": "sample_headers_validated" if structure_ok else "partial",
        "records": records, "schema_consistency": schemas, "errors": errors,
        "full_mock_files_downloaded": False,
        "mock_rows_or_mock_weights_validated": False,
        "physical_pairing_or_covariance_validated": False,
        "odd_sector_observables_computed": False,
        "note": ("A matching filename and a parsed BINTABLE header do not "
                 "establish the joint mock selection or covariance. Validation "
                 "of full paired realizations is a subsequent gate."),
    }
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("MOCK_HEADER_AUDIT", path, result["status"], flush=True)
    return 0 if structure_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
