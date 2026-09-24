#!/usr/bin/env python3
"""Read only public eBOSS DR16 FITS headers; do not inspect catalogue rows.

The output records released FITS columns and header row counts, not valid
galaxy counts, weights, redshift distributions, overlap or an odd data vector.
Only a bounded initial byte range of each publicly listed file is requested.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from inspect_eboss_dr16_catalogues import DEFAULT_URL, list_remote

BLOCK = 2880


def parse_header(data: bytes, start: int) -> tuple[dict, int]:
    if start % BLOCK or start >= len(data):
        raise ValueError("FITS header starts outside available aligned data")
    header: dict[str, object] = {}
    pos = start
    while pos + 80 <= len(data):
        card = data[pos:pos + 80]
        try:
            text = card.decode("ascii")
        except UnicodeDecodeError as exc:
            raise ValueError("Non-ASCII FITS header card") from exc
        key = text[:8].strip()
        pos += 80
        if key == "END":
            return header, ((pos + BLOCK - 1) // BLOCK) * BLOCK
        if not key or text[8:10] != "= ":
            continue
        raw = text[10:80].strip()
        if raw.startswith("'"):
            value = raw[1:].split("'", 1)[0].replace("''", "'").strip()
        else:
            value = raw.split("/", 1)[0].strip()
            if value in ("T", "F"):
                value = value == "T"
            elif re.fullmatch(r"[+-]?\d+", value):
                value = int(value)
        header[key] = value
    raise ValueError("END card not found in bounded FITS header")


def hdu_data_length(header: dict) -> int:
    naxis = int(header.get("NAXIS", 0))
    if naxis == 0:
        return 0
    if str(header.get("XTENSION", "")).upper() == "BINTABLE":
        return (int(header["NAXIS1"]) * int(header["NAXIS2"])
                + int(header.get("PCOUNT", 0))) * int(header.get("GCOUNT", 1))
    size = abs(int(header["BITPIX"])) // 8
    for axis in range(1, naxis + 1):
        size *= int(header[f"NAXIS{axis}"])
    return (size + int(header.get("PCOUNT", 0))) * int(header.get("GCOUNT", 1))


def metadata(data: bytes) -> dict:
    primary, pos = parse_header(data, 0)
    if primary.get("SIMPLE") is not True:
        raise ValueError("The primary FITS header does not declare SIMPLE = T")
    pos += ((hdu_data_length(primary) + BLOCK - 1) // BLOCK) * BLOCK
    for _ in range(6):
        if pos >= len(data):
            break
        header, next_pos = parse_header(data, pos)
        if str(header.get("XTENSION", "")).upper() == "BINTABLE":
            nfields = int(header.get("TFIELDS", 0))
            columns = [str(header.get(f"TTYPE{i}", "")).strip()
                       for i in range(1, nfields + 1)]
            forms = [str(header.get(f"TFORM{i}", "")).strip()
                     for i in range(1, nfields + 1)]
            names = {c.upper() for c in columns}
            return {
                "extension": "BINTABLE",
                "header_rows": int(header["NAXIS2"]),
                "row_bytes": int(header["NAXIS1"]),
                "column_count": nfields,
                "column_names": columns,
                "column_formats": forms,
                "has_ra_dec_z": all(c in names for c in ("RA", "DEC", "Z")),
                "weight_columns": [c for c in columns if "WEIGHT" in c.upper()],
            }
        pos = next_pos + ((hdu_data_length(header) + BLOCK - 1) // BLOCK) * BLOCK
    raise ValueError("No BINTABLE extension is present in the bounded header bytes")


def inspect(url: str, timeout: float, max_bytes: int) -> dict:
    allowed = urlsplit(DEFAULT_URL)
    supplied = urlsplit(url)
    if (supplied.scheme != "https" or supplied.netloc != allowed.netloc
            or not supplied.path.startswith(allowed.path) or supplied.query):
        raise ValueError("Refusing a URL outside the published SDSS DR16 catalogue directory")
    request = Request(url, headers={
        "Range": f"bytes=0-{max_bytes - 1}",
        "User-Agent": "eBOSS-DR16-FITS-header-inventory/1.0",
    })
    with urlopen(request, timeout=timeout) as response:
        final = urlsplit(response.geturl())
        if final.netloc != allowed.netloc or not final.path.startswith(allowed.path):
            raise ValueError("FITS request redirected outside the SDSS release")
        size_header = response.headers.get("Content-Range", "")
        total = re.search(r"/(\d+)$", size_header)
        remote_bytes = int(total.group(1)) if total else (
            int(response.headers["Content-Length"])
            if response.headers.get("Content-Length", "").isdigit()
            and response.status == 200 else None
        )
        record = {
            "http_status": response.status,
            "content_range": size_header or None,
            "remote_bytes_reported": remote_bytes,
            "etag": response.headers.get("ETag"),
            "last_modified": response.headers.get("Last-Modified"),
        }
        payload = response.read(max_bytes)
    record["bounded_bytes_read"] = len(payload)
    record["fits_header"] = metadata(payload)
    return record


def self_test() -> None:
    def card(key: str, value: str = "") -> str:
        return (f"{key:<8}= {value}" if value else key).ljust(80)

    def block(cards: list[str]) -> bytes:
        s = "".join(cards + [card("END")]).encode("ascii")
        return s.ljust((len(s) + BLOCK - 1) // BLOCK * BLOCK, b" ")

    primary = block([card("SIMPLE", "T"), card("BITPIX", "8"), card("NAXIS", "0")])
    table = block([card("XTENSION", "'BINTABLE'"), card("BITPIX", "8"),
                   card("NAXIS", "2"), card("NAXIS1", "32"), card("NAXIS2", "11"),
                   card("PCOUNT", "0"), card("GCOUNT", "1"), card("TFIELDS", "4"),
                   card("TTYPE1", "'RA'"), card("TFORM1", "'1D'"),
                   card("TTYPE2", "'DEC'"), card("TFORM2", "'1D'"),
                   card("TTYPE3", "'Z'"), card("TFORM3", "'1D'"),
                   card("TTYPE4", "'WEIGHT_FKP'"), card("TFORM4", "'1D'")])
    result = metadata(primary + table)
    assert result["header_rows"] == 11
    assert result["has_ra_dec_z"] and result["weight_columns"] == ["WEIGHT_FKP"]
    print("Synthetic FITS header parser test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/fits_header_inventory.json")
    parser.add_argument("--max-bytes", type=int, default=128 * 1024)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.max_bytes < 2 * BLOCK or args.timeout <= 0:
        parser.error("max-bytes must be at least 5760 and timeout positive")
    records, listing_url = list_remote(DEFAULT_URL, args.timeout)
    results = []
    for record in records:
        result = {k: record[k] for k in ("tracer", "cap", "role", "name", "url")}
        try:
            result.update(inspect(record["url"], args.timeout, args.max_bytes))
            result["header_status"] = "parsed"
        except (OSError, ValueError, KeyError, TimeoutError) as exc:
            result["header_status"] = "unavailable"
            result["error"] = str(exc)
        results.append(result)
    complete = len(results) == 8 and all(
        r.get("header_status") == "parsed" and
        r["fits_header"]["has_ra_dec_z"] for r in results
    )
    report = {
        "study": "eBOSS DR16 public clustering FITS header metadata",
        "listing_url_checked": listing_url,
        "status": "headers_parsed" if complete else "partial",
        "records": results,
        "no_catalogue_rows_interpreted": True,
        "no_odd_sector_vector_measured": True,
        "full_file_sha256_checked": False,
        "note": (
            "Header NAXIS2 is a declared row count, not a validated count of "
            "usable tracers. FITS data rows, angular/redshift overlap, weights, "
            "and mock covariance remain unverified."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"FITS header inventory: {output}; status: {report['status']}")
    for r in results:
        if r["header_status"] == "parsed":
            h = r["fits_header"]
            print(r["name"], "rows_in_header=", h["header_rows"],
                  "columns=", h["column_count"], "weights=", h["weight_columns"])
        else:
            print(r["name"], "header unavailable:", r["error"])
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
