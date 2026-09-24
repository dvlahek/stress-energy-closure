#!/usr/bin/env python3
"""Inspect eBOSS catalogue metadata without measuring clustering.

We list release-directory filenames or inspect FITS headers and column names.
No galaxy positions, redshifts, pair counts or odd multipoles are evaluated.
"""
from __future__ import annotations

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urljoin
from urllib.request import Request, urlopen


DEFAULT_INDEX = "https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/"


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)
                break


def read_index(url: str) -> dict:
    req = Request(url, headers={"User-Agent": "scientific-catalogue-metadata-audit/1.0"})
    with urlopen(req, timeout=45) as response:
        effective_url = response.url
        content_type = response.headers.get("Content-Type", "")
        html = response.read(4_000_000).decode("utf-8", errors="replace")
    parser = Links()
    parser.feed(html)
    found = []
    for href in parser.links:
        name = href.rsplit("/", 1)[-1]
        if re.search(r"LRG|ELG|DR16|clustering|\.fits|EZmock", href, re.I):
            found.append({"name": name or href, "url": urljoin(effective_url, href)})
    found.sort(key=lambda item: item["name"].lower())
    return {
        "requested_url": url,
        "effective_url": effective_url,
        "content_type": content_type,
        "matched_links": found,
        "note": (
            "The directory listing is a catalogue-discovery aid. File versions, "
            "tracer samples and columns must be checked before pair counting."
        ),
    }


def inspect_fits(path: str) -> dict:
    from astropy.io import fits

    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(p)
    tables = []
    with fits.open(p, mode="readonly", memmap=True, lazy_load_hdus=True) as hdul:
        for i, hdu in enumerate(hdul):
            if not isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                continue
            names = list(hdu.columns.names)
            tables.append({
                "hdu": i,
                "name": hdu.name,
                "n_rows_header": int(hdu.header.get("NAXIS2", 0)),
                "columns": names,
                "position_columns_present": {
                    name: (name in names) for name in ("RA", "DEC", "Z")
                },
                "weight_columns": [
                    name for name in names
                    if re.search(r"WEIGHT|FKP|SYS|COMP|NOZ|CP", name, re.I)
                ],
            })
    if not tables:
        raise ValueError(f"No FITS table HDU found in {p}")
    return {
        "path": str(p),
        "filename": p.name,
        "size_bytes": p.stat().st_size,
        "table_hdus": tables,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index-url", help="Inspect an official directory listing without downloading FITS files.")
    ap.add_argument("--fits", nargs="+", help="Inspect one or more local FITS files (headers only).")
    ap.add_argument("--out", help="Optional JSON output path.")
    args = ap.parse_args()
    if not args.index_url and not args.fits:
        ap.error("Provide --index-url and/or --fits.")
    summary = {
        "scope": "eBOSS DR16 catalogue metadata",
        "uses_pair_counts": False,
        "uses_odd_data_vector": False,
    }
    if args.index_url:
        summary["release_index"] = read_index(args.index_url)
    if args.fits:
        summary["fits_files"] = [inspect_fits(p) for p in args.fits]
    payload = json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        dest = Path(args.out)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
