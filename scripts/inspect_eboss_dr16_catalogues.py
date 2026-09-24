#!/usr/bin/env python3
"""Inspect the public eBOSS DR16 LSS catalogue index without downloading FITS data.

We enumerate published LRG and ELG clustering files and distinguish the
catalogue hemispheres, data/random roles, and reconstructed BAO products.
The resulting JSON records catalogue locations; it does not validate FITS
columns, selection functions, or the physical overlap of the two tracers.
"""
from __future__ import annotations

import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import sys
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import Request, urlopen


DEFAULT_URL = "https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/"
FITS_SUFFIXES = (".fits", ".fits.gz", ".fits.fz", ".fit", ".fit.gz")


class CatalogueLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.append(href)


def classify(name):
    """Classify a published filename without assuming a particular DR16 suffix."""
    lower = name.lower()
    if not lower.endswith(FITS_SUFFIXES):
        return None
    if "clustering" not in lower or re.search(r"(?:^|[_-])rec(?:[_.-]|$)", lower):
        return None

    tracer = None
    if re.search(r"eboss[_-]lrg[_-]clustering", lower):
        tracer = "LRG"
    elif re.search(r"eboss[_-]elg[_-]clustering", lower):
        tracer = "ELG"
    if tracer is None:
        return None

    cap_match = re.search(r"(?:^|[_-])(NGC|SGC)(?:[_\.-]|$)", name, re.I)
    cap = cap_match.group(1).upper() if cap_match else "unspecified"

    if re.search(r"(?:^|[_\.-])(ran|random|randoms)(?:[_\.-]|$)", name, re.I):
        role = "random"
    elif re.search(r"(?:^|[_\.-])(dat|data)(?:[_\.-]|$)", name, re.I):
        role = "data"
    else:
        role = "unclassified"
    return tracer, cap, role


def list_remote(url, timeout):
    """Read an SDSS directory listing; do not download the FITS targets."""
    req = Request(url, headers={"User-Agent": "eBOSS-DR16-catalogue-index/1.0"})
    with urlopen(req, timeout=timeout) as response:
        final_url = response.geturl()
        content_type = response.headers.get("Content-Type", "")
        if "html" not in content_type.lower():
            raise ValueError(
                f"Expected an HTML directory listing; received {content_type!r}"
            )
        raw = response.read(8 * 1024 * 1024 + 1)
        if len(raw) > 8 * 1024 * 1024:
            raise ValueError("Directory listing exceeds the 8 MiB inspection limit.")
        charset = response.headers.get_content_charset() or "utf-8"
        html = raw.decode(charset, errors="replace")

    parser = CatalogueLinks()
    parser.feed(html)
    base = urlsplit(final_url)
    base_path = base.path.rstrip("/") + "/"
    records = []
    for href in parser.hrefs:
        target = urljoin(final_url, href)
        parts = urlsplit(target)
        if parts.netloc != base.netloc or not parts.path.startswith(base_path):
            continue
        filename = unquote(parts.path.rsplit("/", 1)[-1])
        kind = classify(filename)
        if kind is None:
            continue
        tracer, cap, role = kind
        records.append({
            "tracer": tracer,
            "cap": cap,
            "role": role,
            "name": filename,
            "url": target,
            "source": "SDSS DR16 SAS directory listing",
        })
    return records, final_url


def list_local(directory):
    records = []
    if not directory.is_dir():
        raise FileNotFoundError(f"Local catalogue directory does not exist: {directory}")
    for path in sorted(directory.rglob("*")):
        if not path.is_file():
            continue
        kind = classify(path.name)
        if kind is None:
            continue
        tracer, cap, role = kind
        records.append({
            "tracer": tracer,
            "cap": cap,
            "role": role,
            "name": path.name,
            "path": str(path.resolve()),
            "bytes": path.stat().st_size,
            "source": "local file",
        })
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", default="eboss_workspace/release_index.json",
        help="Output JSON path (default: eboss_workspace/release_index.json)",
    )
    parser.add_argument("--base-url", default=DEFAULT_URL)
    parser.add_argument(
        "--local-dir", type=Path,
        help="Also inspect previously downloaded catalogue filenames in this directory.",
    )
    parser.add_argument(
        "--skip-remote", action="store_true",
        help="Inspect only --local-dir; do not contact the SDSS server.",
    )
    parser.add_argument("--timeout", type=float, default=25.0)
    args = parser.parse_args()

    if args.skip_remote and args.local_dir is None:
        parser.error("--skip-remote requires --local-dir")

    records = []
    errors = []
    checked_url = None
    if not args.skip_remote:
        try:
            records, checked_url = list_remote(args.base_url, args.timeout)
        except (HTTPError, URLError, ValueError, TimeoutError, OSError) as exc:
            errors.append(f"SDSS directory listing unavailable: {exc}")

    if args.local_dir is not None:
        try:
            records.extend(list_local(args.local_dir))
        except OSError as exc:
            errors.append(f"Local catalogue inspection failed: {exc}")

    # Keep remote and local occurrences separately, but discard duplicate entries.
    unique = {(r["source"], r.get("url", r.get("path", ""))): r for r in records}
    records = sorted(
        unique.values(), key=lambda r: (r["tracer"], r["cap"], r["role"], r["name"], r["source"])
    )
    by_tracer = {
        tracer: {
            cap: {role: [] for role in ("data", "random", "unclassified")}
            for cap in ("NGC", "SGC", "unspecified")
        }
        for tracer in ("LRG", "ELG")
    }
    for record in records:
        by_tracer[record["tracer"]][record["cap"]][record["role"]].append(record)

    completeness = {
        tracer: {
            cap: {
                "has_data": bool(by_tracer[tracer][cap]["data"]),
                "has_random": bool(by_tracer[tracer][cap]["random"]),
            }
            for cap in ("NGC", "SGC")
        }
        for tracer in ("LRG", "ELG")
    }
    complete = all(
        completeness[tracer][cap]["has_data"]
        and completeness[tracer][cap]["has_random"]
        for tracer in ("LRG", "ELG")
        for cap in ("NGC", "SGC")
    )
    status = "complete" if complete else ("incomplete" if records else "unavailable")
    result = {
        "study": "eBOSS DR16 LRG–ELG catalogue inventory",
        "release_url": args.base_url,
        "listing_url_checked": checked_url,
        "status": status,
        "downloaded_fits": False,
        "records": records,
        "by_tracer": by_tracer,
        "completeness": completeness,
        "errors": errors,
        "note": (
            "The inventory records filenames and URLs only. Catalogue weights, "
            "FITS columns, survey selection and common LRG–ELG footprint require "
            "separate validation before the odd-sector estimator is applied."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Catalogue index: {output}")
    print(f"Status: {status}; matching files: {len(records)}")
    for tracer in ("LRG", "ELG"):
        for cap in ("NGC", "SGC"):
            ndata = len(by_tracer[tracer][cap]["data"])
            nrandom = len(by_tracer[tracer][cap]["random"])
            print(f"{tracer} {cap}: data={ndata}, random={nrandom}")
    for error in errors:
        print(error, file=sys.stderr)
    if status != "complete":
        print(
            "The index is incomplete. Inspect its errors and records; no catalogue "
            "paths should be inferred from missing listings.",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
