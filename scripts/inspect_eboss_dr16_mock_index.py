#!/usr/bin/env python3
"""Inventory public eBOSS DR16 EZmock directory entries without reading mock data.

The output is a discovery index only. Identical realization IDs, common
initial conditions and usable cross-tracer mock covariance require separate
validation of the released products.
"""
from __future__ import annotations

import argparse
from collections import deque
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urljoin, urlsplit
from urllib.request import Request, urlopen


ROOT_URL = "https://data.sdss.org/sas/dr17/eboss/lss/EZmocks/"
FILE_SUFFIXES = (
    ".fits", ".fits.gz", ".fits.fz", ".fit", ".fit.gz",
    ".csv", ".dat", ".dat.gz", ".txt", ".txt.gz", ".json", ".h5",
    ".hdf5", ".tar", ".tar.gz", ".zip", ".gz", ".parquet",
)


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if href:
                self.hrefs.append(href)


def fetch_listing(url: str, timeout: float) -> tuple[str, list[str]]:
    request = Request(url, headers={"User-Agent": "eBOSS-DR16-EZmock-inventory/1.0"})
    with urlopen(request, timeout=timeout) as response:
        final = response.geturl()
        ctype = response.headers.get("Content-Type", "")
        if "html" not in ctype.lower():
            raise ValueError(f"Expected an HTML directory listing, got {ctype!r}")
        body = response.read(2 * 1024 * 1024 + 1)
        if len(body) > 2 * 1024 * 1024:
            raise ValueError("Listing exceeds the 2 MiB read limit")
        charset = response.headers.get_content_charset() or "utf-8"
    parser = Links()
    parser.feed(body.decode(charset, errors="replace"))
    return final, parser.hrefs


def within_release(target: str, root: str) -> bool:
    a, b = urlsplit(target), urlsplit(root)
    return (a.scheme == b.scheme == "https" and a.netloc == b.netloc
            and a.path.startswith(b.path) and not a.query and not a.fragment)


def sample_type(path: str) -> str:
    lower = path.lower()
    if "elg" in lower:
        return "ELG"
    if "lrg" in lower or "cmass" in lower:
        return "LRG_or_CMASS"
    if "qso" in lower:
        return "QSO"
    return "unspecified"


def run(max_depth: int, max_pages: int, timeout: float) -> dict:
    base = ROOT_URL
    queue = deque([(base, 0)])
    seen: set[str] = set()
    records: dict[str, dict] = {}
    errors: list[str] = []
    page_count = 0

    while queue and page_count < max_pages:
        url, depth = queue.popleft()
        if url in seen:
            continue
        seen.add(url)
        try:
            final, hrefs = fetch_listing(url, timeout)
        except (HTTPError, URLError, ValueError, TimeoutError, OSError) as exc:
            errors.append(f"{url}: {exc}")
            continue
        page_count += 1
        if not within_release(final, base):
            errors.append(f"Listing redirected outside the allowed release: {url}")
            continue
        for href in hrefs:
            target = urljoin(final, href)
            if not within_release(target, base) or target == final:
                continue
            rel = unquote(urlsplit(target).path[len(urlsplit(base).path):])
            if not rel:
                continue
            is_dir = urlsplit(target).path.endswith("/")
            record = {
                "path": rel,
                "url": target,
                "kind": "directory" if is_dir else "file",
                "tracer_hint": sample_type(rel),
                "discovered_at_depth": depth,
            }
            records[target] = record
            if is_dir and depth < max_depth and target not in seen:
                queue.append((target, depth + 1))

    truncated = bool(queue)
    data_files = [
        r for r in records.values()
        if r["kind"] == "file" and r["path"].lower().endswith(FILE_SUFFIXES)
    ]
    ordered = sorted(records.values(), key=lambda r: r["path"])
    return {
        "study": "eBOSS DR16 EZmock public directory inventory",
        "release_url": base,
        "status": "partial" if errors or truncated else ("listed" if page_count else "unavailable"),
        "pages_inspected": page_count,
        "max_depth": max_depth,
        "max_pages": max_pages,
        "truncated": truncated,
        "entries": ordered,
        "potential_data_files": sorted(data_files, key=lambda r: r["path"]),
        "errors": errors,
        "downloaded_mock_catalogues": False,
        "joint_lrg_elg_realizations_verified": False,
        "note": (
            "This directory listing does not establish matching realization IDs, "
            "catalogue contents, tracer selection, per-realization randoms or "
            "eligibility for cross-tracer covariance estimation."
        ),
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="eboss_workspace/mock_release_index.json")
    p.add_argument("--max-depth", type=int, default=2)
    p.add_argument("--max-pages", type=int, default=30)
    p.add_argument("--timeout", type=float, default=30)
    args = p.parse_args()
    if args.max_depth < 0 or args.max_pages < 1 or args.timeout <= 0:
        p.error("Require max-depth >= 0, max-pages >= 1 and timeout > 0")
    result = run(args.max_depth, args.max_pages, args.timeout)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Mock inventory: {path}")
    print(f"Status: {result['status']}; pages: {result['pages_inspected']}; "
          f"entries: {len(result['entries'])}; candidate files: "
          f"{len(result['potential_data_files'])}")
    for entry in result["entries"][:40]:
        print(entry["kind"], entry["path"])
    for error in result["errors"]:
        print("Directory error:", error)
    return 0 if result["pages_inspected"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
