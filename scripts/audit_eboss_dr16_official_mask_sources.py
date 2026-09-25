#!/usr/bin/env python3
"""Discover eBOSS DR16 official mask files without reading any galaxy data.

This is a SOURCE INVENTORY, not an angular-mask implementation. Retrieve
bounded official directory listings; record exact per-ELG-chunk maskbit
filenames and the three published ELG extra-veto polygon names. Report
LRG-related names from the DR16 index without guessing which MANGLE
products or veto bits belong to our tracer-pair selection.

Optional --download-extra pins SHA256 of the three *listed* ELG extra
polygons, but does not apply masks or certify a common footprint.
"""
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urljoin, urlsplit, unquote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_official_mask_source_protocol_2026-09-25.json"
DEFAULT_OUT = Path("eboss_workspace/official_mask_inventory")
ALLOWED_ORIGINS = ("https://portal.nersc.gov",
                   "https://data.sdss.org")
MAX_INDEX_BYTES = 16 * 1024 * 1024
MAX_POLYGON_BYTES = 1024 * 1024 * 1024
CHUNKS = ("eboss21", "eboss22", "eboss23", "eboss25")
EXTRAS = ("ELG_centerpost.ply",
          "ELG_TDSSFES_62arcsec.pix.snap.balk.ply",
          "ebosselg_badphot.26Aug2019.ply")
ALLOWED_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")
MASK_RE = re.compile(r"^mask-(eboss21|eboss22|eboss23|eboss25)-.+\.fits\.gz$")


class IndexLinks(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            href = dict(attrs).get("href")
            if isinstance(href, str):
                self.hrefs.append(href)


def accept_root(root_url):
    split = urlsplit(root_url)
    origin = split.scheme + "://" + split.netloc
    if (origin not in ALLOWED_ORIGINS or split.query or split.fragment
            or not split.path.endswith("/")
            or not split.path.startswith("/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/")
            and not split.path.startswith("/sas/dr17/eboss/lss/catalogs/DR16/")):
        raise ValueError("Unapproved DR17-published DR16 mask directory root")
    return root_url


def listing_names(listing_html, base_url):
    root = urlsplit(accept_root(base_url))
    parser = IndexLinks()
    parser.feed(listing_html)
    result = {}
    for href in parser.hrefs:
        split = urlsplit(href)
        if split.query or split.fragment:
            continue
        name = unquote(split.path).rstrip("/").split("/")[-1]
        if not ALLOWED_NAME.fullmatch(name):
            continue
        url = urljoin(base_url, href)
        full = urlsplit(url)
        if (full.scheme != root.scheme or full.netloc != root.netloc
                or not full.path.startswith(root.path)
                or full.query or full.fragment):
            continue
        is_dir = split.path.endswith("/")
        if name in result and result[name] != (url, is_dir):
            raise ValueError("Conflicting official index item: " + name)
        result[name] = (url, is_dir)
    return result


def fetch_limited(url, *, max_bytes, root_url, timeout=90):
    base = urlsplit(accept_root(root_url))
    target = urlsplit(url)
    if (target.scheme != base.scheme or target.netloc != base.netloc
            or not target.path.startswith(base.path)
            or target.query or target.fragment):
        raise ValueError("Official source URL outside approved DR16 root")
    req = Request(url, headers={
        "User-Agent": "eBOSS-DR16-blinded-official-mask-input/1.0",
        "Accept-Encoding": "identity",
    })
    with urlopen(req, timeout=timeout) as response:
        resolved = urlsplit(response.geturl())
        if (resolved.scheme != base.scheme or resolved.netloc != base.netloc
                or not resolved.path.startswith(base.path)):
            raise ValueError("Source redirected outside approved DR16 release")
        nbytes, digest, chunks = 0, hashlib.sha256(), []
        while True:
            block = response.read(min(1024 * 1024, max_bytes - nbytes + 1))
            if not block:
                break
            nbytes += len(block)
            if nbytes > max_bytes:
                raise ValueError("Official directory/product exceeds declared byte bound")
            digest.update(block)
            chunks.append(block)
    if nbytes == 0:
        raise ValueError("Empty official directory/product")
    return b"".join(chunks), digest.hexdigest(), nbytes


def scan(protocol, out_dir, *, timeout, download_extra):
    source = protocol["source_data_roots"]
    root = accept_root(source["nersc_mirror"])
    base = root + source["elg_mask_directory"]
    root_bytes, root_sha, root_size = fetch_limited(
        root, max_bytes=MAX_INDEX_BYTES, root_url=root, timeout=timeout)
    elg_bytes, elg_sha, elg_size = fetch_limited(
        base, max_bytes=MAX_INDEX_BYTES, root_url=root, timeout=timeout)
    root_entries = listing_names(root_bytes.decode("utf-8"), root)
    elg_entries = listing_names(elg_bytes.decode("utf-8"), base)
    mask_groups = {chunk: [] for chunk in CHUNKS}
    for name, (url, is_dir) in elg_entries.items():
        match = MASK_RE.fullmatch(name)
        if match and not is_dir:
            mask_groups[match.group(1)].append({
                "filename": name, "source_url": url,
            })
    for chunk in CHUNKS:
        mask_groups[chunk].sort(key=lambda x: x["filename"])
    extras = {}
    for name in EXTRAS:
        entry = elg_entries.get(name)
        extras[name] = {
            "listed": bool(entry and not entry[1]),
            "source_url": entry[0] if entry and not entry[1] else None,
        }
    lrg_names = [
        {"filename": n, "directory": is_dir, "source_url": url}
        for n, (url, is_dir) in sorted(root_entries.items())
        if ("lrg" in n.lower() or "mask" in n.lower()
            or "veto" in n.lower() or "sector" in n.lower())
    ]
    if not root_entries or not elg_entries:
        raise ValueError("Official listing has no usable safe links")
    report = {
        "status": "official_mask_index_incomplete",
        "study": "eBOSS DR16 official mask-source inventory only",
        "approved_source_root": root,
        "approved_elg_directory": base,
        "root_index": {"sha256": root_sha, "bytes": root_size,
                       "usable_entries": len(root_entries)},
        "elg_index": {"sha256": elg_sha, "bytes": elg_size,
                      "usable_entries": len(elg_entries)},
        "elg_chunks": {
            k: {"listed_maskbits": len(v),
                "maskbit_filenames": [e["filename"] for e in v]}
            for k, v in mask_groups.items()
        },
        "elg_extra_polygons": extras,
        "lrg_root_candidates": lrg_names,
        "upstream_brickmask": protocol["published_brickmask"],
        "source_filenames_confirmed": False,
        "physical_joint_mask_certified": False,
        "observed_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "errors": [],
        "warning": ("Listing and optional file-byte SHA256 do not specify "
                    "maskbit rejection logic, completeness, continuous joint "
                    "geometry or pair-level window; official LRG veto products "
                    "must be identified separately."),
    }
    if any(not mask_groups[k] for k in CHUNKS):
        report["errors"].append("Missing one or more predeclared ELG maskbits chunk families")
    for name in EXTRAS:
        if not extras[name]["listed"]:
            report["errors"].append("Not listed in ELGmasks: " + name)
    if download_extra and not report["errors"]:
        for name in EXTRAS:
            destination = out_dir / "official_elg_extra" / name
            if destination.exists():
                digest = hashlib.sha256(destination.read_bytes()).hexdigest()
                size = destination.stat().st_size
            else:
                payload, digest, size = fetch_limited(
                    extras[name]["source_url"], max_bytes=MAX_POLYGON_BYTES,
                    root_url=root, timeout=timeout)
                destination.parent.mkdir(parents=True, exist_ok=True)
                temporary = destination.with_suffix(".part")
                temporary.write_bytes(payload)
                temporary.replace(destination)
            if size < 16:
                raise ValueError("Truncated official ELG veto polygon: " + name)
            extras[name]["local_filename"] = str(destination)
            extras[name]["sha256"] = digest
            extras[name]["bytes"] = size
    if not report["errors"]:
        report["status"] = ("official_elg_extra_bytes_pinned_and_chunk_indexed"
                            if download_extra
                            else "official_elg_chunk_index_discovered")
        report["source_filenames_confirmed"] = True
    return report


def self_test():
    base = "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/ELGmasks/"
    html = ('<html><a href="../">parent</a>'
            '<a href="mask-eboss21-0112p345.fits.gz">21</a>'
            '<a href="mask-eboss22-0991m121.fits.gz">22</a>'
            '<a href="ELG_centerpost.ply">centerpost</a>'
            '<a href="https://example.org/evil.fits">evil</a>'
            '<a href="../../outside.ply">escape</a>'
            '<a href="mask-eboss25-bad.fits.gz?x=1">query</a></html>')
    links = listing_names(html, base)
    assert set(links) == {
        "mask-eboss21-0112p345.fits.gz",
        "mask-eboss22-0991m121.fits.gz", "ELG_centerpost.ply"}
    assert MASK_RE.fullmatch("mask-eboss23-1000p020.fits.gz")
    assert not MASK_RE.fullmatch("mask-eboss24-1000p020.fits.gz")
    try:
        accept_root("https://example.com/masks/")
    except ValueError:
        pass
    else:
        raise AssertionError("Non-official origin accepted")
    src = json.loads(PROTOCOL.read_text())
    assert tuple(src["elg_chunks"]) == CHUNKS
    assert tuple(src["elg_extra_polygons"]) == EXTRAS
    assert src["observed_odd_data_vector_read"] is False
    print("EBOSS_OFFICIAL_MASK_SOURCE_INVENTORY_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--download-extra", action="store_true",
                    help="Also SHA-pin the three listed official ELG extra veto polygons")
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--timeout", type=float, default=90)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout <= 0:
        ap.error("timeout must be positive")
    protocol = json.loads(PROTOCOL.read_text())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    try:
        report = scan(protocol, args.out_dir, timeout=args.timeout,
                      download_extra=args.download_extra)
    except Exception as exc:
        report = {
            "status": "official_mask_source_inventory_incomplete",
            "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "physical_joint_mask_certified": False,
        }
    output = args.out_dir / "official_mask_index.json"
    tmp = output.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    tmp.replace(output)
    print("EBOSS_MASK_SOURCE_INVENTORY", report["status"],
          output, flush=True)
    if report["errors"]:
        print("EBOSS_MASK_SOURCE_ERRORS", *report["errors"], sep="\n", flush=True)
        return 2
    for chunk, item in report["elg_chunks"].items():
        print("EBOSS_ELG_MASKBITS_LISTED", chunk, item["listed_maskbits"],
              flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
