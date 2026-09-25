#!/usr/bin/env python3
"""Byte-pin DR16 LRG/ELG official polygon inputs, without applying a mask.

Consumes an index acquired from the published NERSC mirror and the
previously committed official-mask protocol. Downloads exactly seven
published LRG/QSO veto polygons, the published QSO+LRG *noveto*
footprint polygon, and three ELG extra veto polygons. Full SHA256 is
computed while streaming to disk. This is provenance only: presence
does not determine the official LRG veto composition, ELG brick bits,
survey completeness, or the final LRG x ELG selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_official_mask_source_protocol_2026-09-25.json"
LIMIT_PER_POLYGON = 512 * 1024 * 1024
EXPECTED_COUNT = 11


def digest_file(path):
    h = hashlib.sha256()
    total = 0
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
            total += len(block)
    return h.hexdigest(), total


def atomic_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def validate_url(url, root, expected_filename):
    src, candidate = urlsplit(root), urlsplit(url)
    if (src.scheme != "https" or candidate.scheme != src.scheme
            or candidate.netloc != src.netloc or candidate.query
            or candidate.fragment or not candidate.path.startswith(src.path)
            or candidate.path.split("/")[-1] != expected_filename):
        raise ValueError("Official polygon URL does not match pinned source tree")
    if not re.fullmatch(r"[A-Za-z0-9_.+-]+\.ply", expected_filename):
        raise ValueError("Unexpected polygon source filename")


def stream_pinned(url, path, root, filename, timeout):
    validate_url(url, root, filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".partial")
    tmp.unlink(missing_ok=True)
    h, total = hashlib.sha256(), 0
    req = Request(url, headers={
        "User-Agent": "eBOSS-DR16-pinned-official-polygon-byte-audit/1.0",
        "Accept-Encoding": "identity"})
    try:
        with urlopen(req, timeout=timeout) as response, tmp.open("wb") as f:
            redirected = response.geturl()
            validate_url(redirected, root, filename)
            while True:
                data = response.read(1024 * 1024)
                if not data:
                    break
                total += len(data)
                if total > LIMIT_PER_POLYGON:
                    raise ValueError("Official polygon exceeds fixed safety bound")
                h.update(data)
                f.write(data)
        if total < 16:
            raise ValueError("Official polygon is missing or truncated")
        tmp.replace(path)
        return h.hexdigest(), total
    finally:
        tmp.unlink(missing_ok=True)


def expected_products(protocol, inventory):
    roots = protocol["source_data_roots"]
    root = roots["nersc_mirror"]
    lrg_root = root + roots["lrg_mask_directory"]
    elg_root = root + roots["elg_mask_directory"]
    if (inventory.get("status") != "official_elg_chunk_index_discovered"
            or inventory.get("errors")
            or inventory.get("approved_source_root") != root
            or inventory.get("lrg_mask_directory") != lrg_root
            or inventory.get("root_index", {}).get("sha256")
            != protocol["source_index_sha256"]["dr16_root"]
            or inventory.get("lrg_index", {}).get("sha256")
            != protocol["source_index_sha256"]["lrg_directory"]
            or inventory.get("elg_index", {}).get("sha256")
            != protocol["source_index_sha256"]["elg_directory"]
            or inventory.get("observed_odd_data_vector_read") is not False
            or inventory.get("physical_joint_mask_certified") is not False):
        raise ValueError("Official source inventory missing or not SHA-pinned")
    if any(inventory["elg_chunks"][c]["listed_maskbits"] < 1
           for c in protocol["elg_chunks"]):
        raise ValueError("One or more published ELG chunk mask families are empty")
    lrg = {x["filename"]: x["source_url"]
           for x in inventory["lrg_published_products"]}
    extra = inventory["elg_extra_polygons"]
    root_candidates = {x["filename"]: x["source_url"]
                       for x in inventory["lrg_root_candidates"]
                       if not x["directory"]}
    products = []
    for name in protocol["lrg_published_polygons"]:
        if name not in lrg:
            raise ValueError("Missing pinned official LRG/QSO polygon: " + name)
        products.append({"role": "lrg_published_polygon", "filename": name,
                         "source_url": lrg[name], "root": lrg_root})
    root_file = protocol["lrg_root_geometry_polygon"]
    if root_file not in root_candidates:
        raise ValueError("Missing official LRG/QSO noveto geometry polygon")
    products.append({"role": "lrg_noveto_footprint_polygon",
                     "filename": root_file, "source_url": root_candidates[root_file],
                     "root": root})
    for name in protocol["elg_extra_polygons"]:
        item = extra.get(name, {})
        if item.get("listed") is not True or not item.get("source_url"):
            raise ValueError("Missing official ELG extra veto polygon: " + name)
        products.append({"role": "elg_extra_polygon", "filename": name,
                         "source_url": item["source_url"], "root": elg_root})
    if len(products) != EXPECTED_COUNT or len(set(p["filename"] for p in products)) != EXPECTED_COUNT:
        raise ValueError("Fixed official polygon source cohort is incomplete or duplicated")
    # Pin smaller ELG extra masks before the potentially large allsky LRG
    # polygons; a timed-out full transfer then preserves useful, clearly
    # labelled partial provenance instead of losing the whole run.
    priority = {"elg_extra_polygon": 0, "lrg_noveto_footprint_polygon": 1,
                "lrg_published_polygon": 2}
    products.sort(key=lambda p: (
        priority[p["role"]],
        p["filename"] == "allsky_bright_star_mask_pix.ply",
        p["filename"]))
    for item in products:
        validate_url(item["source_url"], item["root"], item["filename"])
    return products


def audit(protocol, inventory, out_dir, timeout, only_elg_extra=False):
    products = expected_products(protocol, inventory)
    if only_elg_extra:
        products = [p for p in products if p["role"] == "elg_extra_polygon"]
        if len(products) != 3:
            raise ValueError("Exactly three official ELG polygon inputs required")
    output = out_dir / "official_polygon_sha_manifest.json"
    if output.exists():
        prior = json.loads(output.read_text(encoding="utf-8"))
        if (prior.get("source_inventory_sha256")
                != hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()
                or prior.get("observed_odd_data_vector_read") is not False):
            raise ValueError("Previous polygon checkpoint refers to a different source index")
        old_records = {p["filename"]: p for p in prior.get("products", [])}
    else:
        old_records = {}
    inv_sha = hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()
    report = {
        "status": "official_polygon_sha_inventory_incomplete",
        "source_inventory_sha256": inv_sha,
        "source_index_sha256": protocol["source_index_sha256"],
        "upstream_brickmask": protocol["published_brickmask"],
        "products": [],
        "errors": [],
        "observed_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "physical_joint_mask_certified": False,
        "official_lrg_veto_semantics_certified": False,
        "elg_brickmask_products_fully_pinned": False,
        "mock_galaxy_covariance_computed": False,
    }
    for item in products:
        file_path = out_dir / "official_polygons" / item["filename"]
        previous = old_records.get(item["filename"])
        if file_path.exists():
            sha, size = digest_file(file_path)
            if size < 16 or size > LIMIT_PER_POLYGON:
                raise ValueError("Invalid cached official polygon byte length")
            print("EBOSS_OFFICIAL_POLYGON_CACHE", item["filename"], sha, size,
                  flush=True)
        else:
            sha, size = stream_pinned(item["source_url"], file_path,
                                      item["root"], item["filename"], timeout)
            print("EBOSS_OFFICIAL_POLYGON_FETCHED", item["filename"], sha, size,
                  flush=True)
        if previous and (
                previous["sha256"] != sha or previous["bytes"] != size
                or previous["source_url"] != item["source_url"]
                or previous["role"] != item["role"]):
            raise ValueError("Pinned official polygon bytes or role changed: "
                             + item["filename"])
        report["products"].append({
            "role": item["role"], "filename": item["filename"],
            "source_url": item["source_url"], "sha256": sha, "bytes": size,
        })
        atomic_json(output, report)
    expected = 3 if only_elg_extra else EXPECTED_COUNT
    if len(report["products"]) != expected:
        raise ValueError("Published fixed polygon SHA inventory incomplete")
    report["status"] = (
        "official_elg_extra_polygon_bytes_pinned_only"
        if only_elg_extra else "official_lrg_elg_polygon_bytes_pinned_only"
    )
    report["note"] = (
        "Byte provenance only. Three extra ELG polygons are pinned in the "
        "limited run; the full run requires eleven official source polygons. "
        "This does not implement any official angular mask, ELG brick "
        "maskbits, LRG veto union, completeness or exact pair selection."
    )
    atomic_json(output, report)
    print("EBOSS_OFFICIAL_POLYGON_SHA_INVENTORY_COMPLETE",
          len(report["products"]), output, flush=True)
    return report


def self_test():
    source = "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/"
    validate_url(source + "ELGmasks/ELG_centerpost.ply",
                 source + "ELGmasks/", "ELG_centerpost.ply")
    for invalid in [source + "../wrong.ply",
                    "https://example.com/polygon.ply",
                    source + "ELGmasks/ELG_centerpost.ply?accept=1"]:
        try:
            validate_url(invalid, source + "ELGmasks/", "ELG_centerpost.ply")
        except ValueError:
            pass
        else:
            raise AssertionError("Unsafe source URL allowed")
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmp:
        p = Path(tmp) / "fake.ply"
        p.write_bytes(b"13 polygons\nsynthetic only\n")
        h, n = digest_file(p)
        assert n == p.stat().st_size
        assert h == hashlib.sha256(p.read_bytes()).hexdigest()
    p = json.loads(PROTOCOL.read_text())
    assert len(p["lrg_published_polygons"]) == 7
    assert len(p["elg_extra_polygons"]) == 3
    assert p["lrg_root_geometry_polygon"].endswith(".ply")
    print("EBOSS_OFFICIAL_POLYGON_SHA_SOURCE_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory-json", type=Path)
    ap.add_argument("--out-dir", type=Path,
                    default=Path("eboss_workspace/official_mask_inventory"))
    ap.add_argument("--timeout", type=float, default=100)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only-elg-extra", action="store_true",
                    help="Pin only the three small required ELG extra-veto polygons; never claim all eleven")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.inventory_json:
        ap.error("--inventory-json is required unless --self-test is set")
    if args.timeout <= 0:
        ap.error("timeout must be positive")
    inv = json.loads(args.inventory_json.read_text(encoding="utf-8"))
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    try:
        audit(protocol, inv, args.out_dir, args.timeout,
              only_elg_extra=args.only_elg_extra)
    except Exception as exc:
        print("EBOSS_OFFICIAL_POLYGON_SHA_INVENTORY_FAILED", type(exc).__name__,
              str(exc), flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
