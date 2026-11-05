#!/usr/bin/env python3
"""Discover *candidate* official DR16 ELG bit-8 reference files by filename.

ONLY read the SHA-pinned public DR16 root directory HTML index; never open
FITS files, catalogues, mask images, coordinates, or odd measurements.
A filename or post-veto clustering catalogue is NOT a verified reference.
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
from pathlib import Path

from audit_eboss_dr16_official_mask_sources import (
    MAX_INDEX_BYTES,
    accept_root,
    fetch_limited,
    listing_names,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json"


def load_protocol():
    p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if (
        p["official_root"]
        != "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/"
        or len(p["root_listing_sha256"]) != 64
        or p["byte_cap"] != MAX_INDEX_BYTES
        or p["official_reference_found"] is not False
        or p["official_bit8_applied"] is not False
        or p["observed_odd_data_vector_read"] is not False
    ):
        raise ValueError("Unexpected frozen reference-discovery protocol")
    accept_root(p["official_root"])
    return p


def classify(entries, p):
    out = []
    for name, (url, is_dir) in sorted(entries.items()):
        if not any(x in name.lower() for x in p["candidate_name_tokens"]):
            continue
        if is_dir:
            kind = "unverified_directory_name_only"
        elif any(
            fnmatch.fnmatchcase(name, pat)
            for pat in p["known_nonreference_globs"]
        ):
            kind = "known_postselection_clustering_not_bit8_reference"
        else:
            kind = "unverified_file_name_only"
        out.append(
            {"filename": name, "source_url": url,
             "is_directory": bool(is_dir), "classification": kind}
        )
    return out


def evaluate(raw, p):
    if len(raw) > p["byte_cap"]:
        raise ValueError("Official root index exceeds byte cap")
    digest = hashlib.sha256(raw).hexdigest()
    if digest != p["root_listing_sha256"]:
        return {
            "status": "PINNED_OFFICIAL_ROOT_INDEX_SHA_MISMATCH_STOP",
            "actual_root_index_sha256": digest,
            "expected_root_index_sha256": p["root_listing_sha256"],
            "source_bytes": len(raw),
            "candidate_files": [],
            "errors": ["No reference candidates accepted; do not silently repin the official root index"],
            "official_pre_veto_bit8_reference_verified": False,
            "observed_odd_data_vector_read": False,
        }
    html = raw.decode("utf-8", errors="strict")
    entries = listing_names(html, p["official_root"])
    if not entries:
        raise ValueError("Empty or unparsable official root index")
    candidates = classify(entries, p)
    return {
        "status": "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED",
        "official_root": p["official_root"],
        "root_index_sha256": digest,
        "source_bytes": len(raw),
        "safe_index_entries": len(entries),
        "candidate_files": candidates,
        "candidate_filename_count": len(candidates),
        "official_pre_veto_bit8_reference_verified": False,
        "official_bit8_applied": False,
        "physical_LRG_ELG_pair_window_certified": False,
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_random_positions_read": False,
        "mock_random_positions_read": False,
        "observed_odd_data_vector_read": False,
        "limitations": (
            "Filename-only inventory; no FITS headers/rows or production bitmap read. "
            "The official reference gate requires authenticated pre-veto bit-coded "
            "output, exact release and SHA, and discriminating positive/negative positions."
        ),
        "errors": [],
    }


def self_test():
    p = load_protocol()
    base = p["official_root"]
    fake = (
        '<a href="../">parent</a>'
        '<a href="eBOSS_ELG_full_ALLdata-vDR16.fits">full</a>'
        '<a href="eBOSS_ELG_clustering_data-NGC-vDR16.fits">data</a>'
        '<a href="ELGmasks/">masks</a>'
        '<a href="https://evil.example.invalid/bit8.fits">evil</a>'
        '<a href="../../escape.fits">escape</a>'
    )
    entries = listing_names(fake, base)
    classes = {i["filename"]: i["classification"] for i in classify(entries, p)}
    assert classes["eBOSS_ELG_full_ALLdata-vDR16.fits"] == "unverified_file_name_only"
    assert classes["eBOSS_ELG_clustering_data-NGC-vDR16.fits"] == (
        "known_postselection_clustering_not_bit8_reference"
    )
    assert classes["ELGmasks"] == "unverified_directory_name_only"
    assert "bit8.fits" not in classes and "escape.fits" not in classes
    print("EBOSS_ELG_BIT8_OFFICIAL_REFERENCE_INDEX_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--index-html", type=Path,
                    help="Optional already downloaded official root HTML; still SHA verified")
    ap.add_argument("--out", type=Path, default=Path(
        "eboss_workspace/official_mask_inventory/bit8_reference_index.json"
    ))
    ap.add_argument("--timeout", type=float, default=90)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    p = load_protocol()
    if args.timeout <= 0:
        ap.error("timeout must be positive")
    try:
        if args.index_html:
            with args.index_html.open("rb") as f:
                raw = f.read(p["byte_cap"] + 1)
        else:
            raw, _, _ = fetch_limited(
                p["official_root"], max_bytes=p["byte_cap"],
                root_url=p["official_root"], timeout=args.timeout,
            )
        report = evaluate(raw, p)
    except Exception as exc:
        report = {
            "status": "OFFICIAL_REFERENCE_INDEX_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "official_pre_veto_bit8_reference_verified": False,
            "observed_odd_data_vector_read": False,
        }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.out)
    print("EBOSS_ELG_BIT8_REFERENCE_INDEX", report["status"], flush=True)
    print("REPORT", args.out, flush=True)
    for item in report.get("candidate_files", []):
        print("CANDIDATE", item["classification"], item["filename"], flush=True)
    if report.get("errors"):
        print("ERRORS", *report["errors"], sep="\n", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
