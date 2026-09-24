#!/usr/bin/env python3
"""Compare released eBOSS LRG and ELG EZmock filenames by realization ID.

This validates published paths and matching IDs only. It does not open mock
FITS files, validate common physical initial conditions or estimate covariance.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

MOCK = re.compile(
    r"^v1_0_0/(?P<mode>complete|realistic)/"
    r"(?P<tracer>eBOSS_LRG|eBOSS_ELG)/(?P<role>dat|ran|ran\.shuf)/"
    r"EZmock_(?P=mode)_(?P=tracer)_(?P<cap>NGC|SGC)_"
    r"(?:(?P<version>v7)_)?(?P<id>\d{4})\."
    r"(?P<suffix>dat|ran|ran\.shuf)\.fits\.gz$"
)


def audit(index: dict, expected: int) -> dict:
    groups = defaultdict(set)
    malformed = []
    for entry in index["potential_data_files"]:
        path = entry["path"]
        if "eBOSS_LRG/" not in path and "eBOSS_ELG/" not in path:
            continue
        m = MOCK.fullmatch(path)
        if m is None:
            # Complete mocks also have shared, non-realization random files.
            if re.fullmatch(
                r"v1_0_0/complete/eBOSS_(?:LRG|ELG)/ran/"
                r"EZmock_complete_eBOSS_(?:LRG|ELG)_(?:NGC|SGC)\.ran\.fits\.gz",
                path,
            ):
                continue
            malformed.append(path)
            continue
        item = m.groupdict()
        mode, role = item["mode"], item["role"]
        version = item["version"]
        if (mode == "realistic" and (version != "v7" or role not in ("dat", "ran"))
                or mode == "complete" and (version is not None or role not in ("dat", "ran.shuf"))
                or role != item["suffix"]):
            malformed.append(path)
            continue
        key = (mode, item["cap"], item["tracer"], role)
        realization = int(item["id"])
        if realization in groups[key]:
            malformed.append("duplicate realization: " + path)
            continue
        groups[key].add(realization)

    expected_ids = set(range(1, expected + 1))
    by_mode_cap = {}
    filename_complete = not index.get("errors") and not index.get("truncated") and not malformed
    for mode in ("complete", "realistic"):
        for cap in ("NGC", "SGC"):
            random_role = "ran.shuf" if mode == "complete" else "ran"
            selections = {}
            streams = []
            for tracer in ("eBOSS_LRG", "eBOSS_ELG"):
                for role in ("dat", random_role):
                    ids = groups[(mode, cap, tracer, role)]
                    streams.append(ids)
                    selections[f"{tracer}_{role}"] = {
                        "count": len(ids),
                        "missing_ids": sorted(expected_ids - ids)[:20],
                        "unexpected_ids": sorted(ids - expected_ids)[:20],
                        "complete": ids == expected_ids,
                    }
                    filename_complete &= ids == expected_ids
            matched = set.intersection(*streams)
            by_mode_cap[f"{mode}_{cap}"] = {
                "tracer_and_random_streams": selections,
                "joint_realization_ids_count": len(matched),
                "joint_expected_ids_complete": matched == expected_ids,
                "first_matched_id": min(matched) if matched else None,
                "last_matched_id": max(matched) if matched else None,
            }
            filename_complete &= matched == expected_ids

    listed_paths = sorted(item["path"] for item in index["potential_data_files"])
    fingerprint = hashlib.sha256(("\n".join(listed_paths) + "\n").encode("utf-8")).hexdigest()
    return {
        "study": "eBOSS DR16 EZmock LRG–ELG filename-level pairing",
        "release_url": index["release_url"],
        "release_version": "v1_0_0",
        "expected_realization_ids": [1, expected],
        "index_file_count": len(listed_paths),
        "listing_path_sha256": fingerprint,
        "status": "filename_sets_complete" if filename_complete else "incomplete",
        "by_mode_cap": by_mode_cap,
        "malformed_or_unexpected_filename_examples": malformed[:20],
        "malformed_or_unexpected_count": len(malformed),
        "fits_contents_inspected": False,
        "common_initial_conditions_verified_from_files": False,
        "selection_and_weight_compatibility_verified": False,
        "cross_tracer_covariance_calculated": False,
        "note": (
            "The complete and realistic releases contain matching LRG and ELG "
            "realization identifiers and matching tracer-specific random files. "
            "Common initial conditions are documented in the release literature "
            "but have not been verified by this filename-level audit."
        ),
    }


def self_test() -> None:
    examples = (
        "v1_0_0/realistic/eBOSS_ELG/dat/"
        "EZmock_realistic_eBOSS_ELG_NGC_v7_0001.dat.fits.gz",
        "v1_0_0/complete/eBOSS_LRG/ran.shuf/"
        "EZmock_complete_eBOSS_LRG_SGC_1000.ran.shuf.fits.gz",
    )
    assert all(MOCK.fullmatch(p) for p in examples)
    assert MOCK.fullmatch(examples[0]).group("version") == "v7"
    assert MOCK.fullmatch(examples[1]).group("version") is None
    print("EZmock filename parser test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--index", default="eboss_workspace/mock_release_index.json")
    parser.add_argument("--out", default="eboss_workspace/mock_pair_summary.json")
    parser.add_argument("--expected-count", type=int, default=1000)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.expected_count <= 0:
        parser.error("expected-count must be positive")
    summary = audit(json.loads(Path(args.index).read_text(encoding="utf-8")),
                    args.expected_count)
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("Mock pair summary:", path, "status:", summary["status"])
    for key, result in summary["by_mode_cap"].items():
        print(key, "matched IDs:", result["joint_realization_ids_count"],
              "complete:", result["joint_expected_ids_complete"])
    return 0 if summary["status"] == "filename_sets_complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
