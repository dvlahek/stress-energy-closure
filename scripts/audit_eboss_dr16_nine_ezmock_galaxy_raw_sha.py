#!/usr/bin/env python3
"""Source-byte-only, exact predeclared nine-ID realistic EZmock DATA inventory.

Reuse four already SHA-pinned 0001 galaxies and acquire 32 exactly fixed
additional mock galaxy gzip files. No gzip decompression, FITS headers,
galaxy/random/observed rows, pair counts or observed odd signal.
Previously completed checkpoints only resume under identical source SHAs.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

from audit_eboss_dr16_ezmock0001_galaxy_raw_sha import (
    ROOT, BASE, ORDER, mock_path, digest, stream_bytes, acquire, load, atomic
)

PROTOCOL = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_protocol_2026-09-26.json"
PASS = "EZMOCK_PREDECLARED_NINE_GALAXY_RAW_SHA_INVENTORY_ONLY"
STOP = "EZMOCK_PREDECLARED_NINE_GALAXY_RAW_SHA_INCOMPLETE_STOP"
NEW_IDS = (125, 250, 375, 500, 625, 750, 875, 1000)
ALL_IDS = (1,) + NEW_IDS
ORDER_KEYS = tuple(cap + "/" + tracer for cap, tracer in ORDER)
NEW_KEYS = tuple(f"{id_:04d}/{cap}/{tracer}" for id_ in NEW_IDS
                 for cap, tracer in ORDER)


def git_blob_sha(raw):
    return hashlib.sha1(
        ("blob " + str(len(raw))).encode("ascii") + bytes((0,)) + raw
    ).hexdigest()


def checked_prior(p):
    base = load(ROOT / p["previously_fixed_window_protocol"])
    random_pilot = load(ROOT / p["parent_ninemock_random_only_snapshot"])
    raw_manifest = load(ROOT / p["original_id0001_sha_only_manifest"])
    odd_manifest = load(ROOT / p["successful_id0001_mock_odd_projection_manifest"])
    original = (ROOT / p["original_id0001_raw_report_archive"]).read_bytes()
    odd_raw = (ROOT / p["successful_id0001_mock_odd_projection_archive"]).read_bytes()
    if (
        tuple(p["previously_fixed_nine_ids"]) != ALL_IDS
        or tuple(base["predeclared_mock_realization_ids"]) != ALL_IDS
        or tuple(random_pilot["mock_ids"]) != ALL_IDS
        or tuple(p["other_eight_ids_in_original_order"]) != NEW_IDS
        or p["caps"] != ["NGC", "SGC"]
        or p["tracers"] != ["eBOSS_LRG", "eBOSS_ELG"]
        or p["role"] != "dat"
        or p["official_release_base_url"] != BASE
        or p["unchanged_existing_frozen_id"] != 1
        or p["strict_source_relpath_template"] !=
           "{tracer}/dat/EZmock_realistic_{tracer}_{cap}_v7_{id:04d}.dat.fits.gz"
        or p["file_count_existing_id0001"] != 4
        or p["file_count_additional_exact_galaxy_sources"] != 32
        or p["file_count_total_expected"] != 36
        or p["max_compressed_bytes_each_additional_source"] != 16777216
        or p["max_combined_additional_compressed_bytes"] != 536870912
        or random_pilot["status"] != "nine_mock_random_only_cross_ls_algebra_complete"
        or random_pilot["observed_odd_data_vector_read"] is not False
        or p["no_mock_galaxy_rows_read"] is not True
        or p["no_mock_random_rows_read"] is not True
        or p["no_observed_galaxy_rows_read"] is not True
        or p["no_observed_odd_data_vector_read"] is not True
        or p["no_new_science_selection_applied"] is not True
        or hashlib.sha256(original).hexdigest() !=
           p["original_id0001_raw_report_sha256"]
        or hashlib.sha256(original).hexdigest() !=
           raw_manifest["uploaded_report_sha256"]
        or len(original) != raw_manifest["uploaded_report_bytes"]
        or git_blob_sha(original) != raw_manifest["uploaded_report_git_blob_sha1"]
        or hashlib.sha256(odd_raw).hexdigest() !=
           p["successful_id0001_mock_odd_projection_sha256"]
        or hashlib.sha256(odd_raw).hexdigest() !=
           odd_manifest["exact_uploaded_SHA256"]
        or len(odd_raw) != odd_manifest["exact_uploaded_bytes"]
        or git_blob_sha(odd_raw) != odd_manifest["exact_uploaded_git_blob_SHA1"]
    ):
        raise ValueError("Previously registered mock-nine cohort or byte-identical 0001 reports changed")
    original_json, odd_json = json.loads(original), json.loads(odd_raw)
    if (
        original_json["status"] != raw_manifest["raw_status"]
        or original_json["mock_realization_id"] != 1
        or original_json["errors"] != []
        or original_json["observed_odd_data_vector_read"] is not False
        or odd_json["status"] != "EZMOCK0001_MOCK_GALAXY_ODD_PROJECTION_DESCRIPTIVE_ONLY"
        or odd_json["errors"] != []
        or odd_json["original_pilot_report_exact_sha256"] !=
           odd_manifest["original_galaxy_cross_ls_report_sha256"]
        or odd_json["source_compressed_SHA_and_prior_pair_histograms_reverified"] is not True
        or odd_json["observed_odd_data_vector_read"] is not False
        or list(original_json["samples"]) != list(ORDER_KEYS)
    ):
        raise ValueError("Successful original source/odd descriptive mock0001 archive is invalid")
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        rec = original_json["samples"][key]
        frozen = raw_manifest["per_source_frozen_from_uploaded_report"][key]
        relative = mock_path(tracer, cap, "dat", 1)
        if (
            rec["source_url"] != BASE + relative
            or rec["returned_url"] != rec["source_url"]
            or rec["first_seen_full_compressed_sha256"] != frozen["full_compressed_sha256"]
            or rec["compressed_bytes"] != frozen["compressed_bytes"]
            or rec["gzip_payload_decompressed"] is not False
            or rec["FITS_headers_or_mock_rows_read"] is not False
            or not rec["local_path"].endswith("/" + relative)
            or frozen["released_filename"] != Path(relative).name
        ):
            raise ValueError("Previously source-pinned 0001 raw galaxy input changed: " + key)
    return original_json, raw_manifest


def existing_0001_sources(p, original, manifest):
    base = load(ROOT / manifest["parent_source_protocol"])
    if base["mock_realization_id"] != 1:
        raise ValueError("Old 0001 original file path protocol changed")
    records = {}
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        relative = mock_path(tracer, cap, "dat", 1)
        dest = ROOT / base["local_quarantine_dir"] / relative
        expected = manifest["per_source_frozen_from_uploaded_report"][key]
        sha, n = digest(dest, p["max_compressed_bytes_each_additional_source"])
        if sha != expected["full_compressed_sha256"] or n != expected["compressed_bytes"]:
            raise ValueError("Previously pinned local mock0001 source bytes changed: " + key)
        records["0001/" + key] = {
            "release_source_url": BASE + relative,
            "full_compressed_sha256_reverified": sha,
            "compressed_bytes": n,
            "from_preexisting_frozen_0001_manifest": True,
            "FITS_headers_or_mock_rows_read": False,
        }
        print("REVERIFIED_EXISTING_0001_GALAXY_SHA", key, n, sha, flush=True)
    return records


def run(p, timeout):
    if timeout <= 0:
        raise ValueError("Source transfer timeout must be positive")
    original, manifest = checked_prior(p)
    protocol_sha = hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    dest = ROOT / p["local_checkpoint_report"]
    if dest.exists():
        result = load(dest)
        if (
            result.get("protocol_sha256") != protocol_sha
            or result.get("status") not in (STOP, PASS)
            or result.get("mock_ids") != list(ALL_IDS)
            or result.get("observed_odd_data_vector_read") is not False
            or result.get("mock_galaxy_or_random_rows_read") is not False
            or not isinstance(result.get("new_galaxy_sources"), dict)
            or not isinstance(result.get("existing_frozen_0001_sources"), dict)
            or not set(result["new_galaxy_sources"]).issubset(NEW_KEYS)
        ):
            raise ValueError("Existing nine-ID source report incompatible; refuse repinning")
        result["status"] = STOP
        result["errors"] = []
    else:
        result = {
            "status": STOP,
            "protocol_sha256": protocol_sha,
            "mock_ids": list(ALL_IDS),
            "existing_frozen_0001_sources": {},
            "new_galaxy_sources": {},
            "mock_galaxy_or_random_rows_read": False,
            "mock_pair_counts_computed": False,
            "observed_galaxy_rows_read": False,
            "observed_odd_data_vector_read": False,
            "new_science_selection_applied": False,
            "not_an_18dimensional_mock_covariance_or_detection": True,
            "errors": [],
        }
    prior_0001 = existing_0001_sources(p, original, manifest)
    if result["existing_frozen_0001_sources"] and (
        result["existing_frozen_0001_sources"] != prior_0001
    ):
        raise ValueError("Previously checkpointed 0001 source identities differ")
    result["existing_frozen_0001_sources"] = prior_0001
    atomic(dest, result)
    sum_additional = 0
    for id_ in NEW_IDS:
        for cap, tracer in ORDER:
            key = f"{id_:04d}/{cap}/{tracer}"
            relative = mock_path(tracer, cap, "dat", id_)
            url = BASE + relative
            path = ROOT / p["local_quarantine_dir"] / relative
            previous = result["new_galaxy_sources"].get(key)
            if previous:
                sha, size = digest(path, p["max_compressed_bytes_each_additional_source"])
                if (
                    previous.get("official_source_url") != url
                    or previous.get("local_path") != str(path.resolve())
                    or previous.get("full_compressed_sha256_first_seen") != sha
                    or previous.get("compressed_bytes") != size
                    or previous.get("gzip_payload_decompressed") is not False
                    or previous.get("FITS_headers_or_mock_rows_read") is not False
                ):
                    raise ValueError("Previously checkpointed new source SHA changed: " + key)
                print("REUSED_NINE_MOCK_GALAXY_SHA", key, size, sha, flush=True)
            else:
                got = acquire(url, path, p["max_compressed_bytes_each_additional_source"], timeout)
                record = {
                    "official_source_url": url,
                    "returned_official_url": got["returned_url"],
                    "local_path": got["local_path"],
                    "compressed_bytes": got["compressed_bytes"],
                    "full_compressed_sha256_first_seen": got["first_seen_full_compressed_sha256"],
                    "HTTP_status": got["HTTP_status"],
                    "HTTP_content_length_if_available": got["HTTP_content_length_if_available"],
                    "HTTP_etag_if_available": got["HTTP_etag_if_available"],
                    "HTTP_last_modified_if_available": got["HTTP_last_modified_if_available"],
                    "gzip_magic_checked": got["gzip_magic_checked"],
                    "gzip_payload_decompressed": False,
                    "FITS_headers_or_mock_rows_read": False,
                    "source_SHA_is_first_seen_not_official_reference": True,
                }
                result["new_galaxy_sources"][key] = record
                atomic(dest, result)
                sha, size = record["full_compressed_sha256_first_seen"], record["compressed_bytes"]
                print("FIRST_SEEN_NINE_MOCK_GALAXY_SHA", key, size, sha, flush=True)
            sum_additional += size
            if sum_additional > p["max_combined_additional_compressed_bytes"]:
                raise ValueError("Frozen total compressed-source byte cap exceeded")
    if set(result["new_galaxy_sources"]) != set(NEW_KEYS):
        raise ValueError("Not all 32 additional source filenames were SHA-pinned")
    if set(result["existing_frozen_0001_sources"]) != set("0001/" + k for k in ORDER_KEYS):
        raise ValueError("Not all four old mock0001 source filenames were SHA-reverified")
    result["new_additional_compressed_bytes"] = sum_additional
    result["previously_frozen_0001_compressed_bytes"] = sum(
        x["compressed_bytes"] for x in prior_0001.values()
    )
    result["status"] = PASS
    result["errors"] = []
    atomic(dest, result)
    return result


def self_test(p):
    old, pins = checked_prior(p)  # CI only archived metadata, no real FITS.
    if (len(old["samples"]) != 4 or len(pins["per_source_frozen_from_uploaded_report"]) != 4
        or len(NEW_KEYS) != 32 or len(set(NEW_KEYS)) != 32):
        raise AssertionError("Frozen nine-ID source count/order changed")
    fake = bytes((31,139,8)) + b"synthetic raw mock source only"
    out = io.BytesIO()
    sha, n = stream_bytes(io.BytesIO(fake), out, 256, len(fake))
    if n != len(fake) or out.getvalue() != fake or sha != hashlib.sha256(fake).hexdigest():
        raise AssertionError("Synthetic full source SHA stream changed")
    for broken, declared, cap in (
        (fake[:-1], len(fake), 256),
        (b"BAD" + fake[3:], len(fake), 256),
        (fake + b"x", len(fake), 256),
        (fake, len(fake), len(fake)-1),
    ):
        try:
            stream_bytes(io.BytesIO(broken), io.BytesIO(), cap, declared)
        except ValueError:
            pass
        else:
            raise AssertionError("Damaged synthetic source or size cap accepted")
    print("EBOSS_NINE_EZMOCK_GALAXY_RAW_SHA_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--timeout", type=float, default=120)
    args = ap.parse_args()
    p = load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    out = ROOT / p["local_checkpoint_report"]
    try:
        result = run(p, args.timeout)
    except Exception as exc:
        if out.exists():
            result = load(out)
            result["status"] = STOP
            result["errors"] = [str(exc)]
        else:
            result = {
                "status": STOP,
                "errors": [str(exc)],
                "mock_galaxy_or_random_rows_read": False,
                "observed_galaxy_rows_read": False,
                "observed_odd_data_vector_read": False,
                "new_science_selection_applied": False,
            }
        atomic(out, result)
        print("EBOSS_NINE_EZMOCK_GALAXY_SOURCE", result["status"], flush=True)
        print("REPORT", out, flush=True)
        print("ERRORS", *result["errors"], sep="\n", flush=True)
        return 2
    print("EBOSS_NINE_EZMOCK_GALAXY_SOURCE", result["status"], flush=True)
    print("REPORT", out, flush=True)
    print("FOUR_FROZEN_0001_SOURCES_REVERIFIED", len(result["existing_frozen_0001_sources"]), flush=True)
    print("NEW_SOURCE_COUNT", len(result["new_galaxy_sources"]), flush=True)
    print("ADDITIONAL_COMPRESSED_BYTES", result["new_additional_compressed_bytes"], flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
