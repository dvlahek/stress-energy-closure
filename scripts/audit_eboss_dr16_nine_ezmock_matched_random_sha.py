#!/usr/bin/env python3
"""SHA-only 36 matched same-ID realistic EZmock random gzip source audit.

Reference SHA256s are from an archived PRIOR successful 2026-09-24
nine-ID ensemble artifact, independently cross-checked on the 12 older
three-ID random pins. Reuse known caches by FULL byte identity; optional
exact official downloads must match the OLD SHA. No FITS decompression,
mock/observed row, pair, covariance or observed odd measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from audit_eboss_dr16_ezmock0001_galaxy_raw_sha import (
    ROOT, BASE, digest, stream_bytes, acquire, load, atomic, mock_path
)

PROTOCOL = ROOT / "source_data/eboss_dr16_nine_ezmock_matched_random_sha_protocol_2026-09-26.json"
PASS = "EZMOCK_PREDECLARED_NINE_MATCHED_RANDOM_FULL_SHA_REVERIFIED_ONLY"
STOP = "EZMOCK_PREDECLARED_NINE_MATCHED_RANDOM_FULL_SHA_INCOMPLETE_STOP"
IDS = (1,125,250,375,500,625,750,875,1000)
CAPS = ("NGC","SGC")
TRACERS = ("eBOSS_LRG","eBOSS_ELG")
ALL = tuple((mid,cap,tracer) for mid in IDS for cap in CAPS for tracer in TRACERS)


def git_blob_sha(raw):
    return hashlib.sha1(("blob " + str(len(raw))).encode("ascii")
                        + bytes([0]) + raw).hexdigest()


def source_key(mid, cap, tracer):
    return f"{mid:04d}/{cap}/{tracer}"


def source_rel(mid, cap, tracer):
    return mock_path(tracer, cap, "ran", mid)


def exact_sha(path, expected, limit):
    if path.is_symlink() or not path.is_file():
        raise ValueError("Expected a regular, non-symlink mock random gzip: " + str(path))
    got, n = digest(path, limit)
    if got != expected:
        raise ValueError("COMPRESSED_SHA_MISMATCH_NO_REPIN " + str(path)
                         + " expected=" + expected + " actual=" + got)
    return n


def source_preflight(p, *, require_local_report):
    manifest_path = ROOT / p["galaxy_uploaded_manifest"]
    reference_path = ROOT / p["random_sha_reference"]
    archive_path = ROOT / p["galaxy_archived_local_report"]
    manifest_raw = manifest_path.read_bytes()
    reference_raw = reference_path.read_bytes()
    archived_raw = archive_path.read_bytes()
    galaxy_manifest = json.loads(manifest_raw)
    reference = json.loads(reference_raw)
    report = json.loads(archived_raw)
    original_protocol = load(ROOT / p["previous_nine_mock_window_protocol"])
    old_ref = load(ROOT / p["old_three_mock_random_source_reference"])
    expected_keys = tuple(source_key(*entry) for entry in ALL)
    if (
        git_blob_sha(manifest_raw) != p["galaxy_uploaded_manifest_git_blob_sha1"]
        or git_blob_sha(reference_raw) != p["random_sha_reference_git_blob_sha1"]
        or p["reference_original_workflow_run"] != 36017670812
        or p["reference_original_artifact_id"] != 10815013322
        or reference["source_workflow_run_id"] != p["reference_original_workflow_run"]
        or reference["source_workflow_artifact_id"] != p["reference_original_artifact_id"]
        or reference["source_json_sha256"] !=
           p["reference_original_artifact_member_sha256"]
        or reference["source_json_revision_commit"] !=
           "054005edc6ca193176b129e1951e4bd3ce8751a7"
        or reference["source_sha_count"] != 36
        or tuple(reference["fixed_mock_ids"]) != IDS
        or tuple(reference["fixed_caps"]) != CAPS
        or tuple(reference["fixed_tracers"]) != TRACERS
        or tuple(reference["full_random_gzip_sha256_by_id_cap_tracer"]) != expected_keys
        or tuple(p["previously_fixed_ids"]) != IDS
        or tuple(original_protocol["predeclared_mock_realization_ids"]) != IDS
        or p["caps"] != list(CAPS)
        or p["tracers"] != list(TRACERS)
        or p["role"] != "ran"
        or p["official_base_url"] != BASE
        or p["filename_rule"] !=
           "{tracer}/ran/EZmock_realistic_{tracer}_{cap}_v7_{id:04d}.ran.fits.gz"
        or p["max_compressed_bytes_each"] != 128 * 1024 * 1024
        or p["max_total_compressed_bytes"] != 4 * 1024 * 1024 * 1024
        or p["expected_total_random_source_files"] != 36
        or p["expected_prior_3_id_sha_matches"] != 12
        or p["no_gzip_decompression_or_FITS_header_row_access"] is not True
        or p["no_mock_galaxy_random_or_observed_rows_read"] is not True
        or p["no_observed_odd_vector_read"] is not True
        or p["no_new_sky_mask_or_science_cuts"] is not True
        or not reference["reference_source_NOT_latest_local_random_binaries"]
        or p["galaxy_uploaded_report_sha256"] !=
           "2c6b55bf5611dc71449ecbbfc2dedf545f2d22250147ed08f6a072536b50f06c"
        or p["galaxy_uploaded_report_git_blob_sha1"] !=
           "43be7aae5b06580867bee4fff08292ebeee468bf"
        or hashlib.sha256(archived_raw).hexdigest() != p["galaxy_uploaded_report_sha256"]
        or len(archived_raw) != p["galaxy_uploaded_report_bytes"]
        or git_blob_sha(archived_raw) != p["galaxy_uploaded_report_git_blob_sha1"]
        or galaxy_manifest["uploaded_exact_sha256"] != p["galaxy_uploaded_report_sha256"]
        or galaxy_manifest["archive_git_blob_sha1"] != p["galaxy_uploaded_report_git_blob_sha1"]
        or report["status"] != "EZMOCK_PREDECLARED_NINE_GALAXY_RAW_SHA_INVENTORY_ONLY"
        or report["errors"] != []
        or report["mock_galaxy_or_random_rows_read"] is not False
        or report["observed_odd_data_vector_read"] is not False
        or tuple(report["mock_ids"]) != IDS
        or len(report["new_galaxy_sources"]) != 32
        or len(report["existing_frozen_0001_sources"]) != 4
        or galaxy_manifest["count_total"] != 36
        or galaxy_manifest["source_protocol_sha256"] != report["protocol_sha256"]
    ):
        raise ValueError("Original nine-ID source report/manifest or 2026-09-24 random SHA provenance changed")
    old = {
        item["id"] + "/" + item["cap"] + "/eBOSS_" + item["tracer"]:
            item["compressed_file_sha256"]
        for item in old_ref["mock_random_catalogues"]
    }
    if (
        len(old) != 12
        or reference["original_three_ID_source_SHA_checks_against_prior_committed_audit"] != 12
        or any(reference["full_random_gzip_sha256_by_id_cap_tracer"].get(k) != sha
               for k,sha in old.items())
    ):
        raise ValueError("Prior independent 0001/0500/1000 random SHA validation changed")
    for mid,cap,tracer in ALL:
        value = reference["full_random_gzip_sha256_by_id_cap_tracer"][
            source_key(mid,cap,tracer)]
        if not isinstance(value,str) or len(value) != 64 or any(
            ch not in "0123456789abcdef" for ch in value
        ):
            raise ValueError("Invalid archived previously pinned 36-random SHA256")
    # CI synthetic preflight uses archived report alone; the real source
    # audit MUST separately verify local uploaded report's exact raw bytes.
    if require_local_report and (
        ROOT / p["galaxy_uploaded_report_local"]
    ).read_bytes() != archived_raw:
        raise ValueError("Local user uploaded galaxy source report differs byte-for-byte from Git archive")
    return reference, report


def cache_candidates(p, mid, cap, tracer):
    rel = Path(source_rel(mid,cap,tracer))
    dest = ROOT / p["new_quarantine_dir"] / rel
    candidates = []
    for directory in p["cache_lookup_relative_dirs"]:
        root = ROOT / directory
        for item in (root / rel.name, root / rel):
            if item not in candidates:
                candidates.append(item)
    if dest not in candidates:
        raise ValueError("Approved quarantine source path omitted from frozen cache discovery")
    return candidates, dest


def source_record(p, ref, mid, cap, tracer, *, download_missing, timeout):
    key = source_key(mid,cap,tracer)
    expected = ref["full_random_gzip_sha256_by_id_cap_tracer"][key]
    rel = source_rel(mid,cap,tracer)
    url = BASE + rel
    candidates, dest = cache_candidates(p, mid,cap,tracer)
    approved_existing = [path for path in candidates if path.exists() or path.is_symlink()]
    if approved_existing:
        source = approved_existing[0]
        n = exact_sha(source, expected, p["max_compressed_bytes_each"])
        mode = "verified_existing_cache"
        http = None
        print("REVERIFIED_MATCHED_NINE_RANDOM_SHA", key, n, expected, flush=True)
    else:
        if not download_missing:
            raise FileNotFoundError(
                "Approved matched mock random not cached: " + key
                + ". Re-run with --download-missing to fetch ONLY its unchanged pinned official URL."
            )
        if dest.exists() or dest.is_symlink():
            raise ValueError("Refusing overwrite of an existing unchecked random quarantine file")
        fetched = acquire(url, dest, p["max_compressed_bytes_each"], timeout)
        if fetched["first_seen_full_compressed_sha256"] != expected:
            dest.unlink(missing_ok=True)
            raise ValueError("OFFICIAL_RANDOM_FULL_SHA_MISMATCH_NO_REPIN " + key
                             + " prior_expected=" + expected
                             + " new_actual=" + fetched["first_seen_full_compressed_sha256"])
        if fetched["returned_url"] != url or fetched["HTTP_status"] != 200:
            dest.unlink(missing_ok=True)
            raise ValueError("Unexpected official mock random redirect or HTTP response")
        source, n, mode = dest, fetched["compressed_bytes"], "downloaded_same_prior_pinned_URL"
        http = {
            "HTTP_status":fetched["HTTP_status"],
            "HTTP_content_length_if_available":fetched["HTTP_content_length_if_available"],
            "HTTP_etag_if_available":fetched["HTTP_etag_if_available"],
            "HTTP_last_modified_if_available":fetched["HTTP_last_modified_if_available"],
        }
        print("DOWNLOAD_MATCHED_NINE_RANDOM_SHA_OK", key, n, expected, flush=True)
    if n <= 3:
        raise ValueError("Unexpected empty official mock random compressed source")
    return {
        "source_key":key,
        "official_source_url":url,
        "local_verified_source_path":str(source.resolve()),
        "relative_released_filename":str(rel),
        "source_mode":mode,
        "compressed_bytes":n,
        "prior_20260924_artifact_expected_sha256":expected,
        "fully_reverified_compressed_sha256":expected,
        "sha_matches_original_artifact_reference":True,
        "gzip_magic_checked":True,
        "gzip_decompressed_or_FITS_rows_read":False,
        "newly_downloaded_HTTP_details_if_any":http,
    }


def run(p, *, download_missing, timeout):
    if timeout <= 0:
        raise ValueError("Source transfer timeout must be positive")
    ref, _ = source_preflight(p,require_local_report=True)
    protocol_sha = hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    out = ROOT / p["local_report"]
    if out.exists():
        result = load(out)
        if (
            result.get("protocol_sha256") != protocol_sha
            or result.get("status") not in (STOP,PASS)
            or result.get("fixed_mock_ids") != list(IDS)
            or result.get("observed_odd_data_vector_read") is not False
            or result.get("no_mock_galaxy_random_or_observed_rows_read") is not True
            or not isinstance(result.get("verified_random_sources"),dict)
            or not set(result["verified_random_sources"]).issubset(
                {source_key(*t) for t in ALL}
            )
        ):
            raise ValueError("Existing random SHA checkpoint is inconsistent; refuse silent repinning")
        result["status"] = STOP
        result["errors"] = []
    else:
        result = {
            "status":STOP,
            "protocol_sha256":protocol_sha,
            "prior_20260924_source_json_sha256":ref["source_json_sha256"],
            "galaxy_uploaded_report_sha256":p["galaxy_uploaded_report_sha256"],
            "fixed_mock_ids":list(IDS),
            "caps":list(CAPS),
            "tracers":list(TRACERS),
            "verified_random_sources":{},
            "no_mock_galaxy_random_or_observed_rows_read":True,
            "mock_pairs_computed":False,
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
            "joint_mock_covariance_computed":False,
            "errors":[],
        }
    total = 0
    for mid,cap,tracer in ALL:
        key = source_key(mid,cap,tracer)
        prev = result["verified_random_sources"].get(key)
        if prev is not None:
            got_path = Path(prev["local_verified_source_path"])
            expected_rel = Path(source_rel(mid,cap,tracer))
            if (
                prev.get("source_key") != key
                or prev.get("official_source_url") != BASE+str(expected_rel)
                or prev.get("relative_released_filename") != str(expected_rel)
                or prev.get("prior_20260924_artifact_expected_sha256") !=
                   ref["full_random_gzip_sha256_by_id_cap_tracer"][key]
                or prev.get("fully_reverified_compressed_sha256") !=
                   ref["full_random_gzip_sha256_by_id_cap_tracer"][key]
                or prev.get("gzip_decompressed_or_FITS_rows_read") is not False
                or not str(got_path).endswith("/"+expected_rel.name)
            ):
                raise ValueError("Previously verified random source checkpoint mutated: "+key)
            if not any(got_path.resolve() == candidate.resolve()
                       for candidate in cache_candidates(p,mid,cap,tracer)[0]):
                raise ValueError("Previously verified random cache moved outside approved locations")
            n = exact_sha(got_path,ref["full_random_gzip_sha256_by_id_cap_tracer"][key],
                          p["max_compressed_bytes_each"])
            if prev["compressed_bytes"] != n:
                raise ValueError("Previously verified random gzip source size changed: "+key)
            print("REUSED_PINNED_MATCHED_NINE_RANDOM_SHA", key,n,
                  prev["fully_reverified_compressed_sha256"], flush=True)
        else:
            row = source_record(p,ref,mid,cap,tracer,
                                download_missing=download_missing,timeout=timeout)
            result["verified_random_sources"][key] = row
            atomic(out,result)
            n = row["compressed_bytes"]
        total += n
        if total > p["max_total_compressed_bytes"]:
            raise ValueError("Previously frozen 36 mock random compressed source total exceeded")
    if (
        len(result["verified_random_sources"]) != 36
        or set(result["verified_random_sources"]) != {source_key(*t) for t in ALL}
    ):
        raise ValueError("Full 36 source matched random audit incomplete")
    result["total_full_compressed_random_bytes"] = total
    result["complete_random_source_count"] = 36
    result["status"] = PASS
    result["errors"] = []
    atomic(out,result)
    return result


def self_test(p):
    ref, original = source_preflight(p,require_local_report=False)
    if (
        len(ref["full_random_gzip_sha256_by_id_cap_tracer"]) != 36
        or original["new_additional_compressed_bytes"] != 64937949
        or original["previously_frozen_0001_compressed_bytes"] != 8307621
    ):
        raise AssertionError("Exact uploaded 36 galaxy and prior 36 random reference not fixed")
    with TemporaryDirectory() as td:
        path = Path(td)/"synthetic.ran.fits.gz"
        raw = bytes((31,139,8)) + b"only synthetic SHA reference, no FITS"
        expected = hashlib.sha256(raw).hexdigest()
        path.write_bytes(raw)
        if exact_sha(path,expected,len(raw)) != len(raw):
            raise AssertionError("Valid synthetic compressed input SHA not verified")
        for changed,limit in ((expected[:-1]+"0",len(raw)),(expected,len(raw)-1)):
            try:
                exact_sha(path,changed,limit)
            except ValueError:
                pass
            else:
                raise AssertionError("Corrupted prior random SHA/size gate accepted")
        bad = b"BAD"+raw[3:]
        try:
            stream_bytes(io.BytesIO(bad),io.BytesIO(),1024,len(bad))
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid full gzip transfer accepted")
    print("EBOSS_NINE_MATCHED_EZMOCK_RANDOM_SHA_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--download-missing",action="store_true",
                    help="Opt in to large, bounded official random gzip downloads; only prior 2026-09-24 SHA refs accepted")
    ap.add_argument("--timeout",type=float,default=120)
    args=ap.parse_args()
    p=load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    out=ROOT/p["local_report"]
    try:
        result=run(p,download_missing=args.download_missing,timeout=args.timeout)
    except Exception as exc:
        if out.exists():
            result=load(out)
            result["status"]=STOP
            result["errors"]=[str(exc)]
        else:
            result={
                "status":STOP,"errors":[str(exc)],
                "no_mock_galaxy_random_or_observed_rows_read":True,
                "observed_odd_data_vector_read":False,
                "new_science_selection_applied":False,
            }
        atomic(out,result)
        print("EBOSS_NINE_MATCHED_EZMOCK_RANDOM_SHA",result["status"],flush=True)
        print("REPORT",out,flush=True)
        print("ERRORS",*result["errors"],sep="\n",flush=True)
        return 2
    print("EBOSS_NINE_MATCHED_EZMOCK_RANDOM_SHA",result["status"],flush=True)
    print("REPORT",out,flush=True)
    print("VERIFIED_RANDOM_SOURCE_COUNT",result["complete_random_source_count"],flush=True)
    print("TOTAL_COMPRESSED_RANDOM_BYTES",result["total_full_compressed_random_bytes"],flush=True)
    print("OBSERVED_ODD_DATA_READ",result["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
