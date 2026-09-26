#!/usr/bin/env python3
"""SHA-only preflight for four predeclared eBOSS realistic EZmock0001 galaxies.

Downloads only exact released LRG/ELG mock DATA gzip bytes in NGC/SGC.
No FITS decompression, headers, galaxy positions, redshift, weights,
random-catalogue rows, pair counts or observed odd-sector measurements.
First-seen source SHAs must be archived in a separate protocol BEFORE
the next pilot is permitted to open galaxy rows.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.error import URLError
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

from inspect_eboss_dr16_mock_headers import BASE, mock_path

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_raw_bytes_protocol_2026-09-26.json"
PASS = "EZMOCK0001_FOUR_GALAXY_RAW_GZIP_SOURCES_FIRST_SEEN_SHA_ONLY"
PARTIAL = "EZMOCK0001_GALAXY_RAW_SOURCE_INCOMPLETE_STOP"
ORDER = ((cap, tracer) for cap in ("NGC", "SGC")
         for tracer in ("eBOSS_LRG", "eBOSS_ELG"))
GZIP_PREFIX = bytes((31, 139, 8))


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(path, cap):
    sha, n, prefix = hashlib.sha256(), 0, b""
    with path.open("rb") as fd:
        while True:
            block = fd.read(1024 * 1024)
            if not block:
                break
            if n + len(block) > cap:
                raise ValueError("Compressed file exceeds frozen cap: " + path.name)
            prefix += block[:max(0, 3 - len(prefix))]
            sha.update(block)
            n += len(block)
    if n <= 3 or prefix != GZIP_PREFIX:
        raise ValueError("Gzip prefix or source length invalid: " + path.name)
    return sha.hexdigest(), n


def stream_bytes(src, dst, cap, declared=None):
    if declared is not None and (declared <= 3 or declared > cap):
        raise ValueError("Official source declared a disallowed file length")
    sha, n, prefix = hashlib.sha256(), 0, b""
    while True:
        block = src.read(1024 * 1024)
        if not block:
            break
        if n + len(block) > cap or (
            declared is not None and n + len(block) > declared
        ):
            raise ValueError("Official source byte cap or Content-Length exceeded")
        prefix += block[:max(0, 3 - len(prefix))]
        dst.write(block)
        sha.update(block)
        n += len(block)
    if n <= 3 or prefix != GZIP_PREFIX or (
        declared is not None and n != declared
    ):
        raise ValueError("Official source truncated, gzip magic or Content-Length mismatch")
    return sha.hexdigest(), n


def preflight(p):
    win = load(ROOT / p["parent_fixed_nine_mock_window_protocol"])
    header = load(ROOT / p["parent_previous_mock_header_audit"])
    sample = load(ROOT / p["parent_previous_3id_mock_data_selection_audit"])
    ran = load(ROOT / p["parent_pinned_0001_random_source_audit"])
    random_algebra = load(ROOT / p["parent_completed_nine_mock_random_only_estimator_audit"])
    if (
        p["mock_realization_id"] != 1
        or p["prior_ensemble_ids"] != win["predeclared_mock_realization_ids"]
        or p["prior_ensemble_ids"] != random_algebra["mock_ids"]
        or p["caps"] != ["NGC", "SGC"]
        or p["tracers"] != ["eBOSS_LRG", "eBOSS_ELG"]
        or p["role"] != "dat"
        or p["official_base_url"] != BASE
        or header["status"] != "sample_headers_validated"
        or sample["status"] != "sample_candidate_selection_validated_with_numerical_zero_convention"
        or 1 not in sample["mock_realization_ids"]
        or 1 not in header["realization_ids"]
        or random_algebra["status"] != "nine_mock_random_only_cross_ls_algebra_complete"
        or ran["status"] != "mock_random_sample_selection_compatible"
        or p["max_compressed_bytes_per_file"] != 16 * 1024 * 1024
        or p["max_total_compressed_bytes"] != 64 * 1024 * 1024
        or p["no_current_FITS_header_or_mock_galaxy_row_read"] is not True
        or p["no_observed_galaxy_or_observed_odd_read"] is not True
        or p["no_new_mask_or_cut"] is not True
        or header["data_header_row_counts"]["0001"] != {
            "LRG_NGC": 129262, "LRG_SGC": 82921,
            "ELG_NGC": 83747, "ELG_SGC": 90650,
        }
    ):
        raise ValueError("Previously frozen mock source, ensemble or raw-only contract changed")
    registered = {}
    for row in ran["mock_random_catalogues"]:
        if row["id"] == "0001":
            registered[row["cap"] + "/eBOSS_" + row["tracer"]] = row["compressed_file_sha256"]
    if set(registered) != {cap + "/" + tracer for cap,tracer in ORDER}:
        raise ValueError("Missing four previously SHA-pinned matched EZmock0001 randoms")
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        label = tracer.removeprefix("eBOSS_") + "_" + cap
        if (
            p["exact_expected_header_rows_by_cap_tracer"][key] !=
            header["data_header_row_counts"]["0001"][label]
            or p["filename_rule"] !=
            "{tracer}/dat/EZmock_realistic_{tracer}_{cap}_v7_0001.dat.fits.gz"
        ):
            raise ValueError("Predeclared mock-galaxy filename or row count changed")
    return registered


def acquire(url, dest, cap, timeout):
    approved = urlsplit(BASE)
    requested = urlsplit(url)
    if (
        requested.scheme != "https"
        or requested.netloc != approved.netloc
        or not requested.path.startswith(approved.path)
        or requested.query or requested.fragment
    ):
        raise ValueError("Unexpected source outside exact released EZmock directory")
    partial = Path(str(dest) + ".part")
    if dest.exists():
        raise ValueError("Unpinned mock file already exists; refusing overwrite")
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={
        "User-Agent": "eBOSS-DR16-fixed-EZmock0001-galaxy-raw-SHA/1.0",
        "Accept-Encoding": "identity",
    })
    last = None
    for attempt in range(2):
        if partial.exists():
            partial.unlink()
        try:
            with urlopen(request, timeout=timeout) as response:
                if response.status != 200:
                    raise ValueError("Official full gzip GET was not HTTP 200")
                actual = urlsplit(response.geturl())
                if (
                    actual.scheme != "https"
                    or actual.netloc != requested.netloc
                    or actual.path != requested.path
                    or actual.query or actual.fragment
                ):
                    raise ValueError("Official source redirected to a different file")
                declared_raw = response.headers.get("Content-Length")
                declared = None
                if declared_raw is not None:
                    if not declared_raw.isascii() or not declared_raw.isdecimal():
                        raise ValueError("Invalid HTTP Content-Length")
                    declared = int(declared_raw)
                with partial.open("xb") as writer:
                    sha, size = stream_bytes(response, writer, cap, declared)
                etag = response.headers.get("ETag")
                modified = response.headers.get("Last-Modified")
            partial.replace(dest)
            return {
                "source_url":url,
                "returned_url":response.geturl(),
                "local_path":str(dest.resolve()),
                "first_seen_full_compressed_sha256":sha,
                "compressed_bytes":size,
                "HTTP_status":200,
                "HTTP_content_length_if_available":declared,
                "HTTP_etag_if_available":etag,
                "HTTP_last_modified_if_available":modified,
                "gzip_magic_checked":True,
                "gzip_payload_decompressed":False,
                "FITS_headers_or_mock_rows_read":False,
            }
        except (OSError, URLError, TimeoutError) as exc:
            last = exc
            if attempt == 1:
                raise
            print("RETRY_SAME_OFFICIAL_MOCK_URL",attempt+1,
                  type(exc).__name__,flush=True)
            time.sleep(2)
        finally:
            if partial.exists():
                partial.unlink()
    raise RuntimeError("Official mock source failed") from last


def atomic(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result, indent=2) + chr(10), encoding="utf-8")
    tmp.replace(path)


def run(p, timeout):
    if timeout <= 0:
        raise ValueError("Timeout must be positive")
    random_sha = preflight(p)
    dest = ROOT / p["local_report"]
    protocol_sha = hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    if dest.exists():
        report = load(dest)
        if (
            report.get("protocol_sha256") != protocol_sha
            or report.get("status") not in (PASS, PARTIAL)
            or report.get("mock_realization_id") != 1
            or report.get("observed_odd_data_vector_read") is not False
            or report.get("mock_galaxy_rows_read") is not False
        ):
            raise ValueError("Prior mock-source checkpoint differs; refusing silent repinning")
        report["errors"] = []
        report["status"] = PARTIAL
    else:
        report = {
            "status": PARTIAL,
            "protocol_sha256":protocol_sha,
            "mock_realization_id":1,
            "previously_pinned_companion_mock_random_SHA256":random_sha,
            "samples":{},
            "mock_galaxy_rows_read":False,
            "mock_random_rows_read":False,
            "observed_galaxy_rows_read":False,
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
            "mock_galaxy_pair_counts_computed":False,
            "errors":[],
        }
    total = 0
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        relative = mock_path(tracer, cap, "dat", 1)
        url = BASE + relative
        path = ROOT / p["local_quarantine_dir"] / relative
        prior = report["samples"].get(key)
        if prior is not None:
            sha, size = digest(path, p["max_compressed_bytes_per_file"])
            if (
                prior.get("source_url") != url
                or prior.get("local_path") != str(path.resolve())
                or prior.get("first_seen_full_compressed_sha256") != sha
                or prior.get("compressed_bytes") != size
                or prior.get("FITS_headers_or_mock_rows_read") is not False
            ):
                raise ValueError("Previously SHA-pinned mock source changed: " + key)
            print("REUSED_PINNED_MOCK_GALAXY_SOURCE",key,size,sha,flush=True)
        else:
            record = acquire(url,path,p["max_compressed_bytes_per_file"],timeout)
            record["previously_declared_header_rows_NOT_READ"] = (
                p["exact_expected_header_rows_by_cap_tracer"][key]
            )
            record["matched_mock_random_sha256_FROM_PRIOR_AUDIT"] = random_sha[key]
            report["samples"][key] = record
            atomic(dest,report)
            size=record["compressed_bytes"]
            print("FIRST_SEEN_MOCK_GALAXY_SHA",key,size,
                  record["first_seen_full_compressed_sha256"],flush=True)
        total += size
        if total > p["max_total_compressed_bytes"]:
            raise ValueError("Four official mock sources exceed frozen total byte cap")
    if set(report["samples"]) != {cap + "/" + tracer for cap,tracer in ORDER}:
        raise ValueError("Four fixed mock sources incomplete")
    report["total_four_galaxy_compressed_bytes"] = total
    report["status"] = PASS
    report["errors"] = []
    atomic(dest,report)
    return report


def self_test(p):
    registered=preflight(p)
    assert len(registered)==4
    assert list(registered)==[
        "NGC/eBOSS_LRG","NGC/eBOSS_ELG",
        "SGC/eBOSS_LRG","SGC/eBOSS_ELG"]
    raw=GZIP_PREFIX + b"only synthetic raw byte data"
    buf=io.BytesIO()
    got,n=stream_bytes(io.BytesIO(raw),buf,128,len(raw))
    assert n==len(raw) and got==hashlib.sha256(raw).hexdigest()
    for blob,declared,cap in (
        (raw[:-1],len(raw),128),
        (b"BAD"+raw[3:],len(raw),128),
        (raw+b"extra",len(raw),128),
        (raw,len(raw),len(raw)-1),
    ):
        try:
            stream_bytes(io.BytesIO(blob),io.BytesIO(),cap,declared)
        except ValueError:
            pass
        else:
            raise AssertionError("Damaged/over-limit synthetic official mock source accepted")
    print("EBOSS_EZMOCK0001_GALAXY_SOURCE_SHA_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--timeout",type=float,default=120)
    args=ap.parse_args()
    p=load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    try:
        res=run(p,args.timeout)
    except Exception as exc:
        dest=ROOT/p["local_report"]
        if dest.exists():
            res=load(dest)
            res["status"]=PARTIAL
            res["errors"]=[str(exc)]
        else:
            res={
                "status":PARTIAL,"errors":[str(exc)],
                "observed_odd_data_vector_read":False,
                "mock_galaxy_rows_read":False,
                "new_science_selection_applied":False,
            }
        atomic(dest,res)
        print("EBOSS_EZMOCK0001_GALAXY_SOURCE",res["status"],flush=True)
        print("REPORT",dest,flush=True)
        print("ERRORS",*res["errors"],sep=chr(10),flush=True)
        return 2
    print("EBOSS_EZMOCK0001_GALAXY_SOURCE",res["status"],flush=True)
    print("REPORT",ROOT/p["local_report"],flush=True)
    print("TOTAL_COMPRESSED_BYTES",res["total_four_galaxy_compressed_bytes"],
          flush=True)
    print("ODD_DATA_READ",res["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
