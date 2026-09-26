#!/usr/bin/env python3
"""SHA-pin EXACT raw gzip bytes of four preselected official ELG FITS images.

No gzip decompression, FITS header/image inspection, galaxy/random reads,
catalogue cuts or observed odd-vector access. One predeclared middle
filename per frozen eBOSS ELG chunk, maximum 64 MiB per image.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

from audit_eboss_dr16_elg_maskbit_index import (
    PIN as INDEX_PIN,
    ROOT,
    audit as audit_index,
    load,
)

PROTOCOL = ROOT / "source_data/eboss_dr16_elg_maskbit_four_sample_raw_sha_protocol_2026-09-26.json"
PARENT_MANIFEST = ROOT / "source_data/eboss_dr16_elg_maskbit_index_uploaded_manifest_2026-09-26.json"
OFFLINE_ARCHIVE = ROOT / "source_data/eboss_dr16_elg_maskbit_index_local_report_2026-09-26.json"
CHUNK_ORDER = ("eboss21", "eboss22", "eboss23", "eboss25")
COMPLETE = "OFFICIAL_ELG_FOUR_SAMPLE_RAW_GZIP_SHA_FIRST_SEEN_ONLY"
PARTIAL = "OFFICIAL_ELG_FOUR_SAMPLE_RAW_GZIP_INCOMPLETE_STOP"


def digest_file(path, cap):
    sha = hashlib.sha256()
    n, prefix = 0, b""
    with path.open("rb") as fd:
        while True:
            data = fd.read(1024 * 1024)
            if not data:
                break
            if n + len(data) > cap:
                raise ValueError("Raw gzip exceeds frozen byte cap: " + str(path))
            prefix += data[:max(0, 3 - len(prefix))]
            sha.update(data)
            n += len(data)
    if prefix != bytes((31, 139, 8)):
        raise ValueError("Expected gzip file magic/compression method missing: " + str(path))
    return sha.hexdigest(), n


def copy_sha(reader, writer, expected, cap):
    if expected < 16 or expected > cap:
        raise ValueError("Official Content-Length outside frozen gzip byte cap")
    sha = hashlib.sha256()
    n, prefix = 0, b""
    while True:
        data = reader.read(1024 * 1024)
        if not data:
            break
        if n + len(data) > expected or n + len(data) > cap:
            raise ValueError("Official response delivered more bytes than declared")
        prefix += data[:max(0, 3 - len(prefix))]
        writer.write(data)
        sha.update(data)
        n += len(data)
    if n != expected or prefix != bytes((31, 139, 8)):
        raise ValueError("Official gzip Content-Length or magic mismatch")
    return sha.hexdigest(), n


def preflight(p):
    manifest, index_pin = load(PARENT_MANIFEST), load(INDEX_PIN)
    local_raw = (ROOT / p["parent_offline_report_local"]).read_bytes()
    archive_raw = OFFLINE_ARCHIVE.read_bytes()
    report_sha = hashlib.sha256(local_raw).hexdigest()
    if (
        local_raw != archive_raw
        or report_sha != p["parent_offline_index_report_sha256"]
        or report_sha != manifest["report_sha256"]
        or p["parent_offline_index_report_git_blob_sha1"] != manifest["report_git_blob_sha1"]
        or p["parent_offline_index_report_git_blob_sha1"] !=
           hashlib.sha1(
               ("blob " + str(len(local_raw))).encode("ascii")
               + bytes((0,)) + local_raw
           ).hexdigest()
        or p["parent_uploaded_official_inventory_sha256"] !=
           manifest["original_official_inventory_sha256"]
        or p["parent_uploaded_official_inventory_sha256"] !=
           index_pin["original_uploaded_report_sha256"]
        or p["approved_source_directory"] != index_pin["published_elg_directory"]
        or p["expected_family_order"] != list(index_pin["ordered_chunk_families"])
        or tuple(p["expected_family_order"]) != CHUNK_ORDER
        or p["max_compressed_bytes_per_selected_image"] != 67108864
        or p["max_total_selected_compressed_bytes"] != 268435456
        or p["observed_odd_data_vector_read"] is not False
        or p["gzip_payload_decompressed"] is not False
        or p["fits_image_pixels_read"] is not False
        or p["new_mask_or_selection_applied"] is not False
    ):
        raise ValueError("Exact archived offline index/parent sample protocol changed")
    local_report = json.loads(local_raw)
    if (
        local_report.get("status") != "OFFICIAL_ELG_FOUR_FAMILY_FILENAME_INDEX_SHA_PINNED_ONLY"
        or local_report.get("total_listed_maskbit_fits") != 19381
        or local_report.get("observed_odd_data_vector_read") is not False
        or local_report.get("maskbit_image_bytes_downloaded_or_verified") is not False
        or local_report.get("errors") != []
    ):
        raise ValueError("Previously archived four-family local report did not pass")
    source_index_file = ROOT / index_pin["local_source_report"]
    new_report = audit_index(index_pin, source_index_file)
    if local_report != new_report:
        raise ValueError("Previously archived offline index result changed")
    official_index = load(source_index_file)
    for chunk in CHUNK_ORDER:
        names = official_index["elg_chunks"][chunk]["maskbit_filenames"]
        chosen = names[len(names) // 2]
        if (
            p["selected_filename_by_chunk"][chunk] != chosen
            or chosen != local_report["families"][chunk]["middle"]
            or chosen != index_pin["ordered_chunk_families"][chunk]["middle"]
        ):
            raise ValueError("Changed predeclared central sample filename: " + chunk)
    return report_sha


def new_report(p, report_sha):
    return {
        "status": PARTIAL,
        "source_offline_index_report_sha256": report_sha,
        "official_inventory_sha256": p["parent_uploaded_official_inventory_sha256"],
        "sampling_rule": p["sampling_rule"],
        "selected_filename_by_chunk": p["selected_filename_by_chunk"],
        "samples": {},
        "sample_gzip_bytes_downloaded_or_verified": True,
        "gzip_or_fits_headers_parsed": False,
        "fits_image_pixels_read": False,
        "observed_catalogue_random_mock_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_mask_or_selection_applied": False,
        "physical_elg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "first_seen_gzip_sha256_NOT_OFFICIAL_PUBLISHED_REFERENCE": True,
        "errors": [],
    }


def atomic_report(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(".tmp.json")
    part.write_text(json.dumps(result, indent=2) + chr(10), encoding="utf-8")
    part.replace(path)


def acquire_one(url, dest, cap, timeout):
    if not url.startswith(
        "https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/ELGmasks/"
    ):
        raise ValueError("Unapproved official selected-image source URL")
    req = Request(url, headers={
        "User-Agent": "eBOSS-DR16-ELG-blinded-four-sample-raw-byte-audit/1.0",
        "Accept-Encoding": "identity",
    })
    part = Path(str(dest) + ".part")
    if dest.exists():
        raise ValueError("Unpinned destination already exists; refusing replacement: " + str(dest))
    if part.exists():
        part.unlink()  # Previous partial download is not a certified source.
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urlopen(req, timeout=timeout) as response:
            if response.status != 200 or response.geturl() != url:
                raise ValueError("Official response status or URL changed")
            declared = response.headers.get("Content-Length")
            if declared is None or not declared.isascii() or not declared.isdecimal():
                raise ValueError("Official gzip response lacks exact Content-Length")
            declared = int(declared)
            with part.open("xb") as fd:
                sha, size = copy_sha(response, fd, declared, cap)
            etag = response.headers.get("ETag")
            modified = response.headers.get("Last-Modified")
        part.replace(dest)
    except Exception:
        if part.exists():
            part.unlink()
        raise
    return {
        "source_url": url,
        "local_path": str(dest.resolve()),
        "compressed_file_bytes": size,
        "full_compressed_file_sha256_first_seen": sha,
        "http_status": 200,
        "http_content_length": declared,
        "http_etag_if_available": etag,
        "http_last_modified_if_available": modified,
        "raw_gzip_magic_verified": True,
        "gzip_decompressed": False,
        "fits_header_or_pixels_read": False,
    }


def run(p, timeout):
    if timeout <= 0:
        raise ValueError("Timeout must be positive")
    report_sha = preflight(p)
    out = ROOT / p["local_report"]
    if out.exists():
        result = load(out)
        if (
            result.get("status") not in (COMPLETE, PARTIAL)
            or result.get("source_offline_index_report_sha256") != report_sha
            or result.get("selected_filename_by_chunk") != p["selected_filename_by_chunk"]
            or result.get("observed_odd_data_vector_read") is not False
            or result.get("gzip_or_fits_headers_parsed") is not False
        ):
            raise ValueError("Existing sample report differs; no silent repinning")
        result["errors"] = []
        result["status"] = PARTIAL
    else:
        result = new_report(p, report_sha)
    cap = p["max_compressed_bytes_per_selected_image"]
    folder = ROOT / p["retained_raw_files_directory"]
    cumulative = 0
    for chunk in CHUNK_ORDER:
        filename = p["selected_filename_by_chunk"][chunk]
        url = p["approved_source_directory"] + filename
        dest = folder / chunk / filename
        prior = result["samples"].get(chunk)
        if prior is not None:
            sha, size = digest_file(dest, cap)
            if (
                prior.get("source_url") != url
                or prior.get("local_path") != str(dest.resolve())
                or prior.get("full_compressed_file_sha256_first_seen") != sha
                or prior.get("compressed_file_bytes") != size
                or prior.get("gzip_decompressed") is not False
                or prior.get("fits_header_or_pixels_read") is not False
            ):
                raise ValueError("Previously downloaded selected gzip bytes changed: " + chunk)
            print("REUSED_PINNED_CHUNK", chunk, size, sha, flush=True)
        else:
            acquired = acquire_one(url, dest, cap, timeout)
            result["samples"][chunk] = acquired
            atomic_report(out, result)
            size = acquired["compressed_file_bytes"]
            print("SHA_PINNED_CHUNK", chunk, size,
                  acquired["full_compressed_file_sha256_first_seen"], flush=True)
        cumulative += size
        if cumulative > p["max_total_selected_compressed_bytes"]:
            raise ValueError("Four selected gzip files exceed combined byte cap")
    if set(result["samples"]) != set(CHUNK_ORDER):
        raise ValueError("Incomplete four-family sample response")
    result["status"] = COMPLETE
    result["total_four_sample_compressed_bytes"] = cumulative
    result["errors"] = []
    atomic_report(out, result)
    return result


def self_test(p):
    if (
        tuple(p["selected_filename_by_chunk"]) != CHUNK_ORDER
        or p["selected_filename_by_chunk"] != {
            "eboss21": "mask-eboss21-3386m005.fits.gz",
            "eboss22": "mask-eboss22-0226m005.fits.gz",
            "eboss23": "mask-eboss23-1415p232.fits.gz",
            "eboss25": "mask-eboss25-1538p312.fits.gz",
        }
        or p["max_compressed_bytes_per_selected_image"] != 67108864
        or p["observed_odd_data_vector_read"] is not False
    ):
        raise AssertionError("Four frozen official sample filenames or bounds changed")
    fake = bytes((31, 139, 8)) + b"synthetic only 123"
    copied = io.BytesIO()
    sha, size = copy_sha(io.BytesIO(fake), copied, len(fake), len(fake))
    assert copied.getvalue() == fake
    assert size == len(fake)
    assert sha == hashlib.sha256(fake).hexdigest()
    for damaged, declared, cap in (
        (fake[:-1], len(fake), len(fake)),
        (b"BAD" + fake[3:], len(fake), len(fake)),
        (fake + b"X", len(fake), len(fake) + 1),
        (fake, len(fake), len(fake) - 1),
    ):
        try:
            copy_sha(io.BytesIO(damaged), io.BytesIO(), declared, cap)
        except ValueError:
            pass
        else:
            raise AssertionError("Corrupted synthetic official gzip response accepted")
    print("EBOSS_ELG_FOUR_GZIP_SAMPLE_RAW_SHA_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--timeout", type=float, default=120)
    args = ap.parse_args()
    p = load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    try:
        result = run(p, args.timeout)
    except Exception as exc:
        out = ROOT / p["local_report"]
        if out.exists():
            result = load(out)
            result["status"] = PARTIAL
            result["errors"] = [str(exc)]
        else:
            result = {
                "status": PARTIAL,
                "errors": [str(exc)],
                "observed_odd_data_vector_read": False,
                "new_mask_or_selection_applied": False,
                "physical_elg_mask_certified": False,
            }
        atomic_report(out, result)
        print("EBOSS_ELG_FOUR_SAMPLE_RAW_SHA", result["status"], flush=True)
        print("REPORT", out, flush=True)
        print("ERRORS", *result["errors"], sep=chr(10), flush=True)
        return 2
    print("EBOSS_ELG_FOUR_SAMPLE_RAW_SHA", result["status"], flush=True)
    print("REPORT", ROOT / p["local_report"], flush=True)
    print("TOTAL_COMPRESSED_BYTES", result["total_four_sample_compressed_bytes"], flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
