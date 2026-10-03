#!/usr/bin/env python3
"""SHA-first four-image eBOSS ELG pixel BYTE histogram audit.

Input: exactly four already preregistered and SHA-pinned FITS.gz maskbit
files. This source-informed QA reads 8-bit mask IMAGE payloads only,
producing 256-bin counts and bit0/bit2 totals (not spatial selections).
It never reads galaxy/random/mock rows, applies new masks or accesses
the observed odd-sector vector. Never generalize 4 samples to 19381.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np

from audit_eboss_dr16_elg_four_sample_fits_header_wcs import (
    CHUNK_ORDER, PROTOCOL as HEADER_PROTOCOL, image_metadata,
    read_single_header, verified_parent,
)
from audit_eboss_dr16_elg_maskbit_index import ROOT, load

PROTOCOL = ROOT / "source_data/eboss_dr16_elg_four_sample_pixel_bitcode_protocol_2026-09-26.json"
HEADER_ARCHIVE = ROOT / "source_data/eboss_dr16_elg_four_sample_fits_header_wcs_report_2026-09-26.json"
HEADER_MANIFEST = ROOT / "source_data/eboss_dr16_elg_four_sample_fits_header_wcs_uploaded_manifest_2026-09-26.json"
COMPLETE = "OFFICIAL_ELG_FOUR_SAMPLE_PIXEL_AGGREGATES_SOURCE_QA_ONLY"
INCOMPLETE = "OFFICIAL_ELG_FOUR_SAMPLE_PIXEL_AGGREGATES_INCOMPLETE_STOP"
CHUNK_BYTES = 1024 * 1024


def prior_gate(p):
    header_protocol, manifest = load(HEADER_PROTOCOL), load(HEADER_MANIFEST)
    archived_bytes = HEADER_ARCHIVE.read_bytes()
    local_bytes = (ROOT / p["parent_exact_header_report_local"]).read_bytes()
    sha = hashlib.sha256(archived_bytes).hexdigest()
    blob = hashlib.sha1(
        ("blob " + str(len(archived_bytes))).encode("ascii")
        + bytes((0,)) + archived_bytes
    ).hexdigest()
    if (
        archived_bytes != local_bytes
        or len(archived_bytes) != manifest["uploaded_report_bytes"]
        or sha != p["parent_exact_header_report_sha256"]
        or sha != manifest["uploaded_report_sha256"]
        or blob != p["parent_exact_header_report_git_blob_sha1"]
        or blob != manifest["uploaded_report_git_blob_sha1"]
        or p["frozen_sample_order"] != list(CHUNK_ORDER)
        or p["frozen_image_header_blocks_sha256"] !=
           manifest["complete_header_SHA256_by_chunk"]
        or p["expected_pixels_per_image"] != 12960000
        or p["max_uncompressed_fits_file_bytes_per_image"] != 20000000
        or p["histogram_bins"] != 256
        or p["bit_values_whitelist"] != list(range(8))
        or p["frozen_image_metadata"] != {
            "BITPIX": 8, "NAXIS": 2, "NAXIS1": 3600, "NAXIS2": 3600,
            "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
            "image_hdu": "PRIMARY_IMAGE", "complete_header_bytes": 5760,
        }
        or p["data_scope_guards"]["observed_odd_data_vector_read"] is not False
        or p["data_scope_guards"]["new_science_selection_applied"] is not False
    ):
        raise ValueError("Exact uploaded FITS header report or pixel protocol changed")
    code = p["source_code_contract"]
    if (
        code["repo"] != "cheng-zhao/brickmask"
        or code["commit"] != "b9eb684a579b56ec3dbdb46549224be7e3fa2830"
        or code["compile_flag"] != "-DEBOSS"
        or code["define_header_git_blob_sha1"] !=
           "a6b85bfbc7745388a076fba822a0dc27d52d9809"
        or code["bitcode_git_blob_sha1"] !=
           "0ab07730c491324fac3222ebac1ae1c763efe980"
        or code["fits_reader_git_blob_sha1"] !=
           "f5262cd33e9b011131915e5a6957a3fc86d37d5a"
    ):
        raise ValueError("Pinned upstream EBOSS bit0/bit2 source contract changed")
    report = json.loads(archived_bytes)
    if (
        report["status"] !=
           "OFFICIAL_ELG_FOUR_SAMPLE_FITS_HEADER_WCS_METADATA_ONLY"
        or report["first_seen_compressed_SHA256_reverified_before_any_FITS_header"]
           is not True
        or report["exact_uploaded_parent_report_sha256"] !=
           header_protocol["parent_local_gzip_report_sha256"]
        or report["no_image_pixel_values_returned_or_decoded"] is not True
        or report["observed_odd_data_vector_read"] is not False
        or report["errors"] != []
        or list(report["images"]) != list(CHUNK_ORDER)
        or list(report["predeclared_sample_metadata"]) != list(CHUNK_ORDER)
    ):
        raise ValueError("Previously uploaded actual four-image FITS header check changed")
    for chunk in CHUNK_ORDER:
        item = report["images"][chunk]
        m = item["image_header_metadata"]
        frozen = p["frozen_image_metadata"]
        if (
            item["first_data_image_hdu"] != frozen["image_hdu"]
            or item["total_decompressed_fits_header_bytes_READ_ONLY"] !=
               frozen["complete_header_bytes"]
            or item["complete_header_blocks_sha256"] !=
               p["frozen_image_header_blocks_sha256"][chunk]
            or item["image_pixel_values_decoded"] is not False
            or m["declared_image_pixel_count_NOT_READ"] != p["expected_pixels_per_image"]
            or any(m[k] != frozen[k] for k in (
                "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "CTYPE1", "CTYPE2"
            ))
            or report["predeclared_sample_metadata"][chunk] !=
               header_protocol["four_predeclared_compressed_source_fingerprints"][chunk]
        ):
            raise ValueError("Already verified image header or compressed SHA changed: " + chunk)
    # This rechecks all four complete compressed-source SHA256 BEFORE gzip opens.
    paths = verified_parent(header_protocol)
    return report, paths


def stream_pixel_histogram(gz, nbytes, max_uncompressed, header_bytes):
    if (
        nbytes < 1 or nbytes != 12960000 and nbytes != 6
        or len(header_bytes) + nbytes > max_uncompressed
    ):
        raise ValueError("Pixel byte count exceeds frozen image size/bound")
    bins = np.zeros(256, dtype=np.int64)
    payload_sha = hashlib.sha256()
    fits_sha = hashlib.sha256()
    fits_sha.update(header_bytes)
    remaining = nbytes
    while remaining:
        block = gz.read(min(CHUNK_BYTES, remaining))
        if not block:
            raise ValueError("Image payload ended before declared FITS NAXIS1*NAXIS2")
        if len(block) > remaining:
            raise ValueError("Read beyond frozen primary image payload")
        arr = np.frombuffer(block, dtype=np.uint8)
        bins += np.bincount(arr, minlength=256)
        payload_sha.update(block)
        fits_sha.update(block)
        remaining -= len(block)
    # For the selected 3600x3600 8-bit images, FITS image data is
    # already a whole number of 2880-byte blocks; no padding is expected.
    if nbytes % 2880 == 0:
        if gz.read(1) != b"":
            raise ValueError("Unexpected extra decompressed bytes beyond frozen image payload")
    else:
        padding = (-nbytes) % 2880
        if len(header_bytes) + nbytes + padding > max_uncompressed:
            raise ValueError("FITS data padding exceeds frozen decompression cap")
        tail = gz.read(padding)
        if tail != bytes(padding) or gz.read(1) != b"":
            raise ValueError("FITS image padding/trailing bytes inconsistent")
        fits_sha.update(tail)
    if int(np.sum(bins)) != nbytes:
        raise ValueError("Pixel histogram does not account for every declared image byte")
    masks = np.arange(256, dtype=np.uint16)
    bit0 = (masks & 1) != 0
    bit2 = (masks & 4) != 0
    return {
        "total_maskbit_image_pixels": nbytes,
        "raw_image_payload_sha256": payload_sha.hexdigest(),
        "decompressed_fits_header_plus_payload_sha256": fits_sha.hexdigest(),
        "full_byte_value_histogram_0_to_255": [int(n) for n in bins],
        "bit0_set_pixels": int(np.sum(bins[bit0])),
        "bit2_set_pixels": int(np.sum(bins[bit2])),
        "bit0_and_bit2_set_pixels": int(np.sum(bins[bit0 & bit2])),
        "bit0_set_bit2_unset_pixels": int(np.sum(bins[bit0 & ~bit2])),
        "bit2_set_bit0_unset_pixels": int(np.sum(bins[~bit0 & bit2])),
        "gzip_EOF_CRC_and_ISIZE_validated": True,
        "no_pixel_positions_or_sky_coordinates_output": True,
    }


def inspect(path, p, expected_header_sha, expected_metadata):
    with path.open("rb") as source:
        with gzip.GzipFile(fileobj=source, mode="rb") as gz:
            image, header_bytes = read_single_header(
                gz, p["frozen_image_metadata"]["complete_header_bytes"]
            )
            if (
                image.get("SIMPLE") is not True
                or hashlib.sha256(header_bytes).hexdigest() != expected_header_sha
            ):
                raise ValueError("Pre-pixel FITS header SHA/SIMPLE changed")
            canonical_meta = image_metadata(image, load(HEADER_PROTOCOL))
            if canonical_meta != expected_metadata:
                raise ValueError("FITS WCS metadata changed from archived prior report")
            result = stream_pixel_histogram(
                gz, p["expected_pixels_per_image"],
                p["max_uncompressed_fits_file_bytes_per_image"], header_bytes
            )
    result["reverified_pre_pixel_header_sha256"] = expected_header_sha
    result["frozen_image_bitpix"] = image["BITPIX"]
    result["observed_odd_data_vector_read"] = False
    result["new_science_selection_applied"] = False
    return result


def self_test(p):
    if (
        p["expected_pixels_per_image"] != 12960000
        or p["bit_values_whitelist"] != list(range(8))
        or p["source_code_contract"]["define_header_git_blob_sha1"] !=
           "a6b85bfbc7745388a076fba822a0dc27d52d9809"
        or list(p["frozen_image_header_blocks_sha256"]) != list(CHUNK_ORDER)
    ):
        raise AssertionError("Frozen pixel aggregate contract changed")
    header = b"H" * 2880
    pixels = bytes((0, 1, 4, 5, 5, 255))
    raw = header + pixels + bytes(2880 - len(pixels))
    zipped = gzip.compress(raw)
    with gzip.GzipFile(fileobj=io.BytesIO(zipped), mode="rb") as gz:
        assert gz.read(len(header)) == header
        outcome = stream_pixel_histogram(gz, 6, len(raw), header)
    assert outcome["total_maskbit_image_pixels"] == 6
    assert outcome["bit0_set_pixels"] == 4
    assert outcome["bit2_set_pixels"] == 4
    assert outcome["bit0_and_bit2_set_pixels"] == 3
    assert outcome["full_byte_value_histogram_0_to_255"][5] == 2
    assert outcome["gzip_EOF_CRC_and_ISIZE_validated"] is True
    corrupted = bytearray(zipped)
    corrupted[-8] ^= 1
    try:
        with gzip.GzipFile(fileobj=io.BytesIO(corrupted), mode="rb") as gz:
            assert gz.read(len(header)) == header
            stream_pixel_histogram(gz, 6, len(raw), header)
    except (gzip.BadGzipFile, EOFError, OSError):
        pass
    else:
        raise AssertionError("Synthetic corrupted gzip CRC was accepted")
    print("EBOSS_ELG_FOUR_SAMPLE_PIXEL_BITCODE_SYNTHETIC_SELF_TEST_OK", flush=True)


def atomic_report(path, result):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp.json")
    temp.write_text(json.dumps(result, indent=2) + chr(10), encoding="utf-8")
    temp.replace(path)


def audit(p):
    header_report, paths = prior_gate(p)
    report = {
        "status": COMPLETE,
        "source_informed_aggregate_only_not_physical_mask": True,
        "parent_exact_header_report_sha256": p["parent_exact_header_report_sha256"],
        "frozen_sample_order": list(CHUNK_ORDER),
        "sample_histograms": {},
        "code_contract": p["source_code_contract"],
        "histogram_definitions_frozen_before_pixel_read": True,
        "only_four_preselected_mask_image_pixels_decoded": True,
        "no_pixel_positions_or_sky_coordinate_mapping": True,
        "observed_catalogue_or_random_or_mock_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_science_selection_applied": False,
        "physical_elg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "errors": [],
    }
    for chunk in CHUNK_ORDER:
        item = header_report["images"][chunk]
        result = inspect(
            paths[chunk], p,
            p["frozen_image_header_blocks_sha256"][chunk],
            item["image_header_metadata"],
        )
        report["sample_histograms"][chunk] = result
        print("ELG_SAMPLE_PIXEL_AGGREGATE", chunk,
              "PIXELS", result["total_maskbit_image_pixels"],
              "BIT0", result["bit0_set_pixels"],
              "BIT2", result["bit2_set_pixels"],
              "BOTH", result["bit0_and_bit2_set_pixels"],
              "PAYLOAD_SHA256", result["raw_image_payload_sha256"],
              flush=True)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    p = load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    out = ROOT / p["report_output"]
    try:
        result = audit(p)
    except Exception as exc:
        result = {
            "status": INCOMPLETE, "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "new_science_selection_applied": False,
            "physical_elg_mask_certified": False,
        }
    atomic_report(out, result)
    print("EBOSS_ELG_FOUR_SAMPLE_PIXEL_AGGREGATE", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep=chr(10), flush=True)
        return 2
    print("IMAGES", len(result["sample_histograms"]), flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
