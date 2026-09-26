#!/usr/bin/env python3
"""Metadata-only FITS/WCS audit of FOUR SHA-pinned official ELG gzip images.

Preflight requires exact byte-identical previously uploaded raw-gzip
report and all four already frozen compressed file SHA256 identities.
Only bounded decompressed FITS header blocks are returned to this code.
gzip may buffer/decompress bytes internally; no image pixel values
are requested, decoded, summarized or used for any mask selection.
No network, catalogue/random access or observed odd-vector access.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

from audit_eboss_dr16_elg_four_sample_raw_sha import (
    CHUNK_ORDER, PROTOCOL as PARENT, digest_file,
)
from audit_eboss_dr16_elg_maskbit_index import ROOT, load
from inspect_eboss_dr16_fits_headers import BLOCK, parse_header

PROTOCOL = ROOT / "source_data/eboss_dr16_elg_four_sample_fits_header_wcs_protocol_2026-09-26.json"
ARCHIVE = ROOT / "source_data/eboss_dr16_elg_four_sample_raw_sha_report_2026-09-26.json"
MANIFEST = ROOT / "source_data/eboss_dr16_elg_four_sample_raw_sha_uploaded_manifest_2026-09-26.json"
COMPLETE = "OFFICIAL_ELG_FOUR_SAMPLE_FITS_HEADER_WCS_METADATA_ONLY"
INCOMPLETE = "OFFICIAL_ELG_FOUR_SAMPLE_FITS_HEADER_WCS_INCOMPLETE_STOP"


def verified_parent(p):
    prior = load(PARENT)
    manifest = load(MANIFEST)
    raw = ARCHIVE.read_bytes()
    local = (ROOT / p["parent_local_gzip_report"]).read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    blob = hashlib.sha1(
        ("blob " + str(len(raw))).encode("ascii") + bytes((0,)) + raw
    ).hexdigest()
    if (
        raw != local
        or sha != p["parent_local_gzip_report_sha256"]
        or sha != manifest["uploaded_report_sha256"]
        or blob != p["parent_local_gzip_report_git_blob_sha1"]
        or blob != manifest["uploaded_report_git_blob_sha1"]
        or len(raw) != manifest["uploaded_report_bytes"]
        or p["parent_index_report_sha256"] !=
           prior["parent_offline_index_report_sha256"]
        or p["expected_family_order"] != list(CHUNK_ORDER)
        or p["four_predeclared_compressed_source_fingerprints"] !=
           manifest["sample_sha256_first_seen_not_official_reference"]
        or prior["selected_filename_by_chunk"] != {
            chunk: value["filename"] for chunk, value in
            p["four_predeclared_compressed_source_fingerprints"].items()
        }
        or p["total_sample_compressed_bytes"] != 676165
        or p["observed_odd_data_vector_read"] is not False
        or p["image_pixel_values_decoded"] is not False
        or p["new_mask_or_selection_applied"] is not False
        or p["max_fits_header_bytes_per_image"] != 57600
    ):
        raise ValueError("Previously archived gzip source report or frozen header protocol differs")
    report = json.loads(raw)
    if (
        report.get("status") != "OFFICIAL_ELG_FOUR_SAMPLE_RAW_GZIP_SHA_FIRST_SEEN_ONLY"
        or report.get("errors") != []
        or report.get("source_offline_index_report_sha256") !=
           p["parent_index_report_sha256"]
        or report.get("total_four_sample_compressed_bytes") !=
           p["total_sample_compressed_bytes"]
        or report.get("gzip_or_fits_headers_parsed") is not False
        or report.get("fits_image_pixels_read") is not False
        or report.get("observed_odd_data_vector_read") is not False
        or report.get("new_mask_or_selection_applied") is not False
        or list(report.get("samples", {})) != list(CHUNK_ORDER)
    ):
        raise ValueError("First-seen local four-gzip report is incomplete or changed")
    cap = prior["max_compressed_bytes_per_selected_image"]
    paths = {}
    for chunk in CHUNK_ORDER:
        item = report["samples"][chunk]
        frozen = p["four_predeclared_compressed_source_fingerprints"][chunk]
        dest = (ROOT / prior["retained_raw_files_directory"] /
                chunk / frozen["filename"])
        if (
            item.get("source_url") != prior["approved_source_directory"] + frozen["filename"]
            or item.get("local_path") != str(dest.resolve())
            or item.get("full_compressed_file_sha256_first_seen") != frozen["sha256"]
            or item.get("compressed_file_bytes") != frozen["bytes"]
            or item.get("http_content_length") != frozen["bytes"]
            or item.get("http_status") != 200
            or item.get("raw_gzip_magic_verified") is not True
            or item.get("gzip_decompressed") is not False
            or item.get("fits_header_or_pixels_read") is not False
        ):
            raise ValueError("Preselected official gzip sample record mismatch: " + chunk)
        # Full compressed-file SHA check happens BEFORE any gzip object opens.
        checked_sha, checked_size = digest_file(dest, cap)
        if checked_sha != frozen["sha256"] or checked_size != frozen["bytes"]:
            raise ValueError("Official selected FITS.gz raw bytes changed: " + chunk)
        paths[chunk] = dest
        print("REVERIFIED_COMPRESSED_SHA", chunk, checked_size,
              checked_sha, flush=True)
    if sum(report["samples"][chunk]["compressed_file_bytes"]
           for chunk in CHUNK_ORDER) != p["total_sample_compressed_bytes"]:
        raise ValueError("Four official compressed sample byte totals changed")
    return paths


def read_single_header(gz, max_bytes):
    if max_bytes < BLOCK:
        raise ValueError("FITS header byte cap reached before header")
    raw = bytearray()
    for _ in range(max_bytes // BLOCK):
        block = gz.read(BLOCK)
        if len(block) != BLOCK:
            raise ValueError("Truncated FITS header block in selected gzip file")
        raw.extend(block)
        if any(block[i:i+8] == b"END     " for i in range(0, BLOCK, 80)):
            parsed, end_offset = parse_header(bytes(raw), 0)
            if end_offset != len(raw):
                raise ValueError("FITS header parser found inconsistent END block")
            return parsed, bytes(raw)
    raise ValueError("FITS END card absent before frozen bounded header cap")


def required_number(header, key):
    if key not in header or isinstance(header[key], bool):
        raise ValueError("Required FITS WCS numeric keyword missing/boolean: " + key)
    try:
        v = float(str(header[key]).replace("D", "E").replace("d", "E"))
    except (ValueError, TypeError) as exc:
        raise ValueError("Required FITS WCS keyword is not numeric: " + key) from exc
    if not math.isfinite(v):
        raise ValueError("Nonfinite required FITS WCS keyword: " + key)
    return v


def image_metadata(image, p):
    if (
        image.get("NAXIS") != p["expected_NAXIS"]
        or image.get("BITPIX") not in p["allowed_BITPIX_values"]
        or image.get("CTYPE1") != p["expected_CTYPE1"]
        or image.get("CTYPE2") != p["expected_CTYPE2"]
    ):
        raise ValueError("Published eBOSS BRICKMASK IMAGE/NAXIS/BITPIX/TAN WCS contract failed")
    n1, n2 = image.get("NAXIS1"), image.get("NAXIS2")
    if (
        not isinstance(n1, int) or isinstance(n1, bool)
        or not isinstance(n2, int) or isinstance(n2, bool)
        or n1 <= 0 or n2 <= 0
        or n1 * n2 > p["max_image_pixel_count_declared_header_only"]
    ):
        raise ValueError("Invalid or unexpectedly large declared FITS image dimensions")
    keys = p["exact_required_image_header_keys"]
    if keys != [
        "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
        "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2",
        "CRPIX1", "CRPIX2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
    ]:
        raise ValueError("Pinned FITS WCS metadata whitelist differs")
    if not all(key in image for key in keys):
        raise ValueError("One or more required upstream FITS/WCS keywords missing")
    numeric = {key: required_number(image, key) for key in
               ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
                "CD1_1", "CD1_2", "CD2_1", "CD2_2")}
    determinant = (numeric["CD1_1"] * numeric["CD2_2"]
                   - numeric["CD1_2"] * numeric["CD2_1"])
    if not math.isfinite(determinant) or determinant == 0:
        raise ValueError("FITS image WCS CD matrix is nonfinite/singular")
    if not (0 <= numeric["CRVAL1"] <= 360
            and -90 <= numeric["CRVAL2"] <= 90):
        raise ValueError("FITS image WCS reference RA/Dec outside standard domain")
    return {
        "BITPIX": image["BITPIX"],
        "NAXIS": image["NAXIS"],
        "NAXIS1": n1,
        "NAXIS2": n2,
        "CTYPE1": image["CTYPE1"],
        "CTYPE2": image["CTYPE2"],
        "CRVAL1": numeric["CRVAL1"],
        "CRVAL2": numeric["CRVAL2"],
        "CRPIX1": numeric["CRPIX1"],
        "CRPIX2": numeric["CRPIX2"],
        "CD1_1": numeric["CD1_1"],
        "CD1_2": numeric["CD1_2"],
        "CD2_1": numeric["CD2_1"],
        "CD2_2": numeric["CD2_2"],
        "CD_matrix_determinant_nonzero": True,
        "declared_image_pixel_count_NOT_READ": n1 * n2,
    }


def inspect_selected(path, p):
    cap = p["max_fits_header_bytes_per_image"]
    with path.open("rb") as compressed:
        with gzip.GzipFile(fileobj=compressed, mode="rb") as gz:
            primary, primary_bytes = read_single_header(gz, cap)
            if primary.get("SIMPLE") is not True:
                raise ValueError("First selected FITS HDU is not SIMPLE=T primary")
            if primary.get("NAXIS") == 2:
                image, image_bytes = primary, primary_bytes
                kind = "PRIMARY_IMAGE"
                total_raw = primary_bytes
            elif primary.get("NAXIS") == 0:
                extension, image_bytes = read_single_header(
                    gz, cap - len(primary_bytes)
                )
                if extension.get("XTENSION") != "IMAGE":
                    raise ValueError("Empty FITS primary not followed by IMAGE HDU")
                image = extension
                kind = "EMPTY_PRIMARY_PLUS_FIRST_IMAGE_EXTENSION"
                total_raw = primary_bytes + image_bytes
            else:
                raise ValueError("First FITS HDU is neither image nor empty primary")
    result = image_metadata(image, p)
    return {
        "first_data_image_hdu": kind,
        "total_decompressed_fits_header_bytes_READ_ONLY": len(total_raw),
        "complete_header_blocks_sha256": hashlib.sha256(total_raw).hexdigest(),
        "primary_header_blocks_sha256": hashlib.sha256(primary_bytes).hexdigest(),
        "image_header_blocks_sha256": hashlib.sha256(image_bytes).hexdigest(),
        "image_header_metadata": result,
        "gzip_may_have_buffered_data_internally": True,
        "image_pixel_values_decoded": False,
    }


def self_test(p):
    def block(cards):
        raw = b"".join(
            (f"{key:<8}= {value}" if value is not None else key)
            .ljust(80).encode("ascii")
            for key, value in cards + [("END", None)]
        )
        return raw.ljust(((len(raw) + BLOCK - 1) // BLOCK) * BLOCK, b" ")
    image = [
        ("BITPIX", "8"), ("NAXIS", "2"), ("NAXIS1", "3"),
        ("NAXIS2", "2"), ("CTYPE1", "'RA---TAN'"),
        ("CTYPE2", "'DEC--TAN'"), ("CRVAL1", "338.6"),
        ("CRVAL2", "-0.5"), ("CRPIX1", "1.0"), ("CRPIX2", "1.0"),
        ("CD1_1", "-0.0001"), ("CD1_2", "0.0"),
        ("CD2_1", "0.0"), ("CD2_2", "0.0001"),
    ]
    primary = block([("SIMPLE", "T")] + image)
    synthetic = gzip.compress(primary + b"synthetic-pixel-data-never-requested")
    tmp = io.BytesIO(synthetic)
    with gzip.GzipFile(fileobj=tmp, mode="rb") as gz:
        header, raw = read_single_header(gz, p["max_fits_header_bytes_per_image"])
        meta = image_metadata(header, p)
    assert raw == primary
    assert meta["declared_image_pixel_count_NOT_READ"] == 6
    assert meta["CD_matrix_determinant_nonzero"] is True
    empty = block([("SIMPLE", "T"), ("BITPIX", "8"), ("NAXIS", "0")])
    ext = block([("XTENSION", "'IMAGE'")] + image)
    with gzip.GzipFile(fileobj=io.BytesIO(gzip.compress(empty + ext + b"pixels")),
                       mode="rb") as gz:
        first, a = read_single_header(gz, p["max_fits_header_bytes_per_image"])
        second, b = read_single_header(gz, p["max_fits_header_bytes_per_image"] - len(a))
    assert first["NAXIS"] == 0 and second["XTENSION"] == "IMAGE"
    assert a == empty and b == ext
    image_metadata(second, p)
    invalid = dict(header, CTYPE1="DEC--TAN")
    try:
        image_metadata(invalid, p)
    except ValueError:
        pass
    else:
        raise AssertionError("Invalid synthetic FITS WCS CTYPE was accepted")
    singular = dict(header, CD1_1="0", CD1_2="0", CD2_1="0", CD2_2="0")
    try:
        image_metadata(singular, p)
    except ValueError:
        pass
    else:
        raise AssertionError("Singular synthetic FITS WCS was accepted")
    print("EBOSS_ELG_FOUR_SAMPLE_FITS_HEADER_WCS_SYNTHETIC_SELF_TEST_OK", flush=True)


def atomic_report(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(".tmp.json")
    temp.write_text(json.dumps(value, indent=2) + chr(10), encoding="utf-8")
    temp.replace(path)


def audit(p):
    sources = verified_parent(p)
    result = {
        "status": COMPLETE,
        "first_seen_compressed_SHA256_reverified_before_any_FITS_header": True,
        "exact_uploaded_parent_report_sha256": p["parent_local_gzip_report_sha256"],
        "predeclared_sample_metadata": {},
        "images": {},
        "gzip_may_buffer_decompress_internally": True,
        "no_image_pixel_values_returned_or_decoded": True,
        "observed_catalogue_random_mock_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_mask_or_selection_applied": False,
        "full_elg_image_cohort_certified": False,
        "physical_elg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "errors": [],
    }
    for chunk in CHUNK_ORDER:
        frozen = p["four_predeclared_compressed_source_fingerprints"][chunk]
        result["predeclared_sample_metadata"][chunk] = frozen
        result["images"][chunk] = inspect_selected(sources[chunk], p)
        print("HEADER_WCS_CHUNK", chunk,
              json.dumps(result["images"][chunk], sort_keys=True), flush=True)
    return result


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
            "status": INCOMPLETE,
            "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "new_mask_or_selection_applied": False,
            "physical_elg_mask_certified": False,
        }
    atomic_report(out, result)
    print("EBOSS_ELG_FOUR_SAMPLE_HEADER_WCS", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep=chr(10), flush=True)
        return 2
    print("HEADER_COUNT", len(result["images"]), flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
