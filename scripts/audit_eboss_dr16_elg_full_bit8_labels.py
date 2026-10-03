#!/usr/bin/env python3
"""Blinded comparison of official ELG full-catalogue bit-8 labels with frozen pixels.

The source file's exact 189-MB SHA256 was obtained in an earlier *byte-only*
run and frozen BEFORE this script was written. Check full SHA and header SHA
again, then decode ONLY three requested observed FITS columns: RA, DEC,
mskbit. No IDs, redshifts, weights, randoms, galaxy pair counts or odd vector.
Output aggregate confusion counts only; do not change any science selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from audit_eboss_dr16_elg_bit8_ra_orientation import (
    NSIDE, candidate_pixels, frozen_inputs,
)
from inspect_eboss_dr16_fits_headers import BLOCK, hdu_data_length, parse_header

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_full_bit8_label_protocol_2026-09-25.json"
PARENT = ROOT / "source_data/eboss_dr16_elg_full_bytes_provenance_protocol_2026-09-25.json"
HEADER_REPORT = ROOT / "eboss_workspace/official_mask_inventory/elg_full_header_only.json"
RAW_REPORT = ROOT / "eboss_workspace/official_mask_inventory/elg_full_bytes_provenance.json"
WIDTHS = {"A": 1, "L": 1, "B": 1, "I": 2, "J": 4, "K": 8,
          "E": 4, "D": 8, "C": 8, "M": 16, "P": 8, "Q": 16}
FORM = re.compile(r"^([1-9][0-9]*)?([A-Z])$")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def verify_contract(p, parent, raw, hdr, path):
    if (p["source_url"] != parent["official_source_url"]
            or p["source_file"] != parent["local_output"]
            or p["expected_source_bytes"] != parent["expected_file_bytes"]
            or p["expected_header_sha256"] != parent["expected_header_sha256"]
            or p["expected_header_bytes"] != parent["expected_header_bytes"]
            or p["first_observed_full_file_sha256_from_local_byte_only_run"]
               != "8806699e14422904ef91efb6f1632184171fb740c2e074456c0534b7d2dc28b3"
            or p["bit8_integer_mask"] != 1 << 8
            or p["allowed_observed_row_columns"] != ["RA", "DEC", "mskbit"]
            or p["candidate_coordinate_conventions"] != [
                "literal_upstream", "units_only_reflected_RA",
                "native_RA_radians", "native_RA_degrees"]
            or p["observed_odd_data_vector_read"] is not False
            or p["mask_selection_changed"] is not False
            or raw.get("status") != "OFFICIAL_ELG_FULL_BYTES_SHA_PINNED_ONLY"
            or raw.get("full_file_sha256")
               != p["first_observed_full_file_sha256_from_local_byte_only_run"]
            or raw.get("file_bytes") != p["expected_source_bytes"]
            or raw.get("expected_header_sha256") != p["expected_header_sha256"]
            or raw.get("observed_galaxy_data_rows_parsed") is not False
            or raw.get("observed_odd_data_vector_read") is not False
            or raw.get("local_file_path") != str(path.resolve())
            or hdr.get("status") != "ELG_FULL_HEADER_ONLY_PARSED_REFERENCE_NOT_VERIFIED"
            or hdr.get("source_url") != p["source_url"]
            or hdr.get("header_range_sha256") != p["expected_header_sha256"]
            or hdr.get("header_bytes_requested") != p["expected_header_bytes"]
            or hdr.get("remote_bytes_reported") != p["expected_source_bytes"]
            or hdr.get("header_declared_rows_NOT_READ") != p["expected_declared_rows"]
            or hdr.get("header_declared_row_bytes") != p["expected_row_bytes"]
            or hdr.get("column_count") != p["expected_column_count"]
            or hdr.get("has_mskbit_column") is not True):
        raise ValueError("Frozen bit8 protocol / local SHA provenance / FITS header mismatch")
    if path.stat().st_size != p["expected_source_bytes"]:
        raise ValueError("Local FITS byte count differs from frozen release")


def hash_full_and_read_header(path, p):
    full = hashlib.sha256()
    first = hashlib.sha256()
    total = 0
    head_size = p["expected_header_bytes"]
    with path.open("rb") as f:
        while True:
            block = f.read(1024 * 1024)
            if not block:
                break
            full.update(block)
            if total < head_size:
                first.update(block[:min(len(block), head_size - total)])
            total += len(block)
            if total > p["expected_source_bytes"]:
                raise ValueError("Local source exceeded exact frozen byte cap")
        if total != p["expected_source_bytes"]:
            raise ValueError("Local source ended before frozen byte count")
        if full.hexdigest() != p["first_observed_full_file_sha256_from_local_byte_only_run"]:
            raise ValueError("Full-file SHA256 mismatch; do not read FITS data rows")
        if first.hexdigest() != p["expected_header_sha256"]:
            raise ValueError("First-header SHA256 mismatch; do not read FITS data rows")
        f.seek(0)
        header_bytes = f.read(head_size)
    return header_bytes, full.hexdigest()


def width(tform):
    match = FORM.fullmatch(tform.strip().upper())
    if not match:
        raise ValueError("Unsupported FITS TFORM in pinned source")
    count = int(match.group(1) or "1")
    kind = match.group(2)
    if kind == "X":
        return (count + 7) // 8
    if kind not in WIDTHS:
        raise ValueError("Unsupported FITS TFORM type in pinned source")
    return count * WIDTHS[kind]


def column_offsets(names, forms, row_bytes):
    if len(names) != len(forms) or len(set(x.upper() for x in names)) != len(names):
        raise ValueError("Invalid/duplicate FITS TTYPE fields")
    offsets = {}
    pos = 0
    for name, form in zip(names, forms):
        offsets[name.upper()] = (pos, form.upper())
        pos += width(form)
    if pos != row_bytes:
        raise ValueError("FITS TFORM widths do not match pinned NAXIS1 stride")
    for name, form in (("RA", "D"), ("DEC", "D"), ("MSKBIT", "I")):
        if name not in offsets or offsets[name][1] != form:
            raise ValueError("Required FITS RA/DEC/mskbit column format differs")
    return offsets


def parse_pinned_header(data, p, hdr):
    primary, first_end = parse_header(data, 0)
    if primary.get("SIMPLE") is not True or hdu_data_length(primary) != 0:
        raise ValueError("Primary FITS HDU is not the pinned empty primary")
    binary, data_begin = parse_header(data, first_end)
    if (str(binary.get("XTENSION", "")).upper() != "BINTABLE"
            or data_begin != p["expected_header_bytes"]
            or int(binary["NAXIS1"]) != p["expected_row_bytes"]
            or int(binary["NAXIS2"]) != p["expected_declared_rows"]
            or int(binary["TFIELDS"]) != p["expected_column_count"]):
        raise ValueError("FITS table geometry differs from previously pinned header")
    names = [str(binary.get(f"TTYPE{i}", "")).strip()
             for i in range(1, p["expected_column_count"] + 1)]
    forms = [str(binary.get(f"TFORM{i}", "")).strip().upper()
             for i in range(1, p["expected_column_count"] + 1)]
    if names != hdr["column_names"] or forms != hdr["column_formats"]:
        raise ValueError("Current FITS column structure differs from frozen header report")
    offsets = column_offsets(names, forms, p["expected_row_bytes"])
    expected_data_end = data_begin + p["expected_declared_rows"] * p["expected_row_bytes"]
    if expected_data_end > p["expected_source_bytes"]:
        raise ValueError("Table bytes exceed pinned official file length")
    return data_begin, offsets


def isolated_fields(buffer, data_begin, nrows, stride, offsets):
    """Only decode these three FITS fields. Never construct a whole FITS table."""
    specs = (("RA", ">f8"), ("DEC", ">f8"), ("MSKBIT", ">i2"))
    result = []
    for name, dtype in specs:
        offset, _ = offsets[name]
        view = np.ndarray(
            shape=(nrows,), dtype=np.dtype(dtype), buffer=buffer,
            offset=data_begin + offset, strides=(stride,),
        )
        result.append(view.astype(np.float64 if name != "MSKBIT" else np.int64))
    return tuple(result)


def compare(ra, dec, mask, pixels, hp, p):
    if (ra.ndim != 1 or ra.shape != dec.shape or ra.shape != mask.shape
            or np.any(~np.isfinite(ra)) or np.any(~np.isfinite(dec))
            or np.any(ra < 0) or np.any(ra >= 360)
            or np.any(dec < -90) or np.any(dec > 90)):
        raise ValueError("Invalid RA/DEC range or shape; no coordinate repair allowed")
    labels = np.bitwise_and(mask, p["bit8_integer_mask"]) != 0
    candidates = candidate_pixels(ra, dec, hp)
    if list(candidates) != p["candidate_coordinate_conventions"]:
        raise ValueError("Candidate mapping list differs from frozen source-only audit")
    if not np.array_equal(
        candidates["native_RA_radians"], candidates["native_RA_degrees"]
    ):
        raise ValueError("Independent native-RA radian/degree implementations disagree")
    membership = {}
    truth = int(np.count_nonzero(labels))
    for name, candidate in candidates.items():
        pred = np.isin(candidate, pixels)
        membership[name] = {
            "label_positive_and_candidate_positive_TP": int(np.count_nonzero(labels & pred)),
            "label_negative_and_candidate_positive_FP": int(np.count_nonzero(~labels & pred)),
            "label_positive_and_candidate_negative_FN": int(np.count_nonzero(labels & ~pred)),
            "label_negative_and_candidate_negative_TN": int(np.count_nonzero(~labels & ~pred)),
            "candidate_positive_total": int(np.count_nonzero(pred)),
        }
        counts = membership[name]
        if sum(counts[k] for k in (
            "label_positive_and_candidate_positive_TP",
            "label_negative_and_candidate_positive_FP",
            "label_positive_and_candidate_negative_FN",
            "label_negative_and_candidate_negative_TN",
        )) != len(ra):
            raise ValueError("Confusion counts do not account for every published row")
    return {
        "status": ("OFFICIAL_ELG_FULL_BIT8_LABELS_PRESENT_GEOMETRY_AUDITED_ONLY"
                   if truth > 0 else
                   "OFFICIAL_ELG_FULL_NO_BIT8_POSITIVE_LABELS_REFERENCE_UNRESOLVED"),
        "catalogue_header_declared_and_read_rows": int(len(ra)),
        "bit8_positive_catalogue_rows": truth,
        "published_removed_targets_diagnostic_only":
            p["published_bit8_removed_targets_for_diagnostic_only"],
        "bit8_count_equals_published_removed_targets_DIAGNOSTIC_ONLY": (
            truth == p["published_bit8_removed_targets_for_diagnostic_only"]
        ),
        "candidate_confusion_counts": membership,
        "native_radians_degrees_identical_all_rows": True,
        "published_source_pixel_count": int(len(pixels)),
        "positive_labels_present": truth > 0,
        "independent_official_production_chain_authenticated": False,
        "official_bit8_production_coordinate_convention_certified": False,
        "official_mask_applied_to_science_selection": False,
        "physical_LRG_ELG_pair_window_certified": False,
        "observed_row_columns_decoded": ["RA", "DEC", "mskbit"],
        "observed_redshift_or_weight_values_decoded": False,
        "observed_target_identifiers_decoded": False,
        "observed_random_or_mock_rows_decoded": False,
        "observed_odd_data_vector_read": False,
        "note": (
            "Only aggregate agreement of the published full-catalogue bit8 labels "
            "against four predeclared published-source geometries. Even perfect "
            "agreement identifies the catalogue mapping at its sampled targets, "
            "NOT an independently authenticated production pipeline, a full "
            "continuous angular mask, a tracer-pair window or an odd-sector detection."
        ),
    }


def audit(p, hp):
    source = ROOT / p["source_file"]
    raw = read_json(RAW_REPORT)
    hdr = read_json(HEADER_REPORT)
    parent = read_json(PARENT)
    verify_contract(p, parent, raw, hdr, source)
    header, digest = hash_full_and_read_header(source, p)
    begin, offsets = parse_pinned_header(header, p, hdr)
    _, pixels, upstream_blob, _ = frozen_inputs()
    expected = p["source_pixels"]
    if not expected.startswith("Exactly 37") or len(pixels) != 37:
        raise ValueError("Pinned upstream 37-pixel source contract changed")
    mmap = np.memmap(source, mode="r", dtype="u1")
    try:
        ra, dec, mskbit = isolated_fields(
            mmap, begin, p["expected_declared_rows"],
            p["expected_row_bytes"], offsets,
        )
        result = compare(ra, dec, mskbit, pixels, hp, p)
    finally:
        del mmap
    result.update({
        "full_file_sha256_reverified": digest,
        "first_header_sha256_reverified": p["expected_header_sha256"],
        "upstream_git_blob_sha": upstream_blob,
        "nside": NSIDE, "ordering": "RING",
        "healpy_version": hp.__version__,
        "errors": [],
    })
    return result


def self_test(hp):
    p = read_json(PROTOCOL)
    assert p["first_observed_full_file_sha256_from_local_byte_only_run"] == (
        "8806699e14422904ef91efb6f1632184171fb740c2e074456c0534b7d2dc28b3"
    )
    assert p["allowed_observed_row_columns"] == ["RA", "DEC", "mskbit"]
    offs = column_offsets(["RA", "DEC", "IGNORED", "MSKBIT"],
                          ["D", "D", "4A", "I"], 22)
    import struct
    raw = (struct.pack(">dd4sh", 0., 0., b"XXXX", 0)
           + struct.pack(">dd4sh", 180., 0., b"YYYY", 256))
    ra, dec, m = isolated_fields(raw, 0, 2, 22, offs)
    assert ra.tolist() == [0., 180.] and dec.tolist() == [0., 0.]
    assert m.tolist() == [0, 256]
    _, pixels, _, _ = frozen_inputs()
    positive_ra, positive_dec = hp.pix2ang(
        NSIDE, pixels[0], nest=False, lonlat=True)
    ra = np.array([float(positive_ra), 180.0], dtype="f8")
    dec = np.array([float(positive_dec), 0.0], dtype="f8")
    mask = np.array([256, 0], dtype="i8")
    res = compare(ra, dec, mask, pixels, hp, p)
    c = res["candidate_confusion_counts"]["native_RA_degrees"]
    assert c["label_positive_and_candidate_positive_TP"] == 1
    assert c["label_negative_and_candidate_negative_TN"] == 1
    assert res["observed_odd_data_vector_read"] is False
    print("EBOSS_ELG_BIT8_THREE_COLUMN_LABEL_AUDIT_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    import healpy as hp

    if args.self_test:
        self_test(hp)
        return 0
    p = read_json(PROTOCOL)
    out = args.out or ROOT / p["output"]
    try:
        result = audit(p, hp)
    except Exception as exc:
        result = {
            "status": "OFFICIAL_ELG_FULL_BIT8_LABEL_AUDIT_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "official_bit8_production_coordinate_convention_certified": False,
            "observed_odd_data_vector_read": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(".tmp.json")
    temp.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temp.replace(out)
    print("EBOSS_ELG_BIT8_FULL_CATALOGUE_AUDIT", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result.get("errors"):
        print("ERRORS", *result["errors"], sep="\n", flush=True)
        return 2
    print("BIT8_POSITIVE_CATALOGUE_ROWS",
          result["bit8_positive_catalogue_rows"], flush=True)
    for name, counts in result["candidate_confusion_counts"].items():
        print("MAPPING", name, json.dumps(counts, sort_keys=True), flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
