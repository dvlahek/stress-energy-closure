#!/usr/bin/env python3
"""Blinded source-only audit of ELG bit-8 RA orientation.

This is a follow-up to the upstream lonlat/radian diagnostic. The previous
37/37 roundtrip inverted the *published* phi = 360 - RA expression and thus
did not test the native celestial-RA interpretation of the listed pixels.
This script tests both readings without accessing any survey catalogue.
It NEVER applies bit 8 or identifies the official production convention.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_eboss_dr16_elg_bit8_coordinate_contract import (
    upstream_pixels_and_identity,
)

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = (
    ROOT
    / "source_data/eboss_dr16_elg_bit8_ra_orientation_protocol_2026-09-25.json"
)
NSIDE = 1024


def frozen_inputs():
    raw = PROTOCOL.read_bytes()
    protocol = json.loads(raw)
    pixels, source_blob = upstream_pixels_and_identity()
    source = protocol["published_upstream"]
    if (
        source_blob != source["git_blob_sha"]
        or source["source_pixels_sha256"]
        != hashlib.sha256(pixels.tobytes()).hexdigest()
        or protocol["nside"] != NSIDE
        or protocol["ordering"] != "RING"
        or len(pixels) != protocol["published_source_pixel_count"]
        or protocol["official_reference_output_compared"] is not False
        or protocol["official_bit8_applied"] is not False
        or protocol["observed_odd_data_vector_read"] is not False
    ):
        raise ValueError("RA-orientation source/protocol contract changed")
    return protocol, pixels, source_blob, hashlib.sha256(raw).hexdigest()


def candidate_pixels(ra, dec, hp):
    theta = np.radians(90.0 - dec)
    phi_reflected = np.radians(360.0 - ra)
    phi_native = np.radians(ra)
    return {
        "literal_upstream": np.asarray(
            hp.ang2pix(NSIDE, theta, phi_reflected, nest=False, lonlat=True),
            dtype="<i8",
        ),
        "units_only_reflected_RA": np.asarray(
            hp.ang2pix(NSIDE, theta, phi_reflected, nest=False, lonlat=False),
            dtype="<i8",
        ),
        "native_RA_radians": np.asarray(
            hp.ang2pix(NSIDE, theta, phi_native, nest=False, lonlat=False),
            dtype="<i8",
        ),
        "native_RA_degrees": np.asarray(
            hp.ang2pix(NSIDE, ra, dec, nest=False, lonlat=True),
            dtype="<i8",
        ),
    }


def rectangle_counts(ra, dec, rectangles):
    n = rectangles["NGC"]
    s = rectangles["SGC"]
    ngc = (
        (ra > n["ra_deg"][0])
        & (ra < n["ra_deg"][1])
        & (dec > n["dec_deg"][0])
        & (dec < n["dec_deg"][1])
    )
    sgc = (
        ((ra > s["ra_deg_wrap"][0]) | (ra < s["ra_deg_wrap"][1]))
        & (dec > s["dec_deg"][0])
        & (dec < s["dec_deg"][1])
    )
    if np.any(ngc & sgc):
        raise ValueError("Approximate cap rectangles unexpectedly overlap")
    return {
        "NGC": int(np.count_nonzero(ngc)),
        "SGC": int(np.count_nonzero(sgc)),
        "outside": int(np.count_nonzero(~(ngc | sgc))),
    }


def audit(hp):
    protocol, pixels, source_blob, protocol_sha = frozen_inputs()
    # healpy lonlat=True returns longitude/latitude in DEGREES, and its
    # longitude is here the native HEALPix angle phi. Interpret it as RA
    # only for the declared source/footprint compatibility diagnostic.
    ra_native, dec = hp.pix2ang(NSIDE, pixels, nest=False, lonlat=True)
    ra_native = np.mod(np.asarray(ra_native, dtype="f8"), 360.0)
    dec = np.asarray(dec, dtype="f8")
    versions = candidate_pixels(ra_native, dec, hp)

    np.testing.assert_array_equal(versions["native_RA_radians"], pixels)
    np.testing.assert_array_equal(versions["native_RA_degrees"], pixels)
    np.testing.assert_array_equal(
        versions["native_RA_radians"], versions["native_RA_degrees"]
    )
    listed = set(int(x) for x in pixels)
    matches = {
        name: {
            "exact_source_pixel_centre_roundtrip": int(
                np.count_nonzero(values == pixels)
            ),
            "returned_pixel_in_37_source_pixels": int(
                sum(int(x) in listed for x in values)
            ),
        }
        for name, values in versions.items()
    }

    rectangles = protocol["published_reference"][
        "approximate_survey_rectangles_for_diagnostic_only"
    ]
    native_rect = rectangle_counts(ra_native, dec, rectangles)
    reflected_rect = rectangle_counts(
        np.mod(360.0 - ra_native, 360.0), dec, rectangles
    )
    exploratory = protocol["exploratory_values_disclosed_before_frozen_followup"]
    if (
        native_rect != exploratory["native_RA_rectangles"]
        or reflected_rect != exploratory["reflected_RA_rectangles"]
        or matches["units_only_reflected_RA"][
            "returned_pixel_in_37_source_pixels"
        ] != exploratory["native_centre_reflected_mapping_source_memberships"]
    ):
        raise ValueError("Frozen follow-up geometry differs from exploratory record")

    area = float(len(pixels) * hp.nside2pixarea(NSIDE, degrees=True))
    if round(area, 1) != protocol["published_reference"][
        "reported_bit8_area_deg2_rounded"
    ]:
        raise ValueError("Listed source pixels disagree with rounded paper area")
    examples = [
        {
            "source_pixel": int(pixels[i]),
            "native_HEALPix_longitude_deg": float(ra_native[i]),
            "latitude_deg": float(dec[i]),
            "reflected_longitude_deg": float(
                np.mod(360.0 - ra_native[i], 360.0)
            ),
            "literal_upstream_pixel": int(versions["literal_upstream"][i]),
            "units_only_reflected_pixel": int(
                versions["units_only_reflected_RA"][i]
            ),
            "native_RA_pixel": int(versions["native_RA_radians"][i]),
        }
        for i in range(min(5, len(pixels)))
    ]
    return {
        "status": "published_source_bit8_ra_orientation_diagnostic_only",
        "upstream_git_blob_sha": source_blob,
        "protocol_sha256": protocol_sha,
        "healpy_version": hp.__version__,
        "nside": NSIDE,
        "ordering": "RING",
        "source_pixel_count": int(len(pixels)),
        "source_pixel_array_sha256": hashlib.sha256(
            pixels.tobytes()
        ).hexdigest(),
        "candidate_counts_at_native_HEALPix_pixel_centres": matches,
        "approximate_published_survey_rectangles_NATIVE_ra": native_rect,
        "approximate_published_survey_rectangles_REFLECTED_ra": reflected_rect,
        "source_37_pixel_area_deg2": area,
        "paper_rounded_bit8_area_deg2": protocol["published_reference"][
            "reported_bit8_area_deg2_rounded"
        ],
        "paper_reported_removed_targets_not_tested": protocol[
            "published_reference"
        ]["reported_removed_targets"],
        "first_five_source_order_examples": examples,
        "independent_official_reference_bitmap_compared": False,
        "official_bit8_production_coordinate_convention_resolved": False,
        "official_bit8_applied": False,
        "physical_ELG_mask_certified": False,
        "physical_LRG_ELG_pair_window_certified": False,
        "observed_galaxy_catalogue_read": False,
        "mock_galaxy_catalogue_read": False,
        "observed_random_positions_read": False,
        "mock_random_positions_read": False,
        "observed_odd_data_vector_read": False,
        "inference_protocol_frozen": False,
        "note": (
            "37/37 native-RA and 37/37 reflected-RA results are different "
            "coordinate roundtrips at DIFFERENT constructed sky positions. "
            "Neither proves the celestial RA convention of the official "
            "DR16 production mask. Published bounding rectangles are "
            "approximate only and must never be used as a joint mask."
        ),
    }


def self_test(hp):
    pixels = np.array([111, 333333, 2981667, 6031493, 10500000], dtype="<i8")
    ra, dec = hp.pix2ang(NSIDE, pixels, nest=False, lonlat=True)
    result = candidate_pixels(
        np.asarray(ra, dtype="f8"), np.asarray(dec, dtype="f8"), hp
    )
    np.testing.assert_array_equal(result["native_RA_radians"], pixels)
    np.testing.assert_array_equal(result["native_RA_degrees"], pixels)
    if not np.any(result["units_only_reflected_RA"] != pixels):
        raise AssertionError("Synthetic reflected-longitude test did not discriminate")
    frozen_inputs()
    print("EBOSS_ELG_BIT8_RA_ORIENTATION_SELF_TEST_OK", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(
            "reports/eboss_dr16_elg_bit8_ra_orientation_followup.json"
        ),
    )
    args = parser.parse_args()
    import healpy as hp

    if args.self_test:
        self_test(hp)
        return 0
    report = audit(hp)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.out)
    print(
        "EBOSS_ELG_BIT8_RA_ORIENTATION_SOURCE_ONLY",
        report["source_pixel_count"],
        report["candidate_counts_at_native_HEALPix_pixel_centres"][
            "native_RA_radians"
        ]["exact_source_pixel_centre_roundtrip"],
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
