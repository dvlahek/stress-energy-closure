#!/usr/bin/env python3
"""Source-pinned HEALPix bit-8 coordinate-contract audit (NO MASK APPLICATION).

The published extra eBOSS ELG script passes radian-transformed angles to
healpy.ang2pix(..., lonlat=True), while healpy documents degree longitude
and latitude for that option. Extract the published mask_pixel array from
the byte-identical upstream Python AST and test all its HEALPix pixel
centres under the published call and two mathematically equivalent
co-latitude/longitude conventions. This is an independent geometry
diagnostic only; it does not select or remove any catalogue rows.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_coordinate_contract_protocol_2026-09-25.json"
UPSTREAM = ROOT / "scripts/eBOSS_ELG_extra.py"
UPSTREAM_GIT_BLOB = "d66c5e6ed6a8260b5e8dcc840795485345392ced"
NSIDE = 1024


def upstream_pixels_and_identity():
    raw = UPSTREAM.read_bytes()
    blob = hashlib.sha1(
        f"blob {len(raw)}\0".encode("ascii") + raw
    ).hexdigest()
    if blob != UPSTREAM_GIT_BLOB:
        raise ValueError("Published upstream ELG extra-mask source changed")
    tree = ast.parse(raw.decode("utf-8"), filename=str(UPSTREAM))
    assignments = [
        node.value for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "mask_pixel"
                for t in node.targets)
    ]
    if len(assignments) != 1:
        raise ValueError("Missing or repeated published mask_pixel assignment")
    source_pixels = ast.literal_eval(assignments[0])
    if (not isinstance(source_pixels, list) or not 0 < len(source_pixels) < 1000
            or len(set(source_pixels)) != len(source_pixels)
            or any(type(x) is not int or not 0 <= x < 12 * NSIDE**2
                   for x in source_pixels)):
        raise ValueError("Malformed published extra bit-8 HEALPix pixel list")
    source = raw.decode("utf-8")
    if ("hp.pixelfunc.ang2pix(1024, theta, phi, nest=False, lonlat=True)"
            not in source
            or "theta = np.radians(90 - dec)" not in source
            or "phi = np.radians(360 - ra)" not in source):
        raise ValueError("Published HEALPix coordinate source contract differs")
    return np.asarray(source_pixels, dtype="<i8"), blob


def mappings(ra_deg, dec_deg, hp):
    ra = np.asarray(ra_deg, dtype="f8")
    dec = np.asarray(dec_deg, dtype="f8")
    if (ra.shape != dec.shape or not np.isfinite(ra).all()
            or not np.isfinite(dec).all()
            or np.any(ra < 0) or np.any(ra >= 360)
            or np.any(dec < -90) or np.any(dec > 90)):
        raise ValueError("Invalid synthetic spherical coordinates")
    theta = np.radians(90.0 - dec)
    phi = np.radians(360.0 - ra)
    as_written = hp.ang2pix(NSIDE, theta, phi, nest=False, lonlat=True)
    colat_rad = hp.ang2pix(NSIDE, theta, phi, nest=False, lonlat=False)
    degrees = hp.ang2pix(
        NSIDE, (360.0 - ra) % 360.0, dec, nest=False, lonlat=True)
    return (np.asarray(as_written,dtype="<i8"),
            np.asarray(colat_rad,dtype="<i8"),
            np.asarray(degrees,dtype="<i8"))


def audit(hp):
    p = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    pixels, blob = upstream_pixels_and_identity()
    if (p.get("source_git_blob_sha") != blob
            or p.get("nside") != NSIDE or p.get("ordering") != "RING"
            or p.get("observed_odd_data_vector_read") is not False
            or p.get("official_bit8_applied") is not False
            or p.get("pixel_list") !=
            "Extract complete literal mask_pixel assignment from byte-identical upstream source with Python AST; no hand-selected pixel subset"):
        raise ValueError("Bit-8 coordinate diagnostic protocol changed")
    theta, phi = hp.pix2ang(NSIDE, pixels, nest=False, lonlat=False)
    ra = np.mod(360.0 - np.degrees(phi), 360.0)
    dec = 90.0 - np.degrees(theta)
    literal, radians, degrees = mappings(ra,dec,hp)
    matches_literal = literal == pixels
    matches_radians = radians == pixels
    matches_degrees = degrees == pixels
    same = radians == degrees
    disagreement = literal != radians
    if not np.all(matches_radians) or not np.all(matches_degrees) or not np.all(same):
        raise ValueError("Internal spherical-coordinate/HEALPix roundtrip failed")
    differing = np.flatnonzero(disagreement)[:10]
    return {
        "status": ("published_elg_bit8_source_coordinate_discrepancy_documented"
                   if np.any(disagreement)
                   else "published_elg_bit8_source_coordinate_conventions_agree_on_probe"),
        "source_upstream_commit": p["source_commit"],
        "source_git_blob_sha": blob,
        "mask_pixel_array_sha256": hashlib.sha256(pixels.tobytes()).hexdigest(),
        "nside": NSIDE, "ordering": "RING",
        "source_mask_pixel_count": int(len(pixels)),
        "fixed_probe": "All source-published RING pixel centres, RA=360-phi, DEC=90-theta",
        "published_literal_lonlat_true_matching_source_pixels": int(np.count_nonzero(matches_literal)),
        "radian_colatitude_longitude_matching_source_pixels":
            int(np.count_nonzero(matches_radians)),
        "equivalent_degree_longitude_latitude_matching_source_pixels":
            int(np.count_nonzero(matches_degrees)),
        "radian_vs_degree_disagreement_count":int(np.count_nonzero(~same)),
        "published_literal_vs_spherical_radians_disagreement_count":
            int(np.count_nonzero(disagreement)),
        "first_source_order_differences": [
            {"source_pixel":int(pixels[i]), "ra_deg":float(ra[i]),
             "dec_deg":float(dec[i]), "literal_lonlat_true_pixel":int(literal[i]),
             "radian_spherical_pixel":int(radians[i]),
             "degree_spherical_pixel":int(degrees[i])}
            for i in differing
        ],
        "mathematical_spherical_coordinate_roundtrip_verified":True,
        "official_reference_bitmap_compared":False,
        "intended_published_reference_mask_convention_identified":False,
        "official_bit8_applied":False,
        "official_physical_elg_mask_certified":False,
        "observed_galaxy_positions_read":False,
        "mock_galaxy_positions_read":False,
        "observed_random_positions_read":False,
        "mock_random_positions_read":False,
        "observed_odd_data_vector_read":False,
        "physical_joint_pair_window_certified":False,
        "note":"Identical published pixel centres expose only source/API coordinate consistency. Do not silently substitute either convention into published DR16 mask selection until checked against an official bitmap or verified output catalogue."
    }


def self_test(hp):
    points = np.asarray([111,333333,2981667,6031493,10500000],dtype="i8")
    theta,phi = hp.pix2ang(NSIDE,points,nest=False,lonlat=False)
    ra = np.mod(360.0-np.degrees(phi),360.0)
    dec = 90.0-np.degrees(theta)
    original,rad,deg = mappings(ra,dec,hp)
    np.testing.assert_array_equal(rad,points)
    np.testing.assert_array_equal(deg,points)
    assert np.count_nonzero(original != rad)>0
    pixels,sha=upstream_pixels_and_identity()
    assert len(pixels)>=30 and sha==UPSTREAM_GIT_BLOB
    print("EBOSS_ELG_BIT8_COORDINATE_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--out",type=Path,
                        default=Path("eboss_workspace/official_mask_inventory/bit8_coordinate_contract.json"))
    args=parser.parse_args()
    import healpy as hp
    if args.self_test:
        self_test(hp)
        return 0
    result=audit(hp)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    temp=args.out.with_suffix(".tmp.json")
    temp.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    temp.replace(args.out)
    print("EBOSS_ELG_BIT8_COORDINATE_CONTRACT",
          result["status"],result["source_mask_pixel_count"],
          result["published_literal_vs_spherical_radians_disagreement_count"],
          flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
