#!/usr/bin/env python3
"""Offline historical eBOSS ELG public extra-mask source-lineage gate.

Only reads pinned vendored Python *text* and already archived JSON
manifests/aggregate bit8 label report. Never imports or executes the
vendor helper; no FITS, galaxy, random, pixel, polygon geometry, odd
measurement, remote download or catalogue selection.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_production_lineage_gate_2026-09-26.json"
BIT8_REPORT = ROOT / "source_data/eboss_dr16_elg_full_bit8_label_geometry_audit_2026-09-25.json"
BIT8_MANIFEST = ROOT / "source_data/eboss_dr16_elg_full_bit8_label_uploaded_report_manifest_2026-09-25.json"
POLYGONS = ROOT / "source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json"
COMPLETE = "ELG_PUBLIC_EXTRA_HELPER_SOURCE_AUDITED_HISTORICAL_PRODUCTION_OPEN"


def git_blob(raw):
    return hashlib.sha1(
        ("blob " + str(len(raw))).encode("ascii") + bytes([0]) + raw
    ).hexdigest()


def extract_assign(function, name):
    found = [
        n.value for n in ast.walk(function)
        if isinstance(n, ast.Assign)
        and len(n.targets) == 1
        and isinstance(n.targets[0], ast.Name)
        and n.targets[0].id == name
    ]
    if len(found) != 1:
        raise ValueError("Missing or duplicate vendor AST assignment: " + name)
    return found[0]


def inspect_source(raw, p):
    upstream = p["public_upstream"]
    if git_blob(raw) != upstream["vendored_extra_source_git_blob_sha1"]:
        raise ValueError("Pinned public eBOSS ELG helper bytes changed")
    try:
        tree = ast.parse(raw)
    except SyntaxError as exc:
        raise ValueError("Published helper source is syntactically invalid") from exc
    functions = [
        x for x in tree.body
        if isinstance(x, ast.FunctionDef) and x.name == "eBOSS_ELG_mask"
    ]
    if len(functions) != 1:
        raise ValueError("Published extra-mask helper function is not unique")
    fun = functions[0]
    arguments = [x.arg for x in fun.args.args]
    if arguments != ["ra", "dec", "mskdir", "mskbit"]:
        raise ValueError("Extra-mask helper input contract changed")
    pixels = ast.literal_eval(extract_assign(fun, "mask_pixel"))
    plys = ast.literal_eval(extract_assign(fun, "plys"))
    codes = ast.literal_eval(extract_assign(fun, "ply_codes"))
    if (
        len(pixels) != upstream["bit8_source_pixel_count"]
        or len(set(pixels)) != len(pixels)
        or not all(isinstance(k, int) and k >= 0 for k in pixels)
        or plys != list(upstream["extra_polygon_map"])
        or codes != list(upstream["extra_polygon_map"].values())
    ):
        raise ValueError("Upstream 37-pixel list or bits 9-11 polygon assignment changed")
    theta = ast.unparse(extract_assign(fun, "theta"))
    phi = ast.unparse(extract_assign(fun, "phi"))
    pix = ast.unparse(extract_assign(fun, "pix"))
    bit8_assign = [
        ast.unparse(n)
        for n in ast.walk(fun)
        if isinstance(n, ast.AugAssign) and isinstance(n.target, ast.Name)
        and n.target.id == "mask"
    ]
    polygon_bit_assignments = [ast.unparse(n.value) for n in ast.walk(fun)
                               if isinstance(n, ast.Assign) and len(n.targets) == 1
                               and isinstance(n.targets[0], ast.Name)
                               and n.targets[0].id == "bit"
                               and isinstance(n.value, ast.Compare)]
    if len(polygon_bit_assignments) != 1:
        raise ValueError("Published helper polygon membership bit assignment changed")
    polygon_assign = polygon_bit_assignments[0]
    loops = [
        ast.unparse(n.iter) for n in ast.walk(fun)
        if isinstance(n, ast.For)
    ]
    expected_source={
        "theta": "np.radians(90 - dec)",
        "phi": "np.radians(360 - ra)",
        "pix": "hp.pixelfunc.ang2pix(1024, theta, phi, nest=False, lonlat=True)",
    }
    if (
        (theta, phi, pix) != tuple(expected_source.values())
        or sorted(bit8_assign) != sorted([
            "mask += bit * 2 ** 8", "mask += bit * 2 ** code"
        ])
        or "zip(plys, ply_codes)" not in loops
        or "m.polyid(ra, dec) != -1" != polygon_assign
    ):
        raise ValueError("Public extra-helper coordinate or bit-setting source changed")
    return {
        "vendor_source_blob_sha1": git_blob(raw),
        "published_bit8_healpix_list_length_NOT_OUTPUT": len(pixels),
        "published_helper_literal_bit8_formula": expected_source,
        "public_extra_polygon_bit_positions": dict(zip(plys,codes)),
        "public_source_unmodified": True,
        "historical_execution_or_RA_convention_authenticated": False,
    }


def check_earlier_metadata(p):
    raw, prior_manifest = BIT8_REPORT.read_bytes(), json.loads(
        BIT8_MANIFEST.read_text(encoding="utf-8"))
    earlier = p["existing_catalogue_only_bit8_evidence"]
    if (
        hashlib.sha256(raw).hexdigest() != earlier["archive_report_sha256"]
        or prior_manifest["uploaded_report_sha256"] != earlier["archive_report_sha256"]
        or git_blob(raw) != prior_manifest["uploaded_report_git_blob_sha1"]
    ):
        raise ValueError("Archived bit8 released-label evidence is not byte-identical")
    report = json.loads(raw)
    native = report["candidate_confusion_counts"]["native_RA_radians"]
    literal = report["candidate_confusion_counts"]["literal_upstream"]
    if (
        report["status"] != "OFFICIAL_ELG_FULL_BIT8_LABELS_PRESENT_GEOMETRY_AUDITED_ONLY"
        or report["catalogue_header_declared_and_read_rows"] !=
           earlier["full_catalogue_rows_checked"]
        or report["bit8_positive_catalogue_rows"] != earlier["bit8_label_positives"]
        or report["full_file_sha256_reverified"] != earlier["full_official_FITS_sha256"]
        or report["observed_odd_data_vector_read"] is not False
        or report["independent_official_production_chain_authenticated"] is not False
        or tuple(native[k] for k in (
            "label_positive_and_candidate_positive_TP",
            "label_negative_and_candidate_positive_FP",
            "label_positive_and_candidate_negative_FN",
            "label_negative_and_candidate_negative_TN"
        )) != tuple(earlier["native_RA_candidate_counts"][k]
                     for k in ("TP","FP","FN","TN"))
        or tuple(literal[k] for k in (
            "label_positive_and_candidate_positive_TP",
            "label_negative_and_candidate_positive_FP",
            "label_positive_and_candidate_negative_FN",
            "label_negative_and_candidate_negative_TN"
        )) != tuple(earlier["literal_upstream_candidate_counts"][k]
                     for k in ("TP","FP","FN","TN"))
    ):
        raise ValueError("Existing full-catalogue released bit8 label audit changed")
    polygons = json.loads(POLYGONS.read_text(encoding="utf-8"))
    published = {
        item["filename"]:item["sha256"]
        for item in polygons["products"]
        if item["role"] == "elg_extra_polygon"
    }
    expected = {
        "ELG_centerpost.ply":"2a9494725c9265894b9b463d6ccc8616d4979ed2c7bbc213764bf81b94e0b4c9",
        "ELG_TDSSFES_62arcsec.pix.snap.balk.ply":"ad1e4c135c351008eec316d8ff7e95f749ffd1b1b207fdcbd5ebae2181bcee0c",
        "ebosselg_badphot.26Aug2019.ply":"80aacd51459be796ccc493bad0450c9f8c26b0de7005dd2e3e5bfaadd0095a8d"
    }
    if published != expected or set(published) != set(
        p["public_upstream"]["extra_polygon_map"]
    ):
        raise ValueError("Already pinned official extra polygon bytes changed")
    return published


def synthetic_test(p):
    vendor=(ROOT/p["public_upstream"]["vendored_extra_source_path"]).read_bytes()
    positive=inspect_source(vendor,p)
    assert positive["published_bit8_healpix_list_length_NOT_OUTPUT"]==37
    assert list(positive["public_extra_polygon_bit_positions"].values())==[9,10,11]
    for invalid in (
        vendor.replace(b"ply_codes = [9, 10, 11]",
                       b"ply_codes = [9, 11, 10]"),
        vendor.replace(b"lonlat=True",b"lonlat=False"),
        vendor.replace(b"np.radians(360 - ra)",b"np.radians(ra)")
    ):
        assert invalid != vendor
        try:
            inspect_source(invalid,p)
        except ValueError:
            pass
        else:
            raise AssertionError("Changed vendor mask-code/RA source was accepted")
    print("EBOSS_ELG_PUBLIC_PRODUCTION_SOURCE_AST_SYNTHETIC_SELF_TEST_OK",
          flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    p=json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if args.self_test:
        synthetic_test(p)
        return 0
    dest=args.out or ROOT/"eboss_workspace/official_mask_inventory/elg_public_source_lineage_gate.json"
    try:
        contract=inspect_source(
            (ROOT/p["public_upstream"]["vendored_extra_source_path"]).read_bytes(),p
        )
        polygon_sha=check_earlier_metadata(p)
        result={
            "status":COMPLETE,
            "public_helper_source_contract":contract,
            "official_three_extra_polygon_sha256_from_prior_manifest":polygon_sha,
            "published_full_catalogue_bit8_label_crosscheck_from_exact_prior_report":{
                "source_report_sha256":
                    p["existing_catalogue_only_bit8_evidence"]["archive_report_sha256"],
                "full_catalogue_rows":
                    p["existing_catalogue_only_bit8_evidence"]["full_catalogue_rows_checked"],
                "bit8_positive_labels":
                    p["existing_catalogue_only_bit8_evidence"]["bit8_label_positives"],
                "native_RA_radians_confusion_counts":
                    p["existing_catalogue_only_bit8_evidence"]["native_RA_candidate_counts"],
            },
            "historical_dr16_production_invocation_authenticated":False,
            "historical_custom_low_quality_eboss22_plate_code_authenticated":False,
            "released_random_full_mask_semantics_authenticated":False,
            "physical_lrg_elg_pair_window_certified":False,
            "no_FITS_pixel_galaxy_random_or_mock_row_read":True,
            "observed_odd_data_vector_read":False,
            "new_mask_or_selection_applied":False,
            "errors":[],
        }
    except Exception as exc:
        result={
            "status":"ELG_PUBLIC_EXTRA_SOURCE_LINEAGE_INCOMPLETE_STOP",
            "errors":[str(exc)],
            "observed_odd_data_vector_read":False,
            "new_mask_or_selection_applied":False,
        }
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result,indent=2)+chr(10),encoding="utf-8")
    tmp.replace(dest)
    print("EBOSS_ELG_PUBLIC_SOURCE_LINEAGE",result["status"],flush=True)
    print("REPORT",dest,flush=True)
    if result["errors"]:
        print("ERRORS",*result["errors"],sep=chr(10),flush=True)
        return 2
    print("PUBLISHED_EXTRA_POLYGON_BITS",
          json.dumps(result["public_helper_source_contract"][
              "public_extra_polygon_bit_positions"],sort_keys=True),flush=True)
    print("PUBLISHED_BIT8_LABELS",
          result["published_full_catalogue_bit8_label_crosscheck_from_exact_prior_report"][
              "bit8_positive_labels"],flush=True)
    print("HISTORICAL_PRODUCTION_CERTIFIED",
          result["historical_dr16_production_invocation_authenticated"],flush=True)
    print("ODD_DATA_READ",result["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
