#!/usr/bin/env python3
"""A-03 published-source provenance stop gate, no FITS or observed odd reads.

This verifies that the exact archived 18/18 A-02 mock-galaxy algebra result
does NOT get silently promoted to a certified historical ELG production
mask or an observational LRG×ELG physical pair window. Read repository
metadata/archived JSON bytes only; never download or open FITS catalogues.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROTOCOL = ROOT / "source_data/eboss_dr16_a03_selection_window_provenance_decision_2026-09-26.json"
STOP = "EBOSS_A03_HISTORICAL_MASK_AND_EMPIRICAL_WINDOW_NOT_YET_CERTIFIED"


def read(path):
    return json.loads((ROOT/path).read_bytes())


def git_blob_sha(path):
    raw=(ROOT/path).read_bytes()
    return hashlib.sha1(("blob "+str(len(raw))).encode()+b"\0"+raw).hexdigest()


def validate_claims(p):
    ev=p["existing_selection_source_evidence"]
    a=p["completed_A02"]
    h=p["two_nonconflated_paths"]["historical_full_mask_reconstruction"]
    e=p["two_nonconflated_paths"]["published_random_empirical_pair_window"]
    if (
        p["master_plan_version"]!="1.2"
        or p["pr_number"]!=1
        or p["pr_is_draft"] is not True
        or p["must_not_update_main"] is not True
        or p["observed_odd_data_vector_read"] is not False
        or p["new_observed_galaxy_rows_read"] is not False
        or p["new_science_selection_applied"] is not False
        or p["physical_ELG_production_mask_certified"] is not False
        or p["physical_LRG_ELG_pair_window_certified"] is not False
        or p["valid_18D_inferential_covariance_certified"] is not False
        or p["unblinding_authorized"] is not False
        or p["historical_source_request"]["sent"] is not False
        or p["historical_source_request"]["transferred_reply_or_files_available"] is not False
        or h["status"]!="NOT_CERTIFIED"
        or e["status"]!="POTENTIAL_LIMITED_PATH_NOT_YET_VALIDATED"
        or a["completed_fixed_mock_id_cap_cases"]!=18
        or a["min_RR_supported_cells_per_case"]!=144
        or a["does_not_certify_production_mask_or_physical_window"] is not True
        or a["does_not_provide_inferential_18D_covariance"] is not True
        or ev["released_full_catalogue_bit8_native_RA_TP_FP_FN_TN"]!=[15,0,0,269163]
        or ev["released_full_catalogue_bit8_label_match_not_production_executable_proof"] is not True
        or ev["four_image_sample_FITS_not_full_ELG_19381_image_certification"] is not True
        or ev["plate_footprint_geometry_and_release_acceptance_not_proven"] is not True
        or ev["published_bad_plate_pairs"]!=[
            {"PLATE":9430,"MJD":58112},{"PLATE":9395,"MJD":58113}]
        or ev["later_public_README_credit"]!=
            "README acknowledges Arnaud de Mattia shared original extra maskbit eBOSS ELG script with public brickmask author. This is a lead for source request, not proof later helper is DR16 production original."
    ):
        raise ValueError("A-03 claim overstated or original frozen evidence changed")


def check_archived_references(p):
    ev=p["existing_selection_source_evidence"]
    a=p["completed_A02"]
    refs=(
        (a["uploaded_manifest"],a["uploaded_manifest_git_blob_sha1"]),
        (ev["production_lineage_gate"],ev["production_lineage_gate_git_blob_sha1"]),
        (ev["source_version_comparison"],ev["source_version_comparison_git_blob_sha1"]),
        (ev["released_ELG_bit8_label_audit_manifest"],
         ev["released_ELG_bit8_label_manifest_git_blob_sha1"]),
        (ev["published_bad_plate_identifier_reference"],
         ev["published_bad_plate_identifier_reference_git_blob_sha1"]),
    )
    for path,expected in refs:
        if git_blob_sha(path)!=expected:
            raise ValueError("Original immutable A-03/A-02 Git source changed: "+path)
    raw=(ROOT/a["report"]).read_bytes()
    if (len(raw)!=a["report_bytes"]
        or hashlib.sha256(raw).hexdigest()!=a["report_SHA256"]):
        raise ValueError("Original uploaded 18/18 mock source report changed")
    report=json.loads(raw)
    manifest=read(a["uploaded_manifest"])
    if (
        report["status"]!="EZMOCK_NINE_GALAXY_CROSS_LS_CODE_TRANSPORT_ALGEBRA_ONLY"
        or report["completed_cases"]!=18
        or report["failed_cases"]!=0
        or report["errors"]!=[]
        or report["all_72_full_compressed_gzip_SHA_checked_before_any_FITS_rows"] is not True
        or report["observed_odd_data_vector_read"] is not False
        or manifest["exact_uploaded_SHA256"]!=a["report_SHA256"]
        or manifest["physical_pair_window_certified"] is not False
        or manifest["inferential_18D_mock_covariance_computed"] is not False
        or any(c["status"]!="complete" or c["cross_ls"]["RR_supported_s_mu_cells"]!=144
               for c in report["cases"].values())
    ):
        raise ValueError("A-02 archived result promoted beyond its algorithmic scope")
    lineage=read(ev["production_lineage_gate"])
    version=read(ev["source_version_comparison"])
    bit8=read(ev["released_ELG_bit8_label_audit_manifest"])
    plates=read(ev["published_bad_plate_identifier_reference"])
    if (
        lineage["guardrails"]["physical_elg_mask_certified"] is not False
        or lineage["guardrails"]["physical_lrg_elg_pair_window_certified"] is not False
        or lineage["guardrails"]["inference_protocol_unblinding_authorized"] is not False
        or version["findings"]["legacy_VETOMASK_to_released_full_mskbit_pipeline_NOT_AUTHENTICATED"] is not True
        or version["findings"]["tagged_release_source_not_evidence_of_exact_historical_execution"] is not True
        or bit8["summary"]["native_RA_radians"]!=
             {"tp":15,"fp":0,"fn":0,"tn":269163}
        or bit8["physical_mask_certified"] is not False
        or plates["published_plate_mjd_pairs"]!=
             ev["published_bad_plate_pairs"]
        or plates["physical_lrg_elg_pair_window_certified"] is not False
    ):
        raise ValueError("Published-source/label/bad-plate provenance changed")
    return report


def self_test(p):
    validate_claims(p)
    check_archived_references(p)
    for path,invalid in (
        (("physical_ELG_production_mask_certified",),True),
        (("unblinding_authorized",),True),
        (("two_nonconflated_paths","historical_full_mask_reconstruction","status"),"CERTIFIED"),
        (("historical_source_request","sent"),True),
    ):
        fake=copy.deepcopy(p)
        target=fake
        for member in path[:-1]:
            target=target[member]
        target[path[-1]]=invalid
        try:
            validate_claims(fake)
        except ValueError:
            pass
        else:
            raise AssertionError("Synthetic unsupported A-03 certification accepted: "+str(path))
    print("EBOSS_A03_SOURCE_ONLY_PROVENANCE_NEGATIVE_CONTROLS_OK",flush=True)
    print(STOP,flush=True)
    print("OBSERVED_ODD_DATA_READ",False,flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    p=json.loads(PROTOCOL.read_bytes())
    if args.self_test:
        self_test(p)
    else:
        validate_claims(p)
        check_archived_references(p)
        print(STOP,flush=True)
        print("REQUEST_SENT",False,flush=True)
        print("PHYSICAL_PAIR_WINDOW_CERTIFIED",False,flush=True)
        print("OBSERVED_ODD_DATA_READ",False,flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
