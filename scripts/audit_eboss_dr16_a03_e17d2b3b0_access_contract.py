#!/usr/bin/env python3
"""E17D2b3b0: OFFLINE Abacus c000 accession inventory and fail-closed dry-run.

This module DOES NOT contain a downloader, HTTP client, ASDF reader or
halo/tree loader. A successful CI result certifies only a provenance-locked
BLOCKED data-access checklist. A user must separately preregister actual
provider item access and independently verify file bytes and actual headers
before claiming real M200c or any physical halo response.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
B="audit/eboss-elg-bit8-ra-orientation-20260925"
PROTO=ROOT/"source_data/eboss_dr16_a03_e17d2b3b0_abacus_real_data_access_dryrun_prereg_2026-09-28.json"
PROTO_BLOB="49df0cd9c110b3ef9db02a3ad0cb9a4ff20554e1"
TEMPLATE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b0_abacus_access_candidate_BLOCKED_template_2026-09-28.json"
PARENT=ROOT/"source_data/eboss_dr16_a03_e17d2b3a_archived_CI_2026_09_28"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
B2=ROOT/"source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28/e17d2b2_original_abacus_c000_documented_metadata_bridge.json"
B3A=PARENT/"e17d2b3a_original_l1_200critical_nonidentifiability.json"
B3AI=PARENT/"e17d2b3a_independent_decimal_bisection_200critical_replay.json"
B3AMAN=PARENT/"archive_manifest.json"
B3APRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3a_mass_reference_nonidentifiability_prereg_2026-09-28.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3b0_abacus_access_dryrun"
REPORT="e17d2b3b0_offline_real_abacus_access_contract_BLOCKED.json"

def require(ok,message):
    if not ok:raise ValueError("E17D2B3B0_FAIL_CLOSED: "+message)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def write_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"refuse to overwrite a changed report: "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3b0_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B3B0_OFFLINE_REPORT",str(path),"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PROTO.read_bytes())==PROTO_BLOB,"prospective no-download protocol Git blob")
    p=json.loads(PROTO.read_bytes());parents=p["immutable_science_parents"]
    for path,key in ((E8,"E8_frozen_4000q_sha256"),
                     (E16,"E16_frozen_576_triangle_sha256"),
                     (E17A,"E17A_CLASS_original_joint_sha256"),
                     (B2,"E17D2b2_official_metadata_original_sha256"),
                     (B3A,"E17D2b3a_original_l1_m200c_proof_sha256"),
                     (B3AI,"E17D2b3a_independent_proof_sha256")):
        require(sha(path.read_bytes())==parents[key],"immutable original SHA "+key)
    for path,key in ((B3AMAN,"E17D2b3a_SHA_archive_manifest_git_blob"),
                     (B3APRE,"E17D2b3a_prereg_git_blob")):
        require(blob(path.read_bytes())==parents[key],"immutable parent Git blob "+key)
    z16=json.loads(E16.read_bytes())
    za=json.loads(E17A.read_bytes())
    prior=json.loads(B3A.read_bytes())
    b2=json.loads(B2.read_bytes())
    require(z16["QA"]["original_72_cases"]==72
            and z16["QA"]["original_24_contrasts"]==24
            and za["geometry_count"]==576
            and za["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED"
            and za["eBOSS_observed_odd_read"] is False,
            "frozen E16/E17A source/contrast/observed science scope")
    require(prior["actual_Abacus_halo_particle_or_merger_tree_read"] is False
            and prior["any_real_halo_M200c_remeasured"] is False
            and prior["observed_odd_SEALED"] is True
            and prior["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and b2["real_ASDF_header_or_halo_or_merger_tree_loaded"] is False
            and b2["observed_odd_SEALED"] is True
            and b2["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED",
            "immutable E17D2b2/E17D2b3a real-halo hard STOP")
    require(all(p["absolute_STOP"].values()),"future data prereg absolute no-read STOP")
    expected=p["prospective_locked_dry_run_candidate_template"]
    t=json.loads(TEMPLATE.read_bytes())
    require(t==expected,"exact prospective BLOCKED template must equal preregistered object")
    require(t["real_data_access_authorized"] is False
            and t["release_selection"]=="UNSELECTED_RAW_OR_CLEANED"
            and t["actual_halo_ASDF_header_read"] is False
            and t["actual_ASDF_redshift"] is None
            and t["observed_odd_SEALED"] is True
            and t["any_real_ASDF_or_halo_catalogue_or_merger_tree_downloaded_in_this_stage"] is False
            and t["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED",
            "no fabricated actual provider item, header or physical B")
    return p,t

def assess(template):
    """Only report unmet conditions. This version CAN NEVER claim REAL_DATA_READY."""
    require(isinstance(template,dict),"invalid candidate manifest")
    require(template.get("observed_odd_SEALED") is True,
            "OBSERVED_ODD_WOULD_BE_UNSEALED")
    require(template.get("any_real_ASDF_or_halo_catalogue_or_merger_tree_downloaded_in_this_stage")
            is False,"ACTUAL_ABACUS_DATA_ACCESSED_IN_NO_READ_STAGE")
    require(template.get("full_physical_finiteK_galaxy_bispectrum")=="BLOCKED",
            "FULL_PHYSICAL_GALAXY_B_CANNOT_BE_DECLARED")
    require(template.get("provider_checksum_method") not in ("zlib_crc32","ZIP_CRC32"),
            "ZLIB_CRC32_IS_NOT_GNU_POSIX_CKSUM")
    require(template.get("source_mass_definition")!="M200c_FROM_COMPASO_L1_ONLY",
            "COMPASO_L1_N_OR_SO_RADIUS_CANNOT_BE_RELABELED_M200C")
    reasons=[]
    def check_flag(key,message):
        if template.get(key) is not True:reasons.append(message)
    def check_value(key,ok,message):
        if not ok:reasons.append(message)
    check_flag("real_data_access_authorized","NO_SEPARATE_APPROVED_REAL_DATA_ACCESS_STAGE")
    check_flag("phase_and_box_independently_verified","C000_BOX_PHASE_NOT_VERIFIED")
    check_value("release_selection",
        template.get("release_selection") in ("RAW","CLEANED"),
        "RAW_OR_CLEANED_RELEASE_NOT_SELECTED")
    check_value("exact_provider_item_or_globus_path",
        bool(template.get("exact_provider_item_or_globus_path")),
        "EXACT_PROVIDER_FILE_AND_TREE_ITEM_PATHS_UNVERIFIED")
    check_value("selected_products_inventory",bool(template.get("selected_products_inventory")),
        "SMALL_EXACT_SUPERSLAB_AND_TREE_PRODUCTS_NOT_INVENTORIED")
    check_value("selected_transfer_bytes_verified",
        isinstance(template.get("selected_transfer_bytes_verified"),int)
        and not isinstance(template.get("selected_transfer_bytes_verified"),bool)
        and template["selected_transfer_bytes_verified"]>0,
        "SELECTED_FILE_BYTE_SIZE_NOT_VERIFIED_NO_BULK_TRANSFER")
    check_value("local_storage_bytes_verified",
        isinstance(template.get("local_storage_bytes_verified"),int)
        and not isinstance(template.get("local_storage_bytes_verified"),bool)
        and template["local_storage_bytes_verified"]>0,
        "LOCAL_STORAGE_CAPACITY_NOT_VERIFIED")
    check_flag("actual_halo_ASDF_header_read","ACTUAL_ASDF_HEADER_NOT_READ")
    check_value("actual_ASDF_redshift",
        isinstance(template.get("actual_ASDF_redshift"),(int,float))
        and not isinstance(template.get("actual_ASDF_redshift"),bool)
        and template["actual_ASDF_redshift"]>=0,
        "ACTUAL_ASDF_Z_UNVERIFIED_DIRECTORY_Z0P950_NOT_ACCEPTED")
    check_flag("provider_POSIX_cksum_verified","OFFICIAL_GNU_POSIX_CKSUM_NOT_VERIFIED")
    check_flag("local_file_sha256_verified","LOCAL_REAL_FILE_SHA256_NOT_VERIFIED")
    check_flag("true_particle_mass_and_units_from_header_verified",
               "PARTICLE_MASS_AND_UNITS_FROM_TRUE_HEADER_UNVERIFIED")
    check_flag("actual_SODensityL1_reference_species_and_units_verified",
               "ACTUAL_SO_DENSITY_REFERENCE_SPECIES_AND_UNITS_UNVERIFIED")
    check_flag("exact_real_M200c_200critical_provider_or_full_particle_remeasurement_verified",
               "REAL_200CRITICAL_M200C_NOT_VALIDATED")
    check_flag("tree_product_actual_availability_and_branch_verified",
               "EXACT_TREE_PRODUCT_AND_PROGENITOR_BRANCH_UNVERIFIED")
    check_flag("raw_cleaned_tree_mass_definition_and_epochs_consistent",
               "RAW_CLEANED_MASS_CONVENTION_OR_TRUE_EPOCHS_UNVERIFIED")
    check_flag("original_CLASS_As_tau_and_Fplus_Fminus_neutrino_response_consistent",
               "ORIGINAL_CLASS_AS_TAU_FPLUS_FMINUS_NEUTRINO_WAKE_NOT_MATCHED")
    check_flag("actual_LRG_ELG_HOD_selection_verified","ACTUAL_LRG_ELG_HOD_AND_SELECTION_UNVERIFIED")
    # Critical: no actual-file reader or provider file fetch exists in this stage,
    # and a user-edited JSON boolean is NOT a file integrity or physics certificate.
    reasons.append("NO_ACTUAL_ASDF_BYTES_HEADER_POSIX_CKSUM_M200C_VERIFIER_IN_THIS_STAGE")
    reasons.append("NO_ACTUAL_PROVIDER_MERGER_TREE_AND_LRG_ELG_HALO_SOURCE_ACCEPTED")
    return {
        "decision":"BLOCKED_DOCUMENTATION_ONLY_NO_REAL_ABACUS_DATA",
        "real_ABACUS_data_ready":False,
        "actual_halo_M200c_branch_accepted":False,
        "physical_finiteK_galaxy_bispectrum":"BLOCKED",
        "observed_odd_SEALED":True,
        "blockers":reasons,
        "blocker_count":len(reasons)
    }

def internal_tests(p,template):
    baseline=assess(template)
    require(baseline["real_ABACUS_data_ready"] is False
            and baseline["blocker_count"]>=18
            and "NO_ACTUAL_ASDF_BYTES_HEADER_POSIX_CKSUM_M200C_VERIFIER_IN_THIS_STAGE"
                in baseline["blockers"],
            "empty candidate could be treated as real-data ready")
    # Deliberately false z value entered from directory, no real header bytes:
    fake=copy.deepcopy(template)
    fake["actual_ASDF_redshift"]=.95
    res=assess(fake)
    require("ACTUAL_ASDF_HEADER_NOT_READ" in res["blockers"]
            and res["real_ABACUS_data_ready"] is False,
            "directory-only secondary redshift incorrectly accepted")
    # A naive textual statement of a perfectly full positive manifest must
    # NOT substitute for an independently verified provider ASDF file.
    forced=copy.deepcopy(template)
    for k,v in list(forced.items()):
        if v is False and k not in (
           "any_real_ASDF_or_halo_catalogue_or_merger_tree_downloaded_in_this_stage",):
            forced[k]=True
    forced["real_data_access_authorized"]=True
    forced["release_selection"]="RAW"
    forced["exact_provider_item_or_globus_path"]="FAKE://unverified"
    forced["selected_products_inventory"]=["FAKE_UNVERIFIED_FILE"]
    forced["selected_transfer_bytes_verified"]=10
    forced["local_storage_bytes_verified"]=10
    forced["actual_ASDF_redshift"]=.95
    require(assess(forced)["real_ABACUS_data_ready"] is False
            and "NO_ACTUAL_ASDF_BYTES_HEADER_POSIX_CKSUM_M200C_VERIFIER_IN_THIS_STAGE"
                in assess(forced)["blockers"],
            "all user-supplied self-reported flags cannot constitute file validation")
    for patch,expected in [
        ({"observed_odd_SEALED":False},"OBSERVED_ODD_WOULD_BE_UNSEALED"),
        ({"any_real_ASDF_or_halo_catalogue_or_merger_tree_downloaded_in_this_stage":True},
         "ACTUAL_ABACUS_DATA_ACCESSED_IN_NO_READ_STAGE"),
        ({"full_physical_finiteK_galaxy_bispectrum":"DETECTED"},
         "FULL_PHYSICAL_GALAXY_B_CANNOT_BE_DECLARED"),
        ({"provider_checksum_method":"zlib_crc32"},"ZLIB_CRC32_IS_NOT_GNU_POSIX_CKSUM"),
        ({"source_mass_definition":"M200c_FROM_COMPASO_L1_ONLY"},
         "COMPASO_L1_N_OR_SO_RADIUS_CANNOT_BE_RELABELED_M200C")]:
        fake=copy.deepcopy(template);fake.update(patch)
        try:assess(fake)
        except ValueError as e:
            require(expected in str(e),"negative control failed at wrong gate: "+expected)
        else:raise AssertionError("unaccepted real data or physics claim passed: "+expected)
    print("E17D2B3B0_ALL_7_OFFLINE_NEGATIVE_CONTROLS_PASS",
          "BASELINE_BLOCKERS",baseline["blocker_count"],
          "REAL_ABACUS_AND_B_BLOCKED",flush=True)
    return baseline

def report(p,template):
    res=internal_tests(p,template)
    out={
      "date":"2026-09-28",
      "status":"E17D2B3B0_OFFLINE_PROVIDER_ACCESS_CONTRACT_PREREG_TEMPLATE_SHA_GATES_PASS_REAL_ASDF_UNREAD",
      "prospective_protocol_git_blob":PROTO_BLOB,
      "exact_preregistered_candidate_template_sha256":sha(TEMPLATE.read_bytes()),
      "immutable_E8_E16_E17A_E17D2b2_E17D2b3a_source_ids":p["immutable_science_parents"],
      "official_source_URLs":[x["url"] for x in p["primary_documented_official_sources"]],
      "provider_documented_raw_doi":"10.13139/OLCF/1811689",
      "provider_documented_cleaned_doi":"10.13139/OLCF/1828535",
      "candidate_name_UNVERIFIED":"AbacusSummit_base_c000_ph000",
      "potential_narrow_access":"One provider-verified halo_info_XXX.asdf superslab and selected fields only after new authorized separate data stage; cleaned provenance and matching tree product separately verified.",
      "accession_and_bytes_not_assumed":True,
      "no_actual_ASDF_z_or_product_tree_mass_measurement":True,
      "no_new_science_redshift_or_mock_or_odd_or_CLASS":True,
      "no_actual_ASDF_halo_or_tree_downloaded":True,
      "offline_decision":res,
      "all_negatives_passed":True,
      "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
      "observed_odd_SEALED":True,
      "main_untouched_PR_draft":True}
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--self-test",action="store_true")
    modes.add_argument("--report",action="store_true")
    modes.add_argument("--local-steps",action="store_true")
    parser.add_argument("--output-dir",type=Path,default=OUT)
    args=parser.parse_args()
    p,t=gate()
    if args.local_steps:
        print("E17D2B3B0_WSL_NO_SCIENCE_DATA_DOWNLOAD",flush=True)
        print("cd ~/stress-energy-closure",flush=True)
        print("git fetch origin "+B,flush=True)
        print("git switch "+B,flush=True)
        print("git pull --ff-only origin "+B,flush=True)
        print("python -u scripts/audit_eboss_dr16_a03_e17d2b3b0_access_contract.py --self-test",flush=True)
        print("NO_GLOBUS_NO_ASDF_NO_MOCK_NO_CLASS_IN_THIS_COMMAND",flush=True)
        return
    result=report(p,t)
    if args.self_test:
        print("E17D2B3B0_OFFLINE_SHA_PINNED_PROVIDER_ACCESS_SELF_TEST_OK",
              "STATUS",result["offline_decision"]["decision"],
              "BLOCKERS",result["offline_decision"]["blocker_count"],
              "OBSERVED_ODD_SEALED_TRUE",flush=True)
    else:
        write_once(args.output_dir/REPORT,result)
        for i,message in enumerate(result["offline_decision"]["blockers"],1):
            print("E17D2B3B0_REAL_DATA_GATE_BLOCKER",i,message,flush=True)
        print("E17D2B3B0_DRY_RUN_COMPLETED_NO_REAL_ABACUS_NO_LOCAL_DOWNLOAD",
              "HALO_M200C_AND_B_BLOCKED",flush=True)
if __name__=="__main__":main()
