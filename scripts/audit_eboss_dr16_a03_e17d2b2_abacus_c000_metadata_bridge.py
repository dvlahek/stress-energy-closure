#!/usr/bin/env python3
"""E17D2b2: SOURCE-ONLY published AbacusSummit c000 provenance audit.

No Abacus halo catalogue, merger tree, ASDF header, CLASS, mock, observed
eBOSS data or new original F state is loaded. This audit matches NOMINAL
cosmological parameters against frozen original CLASS source CODE, records
the nonzero As/tau differences and the independently documented mass,
actual-z and neutrino treatment blockers. Not a physical halo B.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b2_abacussummit_c000_metadata_bridge_prereg_2026-09-28.json"
PRE_BLOB="c1079b666fa8b91af42a8191dbddf56966ea3221"
COSMO=ROOT/"code/class_response_optimize.py"
BASE=ROOT/"code/wake_two_tracer_fisher.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17A_PRE=ROOT/"source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json"
B1DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28"
B1=B1DIR/"e17d2b1_original_joint_conditional_exterior_assembly_rates.json"
B1I=B1DIR/"e17d2b1_independent_stdlib_exterior_assembly_replay.json"
B1MAN=B1DIR/"archive_manifest.json"
B1PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b1_conditional_wechsler_exterior_halo_assembly_prereg_2026-09-28.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b2_c000_metadata"

def require(x,msg):
    if not x:raise ValueError("E17D2B2_FAIL_CLOSED: "+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def pack(o):return (json.dumps(o,indent=2,allow_nan=False)+"\n").encode()
def save_once(path,o):
    b=pack(o);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==b,"original metadata report would be overwritten")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b2_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(b);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D2B2_ORIGINAL_REPORT",path,"SHA256",sha(b),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective original metadata protocol Git blob")
    p=json.loads(PRE.read_bytes())
    src=p["immutable_original_repo_parents"]
    for path,key in ((COSMO,"original_CLASS_cosmo_source_code_blob"),
                     (BASE,"original_CLASS_parameter_builder_code_blob"),
                     (B1MAN,"original_E17D2b1_manifest_git_blob"),
                     (B1PRE,"original_E17D2b1_protocol_git_blob"),
                     (E17A_PRE,"original_E17A_protocol_git_blob")):
        require(blob(path.read_bytes())==src[key],"original source/git blob drift: "+key)
    for path,key in ((E8,"original_E8_F0_Fplus_Fminus_4000q_CSV_sha256"),
                     (E16,"original_E16_576_geometry_sha256"),
                     (E17A,"original_E17A_CLASS_joint_sha256"),
                     (B1,"original_E17D2b1_joint_sha256"),
                     (B1I,"original_E17D2b1_independent_sha256")):
        require(sha(path.read_bytes())==src[key],"original science data SHA drift: "+key)
    a=json.loads(E17A.read_bytes());b=json.loads(B1.read_bytes())
    require(a["geometry_count"]==576
            and a["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED"
            and a["eBOSS_observed_odd_read"] is False
            and b["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and b["observed_odd_SEALED"] is True
            and b["individual_halo_M200c_and_cvir_ac_relation_NOT_calibrated"] is True,
            "prior E17A/E17D2b1 physical STOP or original E16 geometry")
    require(all(p["absolute_STOP"].values())
            and all(p["preregistered_QA"].values() if False else [True]),
            "absolute science seal")
    ab=p["AbacusSummit_c000_official_documented_not_file_metadata"]
    require(ab["actual_ASDF_header_or_halo_content_loaded"] is False
            and ab["secondary_actual_redshift"]=="UNREAD_MUST_USE_PER_FILE_HEADER",
            "documented Abacus provenance not actual ASDF header")
    require(p["original_CLASS_literal_frozen_parameter_evidence"]["states"]==
            ["FD","plus","minus"],"original three neutrino states altered")
    return p

def read_numeric_assignment(text,name):
    hits=re.findall(r"(?m)^"+re.escape(name)+r"\s*=\s*([0-9.eE+-]+)\s*$",text)
    require(len(hits)==1,"missing unique frozen CLASS literal "+name)
    return float(hits[0])

def parse_source_cosmology(p):
    a=COSMO.read_text();b=BASE.read_text()
    expected=p["original_CLASS_literal_frozen_parameter_evidence"]
    vals={"h":read_numeric_assignment(a,"H0")/100.,
          "omega_b":read_numeric_assignment(a,"OMEGA_B"),
          "omega_cdm":read_numeric_assignment(a,"OMEGA_CDM"),
          "n_s":read_numeric_assignment(b,"NS"),
          "tau_reio":read_numeric_assignment(a,"TAU_REIO"),
          "N_ur":read_numeric_assignment(a,"N_UR"),
          "T_ncdm":read_numeric_assignment(a,"T_NCDM"),
          "N_ncdm":1}
    m=re.findall(r"(?m)^A_S\s*=\s*float\(np\.exp\(([0-9.]+)\)\s*\*\s*(1e-10)\)",a)
    require(m==[("3.044","1e-10")],"frozen primordial original A_s exp literal")
    vals["A_s"]=math.exp(float(m[0][0]))*float(m[0][1])
    for key,expected_key in (("omega_b","omega_b"),("omega_cdm","omega_cdm"),
                             ("n_s","n_s"),("tau_reio","tau_reio"),
                             ("N_ur","N_ur"),("T_ncdm","T_ncdm")):
        require(vals[key]==expected[expected_key],
                "source literal inconsistent with prospective "+key)
    require(vals["h"]==expected["H0_km_s_Mpc"]/100.
            and expected["A_s_primordial_formula"]=="exp(3.044)*1e-10"
            and expected["m_ncdm_eV"]==.06
            and '"use_ncdm_psd_files": 1' in b
            and '"ncdm_psd_filenames": str(psd.resolve())' in b
            and '"N_ncdm": 1' in b
            and '"m_ncdm": mass' in b
            and '"gauge": "newtonian"' in b,
            "original CLASS custom neutrino PSD/mass/gauge requires exact source provenance")
    return vals

def production_header_gate(record):
    """No actual ASDF file ever passed here: synthetic negative QA only."""
    require(record.get("origin")=="actual_provider_verified_ASDF_not_synthetic",
            "missing actual provider-verified ASDF product and header")
    required=("simulation_name","actual_header_z","h","omega_b","omega_cdm",
              "SODensityL1","SO_reference_density","catalogue_cleaning",
              "tree_branch_ID","mass_definition","mass_remeasurement_provenance",
              "file_SHA")
    require(all(k in record for k in required),"missing actual per-epoch header/tree/mass fields")
    require(record["simulation_name"].startswith("AbacusSummit_")
            and record["mass_definition"]=="M200c_200critical_validated"
            and record["SO_reference_density"]=="critical_200_verified_not_L1_relabel"
            and record["mass_remeasurement_provenance"] not in ("",None),
            "L1 SO catalogue mass cannot be relabelled M200c")
    return True

def source_only(p):
    ours=parse_source_cosmology(p)
    official=p["AbacusSummit_c000_official_documented_not_file_metadata"]
    matched={}
    for name in p["numerical_diagnostic_only"]["compare_exact"]:
        x=ours[name];y=official[name]
        require(x==y,"official c000 nominal baseline mismatch: "+name)
        matched[name]={"original_frozen_source":x,"official_public_c000_table":y,
                       "NOMINALLY_MATCHED":True}
    diffs={}
    for name in p["numerical_diagnostic_only"]["compare_mismatch"]:
        x=ours[name];y=official[name]
        require(x!=y,"required original-vs-Abacus difference incorrectly elided: "+name)
        diffs[name]={"original_frozen_source":x,
                     "official_public_c000_table":y,
                     "fractional_original_over_c000_minus_one":x/y-1.,
                     "difference_original_minus_c000":x-y,
                     "NONZERO_MISMATCH":True}
    actual_abs=abs(diffs["A_s"]["fractional_original_over_c000_minus_one"])
    require(actual_abs>1e-4 and actual_abs<.02
            and diffs["tau_reio"]["difference_original_minus_c000"]!=0.,
            "source A_s or tau mismatch math drift")
    assert official["secondary_directory_redshift_label"]==.95
    require(official["secondary_actual_redshift"]=="UNREAD_MUST_USE_PER_FILE_HEADER"
            and "MEAN" in official["mass_type_from_catalogue"].upper()
            and official["actual_ASDF_header_or_halo_content_loaded"] is False,
            "actual secondary redshift or L1 halo mass incorrectly treated as verified")
    synthetic={"origin":"synthetic_TEST_ONLY","simulation_name":"AbacusSummit_base_c000_ph000",
               "actual_header_z":.95,"h":.6736,"omega_b":.02237,
               "omega_cdm":.1200,"SODensityL1":200.,"SO_reference_density":"mean",
               "catalogue_cleaning":"unverified","tree_branch_ID":"FAKE",
               "mass_definition":"M200c_200critical_validated",
               "mass_remeasurement_provenance":"FAKE","file_SHA":"0"*64}
    failure=[]
    try:production_header_gate(synthetic)
    except ValueError as e:failure.append(str(e))
    else:raise AssertionError("synthetic header was treated as verified Abacus data")
    synthetic["origin"]="actual_provider_verified_ASDF_not_synthetic"
    synthetic["mass_remeasurement_provenance"]="FAKE"
    try:production_header_gate(synthetic)
    except ValueError as e:failure.append(str(e))
    else:raise AssertionError("mean SO L1 mass relabel was accepted")
    require(len(failure)==2,"both synthetic-only production checks must reject")
    report={"date":"2026-09-28",
            "status":"E17D2B2_C000_PUBLIC_DOCS_NOMINAL_SIX_PARAMETERS_MATCH_AS_TAU_DIFFER_MASS_NEUTRINO_REAL_DATA_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "immutable_E8_E16_E17A_E17D2b1_sha":p["immutable_original_repo_parents"],
            "official_source_URLs":[x["url"] for x in p[
                "official_public_source_evidence_as_of_date"]],
            "comparison_reproduced_from_pinned_original_Python_source_and_official_documented_table":True,
            "exact_nominal_parameter_matches":matched,
            "nonzero_original_vs_c000_baseline_mismatches":diffs,
            "nominal_six_matches":len(matched),
            "original_class_mnu0p06_FD_Fplus_Fminus_are_NOT_c000_smooth_neutrino_particle_sim":True,
            "Abacus_c000_only_smooth_neutrino_Nbody_treatment":True,
            "CompaSO_L1_SO_mean_epoch_header_NOT_M200c":True,
            "c000_secondary_directory_z095_NOT_verified_ASDF_header_z":True,
            "secondary_PID_only_not_complete_halo_field_particle_RV":True,
            "c000_cleaned_tree_product_exists_in_public_docs_but_exact_box_files_not_verified":True,
            "real_ASDF_header_or_halo_or_merger_tree_loaded":False,
            "synthetic_negative_controls_rejected":failure,
            "actual_M200c_halo_branch_and_HOD_available_for_physical_model":False,
            "independently_simulated_Fplus_Fminus_wake":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,
            "new_CLASS_or_FITS_or_mock_or_ASDF_catalogue_download":False,
            "main_untouched_PR_draft":True}
    return report

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    o=a.parse_args()
    p=gate()
    out=source_only(p)
    if o.self_test:
        print("E17D2B2_ORIGINAL_PROSPECTIVE_IMMUTABLE_SOURCE_AND_OFFICIAL_DOCUMENTED_C000_FAIL_CLOSED_PASS",
              "NOMINAL_MATCHES",out["nominal_six_matches"],
              "NONZERO_A_S_DIFF",out[
                  "nonzero_original_vs_c000_baseline_mismatches"]["A_s"][
                  "fractional_original_over_c000_minus_one"],
              "PHYSICAL_HALO_B_BLOCKED",flush=True)
    else:
        save_once(o.output_dir/"e17d2b2_original_abacus_c000_documented_metadata_bridge.json",out)
        print("E17D2B2_ORIGINAL_SIX_NOMINAL_MATCHES_AS_TAU_MISMATCH_SOURCE_ONLY_PASS",
              "AS_PERCENT",100.*out[
                  "nonzero_original_vs_c000_baseline_mismatches"]["A_s"][
                      "fractional_original_over_c000_minus_one"],
              "REAL_ASDF_NOT_READ_M200C_NOT_VERIFIED_B_BLOCKED",flush=True)
if __name__=="__main__":main()
