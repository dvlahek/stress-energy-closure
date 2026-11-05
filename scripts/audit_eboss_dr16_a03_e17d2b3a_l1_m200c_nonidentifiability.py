#!/usr/bin/env python3
"""E17D2b3a: exact positive-profile nonidentifiability of M200c from L1 mass/R.

Only synthetic dimensionless mathematical models, the public AbacusSummit
documentation and SHA-locked ORIGINAL E8/E16/E17A/E17D2b1/E17D2b2 reports.
NO Abacus ASDF, halo catalog/tree/particles, mock, CLASS, odd or galaxy input.
Neither illustrative Delta_L1=200 nor z=.95 is an observed ASDF header.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3a_mass_reference_nonidentifiability_prereg_2026-09-28.json"
PRE_BLOB="81554b0f7cdefa599f302d828ab3529812bb3a59"
B2DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28"
B2=B2DIR/"e17d2b2_original_abacus_c000_documented_metadata_bridge.json"
B2I=B2DIR/"e17d2b2_independent_ast_decimal_c000_metadata_replay.json"
B2MAN=B2DIR/"archive_manifest.json"
B2PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b2_abacussummit_c000_metadata_bridge_prereg_2026-09-28.json"
B2RUN=ROOT/"scripts/audit_eboss_dr16_a03_e17d2b2_abacus_c000_metadata_bridge.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
B1=ROOT/"source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28/e17d2b1_original_joint_conditional_exterior_assembly_rates.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3a_l1_vs_200c"
H=.6736;WB=.02237;WC=.1200;WNU=.00064420;Z=.95;DELTA_L1=200.
P=(1,2)

def require(ok,msg):
    if not ok:raise ValueError("E17D2B3A_FAIL_CLOSED: "+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def write_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"original synthetic output would change")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3a_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B3A_ORIGINAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective source-only protocol Git blob")
    p=json.loads(PRE.read_bytes());pa=p["original_pinned"]
    for path,key in ((B2MAN,"E17D2b2_archive_manifest_blob"),
                     (B2PRE,"E17D2b2_protocol_blob"),
                     (B2RUN,"E17D2b2_original_runner_blob")):
        require(blob(path.read_bytes())==pa[key],"immutable source Git blob: "+key)
    for path,key in ((E8,"E8_csv_sha256"),(E16,"E16_full_sha256"),
                     (E17A,"E17A_joint_sha256"),(B1,"E17D2b1_joint_sha256"),
                     (B2,"E17D2b2_original_sha256"),
                     (B2I,"E17D2b2_independent_sha256")):
        require(sha(path.read_bytes())==pa[key],"immutable source SHA: "+key)
    old=json.loads(B2.read_bytes())
    geometry=json.loads(E16.read_bytes())
    require(old["real_ASDF_header_or_halo_or_merger_tree_loaded"] is False
            and old["CompaSO_L1_SO_mean_epoch_header_NOT_M200c"] is True
            and old["c000_secondary_directory_z095_NOT_verified_ASDF_header_z"] is True
            and old["Abacus_c000_only_smooth_neutrino_Nbody_treatment"] is True
            and old["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and old["observed_odd_SEALED"] is True
            and geometry["QA"]["original_72_cases"]==72
            and geometry["QA"]["original_24_contrasts"]==24,
            "E17D2b2 documented original data/physics STOP")
    s=p["synthetic_scope"]
    require(s["h_c000_documented"]==H and s["omega_b_documented"]==WB
            and s["omega_cdm_documented"]==WC
            and s["omega_ncdm_documented"]==WNU
            and s["z_illustrative_ONLY_not_real_header"]==Z
            and s["delta_L1_mean_illustrative_ONLY_not_real_header"]==DELTA_L1
            and s["reference_radius_R_L1_normalized"]==1
            and s["reference_mass_M_L1_normalized"]==1
            and s["positive_enclosed_mass_exponents"]==list(P)
            and all(p["absolute_STOP"].values()),
            "locked synthetic shape/cosmology/no real data scope changed")
    return p

def exact_background():
    omega_m0=(WB+WC+WNU)/H**2
    e2=omega_m0*(1.+Z)**3+1.-omega_m0
    omega_m_z=omega_m0*(1.+Z)**3/e2
    delta_200c_mean=200./omega_m_z
    require(0.<omega_m0<1. and 0.<omega_m_z<1.
            and delta_200c_mean>DELTA_L1,
            "synthetic flat matter-Lambda positivity / two-root interior domain")
    return {"omega_m0_flat_matterLambda_TOY":omega_m0,
            "E2_z_flat_matterLambda_TOY":e2,
            "omega_m_z_flat_matterLambda_TOY":omega_m_z,
            "synthetic_200crit_over_total_matter_mean_threshold":delta_200c_mean,
            "synthetic_L1_SO_threshold_over_total_matter_mean_NOT_ABACUS_HEADER":DELTA_L1,
            "synthetic_ratio_threshold_200c_to_L1":delta_200c_mean/DELTA_L1}

def synthetic_model(power,b):
    require(power in P,"non-preregistered radial exponent")
    target=b["synthetic_200crit_over_total_matter_mean_threshold"]
    x=(DELTA_L1/target)**(1./(3.-power))
    m=x**power
    # Mean enclosed density is Delta_L1*(r/R_L1)^(power-3);
    # physical reference normalization uses rho_mean chosen from M_L1/R_L1^3.
    reclose=DELTA_L1*m/x**3/target
    require(0.<x<1. and 0.<m<1. and abs(reclose-1.)<1e-12,
            "synthetic exact 200critical SO root/mass reclosure")
    def radial_mass(radius):
        require(radius>=0.,"radius must be nonnegative")
        return min(radius,1.)**power
    qa_radii=(0.,x/2.,x,1.)
    masses=[radial_mass(r) for r in qa_radii]
    require(masses[0]==0. and masses[-1]==1.
            and all(masses[i]<masses[i+1] for i in range(3)),
            "positive integrable radial profile and monotone enclosed mass")
    return {"enclosed_mass_exponent_synthetic_ONLY":power,
            "normalization_shared_L1_mass":1.,
            "normalization_shared_L1_radius":1.,
            "shared_exact_L1_SO_delta_mean_idealized":DELTA_L1,
            "mass_profile":"M_p(<r)/M_L1=(r/R_L1)^p for 0<=r<=R_L1, constant outside",
            "positive_density_at_r_gt_zero":
                "rho_p(r)=p*M_L1*r^(p-3)/(4*pi*R_L1^p)",
            "synthetic_QA_radii_normalized":list(qa_radii),
            "positive_monotone_enclosed_mass_at_QA_radii":masses,
            "r200c_over_shared_R_L1_TOY":x,
            "M200c_over_shared_M_L1_TOY":m,
            "200critical_mean_density_reclosure_relative":abs(reclose-1.)}

def validate_actual_header(record):
    """Fail closed for synthetic negative controls; this stage never reads ASDF."""
    require(record.get("provenance")=="REAL_PROVIDER_VERIFIED_ASDF_FILE_WITH_LOCAL_SHA_AND_POSIX_CKSUM",
            "REAL_ASDF_PROVENANCE_MISSING")
    require("actual_header_redshift" in record
            and isinstance(record["actual_header_redshift"],(int,float))
            and not isinstance(record["actual_header_redshift"],bool)
            and math.isfinite(record["actual_header_redshift"]),
            "REAL_ASDF_HEADER_Z_MISSING_DIRECTORY_Z_NOT_ACCEPTED")
    require(record.get("mass_definition")=="INDEPENDENTLY_VALIDATED_200_CRITICAL_MASS"
            and record.get("reference_density")=="200_CRITICAL_AT_ACTUAL_HEADER_Z"
            and record.get("radial_remeasurement_or_provider_validation"),
            "COMPASO_L1_MASS_OR_SO_RADIUS_NOT_M200C")
    require(record.get("particle_positions_or_independent_M200c_provenance_verified") is True
            and record.get("tree_branch_and_cleaning_verified") is True,
            "EXACT_200CRITICAL_MASS_HISTORY_AND_TREE_NOT_VERIFIED")
    return True

def original():
    p=gate();b=exact_background()
    rows=[synthetic_model(k,b) for k in P]
    difference=abs(rows[0]["M200c_over_shared_M_L1_TOY"]-
                   rows[1]["M200c_over_shared_M_L1_TOY"])
    difference_radius=abs(rows[0]["r200c_over_shared_R_L1_TOY"]-
                          rows[1]["r200c_over_shared_R_L1_TOY"])
    require(difference>.05 and difference_radius>0.,
            "distinct non-identical M200c from same idealized L1 M,R")
    # These are controlled negative cases, never real provider metadata.
    fake={"provenance":"SYNTHETIC_HEADER_TEST_ONLY",
          "directory_label_redshift":Z,
          "mass_definition":"COMPASO_L1_N_TIMES_PARTICLE_MASS",
          "reference_density":"L1_SO_RELATIVE_MEAN"}
    negative=[]
    try:validate_actual_header(fake)
    except ValueError as e:
        require("REAL_ASDF_PROVENANCE_MISSING" in str(e),
                "fake header was rejected at incorrect gate")
        negative.append("synthetic ASDF not accepted as actual provider file")
    else:raise AssertionError("FAKE_ASDF_WAS_ACCEPTED")
    fake["provenance"]="REAL_PROVIDER_VERIFIED_ASDF_FILE_WITH_LOCAL_SHA_AND_POSIX_CKSUM"
    try:validate_actual_header(fake)
    except ValueError as e:
        require("REAL_ASDF_HEADER_Z_MISSING" in str(e),
                "directory redshift incorrectly accepted as actual header z")
        negative.append("directory-only z0.95 rejected")
    else:raise AssertionError("DIRECTORY_Z_WAS_ACCEPTED")
    fake["actual_header_redshift"]=Z
    try:validate_actual_header(fake)
    except ValueError as e:
        require("COMPASO_L1_MASS_OR_SO_RADIUS_NOT_M200C" in str(e),
                "L1 mass illegally relabelled M200c")
        negative.append("CompaSO L1 N or SO_radius rejected as 200crit")
    else:raise AssertionError("L1_MASS_AS_M200C_ACCEPTED")
    return {"date":"2026-09-28",
      "status":"E17D2B3A_IDEALIZED_TWO_POSITIVE_SAME_L1_MASS_RADIUS_PROFILES_DIFFERENT_M200C_CERTIFIED_REAL_ABACUS_DATA_UNREAD",
      "prospective_protocol_git_blob":PRE_BLOB,
      "immutable_original_parent_IDs":p["original_pinned"],
      "official_public_documentation_URLs":[s["url"] for s in p["official_sources"]],
      "source_only_analytic_flat_matterLambda_TOY":b,
      "two_predeclared_positive_normalized_radial_profiles":rows,
      "same_L1_mass_radius_density_and_synthetic_z":True,
      "distinct_M200c_mass_fraction_absolute_gap_synthetic_only":difference,
      "distinct_r200c_over_L1_radius_gap_synthetic_only":difference_radius,
      "synthetic_header_directory_z_and_L1_mass_negative_controls":negative,
      "idealized_SO_crossing_is_stronger_than_real_SO_radius_metadata":True,
      "actual_Abacus_SODensityL1_header_read":False,
      "actual_Abacus_ASDF_redshift_header_read":False,
      "actual_Abacus_halo_particle_or_merger_tree_read":False,
      "any_real_halo_M200c_remeasured":False,
      "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
      "observed_odd_SEALED":True,
      "new_CLASS_ASDF_FITS_halo_catalogue_mock_download":False,
      "main_untouched_PR_draft":True}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    a=ap.parse_args()
    r=original()
    if a.self_test:
        print("E17D2B3A_ORIGINAL_SHA_AND_TOY_POSITIVE_DENSITY_PREFLIGHT_PASS",
              "L1_MASS_SAME_BUT_M200C_GAP",
              r["distinct_M200c_mass_fraction_absolute_gap_synthetic_only"],
              "REAL_ASDF_NOT_READ",flush=True)
    else:
        write_once(a.output_dir/"e17d2b3a_original_l1_200critical_nonidentifiability.json",r)
        print("E17D2B3A_ORIGINAL_SAME_L1_DIFFERENT_200CRITICAL_MASS_PASS",
              "MASS_GAP",r["distinct_M200c_mass_fraction_absolute_gap_synthetic_only"],
              "REAL_M200C_NOT_MEASURED_GALAXY_B_BLOCKED",flush=True)
if __name__=="__main__":main()
