#!/usr/bin/env python3
"""E17D2b1: conditional Wechsler-type M200c-shape / exterior halo monopole.

The Wechsler functional form was fitted for M_vir; its use for the original
DM14 M200c snapshots is an EXPLICIT UNCALIBRATED test ansatz. Two chosen
growth parameters are NOT inferred halo formation times, LRG/ELG HOD or
simulation merger trees. Time is ordinary DM14 cosmic-time background, NOT
E17D1 arbitrary u and NOT original massive-neutrino CLASS. NO GALAXY B.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b1_conditional_wechsler_exterior_halo_assembly_prereg_2026-09-28.json"
PRE_BLOB="28e33fad09db24b8b19739f50322169180442988"
B0DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28"
B0_JOINT=B0DIR/"e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json"
B0_INDEPENDENT=B0DIR/"e17d2b0_independent_stdlib_dm14_full_radial_replay.json"
B0_MAN=B0DIR/"archive_manifest.json"
B0_PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b0_dm14_population_halo_mass_concentration_prereg_2026-09-28.json"
B0_RUNNER=ROOT/"scripts/audit_eboss_dr16_a03_e17d2b0_dm14_halo_snapshots.py"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
D1=ROOT/"source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_original_joint_causal_history_nonidentifiability.json"
D2A=ROOT/"source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28/e17d2a_original_joint_nfw_static_spatial_poisson_family.json"
Z0=.95
ZS=(.945,.95,.955)
ALPHA=(.4,.8)
MASS_HINV=(1000000000000,10000000000000)
ANCHOR_INDICES=(4,13)
HDM=.671;HCLASS=.6736
OM=.3175;OL=.6825
G=4.30091e-9
MPC_KM=3.0856775814913673e19
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b1_wechsler_exterior"

def require(ok,msg):
    if not ok:raise ValueError("E17D2B1_FAIL_CLOSED: "+msg)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def pack(obj):return (json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
def save_once(path,obj):
    raw=pack(obj);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"existing archive would change: "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b1_",delete=False) as fh:
            temp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17D2B1_ORIGINAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective Wechsler source/history Git blob")
    p=json.loads(PRE.read_bytes());parent=p["immutable_parents"]
    paths=((E8,"E8_original_4000q_CSV_sha256"),
           (E16,"E16_original_full_sha256"),
           (D1,"E17D1_original_joint_sha256"),
           (D2A,"E17D2a_original_joint_sha256"),
           (B0_JOINT,"E17D2b0_original_joint_sha256"),
           (B0_INDEPENDENT,"E17D2b0_independent_sha256"))
    for path,key in paths:
        require(sha(path.read_bytes())==parent[key],"frozen source SHA "+key)
    require(blob(B0_MAN.read_bytes())==parent["E17D2b0_archive_manifest_git_blob"]
            and blob(B0_PRE.read_bytes())==parent["E17D2b0_protocol_git_blob"]
            and blob(B0_RUNNER.read_bytes())==parent["E17D2b0_original_runner_git_blob"],
            "original B0 archive/protocol/runner Git blob drift")
    previous=json.loads(B0_JOINT.read_bytes())
    require(previous["number_of_frozen_conditional_halo_snapshots"]==18
            and previous["total_conditional_spatial_Fourier_modes"]==31104
            and previous["individual_LRG_ELG_halo_M_HOD_time_history_known"] is False
            and previous["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and previous["observed_odd_SEALED"] is True,
            "B0 source-only absolute STOP or original QA drift")
    config=p["frozen"]
    require(config["original_mass_reference_hDM_inverse_Msun"]==list(MASS_HINV)
            and config["original_DM14_anchor_case_indices"]==list(ANCHOR_INDICES)
            and config["z0_anchor_only"]==Z0
            and config["fixed_original_E17A_z_math_nodes_NOT_real_halo_formation_samples"]==list(ZS)
            and config["Wechsler_shape_alpha_equals_2_a_c_QA_only"]==list(ALPHA)
            and config["physical_exterior_test_radius_ratio_to_anchor_r200c"]==3
            and config["DM14_h"]==HDM
            and config["DM14_Omega_m"]==OM
            and config["DM14_Omega_lambda"]==OL
            and config["original_CLASS_h"]==HCLASS
            and config["G_Mpc_km2_s2_Msun"]==G
            and config["original_tracer_order"]==["LRG","ELG"]
            and all(p["absolute_STOP"].values()),
            "original geometry/scoped shape family or mandatory STOP drift")
    require(json.loads(E16.read_bytes())["QA"]["original_72_cases"]==72
            and json.loads(E16.read_bytes())["QA"]["original_24_contrasts"]==24,
            "original E14/E15 science source counts")
    return p,previous

def original_anchor(i,p,previous):
    require(i in range(2),"only prospectively anchored original mass cases")
    index=ANCHOR_INDICES[i]
    path=B0DIR/("e17d2b0_original_DM14_conditional_snapshot_%02d.json"%index)
    raw=path.read_bytes()
    key=("E17D2b0_anchor_mass_1e12_z095_original_SHA" if i==0 else
         "E17D2b0_anchor_mass_1e13_z095_original_SHA")
    require(sha(raw)==p["immutable_parents"][key],
            "exact immutable original DM14 median anchor SHA")
    require(sha(raw)==previous["all_original_case_SHA"][path.name]["sha256"],
            "original joint SHA does not match published median anchor")
    data=json.loads(raw)
    snap=data["published_population_halo_snapshot"]
    require(data["case_index"]==index
            and data["original_CLASS_and_DM14_h_equal"] is False
            and data["same_mass_across_z_is_NOT_one_halo_mass_accretion_history"] is True
            and snap["M200c_hDM_inverse_Msun_QA_anchor"]==MASS_HINV[i]
            and snap["z_original_CLASS_math_node_only"]==Z0
            and snap["illustrative_log10_concentration_offset_dex_NOT_z095_posterior"]==0.
            and data["observed_odd_SEALED"] is True,
            "original halo median anchor or physical provenance changed")
    return path,raw,snap

def H_at_z(z):
    return 100.*HDM*math.sqrt(OM*(1.+z)**3+OL)

def mass_at_z(M0,alpha,z):
    return M0*math.exp(-alpha*(z-Z0)/(1.+Z0))

def snapshot(M0,alpha,z,R):
    M=mass_at_z(M0,alpha,z)
    H=H_at_z(z)
    rho=3.*H*H/(8.*math.pi*G)
    radius=(3.*M/(4.*math.pi*200.*rho))**(1./3.)
    Phi=-G*M/R
    grad=G*M/(R*R)  # dPhi/dr positive; acceleration inward.
    D=-alpha/(1.+Z0)
    gamma=alpha*(1.+z)/(1.+Z0)
    growth_rate_s=gamma*(H/MPC_KM)
    mass_reclosure=(4.*math.pi/3.)*200.*rho*radius**3/M
    phi_reclosure=(-Phi*R)/(G*M)
    force_reclosure=grad*R*R/(G*M)
    return {"z_frozen_original_math_node":z,
            "mass_M200c_physical_Msun_conditional_ansatz":M,
            "mass_ratio_to_same_frozen_anchor":M/M0,
            "DM14_H_km_s_Mpc":H,
            "DM14_rho_crit_Msun_Mpc3":rho,
            "r200c_phys_Mpc_conditional_ansatz":radius,
            "fixed_phys_exterior_R_Mpc_QA_only":R,
            "exterior_margin_R_over_r200c":R/radius,
            "exterior_Phi_phys_km2_s2":Phi,
            "exterior_inward_acceleration_magnitude_km2_s2_per_Mpc":grad,
            "exact_dlnM_dz_and_dlnAbsPhi_dz":D,
            "exact_dlnM_dln_a_and_dlnAbsPhi_dln_a":gamma,
            "exact_dlnM_dt_seconds_inverse_DM14_only":growth_rate_s,
            "exact_dPhi_dz_at_fixed_physical_R_km2_s2":Phi*D,
            "mass_radius_reclosure_relative":abs(mass_reclosure-1.),
            "exterior_Phi_mass_reclosure_relative":abs(phi_reclosure-1.),
            "exterior_force_mass_reclosure_relative":abs(force_reclosure-1.)}

def filename(index):
    return "e17d2b1_original_conditional_wechsler_exterior_%02d.json"%index

def run_case(index,out):
    p,previous=gate()
    require(index in range(4),"only two original masses and two prereg rate anchors")
    i=index//2;alpha=ALPHA[index%2]
    anchor_path,anchor_raw,anchor=original_anchor(i,p,previous)
    M0=anchor["M200c_physical_Msun_DM14"]
    radius0=anchor["r200c_physical_Mpc_DM14"]
    R=3.*radius0
    snaps=[snapshot(M0,alpha,z,R) for z in ZS]
    q=p["prospectively_locked_QA"]
    require(len(snaps)==3
            and snaps[0]["mass_M200c_physical_Msun_conditional_ansatz"]>M0
            and snaps[2]["mass_M200c_physical_Msun_conditional_ansatz"]<M0,
            "monotone positive mass accretion towards later a")
    require(abs(snaps[1]["mass_M200c_physical_Msun_conditional_ansatz"]/M0-1.)<
            q["exact_anchor_mass_and_reference_r200c_relative_max"]
            and abs(snaps[1]["r200c_phys_Mpc_conditional_ansatz"]/radius0-1.)<
            q["exact_anchor_mass_and_reference_r200c_relative_max"],
            "cannot reclose exact immutable DM14 anchor at z0")
    for obj in snaps:
        require(obj["mass_M200c_physical_Msun_conditional_ansatz"]>0
                and obj["r200c_phys_Mpc_conditional_ansatz"]>0
                and obj["exterior_margin_R_over_r200c"]>1.
                and obj["exterior_Phi_phys_km2_s2"]<0
                and obj["exterior_inward_acceleration_magnitude_km2_s2_per_Mpc"]>0,
                "conditional exterior NFW monopole positivity/exteriority")
        for key in ("mass_radius_reclosure_relative",
                    "exterior_Phi_mass_reclosure_relative",
                    "exterior_force_mass_reclosure_relative"):
            require(obj[key]<q["mass_radius_200critical_reclosure_relative_max"],
                    "exterior shell theorem or 200-critical mass-radial closure: "+key)
    dZ=ZS[2]-ZS[0]
    analytical=-alpha/(1.+Z0)
    fd_log=(math.log(snaps[2]["mass_M200c_physical_Msun_conditional_ansatz"])-
            math.log(snaps[0]["mass_M200c_physical_Msun_conditional_ansatz"]))/dZ
    fd_mass=(snaps[2]["mass_M200c_physical_Msun_conditional_ansatz"]-
             snaps[0]["mass_M200c_physical_Msun_conditional_ansatz"])/dZ
    exact_mass=analytical*M0
    err_fd_mass=abs(fd_mass-exact_mass)/abs(exact_mass)
    err_fd_log=abs(fd_log-analytical)
    require(err_fd_log<q["exact_analytic_log_derivative_vs_two_sided_original_z_fd_abs_max"]
            and err_fd_mass<q["exact_analytic_mass_derivative_vs_two_sided_original_z_fd_relative_max"],
            "frozen original z-grid derivative QA of analytic halo mass history")
    outj={"date":"2026-09-28",
          "status":"E17D2B1_ORIGINAL_WECHSLER_TYPE_CONDITIONAL_EXTERIOR_MONOPOLE_HISTORY_QA_ONLY_FULL_PHYSICAL_B_BLOCKED",
          "prospective_protocol_git_blob":PRE_BLOB,
          "original_source_mass_reference_hDM_inverse_Msun":MASS_HINV[i],
          "original_E17D2b0_anchor_index":ANCHOR_INDICES[i],
          "original_E17D2b0_anchor_SHA256":sha(anchor_raw),
          "original_DM14_median_c200c_at_z0_ONLY_not_evolved":anchor[
              "c200c_population_median_DM14"],
          "original_DM14_anchor_r200c_phys_Mpc":radius0,
          "original_DM14_anchor_mass_physical_Msun":M0,
          "Wechsler_form_alpha_equals_two_ac_math_QA_NOT_fitted_to_halo":alpha,
          "redshift_anchor_z0":Z0,
          "fixed_physical_exterior_R_equals_three_times_anchor_r200c_Mpc":R,
          "snapshots_sorted_z_increasing":snaps,
          "QA":{"max_mass_radius_and_exterior_reclosure_relative":max(
                    max(o["mass_radius_reclosure_relative"],
                        o["exterior_Phi_mass_reclosure_relative"],
                        o["exterior_force_mass_reclosure_relative"]) for o in snaps),
                "fd_3node_logmass_derivative_abs_gap":err_fd_log,
                "fd_3node_mass_derivative_relative_gap":err_fd_mass,
                "min_exterior_radius_margin_R_over_r200c":min(
                    o["exterior_margin_R_over_r200c"] for o in snaps)},
          "same_DM14_median_c200c_at_z0_but_no_cvir_conversion_or_evolved_c":True,
          "halo_M200c_use_of_Wechsler_virial_M_history_UNCALIBRATED_ansatz":True,
          "physical_CLASS_plus_DM14_cosmology_consistent":False,
          "actual_LRG_ELG_halo_assembly_HOD_known":False,
          "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
          "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
          "main_untouched_PR_draft":True}
    save_once(out/filename(index),outj)
    print("E17D2B1_ORIGINAL_CASE",index,"M_hDM",MASS_HINV[i],
          "ALPHA_SHAPE_ONLY",alpha,
          "ANCHOR_C",outj["original_DM14_median_c200c_at_z0_ONLY_not_evolved"],
          "M_GROWTH_RATIO_OLD_VS_NEW",
          snaps[2]["mass_ratio_to_same_frozen_anchor"],
          snaps[0]["mass_ratio_to_same_frozen_anchor"],
          "FD_MASS_REL",err_fd_mass,
          "NO_TRUE_HALO_HISTORY_NO_B",flush=True)

def aggregate(out):
    p,previous=gate()
    reports={}
    for i in range(4):
        raw=(out/filename(i)).read_bytes();r=json.loads(raw)
        require(r["original_source_mass_reference_hDM_inverse_Msun"]==MASS_HINV[i//2]
                and r["Wechsler_form_alpha_equals_two_ac_math_QA_NOT_fitted_to_halo"]==
                ALPHA[i%2] and len(r["snapshots_sorted_z_increasing"])==3
                and r["halo_M200c_use_of_Wechsler_virial_M_history_UNCALIBRATED_ansatz"] is True
                and r["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
                and r["observed_odd_SEALED"] is True,
                "original 4 conditional histories or mandatory stop")
        reports[filename(i)]={"sha256":sha(raw),"bytes":len(raw),
                              "mass_anchor_hDM_Msun":MASS_HINV[i//2],
                              "alpha_shape_only":ALPHA[i%2],
                              "QA":r["QA"]}
    differences={}
    for i in range(2):
        x=json.loads((out/filename(i*2)).read_bytes())
        y=json.loads((out/filename(i*2+1)).read_bytes())
        dx=x["snapshots_sorted_z_increasing"][1]
        dy=y["snapshots_sorted_z_increasing"][1]
        gapmass=abs(dx["mass_M200c_physical_Msun_conditional_ansatz"]-
                    dy["mass_M200c_physical_Msun_conditional_ansatz"])/dx[
                        "mass_M200c_physical_Msun_conditional_ansatz"]
        gapphi=abs(dx["exterior_Phi_phys_km2_s2"]-
                   dy["exterior_Phi_phys_km2_s2"])/abs(dx[
                       "exterior_Phi_phys_km2_s2"])
        gapc=abs(x["original_DM14_median_c200c_at_z0_ONLY_not_evolved"]-
                 y["original_DM14_median_c200c_at_z0_ONLY_not_evolved"])
        gaprate=abs(dx["exact_dlnM_dln_a_and_dlnAbsPhi_dln_a"]-
                    dy["exact_dlnM_dln_a_and_dlnAbsPhi_dln_a"])
        q=p["prospectively_locked_QA"]
        require(gapmass<q["distinct_alpha_fixed_anchor_at_z095_zero_mass_gap_max"]
                and gapphi<q["distinct_alpha_fixed_anchor_at_z095_zero_mass_gap_max"]
                and gapc<q["distinct_alpha_fixed_anchor_at_z095_zero_mass_gap_max"]
                and gaprate>q["distinct_alpha_growth_rate_gap_at_z095_min"],
                "two rate histories do not share exact frozen mass/c/Phi anchor or derivative witness")
        differences[str(MASS_HINV[i])]={"same_z095_M_c200c_and_exterior_Phi_anchor_gap_max":max(
                                           gapmass,gapphi,gapc),
                                       "dlnM_dln_a_at_z095_difference_abs":gaprate,
                                       "rate_alpha_frozen_choices":list(ALPHA)}
    result={"date":"2026-09-28",
            "status":"E17D2B1_TWO_WECHSLER_TYPE_CONDITIONAL_HALO_MASS_GROWTH_AND_EXTERIOR_POISSON_HISTORY_QA_CERTIFIED_FULL_PHYSICAL_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "all_immutable_parents":p["immutable_parents"],
            "original_four_trajectories_3_math_nodes_each":reports,
            "conditional_mass_references":list(MASS_HINV),
            "mathematical_rate_parameters_NOT_empirically_fitted":list(ALPHA),
            "same_exact_anchor_two_distinct_exterior_growth_rate_witnesses":differences,
            "individual_halo_M200c_and_cvir_ac_relation_NOT_calibrated":True,
            "actual_neutrino_halo_tracer_dynamics_calculated":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(out/"e17d2b1_original_joint_conditional_exterior_assembly_rates.json",result)
    print("E17D2B1_ORIGINAL_CONDITIONAL_ANCHORED_EXTERIOR_DERIVATIVE_WITNESS_PASS",
          differences,"FULL_PHYSICAL_B_BLOCKED",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--case-index",type=int)
    g.add_argument("--aggregate",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    x=ap.parse_args()
    if x.self_test:
        p,previous=gate()
        for j in range(2):
            _,_,a=original_anchor(j,p,previous)
            require(abs(mass_at_z(a["M200c_physical_Msun_DM14"],ALPHA[0],Z0)/
                        a["M200c_physical_Msun_DM14"]-1.)<1e-12,
                    "original halo anchor no longer closes")
        print("E17D2B1_ORIGINAL_PROSPECTIVE_SHA_WECHSLER_SHAPE_AND_ORIGINAL_MEDIAN_ANCHOR_GATE_PASS",
              "PHYSICAL_HALO_ASSEMBLY_AND_B_BLOCKED",flush=True)
    elif x.aggregate:aggregate(x.output_dir)
    else:run_case(x.case_index,x.output_dir)
if __name__=="__main__":main()
