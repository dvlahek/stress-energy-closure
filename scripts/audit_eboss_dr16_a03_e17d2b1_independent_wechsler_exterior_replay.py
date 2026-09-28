#!/usr/bin/env python3
"""Independent pure-stdlib E17D2b1 conditional exterior source-history audit.

No original E17D2b1 runner imports, NumPy, SciPy or CLASS; independently
evaluate original E17D2b0 SHA-pinned 200c anchors and algebraic derivatives
of conditional Wechsler-type M(a) on only original z=.945,.95,.955 nodes.
Does not interpret c_vir as c_200c or treat test parameters as halo priors.
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
B0=ROOT/"source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28"
B0_JOINT=B0/"e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json"
B0_INDEP=B0/"e17d2b0_independent_stdlib_dm14_full_radial_replay.json"
B0_MAN=B0/"archive_manifest.json"
B0_PROTO=ROOT/"source_data/eboss_dr16_a03_e17d2b0_dm14_population_halo_mass_concentration_prereg_2026-09-28.json"
B0_RUNNER=ROOT/"scripts/audit_eboss_dr16_a03_e17d2b0_dm14_halo_snapshots.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
D1=ROOT/"source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_original_joint_causal_history_nonidentifiability.json"
D2A=ROOT/"source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28/e17d2a_original_joint_nfw_static_spatial_poisson_family.json"
MASSES=(1000000000000,10000000000000)
INDICES=(4,13)
GROWTH=(.4,.8)
Z=(.945,.95,.955);Z0=.95
H0=67.1;OM=.3175;OL=.6825
G=4.30091e-9
MPC_KM=3.0856775814913673e19
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b1_wechsler_exterior"

def check(ok,why):
    if not ok:raise ValueError("E17D2B1_INDEPENDENT_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"prospective E17D2b1 physical STOP Git blob")
    p=json.loads(PRE.read_bytes());parent=p["immutable_parents"]
    for path,key in ((E8,"E8_original_4000q_CSV_sha256"),
                     (E16,"E16_original_full_sha256"),
                     (D1,"E17D1_original_joint_sha256"),
                     (D2A,"E17D2a_original_joint_sha256"),
                     (B0_JOINT,"E17D2b0_original_joint_sha256"),
                     (B0_INDEP,"E17D2b0_independent_sha256")):
        check(sha(path.read_bytes())==parent[key],"immutable source SHA "+key)
    check(blob(B0_MAN.read_bytes())==parent["E17D2b0_archive_manifest_git_blob"]
          and blob(B0_PROTO.read_bytes())==parent["E17D2b0_protocol_git_blob"]
          and blob(B0_RUNNER.read_bytes())==parent["E17D2b0_original_runner_git_blob"],
          "immutable original DM14 archive/protocol/runner Git blob")
    original=json.loads(B0_JOINT.read_bytes())
    check(original["number_of_frozen_conditional_halo_snapshots"]==18
          and original["total_conditional_spatial_Fourier_modes"]==31104
          and original["individual_LRG_ELG_halo_M_HOD_time_history_known"] is False
          and original["observed_odd_SEALED"] is True
          and original["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED",
          "original DM14 POPULATION source-only parent identity")
    f=p["frozen"]
    check(f["original_mass_reference_hDM_inverse_Msun"]==list(MASSES)
          and f["original_DM14_anchor_case_indices"]==list(INDICES)
          and f["Wechsler_shape_alpha_equals_2_a_c_QA_only"]==list(GROWTH)
          and f["z0_anchor_only"]==Z0
          and f["fixed_original_E17A_z_math_nodes_NOT_real_halo_formation_samples"]==list(Z)
          and f["physical_exterior_test_radius_ratio_to_anchor_r200c"]==3
          and math.isclose(f["DM14_h"],H0/100.,rel_tol=0.,abs_tol=1e-15)
          and f["DM14_Omega_m"]==OM
          and f["DM14_Omega_lambda"]==OL
          and f["G_Mpc_km2_s2_Msun"]==G
          and all(p["absolute_STOP"].values()),
          "prospective time/halo source/growth QA or physical STOP changed")
    return p,original

def anchormass(mass_index,p,b0):
    index=INDICES[mass_index]
    path=B0/("e17d2b0_original_DM14_conditional_snapshot_%02d.json"%index)
    raw=path.read_bytes()
    key=("E17D2b0_anchor_mass_1e12_z095_original_SHA" if mass_index==0
         else "E17D2b0_anchor_mass_1e13_z095_original_SHA")
    check(sha(raw)==p["immutable_parents"][key]
          and sha(raw)==b0["all_original_case_SHA"][path.name]["sha256"],
          "original E17D2b0 median case SHA")
    a=json.loads(raw)
    check(a["case_index"]==index and a["observed_odd_SEALED"] is True
          and a["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED",
          "original anchored POPULATION halo source scope")
    b=a["published_population_halo_snapshot"]
    check(b["M200c_hDM_inverse_Msun_QA_anchor"]==MASSES[mass_index]
          and b["z_original_CLASS_math_node_only"]==Z0
          and b["illustrative_log10_concentration_offset_dex_NOT_z095_posterior"]==0.,
          "original median mass reference or z snapshot changed")
    return sha(raw),b

def independent_snapshot(M0,parameter,z,R):
    """Independently derived Wechsler form and critical-density shell theorem."""
    growth_exp=-parameter*(z-Z0)/(1.+Z0)
    mass=M0*math.exp(growth_exp)
    H=H0*math.sqrt(OM*(1.+z)**3+OL)
    critical=(3./(8.*math.pi*G))*H*H
    radius=(mass/(200.*critical*(4.*math.pi/3.)))**(1./3.)
    potential=-(G/R)*mass
    radial_gradient=G/(R*R)*mass
    dlogmass_dz=-parameter/(1.+Z0)
    dloga=parameter*(1.+z)/(1.+Z0)
    return {"mass_M200c_physical_Msun_conditional_ansatz":mass,
            "mass_ratio_to_same_frozen_anchor":math.exp(growth_exp),
            "DM14_H_km_s_Mpc":H,
            "DM14_rho_crit_Msun_Mpc3":critical,
            "r200c_phys_Mpc_conditional_ansatz":radius,
            "fixed_phys_exterior_R_Mpc_QA_only":R,
            "exterior_margin_R_over_r200c":R/radius,
            "exterior_Phi_phys_km2_s2":potential,
            "exterior_inward_acceleration_magnitude_km2_s2_per_Mpc":radial_gradient,
            "exact_dlnM_dz_and_dlnAbsPhi_dz":dlogmass_dz,
            "exact_dlnM_dln_a_and_dlnAbsPhi_dln_a":dloga,
            "exact_dlnM_dt_seconds_inverse_DM14_only":dloga*(H/MPC_KM),
            "exact_dPhi_dz_at_fixed_physical_R_km2_s2":potential*dlogmass_dz}

def full_replay(out,joint,p,b0):
    check(joint["prospective_protocol_git_blob"]==PRE_BLOB
          and len(joint["original_four_trajectories_3_math_nodes_each"])==4
          and joint["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and joint["observed_odd_SEALED"] is True
          and joint["individual_halo_M200c_and_cvir_ac_relation_NOT_calibrated"] is True,
          "joint original four source-only histories and physical STOP drift")
    mx=0.;reports={};n=0
    q=p["prospectively_locked_QA"]
    for i in range(4):
        name="e17d2b1_original_conditional_wechsler_exterior_%02d.json"%i
        raw=(out/name).read_bytes()
        check(sha(raw)==joint["original_four_trajectories_3_math_nodes_each"][name]["sha256"],
              "immutable original E17D2b1 report SHA mismatch: "+name)
        row=json.loads(raw)
        parsha,anchor=anchormass(i//2,p,b0)
        alpha=GROWTH[i%2]
        R=3.*anchor["r200c_physical_Mpc_DM14"]
        check(row["original_E17D2b0_anchor_SHA256"]==parsha
              and row["original_source_mass_reference_hDM_inverse_Msun"]==MASSES[i//2]
              and row["Wechsler_form_alpha_equals_two_ac_math_QA_NOT_fitted_to_halo"]==alpha
              and row["original_DM14_median_c200c_at_z0_ONLY_not_evolved"]==
                  anchor["c200c_population_median_DM14"]
              and row["halo_M200c_use_of_Wechsler_virial_M_history_UNCALIBRATED_ansatz"] is True
              and row["same_DM14_median_c200c_at_z0_but_no_cvir_conversion_or_evolved_c"] is True
              and row["physical_CLASS_plus_DM14_cosmology_consistent"] is False
              and row["actual_LRG_ELG_halo_assembly_HOD_known"] is False
              and row["observed_odd_SEALED"] is True
              and row["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED",
              "original conditional case provenance/independent parent or absolute STOP")
        check(len(row["snapshots_sorted_z_increasing"])==3,
              "original frozen three z nodes")
        for z,stored in zip(Z,row["snapshots_sorted_z_increasing"]):
            check(stored["z_frozen_original_math_node"]==z,
                  "original redshift is not originally frozen E17A math node")
            independent=independent_snapshot(
                anchor["M200c_physical_Msun_DM14"],alpha,z,R)
            for key,v in independent.items():
                old=stored[key]
                mx=max(mx,abs(old-v)/max(1.,abs(old),abs(v)))
            check(stored["exterior_margin_R_over_r200c"]>1.
                  and stored["exterior_Phi_phys_km2_s2"]<0,
                  "exterior of shell theorem/positive mass condition")
            close=stored["r200c_phys_Mpc_conditional_ansatz"]
            density=stored["DM14_rho_crit_Msun_Mpc3"]
            mass=stored["mass_M200c_physical_Msun_conditional_ansatz"]
            reconstruction=200.*density*(4.*math.pi/3.)*close**3
            check(abs(reconstruction/mass-1.)<
                  q["mass_radius_200critical_reclosure_relative_max"],
                  "independent 200critical mass closure")
            n+=1
        records=row["snapshots_sorted_z_increasing"]
        dlog=(math.log(records[2]["mass_M200c_physical_Msun_conditional_ansatz"])-
              math.log(records[0]["mass_M200c_physical_Msun_conditional_ansatz"]))/(Z[2]-Z[0])
        rate=-alpha/(1.+Z0)
        mx=max(mx,abs(dlog-rate))
        check(abs(dlog-rate)<q[
            "exact_analytic_log_derivative_vs_two_sided_original_z_fd_abs_max"],
            "original three-node central logarithmic mass derivative")
        reports[name]={"original_sha256":sha(raw),"mass_anchor_hDM_inverse_Msun":MASSES[i//2],
                       "alpha_shape_math_only":alpha,"three_original_math_z_samples":3}
        print("E17D2B1_INDEPENDENT_ORIGINAL_CASE",i,
              "REPLAYED_MATH_NODES",3,
              "MAX_SCALED_GAP",mx,flush=True)
    check(n==12,"full four conditional histories × original three math nodes")
    check(mx<q["full_independent_stdlib_four_trajectories_max_scaled_gap"],
          "independent full four history + exterior Poisson closure mismatch")
    for j in range(2):
        x=json.loads((out/("e17d2b1_original_conditional_wechsler_exterior_%02d.json"%(2*j))).read_bytes())
        y=json.loads((out/("e17d2b1_original_conditional_wechsler_exterior_%02d.json"%(2*j+1))).read_bytes())
        ax=x["snapshots_sorted_z_increasing"][1]
        ay=y["snapshots_sorted_z_increasing"][1]
        gap=abs(ax["exact_dlnM_dln_a_and_dlnAbsPhi_dln_a"]-
                ay["exact_dlnM_dln_a_and_dlnAbsPhi_dln_a"])
        check(gap>q["distinct_alpha_growth_rate_gap_at_z095_min"]
              and ax["mass_M200c_physical_Msun_conditional_ansatz"]==
                  ay["mass_M200c_physical_Msun_conditional_ansatz"]
              and ax["exterior_Phi_phys_km2_s2"]==
                  ay["exterior_Phi_phys_km2_s2"],
              "independent same anchored M/Phi and distinct model growth rates")
    return mx,n,reports

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"independent immutable certificate collision")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b1i_",delete=False) as fh:
            temp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17D2B1_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--negative-control",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    x=ap.parse_args()
    p,b0=gate()
    if x.self_test:
        for i in range(2):
            sha_value,anchor=anchormass(i,p,b0)
            s=independent_snapshot(anchor["M200c_physical_Msun_DM14"],GROWTH[0],Z0,
                                   3.*anchor["r200c_physical_Mpc_DM14"])
            check(abs(s["r200c_phys_Mpc_conditional_ansatz"]/
                      anchor["r200c_physical_Mpc_DM14"]-1.)<1e-12,
                  "original DM14 reference not recovered")
        print("E17D2B1_INDEPENDENT_PROSPECTIVE_SHA_AND_EXACT_ANCHOR_GATE_PASS",
              "NO_CVIR_C200C_INTERCHANGE_NO_PHYSICAL_B",flush=True)
        return
    path=x.output_dir/"e17d2b1_original_joint_conditional_exterior_assembly_rates.json"
    raw=path.read_bytes();joint=json.loads(raw)
    if x.negative_control:
        fake=json.loads(raw)
        fake["original_four_trajectories_3_math_nodes_each"][
            "e17d2b1_original_conditional_wechsler_exterior_00.json"]["sha256"]="0"*64
        try:full_replay(x.output_dir,fake,p,b0)
        except ValueError as err:
            check("immutable original E17D2B1 report SHA mismatch" in str(err),
                  "negative control rejected at wrong gate")
            print("E17D2B1_INDEPENDENT_TAMPERED_ORIGINAL_SHA_REJECTED",flush=True)
            return
        raise AssertionError("deliberately mutated full original SHA was accepted")
    mx,n,reports=full_replay(x.output_dir,joint,p,b0)
    result={"date":"2026-09-28",
            "status":"E17D2B1_INDEPENDENT_PURE_STDLIB_FOUR_CONDITIONAL_WECHSLER_EXTERIOR_HISTORY_AND_FROZEN_PARENT_REPLAY_PASS_FULL_PHYSICAL_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_joint_sha256":sha(raw),
            "original_case_SHA_and_frozen_node_counts":reports,
            "mass_references":list(MASSES),"math_growth_rate_anchors":list(GROWTH),
            "full_independent_conditional_exterior_redshift_snapshots":n,
            "maximum_scaled_independent_source_history_radius_force_gap":mx,
            "method":"Pure standard-library separate E17D2b0 SHA-pinned anchor reader, independent exponential Wechsler-type CONDITIONAL form, DM14 critical mass/radius, same fixed exterior Newtonian shell theorem and exact/fd derivative; no original E17D2b1 import, SciPy/NumPy/CLASS, merger tree or observed odd.",
            "Wechsler_virial_mass_to_DM14_M200c_UNCALIBRATED_shape_only":True,
            "CLASS_DM14_background_mismatch_preserved":True,
            "individual_halo_assembly_or_LRG_ELG_HOD_inferred":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(x.output_dir/"e17d2b1_independent_stdlib_exterior_assembly_replay.json",result)
    print("E17D2B1_INDEPENDENT_ALL_4x3_CONDITIONAL_OUTER_MONOPOLE_HISTORY_PASS",
          "MAX",mx,"PHYSICAL_HALO_HISTORY_AND_B_BLOCKED",flush=True)
if __name__=="__main__":main()
