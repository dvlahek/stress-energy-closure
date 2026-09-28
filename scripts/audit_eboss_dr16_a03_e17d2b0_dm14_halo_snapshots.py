#!/usr/bin/env python3
"""E17D2b0 DM14 population-median c200c(M,z) conditional NFW snapshots.

This is not actual LRG/ELG halo mass or a single halo mass-accretion history.
Dutton-Maccio 2014 is Planck2013 DM-only, NOT exact original massive-neutrino
CLASS cosmology. Both h values are explicit. No physical B or odd reading.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e17d2a_nfw_spatial_poisson_halo as nfw

PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b0_dm14_population_halo_mass_concentration_prereg_2026-09-28.json"
PRE_BLOB="8faaf167403e09f321a6c62f452120f8b96898a2"
NFW_RUNNER_BLOB="5d7248095273f8ab0dc74cb9f95e7f78200044fa"
NFW_ARCH=ROOT/"source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28"
NFW_JOINT=NFW_ARCH/"e17d2a_original_joint_nfw_static_spatial_poisson_family.json"
NFW_IND=NFW_ARCH/"e17d2a_independent_stdlib_direct_radial_fourier_replay.json"
NFW_MAN=NFW_ARCH/"archive_manifest.json"
E17A_ARCH=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27"
E17A=E17A_ARCH/"e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17A_MAN=E17A_ARCH/"archive_manifest.json"
STATES=("FD","plus","minus")
MASSR=(1.,10.)
ZS=(.945,.95,.955)
DEX=(-.11,0.,.11)
HDM=.671
HCLASS=.6736
OM=.3175
OL=.6825
G=4.30091e-9 # Mpc*(km/s)^2/Msun
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b0_dm14_population"

def require(ok,why):
    if not ok:raise ValueError("E17D2B0_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"immutable conditional halo report changed")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b0_",delete=False) as fh:
            t=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B0_ORIGINAL_CONDITIONAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective E17D2b0 protocol Git blob")
    p=json.loads(PRE.read_bytes());pa=p["immutable_source_parents"]
    require(blob(Path(nfw.__file__).read_bytes())==NFW_RUNNER_BLOB
            and blob(NFW_MAN.read_bytes())==pa["E17D2a_archive_manifest_git_blob"]
            and blob(nfw.PRE.read_bytes())==pa["E17D2a_protocol_git_blob"]
            and blob(E17A_MAN.read_bytes())==pa["E17A_manifest_git_blob"]
            and blob(nfw.ROOT.joinpath(
                "source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json"
                ).read_bytes())==pa["E17A_protocol_git_blob"],
            "original source runners/manifests/protocols not immutable")
    for path,name in ((nfw.E8,"E8_original_4000q_CSV_sha256"),
                      (nfw.E16,"E16_original_full_sha256"),
                      (nfw.D1,"E17D1_original_joint_sha256"),
                      (NFW_JOINT,"E17D2a_original_joint_sha256"),
                      (NFW_IND,"E17D2a_independent_sha256"),
                      (E17A,"E17A_original_joint_sha256")):
        require(sha(path.read_bytes())==pa[name],"original science source SHA "+name)
    a=json.loads(E17A.read_bytes())
    d=json.loads(NFW_JOINT.read_bytes())
    require(a["geometry_count"]==576
            and a["full_physical_finite_K_bispectrum"]=="BLOCKED"
            and a["observed_odd_read"] is False
            and d["original_E16_geometries_per_case"]==576
            and d["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and d["observed_odd_SEALED"] is True,
            "original science STOP/source geometry altered")
    f=p["frozen_geometry"]
    config=p["precommitted_evaluation_only"]
    dm=p["halo_background_units"]
    require(f["original_E16_closed_geometries"]==576
            and f["k_short_original_CLASS_h_per_Mpc"]==list(nfw.KS)
            and f["K_long_original_CLASS_h_per_Mpc"]==list(nfw.KLS)
            and f["mu_short"]==list(nfw.MS)
            and f["mu_long"]==list(nfw.ML)
            and f["redshift_math_nodes_NOT_eBOSS_z_eff"]==list(ZS)
            and f["original_states"]==list(STATES)
            and f["original_tracer_order"]==["LRG","ELG"]
            and f["original_CLASS_h"]==HCLASS
            and config["halo_mass_ratio_to_published_1e12_hDM_inverse_Msun"]==list(MASSR)
            and config["redshift_nodes"]==list(ZS)
            and config["log10_concentration_offset_dex_illustrative"]==list(DEX)
            and config["cases"]==18 and config["mode_values"]==31104
            and dm["h_DM14"]==HDM
            and dm["Omega_m_DM14"]==OM
            and dm["Omega_Lambda_flat_DM14"]==OL
            and dm["physical_G_Mpc_kms2_per_Msun"]==G
            and all(p["absolute_STOP"].values()),
            "precommitted physical calibration or governance STOP changed")
    nfw.gate()
    return p

def evaluations():
    return [(mass,z,dex) for mass in MASSR for z in ZS for dex in DEX]

def case_filename(index):
    return "e17d2b0_original_DM14_conditional_snapshot_%02d.json"%index

def coefficients(z):
    require(0.<=z<=5.,"DM14 published redshift range")
    return (.520+(.905-.520)*math.exp(-.617*z**1.21),
            -.101+.026*z)

def concentration_median(mass_ratio,z):
    require(mass_ratio in MASSR,"precommitted mass anchors only")
    a,b=coefficients(z)
    return 10.**(a+b*math.log10(mass_ratio))

def halo_snapshot(mass_ratio,z,dex):
    require(dex in DEX,"precommitted dex bracket only")
    a,b=coefficients(z)
    median=concentration_median(mass_ratio,z)
    c=median*10.**dex
    mass_hinv=1.e12*mass_ratio # Msun/h_DM
    mass_phys=mass_hinv/HDM
    E2=OM*(1.+z)**3+OL
    H=100.*HDM*math.sqrt(E2)  # km/s/Mpc
    rho=3.*H*H/(8.*math.pi*G) # Msun/Mpc^3
    r200=(3.*mass_phys/(4.*math.pi*200.*rho))**(1./3.) # Mpc physical
    rs_phys=r200/c
    rs_com=rs_phys*(1.+z)
    reclosure=(4.*math.pi/3.)*200.*rho*r200**3/mass_phys
    grav=G*mass_phys/rs_phys # (km/s)^2
    require(math.isfinite(c) and c>0 and r200>0 and rs_phys>0
            and abs(reclosure-1.)<1e-12,
            "published c200c and critical density radius mass closure")
    return {"M200c_hDM_inverse_Msun_QA_anchor":mass_hinv,
            "M200c_physical_Msun_DM14":mass_phys,
            "z_original_CLASS_math_node_only":z,
            "log10_c200c_population_median_DM14":math.log10(median),
            "fit_a_of_z":a,"fit_b_of_z":b,
            "c200c_population_median_DM14":median,
            "illustrative_log10_concentration_offset_dex_NOT_z095_posterior":dex,
            "c200c_conditional_illustrative":c,
            "E_DM14_z_squared":E2,
            "H_DM14_km_s_Mpc":H,
            "rho_crit_DM14_Msun_Mpc3":rho,
            "r200c_physical_Mpc_DM14":r200,
            "r_s_physical_Mpc_DM14":rs_phys,
            "r_s_comoving_Mpc_DM14":rs_com,
            "static_GM_over_rs_km2_s2_ONLY_not_actual_eBOSS_halo":grav,
            "r200_critical_mass_reclosure_rel":abs(reclosure-1.),
            "h_original_CLASS_over_h_DM14_unit_ratio":HCLASS/HDM}

def source_mode(mod_hCLASS,z,rs_com,c):
    require(mod_hCLASS>0,"strict nonzero original Fourier sample")
    kphys=mod_hCLASS*HCLASS # CLASS original comoving physical 1/Mpc
    x=kphys*rs_com
    u=nfw.analytic_u(x,c)
    phi=-u/(x*x)
    return {"k_com_original_CLASS_h_per_Mpc":mod_hCLASS,
            "k_com_physical_inv_Mpc_from_original_CLASS_h":kphys,
            "dimensionless_x_mixed_cosmology_QA_only":x,
            "normalized_static_NFW_density_Fourier_u_only":u,
            "dimensionless_Poisson_Phi_shape_minus_u_over_x2_ONLY":phi,
            "conditional_point_mass_shape_deficit_abs":abs(1.-u)}

def run_case(index,out):
    p=gate();cases=evaluations()
    require(index in range(len(cases)),"case not preregistered")
    mass,z,dex=cases[index]
    halo=halo_snapshot(mass,z,dex)
    c=halo["c200c_conditional_illustrative"]
    rs=halo["r_s_comoving_Mpc_DM14"]
    maximum_deficit=0.
    maximum_abs_u=0.
    maximum_closure=0.
    max_reverse=0.
    geometry=nfw.geometry()
    records=[]
    for g in geometry:
        legs=[source_mode(mod,z,rs,c)
              for mod in g["original_short_leg_moduli_h_Mpc"]]
        long=source_mode(g["original_long_K_h_Mpc"],z,rs,c)
        for o in legs+[long]:
            maximum_abs_u=max(maximum_abs_u,
                              abs(o["normalized_static_NFW_density_Fourier_u_only"]))
            x=o["dimensionless_x_mixed_cosmology_QA_only"]
            max_reverse=max(max_reverse,abs(
                o["normalized_static_NFW_density_Fourier_u_only"]-
                nfw.analytic_u(-x,c)))
            maximum_deficit=max(maximum_deficit,
                                o["conditional_point_mass_shape_deficit_abs"])
        maximum_closure=max(maximum_closure,
                            g["original_triangle_closure_h_Mpc"])
        records.append({"k":g["k"],"K":g["K"],"mu_s":g["mu_s"],
                        "mu_L":g["mu_L"],"phi_index":g["phi_index"],
                        "original_short_leg_moduli_h_Mpc":g[
                            "original_short_leg_moduli_h_Mpc"],
                        "original_long_K_h_Mpc":g["original_long_K_h_Mpc"],
                        "original_triangle_closure_h_Mpc":g[
                            "original_triangle_closure_h_Mpc"],
                        "actual_two_short_legs_static_spatial_form_factors":legs,
                        "original_long_K_static_spatial_form_factor":long})
    qa=p["preregistered_QA"]
    require(len(records)==576
            and maximum_closure<qa["triangle_closure_abs_h_per_Mpc_max"]
            and maximum_abs_u<qa["model_nfw_Fourier_u_positive_mass_abs_bound"]
            and max_reverse<qa["fourier_reversal_abs_gap_max"],
            "frozen geometry or positive-density static Fourier QA")
    result={"date":"2026-09-28",
            "status":"E17D2B0_ORIGINAL_DM14_PUBLISHED_CONDITIONAL_MASS_CONCENTRATION_RADIUS_AND_NFW_SPATIAL_SNAPSHOTS_ONLY_FULL_B_BLOCKED",
            "case_index":index,
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_E17D2a_joint_sha256":p["immutable_source_parents"][
                "E17D2a_original_joint_sha256"],
            "external_DM14_cosmology_NOT_original_CLASS_cosmology":True,
            "original_CLASS_and_DM14_h_equal":False,
            "actual_eBOSS_LRG_ELG_halo_mass_or_HOD_selected":False,
            "same_mass_across_z_is_NOT_one_halo_mass_accretion_history":True,
            "no_true_halo_formation_or_potential_time_history":True,
            "published_population_halo_snapshot":halo,
            "original_576_E16_geometry_and_two_legs_plus_long_QA":records,
            "original_E16_geometry_count":len(records),
            "actual_short_short_long_values":3*len(records),
            "QA":{"max_original_triangle_closure_h_per_Mpc":maximum_closure,
                  "max_abs_normalized_positive_NFW_u":maximum_abs_u,
                  "max_reverse_abs":max_reverse,
                  "max_conditional_static_point_mass_shape_deficit_report_only":maximum_deficit,
                  "critical_mass_reclosure_rel":halo["r200_critical_mass_reclosure_rel"]},
            "physical_consistent_CLASS_plus_DM14_neutrino_halo_coupling":"BLOCKED",
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,
            "no_new_CLASS_FITS_mock":True,
            "main_untouched_PR_draft":True}
    save_once(out/case_filename(index),result)
    print("E17D2B0_ORIGINAL_CONDITIONAL_DM14_CASE",index,
          "MASS_RATIO",mass,"Z",z,"DEX_ILLUSTRATIVE",dex,
          "C",c,"RS_COM_Mpc",rs,
          "N",len(records)*3,
          "MAX_POINT_MASS_SHAPE_DEFICIT",maximum_deficit,
          "HALO_HISTORY_AND_GALAXY_B_BLOCKED",flush=True)

def aggregate(out):
    p=gate();cases=evaluations()
    reports={};maxdiff=0.
    for i,(mass,z,dex) in enumerate(cases):
        raw=(out/case_filename(i)).read_bytes();r=json.loads(raw)
        h=r["published_population_halo_snapshot"]
        require(r["case_index"]==i and r["original_E16_geometry_count"]==576
                and r["actual_short_short_long_values"]==1728
                and h["z_original_CLASS_math_node_only"]==z
                and h["illustrative_log10_concentration_offset_dex_NOT_z095_posterior"]==dex
                and h["M200c_hDM_inverse_Msun_QA_anchor"]==1.e12*mass
                and r["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
                and r["observed_odd_SEALED"] is True,
                "missing original case, identity or mandatory physical STOP")
        reports[case_filename(i)]={"sha256":sha(raw),"size_bytes":len(raw),
                                   "M_hinv_Msun_QA_only":h["M200c_hDM_inverse_Msun_QA_anchor"],
                                   "z_math_only":z,
                                   "dex_sensitivity_only":dex,
                                   "c":h["c200c_conditional_illustrative"],
                                   "r_s_comoving_Mpc_DM14":h["r_s_comoving_Mpc_DM14"],
                                   "QA":r["QA"]}
        maxdiff=max(maxdiff,r["QA"][
            "max_conditional_static_point_mass_shape_deficit_report_only"])
    result={"date":"2026-09-28",
            "status":"E17D2B0_PUBLISHED_DM14_POPULATION_C200C_MASS_RADIUS_CONDITIONAL_18_SNAPSHOTS_DONE_FULL_PHYSICAL_HALO_HISTORY_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "immutable_original_parent_identifiers":p["immutable_source_parents"],
            "published_DM14_planck2013_cosmology":{
                "Omega_m":OM,"Omega_Lambda_flat":OL,"h":HDM},
            "original_CLASS_h":HCLASS,
            "cosmology_h_fractional_label_mismatch":abs(HCLASS/HDM-1.),
            "number_of_frozen_conditional_halo_snapshots":len(cases),
            "original_E16_geometry_per_snapshot":576,
            "two_actual_short_and_one_long_mode_per_geometry":3,
            "total_conditional_spatial_Fourier_modes":len(cases)*576*3,
            "maximum_conditional_point_mass_spatial_u_deficit_in_declared_mass_QA_only_cases":maxdiff,
            "all_original_case_SHA":reports,
            "cross_redshift_same_mass_snapshots_NOT_single_halo_assembly":True,
            "individual_LRG_ELG_halo_M_HOD_time_history_known":False,
            "physically_consistent_neutrino_CLASS_plus_DM14_halo_combination":"BLOCKED",
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(out/"e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json",result)
    print("E17D2B0_ORIGINAL_DM14_CONDITIONAL_POPULATION_SNAPSHOT_DONE",
          len(cases),len(cases)*576*3,"MAX_U_DEFICIT",maxdiff,
          "FULL_PHYSICAL_HALO_HISTORY_AND_B_BLOCKED",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--case-index",type=int)
    g.add_argument("--aggregate",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    args=ap.parse_args()
    if args.self_test:
        p=gate()
        require(len(evaluations())==18
                and len(nfw.geometry())==576
                and abs(halo_snapshot(1.,.95,0.)["r200_critical_mass_reclosure_rel"])<1e-12,
                "original 18-case, frozen geometry or DM14 200 critical radius")
        print("E17D2B0_ORIGINAL_PUBLISHED_DM14_AND_IMMUTABLE_SOURCE_PREFLIGHT_PASS",
              "NO_GALAXY_HOD_NO_HALO_TIME_NO_PHYSICAL_B",flush=True)
    elif args.aggregate:aggregate(args.output_dir)
    else:run_case(args.case_index,args.output_dir)
if __name__=="__main__":
    main()
