#!/usr/bin/env python3
"""Independent pure-stdlib DM14 conditional snapshot and direct NFW Fourier audit.

No imports of E17D2b0 original worker, SciPy, NumPy or CLASS. Reconstruct
published Dutton-Maccio c200c fit, critical-density radius, original E16
triangle algebra and full positive radial NFW Fourier profile independently.
Conditional population mass samples are NOT actual LRG/ELG HOD/halo histories.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b0_dm14_population_halo_mass_concentration_prereg_2026-09-28.json"
PRE_BLOB="8faaf167403e09f321a6c62f452120f8b96898a2"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A_ARCH=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27"
E17A=E17A_ARCH/"e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17A_MAN=E17A_ARCH/"archive_manifest.json"
NFW_ARCH=ROOT/"source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28"
NFW=NFW_ARCH/"e17d2a_original_joint_nfw_static_spatial_poisson_family.json"
NFW_IND=NFW_ARCH/"e17d2a_independent_stdlib_direct_radial_fourier_replay.json"
NFW_MAN=NFW_ARCH/"archive_manifest.json"
D1=ROOT/"source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_original_joint_causal_history_nonidentifiability.json"
MASS=(1.,10.)
Z=(.945,.95,.955)
DEX=(-.11,0.,.11)
H=.671;HC=.6736
OM=.3175;OL=1.-OM
G=4.30091e-9
KS=(.05,.075,.1);KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.);ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi)
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b0_dm14_population"

def check(ok,why):
    if not ok:raise ValueError("E17D2B0_INDEPENDENT_FAIL_CLOSED: "+why)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"prospective DM14 protocol Git blob")
    p=json.loads(PRE.read_bytes());par=p["immutable_source_parents"]
    for path,key in ((E8,"E8_original_4000q_CSV_sha256"),
                     (E16,"E16_original_full_sha256"),
                     (E17A,"E17A_original_joint_sha256"),
                     (D1,"E17D1_original_joint_sha256"),
                     (NFW,"E17D2a_original_joint_sha256"),
                     (NFW_IND,"E17D2a_independent_sha256")):
        check(sha(path.read_bytes())==par[key],"frozen original source SHA: "+key)
    check(blob(NFW_MAN.read_bytes())==par["E17D2a_archive_manifest_git_blob"]
          and blob(E17A_MAN.read_bytes())==par["E17A_manifest_git_blob"],
          "original E17D2a or E17A archive Git blob")
    a=json.loads(E17A.read_bytes());n=json.loads(NFW.read_bytes())
    check(a["geometry_count"]==576
          and a["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED"
          and a["eBOSS_observed_odd_read"] is False
          and n["original_E16_geometries_per_case"]==576
          and n["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and n["observed_odd_SEALED"] is True,
          "source-only parent actual science STOP")
    f=p["frozen_geometry"];m=p["precommitted_evaluation_only"]
    c=p["halo_background_units"]
    check(f["k_short_original_CLASS_h_per_Mpc"]==list(KS)
          and f["K_long_original_CLASS_h_per_Mpc"]==list(KLS)
          and f["mu_short"]==list(MS) and f["mu_long"]==list(ML)
          and f["original_E16_closed_geometries"]==576
          and f["redshift_math_nodes_NOT_eBOSS_z_eff"]==list(Z)
          and f["original_CLASS_h"]==HC
          and f["original_tracer_order"]==["LRG","ELG"]
          and m["halo_mass_ratio_to_published_1e12_hDM_inverse_Msun"]==list(MASS)
          and m["redshift_nodes"]==list(Z)
          and m["log10_concentration_offset_dex_illustrative"]==list(DEX)
          and m["cases"]==18 and m["mode_values"]==31104
          and c["h_DM14"]==H and c["Omega_m_DM14"]==OM
          and c["Omega_Lambda_flat_DM14"]==OL
          and c["physical_G_Mpc_kms2_per_Msun"]==G
          and all(p["absolute_STOP"].values()),
          "preregistered mass/redshift/scatter/Hubble/STOP scope changed")
    return p

def snapshot(mass,z,dex):
    # Independent algebraic evaluation, no import of original fit.
    pivotlog=math.log10(mass)
    az=.520+.385*math.exp(-.617*math.pow(z,1.21))
    bz=.026*z-.101
    median=math.pow(10.,az+bz*pivotlog)
    conc=math.pow(10.,az+bz*pivotlog+dex)
    mass_h_inv=1.e12*mass
    mass_physical=mass_h_inv/H
    hubble_z=100.*H*math.sqrt(OM*(1.+z)**3+OL)
    rho_crit0=3.*(100.*H)**2/(8.*math.pi*G)
    rho=rho_crit0*(OM*(1.+z)**3+OL)
    radius=(mass_physical/((4.*math.pi/3.)*200.*rho))**(1./3.)
    rs=radius/conc
    return {"a":az,"b":bz,"median":median,"c":conc,"mass":mass_physical,
            "H":hubble_z,"rho":rho,"r200":radius,"rs":rs,
            "rs_com":(1.+z)*rs,"GMr":G*mass_physical/rs}

def geometry():
    for k in KS:
        for K in KLS:
            ratio=K/k
            for ms in MS:
                for ml in ML:
                    for index,phi in enumerate(PH):
                        dot=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(
                            max(0.,1.-ml*ml))*math.cos(phi)
                        yield (k,K,ms,ml,index,
                               k*math.sqrt(1.+ratio*ratio/4.-ratio*dot),
                               k*math.sqrt(1.+ratio*ratio/4.+ratio*dot))

def nfw_direct_radial(x,c,panels):
    """Independent positive-density angular transform via even Simpson panels."""
    check(panels%2==0 and panels>=1024,"radial Simpson independent panel count")
    acc=math.log1p(c)-c/(1.+c)
    dx=c/panels
    def f(y):
        if y==0.:return 0.
        t=x*y
        sinc=1. if t==0. else math.sin(t)/t
        return y/(1.+y)**2*sinc
    odd=[]
    even=[]
    for j in range(1,panels):
        sample=f(dx*j)
        (odd if j%2 else even).append(sample)
    return dx/3.*(f(c)+4.*math.fsum(odd)+2.*math.fsum(even))/acc

def full_replay(original,joint,p):
    check(joint["number_of_frozen_conditional_halo_snapshots"]==18
          and joint["original_E16_geometry_per_snapshot"]==576
          and joint["total_conditional_spatial_Fourier_modes"]==31104
          and joint["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and joint["observed_odd_SEALED"] is True
          and joint["individual_LRG_ELG_halo_M_HOD_time_history_known"] is False,
          "original DM14 joint scope changed")
    g=list(geometry())
    check(len(g)==576,"independent original E16 original 576 angles")
    tol=p["preregistered_QA"]
    maxdiff=0.;n=0;reports={}
    ix=0
    for mass in MASS:
        for z in Z:
            for dex in DEX:
                filename="e17d2b0_original_DM14_conditional_snapshot_%02d.json"%ix
                raw=(original/filename).read_bytes()
                check(sha(raw)==joint["all_original_case_SHA"][filename]["sha256"],
                      "immutable full original DM14 case SHA "+filename)
                o=json.loads(raw)
                check(o["case_index"]==ix
                      and o["original_E16_geometry_count"]==576
                      and o["actual_short_short_long_values"]==1728
                      and o["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
                      and o["observed_odd_SEALED"] is True
                      and o["same_mass_across_z_is_NOT_one_halo_mass_accretion_history"] is True
                      and o["original_CLASS_and_DM14_h_equal"] is False,
                      "original case labels, size or hard physics stop")
                sn=snapshot(mass,z,dex)
                data=o["published_population_halo_snapshot"]
                for name,val in (("fit_a_of_z",sn["a"]),("fit_b_of_z",sn["b"]),
                                 ("c200c_population_median_DM14",sn["median"]),
                                 ("c200c_conditional_illustrative",sn["c"]),
                                 ("M200c_physical_Msun_DM14",sn["mass"]),
                                 ("H_DM14_km_s_Mpc",sn["H"]),
                                 ("rho_crit_DM14_Msun_Mpc3",sn["rho"]),
                                 ("r200c_physical_Mpc_DM14",sn["r200"]),
                                 ("r_s_physical_Mpc_DM14",sn["rs"]),
                                 ("r_s_comoving_Mpc_DM14",sn["rs_com"]),
                                 ("static_GM_over_rs_km2_s2_ONLY_not_actual_eBOSS_halo",sn["GMr"])):
                    observed=data[name]
                    maxdiff=max(maxdiff,abs(observed-val)/max(1.,abs(observed),abs(val)))
                check(data["z_original_CLASS_math_node_only"]==z
                      and data["M200c_hDM_inverse_Msun_QA_anchor"]==mass*1.e12
                      and data["illustrative_log10_concentration_offset_dex_NOT_z095_posterior"]==dex,
                      "mass/z/dex are only preregistered QA anchors")
                reclose=((4.*math.pi/3.)*200.*sn["rho"]*
                         sn["r200"]**3)/sn["mass"]
                check(abs(reclose-1.)<tol[
                    "r200_200critical_mass_reclosure_relative_gap_max"],
                    "independent radius does not enclose M200c")
                # Directly integrate F_NFW from mass density. No original SciPy/SiCi.
                cache={}
                def u(x):
                    key=round(x,13)
                    if key not in cache:
                        hi=nfw_direct_radial(x,sn["c"],4096)
                        lo=nfw_direct_radial(x,sn["c"],2048)
                        check(abs(hi-lo)/15.<tol[
                            "independent_2048_4096_radial_convergence_max"],
                            "independent 2048/4096 radial Fourier not converged")
                        cache[key]=(x,hi)
                    first,value=cache[key]
                    check(abs(first-x)<1e-12,"radial cache aliases distinct k")
                    return value
                rows=o["original_576_E16_geometry_and_two_legs_plus_long_QA"]
                check(len(rows)==576,"original 576 E16 geometry rows")
                for expected,row in zip(g,rows):
                    k,K,ms,ml,index,k1,k2=expected
                    check((row["k"],row["K"],row["mu_s"],row["mu_L"],
                           row["phi_index"])==(k,K,ms,ml,index),
                          "independently reconstructed E16 geometry labels")
                    for leg,km in enumerate((k1,k2)):
                        om=row["original_short_leg_moduli_h_Mpc"][leg]
                        maxdiff=max(maxdiff,abs(om-km)/max(1.,km))
                        stored=row["actual_two_short_legs_static_spatial_form_factors"][leg]
                        x=km*HC*sn["rs_com"]
                        val=u(x)
                        maxdiff=max(maxdiff,
                                    abs(stored["k_com_original_CLASS_h_per_Mpc"]-km)/
                                    max(1.,km),
                                    abs(stored["dimensionless_x_mixed_cosmology_QA_only"]-x)/
                                    max(1.,abs(x)),
                                    abs(stored["normalized_static_NFW_density_Fourier_u_only"]-val),
                                    abs(stored["dimensionless_Poisson_Phi_shape_minus_u_over_x2_ONLY"]+
                                        val/(x*x))/max(1.,abs(val/(x*x))))
                        n+=1
                    stored=row["original_long_K_static_spatial_form_factor"]
                    x=K*HC*sn["rs_com"]
                    val=u(x)
                    maxdiff=max(maxdiff,
                                abs(stored["dimensionless_x_mixed_cosmology_QA_only"]-x)/
                                max(1.,abs(x)),
                                abs(stored["normalized_static_NFW_density_Fourier_u_only"]-val),
                                abs(stored["dimensionless_Poisson_Phi_shape_minus_u_over_x2_ONLY"]+
                                    val/(x*x))/max(1.,abs(val/(x*x))))
                    n+=1
                reports[filename]={"sha256":sha(raw),
                                   "mass_ratio_to_published_pivot_QA_only":mass,
                                   "z_math_only":z,"dex_sensitivity_only":dex,
                                   "independent_unique_radial_integrals":len(cache),
                                   "replayed_frozen_E16_short_short_long":1728}
                print("E17D2B0_INDEPENDENT_CASE",ix,"N",1728,
                      "RADIAL_INTEGRALS",len(cache),"MAX_SCALED_GAP",maxdiff,
                      flush=True)
                ix+=1
    check(ix==18 and n==31104,"independent 18x576x3 original mode closure")
    check(maxdiff<tol["independent_stdlib_direct_radial_Simpson_full_31104_mode_gap_max"],
          "full independent DM14+NFW snapshot vs original analytic SiCi")
    return maxdiff,n,reports

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"independent full report differs")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b0i_",delete=False) as fh:
            t=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B0_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--output-dir",type=Path,default=OUT)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    args=a.parse_args();p=gate()
    check(len(list(geometry()))==576,"independent E16 analytic geometry")
    if args.self_test:
        sn=snapshot(1.,.95,0.)
        check(abs(nfw_direct_radial(0.,sn["c"],4096)-1.)<1e-11,
              "positive-mass independently radial NFW normalization")
        check(abs(sn["r200"]**3*(4.*math.pi/3.)*
                  200.*sn["rho"]/sn["mass"]-1.)<1e-12,
              "independent Dutton Planck 200crit unit reclosure")
        print("E17D2B0_INDEPENDENT_DM14_200CRIT_AND_ORIGINAL_SHA_GATE_PASS",
              "NO_HALO_TIME_OR_GALAXY_B",flush=True)
        return
    path=args.output_dir/"e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json"
    raw=path.read_bytes();joint=json.loads(raw)
    if args.negative_control:
        bad=json.loads(raw)
        bad["all_original_case_SHA"][
            "e17d2b0_original_DM14_conditional_snapshot_00.json"]["sha256"]="0"*64
        try:full_replay(args.output_dir,bad,p)
        except ValueError as exc:
            check("immutable full original DM14 case SHA" in str(exc),
                  "original report SHA tamper failed at unexpected gate")
            print("E17D2B0_INDEPENDENT_TAMPERED_ORIGINAL_SHA_REJECTED",flush=True)
            return
        raise AssertionError("tampered original SHA accepted")
    gap,n,reports=full_replay(args.output_dir,joint,p)
    result={"date":"2026-09-28",
            "status":"E17D2B0_INDEPENDENT_PUBLISHED_DM14_CONDITIONAL_18_POPULATION_SNAPSHOTS_FULL_ORIGINAL_31104_RADIAL_REPLAY_PASS_PHYSICAL_B_BLOCKED",
            "original_joint_sha256":sha(raw),
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_E16_full_sha256":p["immutable_source_parents"][
                "E16_original_full_sha256"],
            "all_original_case_sha_and_mode_counts":reports,
            "original_DM14_population_mass_anchor_redshift_concentration_cases":18,
            "full_original_E16_geometries_per_case":576,
            "independent_replayed_conditional_short_short_long_shapes":n,
            "max_scaled_DM14_200crit_radius_and_independent_NFW_radial_gap":gap,
            "method":"Pure stdlib independently evaluate primary published DM14 c200c fit and physical H/rho_crit200/radii in original DM14 Planck2013 background; independent Cartesian-free E16 analytic magnitudes; direct 4096-panel positive NFW radial Simpson for all 31104 original conditional modes; SHA tamper rejection. No original runner import, NumPy, SciPy, CLASS or observed galaxies.",
            "DM14_population_prior_is_NOT_a_single_halo_mass_assembly":True,
            "original_CLASS_and_DM14_cosmology_not_equal":True,
            "real_LRG_ELG_halo_mass_concentration_HOD_time_history_known":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"no_new_CLASS_FITS_mock":True,
            "main_untouched_PR_draft":True}
    save_once(args.output_dir/"e17d2b0_independent_stdlib_dm14_full_radial_replay.json",result)
    print("E17D2B0_INDEPENDENT_FULL_31104_CONDITIONAL_SPATIAL_SHAPES_PASS",
          "MAX",gap,"PHYSICAL_HALO_HISTORY_AND_B_BLOCKED",flush=True)
if __name__=="__main__":
    main()
