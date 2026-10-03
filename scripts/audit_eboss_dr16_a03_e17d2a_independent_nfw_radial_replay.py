#!/usr/bin/env python3
"""Independent E17D2a stdlib radial NFW integration and E16 scalar replay.

Uses direct 4096-panel Simpson mass-weighted Fourier quadrature, not original
SciPy Si/Ci formula, independently reconstructs every exact E16 triangle,
radial mass/potential/force boundary, and all original 6x576x3 mode shapes.
M_200/r_s/time, high-z tracer HOD and the physical galaxy B remain UNKNOWN.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2a_truncated_nfw_spatial_poisson_halo_family_prereg_2026-09-28.json"
PRE_BLOB="e7e0281f0b6b205214db5faf03a7c364521b3f2c"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
D0=ROOT/"source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27/e17d0_original_joint_external_potential_source_only.json"
D1DIR=ROOT/"source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27"
D1=D1DIR/"e17d1_original_joint_causal_history_nonidentifiability.json"
D1IND=D1DIR/"e17d1_independent_full_finite_q_history_replay.json"
D1MAN=D1DIR/"archive_manifest.json"
D1PRE=ROOT/"source_data/eboss_dr16_a03_e17d1_causal_external_history_nonidentifiability_prereg_2026-09-27.json"
KS=(.05,.075,.1);KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.);ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi)
CS=(2,4,8);ALPHAS=(.1,1.);K0=.05
DEFAULT=ROOT/"eboss_workspace/a03_physics_source/e17d2a_nfw_spatial"

def check(ok,why):
    if not ok:raise ValueError("E17D2A_INDEPENDENT_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"prospective E17D2a protocol Git blob")
    p=json.loads(PRE.read_bytes());par=p["immutable_parents"]
    for path,key in ((E8,"E8_original_4000q_CSV_sha256"),
                     (E16,"E16_original_full_sha256"),
                     (D0,"E17D0_original_joint_sha256"),
                     (D1,"E17D1_original_joint_sha256"),
                     (D1IND,"E17D1_independent_sha256")):
        check(sha(path.read_bytes())==par[key],"immutable original SHA: "+key)
    check(blob(D1MAN.read_bytes())==par["E17D1_archive_manifest_git_blob"]
          and blob(D1PRE.read_bytes())==par["E17D1_protocol_git_blob"],
          "original E17D1 Git manifest/protocol identity")
    check(json.loads(D1.read_bytes())["observed_odd_SEALED"] is True
          and json.loads(D1.read_bytes())["full_physical_finiteK_bispectrum"]=="BLOCKED",
          "original E17D1 physical science STOP")
    fr=p["frozen_E16"];m=p["physical_model_and_exact_equations"]
    check(fr["original_short_k_comoving_h_per_Mpc"]==list(KS)
          and fr["original_long_K_comoving_h_per_Mpc"]==list(KLS)
          and fr["mu_short"]==list(MS) and fr["mu_long"]==list(ML)
          and fr["original_closed_geometries"]==576
          and fr["original_neutrino_states"]==["FD","plus","minus"]
          and fr["tracer_order"]==["LRG","ELG"]
          and m["fixed_numeric_QA_only_concentrations"]==list(CS)
          and m["fixed_numeric_QA_only_alpha_k0_rs_comoving"]==list(ALPHAS)
          and m["k0_comoving_h_per_Mpc"]==K0
          and all(p["absolute_STOP"].values()),
          "frozen source, NFW QA or physics STOP drift")
    return p

def A(c):
    return math.log1p(c)-c/(1.+c)

def radial_mass(y,c):
    if y>=c:return 1.
    if y==0:return 0.
    return (math.log1p(y)-y/(1.+y))/A(c)

def radial_potential(y,c):
    if y>=c:return -1./y
    if y==0:return -(1.-1./(1.+c))/A(c)
    return -(math.log1p(y)/y-1./(1.+c))/A(c)

def radial_force_derivative(y,c):
    if y<=0:raise ValueError("strict positive radial derivative")
    return radial_mass(y,c)/(y*y)

def analytic_geometry():
    for k in KS:
        for K in KLS:
            r=K/k
            for ms in MS:
                for ml in ML:
                    for i,phi in enumerate(PH):
                        dot=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(
                            max(0.,1.-ml*ml))*math.cos(phi)
                        yield (k,K,ms,ml,i,
                               k*math.sqrt(1.+r*r/4.-r*dot),
                               k*math.sqrt(1.+r*r/4.+r*dot))

def density_mass_integrand(y,x):
    if y==0:return 0.
    v=x*y
    sinc=(1. if v==0. else math.sin(v)/v)
    return y/(1.+y)**2*sinc

def direct_simpson(x,c,n):
    """Direct radial integration, independent of SciPy's Si/Ci formula."""
    check(n%2==0 and n>=128,"even independent quadrature panel count")
    x=abs(x)
    step=c/n
    endpoint=density_mass_integrand(c,x)
    odd=[];even=[]
    for j in range(1,n):
        val=density_mass_integrand(step*j,x)
        if j%2:odd.append(val)
        else:even.append(val)
    total=step/3.*(endpoint+4.*math.fsum(odd)+2.*math.fsum(even))
    return total/A(c)

def independent_radial_qa(c,p):
    ac=A(c)
    check(ac>0 and math.isfinite(ac),"NFW positive finite mass normalizer")
    y_values=[0.,.25*c,.5*c,c,2.*c]
    mass=[radial_mass(y,c) for y in y_values]
    check(all(0.<=x<=1.000000000001 for x in mass)
          and all(mass[j]<=mass[j+1] for j in range(4))
          and mass[-1]==1.,"NFW mass monotonicity")
    inphi=-(math.log1p(c)/c-1./(1.+c))/ac
    outphi=-1./c
    d_in=(math.log1p(c)-c/(1.+c))/(ac*c*c)
    d_out=1./(c*c)
    check(abs(inphi-outphi)<p["preregistered_QA"][
        "static_potential_at_truncation_continuity_abs_max"]
          and abs(d_in-d_out)<p["preregistered_QA"][
        "static_radial_force_at_truncation_continuity_abs_max"],
        "independent NFW radial potential/force Poisson boundary")
    return {"A_c":ac,"radial_y":y_values,"mass_fractions":mass,
            "normalized_center_potential":radial_potential(0.,c),
            "boundary_potential_abs_gap":abs(inphi-outphi),
            "boundary_force_abs_gap":abs(d_in-d_out)}

def replay(original,joint,p):
    check(joint["original_numerical_QA_cases"]==6
          and joint["original_E16_geometries_per_case"]==576
          and joint["actual_short_legs_plus_long_per_geometry"]==3
          and joint["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and joint["observed_odd_SEALED"] is True,
          "original NFW joint report case scope or physical STOP")
    geo=list(analytic_geometry())
    check(len(geo)==576,"independent original E16 576 analytic geometry")
    maxdiff=0.;n=0;checks={}
    qa_tol=p["preregistered_QA"][
        "analytic_Si_Ci_vs_independent_stdlib_Simpson_scaled_gap_max"]
    for c in CS:
        rqa=independent_radial_qa(c,p)
        for alpha in ALPHAS:
            name="e17d2a_original_nfw_c"+str(c)+"_alpha"+str(alpha).replace(".","p")+".json"
            path=original/name
            raw=path.read_bytes()
            check(sha(raw)==joint["original_case_reports"][name]["sha256"],
                  "original full NFW case SHA changed: "+name)
            report=json.loads(raw)
            check(report["original_E16_geometries"]==576
                  and report["c_numeric_QA_only_NOT_fitted_to_eBOSS"]==c
                  and report["alpha_k0_rs_comoving_numeric_QA_only_NOT_real_halo_scale"]==alpha
                  and report["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
                  and report["observed_odd_SEALED"] is True,
                  "original NFW science/QA drift: "+name)
            rq=report["radial_mass_potential_Poisson_QA"]
            maxdiff=max(maxdiff,abs(rqa["A_c"]-rq["A_c"]))
            for x,y in zip(rqa["mass_fractions"],rq["mass_fractions_at_QA_radii"]):
                maxdiff=max(maxdiff,abs(x-y))
            maxdiff=max(maxdiff,abs(rqa["normalized_center_potential"]-
                                    rq["center_Phi_over_GM_div_rs"]))
            cache={}
            def u(x):
                key=round(x,13)
                if key not in cache:
                    high=direct_simpson(x,c,4096)
                    low=direct_simpson(x,c,2048)
                    check(abs(high-low)/15.<qa_tol/2.,
                          "independent 2048/4096 Simpson convergence: "+name)
                    cache[key]=(x,high)
                first_x,val=cache[key]
                check(abs(first_x-x)<1e-12,"rounded cache changes physical x")
                return val
            data=report["original_E16_two_legs_and_long_static_spatial_Fourier_QA"]
            check(len(data)==576,"original E16 length")
            for expected,row in zip(geo,data):
                k,K,ms,ml,i,k1,k2=expected
                check((row["k"],row["K"],row["mu_s"],row["mu_L"],row["phi_index"])==
                      (k,K,ms,ml,i),"independent frozen geometry labels")
                for j,mod in enumerate((k1,k2)):
                    original_mod=row["original_short_leg_moduli_h_Mpc"][j]
                    maxdiff=max(maxdiff,abs(mod-original_mod))
                    r=row["short_leg_spatial_shape_only"][j]
                    x=alpha*mod/K0
                    check(x>0,"strict positive original E16 short leg")
                    val=u(x)
                    maxdiff=max(maxdiff,
                                abs(x-r["x_k_rs_dimensionless"])/max(1.,x),
                                abs(val-r["u_normalized_mass_Fourier"]),
                                abs(-val/(x*x)-
                                    r["Phi_Fourier_over_4pi_GM_rs2"])/
                                    max(1.,abs(val/(x*x))))
                    n+=1
                # Original real long K is shared across orientations; no 3-point B.
                r=row["long_mode_spatial_shape_only"]
                x=alpha*K/K0
                val=u(x)
                maxdiff=max(maxdiff,
                            abs(x-r["x_k_rs_dimensionless"])/max(1.,x),
                            abs(val-r["u_normalized_mass_Fourier"]),
                            abs(-val/(x*x)-
                                r["Phi_Fourier_over_4pi_GM_rs2"])/
                                max(1.,abs(val/(x*x))))
                n+=1
            checks[name]={"original_sha256":sha(raw),
                          "replayed_original_geometries":len(data),
                          "replayed_two_short_plus_long_modes":len(data)*3,
                          "independent_unique_radial_Fourier_integrals":len(cache)}
            print("E17D2A_INDEPENDENT_CASE",name,"REPLAYED",len(data)*3,
                  "UNIQUE_DIRECT_RADIAL_INTEGRALS",len(cache),
                  "MAX_GAP_SO_FAR",maxdiff,flush=True)
    check(n==6*576*3,"original NFW full modes replay count")
    check(maxdiff<qa_tol,"analytic Si/Ci vs independent direct radial Fourier")
    return maxdiff,n,checks

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"independent NFW audit certificate changed")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2ai_",delete=False) as fh:
            temp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17D2A_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--output-dir",type=Path,default=DEFAULT)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    args=a.parse_args()
    p=gate()
    check(len(list(analytic_geometry()))==576,"independent E16 angular/leg closure")
    if args.self_test:
        for c in CS:
            independent_radial_qa(c,p)
        check(abs(direct_simpson(0.,4,4096)-1.)<1e-11,
              "positive-mass NFW radial Fourier normalization")
        print("E17D2A_INDEPENDENT_ORIGINAL_SHA_SPATIAL_POISSON_PREFLIGHT_PASS",
              "NO_MASS_NO_HALO_TIME_NO_B",flush=True)
        return
    path=args.output_dir/"e17d2a_original_joint_nfw_static_spatial_poisson_family.json"
    raw=path.read_bytes()
    joint=json.loads(raw)
    if args.negative_control:
        bad=json.loads(raw)
        name="e17d2a_original_nfw_c2_alpha0p1.json"
        bad["original_case_reports"][name]["sha256"]="0"*64
        try:replay(args.output_dir,bad,p)
        except ValueError as exc:
            check("original full NFW case SHA changed" in str(exc),
                  "negative SHA control failed at unexpected gate")
            print("E17D2A_INDEPENDENT_MUTATED_SHA_REJECTED",flush=True)
            return
        raise AssertionError("tampered NFW report SHA was accepted")
    gap,n,cases=replay(args.output_dir,joint,p)
    result={"date":"2026-09-28",
            "status":"E17D2A_INDEPENDENT_STDLIB_ALL_FROZEN_NFW_SPATIAL_POISSON_SHAPES_PASS_PHYSICAL_TIME_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_joint_sha256":sha(raw),
            "original_E16_full_sha256":p["immutable_parents"]["E16_original_full_sha256"],
            "original_case_SHA_and_replay_counts":cases,
            "original_E16_geometries_per_case":576,
            "preregistered_c_alpha_QA_cases":6,
            "independently_replayed_actual_two_short_plus_long_modes":n,
            "max_scaled_analytic_SiCi_vs_direct_radial_Simpson_gap":gap,
            "method":"Independent pure-stdlib frozen E16 analytic triangle and 4096-panel mass-weighted direct radial Simpson against original analytic SciPy Si/Ci NFW transform, independently verified radial mass and potential/force boundary. No original NFW runner, SciPy, NumPy, halo mass/HOD or observed eBOSS input.",
            "real_halo_M200_c_rs_time_or_HOD_calibrated":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(args.output_dir/"e17d2a_independent_stdlib_direct_radial_fourier_replay.json",result)
    print("E17D2A_INDEPENDENT_10368_STATIC_SPATIAL_SHAPES_PASS",
          "MAX",gap,"FULL_PHYSICAL_B_BLOCKED",flush=True)

if __name__=="__main__":
    main()
