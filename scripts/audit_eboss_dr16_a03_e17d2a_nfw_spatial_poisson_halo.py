#!/usr/bin/env python3
"""E17D2a: published finite-mass spherical NFW spatial Poisson halo FAMILY.

Independently source a positive-density truncated NFW radial halo profile,
verify mass, real-space Poisson potential and analytic Si/Ci Fourier shapes
for BOTH original E16 short legs and original long K. All M_200, r_s and c
for actual eBOSS galaxies are UNKNOWN. The c/alpha grids are mathematical
QA only, not a halo abundance/HOD or physical retarded neutrino model.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
from scipy.special import sici

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2a_truncated_nfw_spatial_poisson_halo_family_prereg_2026-09-28.json"
PRE_BLOB="e7e0281f0b6b205214db5faf03a7c364521b3f2c"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
D0=ROOT/"source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27/e17d0_original_joint_external_potential_source_only.json"
D1DIR=ROOT/"source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27"
D1=D1DIR/"e17d1_original_joint_causal_history_nonidentifiability.json"
D1I=D1DIR/"e17d1_independent_full_finite_q_history_replay.json"
D1MAN=D1DIR/"archive_manifest.json"
D1PRE=ROOT/"source_data/eboss_dr16_a03_e17d1_causal_external_history_nonidentifiability_prereg_2026-09-27.json"
KS=(.05,.075,.1)
KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.)
ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi)
CS=(2,4,8)
ALPHAS=(.1,1.)
K0=.05
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2a_nfw_spatial"

def check(ok,msg):
    if not ok:raise ValueError("E17D2A_FAIL_CLOSED: "+msg)

def sha(b):return hashlib.sha256(b).hexdigest()

def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"existing source-only report differs: "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2a_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D2A_ORIGINAL_SPATIAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"original prospective NFW protocol Git blob drift")
    p=json.loads(PRE.read_bytes())
    parent=p["immutable_parents"]
    for path,key in ((E8,"E8_original_4000q_CSV_sha256"),
                     (E16,"E16_original_full_sha256"),
                     (D0,"E17D0_original_joint_sha256"),
                     (D1,"E17D1_original_joint_sha256"),
                     (D1I,"E17D1_independent_sha256")):
        check(sha(path.read_bytes())==parent[key],"original parent SHA missing: "+key)
    check(blob(D1MAN.read_bytes())==parent["E17D1_archive_manifest_git_blob"]
          and blob(D1PRE.read_bytes())==parent["E17D1_protocol_git_blob"],
          "E17D1 original SHA manifest/protocol Git blob drift")
    old=json.loads(D1.read_bytes())
    check(old["original_E16_geometries"]==576
          and old["original_state_leg_phase_samples"]==10368
          and old["full_physical_finiteK_bispectrum"]=="BLOCKED"
          and old["observed_odd_SEALED"] is True,
          "E17D1 frozen original physical STOP")
    fr=p["frozen_E16"]
    check(fr["original_short_k_comoving_h_per_Mpc"]==list(KS)
          and fr["original_long_K_comoving_h_per_Mpc"]==list(KLS)
          and fr["mu_short"]==list(MS) and fr["mu_long"]==list(ML)
          and fr["original_closed_geometries"]==576
          and fr["original_neutrino_states"]==["FD","plus","minus"]
          and fr["tracer_order"]==["LRG","ELG"],
          "frozen E16 state or geometry drift")
    m=p["physical_model_and_exact_equations"]
    check(m["fixed_numeric_QA_only_concentrations"]==list(CS)
          and m["fixed_numeric_QA_only_alpha_k0_rs_comoving"]==list(ALPHAS)
          and m["k0_comoving_h_per_Mpc"]==K0
          and all(p["absolute_STOP"].values()),
          "conditional halo QA parameters or science STOP changed")
    check(json.loads(E16.read_bytes())["QA"]["original_72_cases"]==72,
          "original E16 E14/E15 source counts drift")
    return p

def A(c):
    check(c>0,"positive NFW concentration")
    return math.log1p(c)-c/(1.+c)

def mass_ratio(y,c):
    check(y>=0 and c>0,"nonphysical radial mass")
    # The enclosed mass at the origin vanishes; A(0)=0 is a valid
    # radial limit and must not be passed through A's c>0 guard.
    if y==0.:
        return 0.
    return A(min(y,c))/A(c)

def potential_normalized(y,c):
    """Phi/(G*M/r_s), physical static isolated truncated NFW."""
    check(y>=0 and c>0,"nonphysical NFW potential radius")
    if y>=c:return -1./y
    core=1. if y==0. else math.log1p(y)/y
    return -(core-1./(1.+c))/A(c)

def radial_potential_gradient(y,c):
    """d[Phi/(G*M/r_s)]/d(r/r_s); positive, gravity points inward."""
    check(y>0 and c>0,"radial gradient at strictly positive radius")
    return mass_ratio(y,c)/(y*y)

def analytic_u(x,c):
    """Cooray-Sheth Si/Ci Fourier shape; x=k_phys*r_s_phys dimensionless."""
    check(c>0 and math.isfinite(x),"finite Fourier wavenumber and c")
    x=abs(x)
    if x==0.:return 1.
    si1,ci1=sici(x)
    si2,ci2=sici((1.+c)*x)
    out=(math.sin(x)*(si2-si1)+math.cos(x)*(ci2-ci1)
         -math.sin(c*x)/((1.+c)*x))/A(c)
    check(math.isfinite(out),"nonfinite analytic NFW Si/Ci Fourier shape")
    return float(out)

def geometry():
    rows=[]
    for k in KS:
        for K in KLS:
            for ms in MS:
                for ml in ML:
                    for i,phi in enumerate(PH):
                        u=(math.sqrt(max(0.,1.-ms*ms)),0.,ms)
                        v=(math.sqrt(max(0.,1.-ml*ml))*math.cos(phi),
                           math.sqrt(max(0.,1.-ml*ml))*math.sin(phi),ml)
                        L=tuple(K*item for item in v)
                        left=tuple(k*u[j]-.5*L[j] for j in range(3))
                        right=tuple(-k*u[j]-.5*L[j] for j in range(3))
                        k1=math.sqrt(sum(t*t for t in left))
                        k2=math.sqrt(sum(t*t for t in right))
                        closure=math.sqrt(sum((left[j]+right[j]+L[j])**2
                                              for j in range(3)))
                        check(closure<1e-13,"original E16 closed short legs")
                        rows.append({"k":k,"K":K,"mu_s":ms,"mu_L":ml,
                                     "phi_index":i,"phi_rad":phi,
                                     "original_short_leg_moduli_h_Mpc":[k1,k2],
                                     "original_long_K_h_Mpc":K,
                                     "original_triangle_closure_h_Mpc":closure})
    check(len(rows)==576,"E16 12×48 original geometry count")
    return rows

def radial_QA(c,p):
    ac=A(c)
    check(ac>0 and math.isfinite(ac),"positive A(c)")
    ys=[(0. if w==0 else c*w) for w in p[
        "physical_model_and_exact_equations"]["radial_QA_y_over_c"]]
    masses=[mass_ratio(y,c) for y in ys]
    check(all(0.<=m<=1.0000000000001 for m in masses)
          and all(masses[i]<=masses[i+1] for i in range(len(masses)-1))
          and masses[-1]==1.,"positive monotone finite halo mass")
    left_p=-(math.log1p(c)/c-1./(1.+c))/ac
    right_p=-1./c
    potgap=abs(left_p-right_p)
    gradgap=abs((A(c)/ac)/(c*c)-1./(c*c))
    check(potgap<p["preregistered_QA"][
        "static_potential_at_truncation_continuity_abs_max"]
        and gradgap<p["preregistered_QA"][
        "static_radial_force_at_truncation_continuity_abs_max"],
        "boundary potential/force continuity of truncated NFW")
    center=potential_normalized(0.,c)
    check(math.isfinite(center),"finite NFW center potential")
    for y in ys[1:]:
        check(math.isfinite(potential_normalized(y,c))
              and radial_potential_gradient(y,c)>=0.,"radial Poisson shape")
    B=.5*c*c-2.*c-1.+1./(1.+c)+3.*math.log1p(c)
    check(B>0,"second radial moment must be positive")
    x=p["physical_model_and_exact_equations"]["small_x_test"]
    tiny=analytic_u(x,c)
    approx=1.-x*x*B/(6.*ac)
    smallgap=abs(tiny-approx)
    check(smallgap<p["preregistered_QA"][
        "small_x_second_moment_abs_gap_max"],
        "NFW small-x exact finite-radial-moment QA")
    return {"c_QA_only":c,"A_c":ac,"mass_fractions_at_QA_radii":masses,
            "radial_y_QA_only":ys,"center_Phi_over_GM_div_rs":center,
            "max_boundary_Phi_continuity_abs":potgap,
            "max_boundary_gravitational_gradient_continuity_abs":gradgap,
            "B_c_second_radial_moment":B,
            "small_x_second_moment_abs_gap":smallgap}

def run_case(c,alpha,out):
    p=gate()
    check(c in CS and alpha in ALPHAS,"only preregistered dimensionless NFW QA")
    radial=radial_QA(c,p)
    maxabs=0.;maxreversal=0.;maxswap=0.;maxclosure=0.;records=[]
    def mode(k):
        x=alpha*k/K0
        u=analytic_u(x,c)
        return {"x_k_rs_dimensionless":x,
                "u_normalized_mass_Fourier":u,
                "Phi_Fourier_over_4pi_GM_rs2":-u/(x*x)}
    for geo in geometry():
        legs=[mode(k) for k in geo["original_short_leg_moduli_h_Mpc"]]
        long=mode(geo["original_long_K_h_Mpc"])
        for sample in legs+[long]:
            x=sample["x_k_rs_dimensionless"]
            u=sample["u_normalized_mass_Fourier"]
            maxabs=max(maxabs,abs(u))
            maxreversal=max(maxreversal,abs(u-analytic_u(-x,c)))
        maxclosure=max(maxclosure,geo["original_triangle_closure_h_Mpc"])
        geo["short_leg_spatial_shape_only"]=legs
        geo["long_mode_spatial_shape_only"]=long
        records.append(geo)
    checks=p["preregistered_QA"]
    check(maxabs<checks["u_positive_mass_characteristic_abs_upper_bound"]
          and maxreversal<checks["u_Fourier_reversal_abs_max"],
          "NFW real positive mass characteristic Fourier QA")
    lookup={(r["k"],r["K"],r["mu_s"],r["mu_L"],r["phi_index"]):r
            for r in records}
    for row in records:
        other=lookup.get((row["k"],row["K"],-row["mu_s"],row["mu_L"],
                          2-row["phi_index"]))
        if other is None:
            ms,ml,phi=row["mu_s"],row["mu_L"],row["phi_rad"]
            dot=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(
                max(0.,1.-ml*ml))*math.cos(phi)
            z=row["K"]/row["k"]
            reflected=(row["k"]*math.sqrt(1.+z*z/4.+z*dot),
                       row["k"]*math.sqrt(1.+z*z/4.-z*dot))
            for j in range(2):
                maxswap=max(maxswap,abs(row["original_short_leg_moduli_h_Mpc"][j]-
                                        reflected[1-j]))
                x=alpha*reflected[1-j]/K0
                maxswap=max(maxswap,abs(
                    row["short_leg_spatial_shape_only"][j]["u_normalized_mass_Fourier"]-
                    analytic_u(x,c)))
        else:
            for j in range(2):
                r= row["short_leg_spatial_shape_only"][j]
                t=other["short_leg_spatial_shape_only"][1-j]
                maxswap=max(maxswap,abs(r["u_normalized_mass_Fourier"]-
                                        t["u_normalized_mass_Fourier"]),
                             abs(r["x_k_rs_dimensionless"]-
                                 t["x_k_rs_dimensionless"]))
    check(maxswap<checks["same_halo_profile_short_leg_swap_abs_max"]
          and maxclosure<checks["original_E16_576_triangle_closure_abs_h_per_Mpc_max"],
          "original E16 both-leg geometry or model exchange")
    eps=checks["center_squeezed_geometric_probe_epsilon_K_over_k"]
    maxsqueezed=0.
    for k in KS:
        xc=alpha*k/K0
        ref=analytic_u(xc,c)
        for leg in (-.5,+.5):
            x=alpha*k*(1.+leg*eps)/K0
            maxsqueezed=max(maxsqueezed,abs(analytic_u(x,c)-ref))
    check(maxsqueezed<checks["center_squeezed_u_max_abs_difference"],
          "analytic finite-K->0 conditional halo spatial Fourier shape")
    result={"date":"2026-09-28",
            "status":"E17D2A_CONDITIONAL_FINITE_MASS_NFW_STATIC_POISSON_SPATIAL_HALO_FAMILY_ONLY_FULL_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_E16_sha256":p["immutable_parents"]["E16_original_full_sha256"],
            "original_E17D1_joint_sha256":p["immutable_parents"]["E17D1_original_joint_sha256"],
            "c_numeric_QA_only_NOT_fitted_to_eBOSS":c,
            "alpha_k0_rs_comoving_numeric_QA_only_NOT_real_halo_scale":alpha,
            "M200_actual_eBOSS_halo_mass_selected":False,
            "concentration_actual_eBOSS_halo_selected":False,
            "real_halo_formation_time_or_potential_history_selected":False,
            "radial_mass_potential_Poisson_QA":radial,
            "original_E16_geometries":len(records),
            "original_E16_exact_short_legs_per_geometry":2,
            "QA":{"max_abs_normalized_NFW_Fourier":maxabs,
                  "max_Fourier_reversal_abs":maxreversal,
                  "max_short_leg_exchange_abs":maxswap,
                  "max_original_triangle_closure_h_Mpc":maxclosure,
                  "max_center_squeezed_u_abs_gap":maxsqueezed},
            "original_E16_two_legs_and_long_static_spatial_Fourier_QA":records,
            "physical_halo_temporal_response_and_tracer_coupling":"BLOCKED",
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"no_new_CLASS_or_FITS_or_mock":True,
            "main_untouched_PR_draft":True}
    filename="e17d2a_original_nfw_c"+str(c)+"_alpha"+str(alpha).replace(".","p")+".json"
    save_once(out/filename,result)
    print("E17D2A_ORIGINAL_SPATIAL_CASE",c,alpha,"N",len(records),
          "MAX_U",maxabs,"SQUEEZED",maxsqueezed,
          "PHYSICAL_HALO_TIME_AND_GALAXY_B_BLOCKED",flush=True)

def aggregate(out):
    p=gate()
    reports={}
    for c in CS:
        for alpha in ALPHAS:
            name="e17d2a_original_nfw_c"+str(c)+"_alpha"+str(alpha).replace(".","p")+".json"
            raw=(out/name).read_bytes()
            row=json.loads(raw)
            check(row["c_numeric_QA_only_NOT_fitted_to_eBOSS"]==c
                  and row["alpha_k0_rs_comoving_numeric_QA_only_NOT_real_halo_scale"]==alpha
                  and row["original_E16_geometries"]==576
                  and row["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
                  and row["observed_odd_SEALED"] is True,
                  "original QA case report missing or physical STOP drift")
            reports[name]={"sha256":sha(raw),"size_bytes":len(raw),"QA":row["QA"],
                           "radial":row["radial_mass_potential_Poisson_QA"]}
    result={"date":"2026-09-28",
            "status":"E17D2A_THREE_C_TWO_ALPHA_PUBLISHED_STATIC_NFW_POISSON_SPATIAL_FAMILY_ORIGINAL_DONE_PHYSICAL_HALO_TIME_AND_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "immutable_parents":p["immutable_parents"],
            "original_numerical_QA_cases":len(reports),
            "original_E16_geometries_per_case":576,
            "actual_short_legs_plus_long_per_geometry":3,
            "original_case_reports":reports,
            "actual_M200_c_rs_eBOSS_halo_formation_or_HOD_calibrated":False,
            "time_dependent_retarded_Einstein_Vlasov_response_calculated":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(out/"e17d2a_original_joint_nfw_static_spatial_poisson_family.json",result)
    print("E17D2A_ORIGINAL_NFW_SPATIAL_FAMILY_ONLY_DONE",
          "CASES",len(reports),"PHYSICAL_HALO_TIME_AND_GALAXY_B_BLOCKED",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--case",nargs=2,metavar=("C","ALPHA"))
    g.add_argument("--aggregate",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    args=ap.parse_args()
    if args.self_test:
        p=gate()
        check(len(geometry())==576,"frozen original E16 triangles")
        check(all(radial_QA(c,p) for c in CS),"physical NFW radial control")
        print("E17D2A_PROSPECTIVE_NFW_ORIGINAL_SHA_AND_RADIAL_POISSON_GATE_OK",
              "SIX_DIMENSIONLESS_QA_CASES_NO_HALO_MASS",flush=True)
    elif args.case:
        c=int(args.case[0]);alpha=float(args.case[1])
        run_case(c,alpha,args.output_dir)
    else:aggregate(args.output_dir)

if __name__=="__main__":
    main()
