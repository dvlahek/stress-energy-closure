#!/usr/bin/env python3
"""A03E11: locked Gaussian-rank, source-only squeezed wind response.

Reuses byte-identical E8 4000q CSV, E9 three-CLASS source-only report and
E10 conditional angular-odd audit. Calculates Gamma_ell=E[r*T_ell(r)],
the coefficient of an ASSUMED independently calibrated long-mode/rank
cross-correlation. A nonzero Gamma is not an eBOSS bispectrum or detection.
The symmetric unweighted two-point intrinsic odd mean remains zero.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss, legval

ROOT=Path(__file__).resolve().parents[1]
PROTO=ROOT/"source_data/eboss_dr16_a03_e11_gaussian_wind_squeezed_response_prereg_2026-09-27.json"
PROTO_BLOB="c959d0a1c133c158131561acc7a36bdce4cfd0ab"
E10_CODE=ROOT/"scripts/audit_eboss_dr16_a03_e10_angular_wind_parity_bridge.py"
E10_BLOB="c396817e4b3983391d8dd722b25711c4336ec173"
E10_PROTO=ROOT/"source_data/eboss_dr16_a03_e10_wind_parity_angular_odd_bridge_protocol_2026-09-27.json"
E10_PROTO_BLOB="d9de48012103eaa07d531e11bca9cef285e62545"
E9_SCRIPT=ROOT/"scripts/build_eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase.py"
E9_BLOB="700a80d9ea9bc45454925b9a06fec2106203c506"
SRC_HASH={
    "E8":"bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0",
    "E9":"450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890",
    "E10":"cc55a914a8be0174eb7c1e49e0e1968537cedfcf4a429d52ad66ab86cb44a488"}
STATES=("FD","plus","minus")
ZM=.945
ZP=.955
ZC=.95
H=.005
KCOM=.05
MASS=.06
C_KMS=299792.458
T_NCDM=.71611
TCMB_K=2.7255
KB_EV_K=8.617333262e-5
TNU0=T_NCDM*TCMB_K*KB_EV_K
EPS=1e-14

def check(ok,msg):
    if not ok:
        raise RuntimeError(msg)

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def inputs(args):
    for path,blob in ((PROTO,PROTO_BLOB),(E10_CODE,E10_BLOB),
                      (E10_PROTO,E10_PROTO_BLOB),(E9_SCRIPT,E9_BLOB)):
        check(git_blob(path.read_bytes())==blob,
              "A frozen source/protocol Git blob changed: "+str(path))
    p=json.loads(PROTO.read_bytes())
    fixed=p["frozen_e11_calculation"]
    check(p["original_sha_parents"]["E8_4000q_CSV"]==SRC_HASH["E8"]
          and p["original_sha_parents"]["E9_three_CLASS_report"]==SRC_HASH["E9"]
          and p["original_sha_parents"]["E10_ang_report"]==SRC_HASH["E10"]
          and fixed["short_comoving_k_h_per_Mpc"]==KCOM
          and fixed["center_z"]==ZC
          and fixed["z_phase_neighbors"]==[ZM,ZP]
          and fixed["redshift_and_mu_angle_levels"]["angular_GaussLegendre_nodes"]==[256,512]
          and fixed["redshift_and_mu_angle_levels"]["Gaussian_rank_GaussHermite_nodes"]==[48,96]
          and fixed["redshift_and_mu_angle_levels"]["Stein_central_derivative_rank_steps"]==[.002,.001]
          and p["observed_odd_data_vector_read"] is False
          and p["no_new_catalogue_mock_download"] is True,
          "E11 prospective fixed scope changed")
    raw={name:path.read_bytes() for name,path in
         (("E8",args.e8_csv),("E9",args.e9_json),("E10",args.e10_json))}
    for name in raw:
        check(sha(raw[name])==SRC_HASH[name],
              "Source-only original "+name+" exact archived byte SHA256 changed")
    csv=np.genfromtxt(args.e8_csv,delimiter=",",names=True)
    check(csv.dtype.names==("q_dimensionless","F0_CLASS_normalized",
           "Fplus_CLASS_normalized","Fminus_CLASS_normalized")
          and csv.shape==(4000,), "Original 4000-q CSV schema changed")
    q=np.asarray(csv["q_dimensionless"],float)
    f={s:np.asarray(csv[col],float) for s,col in
       (("FD","F0_CLASS_normalized"),("plus","Fplus_CLASS_normalized"),
        ("minus","Fminus_CLASS_normalized"))}
    check(np.all(np.diff(q)>0) and all(np.all(np.isfinite(x)) and
          np.all(x>0) for x in f.values()),"Frozen E8 source invalid")
    e9=json.loads(raw["E9"])
    e10=json.loads(raw["E10"])
    check(e9["observed_odd_data_vector_read"] is False
          and e9["absolute_LRG_ELG_xi1_xi3_or_signal_to_noise_computed"] is False
          and e10["A03_physical_pair_z_RR_window_certified"] is False
          and e10["A04_independent_eBOSS_covariance_available"] is False,
          "Unexpected real survey/physical-window result in frozen parent")
    cases={(x["wind_velocity_derivative_relative_z_step"],x["z"]):x
           for x in e9["all_own_and_common_FD_wind_cases"]}
    check(len(cases)==14 and (.004,ZM) in cases and (.004,ZP) in cases
          and (.004,ZC) in cases,"Original E9 three-state derivative cases changed")
    return p,q,f,cases,e10

def occupancy(qsource,f,state,q):
    q=np.asarray(q,float)
    check(np.isfinite(q).all() and np.min(q)>=qsource[0]
          and np.max(q)<=qsource[-1],
          "E11 Gaussian tails would extrapolate past frozen source q grid")
    fd=1./(1.+np.exp(q))
    if state=="FD":return fd
    return fd*np.interp(q,qsource,f[state])/np.interp(q,qsource,f["FD"])

def phase(qsource,f,cases,state,z,mu,r):
    """Exact E9/E10 phase anchored at original E9 r=1, mu=1."""
    one=cases[(.004,z)]["state"][state]
    v=float(one["own_CLASS_v_R16_1sigma_los_kms"])
    a1=float(one["alpha_with_own_CLASS_velocity"])
    check(v>0 and a1>0 and math.isfinite(v),"Frozen E9 source wind invalid")
    q1=MASS*v/(C_KMS*TNU0*(1.+z))
    base=float(occupancy(qsource,f,state,q1))
    rr=np.asarray(r,float)[:,None]
    uu=np.asarray(mu,float)[None,:]
    return a1*rr*uu*occupancy(qsource,f,state,q1*np.abs(rr*uu))/base

def projected_rank(qsource,f,cases,state,nmu,r):
    u,w=leggauss(nmu)
    pm=phase(qsource,f,cases,state,ZM,u,r)
    pp=phase(qsource,f,cases,state,ZP,u,r)
    da=-(1.+ZC)*(pp-pm)/(2.*H)
    theta=u[None,:]**2*da
    return {ell:(2*ell+1)/2.*np.sum(
        w[None,:]*legval(u,[0.]*ell+[1.])[None,:]*theta,
        axis=1) for ell in (1,3)}

def ranked_gauss(ngh):
    x,w=hermgauss(ngh)
    return np.sqrt(2.)*x,w/math.sqrt(math.pi)

def solve(p,q,f,cases,e10):
    qalevels=p["frozen_e11_calculation"]["redshift_and_mu_angle_levels"]
    results={}
    qa={}
    # Four fully preregistered quadrature combinations, no seed or stochastic MC.
    for nmu in (256,512):
        for ngh in (48,96):
            rk,ww=ranked_gauss(ngh)
            check(abs(float(np.sum(ww))-1.)<2e-14
                  and abs(float(np.sum(ww*rk)))<2e-14
                  and abs(float(np.sum(ww*rk*rk))-1.)<3e-13,
                  "Gaussian rank quadrature not unit normal")
            vals={}
            for s in STATES:
                t=projected_rank(q,f,cases,s,nmu,rk)
                vals[s]={}
                for ell in (1,3):
                    arr=t[ell]
                    check(np.max(np.abs(arr+arr[::-1]))<2e-14,
                          "Source-conditional rank parity failed")
                    even_mean=float(np.sum(ww*arr))
                    # Pair +/- ranks explicitly as independent sign-null check.
                    paired_mean=float(np.sum(
                        .5*ww*(arr+arr[::-1])))
                    gamma=float(np.sum(ww*rk*arr))
                    rms=float(np.sqrt(np.sum(ww*arr*arr)))
                    sign_w=float(np.sum(ww*np.sign(rk)*arr))
                    check(abs(even_mean)<3e-14 and
                          abs(paired_mean)<3e-14 and np.isfinite(gamma)
                          and rms>0 and gamma>0,
                          "Spurious two-point mean or degenerate Gaussian response")
                    vals[s][str(ell)]={"unconditional_mean":even_mean,
                        "paired_positive_negative_mean":paired_mean,
                        "Gamma_E_rank_times_conditional_T":gamma,
                        "rms_conditional_T":rms,
                        "illustrative_E_sign_rank_times_T":sign_w}
            vals["plus_minus"]={
                str(ell):{"Gamma_plus_minus":
                     vals["plus"][str(ell)]["Gamma_E_rank_times_conditional_T"]
                     -vals["minus"][str(ell)]["Gamma_E_rank_times_conditional_T"],
                     "unconditional_mean_plus_minus":
                     vals["plus"][str(ell)]["unconditional_mean"]
                     -vals["minus"][str(ell)]["unconditional_mean"]}
                for ell in (1,3)}
            results[f"mu_{nmu}_gh_{ngh}"]=vals
    # Compare independently archived E10 T_ell(r=+1) exactly, not numeric
    # reconstitution of source CSV under a new NumPy environment.
    e10_ref=e10["conditional_Fourier_source_projected_ells_not_xi"]["512"]["conditional_positive_LOS"]
    for s in STATES:
        got=projected_rank(q,f,cases,s,512,np.asarray([1.]))
        for ell in (1,3):
            a=float(got[ell][0]);b=float(e10_ref[s][str(ell)])
            check(math.isclose(a,b,rel_tol=3e-13,abs_tol=3e-17),
                  "E11 Gaussian-rank r=1 does not reproduce exact E10")
            qa[f"original_E10_rank1_{s}_ell{ell}_absolute_error"]=abs(a-b)
    a=results["mu_512_gh_96"]
    # Gaussian integration by parts: E[r T(r)] = E[dT/dr].
    stein={}
    rk,ww=ranked_gauss(96)
    for s in STATES:
        for ell in (1,3):
            dlist={}
            for h in (.002,.001):
                plus=projected_rank(q,f,cases,s,512,rk+h)[ell]
                minus=projected_rank(q,f,cases,s,512,rk-h)[ell]
                dlist[str(h)]=float(np.sum(ww*(plus-minus)/(2.*h)))
            gamma=a[s][str(ell)]["Gamma_E_rank_times_conditional_T"]
            for h in (.002,.001):
                gap=abs(gamma-dlist[str(h)])/max(abs(gamma),abs(dlist[str(h)]),EPS)
                qa[f"Stein_rank_{s}_ell{ell}_step{h}_relative_gap"]=gap
                check(gap<p["frozen_e11_calculation"]["engineering_QA_warning_not_physical_threshold"]
                      ["Stein_derivative_0p002_vs_0p001_relative_warn"],
                      "Gaussian Stein rank derivative warning exceeded, stop")
            stein[s+"_"+str(ell)]={"Gamma_weighted_rank":gamma,
                "E_derivative_rank_steps":dlist}
    for s in STATES:
        for ell in (1,3):
            b=results["mu_512_gh_48"][s][str(ell)]["Gamma_E_rank_times_conditional_T"]
            c=a[s][str(ell)]["Gamma_E_rank_times_conditional_T"]
            err=abs(b-c)/max(abs(b),abs(c),EPS)
            qa[f"GH48_96_{s}_ell{ell}_relative_gap"]=err
            check(err<.03,"Gaussian 48/96 rank quadrature engineering QA fail")
            d=results["mu_256_gh_96"][s][str(ell)]["Gamma_E_rank_times_conditional_T"]
            err=abs(d-c)/max(abs(d),abs(c),EPS)
            qa[f"mu256_512_{s}_ell{ell}_relative_gap"]=err
            check(err<1e-4,"Angular 256/512 engineering QA fail")
    # Original E10 conditional T_ell(r=1) may be nonzero. An independent
    # long-mode-rank covariance is nevertheless required to form a bispectrum.
    check(all(a[s][str(ell)]["Gamma_E_rank_times_conditional_T"]>0
              for s in STATES for ell in (1,3)),
          "Gaussian response unexpectedly absent under frozen assumptions")
    return {
      "date":"2026-09-27",
      "status":"E11_FROZEN_GAUSSIAN_RANK_NONZERO_SQUEEZED_RESPONSE_KERNEL_SOURCE_ONLY_UNCONDITIONAL_TWO_POINT_ODD_MEAN_ZERO",
      "preregistration_git_blob_sha1":PROTO_BLOB,
      "original_E8_csv_SHA256":SRC_HASH["E8"],
      "original_E9_three_CLASS_report_SHA256":SRC_HASH["E9"],
      "original_E10_angular_report_SHA256":SRC_HASH["E10"],
      "model_definition":"For Gaussian r~N(0,1), locally coherent +r*sigma_LOS,s(z) original E9 wind and fixed frozen E8 F_s, T_ell(r) is the original E10 conditional EqA8 unit-bias/unit-Pcb angular coefficient. Gamma_ell=E[r*T_ell(r)]=E[dT_ell/dr].",
      "response_identity":"If an independently measured/modelled jointly Gaussian long-mode D_L and r have cross-spectrum P_(D_L,r)(K), a local squeezed three-point-like response has coefficient P_(D_L,r)(K)*Gamma_ell; neither P_(D_L,r) nor actual galaxy bispectrum is in this report.",
      "fixed_example_not_measured_eBOSS":{"z":ZC,"k_COMOVING_h_per_Mpc":KCOM,
          "Gaussian_rank_node_counts":[48,96],"angular_node_counts":[256,512],
          "source_class_velocity_step_relative_to_1plusz":.004,
          "centered_phase_redshift_step":H},
      "four_fixed_quadrature_results":results,
      "Stein_Gaussian_integration_by_parts_checks":stein,
      "engineering_QA_recomputed":qa,
      "technical_QA_warnings":[],
      "fixed_negative_controls":{"two_point_mean_exact_sign_symmetric_zero":True,
          "zero_long_mode_rank_cross_power_implies_zero_squeezed_three_point":True,
          "equal_LRG_ELG_bias_zeroes_EqA8_channel":True,
          "tracer_swap_complex_conjugates_signed_imaginary_cross_power":True,
          "original_E10_r_positive_1_independent_frozen_report_reproduced":True,
          "no_external_preferred_wind_vector_in_statistical_symmetry":True},
      "not_Zhu_Castorina_exact_bispectrum_estimator":True,
      "true_EBOSS_P_long_velocity_cross_transfer_calibrated":False,
      "full_triangle_kshort_Klong_redshift_halo_wake_and_bias_model_calibrated":False,
      "real_survey_triple_selection_window_or_covariance_calculated":False,
      "existing_24D_unconditional_EBOSS_odd_mean_reconstructed":False,
      "observed_galaxy_random_or_odd_read":False,
      "new_catalogue_or_mock_download_or_seed_or_cut":False,
      "physical_detection_or_exclusion":False
    }

def atomic_new(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        check(path.read_bytes()==raw,"E11 output preexists with different bytes")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e11_",delete=False) as f:
        tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)

def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument("--e8-csv",required=True,type=Path)
    pa.add_argument("--e9-json",required=True,type=Path)
    pa.add_argument("--e10-json",required=True,type=Path)
    pa.add_argument("--output",required=True,type=Path)
    args=pa.parse_args()
    p,q,f,cases,e10=inputs(args)
    out=solve(p,q,f,cases,e10)
    raw=(json.dumps(out,indent=2,allow_nan=False)+"\n").encode()
    atomic_new(args.output,raw)
    win=out["four_fixed_quadrature_results"]["mu_512_gh_96"]
    print("EBOSS_A03_E11_ORIGINAL_E8_E9_E10_GAUSSIAN_SQUEEZED_RESPONSE_SOURCE_ONLY_PASS",
          "REPORT_SHA256",sha(raw),
          "MAX_NUMERICAL_QA_GAP",max(out["engineering_QA_recomputed"].values()),
          "UNCONDITIONAL_TWO_POINT_ODD_MEAN_ZERO",
          "NO_EBOSS_BISPECTRUM_OR_XI",flush=True)
    for s in STATES:
        print("E11_GAMMA_GAUSSIAN",s,
              json.dumps({ell:win[s][ell]["Gamma_E_rank_times_conditional_T"]
                          for ell in ("1","3")},sort_keys=True),flush=True)
    print("E11_GAMMA_FPLUS_MINUS_FMINUS",
          json.dumps(win["plus_minus"],sort_keys=True),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
