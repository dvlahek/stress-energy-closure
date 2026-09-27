#!/usr/bin/env python3
"""E11 source-only Gaussian/Wick squeezed response of frozen E8/E9/E10.

Compute exact Gaussian-rank averages of ORIGINAL E10 conditional angular
Fourier odd coefficients. The unconditional odd two-point mean cancels for a
symmetric independent wind; a THREE-POINT response to a long mode correlated
with the relative wind is proportional to E[r*T_ell(r)].

This is a local squeezed-response kernel, not an eBOSS bispectrum measurement,
3D full-k bispectrum, halo velocity reconstruction or detection. Observed odd,
galaxies, randoms and survey selection are not opened.
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

import audit_eboss_dr16_a03_e10_angular_wind_parity_bridge as e10

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e11_gaussian_wind_squeezed_response_prereg_2026-09-27.json"
PROTOCOL_BLOB="c959d0a1c133c158131561acc7a36bdce4cfd0ab"
E10_ORIGINAL_SHA="cc55a914a8be0174eb7c1e49e0e1968537cedfcf4a429d52ad66ab86cb44a488"
ZS=(.945,.955)
STATES=("FD","plus","minus")
KCOM=.05
ZCENTER=.95
HSTEP=.005

def require(x,why):
    if not x:raise ValueError(why)

def protocol_gate():
    require(e10.e8.git_blob(PROTOCOL.read_bytes())==PROTOCOL_BLOB,
            "Original E11 precomputed-source protocol Git blob mismatch")
    p=json.loads(PROTOCOL.read_bytes())
    g=p["frozen_e11_calculation"]
    require(g["short_comoving_k_h_per_Mpc"]==KCOM and
            g["center_z"]==ZCENTER and
            g["z_phase_neighbors"]==list(ZS) and
            g["redshift_and_mu_angle_levels"]["angular_GaussLegendre_nodes"]==[256,512] and
            g["redshift_and_mu_angle_levels"]["Gaussian_rank_GaussHermite_nodes"]==[48,96] and
            g["redshift_and_mu_angle_levels"]["Stein_central_derivative_rank_steps"]==[.002,.001] and
            p["original_sha_parents"]["E10_ang_report"]==E10_ORIGINAL_SHA and
            p["observed_odd_data_vector_read"] is False and
            p["no_new_catalogue_mock_download"] is True and
            p["original_PR_remains_draft"] is True,
            "Pre-registered E11 analytical scope, input or observed odd guard changed")
    return p

def source_gate(e8_csv,e9_json,e10_json):
    proto=protocol_gate()
    _,q,fs,cases=e10.data_gate(e8_csv,e9_json)
    require(e10.e8.sha(e10_json.read_bytes())==E10_ORIGINAL_SHA,
            "Original E10 report's immutable JSON SHA mismatch")
    j=json.loads(e10_json.read_bytes())
    require(j["status"].startswith("E10_CONDITIONAL_ANGULAR_FOURIER_SOURCE_TESTS_PASS_") and
            j["cannot_directly_propagate_E9_positive_LOS_into_unconditional_eBOSS_24D"] is True and
            j["observed_galaxies_or_randoms_or_odd_read"] is False,
            "Original E10 source report contradicts current physical estimand gate")
    return proto,q,fs,cases,j

def phase_for_rank(q,fs,case,state,mu,rank):
    """Angular alpha using ORIGINAL E9 positive-rank alpha as sole scale."""
    require(state in STATES and abs(float(case["z"])-float(case["z"]))<1e-10,
            "Unknown E9 state/redshift")
    z=float(case["z"])
    own=case["state"][state]
    v1=float(own["own_CLASS_v_R16_1sigma_los_kms"])
    alpha1=float(own["alpha_with_own_CLASS_velocity"])
    q1=e10.e8.MASS_EV*v1/(e10.e8.C_KMS*e10.e8.cro.T_NCDM*
                  e10.e8.cro.TCMB_K*e10.e8.cro.KB_EV_K*(1.+z))
    require(v1>0 and alpha1>0 and q1>0,"Original E9 positive wind reference invalid")
    mu=np.asarray(mu,dtype="f8")
    rank=np.asarray(rank,dtype="f8")
    fqref=float(e10.occupation(q,fs,state,np.asarray(q1)))
    qres=q1*np.abs(rank*mu)
    require(np.max(qres)<=q[-1] and np.min(qres)>=q[0],
            "Original frozen q grid does not cover Gauss-Hermite wind rank")
    return alpha1*rank*mu*e10.occupation(q,fs,state,qres)/fqref

def transfer_Tells(q,fs,cases,state,ranks,nmu):
    """E10 angular exact-occupancy mode dipole and octupole for signed ranks."""
    mu,w=leggauss(nmu)
    r=np.asarray(ranks,dtype="f8").reshape(-1,1)
    uu=mu.reshape(1,-1)
    phlo=phase_for_rank(q,fs,cases[(e10.VELSTEP,ZS[0])],state,uu,r)
    phhi=phase_for_rank(q,fs,cases[(e10.VELSTEP,ZS[1])],state,uu,r)
    theta=uu**2 * (-(1.+ZCENTER)*(phhi-phlo)/(2.*HSTEP))
    output={}
    for ell in (1,3):
        coef=np.zeros(ell+1);coef[-1]=1.
        output[str(ell)]=(2.*ell+1.)/2. *np.sum(
              theta*(w*legval(mu,coef))[None,:],axis=1)
    require(all(np.isfinite(vals).all() for vals in output.values()),
            "Nonfinite Gaussian E11 conditional E10 projection")
    return output

def expected_gamma(q,fs,cases,state,ghn,nmu):
    x,hw=hermgauss(ghn)
    ranks=np.sqrt(2.)*x
    weights=hw/np.sqrt(np.pi)
    ts=transfer_Tells(q,fs,cases,state,ranks,nmu)
    return {ell:{"unconditional_mean":float(np.dot(weights,T)),
                 "conditional_fixed_rank_gaussian_RMS":float(
                      np.sqrt(np.dot(weights,T*T))),
                 "Gaussian_squeezed_response_E_rank_times_T":float(
                      np.dot(weights,ranks*T))}
                 for ell,T in ts.items()},(ranks,weights,ts)

def stein_gaussian(q,fs,cases,state,ghn,nmu,eps):
    x,hw=hermgauss(ghn)
    ranks=np.sqrt(2.)*x
    weights=hw/np.sqrt(np.pi)
    tp=transfer_Tells(q,fs,cases,state,ranks+eps,nmu)
    tm=transfer_Tells(q,fs,cases,state,ranks-eps,nmu)
    return {ell:float(np.dot(weights,(tp[ell]-tm[ell])/(2.*eps)))
            for ell in ("1","3")}

def output_new(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,"E11 output exists with different source bytes")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e11_",delete=False) as f:
        t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(t,path)
    finally:t.unlink(missing_ok=True)

def analyze(e8_csv,e9_json,e10_json):
    proto,q,fs,cases,original=source_gate(e8_csv,e9_json,e10_json)
    g=proto["frozen_e11_calculation"]
    levels=g["redshift_and_mu_angle_levels"]
    angular=[256,512]
    hermite=[48,96]
    results={}
    qa={}
    warnings=[]
    for nmu in angular:
        for ng in hermite:
            tag=f"mu_{nmu}_GH_{ng}"
            results[tag]={}
            for state in STATES:
                values,_=expected_gamma(q,fs,cases,state,ng,nmu)
                results[tag][state]=values
                for ell in ("1","3"):
                    avg=values[ell]["unconditional_mean"]
                    require(abs(avg)<2e-15,
                            "Symmetric Gaussian wind produced false mean two-point odd")
                    gamma=values[ell]["Gaussian_squeezed_response_E_rank_times_T"]
                    rms=values[ell]["conditional_fixed_rank_gaussian_RMS"]
                    require(math.isfinite(gamma) and rms>0.,
                            "Nonfinite/zero previously fixed E9 wind squeezed response")
            results[tag]["delta_plus_minus"]={
                ell:results[tag]["plus"][ell]["Gaussian_squeezed_response_E_rank_times_T"]-
                    results[tag]["minus"][ell]["Gaussian_squeezed_response_E_rank_times_T"]
                for ell in ("1","3")}
    # E10 exact rank-one unit (bL-bE)P source closes before any ensemble mean.
    t1={s:transfer_Tells(q,fs,cases,s,np.asarray([1.]),512)
        for s in STATES}
    for state in STATES:
        for ell in ("1","3"):
            found=float(t1[state][ell][0])
            original_n=float(original["conditional_Fourier_source_projected_ells_not_xi"]
                             ["512"]["conditional_positive_LOS"][state][ell])
            rel=abs(found-original_n)/max(abs(original_n),1e-14)
            require(rel<levels["engineering_QA_warning_not_physical_threshold"]
                    ["expected_fixed_E10_rank1_tolerance"],
                    "Original E10 rank1 angular result did not reproduce")
            qa[f"rank1_original_E10_{state}_ell{ell}"]=rel
    for state in STATES:
        for ell in ("1","3"):
            g48=results["mu_512_GH_48"][state][ell]["Gaussian_squeezed_response_E_rank_times_T"]
            g96=results["mu_512_GH_96"][state][ell]["Gaussian_squeezed_response_E_rank_times_T"]
            g256=results["mu_256_GH_96"][state][ell]["Gaussian_squeezed_response_E_rank_times_T"]
            qa[f"GH_48_vs_96_{state}_ell{ell}"]=abs(g48-g96)/max(abs(g48),abs(g96),1e-14)
            qa[f"mu_256_vs_512_{state}_ell{ell}"]=abs(g256-g96)/max(abs(g256),abs(g96),1e-14)
    stein={}
    for step in (.002,.001):
        stein[str(step)]={}
        for state in STATES:
            stein[str(step)][state]=stein_gaussian(q,fs,cases,state,96,512,step)
            for ell in ("1","3"):
                measured=results["mu_512_GH_96"][state][ell][
                    "Gaussian_squeezed_response_E_rank_times_T"]
                derived=stein[str(step)][state][ell]
                qa[f"Stein_integrate_{step}_{state}_ell{ell}"]=abs(measured-derived)/max(
                    abs(measured),abs(derived),1e-14)
    for state in STATES:
        for ell in ("1","3"):
            a=stein["0.002"][state][ell];b=stein["0.001"][state][ell]
            qa[f"Stein_steps_{state}_ell{ell}"]=abs(a-b)/max(abs(a),abs(b),1e-14)
    qs=levels["engineering_QA_warning_not_physical_threshold"]
    for key,val in qa.items():
        if key.startswith("GH_48_vs_96_") and val>qs["rank_48_vs_96_relative_warn"]:
            warnings.append(key)
        if key.startswith("mu_256_vs_512_") and val>qs["angular_256_vs_512_relative_warn"]:
            warnings.append(key)
        if key.startswith("Stein_integrate_") and val>qs["Stein_derivative_0p002_vs_0p001_relative_warn"]:
            warnings.append(key)
        if key.startswith("Stein_steps_") and val>qs["Stein_derivative_0p002_vs_0p001_relative_warn"]:
            warnings.append(key)
    require(not warnings,"Preregistered E11 Gaussian/wind source QA warning: "+",".join(warnings))
    gamma=results["mu_512_GH_96"]["delta_plus_minus"]
    # "nonzero" is numerically diagnosed, not a pre-registered minimum physical
    # significance nor a model-independent guarantee of actual eBOSS detectability.
    result={
      "date":"2026-09-27",
      "status":"E11_PINNED_SOURCE_ONLY_GAUSSIAN_WIND_SQUEEZED_RESPONSE_COMPLETE_TWO_POINT_INTRINSIC_MEAN_ZERO_THREE_POINT_REQUIRES_INDEPENDENT_P_DELTA_V",
      "precomputed_prereg_protocol_git_blob_sha1":PROTOCOL_BLOB,
      "frozen_E8_CSV_sha256":e10.ORIG_CSV_SHA,
      "original_E9_three_CLASS_report_sha256":e10.ORIG_E9_SHA,
      "original_E10_conditional_angular_report_sha256":E10_ORIGINAL_SHA,
      "original_CLASS_commit":e10.e8.cro.CLASS_COMMIT,
      "definition":"Same original E9 Gaussian +1sigma LOS wind rank r~N(0,1) correlated between ALL states/z nodes. Exact original Eq20 F_s angular occupancy and original Eq A8 locally uniform wind two-tracer C_LE orientation. T_ell(r) is unit DeltaBias*Pcb conditional Fourier multipole at one short kCOM=.05 h/Mpc z=.95. Gamma_ell=E[r*T_ell(r)]=E[dT_ell/dr].",
      "gaussian_long_density_3point_formula":"For jointly mean-zero Gaussian long scalar D_L(K) and standardized long LOS relative wind r(x), E[D_L(K)*T_ell(r(x))]=exp(-iK dot x)*P_(D_L,r)(K)*Gamma_ell under leading squeezed separate-background short-long independence; requires separately computed/calibrated genuinely LONG K and actual survey tracer mapping. Is NOT an unconditional 24D eBOSS odd two-point prediction.",
      "physical_wind_density_cross_power_role":"P_(D_L,r)(K)=i*(Khat dot LOS)*T_D(K)*V_rel(K)*P_zeta(K)/sigma_LOS up to Fourier and signed-velocity conventions. Not computed here; do not set it to unity or infer its value from single-mode E9 k=.05 phase.",
      "all_angular_and_rank_quadrature_results":results,
      "Stein_identity_independent_centered_derivative_numerical_checks":stein,
      "engineering_QA_relative_gaps":qa,
      "engineering_QA_warnings":warnings,
      "preferred_unit_kernel":"mu_512_GH_96; must still multiply by actual long P_delta_r, DeltaBias, P_cb and physical galaxy response",
      "preferred_plus_minus_Gamma_difference":gamma,
      "zero_intrinsic_unconditional_odd_mean_under_symmetric_model":True,
      "nonzero_Gamma_not_24D_eBOSS_detection_or_exclusion":True,
      "Gaussian_Wick_3pt_cannot_be_predicted_from_9_eBOSS_mocks":True,
      "actual_eBOSS_long_velocity_density_transfer_and_halo_bias_calibrated":False,
      "full_short_k_long_K_redshift_and_survey_3point_window_available":False,
      "A04_valid_independent_bispectrum_covariance_available":False,
      "observed_galaxy_random_or_odd_data_read":False,
      "new_catalogue_or_mock_download":False,
      "main_mutation":False
    }
    return result

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--e8-csv",type=Path,required=True)
    ap.add_argument("--e9-json",type=Path,required=True)
    ap.add_argument("--e10-json",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    result=analyze(args.e8_csv,args.e9_json,args.e10_json)
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    output_new(args.output,raw)
    print("EBOSS_A03_E11_FROZEN_GAUSSIAN_WIND_SQUEEZED_3POINT_RESPONSE_COMPLETE",
          "SHA256",hashlib.sha256(raw).hexdigest(),
          "WARNINGS",len(result["engineering_QA_warnings"]),
          "TWO_POINT_ODD_MEAN_ZERO OBSERVED_ODD_SEALED",flush=True)
    print("E11_UNIT_SQUEEZED_RESPONSE_DELTA_GAMMA",json.dumps(
          result["preferred_plus_minus_Gamma_difference"],sort_keys=True),flush=True)
    print("E11_MAX_NUMERIC_QA_GAP",max(result["engineering_QA_relative_gaps"].values()),
          flush=True)
if __name__=="__main__":
    main()
