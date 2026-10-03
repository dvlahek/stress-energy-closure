#!/usr/bin/env python3
"""E15: exact archived E14 reduced source -> long-LOS parity / toy selection audit.

E14 archived S=B_source/(i mu_long DeltaBias) is NOT a complete triangle
bispectrum or a survey window. This script restores i*mu_long*S only as an
explicit restricted algebraic model, measures necessary angular parity and
demonstrates an entirely synthetic odd angular-selection monopole leakage.
It NEVER calls CLASS, opens FITS, reads observed odd, changes the frozen F states,
or claims a galaxy bispectrum, covariance, S/N or true survey acceptance.
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
from numpy.polynomial.legendre import leggauss, legval

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"source_data/eboss_dr16_a03_e15_long_los_parity_and_selection_leakage_prereg_2026-09-27.json"
E14=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json"
E14MANIFEST=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/archive_manifest.json"
E14PRE=ROOT/"source_data/eboss_dr16_a03_e14_pre_registered_multik_direct_vTk_source_shape_2026-09-27.json"
E14INDEPENDENT=ROOT/"source_data/eboss_dr16_a03_e14_independent_original_source_replay_2026_09_27.json"
BLOBS={
 P:"71a8cb29314036f5899af5fd265d3a8779e46be5",
 E14:"6ef0b2c3268a5bae9564b58abdd914328e861cba",
 E14MANIFEST:"ba747f413988693278a6322bffe2da69d4c52d68",
 E14PRE:"688c56b5486395911ff39ff6f862b405c1f446d7",
 E14INDEPENDENT:"02a8aed04da1cc692a5ac36982b63a10013e3594",
}
E14_SHA="0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8"
STATES=("FD","plus","minus")
KS=("0.05","0.075","0.1")
KLS=("0.001","0.002","0.003","0.005")
ELL=("1","3")
LS=tuple(range(5))
NODES=(32,64)
EPS=.1
OUT=ROOT/"eboss_workspace/a03_physics_source/e15_long_los_parity_source_only.json"

def require(ok,msg):
    if not ok:raise ValueError(msg)

def sha(data):return hashlib.sha256(data).hexdigest()

def git_blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def source_gate():
    for path,blob in BLOBS.items():
        require(git_blob(path.read_bytes())==blob,"Frozen original/preregistered Git blob changed: "+str(path))
    pre=json.loads(P.read_bytes())
    require(pre["source_pins"]["E14_original_full_json_sha256"]==E14_SHA
            and pre["angular_model"]["GL_nodes"]==list(NODES)
            and pre["selection_negative_controls"]["epsilon_synthetic"]==EPS
            and pre["fixed_original_axes"]["states"]==list(STATES)
            and pre["fixed_original_axes"]["n_original_reduced_state_source"]==72
            and pre["fixed_original_axes"]["n_original_Fplus_minus_Fminus_contrasts"]==24
            and pre["prohibitions"]["no_observed_odd_vector_unsealing"] is True
            and pre["prohibitions"]["no_main_mutation"] is True,
            "Prospective E15 source, scope, angular bins, toy epsilon or observed-odd seal changed")
    raw=E14.read_bytes()
    require(sha(raw)==E14_SHA,"Archived ORIGINAL E14 full source SHA changed")
    manifest=json.loads(E14MANIFEST.read_bytes())
    expected=next((a for a in manifest["files"] if
              a["file"]=="e14_frozen_multik_direct_rank_unitbias_source_only.json"),None)
    require(expected is not None and expected["sha256"]==E14_SHA and
            expected["size_bytes"]==len(raw),
            "Original E14 archive manifest is missing correct exact bytes")
    e14=json.loads(raw)
    independent=json.loads(E14INDEPENDENT.read_bytes())
    require(e14["status"].startswith("E14_THREE_FROZEN_CLASS_SHORT_K_EXACT_STATIC_E13_SOURCE_SHAPE_PASS")
            and e14["original_E13_full_JSON_sha256"]==
             "537365b2e758ec215810bb687557771440851cf2869d02578569b1d47d4ec98d"
            and independent["original_E14_sha256"]==E14_SHA
            and independent["numerical_case_count"]==72
            and independent["Fplus_minus_Fminus_contrast_case_count"]==24
            and independent["max_relative_independent_source_replay_gap"]==0
            and e14["scientific_scope"]["observed_galaxies_randoms_odd_sealed"] is True
            and e14["scientific_scope"]["main_untouched"] is True,
            "Original E14 case count/independent numerical closure/observed seal changed")
    require([str(k) for k in e14["short_k_comoving_h_Mpc"]]==list(KS)
            and [str(k) for k in e14["four_original_long_K_comoving_h_Mpc"]]==list(KLS),
            "Original E14 preregistered k axes changed")
    return pre,e14,independent

def Pn(mu,l):
    c=[0.]*l+[1.]
    return legval(mu,c)

def angular_moments(S,n):
    nodes,weights=leggauss(n)
    W_even=1.+EPS*Pn(nodes,2)
    W_odd=1.+EPS*Pn(nodes,1)
    require(np.min(W_even)>0 and np.min(W_odd)>0,"Negative synthetic window")
    mu=nodes
    m={str(L):float((2*L+1)/2 * S * np.sum(weights*Pn(mu,L)*mu)) for L in LS}
    # All returned coefficients are IMAGINARY coefficients of i*mu*S.
    mon_even=float(S*np.sum(weights*W_even*mu)/np.sum(weights*W_even))
    mon_odd=float(S*np.sum(weights*W_odd*mu)/np.sum(weights*W_odd))
    dip_even=float(1.5*S*np.sum(weights*Pn(mu,1)*W_even*mu))
    dip_odd=float(1.5*S*np.sum(weights*Pn(mu,1)*W_odd*mu))
    return {"imag_long_multipoles_over_unit_bias_Mpc6":m,
            "synthetic_weighted_imag_monopole_over_unit_bias_Mpc6":{
              "even_W_1_plus_0p1_P2":mon_even,
              "odd_W_1_plus_0p1_P1":mon_odd},
            "synthetic_weighted_imag_dipole_over_unit_bias_Mpc6":{
              "even_W_1_plus_0p1_P2":dip_even,
              "odd_W_1_plus_0p1_P1":dip_odd}}

def compare(S,calc):
    scale=max(1.,abs(S))
    for key,val in calc["imag_long_multipoles_over_unit_bias_Mpc6"].items():
        target=S if key=="1" else 0.
        require(abs(val-target)/scale < 1e-12,
                "E15 original source long angular Legendre parity/dipole closure failed: "+key)
    m=calc["synthetic_weighted_imag_monopole_over_unit_bias_Mpc6"]
    d=calc["synthetic_weighted_imag_dipole_over_unit_bias_Mpc6"]
    require(abs(m["even_W_1_plus_0p1_P2"])/scale < 1e-12 and
            abs(m["odd_W_1_plus_0p1_P1"]-EPS*S/3.)/scale<1e-12 and
            abs(d["even_W_1_plus_0p1_P2"]-S*(1.+2.*EPS/5.))/scale<1e-12 and
            abs(d["odd_W_1_plus_0p1_P1"]-S)/scale<1e-12,
            "E15 purely synthetic even/odd selection geometry failed")
    # Hermitian reality of Fourier-space real fields for the restricted i*mu*S.
    for mu in (-1.,-.8,-.3,0.,.3,.8,1.):
        B=complex(0.,S*mu)
        Breverse=complex(0.,-S*mu)
        require(Breverse==B.conjugate() and B+Breverse==0j and
                -B==complex(0.,(-S)*mu),
                "E15 long-K vector reversal, conjugate reality or tracer exchange failed")

def compute():
    pre,e14,independent=source_gate()
    all_cases={}
    contrast_cases={}
    maxgap=0.
    seen=0
    for state in STATES:
        all_cases[state]={}
        for k in KS:
            all_cases[state][k]={}
            row=e14["all_original_state_k_short_and_long_results"][state][k]
            for K in KLS:
                all_cases[state][k][K]={}
                krow=row["four_original_filtered_long_modes"][K]
                require(float(krow["K_over_k_short"])<=.1+1e-14,
                        "E14 triangle squeeze ratio outside prospectively locked bound")
                for ell in ELL:
                    S=float(krow["ells"][ell]["Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"])
                    require(math.isfinite(S),"Nonfinite original reduced E14 source")
                    proj={}
                    for n in NODES:
                        a=angular_moments(S,n)
                        compare(S,a)
                        proj[str(n)]=a
                    m32=proj["32"]["imag_long_multipoles_over_unit_bias_Mpc6"]
                    m64=proj["64"]["imag_long_multipoles_over_unit_bias_Mpc6"]
                    gap=max(abs(m32[str(L)]-m64[str(L)])/max(1.,abs(S)) for L in LS)
                    maxgap=max(maxgap,gap)
                    require(gap<1e-12,"Prospectively fixed 32/64 GL source angular QA failed")
                    all_cases[state][k][K][ell]={"original_reduced_E14_S_Mpc6":S,
                        "quadrature_imaginary_coefficient":proj,
                        "GL32_64_max_scaled_gap":gap,
                        "actual_complex_B_at_muLong_plus1_unitBias":{
                           "real":0.,"imag":S},
                        "actual_complex_B_at_muLong_minus1_unitBias":{
                           "real":0.,"imag":-S},
                        "full_physical_eBOSS_triangle_and_window":False}
                    seen+=1
    require(seen==72,"Original E14 72 state/short/long/ell case count lost")
    contrasts=0
    for k in KS:
        contrast_cases[k]={}
        for K in KLS:
            contrast_cases[k][K]={}
            for ell in ELL:
                plus=all_cases["plus"][k][K][ell]["original_reduced_E14_S_Mpc6"]
                minus=all_cases["minus"][k][K][ell]["original_reduced_E14_S_Mpc6"]
                S=plus-minus
                original=float(e14["Fplus_minus_Fminus_reduced_source_difference"][k][K][ell][
                    "Fplus_minus_Fminus_reduced_Bsource_Mpc6"])
                scale=max(1.,abs(plus),abs(minus),abs(original))
                require(abs(S-original)/scale<1e-14,
                        "Original E14 Fplus-minus archived state-wise contrast changed")
                a=angular_moments(S,64);compare(S,a)
                contrast_cases[k][K][ell]={"E14_Fplus_minus_Fminus_reduced_S_Mpc6":original,
                    "source_reconstructed_from_state_values_Mpc6":S,
                    "imag_long_monopole_unweighted_Mpc6":a[
                      "imag_long_multipoles_over_unit_bias_Mpc6"]["0"],
                    "imag_long_dipole_Mpc6":a[
                      "imag_long_multipoles_over_unit_bias_Mpc6"]["1"],
                    "synthetic_odd_selection_imag_monopole_Mpc6":a[
                      "synthetic_weighted_imag_monopole_over_unit_bias_Mpc6"][
                      "odd_W_1_plus_0p1_P1"]}
                contrasts+=1
    require(contrasts==24,"Original E14 24 state-wise Fplus-minus contrast count lost")
    sample=contrast_cases["0.05"]["0.005"]["1"]
    return {"date":"2026-09-27",
       "status":"E15_ORIGINAL_E14_72_LONG_LOS_PARITY_NULL_AND_SYNTHETIC_ODD_SELECTION_SOURCE_ONLY_PASS_NOT_FULL_BISPECTRUM",
       "original_E14_full_SHA256":E14_SHA,
       "E15_prospective_protocol_git_blob":"71a8cb29314036f5899af5fd265d3a8779e46be5",
       "original_E14_independent_replay_git_blob":"02a8aed04da1cc692a5ac36982b63a10013e3594",
       "physics":"Restore only i*mu_long*S(kshort,Klong,short_ell) to the previously reduced E14 unit-DeltaBias source, then perform long-LOS Legendre projections. The short_ell labels are existing original E13 angular coefficients, NOT long-LOS multipoles.",
       "all_72_original_state_short_k_long_K_short_ell_cases":all_cases,
       "all_24_original_Fplus_minus_Fminus_parity_contrasts":contrast_cases,
       "max_32_64_GaussLegendre_scaled_gap":maxgap,
       "analytic_parity_controls":{
          "unweighted_long_LOS_L0_2_3_4_all_zero":True,
          "long_LOS_L1_equals_reduced_E14_source":True,
          "long_K_sign_reversal_imaginary_odd":True,
          "fourier_conjugate_reality_k_to_minus_k":True,
          "tracer_order_reversal_changes_imaginary_sign":True,
          "synthetic_even_angular_window_0p1P2_preserves_monopole_zero":True,
          "synthetic_odd_angular_window_0p1P1_creates_monopole_epsilon_S_over_3":True,
          "synthetic_even_window_long_dipole_multiplies_1_plus_2epsilon_over_5":True,
          "synthetic_odd_window_long_dipole_unchanged":True},
       "example_original_frozen_Fplus_minus_Fminus_kshort_0p05_Klong_0p005_short_ell1":{
          "reduced_E14_source_Mpc6":sample["E14_Fplus_minus_Fminus_reduced_S_Mpc6"],
          "unweighted_long_LOS_monopole_zero":True,
          "unit_bias_imag_long_LOS_dipole_Mpc6":sample["imag_long_dipole_Mpc6"],
          "synthetic_odd_window_leakage_imag_long_LOS_monopole_Mpc6":
            sample["synthetic_odd_selection_imag_monopole_Mpc6"],
          "NOT_eBOSS_selection_or_observation":True},
       "technical_QA_warnings":[],
       "physical_limits":{
           "conditional_source_only_not_true_full_triangle_bispectrum":True,
           "finite_K_over_k_angle_triangle_closure_not_computed":True,
           "true_halo_retarded_Einstein_Vlasov_and_LRG_ELG_bias_HOD_uncalibrated":True,
           "real_eBOSS_cap_chunk_depth_triple_window_not_calculated":True,
           "independent_eBOSS_bispectrum_covariance_absent":True,
           "toy_angular_selection_weight_not_a_measured_eBOSS_window":True,
           "original_unconditional_24D_odd_not_reinterpreted":True,
           "original_observed_odd_sealed":True,
           "new_FITS_mock_download_CLASS_run_or_science_seed_cut":False,
           "main_unchanged":True,
           "PR_remains_draft":True}}

def atomic(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,"Original E15 output exists with different bytes")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e15_",delete=False) as f:
        tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUT)
    args=parser.parse_args()
    result=compute()
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    atomic(args.output,raw)
    x=result["example_original_frozen_Fplus_minus_Fminus_kshort_0p05_Klong_0p005_short_ell1"]
    print("EBOSS_A03_E15_ORIGINAL_E14_72_SOURCE_LONG_LOS_PARITY_ZERO_MONO_AND_SYNTHETIC_SELECTION_QA_PASS",
          "JSON_SHA256",sha(raw),"MAX_GL32_64_SCALED_GAP",
          result["max_32_64_GaussLegendre_scaled_gap"],
          "FPLUS_MINUS_FMINUS_LONG_DIPOLE_MPC6",x["unit_bias_imag_long_LOS_dipole_Mpc6"],
          "TOY_ODD_WINDOW_MONO_MPC6",
          x["synthetic_odd_window_leakage_imag_long_LOS_monopole_Mpc6"],
          "NO_REAL_EBOSS_BISPECTRUM NO_OBSERVED_ODD",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
