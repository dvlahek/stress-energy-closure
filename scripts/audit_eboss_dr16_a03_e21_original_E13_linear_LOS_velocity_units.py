#!/usr/bin/env python3
"""E21 independent exact-original-E13 linear relative LOS velocity units QA.

POSTHOC original E13 discovery, E21 prospective QA. No new CLASS,
FITS, Abacus, WSL, observed galaxy/odd, mock, physical tracer response.
"""
from __future__ import annotations
import copy
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PINS={
"protocol":("source_data/eboss_dr16_a03_e21_linear_LOS_velocity_two_tracer_response_normalization_protocol_2026-09-29.json","7fc84392a915bf3b242bbc12a27e78267dd643f1"),
"result":("source_data/eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json","cc3312c5313cce60fad4a88c3a56ef3d4e675e38"),
"summary":("source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json","d571ef664a4f4c94a0f24e1ba535fbf7c9294f76"),
"FD":("source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_FD.json","d8be18c960bd9a42792e77523a551b70a6aece50"),
"Fplus":("source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_plus.json","60650738ebbd0ca7b4f14780021e501459213f99"),
"Fminus":("source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_minus.json","2e76cf760af7236b1dc4ef5a1260e9b66c468ecf"),
"E12code":("scripts/audit_eboss_dr16_a03_e12_direct_class_vTk_long_cross.py","0f516a1525fc5df44e71c400ea0cfb30664a1ce5"),
"E20":("source_data/eboss_dr16_a03_e20_archived_exact_moment_and_cross_power_structural_QA_2026-09-29.json","c736ea3bca7bd0560ab07225294e1dd0930e0561"),
}
K=("0.001","0.002","0.003","0.005")
STATES=("FD","Fplus","Fminus")

def need(good,msg):
    if not good:raise ValueError("E21_SOURCE_ONLY_STOP: "+msg)
def close(a,b,msg):
    need(isinstance(a,(float,int)) and isinstance(b,(float,int)) and
         math.isfinite(a) and math.isfinite(b) and
         math.isclose(a,b,rel_tol=5e-14,abs_tol=1e-12),msg)
def read_exact():
    out={}
    for tag,(rel,sha) in PINS.items():
        p=ROOT/rel
        need(p.is_file() and not p.is_symlink(),"missing original source "+tag)
        raw=p.read_bytes()
        need(hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()==sha,
             "Git source bytes changed "+tag)
        out[tag]=raw if tag=="E12code" else json.loads(raw)
    p=out["protocol"];res=out["result"];e13=out["summary"]
    need(p["registration_honesty"].startswith("E21 is explicitly POSTHOC") and
         all(p["STOP"].values()),"prospective math-only QA contract violated")
    need(b"P_<delta_cb,v_LOS/sigma> = +i " in out["E12code"] and
         e13["observed_galaxy_random_or_odd_read"] is False and
         e13["new_catalogue_or_mock_download_or_science_cut_seed"] is False and
         e13["status"].startswith("E13_RANK_MATCHED_DIRECT_VTK") and
         out["E20"]["no_observed_odd_access"] is True,
         "original E12 Fourier sign or E13/E20 original data seal changed")
    need(res["no_observed_eBOSS_odd_access"] and res["no_user_WSL"] and
         res["new_physical_halo_galaxy_response"] is False and
         res["no_covariance_snr_or_pvalues"] and
         res["source_git_blobs"]["summary"]==PINS["summary"][1],
         "new E21 physical/result scope changed")
    for s in STATES:
        q=out[s]
        expected="FD" if s=="FD" else s[1:].lower()
        need(q["state"]==expected and
             q["observed_galaxy_random_or_odd_read"] is False and
             q["original_4000q_E8_csv_sha256"]==
             "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0",
             "original E13 state/cosmology/seal drift "+s)
    need(list(res["original_long_modes"])==list(K),"original four K modified")
    return out

def replay(original):
    res=original["result"]
    agg=original["summary"]["conditional_Bsource_unit_DeltaBias_per_E8_state"]
    sigma={}
    for state in STATES:
        q=original[state]
        sigma[state]=q["direct_R16_sigma_LOS_kms_three_fixed_z"]["0.95"]
        need(sigma[state]>0 and math.isfinite(sigma[state]),"sigma invalid "+state)
        need(list(q["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"])==list(K),
             "E13 original K subset/expansion "+state)
    for k in K:
        row=res["original_long_modes"][k]
        need(row["k_h_per_Mpc"]==float(k) and
             set(row["states"])==set(STATES),"E21 K/state pairing failed")
        for state in STATES:
            source=original[state]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"][k]
            reported=row["states"][state]
            C=source["P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3"]
            a="FD" if state=="FD" else state[1:].lower()
            close(C,agg[a]["long_K_fixed"][k]["P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"],
                  "original independent E13 aggregate C disagrees "+state+"/"+k)
            need(source["k_com_h_Mpc"]==float(k) and C>0,
                 "original C or K invalid "+state+"/"+k)
            close(reported["sigma_R16_LOS_km_per_s"],sigma[state],"original state sigma")
            close(reported["C_normalized_P_delta_r_over_i_mu_Mpc3"],C,"original normalized C")
            close(reported["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"],C*sigma[state],
                  "unscaled physical-unit velocity C")
        plus=row["states"]["Fplus"];minus=row["states"]["Fminus"]
        signed=row["signed_Fplus_minus_Fminus"]
        dx=plus["C_normalized_P_delta_r_over_i_mu_Mpc3"]-minus["C_normalized_P_delta_r_over_i_mu_Mpc3"]
        dy=plus["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]-minus["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        close(signed["normalized_Mpc3"],dx,"normalized contrast")
        close(signed["unnormalized_km_s_Mpc3"],dy,"physical velocity contrast")
        need(row["opposite_signed_contrasts"]==(math.copysign(1,dx)!=math.copysign(1,dy)),
             "contrast sign bookkeeping")
    need(res["contrasts_opposite_sign_K"]==["0.003","0.005"],
         "original fixed cohort sign inversion was postselected or altered")

def cross_real_imag(A_L,A_E,c_L,c_E,C,mu):
    # Original E12 E13 P_delta,v = + i*mu*C, P_v,delta=-i*mu*C.
    pdelta_v=complex(0,mu*C)
    pv_delta=pdelta_v.conjugate()
    p=A_L*c_E*pdelta_v+c_L*A_E*pv_delta
    return p.imag

def signs(original):
    for state in STATES:
        C=original["result"]["original_long_modes"]["0.002"]["states"][state]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        sigma=original[state]["direct_R16_sigma_LOS_kms_three_fixed_z"]["0.95"]
        cr=1.5*sigma
        a=cross_real_imag(2,3,1.5,1.5,C,.4)
        b=cross_real_imag(3,2,1.5,1.5,C,.4)
        need(math.isclose(a,-b,abs_tol=1e-8,rel_tol=1e-14),
             "original LRG/ELG conjugation")
        need(math.isclose(a,-cross_real_imag(2,3,1.5,1.5,C,-.4),
                          abs_tol=1e-8,rel_tol=1e-14),"mu orientation parity")
        rn=original["result"]["original_long_modes"]["0.002"]["states"][state]["C_normalized_P_delta_r_over_i_mu_Mpc3"]
        ar=cross_real_imag(2,3,cr,cr,rn,.4)
        need(math.isclose(a,ar,abs_tol=1e-8,rel_tol=1e-14),
             "invalid cr=cv*sigma field reparametrization")
        need(cross_real_imag(2,2,1,1,C,.4)==0 and
             cross_real_imag(2,2,1,2,C,.4)!=0 and
             cross_real_imag(2,3,0,0,C,.4)==0,
             "false universal equal bias or forced nonzero galaxy response")

def reject(name,fn):
    try:fn()
    except (ValueError,KeyError,TypeError):print("E21_NEGATIVE_REJECT",name,flush=True)
    else:raise AssertionError("E21 accepted invalid synthetic perturbation "+name)
def main():
    x=read_exact()
    replay(x);signs(x)
    tamper=copy.deepcopy(x);tamper["Fplus"]["direct_R16_sigma_LOS_kms_three_fixed_z"]["0.95"]*=1.01
    reject("STATE_SIGMA_TAMPER",lambda:replay(tamper))
    tamper=copy.deepcopy(x);tamper["Fminus"]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"]["0.003"]["P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3"]*=1.01
    reject("ORIGINAL_LONG_K_TAMPER",lambda:replay(tamper))
    tamper=copy.deepcopy(x);tamper["summary"]["conditional_Bsource_unit_DeltaBias_per_E8_state"]["plus"]["long_K_fixed"]["0.002"]["P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"]+=2
    reject("INDEPENDENT_E13_SUMMARY_TAMPER",lambda:replay(tamper))
    tamper=copy.deepcopy(x);tamper["result"]["original_long_modes"]["0.005"]["signed_Fplus_minus_Fminus"]["unnormalized_km_s_Mpc3"]*=-1
    reject("CONTRAST_SIGN_TAMPER",lambda:replay(tamper))
    C=x["result"]["original_long_modes"]["0.001"]["states"]["FD"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
    reject("FALSE_EQUAL_BIAS_UNIVERSAL_ZERO",
           lambda:need(cross_real_imag(2,2,1,2,C,.4)==0,
                       "different response can produce linear imaginary cross"))
    print("E21_FIVE_SOURCE_ONLY_SYNTHETIC_NEGATIVES_PASS",flush=True)
    for k in K:
        z=x["result"]["original_long_modes"][k]["signed_Fplus_minus_Fminus"]
        print("E21_ORIGINAL_CLASS_LONG_K",k,
              "NORMALIZED_FPLUS_MINUS_FMINUS_MPC3",format(z["normalized_Mpc3"],".14g"),
              "VELOCITY_FPLUS_MINUS_FMINUS_KMS_MPC3",format(z["unnormalized_km_s_Mpc3"],".14g"),
              flush=True)
    print("E21_ORIGINAL_E13_3_STATES_4_LONG_K_2_UNITS_AND_TRACER_PARITY_PASS",flush=True)
    print("E21_NO_PHYSICAL_TRACER_COUPLING_NO_NEW_CLASS_NO_OBSERVED_ODD_NO_WSL",flush=True)
if __name__=="__main__":main()
