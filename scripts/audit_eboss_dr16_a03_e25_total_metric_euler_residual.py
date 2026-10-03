#!/usr/bin/env python3
"""E25: exact finite toy algebra for total-metric Euler residual and LOS odd.

NO physical neutrino-halo response, observed data, new CLASS, WSL or HOD fit.
All model-dependent source papers and E21/E24 outcomes precede registration.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PINS={
 "protocol":("source_data/eboss_dr16_a03_e25_total_metric_euler_residual_and_linear_los_response_protocol_2026-09-29.json","01c8449bf7388b2391b53dd5c7ba069a23c2186c"),
 "E24_note":("docs/EBOSS_DR16_A03_E24_EXACT_EVEN_SECTOR_SIGN_DUALITY_AND_MINIMAL_NEUTRINO_TRACER_RESPONSE_2026-09-29.md","b3a87f1ed42133e164b440e89f45dcbf6f79dc8a"),
 "E24_report":("source_data/eboss_dr16_a03_e24_source_only_even_odd_linear_response_identifiability_summary_2026-09-29.json","ccb64f189941ed0b89efe1462cf8ca86590b17fc"),
 "E22_note":("docs/EBOSS_DR16_A03_E22_TWO_PHYSICAL_CHANNELS_OKOLI_PHASE_SHARED_WAKE_AND_SIGNED_ESTIMAND_2026-09-29.md","d2d08ec43cd2604e9a6207f40a5ac1732699ad8e"),
 "E23_note":("docs/EBOSS_DR16_A03_E23_HIGHZ_PAIR_REDSHIFT_AND_PUBLISHED_HOD_TRANSPORT_GATE_2026-09-29.md","8fedf17c15d97d49ca45289dc1c0d6cb87b74727"),
 "E17D1_note":("docs/EBOSS_DR16_A03_E17D1_RESTRICTED_CAUSAL_SOURCE_HISTORY_NONIDENTIFIABILITY_2026-09-27.md","3a536aa384eec1a202e77c532beb3a55c7333cb2"),
 "E21_Cv":("source_data/eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json","cc3312c5313cce60fad4a88c3a56ef3d4e675e38"),
}
C_LIGHT_KM_S=F(299792458,1000)  # SI exact c m/s -> km/s, NOT a fitted speed.
K=("0.001","0.002","0.003","0.005")
STATES=("FD","Fplus","Fminus")

def need(ok,msg):
    if not ok:
        raise ValueError("E25_ORIGINAL_SOURCE_AND_PHYSICS_SCOPE_STOP: "+msg)

def read_sources():
    x={}
    for key,(rel,pin) in PINS.items():
        path=ROOT/rel
        need(path.is_file() and not path.is_symlink(),"missing pinned source "+key)
        raw=path.read_bytes()
        actual=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
        need(actual==pin,"pinned source Git blob changed "+key)
        x[key]=json.loads(raw) if key in ("protocol","E24_report","E21_Cv") else raw
    p=x["protocol"]
    need(p["registration_status"].startswith("POSTHOC_THEORY_AFTER_E24") and
         p["original_branch_head_at_decision"]=="fb02886a71af2313c9359c317aaf64e2b0c09d95" and
         all(p["STOP"].values()) and
         p["original_sources"]["E17D1_note"]["git_blob"]==PINS["E17D1_note"][1],
         "new E25 preregistration/safety scope drift")
    need(x["E21_Cv"]["no_observed_eBOSS_odd_access"] and
         x["E21_Cv"]["new_physical_halo_galaxy_response"] is False and
         x["E24_report"]["physical_neutrino_to_galaxy_c_calibrated"] is False and
         b"[0.9,1.0)" in x["E23_note"] and
         b"BLOCKED" in x["E17D1_note"],
         "frozen E21/E24/E23/E17D1 source physics scope drift")
    return x

def replay_E21(x):
    e=x["E21_Cv"]
    need(tuple(e["original_long_modes"])==K,"original four long K changed")
    for k in K:
        row=e["original_long_modes"][k]
        need(tuple(row["states"])==STATES and row["k_h_per_Mpc"]==float(k),
             "E21 state or K changed "+k)
        for state in STATES:
            q=row["states"][state]
            cv=q["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
            sigma=q["sigma_R16_LOS_km_per_s"]
            cr=q["C_normalized_P_delta_r_over_i_mu_Mpc3"]
            need(all(math.isfinite(z) and z>0 for z in (cv,sigma,cr)) and
                 math.isclose(cv,sigma*cr,rel_tol=5e-14,abs_tol=1e-8),
                 "E21 source velocity units or state sigma changed "+state+"/"+k)
        plus=row["states"]["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        minus=row["states"]["Fminus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        need(math.isclose(plus-minus,row["signed_Fplus_minus_Fminus"]["unnormalized_km_s_Mpc3"],
                          rel_tol=2e-12,abs_tol=1e-6),
             "E21 signed source changed "+k)
    return True

def B_coefficient(hprime_over_h2,s_magnification,rH,f_evo):
    """ONLY Bonvin/Lepori 2023 number-count Eq (2) local V convention."""
    need(rH>0,"number-count geometric rH must be positive")
    return (1-hprime_over_h2+(5*s_magnification-2)/rH
            -5*s_magnification+f_evo)

def euler_residual(V,Vprime,H,grad_total):
    return Vprime+H*V+grad_total

def local_count(B,V,Vprime,H,grad_total):
    need(H>0,"conformal H must be positive")
    return B*V+(Vprime+grad_total)/H

def local_reparam(B,V,Vprime,H,grad_total):
    return (B-1)*V+euler_residual(V,Vprime,H,grad_total)/H

def geodesic_and_residual_tests():
    H,V,grad_bg,grad_w=F(2),F(1,100),F(3,100),F(1,50)
    B=B_coefficient(F(-1,2),F(2,5),F(3),F(1))
    need(B==F(1,2),"published-number-count coefficient toy sanity")
    grad_total=grad_bg+grad_w
    Vprime_geo=-H*V-grad_total
    total=euler_residual(V,Vprime_geo,H,grad_total)
    background=euler_residual(V,Vprime_geo,H,grad_bg)
    actual=local_count(B,V,Vprime_geo,H,grad_total)
    need(total==0 and background==-grad_w and grad_w!=0 and
         actual==local_reparam(B,V,Vprime_geo,H,grad_total)==(B-1)*V,
         "total-metric geodesic or corrected local observed number counts")
    # A false, background-only residual is exactly canceled by the
    # missing wake potential in the ACTUAL total-metric number counts.
    wrong=(B-1)*V+background/H
    need(wrong!=actual and
         wrong+grad_w/H==actual,
         "background-potential fake residual was counted as a new force")
    R=F(3,1000)  # SYNTHETIC effective non-geodesic residual, not halo data.
    Vprime_non_geo=Vprime_geo+R
    observed=local_count(B,V,Vprime_non_geo,H,grad_total)
    need(euler_residual(V,Vprime_non_geo,H,grad_total)==R and
         observed==actual+R/H and
         observed==local_reparam(B,V,Vprime_non_geo,H,grad_total),
         "genuine additional toy Euler residual missing")
    return {"B":B,"H":H,"V":V,"gradient_background":grad_bg,
            "gradient_wake":grad_w,"Vprime_geodesic":Vprime_geo,
            "geodesic_E_total":total,"geodesic_E_background":background,
            "geodesic_local_number_count":actual,"wrong_background_only":wrong,
            "synthetic_additional_Euler_residual":R,
            "nongeodesic_local_number_count":observed}

def d_response(B,beta,epsilon,selection_wind):
    """delta Delta = d*(v_rel/c_light) under an ACTUAL physical closure."""
    return (B-1)*beta+epsilon+selection_wind

def odd_LE(A_L,A_E,d_L,d_E,mu,C_v_km_s):
    return mu*(A_L*d_E-A_E*d_L)*C_v_km_s/C_LIGHT_KM_S

def physical_units_and_parity(x):
    B_L,B_E=F(1,2),F(4,3)
    beta_L,beta_E=F(2),F(3)
    epsilon_L,epsilon_E=F(0),F(1,4)
    selection_L,selection_E=F(1,5),F(0)
    dL=d_response(B_L,beta_L,epsilon_L,selection_L)
    dE=d_response(B_E,beta_E,epsilon_E,selection_E)
    need(dL==F(-4,5) and dE==F(5,4),"synthetic response decomposition")
    A_L,A_E,mu=F(2),F(3),F(1,2)
    chi=A_L*dE-A_E*dL
    need(chi==F(49,10),"exact unknown-responses-only contrast")
    originalCv=x["E21_Cv"]["original_long_modes"]["0.005"]["states"]["FD"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
    C=F(str(originalCv))  # original SOURCE physical-unit linear coefficient, not galaxy prediction
    im=odd_LE(A_L,A_E,dL,dE,mu,C)
    e21_cL,e21_cE=dL/C_LIGHT_KM_S,dE/C_LIGHT_KM_S
    need(im==mu*(A_L*e21_cE-A_E*e21_cL)*C and
         im==-odd_LE(A_E,A_L,dE,dL,mu,C) and
         im==-odd_LE(A_L,A_E,dL,dE,-mu,C),
         "physical km/s to dimensionless galaxy number count conversion or parity")
    # Same beta for both tracers does not ensure nonzero chi:
    # the aligned response d_a=t*A_a cancels the cross, even with beta!=0.
    same_A_L,same_A_E=F(2),F(3)
    alignedL,alignedE=F(2),F(3)
    need(alignedL!=0 and alignedE!=0 and
         odd_LE(same_A_L,same_A_E,alignedL,alignedE,mu,C)==0,
         "nonzero halo velocity response does not force nonzero observed chi")
    # Kaiser acts via -i*k*mu*v/H, P_delta,v=+i*mu*C.
    # (i*k*mu/H)*(i*mu*C)=-k*mu^2*C/H is REAL and EVEN.
    k,H=F(1,20),F(2)
    kaiser_cross_real=-k*mu*mu*C/H
    need(kaiser_cross_real==(-k*(-mu)*(-mu)*C/H) and
         kaiser_cross_real!=0 and im!=0,
         "standard Kaiser velocity derivative parity")
    return {"dL":dL,"dE":dE,"chi":chi,"actual_E21_Cv_FD_K_0p005":C,
            "synthetic_Im_LE_source_weighted":im,
            "cL_per_km_s":e21_cL,"cE_per_km_s":e21_cE,
            "Kaiser_density_cross_REAL_even":kaiser_cross_real}

def reject(name,callback):
    try:callback()
    except (ValueError,KeyError,TypeError):
        print("E25_NEGATIVE_REJECT",name,flush=True)
    else:
        raise AssertionError("E25 accepted intentionally false claim "+name)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    x=read_sources()
    replay_E21(x)
    g=geodesic_and_residual_tests()
    u=physical_units_and_parity(x)
    # New *exact* negative tests only. None touches real eBOSS data.
    bad=copy.deepcopy(x)
    bad["protocol"]["STOP"]["observed_odd_24D_SEALED"]=False
    reject("OBSERVED_ODD_SEAL_TAMPER",
           lambda:need(all(bad["protocol"]["STOP"].values()),"observed seal"))
    bad=copy.deepcopy(x)
    bad["E21_Cv"]["original_long_modes"]["0.003"]["states"]["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]*=1.01
    reject("ORIGINAL_E21_CV_TAMPER",lambda:replay_E21(bad))
    reject("FALSE_BACKGROUND_EULER_RESIDUAL_AS_TOTAL",
           lambda:need(g["geodesic_E_background"]==g["geodesic_E_total"],
                       "wake metric gradient omitted"))
    reject("FALSE_GEODESIC_WAKE_EXTRA_FORCE",
           lambda:need(g["wrong_background_only"]==g["geodesic_local_number_count"],
                       "double-counted wake metric"))
    reject("TRUE_NON_GEODESIC_RESIDUAL_DROPPED",
           lambda:need(g["nongeodesic_local_number_count"]==
                       g["geodesic_local_number_count"],
                       "nonzero toy residual deleted"))
    reject("FALSE_KAISER_ODD",
           lambda:need(-F(1,20)*F(1,2)**2*u["actual_E21_Cv_FD_K_0p005"]/F(2)!=
                        -F(1,20)*(-F(1,2))**2*u["actual_E21_Cv_FD_K_0p005"]/F(2),
                       "Kaiser two LOS signs must have same real coefficient"))
    reject("MISSING_SPEED_OF_LIGHT_UNIT_CONVERSION",
           lambda:need(u["synthetic_Im_LE_source_weighted"]==
                       F(1,2)*u["chi"]*u["actual_E21_Cv_FD_K_0p005"],
                       "E21 v physically km/s, observed V dimensionless"))
    reject("FALSE_LE_EL_SAME_SIGN",
           lambda:need(u["synthetic_Im_LE_source_weighted"]==
                       odd_LE(F(3),F(2),u["dE"],u["dL"],F(1,2),
                              u["actual_E21_Cv_FD_K_0p005"]),
                       "LRG/ELG orientation must reverse odd sign"))
    reject("NONZERO_VELOCITY_RESPONSE_IMPLIES_NONZERO_CHI",
           lambda:need(odd_LE(F(2),F(3),F(2),F(3),F(1,2),
                              u["actual_E21_Cv_FD_K_0p005"])!=0,
                       "two aligned nonzero responses can cancel in odd contrast"))
    reject("WRONG_FULL_LOCAL_DOPPLER_B_INSTEAD_OF_B_MINUS_ONE",
           lambda:need(d_response(F(1,2),F(2),F(0),F(0))==F(1),
                       "Euler-reduced D=B-1"))
    print("E25_TEN_SYNTHETIC_NEGATIVE_CONTROLS_PASS",flush=True)
    print("E25_FULL_TOTAL_METRIC_GEODESIC_WAKE_GRADIENT_CANCELS_AS_SEPARATE_FORCE",flush=True)
    print("E25_NONZERO_COARSEGRAINED_TOY_EULER_RESIDUAL_RETAINED",flush=True)
    print("E25_UNIT_CONVERSION_299792p458_KM_S_EXACT_PASS",flush=True)
    print("E25_PINNED_ALL_ORIGINAL_E21_F_STATES_FOUR_LONG_K_PASS",flush=True)
    print("E25_KAISER_EVEN_UNDIFFERENTIATED_DOPPLER_POTENTIAL_ODD_PASS",flush=True)
    report={
      "stage":x["protocol"]["stage"],
      "status":"SOURCE_ONLY_TOTAL_METRIC_EULER_RESIDUAL_AND_LINEAR_LOS_UNIT_IDENTITIES_PASS",
      "original_E25_protocol_git_blob":PINS["protocol"][1],
      "all_parent_git_blobs":{tag:sha for tag,(_path,sha) in PINS.items() if tag!="protocol"},
      "geodesic_and_non_geodesic_exact_toy":{k:str(v) for k,v in g.items()},
      "dimensionless_number_count_and_original_E21_Cv_units_exact_toy":
           {k:str(v) for k,v in u.items()},
      "negative_controls_passed":10,
      "all_3_original_F_four_long_K_replayed":True,
      "E25_not_fully_physical_Einstein_Vlasov_halo_closure":True,
      "physical_beta_epsilon_s_evo_D_halo_selection_all_unmeasured":True,
      "no_observed_galaxy_odd_real_RR_mock_or_new_CLASS_WSL":True,
      "no_eBOSS_significance_or_24D_physical_prediction":True
    }
    if args.out:
        need(args.out.suffix==".json" and not args.out.exists(),
             "output must be a new JSON")
        with args.out.open("x",encoding="utf-8") as f:
            json.dump(report,f,indent=2,allow_nan=False);f.write("\n")
        print("E25_SMALL_SOURCE_ONLY_OUTPUT",args.out,flush=True)
    print("E25_SOURCE_ONLY_PHYSICS_LIMITED_CI_SUCCESS",flush=True)

if __name__=="__main__":main()
