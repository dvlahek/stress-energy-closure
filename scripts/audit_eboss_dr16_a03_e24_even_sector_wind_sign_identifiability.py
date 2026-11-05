#!/usr/bin/env python3
"""E24 restricted linear LOS neutrino response: exact even/odd identifiability.

Prospectively test ONLY new algebraic and synthetic witnesses on original
E21/E22/E23 source bytes. Not galaxy HOD/halo physics or an eBOSS odd read.
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
"protocol":("source_data/eboss_dr16_a03_e24_even_sector_linear_wind_response_sign_identifiability_protocol_2026-09-29.json","3fd79b8cfa116b05c5ab55a248f0f95f2ec7a45b"),
"E21":("source_data/eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json","cc3312c5313cce60fad4a88c3a56ef3d4e675e38"),
"E22":("source_data/eboss_dr16_a03_e22_original_E21_four_K_dual_channel_physical_gate_source_only_summary_2026-09-29.json","a7473fb8b052f37d65e13456d04ad3a885920684"),
"E23":("docs/EBOSS_DR16_A03_E23_HIGHZ_PAIR_REDSHIFT_AND_PUBLISHED_HOD_TRANSPORT_GATE_2026-09-29.md","8fedf17c15d97d49ca45289dc1c0d6cb87b74727")
}
K=("0.001","0.002","0.003","0.005")
STATES=("FD","Fplus","Fminus")

def need(ok,msg):
    if not ok:raise ValueError("E24_RESTRICTED_SOURCE_ONLY_STOP: "+msg)

def load():
    x={}
    for name,(relative,pin) in PINS.items():
        p=ROOT/relative
        need(p.is_file() and not p.is_symlink(),"missing parent "+name)
        raw=p.read_bytes()
        gitsha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
        need(gitsha==pin,"immutable Git SHA mismatch "+name)
        x[name]=raw if name=="E23" else json.loads(raw)
    p=x["protocol"];e=x["E21"];f=x["E22"]
    need(p["research_boundary"].startswith("Original E21 linear LOS response") and
         p["original_audit_head"]=="c57fc17fe1752f5372cb099c680099dd2dc4fb2a" and
         all(p["STOP"].values()) and
         p["frozen_parents"]["E23_high_z_support"]["git_blob"]==PINS["E23"][1],
         "E24 registration/STOP violation")
    need(e["no_observed_eBOSS_odd_access"] and
         e["new_physical_halo_galaxy_response"] is False and
         e["source_git_blobs"]["Fplus"]==
         "60650738ebbd0ca7b4f14780021e501459213f99" and
         f["physics_status"]["observed_odd_SEALED"] and
         f["physics_status"]["real_physical_beta_and_Doppler_tracer_c_missing"] and
         f["no_new_CLASS_Abacus_eBOSS_FITS_WSL"],
         "E21/E22 original physical STOP drift")
    need(b"[0.9,1.0)" in x["E23"] and
         b"M_{200c}" in x["E23"],"E23 high-z/HOD source scope drift")
    need(list(e["original_long_modes"])==list(K) and
         list(f["original_E21_physical_velocity_units_ratios"])==list(K),
         "fixed four E21 original K set drift")
    return x

def replay_original(x):
    e=x["E21"];f=x["E22"]
    for key in K:
        row=e["original_long_modes"][key];r=f["original_E21_physical_velocity_units_ratios"][key]
        need(set(row["states"])==set(STATES) and row["k_h_per_Mpc"]==float(key),
             "original E21 K or state drift "+key)
        z={}
        for s in STATES:
            q=row["states"][s]
            C=q["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
            sigma=q["sigma_R16_LOS_km_per_s"]
            normalized=q["C_normalized_P_delta_r_over_i_mu_Mpc3"]
            need(all(math.isfinite(v) and v>0 for v in (C,sigma,normalized)),
                 "nonpositive original E21 physical cross "+s+"/"+key)
            need(math.isclose(C,sigma*normalized,rel_tol=5e-14,abs_tol=1e-8),
                 "E21 state-wise sigma conversion "+s+"/"+key)
            z[s]=C
        diff=z["Fplus"]-z["Fminus"]
        need(math.isclose(diff,row["signed_Fplus_minus_Fminus"]["unnormalized_km_s_Mpc3"],
                          rel_tol=2e-12,abs_tol=1e-6),
             "original E21 signed Fpm velocity contrast "+key)
        mean=(z["Fplus"]+z["Fminus"])/2
        for k,actual in (("Fplus_over_Fminus",z["Fplus"]/z["Fminus"]),
                         ("Fplus_minus_Fminus_over_two_state_mean",diff/mean),
                         ("Fplus_minus_Fminus_over_FD",diff/z["FD"])):
            need(math.isclose(r[k],actual,rel_tol=1e-12,abs_tol=1e-12),
                 "original E22 ratio changed "+key+"/"+k)
    return True

def spectrum(A_L,A_E,c_L,c_E,mu,Pdd,Pvv,C):
    """Exact real auto/even-cross and imag LE under P_delta,v=+i*mu*C.

    Inputs deliberately rational; Pvv is an even positive LOS-velocity power.
    """
    args=(A_L,A_E,c_L,c_E,mu,Pdd,Pvv,C)
    need(all(isinstance(v,F) for v in args),"exact arithmetic only")
    need(Pdd>=0 and Pvv>=0 and C>=0,"nonphysical toy auto power/cross coefficient")
    return {
      "LL_even":A_L*A_L*Pdd+c_L*c_L*Pvv,
      "EE_even":A_E*A_E*Pdd+c_E*c_E*Pvv,
      "LE_even":A_L*A_E*Pdd+c_L*c_E*Pvv,
      "LE_odd_im":mu*(A_L*c_E-A_E*c_L)*C
    }

def exact_witnesses():
    A_L,A_E,mu,Pdd,Pvv,C=F(2),F(3),F(1,2),F(11),F(5),F(7)
    c_L,c_E=F(1,10),F(2,10)
    base=spectrum(A_L,A_E,c_L,c_E,mu,Pdd,Pvv,C)
    flip=spectrum(A_L,A_E,-c_L,-c_E,mu,Pdd,Pvv,C)
    even=("LL_even","EE_even","LE_even")
    need(all(base[k]==flip[k] for k in even) and
         base["LE_odd_im"]== -flip["LE_odd_im"] and
         base["LE_odd_im"]!=0,"exact global response sign duality")
    swapped=spectrum(A_E,A_L,c_E,c_L,mu,Pdd,Pvv,C)
    reversed_mu=spectrum(A_L,A_E,c_L,c_E,-mu,Pdd,Pvv,C)
    need(swapped["LE_even"]==base["LE_even"] and
         swapped["LE_odd_im"]==-base["LE_odd_im"] and
         all(reversed_mu[k]==base[k] for k in even) and
         reversed_mu["LE_odd_im"]==-base["LE_odd_im"],
         "exact LE/EL and LOS orientation")
    zero=spectrum(A_L,A_E,F(0),F(0),mu,Pdd,Pvv,C)
    eps=F(1,10)
    perturbed=spectrum(A_L,A_E,F(0),eps,mu,Pdd,Pvv,C)
    neg=spectrum(A_L,A_E,F(0),-eps,mu,Pdd,Pvv,C)
    need(all(perturbed[k]-neg[k]==0 for k in even) and
         perturbed["EE_even"]-zero["EE_even"]==eps*eps*Pvv and
         perturbed["LE_odd_im"]-zero["LE_odd_im"]==mu*A_L*eps*C and
         neg["LE_odd_im"]==-perturbed["LE_odd_im"],
         "first-order even null or second-order even correction")
    nullshift=F(3,10)
    sameodd=spectrum(A_L,A_E,c_L+nullshift*A_L,c_E+nullshift*A_E,
                     mu,Pdd,Pvv,C)
    need(sameodd["LE_odd_im"]==base["LE_odd_im"] and
         sameodd["LE_even"]!=base["LE_even"],
         "odd-only tracer rank-null without claiming full even null")
    noresponse=spectrum(A_L,A_E,F(0),F(0),mu,Pdd,Pvv,C)
    need(noresponse["LE_odd_im"]==0,"spurious nonzero without physical c")
    return base,flip,zero,perturbed,swapped,sameodd

def mass_pair_marginals():
    # Same one-tracer 50/50 two-bin occupation weights for BOTH models.
    # x,y label an arbitrary latent host/environment class +-1,
    # not actual LRG/ELG M200c or a physical spatial pairing.
    diag=((F(1,2),F(0)),(F(0),F(1,2)))
    anti=((F(0),F(1,2)),(F(1,2),F(0)))
    def examine(j):
        rows=[sum(j[i]) for i in range(2)]
        cols=[sum(j[i][k] for i in range(2)) for k in range(2)]
        pair=sum(j[i][k]*F((-1 if i==0 else 1)*(-1 if k==0 else 1))
                 for i in range(2) for k in range(2))
        return rows,cols,pair
    m1=examine(diag);m2=examine(anti)
    need(m1[:2]==m2[:2]==([F(1,2),F(1,2)],[F(1,2),F(1,2)]) and
         m1[2]==1 and m2[2]==-1,
         "independent one-tracer HOD marginals do not identify paired environment")
    return m1,m2

def reject(label,call):
    try:call()
    except (ValueError,KeyError,TypeError):
        print("E24_SYNTHETIC_NEGATIVE_REJECT",label,flush=True)
    else:raise AssertionError("E24 accepted deliberately false claim "+label)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    data=load()
    replay_original(data)
    base,flip,zero,eps,swapped,sameodd=exact_witnesses()
    pair1,pair2=mass_pair_marginals()
    # Eight fail-closed negative tests, no science thresholds or physical fit.
    bad=copy.deepcopy(data)
    bad["protocol"]["STOP"]["no_observed_LRG_ELG_galaxies_or_24D_odd"]=False
    reject("OBSERVED_SEAL_TAMPER",
           lambda:need(all(bad["protocol"]["STOP"].values()),"observed seal broken"))
    bad=copy.deepcopy(data)
    bad["E21"]["original_long_modes"]["0.003"]["states"]["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]*=1.01
    reject("ORIGINAL_E21_C_TAMPER",lambda:replay_original(bad))
    bad=copy.deepcopy(data)
    del bad["E21"]["original_long_modes"]["0.005"]
    reject("ORIGINAL_LONG_K_REMOVAL",lambda:replay_original(bad))
    reject("FALSE_EVEN_GLOBAL_SIGN_CHANGE",
           lambda:need(base["LE_even"]!=flip["LE_even"],
                       "even LE cross invariant under both c signs"))
    reject("FALSE_ODD_GLOBAL_SIGN_INVARIANCE",
           lambda:need(base["LE_odd_im"]==flip["LE_odd_im"],
                       "odd LE reverses sign under both c signs"))
    reject("FALSE_FIRST_ORDER_EVEN_DERIVATIVE",
           lambda:need(eps["EE_even"]-zero["EE_even"]!=
                       spectrum(F(2),F(3),F(0),-F(1,10),F(1,2),F(11),F(5),F(7))["EE_even"]-zero["EE_even"],
                       "even is quadratic in c in this exact model"))
    reject("FALSE_TRACER_ODD_SAME_SIGN",
           lambda:need(base["LE_odd_im"]==swapped["LE_odd_im"],
                       "LE/EL must reverse odd orientation"))
    reject("FALSE_HOD_MARGINAL_PAIR_IDENTIFICATION",
           lambda:need(pair1[2]==pair2[2],
                       "same one-tracer marginals do not identify two-host joint"))
    reject("FALSE_ODD_ZERO_RESPONSE_DETECTION",
           lambda:need(zero["LE_odd_im"]!=0,"no neutrino response c, no template odd"))
    print("E24_NINE_SYNTHETIC_NEGATIVE_CONTROLS_PASS",flush=True)
    for key in K:
        z=data["E21"]["original_long_modes"][key]["states"]
        print("E24_PINNED_E21_LINEAR_CV_K",key,
              "FD",format(z["FD"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"],".14g"),
              "Fplus",format(z["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"],".14g"),
              "Fminus",format(z["Fminus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"],".14g"),
              flush=True)
    out={
      "stage":data["protocol"]["stage"],
      "status":"EXACT_RATIONAL_RESTRICTED_LINEAR_TWO_TRACER_EVEN_SIGN_DEGENERACY_PASS",
      "protocol_git_blob":PINS["protocol"][1],
      "original_E21_git_blob":PINS["E21"][1],
      "original_E22_git_blob":PINS["E22"][1],
      "original_E23_git_blob":PINS["E23"][1],
      "rational_witness":{"baseline":{k:str(v) for k,v in base.items()},
                          "simultaneous_c_sign_flip":{k:str(v) for k,v in flip.items()},
                          "zero_c":{k:str(v) for k,v in zero.items()},
                          "epsilon_E_only":{k:str(v) for k,v in eps.items()},
                          "LE_EL_tracer_swap":{k:str(v) for k,v in swapped.items()},
                          "odd_only_null_shift":{k:str(v) for k,v in sameodd.items()}},
      "same_one_tracer_marginals_different_pair_structure":
      {"marginals_both_cases":["1/2","1/2"],"diagonal_pair_moment":str(pair1[2]),
       "antidiagonal_pair_moment":str(pair2[2]),"not_eBOSS_or_physical_halo_pair_data":True},
      "all_3_original_E21_F_and_4_long_K_physically_consistent_source_replayed":True,
      "synthetic_negative_controls_passed":9,
      "no_physical_nu_wind_galaxy_c_calibrated":True,
      "no_actual_highz_HOD_or_pair_response_calibrated":True,
      "no_eBOSS_odd_or_mock_or_new_CLASS_FITS_ASDF_Abacus_WSL":True
    }
    if args.out:
        need(args.out.suffix==".json" and not args.out.exists(),
             "output path must be new JSON")
        with args.out.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
        print("E24_SMALL_SOURCE_ONLY_QA_OUTPUT",args.out,flush=True)
    print("E24_EXACT_ALL_ORDER_RESTRICTED_EVEN_SIGN_DEGENERACY_PASS",flush=True)
    print("E24_FIRST_ORDER_EVEN_ZERO_ODD_LINEAR_RANK_ONE_PASS",flush=True)
    print("E24_NO_PHYSICAL_GALAXY_C_NO_OBSERVED_ODD_NO_WSL",flush=True)

if __name__=="__main__":main()
