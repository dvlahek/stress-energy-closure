#!/usr/bin/env python3
"""E26: original E12/E13/E21 single-adiabatic-mode source covariance rank.

No CLASS execution, original FITS, galaxy, random, odd, mock, or halo response.
"""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASE="source_data/"
PINS={
"protocol":("source_data/eboss_dr16_a03_e26_original_linear_source_rank_protocol_2026-09-29.json","4b09681d188bd324b9d9f3ba459bd5608ba41195"),
"e12code":("scripts/audit_eboss_dr16_a03_e12_direct_class_vTk_long_cross.py","0f516a1525fc5df44e71c400ea0cfb30664a1ce5"),
"e12_FD":(BASE+"eboss_dr16_a03_e12_archived_CI_2026_09_27/e12_direct_vTk_long_cross_FD.json","5977fe98df367d1e3f595b06b7c0ded753a2ad48"),
"e12_plus":(BASE+"eboss_dr16_a03_e12_archived_CI_2026_09_27/e12_direct_vTk_long_cross_plus.json","9c909d7927cada790bcba4cce573a5b874921f91"),
"e12_minus":(BASE+"eboss_dr16_a03_e12_archived_CI_2026_09_27/e12_direct_vTk_long_cross_minus.json","96cecd44c12ddb02174ad72e1c484d9d701fcee2"),
"e13_FD":(BASE+"eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_FD.json","d8be18c960bd9a42792e77523a551b70a6aece50"),
"e13_plus":(BASE+"eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_plus.json","60650738ebbd0ca7b4f14780021e501459213f99"),
"e13_minus":(BASE+"eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_direct_wind_minus.json","2e76cf760af7236b1dc4ef5a1260e9b66c468ecf"),
"e21":(BASE+"eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json","cc3312c5313cce60fad4a88c3a56ef3d4e675e38")
}
STATES=("FD","plus","minus")
KEYS=("FD","Fplus","Fminus")
K=("0.001","0.002","0.003","0.005")
C_LIGHT=299792.458
TOL=5e-11

def require(ok,msg):
    if not ok:
        raise ValueError("E26_SOURCE_ONLY_STOP "+msg)

def load():
    out={}
    for key,(path,pin) in PINS.items():
        p=ROOT/path
        require(p.is_file() and not p.is_symlink(),"missing original "+key)
        raw=p.read_bytes()
        sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
        require(sha==pin,"original Git byte fingerprint changed "+key)
        out[key]=raw if key=="e12code" else json.loads(raw)
    require(out["protocol"]["negative_controls"]==8 and
            out["protocol"]["no_observed_odd"] and
            out["protocol"]["no_new_CLASS"] and
            out["protocol"]["no_new_WSL"],
            "protocol safety scope")
    require(b'"t_cdm"' in out["e12code"] and b'"t_ncdm[0]"' in out["e12code"],
            "E12 original transfer source provenance")
    require(out["e21"]["no_observed_eBOSS_odd_access"] and
            out["e21"]["new_physical_halo_galaxy_response"] is False,
            "E21 observation seal / physical nonclaim")
    for state in STATES:
        require(out["e12_"+state]["state"]==state and
                out["e13_"+state]["state"]==state and
                out["e13_"+state]["observed_galaxy_random_or_odd_read"] is False,
                "original E12/E13 state or observed seal "+state)
    return out

def replay(x,c_light=C_LIGHT,conjugation=1):
    require(tuple(x["e21"]["original_long_modes"])==K,"E21 four K changed")
    report={}
    for state,target in zip(STATES,KEYS):
        e12=x["e12_"+state];e13=x["e13_"+state]
        s=e12["density_proxy_two_fixed_steps_and_four_long_K"]["0.002"]["long_K"]
        require(len(s)==4 and
                tuple(format(q["K_comoving_h_per_Mpc"],".3f") for q in s)==K,
                "original E12 step 0.002 K cohort "+state)
        z=e13["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"]
        require(tuple(z)==K,"original E13 K cohort "+state)
        out={}
        for key,old in zip(K,s):
            fresh=z[key]
            require(abs(fresh["k_com_h_Mpc"]-float(key))<1e-14 and
                    old["K_Mpc_inv"]>0 and
                    math.isclose(old["K_Mpc_inv"],float(key)*.6736,rel_tol=2e-12),
                    "K unit conversion "+state+"/"+key)
            D=fresh["E13_delta_cb_Newtonian_transfer_per_R"]
            theta=fresh["E13_theta_rel_direct_per_R_Mpc_inv"]
            W=fresh["W_R16_original_locked_top_hat"]
            require(0<W<=1 and
                    math.isclose(D,old["delta_cb_Newtonian_per_primordial_curvature"],rel_tol=TOL) and
                    math.isclose(theta,old["theta_relative_direct_vTk_per_primordial_curvature_Mpc_inv"],rel_tol=TOL),
                    "original E12 E13 direct transfer consistency "+state+"/"+key)
            d2=old["primordial_curvature_Delta2"]
            PR=2*math.pi**2*d2/old["K_Mpc_inv"]**3
            V=-c_light*theta*W/old["K_Mpc_inv"]
            Pdd=D*D*PR
            Pvv=V*V*PR
            Cv_theory=conjugation*(-D*V*PR)
            Cv_original=x["e21"]["original_long_modes"][key]["states"][target]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
            require(Pdd>0 and Pvv>0 and Cv_original>0 and
                    math.isclose(Cv_theory,Cv_original,rel_tol=TOL),
                    "independent original E12 E13 transfer vs E21 cross "+state+"/"+key)
            require("theta_cdm" not in old and
                    "t_cdm" not in old,
                    "original E12 archive is NOT a saved raw t_cdm transfer "+state+"/"+key)
            eps=0.0
            for mu in (1.,.5,0.):
                power_v=mu*mu*Pvv
                cross=mu*Cv_original
                residual=abs(Pdd*power_v-cross*cross)
                scale=Pdd*power_v if mu else 1.0
                require((residual/scale if mu else residual)<=TOL,
                        "non-rank-one original source covariance "+state+"/"+key)
                if mu:
                    eps=max(eps,residual/scale)
            out[key]={"C_v_km_s_Mpc3":Cv_original,"P_dd_Mpc3":Pdd,
                      "P_vv_km_s2_Mpc3_at_abs_mu_1":Pvv,
                      "maximum_covariance_determinant_relative_residual":eps}
        report[state]=out
    return report

def reject(name,func):
    try:
        func()
    except (ValueError,KeyError,TypeError):
        print("E26_NEGATIVE_REJECT",name,flush=True)
    else:
        raise AssertionError("E26 accepted incorrect source/claim "+name)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    x=load()
    report=replay(x)
    bad=copy.deepcopy(x);bad["e21"]["original_long_modes"]["0.003"]["states"]["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]*=1.001
    reject("E21_CV_TAMPER",lambda:replay(bad))
    bad=copy.deepcopy(x);bad["e13_FD"]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"]["0.005"]["W_R16_original_locked_top_hat"]*=.9
    reject("E13_FILTER_TAMPER",lambda:replay(bad))
    bad=copy.deepcopy(x);bad["e12_plus"]["density_proxy_two_fixed_steps_and_four_long_K"]["0.002"]["long_K"][1]["primordial_curvature_Delta2"]*=1.01
    reject("PRIMORDIAL_SOURCE_TAMPER",lambda:replay(bad))
    reject("WRONG_C_LIGHT",lambda:replay(x,c_light=1.))
    reject("WRONG_CONJUGATION_SIGN",lambda:replay(x,conjugation=-1))
    bad=copy.deepcopy(x);del bad["e13_minus"]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"]["0.005"]
    reject("MISSING_ORIGINAL_LONG_K",lambda:replay(bad))
    reject("FAKE_ARCHIVED_T_CDM",lambda:require(
        "t_cdm" in x["e12_FD"]["density_proxy_two_fixed_steps_and_four_long_K"]["0.002"]["long_K"][0],
        "E12 does not retain t_cdm long-K numerical transfer"))
    reject("FALSE_RECONSTRUCTED_INDEPENDENT_MODE",lambda:require(
        abs(report["FD"]["0.002"]["maximum_covariance_determinant_relative_residual"])>TOL,
        "same linear adiabatic source covariance is rank one"))
    print("E26_EIGHT_NEW_NEGATIVE_CONTROLS_PASS",flush=True)
    print("E26_ORIGINAL_E12_E13_E21_THREE_STATES_FOUR_K_SOURCE_RANK_PASS",flush=True)
    for s in STATES:
        print("E26_ORIGINAL_STATE",s,"MAX_RELATIVE_DET",
              max(v["maximum_covariance_determinant_relative_residual"]
                  for v in report[s].values()),flush=True)
    out={"stage":"E26_SINGLE_ADIABATIC_SOURCE_RANK",
         "status":"ORIGINAL_SOURCE_SINGLE_ADIABATIC_DENSITY_RELATIVE_VELOCITY_COVARIANCE_RANK_ONE_PASS",
         "protocol_git_blob":PINS["protocol"][1],
         "parent_source_blobs":{k:v[1] for k,v in PINS.items() if k!="protocol"},
         "original_state_K_covariance":report,
         "negative_controls_passed":8,
         "original_t_cdm_individual_long_K_transfer_not_archived":True,
         "no_standard_Doppler_vs_neutrino_template_numeric_rank":True,
         "not_full_eBOSS_galaxy_wake_or_24D_covariance":True,
         "observed_odd_sealed_no_new_CLASS_or_WSL":True}
    if args.out:
        require(args.out.suffix==".json" and not args.out.exists(),
                "output only new JSON")
        with args.out.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
        print("E26_SMALL_SOURCE_ONLY_JSON",args.out,flush=True)
    print("E26_NO_EXTRA_INDEPENDENT_PHASE_FROM_DETERMINISTIC_DENSITY_RECONSTRUCTION",flush=True)

if __name__=="__main__":main()
