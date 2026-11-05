#!/usr/bin/env python3
"""Retrospective mathematical E19 QA on pinned original E10 JSON, no data IO."""
from __future__ import annotations
import copy,hashlib,json,math
from pathlib import Path
R=Path(__file__).resolve().parents[1]/"source_data"
NAMES=[
("eboss_dr16_a03_e10_frozen_conditional_angular_wind_parity_summary_2026-09-27.json","a9e0f49a58b9fd5e0ad975d4d14ba6521a648c28"),
("eboss_dr16_a03_e10_wind_parity_angular_odd_bridge_protocol_2026-09-27.json","d9de48012103eaa07d531e11bca9cef285e62545"),
("eboss_dr16_a03_e6_physical_template_to_empirical_window_bridge_protocol_2026-09-27.json","4df20872cff5b3d8aa09831ab50b22fa5f8b6581"),
("eboss_dr16_a03_e19_posthoc_conditional_Fourier_shape_and_marginal_RR_nonidentifiability_2026-09-29.json","fc5044684b84ff3de86c236988617ead79d8e7e6")
]
def need(ok,msg):
    if not ok:raise ValueError("E19_FAIL_CLOSED: "+msg)
def close(a,b,msg):
    need(isinstance(a,(int,float)) and isinstance(b,(int,float)) and
         math.isfinite(a) and math.isfinite(b) and
         math.isclose(a,b,rel_tol=3e-13,abs_tol=3e-15),msg)
def load():
    out=[]
    for name,sha in NAMES:
        p=R/name
        need(p.is_file() and not p.is_symlink(),"missing exact original JSON")
        raw=p.read_bytes()
        need(hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()==sha,
             "immutable original Git blob drift: "+name)
        out.append(json.loads(raw))
    e,ep,b,r=out
    need(e["registered_protocol"]["git_blob_sha1"]==NAMES[1][1] and
         e["observed_odd_sealed"] and ep["guards"]["not_abs_eBOSS_xi"] and
         ep["guards"]["observed_odd_sealed"] and
         b["math_contract"]["pilot_joint_dimension"]==24 and
         b["physics_prediction_status"].startswith("NO_EBOSS_ABSOLUTE_PHYSICAL_WAKE") and
         r["classification"].startswith("EXPLICITLY_POSTHOC") and
         r["no_observed_odd_access"] and r["no_covariance_inverse_or_snr"],
         "original E10/E6 physics and E19 STOP mismatch")
    return e,r
def check(e,r):
    original=e["results_512_gauss_nodes"]
    t=original["conditional_positive_LOS_per_unit_delta_bias_P_cb"]
    v={}
    for label,upstream in (("FD","FD"),("Fplus","plus"),("Fminus","minus")):
        a=t[upstream]["ell1"];b=t[upstream]["ell3"]
        need(a!=0 and math.isfinite(a) and math.isfinite(b),"invalid original E10 Fourier source")
        reported=r["conditional_per_unit_deltaBias_Pcb"][label]
        close(reported["T1"],a,"T1 changed");close(reported["T3"],b,"T3 changed")
        close(reported["T3_over_T1"],b/a,"conditional Fourier shape ratio")
        v[label]=(a,b)
    for ell,j in (("ell1",0),("ell3",1)):
        close(original["plus_minus_conditional_difference"][ell],
              v["Fplus"][j]-v["Fminus"][j],"E10 original signed difference drift")
    need(original["equal_probability_plus_and_minus_wind_intrinsic_odd_mean_all_3_states_both_ells"]==0,
         "original E10 symmetric intrinsic mean changed")
    p,m=v["Fplus"],v["Fminus"]
    det=p[0]*m[1]-p[1]*m[0]
    angle=math.atan2(abs(det),p[0]*m[0]+p[1]*m[1])
    z=r["conditional_shape_algebra"]
    close(z["determinant_Fplus_Fminus_T1_T3"],det,"determinant")
    close(z["angle_between_conditional_2D_vectors_radians"],angle,"shape angle")
    close(z["absolute_ratio_difference_Fminus_minus_Fplus"],
          m[1]/m[0]-p[1]/p[0],"ratio gap")
    close(z["relative_ratio_difference_over_Fplus"],
          (m[1]/m[0]-p[1]/p[0])/(p[1]/p[0]),"relative ratio gap")
    close(z["quadrature_256_vs_512_max_relative_gap_from_frozen_E10"],
          original["max_relative_quadrature_256_vs_512_gap"],"frozen quadrature")
    need(det>0 and angle>0 and z["not_eBOSS_measurable_shape_or_sigma"],
         "mathematical conditional shape was recast as an observation")
    return v,det,angle
def witness(v,eta):
    need(-1<=eta<=1 and math.isfinite(eta),"invalid sign-weight/negative W")
    R0=7.25 # Positive SYNTHETIC marginal, never observed eBOSS randoms.
    wp,wm=R0*(1+eta),R0*(1-eta)
    need(wp>=0 and wm>=0,"latent selection weight is negative")
    close((wp+wm)/2,R0,"same fixed marginal RR")
    for t in v.values():
        for x in t:close((wp*x-wm*x)/(wp+wm),eta*x,"source mean not eta*T")
def reject(name,f):
    try:f()
    except (ValueError,KeyError,TypeError):print("E19_NEGATIVE_REJECT",name,flush=True)
    else:raise AssertionError("synthetic negative accepted "+name)
def main():
    e,r=load();v,det,angle=check(e,r)
    for eta in (-1,-.8,-.3,0,.3,.8,1):witness(v,eta)
    # Same fixed linear toy window maps opposite conditional signs to zero.
    M=((.7,.1),(.2,.6))
    for t in v.values():
        for row in M:
            a=sum(row[j]*t[j] for j in range(2))
            b=sum(row[j]*(-t[j]) for j in range(2))
            close((a+b)/2,0,"shared linear M cannot create signed mean")
    bad=copy.deepcopy(r);bad["conditional_shape_algebra"]["angle_between_conditional_2D_vectors_radians"]=1
    reject("TAMPERED_ANGLE",lambda:check(e,bad))
    bad=copy.deepcopy(e)
    bad["results_512_gauss_nodes"]["plus_minus_conditional_difference"]["ell1"]=0
    reject("TAMPERED_FROZEN_E10_DIFF",lambda:check(bad,r))
    bad=copy.deepcopy(e)
    bad["results_512_gauss_nodes"]["equal_probability_plus_and_minus_wind_intrinsic_odd_mean_all_3_states_both_ells"]=1
    reject("FORGED_UNCONDITIONAL_SOURCE",lambda:check(bad,r))
    reject("NEGATIVE_CONDITIONAL_WEIGHT",lambda:witness(v,1.1))
    print("E19_FOUR_SYNTHETIC_NEGATIVES_PASS",flush=True)
    print("E19_FPLUS_T3_T1",format(v["Fplus"][1]/v["Fplus"][0],".16g"),flush=True)
    print("E19_FMINUS_T3_T1",format(v["Fminus"][1]/v["Fminus"][0],".16g"),flush=True)
    print("E19_FPLUS_FMINUS_DET",format(det,".16g"),"ANGLE_RAD",format(angle,".16g"),flush=True)
    print("E19_SAME_MARGINAL_RR_DIFFERENT_LATENT_EV_MEANS_PROVEN_TOY",flush=True)
    print("E19_NO_NEW_DATA_NO_PHYSICAL_EBOSS_CLAIM_OBSERVED_ODD_SEALED",flush=True)
if __name__=="__main__":main()
