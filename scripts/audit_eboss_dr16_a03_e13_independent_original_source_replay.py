#!/usr/bin/env python3
"""Independent, retrospective replay of E13 original-source Gamma and product.

No import of E13/E12/E11 executable functions. Uses original archived byte
SHA sources and separate formulas for Gaussian/Legendre quadrature, R16
top-hat numerator, E12 per-state normalization and B-source multiplication.
Pure source-only certificate, not independent N-body/tracer/survey simulation.
"""
from __future__ import annotations
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
HISTORY=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27"
E12=ROOT/"source_data/eboss_dr16_a03_e12_archived_CI_2026_09_27"
E13=ROOT/"source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27"
PARENT={
 "E8":"bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0",
 "E9":"450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890",
 "E13":"537365b2e758ec215810bb687557771440851cf2869d02578569b1d47d4ec98d",
 "FD_E12":"b1c679e8ccdf2e901abbb2ffc200cfcce765bf8ee9492364a280625835d42414",
 "plus_E12":"aedaf4381065ba807adced01b7740ade91db5a83d49113315e7a596ac2652655",
 "minus_E12":"4a12370b578a5482b67cd3d484504bf672396df35d8c252ff5d1150875b04777",
 "FD_E13":"f32527307024480e031e3d7952caf681b577f528d70b86fae1efa4ee7828bfde",
 "plus_E13":"40d0da5e772ee99dce809a48179064a964ae20f46876dacde8d529158a8c7dd2",
 "minus_E13":"bf94d3b701867492f91768b182839b149de68cce147de2d5b21a813c2d907a7a"}
STATES=("FD","plus","minus")
LONG=(.001,.002,.003,.005)
C=299792.458
M=.06
TN=.71611*2.7255*8.617333262e-5

def require(x,msg):
    if not x:raise ValueError(msg)
def sha(data):return hashlib.sha256(data).hexdigest()
def load(path,key):
    raw=path.read_bytes()
    require(sha(raw)==PARENT[key],"Original independent-audit source SHA fail "+str(path))
    if path.suffix==".csv":return np.genfromtxt(path,delimiter=",",names=True)
    return json.loads(raw)

def occupation(src,s,x):
    x=np.asarray(x,float)
    q=np.asarray(src["q_dimensionless"],float)
    require(np.isfinite(x).all() and np.min(x)>=q[0] and np.max(x)<=q[-1],
            "E8 original source extrapolation in independent replay")
    f0=1./(1.+np.exp(x))
    if s=="FD":return f0
    cols={"plus":"Fplus_CLASS_normalized","minus":"Fminus_CLASS_normalized"}
    return (f0*np.interp(x,q,np.asarray(src[cols[s]],float))/
               np.interp(x,q,np.asarray(src["F0_CLASS_normalized"],float)))

def direct_gamma_replay(src,old9,direct,s):
    r,w=hermgauss(96);r=np.sqrt(2.)*r;w=w/math.sqrt(math.pi)
    mu,wm=leggauss(512)
    accum={}
    for z in (.945,.955):
        original=old9[(.004,z)]["state"][s]
        vold=original["own_CLASS_v_R16_1sigma_los_kms"]
        vnew=direct[s]["direct_R16_sigma_LOS_kms_three_fixed_z"][str(z)]
        qa=M*vnew/(C*TN*(1+z))
        qb=M*vold/(C*TN*(1+z))
        Aold=original["alpha_with_own_CLASS_velocity"]
        Anew=Aold*(vnew/vold)*float(occupation(src,s,qa)/occupation(src,s,qb))
        rmu=r[:,None]*mu[None,:]
        # The rank/angle-dependence is independently reconstructed here
        # without importing the E13 phase or projection implementation.
        accum[z]=(Anew*rmu*occupation(src,s,qa*np.abs(rmu))/
                  float(occupation(src,s,qa)))
    dlog=-(1+.95)*(accum[.955]-accum[.945])/.01
    z=mu[None,:]**2*dlog
    for ell in (1,3):
        pp=mu if ell==1 else .5*(5*mu**3-3*mu)
        t=(2*ell+1)/2.*np.sum(wm[None,:]*pp[None,:]*z,axis=1)
        gamma=float(np.sum(w*r*t))
        mean=float(np.sum(w*t))
        require(abs(mean)<4e-14,"A spurious unconditional odd mean survived sign symmetry")
        accum[str(ell)]={"Gamma_independent":gamma,
                         "mean_odd":mean,
                         "source_rank_96_angular_512":True}
    return {k:accum[k] for k in ("1","3")}

def Wtop(k):
    x=k*16.
    if abs(x)<1e-5:return 1.
    return 3*(math.sin(x)-x*math.cos(x))/x**3

def main():
    src=load(HISTORY/"E8/frozen_matched_distributions_4000q.csv","E8")
    e9raw=load(HISTORY/"E9/quasistatic_alpha_dln_a_source_only_plus_CLASS.json","E9")
    e13=load(E13/"e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json","E13")
    assert e13["status"].startswith("E13_RANK_MATCHED_DIRECT_VTK_")
    assert e13["observed_galaxy_random_or_odd_read"] is False
    assert e13["original_E11_Gamma_IS_NOT_USED_IN_FINAL_Bsource"] is True
    original={(q["wind_velocity_derivative_relative_z_step"],q["z"]):q
              for q in e9raw["all_own_and_common_FD_wind_cases"]}
    d12={s:load(E12/("e12_direct_vTk_long_cross_"+s+".json"),s+"_E12")
         for s in STATES}
    d13={s:load(E13/("e13_direct_wind_"+s+".json"),s+"_E13")
         for s in STATES}
    old=e13["conditional_Bsource_unit_DeltaBias_per_E8_state"]
    numerical={}
    all_checks=[]
    for s in STATES:
        independent=direct_gamma_replay(src,original,d13,s)
        numerical[s]={"independent_Gamma":independent,"longK":{}}
        for ell in ("1","3"):
            measured=old[s]["Gamma_direct_ell"+ell]
            derived=independent[ell]["Gamma_independent"]
            err=abs(measured-derived)
            require(err/max(abs(measured),abs(derived),1e-14)<1e-10,
                    "Independent original E8/E9/E13 Gamma failed "+s+ell)
            all_checks.append(err)
        for index,K in enumerate(LONG):
            tag=str(K)
            l13=d13[s]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"][tag]
            l12=d12[s]["density_proxy_two_fixed_steps_and_four_long_K"]["0.004"]["long_K"][index]
            require(l12["K_comoving_h_per_Mpc"]==K,
                    "Original E12 long-K geometry changed")
            sigma_old=d12[s]["direct_R16_sigma_LOS_kms"]
            sigma_new=d13[s]["direct_R16_sigma_LOS_kms_three_fixed_z"]["0.95"]
            corrected=Wtop(K)*l12["P_delta_cb_rdirect_LOS_over_i_mu_Mpc3"]*sigma_old/sigma_new
            reported=l13["P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3"]
            require(abs(corrected-reported)/max(abs(corrected),abs(reported),1e-14)<1e-3,
                    "Original E12b-R16-filtered long direct rank mismatch")
            pk=d13[s]["CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95"]
            perK={}
            for ell in ("1","3"):
                B=independent[ell]["Gamma_independent"]*corrected*pk
                key="ell"+ell+"_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"
                reportedB=old[s]["long_K_fixed"][tag][key]
                err=abs(B-reportedB)/max(abs(B),abs(reportedB),1e-14)
                require(err<1e-3,"Independent Bsource factor replay failed "+s+tag+ell)
                perK[ell]={"Bsource_unit_bias_independent_Mpc6":B,
                           "original_E13_Bsource_Mpc6":reportedB,
                           "relative_product_replay_gap":err}
                all_checks.append(err)
            numerical[s]["longK"][tag]={"W16_direct_calculated":Wtop(K),
                "corrected_P_Dr_independent_Mpc3":corrected,
                "original_E13_P_Dr_Mpc3":reported,
                "unit_bias_Bsrc":perK}
    differences={}
    for K in LONG:
        tag=str(K)
        for ell in ("1","3"):
            v=(numerical["plus"]["longK"][tag]["unit_bias_Bsrc"][ell]
              ["Bsource_unit_bias_independent_Mpc6"]-
              numerical["minus"]["longK"][tag]["unit_bias_Bsrc"][ell]
              ["Bsource_unit_bias_independent_Mpc6"])
            orig=e13["conditional_Fplus_minus_Fminus_Bsource_unit_DeltaBias"][tag][
                 "delta_ell"+ell+"_unit_DeltaBias_Bsource_over_i_muLong_Mpc6"]
            require(abs(v-orig)/max(abs(v),abs(orig),1e-14)<1e-3,
                    "Independent Fplus-Fminus cross contrast does not match original")
            differences[tag+"_ell"+ell]={"original":orig,"independent":v}
    output={
      "date":"2026-09-27",
      "status":"INDEPENDENT_FROM_E13_EXECUTABLE_EXACT_ORIGINAL_E8_E9_E12_E13_GAMMA_W16_LONG_CROSS_AND_BSOURCE_REPLAY_PASS",
      "original_E13_SHA256":PARENT["E13"],
      "original_E8_E9_E12_E13_source_SHA256":PARENT,
      "replay_scope":"Retrospective separate-script 96GHx512GL Gamma from archived original E8+E9 plus direct vTk sigma z from E13, exact long W16 correction of original E12 raw cross, short CLASS P_cb and fixed-K plus/minus Bsource arithmetic; NOT a second CLASS ensemble nor eBOSS sky bispectrum.",
      "original_z_and_k":{"z":.95,"k_short_COM_h_Mpc":.05,
                         "long_K_COM_h_Mpc":list(LONG),"R_Mpc_over_h":16},
      "states":numerical,"original_plus_minus_fixed_K_contrasts":differences,
      "max_absolute_Gamma_or_fractional_product_replay_error":max(all_checks),
      "technical_replay_full_source_no_new_CLASS_or_FITS":True,
      "observed_odd_sealed":True,"original_parent_reports_unchanged":True,
      "main_untouched":True,
      "not_full_observable_eBOSS_bispectrum_or_significance":True}
    dest=ROOT/"eboss_workspace/a03_physics_source/e13_independent_source_replay.json"
    dest.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(output,indent=2,allow_nan=False)+"\n").encode()
    if dest.exists():require(dest.read_bytes()==raw,"Existing E13 independent replay changed")
    else:
        with tempfile.NamedTemporaryFile(dir=dest.parent,prefix=".e13_replay_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,dest)
        finally:t.unlink(missing_ok=True)
    print("E13_INDEPENDENT_DIRECT_GAMMA_AND_FILTERED_LONG_BSOURCE_REPLAY_PASS",
          "REPORT_SHA256",sha(raw),
          "MAX_ABS_GAMMA_OR_REL_PRODUCT_DIFF",max(all_checks),
          "NO_SECOND_CLASS NO_OBSERVED NO_EBOSS_BISPECTRUM",flush=True)

if __name__=="__main__":
    main()
