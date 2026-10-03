#!/usr/bin/env python3
"""E13 original-F± direct-vTk rank-matched Gaussian long-short source.

Three separately evolved pinned CLASS states yield direct neutrino-CDM theta
R16-smoothed sigma(z=.945,.95,.955), NOT old density-derived E9 proxy sigma.
New direct-Gaussian Gamma is calculated under the same filtered LOS rank as
E12b-corrected long P_delta_cb,r. Then multiply original CLASS P_cb(short),
unit Delta_b, and the restricted Gaussian local EqA8 source. Not a complete
or observable eBOSS galaxy bispectrum, not the old 24D odd mean, no FITS.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss,legval

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e12_direct_class_vTk_long_cross as e12

PROTO=ROOT/"source_data/eboss_dr16_a03_e13_direct_vTk_rank_matched_squeezed_source_prereg_2026-09-27.json"
P_BLOB="3c0e2031a14b5bb07e41dce29005e07dc1439deb"
E12B=ROOT/"source_data/eboss_dr16_a03_e12b_filtered_long_rank_cross_window_erratum_2026-09-27.json"
E12B_BLOB="09bac4f72e325578cbed7e790c50114150087e89"
E12_ARCH=ROOT/"source_data/eboss_dr16_a03_e12_archived_CI_2026_09_27"
E12_MANIFEST_SHA="6022a4df6b7f54a4f01bf007d200897769cef680"
STATES=e12.STATES
ZS=(.945,.95,.955)
KSHORT=.05
K_LONG=e12.LONGK
CLASS_SHA=e12.CLASS_SHA
R16=e12.R16
KVEL=e12.KVEL
C=e12.C
TNU0=.71611*2.7255*8.617333262e-5
GH=(48,96)
GL=(256,512)
E12_EXPECT={
 "FD":"b1c679e8ccdf2e901abbb2ffc200cfcce765bf8ee9492364a280625835d42414",
 "plus":"aedaf4381065ba807adced01b7740ade91db5a83d49113315e7a596ac2652655",
 "minus":"4a12370b578a5482b67cd3d484504bf672396df35d8c252ff5d1150875b04777"}
OUT=ROOT/"eboss_workspace/a03_physics_source/e13_direct_rank_matched_source"

def check(condition,why):
    if not condition:raise ValueError(why)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def stable(data):
    return (json.dumps(data,indent=2,allow_nan=False)+"\n").encode()

def write_once(path,data):
    raw=stable(data)
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        check(path.read_bytes()==raw,"E13 result file exists with DIFFERENT bytes")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e13_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E13_SAVED",path,"SHA256",sha(raw),flush=True)

def gate():
    p0,c,e9,_=e12.original_inputs()
    check(e12.blob(PROTO.read_bytes())==P_BLOB and
          e12.blob(E12B.read_bytes())==E12B_BLOB and
          e12.blob((E12_ARCH/"archive_manifest.json").read_bytes())==E12_MANIFEST_SHA,
          "Frozen E13 prospective scope, original E12, or E12b analytic filter erratum changed")
    p=json.loads(PROTO.read_bytes())
    fix=p["fixed_inputs"]
    check(fix["short_k_COMOVING_h_Mpc"]==KSHORT
          and np.array_equal(K_LONG,np.array(fix["long_K_COMOVING_h_Mpc"],float))
          and fix["source_phase_z_nodes"]==list(ZS)
          and fix["R16_filter_comoving_Mpc_over_h"]==R16
          and p["physical_limits"]["observed_galaxy_random_odd_sealed"] is True
          and p["corrected_long_cross"]["mandatory_E12b_filter"].startswith("P_"),
          "E13 K/z/W observation scope changed")
    e12m=json.loads((E12_ARCH/"archive_manifest.json").read_bytes())
    names={x["file"]:x["sha256"] for x in e12m["files"]}
    old={}
    for state in STATES:
        path=E12_ARCH/("e12_direct_vTk_long_cross_"+state+".json")
        required=E12_EXPECT[state]
        check(names.get(path.name)==required and sha(path.read_bytes())==required,
              "Original E12 direct-vTk per-state original SHA changed")
        obj=json.loads(path.read_bytes())
        check(obj["CLASS_commit"]==CLASS_SHA
              and obj["state"]==state
              and obj["E11_Gamma_multiplied_into_real_P_Dr"] is False,
              "E12 source-only physical STOP changed")
        old[state]=obj
    return p,c,e9,old

def source_one(c,state,q):
    """Exact E8 relative occupation, WITHOUT CLASS 2/(2pi)^3 prefactor."""
    q=np.asarray(q,float)
    grid=np.asarray(c["q_dimensionless"],float)
    check(np.isfinite(q).all() and
          np.min(q)>=grid[0] and np.max(q)<=grid[-1],
          "Frozen E8 4000-q source extrapolation forbidden")
    fd=1./(1.+np.exp(q))
    if state=="FD":return fd
    return (fd*np.interp(q,grid,np.asarray(c[e12.COLS[state]],float))/
                np.interp(q,grid,np.asarray(c[e12.COLS["FD"]],float)))

def rescale_original_E9_to_direct_alpha(alpha_old,proxy_v,direct_v,state,z,c):
    check(alpha_old>0 and proxy_v>0 and direct_v>0,
          "Original E9 model or direct CLASS v became degenerate")
    qp=.06*proxy_v/(C*TNU0*(1.+z))
    qd=.06*direct_v/(C*TNU0*(1.+z))
    f0=float(source_one(c,state,np.asarray(qp)))
    fd=float(source_one(c,state,np.asarray(qd)))
    return float(alpha_old*(direct_v/proxy_v)*(fd/f0))

def run_state(state,out):
    check(state in STATES,"Only original frozen F0/F±")
    p,c,e9,old=gate()
    e12.algebraic_test()
    from classy import Class
    import class_response_optimize as cro
    import wake_two_tracer_fisher as base
    check(cro.CLASS_COMMIT==CLASS_SHA and
          np.array_equal(base.KVEL,KVEL) and
          base.R_FILTER==R16,"Original CLASS or R16 changed")
    psd=out/("frozen_E13_"+state+".dat")
    psd.parent.mkdir(parents=True,exist_ok=True)
    q=np.asarray(c["q_dimensionless"],float)
    f=np.asarray(c[e12.COLS[state]],float)
    txt="\n".join(f"{a:.14e} {bb:.14e}" for a,bb in zip(q,f))+"\n"
    if psd.exists():check(psd.read_text()==txt,"E13 original PSD file overwritten")
    else:psd.write_text(txt)
    base.ZBINS=np.array([.9,1.],float)
    params=base.class_params(psd,.06)
    params["output"]="mPk,dTk,vTk"
    refparams=old[state]["CLASS_output_direct_theta"]
    check(refparams["output"]=="mPk,dTk,vTk"
          and refparams["gauge"]=="newtonian"
          and params["gauge"]=="newtonian"
          and math.isclose(params["H0"],67.36)
          and params["z_max_pk"]>=1.05,
          "Original E12/E9 CLASS parameter conventions changed")
    cosm=Class();cosm.set(params);cosm.compute()
    try:
        hh=float(cosm.h())
        check(math.isclose(hh,.6736,rel_tol=0.,abs_tol=1e-14),
              "Original F± h changed")
        kin=None
        sig={}
        tcenter=None
        short_p=None
        for z in ZS:
            tk=cosm.get_transfer(z=float(z),output_format="class")
            if kin is None:
                kin=np.asarray(tk["k (h/Mpc)"],float)
                check(np.all(np.diff(kin)>0)
                      and kin[0]<=KVEL[0] and kin[-1]>=KVEL[-1]
                      and {"t_ncdm[0]","t_cdm","d_b","d_cdm"}.issubset(tk),
                      "Original CLASS direct theta + cb density not available on full K grid")
            check(np.array_equal(kin,np.asarray(tk["k (h/Mpc)"],float)),
                  "CLASS K grid changed between original phase nodes")
            th=np.interp(KVEL,kin,np.asarray(tk["t_ncdm[0]"],float)-
                         np.asarray(tk["t_cdm"],float))
            kt=KVEL*hh
            vel=-C*th/kt*np.sqrt(e12.primordial_delta2(kt,cro.A_S,base.NS))
            val=float(math.sqrt(e12.var_to_cutoff(vel)/3.))
            check(val>0 and math.isfinite(val),"Direct R16 vTk variance invalid")
            sig[str(z)]=val
            if z==.95:
                tcenter=tk
                short_p=float(cosm.pk_cb_lin(KSHORT*hh,.95))
                check(short_p>0 and math.isfinite(short_p),
                      "Original state CLASS P_cb_short nonphysical")
        old12=old[state]
        rel=abs(sig["0.95"]-old12["direct_R16_sigma_LOS_kms"])/old12["direct_R16_sigma_LOS_kms"]
        tol=p["prospective_engineering_QA_only"]["original_E12_center_sigma_relative_tolerance"]
        check(rel<tol,
              "E13 fresh exact direct vTk sigma differs from E12 beyond preregistration")
        newlong={}
        kval=np.asarray(tcenter["k (h/Mpc)"],float)
        td=np.asarray(tcenter["t_ncdm[0]"],float)-np.asarray(tcenter["t_cdm"],float)
        cb=((cro.OMEGA_B*np.asarray(tcenter["d_b"],float)+
             cro.OMEGA_CDM*np.asarray(tcenter["d_cdm"],float))/
             (cro.OMEGA_B+cro.OMEGA_CDM))
        oldlong=old12["density_proxy_two_fixed_steps_and_four_long_K"]["0.004"]["long_K"]
        for i,k in enumerate(K_LONG):
            th=float(np.interp(k,kval,td))
            dcb=float(np.interp(k,kval,cb))
            orig=oldlong[i]
            check(orig["K_comoving_h_per_Mpc"]==k,"E12 original K order changed")
            tdiff=abs(th-orig["theta_relative_direct_vTk_per_primordial_curvature_Mpc_inv"])/max(
                abs(th),abs(orig["theta_relative_direct_vTk_per_primordial_curvature_Mpc_inv"]),1e-14)
            ddiff=abs(dcb-orig["delta_cb_Newtonian_per_primordial_curvature"])/max(
                abs(dcb),abs(orig["delta_cb_Newtonian_per_primordial_curvature"]),1e-14)
            check(max(tdiff,ddiff)<p["prospective_engineering_QA_only"]["original_E12_center_direct_theta_deltaCB_perK_relative_tolerance"],
                  "Fresh E13 CLASS long transfer differs from E12 original beyond locked QA")
            W=float(e12.tophat(np.asarray([k*R16]))[0])
            Krealk=float(k*hh)
            unfiltered=e12.cross_i_mu_real(dcb,th,Krealk,sig["0.95"],cro.A_S,base.NS)
            corrected=W*unfiltered
            reconstructed_from_E12=W*orig["P_delta_cb_rdirect_LOS_over_i_mu_Mpc3"]*(
                old12["direct_R16_sigma_LOS_kms"]/sig["0.95"])
            closediff=abs(corrected-reconstructed_from_E12)/max(
                abs(corrected),abs(reconstructed_from_E12),1e-14)
            check(closediff<p["prospective_engineering_QA_only"]["original_E12_center_direct_theta_deltaCB_perK_relative_tolerance"],
                  "New direct filtered rank long cross fails locked old E12b validation")
            newlong[str(k)]={"k_com_h_Mpc":float(k),
                "W_R16_original_locked_top_hat":W,
                "E13_delta_cb_Newtonian_transfer_per_R":dcb,
                "E13_theta_rel_direct_per_R_Mpc_inv":th,
                "P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3":corrected,
                "old_E12_raw_P_delta_cb_r_over_i_mu_Mpc3":orig["P_delta_cb_rdirect_LOS_over_i_mu_Mpc3"],
                "filter_and_sigma_corrected_original_E12_comparison":reconstructed_from_E12,
                "fresh_theta_vs_E12_relative_gap":tdiff,
                "fresh_dcb_vs_E12_relative_gap":ddiff,
                "corrected_Plong_relative_closure_gap":closediff}
        result={
          "status":"E13_DIRECT_CLASS_VTK_R16_WIND_THREE_REDSHIFT_NODES_AND_FILTERED_LONG_CROSS_COMPLETE",
          "state":state,"CLASS_commit":CLASS_SHA,"preregistration_git_blob_sha1":P_BLOB,
          "original_4000q_E8_csv_sha256":e12.ORIG_E8,
          "original_E9_three_CLASS_sha256":e12.ORIG_E9,
          "original_E12_per_state_sha256":E12_EXPECT[state],
          "original_E12b_filter_erratum_git_blob":E12B_BLOB,
          "reconstructed_PSD_from_original_archived_E8_SHA256":sha(psd.read_bytes()),
          "direct_R16_sigma_LOS_kms_three_fixed_z":sig,
          "CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95":short_p,
          "sigma_center_vs_original_E12_relative_gap":rel,
          "fresh_direct_filtered_Plong_over_i_mu_fixed_long_K":newlong,
          "observed_galaxy_random_or_odd_read":False,
          "actual_complete_eBOSS_LRG_ELG_bispectrum_calculated":False}
        write_once(out/("e13_direct_wind_"+state+".json"),result)
        print("E13_DIRECT_VTK_STATE",state,"SIGMA_0945_095_0955",
              json.dumps(sig,sort_keys=True),"P_CB_SHORT",short_p,
              "E12_CENTER_SIGMA_GAP",rel,flush=True)
    finally:
        cosm.struct_cleanup();cosm.empty()
    return 0

def gamma_for(s,state,q,c,e9,wind,ell,mu,r):
    """Exact original E8/E9 occupancy, now with vTk rank from E13."""
    vals={}
    for z in (.945,.955):
        key=str(z)
        vnew=wind[state]["direct_R16_sigma_LOS_kms_three_fixed_z"][key]
        original=e9[(.004,z)]["state"][state]
        vold=original["own_CLASS_v_R16_1sigma_los_kms"]
        alpha_original=original["alpha_with_own_CLASS_velocity"]
        alpha_direct=rescale_original_E9_to_direct_alpha(
            alpha_original,vold,vnew,state,z,c)
        qq=.06*vnew/(C*TNU0*(1.+z))
        den=float(source_one(c,state,np.asarray(qq)))
        phase=alpha_direct*np.asarray(r)[:,None]*np.asarray(mu)[None,:]*(
            source_one(c,state,qq*np.abs(np.asarray(r)[:,None]*np.asarray(mu)[None,:]))/den)
        vals[key]=phase
    dlog=-(1.+.95)*(vals["0.955"]-vals["0.945"])/(.01)
    theta=np.asarray(mu)[None,:]**2*dlog
    w=s["quad_weights"]
    p=np.zeros(ell+1);p[-1]=1.
    return (2*ell+1)/2.*np.sum(
        w[None,:]*legval(mu,p)[None,:]*theta,axis=1)

def direct_gamma_for(s,state,c,e9,wind,nmu,ngh):
    mu,muw=leggauss(nmu)
    rr,rw=hermgauss(ngh)
    rr=np.sqrt(2.)*rr
    rw=rw/math.sqrt(math.pi)
    check(abs(float(np.sum(rw))-1.)<3e-14 and
          abs(float(np.sum(rw*rr)))<3e-14 and
          abs(float(np.sum(rw*rr*rr))-1.)<3e-13,
          "Original unit-normal Gaussian rank quadrature invalid")
    s={"quad_weights":muw}
    result={}
    for ell in (1,3):
        t=gamma_for(s,state,None,c,e9,wind,ell,mu,rr)
        check(np.max(np.abs(t+t[::-1]))<3e-14,
              "Direct rank conditional source lost original odd symmetry")
        mean=float(np.sum(rw*t))
        gamma=float(np.sum(rw*rr*t))
        rms=float(np.sqrt(np.sum(rw*t*t)))
        check(abs(mean)<5e-14 and math.isfinite(gamma) and rms>0,
              "Spurious symmetric unconditional 2point mean or degenerate Gamma")
        result[str(ell)]={"Gamma_direct_R16_rank":gamma,
            "unconditional_intrinsic_mean_odd_symmetric_rank":mean,
            "conditional_source_rms":rms}
    return result

def aggregate(out):
    p,c,e9full,e12old=gate()
    e9={(row["wind_velocity_derivative_relative_z_step"],row["z"]):row
        for row in e9full["all_own_and_common_FD_wind_cases"]}
    check(all((.004,z) in e9 for z in ZS),"Original E9 phase stencil incomplete")
    wind={}
    reports={}
    for state in STATES:
        fp=out/("e13_direct_wind_"+state+".json")
        raw=fp.read_bytes();j=json.loads(raw)
        check(j["state"]==state and j["preregistration_git_blob_sha1"]==P_BLOB
              and j["observed_galaxy_random_or_odd_read"] is False,
              "E13 original state worker missing or observed seal changed")
        wind[state]=j;reports[state]=sha(raw)
    results={}
    for nmu in GL:
        for ngh in GH:
            values={}
            for state in STATES:
                values[state]=direct_gamma_for({},state,c,e9,wind,nmu,ngh)
            results[f"GL_{nmu}_GH_{ngh}"]=values
    best=results["GL_512_GH_96"]
    qa={}
    for state in STATES:
        for ell in ("1","3"):
            ga=best[state][ell]["Gamma_direct_R16_rank"]
            g48=results["GL_512_GH_48"][state][ell]["Gamma_direct_R16_rank"]
            g256=results["GL_256_GH_96"][state][ell]["Gamma_direct_R16_rank"]
            def gap(a,b):return float(abs(a-b)/max(abs(a),abs(b),1e-14))
            eg=gap(ga,g48);em=gap(ga,g256)
            qa[state+"_ell"+ell+"_GH48_96"]=eg
            qa[state+"_ell"+ell+"_GL256_512"]=em
            check(eg<p["prospective_engineering_QA_only"]["max_GH48_96_relative_gap"]
                  and em<p["prospective_engineering_QA_only"]["max_GL256_512_relative_gap"],
                  "Original preregistered E13 Hermite/Legendre QA failed")
    # Independent algebraic normal-Gaussian Stein identity via direct
    # derivative at two original E11 rank steps, all original F states.
    mu,muw=leggauss(512)
    r,w=hermgauss(96);r=np.sqrt(2.)*r;w=w/math.sqrt(math.pi)
    for state in STATES:
        for ell in (1,3):
            ga=best[state][str(ell)]["Gamma_direct_R16_rank"]
            for h in (.002,.001):
                a=gamma_for({"quad_weights":muw},state,None,c,e9,wind,ell,mu,r+h)
                b=gamma_for({"quad_weights":muw},state,None,c,e9,wind,ell,mu,r-h)
                eq=float(np.sum(w*(a-b)/(2*h)))
                err=abs(ga-eq)/max(abs(ga),abs(eq),1e-14)
                qa[f"Stein_{state}_ell{ell}_step{h}"]=err
                check(err<p["prospective_engineering_QA_only"]["max_Stein_step_Gaussian_identity_relative_gap"],
                      "Predeclared E13 Gaussian rank Stein check failed")
    old11=json.loads((e12.ARCH/"E11/a03_e11_gaussian_squeezed_source_only_report.json").read_bytes())
    oldgamma=old11["four_fixed_quadrature_results"]["mu_512_gh_96"]
    kernel={}
    for state in STATES:
        pk=wind[state]["CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95"]
        allK={}
        for k in K_LONG:
            long=wind[state]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"][str(float(k))]
            pd=long["P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3"]
            allK[str(float(k))]={
                "K_comoving_h_per_Mpc":float(k),
                "P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3":pd,
                "P_cb_short_CLASS_Mpc3":pk,
                "ell1_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6":
                    best[state]["1"]["Gamma_direct_R16_rank"]*pd*pk,
                "ell3_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6":
                    best[state]["3"]["Gamma_direct_R16_rank"]*pd*pk,
                "full_eBOSS_galaxy_bispectrum_inferred":False}
        kernel[state]={"Gamma_direct_ell1":best[state]["1"]["Gamma_direct_R16_rank"],
                       "Gamma_direct_ell3":best[state]["3"]["Gamma_direct_R16_rank"],
                       "old_E11_density_proxy_Gamma_ell1":oldgamma[state]["1"]["Gamma_E_rank_times_conditional_T"],
                       "old_E11_density_proxy_Gamma_ell3":oldgamma[state]["3"]["Gamma_E_rank_times_conditional_T"],
                       "Pcb_short_Mpc3":pk,"long_K_fixed":allK}
    differences={}
    for k in K_LONG:
        tag=str(float(k))
        a=kernel["plus"]["long_K_fixed"][tag]
        b=kernel["minus"]["long_K_fixed"][tag]
        differences[tag]={
          "delta_ell1_unit_DeltaBias_Bsource_over_i_muLong_Mpc6":
            a["ell1_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"]-
            b["ell1_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"],
          "delta_ell3_unit_DeltaBias_Bsource_over_i_muLong_Mpc6":
            a["ell3_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"]-
            b["ell3_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"]}
    outjson={
      "date":"2026-09-27",
      "status":"E13_RANK_MATCHED_DIRECT_VTK_GAUSSIAN_WICK_FIXED_K_SOURCE_ONLY_BISPECTRUM_KERNEL_COMPLETE_NOT_EBOSS_OBSERVABLE",
      "prospectively_registered_E13_protocol_git_blob":P_BLOB,
      "original_E12b_R16_long_mode_filter_erratum_git_blob":E12B_BLOB,
      "original_E8_CSV_SHA256":e12.ORIG_E8,
      "original_E9_three_CLASS_SHA256":e12.ORIG_E9,
      "original_E12_per_state_original_sha256":E12_EXPECT,
      "fresh_E13_direct_vTk_per_state_source_sha256":reports,
      "original_E11_Gamma_IS_NOT_USED_IN_FINAL_Bsource":True,
      "original_E11_Gamma_shown_only_as_density_proxy_comparison":True,
      "original_short_k_COMOVING_h_Mpc":KSHORT,
      "original_long_K_COMOVING_h_Mpc":K_LONG.tolist(),
      "math_midpoint_z_not_eBOSS_measured_z_eff":.95,
      "model":"Gaussian r~N(0,1) from R16-filtered DIRECT CLASS vTk LOS wind; exact frozen F± Eq20 occupation; same rank in original EqA8 angular Gaussian Gamma and in corrected long P_delta_cb,r with SAME W_R16(K); separate CLASS P_cb_short for each frozen F state, unit DeltaBias=1. Units Mpc^6, coefficient of +i mu_long.",
      "all_predeclared_four_quadrature_combinations_Gamma":results,
      "rank_match_and_technical_QA":qa,
      "technical_QA_warnings":[],
      "conditional_Bsource_unit_DeltaBias_per_E8_state":kernel,
      "conditional_Fplus_minus_Fminus_Bsource_unit_DeltaBias":differences,
      "physical_limits":["Not Zhu-Castorina exact 3point estimator or full Einstein-Vlasov retarded halo wake",
         "Only fixed short mode k=.05, 4 fixed long modes and redshift z=.95, Newtonian cb density. No full squeezed triangle angular orientation, unequal-time, gauge-invariant observed galaxy model or full-k/z integration.",
         "Published Eq20 exact FD occupancy versus approximate mu=1 normalization remains unresolved; no physical eBOSS LRG/ELG DeltaBias, HOD, magnification/evolution selection, or full GR bispectrum nuisance.",
         "No validated eBOSS triple window, independent bispectrum covariance or 24D odd observation. Conditional source kernel is NOT survey detection or exclusion."],
      "observed_galaxy_random_or_odd_read":False,
      "new_catalogue_or_mock_download_or_science_cut_seed":False,
      "main_untouched":True,"PR_remains_draft":True}
    write_once(out/"e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json",outjson)
    print("E13_FROZEN_DIRECT_VTK_RANK_MATCHED_GAUSSIAN_SOURCE_KERNEL_COMPLETE",
          "MAX_TECH_QA_GAP",max(qa.values()),
          "GAMMA",json.dumps({s:{ell:best[s][ell]["Gamma_direct_R16_rank"]
                                   for ell in ("1","3")} for s in STATES},sort_keys=True),
          "NO_REAL_EBOSS_BISPECTRUM",flush=True)
    print("E13_FPLUS_MINUS_FMINUS_FIXED_LONG_K_UNITBIAS_BSOURCE",
          json.dumps(differences,sort_keys=True),flush=True)
    return 0

def main():
    a=argparse.ArgumentParser(description=__doc__)
    g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--state",choices=STATES)
    g.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    args=a.parse_args()
    if args.self_test:
        gate();e12.algebraic_test();return 0
    if args.state:return run_state(args.state,args.output_dir)
    return aggregate(args.output_dir)

if __name__=="__main__":
    raise SystemExit(main())
