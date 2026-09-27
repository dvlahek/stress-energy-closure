#!/usr/bin/env python3
"""E12 original frozen F0/F+/- direct CLASS vTk vs archived E9 dTk velocity.

Source-only, two separate types of velocity:
 * direct CLASS t_ncdm[0]-t_cdm (theta divergence transfer in 1/Mpc);
 * prior E9 -H d(d_ncdm[0]-d_cdm)/dz density-derived proxy.
No fitted normalization, preferred wind sign, eBOSS galaxies, survey
xi/bispectrum, observed odd or new catalogue download. Outputs direct
Newtonian-gauge long P_(delta_cb,r) / (i mu), but DOES NOT multiply by
E11 Gamma because E11 rank uses the original density-derived R16 sigma.
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

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
PROTO=ROOT/"source_data/eboss_dr16_a03_e12_direct_CLASS_vTk_long_mode_cross_protocol_2026-09-27.json"
PROTO_BLOB="4805ce5f70d0908ccbf26ccc0a3ec5812f4ce757"
ARCH=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27"
MANIFEST_SHA="83748b452f26b31da34a51c2c5de509d6b7ae3af"
BASE=ROOT/"code/wake_two_tracer_fisher.py"
BASE_SHA="473e6d9b5941d8f67d1ecc1f5f170124d9ef8ea2"
CLASS_SHA="e85808324f51fc694d12e3ed7439552a3c3f9540"
ORIG_E8="bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
ORIG_E9="450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890"
ORIG_E11="cade0959c112ca57e9cd5fd4433a49cbec3be0b74eb1f0a1d6eab9350a756ebd"
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
      "minus":"Fminus_CLASS_normalized"}
LONGK=np.array([.001,.002,.003,.005],dtype=float)
KVEL=np.geomspace(1e-4,.15,180,dtype=float)
ZC=.95
HSTEP=(.004,.002)
C=299792.458
K_PIVOT=.05
CUTOFF=.1
R16=16.
OUTPUT=ROOT/"eboss_workspace/a03_physics_source/e12_direct_vTk_long_cross"

def check(x,msg):
    if not x:raise ValueError(msg)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def deterministic_json(o):
    return (json.dumps(o,indent=2,allow_nan=False)+"\n").encode()

def save_new(path,data):
    raw=deterministic_json(data)
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        check(path.read_bytes()==raw,"Original E12 output already exists with DIFFERENT bytes")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e12_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E12_OUTPUT",path,"SHA256",sha(raw),flush=True)

def original_inputs():
    check(blob(PROTO.read_bytes())==PROTO_BLOB and
          blob((ARCH/"archive_manifest.json").read_bytes())==MANIFEST_SHA and
          blob(BASE.read_bytes())==BASE_SHA,
          "Original E12 prereg/old E8-E11 archive/original velocity source changed")
    p=json.loads(PROTO.read_bytes())
    archived=json.loads((ARCH/"archive_manifest.json").read_bytes())
    original={x["relative_path"]:x["sha256"] for x in archived["files"]}
    needed={"E8/frozen_matched_distributions_4000q.csv":ORIG_E8,
            "E9/quasistatic_alpha_dln_a_source_only_plus_CLASS.json":ORIG_E9,
            "E11/a03_e11_gaussian_squeezed_source_only_report.json":ORIG_E11}
    for file,digest in needed.items():
        check(original.get(file)==digest and sha((ARCH/file).read_bytes())==digest,
              "Frozen parent SHA mismatch "+file)
    check(np.array_equal(np.array(p["locked_E12_geometry"]["K_LONG_COMOVING_hMpc"],float),LONGK)
          and p["locked_E12_geometry"]["redshift_center"]==ZC
          and p["locked_E12_geometry"]["redshift_density_derivative_relative_steps"]==list(HSTEP)
          and p["locked_E12_geometry"]["original_R16_comoving_Mpc_over_h"]==R16
          and p["locked_E12_geometry"]["integrated_K_upper_hMpc"]==CUTOFF
          and p["scientific_STOP"]["observed_galaxy_random_or_odd_read"] is False
          and p["scientific_STOP"]["do_not_multiply_E11_Gamma_by_new_direct_P_Dr_until_same_rank_sigma_compatible_or_Gamma_recomputed"] is True,
          "Preregistered E12 science and observed seal changed")
    c=np.genfromtxt(ARCH/"E8/frozen_matched_distributions_4000q.csv",
                    delimiter=",",names=True)
    check(c.shape==(4000,) and
          set(COLS.values()).issubset(set(c.dtype.names or ())) and
          np.all(np.diff(c["q_dimensionless"])>0),
          "E8 source CSV frozen 4000q invalid")
    e9=json.loads((ARCH/"E9/quasistatic_alpha_dln_a_source_only_plus_CLASS.json").read_bytes())
    check(e9["original_Eq20_exact_FD_vs_printed_mu_approx_normalization_unresolved"] is True
          and e9["observed_odd_data_vector_read"] is False,
          "Original E9 source-only/preregistered physical stop changed")
    return p,c,e9,original

def tophat(x):
    x=np.asarray(x)
    w=np.ones_like(x)
    m=np.abs(x)>1e-5
    w[m]=3*(np.sin(x[m])-x[m]*np.cos(x[m]))/x[m]**3
    return w

def cumulative(y,k):
    out=np.zeros_like(y)
    out[1:]=np.cumsum(.5*(y[1:]+y[:-1])*np.diff(np.log(k)))
    return out

def var_to_cutoff(vel):
    w=vel*vel*tophat(KVEL*R16)**2
    return float(np.interp(CUTOFF,KVEL,cumulative(w,KVEL)))

def primordial_delta2(k_mpc,As,ns):
    return As*(k_mpc/K_PIVOT)**(ns-1.)

def cross_i_mu_real(delta_cb,theta_rel,k_mpc,sigma_los,As,ns):
    """<delta_cb (v_LOS/sigma)^*> = i mu times this Mpc³ coefficient."""
    check(k_mpc>0 and sigma_los>0,"Nonphysical k/sigma")
    P=2*math.pi**2/k_mpc**3*primordial_delta2(k_mpc,As,ns)
    return float(P*delta_cb*C*theta_rel/(k_mpc*sigma_los))

def algebraic_test():
    k=.002*.6736
    z=cross_i_mu_real(0.,8.,k,90.,2.1e-9,.9649)
    t=cross_i_mu_real(25.,0.,k,90.,2.1e-9,.9649)
    p=cross_i_mu_real(25.,8.,k,90.,2.1e-9,.9649)
    check(z==0. and t==0. and p>0
          and math.isclose(cross_i_mu_real(-25.,8.,k,90.,2.1e-9,.9649),-p,rel_tol=1e-15)
          and np.allclose(cumulative(np.ones_like(KVEL),KVEL)[-1],
                          np.log(KVEL[-1]/KVEL[0]),rtol=5e-14),
          "E12 Fourier sign/cross null or R16 log-integral control failed")
    print("E12_LOCKED_ORIGINAL_INPUTS_DIRECT_THETA_VS_DENSITY_PROXY_FOURIER_PARITY_SELF_TEST_OK",flush=True)

def state_calc(state,out):
    check(state in STATES,"Nonfrozen source state")
    p,c,e9,files=original_inputs()
    algebraic_test()
    import wake_two_tracer_fisher as base
    import class_response_optimize as cro
    from classy import Class
    check(cro.CLASS_COMMIT==CLASS_SHA and
          np.array_equal(base.KVEL,KVEL) and
          base.R_FILTER==R16 and
          base.NS==.9649 and
          float(cro.H0)==67.36 and
          math.isclose(base.C_KMS,C),
          "Original CLASS/model source values drifted")
    a=e9["CLASS_state_JSON_sha256"][state]
    source=(ARCH/"E9"/("class_R16_wind_"+state+".json"))
    check(sha(source.read_bytes())==a,"Original E9 per-state source SHA changed")
    old=json.loads(source.read_bytes())
    check(old["state"]==state and old["CLASS_commit_expected"]==CLASS_SHA
          and old["one_sigma_positive_los_rank_not_measured_halo_wind"] is True,
          "Frozen E9 source run status changed")
    workdir=out/("state_"+state)
    workdir.mkdir(parents=True,exist_ok=True)
    q=np.asarray(c["q_dimensionless"],float)
    f=np.asarray(c[COLS[state]],float)
    check(np.min(f)>0 and np.all(np.isfinite(f)),
          "Nonpositive original frozen source")
    psd=workdir/("original_F_"+state+".dat")
    psdraw="\n".join(f"{qi:.14e} {fi:.14e}" for qi,fi in zip(q,f))+"\n"
    if psd.exists():check(psd.read_text()==psdraw,"Original generated PSD was replaced")
    else:psd.write_text(psdraw)
    orig_psd_hash=old["frozen_state_PSD_SHA256"]
    psd_match=(sha(psd.read_bytes())==orig_psd_hash)
    # E12 originates from BYTE-EXACT archived original 4000q CSV. If a different
    # Python/NumPy release rounds the new two-column PSD differently, report
    # it; do NOT patch E9 SHA or relax its archival parent.
    base.ZBINS=np.array([.9,1.],dtype=float)
    params=base.class_params(psd,.06)
    params["output"]="mPk,dTk,vTk"
    oldp=old["CLASS_params_from_original_wake_two_tracer_fisher"]
    for key in oldp:
        if key!="ncdm_psd_filenames":
            check(params.get(key)==oldp[key],
                  "CLASS parameter changed from completed E9: "+key)
    check(params["gauge"]=="newtonian" and params["z_max_pk"]>=1.05,
          "Original CLASS gauge/z range changed")
    cosmo=Class()
    cosmo.set(params);cosmo.compute()
    try:
        h=float(cosmo.h())
        check(math.isclose(h,.6736,rel_tol=0,abs_tol=2e-14),
              "CLASS original h drift")
        zm={(step,sgn):ZC+sgn*step*(1.+ZC) for step in HSTEP for sgn in (-1,1)}
        zs={ZC,*zm.values()}
        tcache={}
        kin=None
        for z in sorted(zs):
            t=cosmo.get_transfer(z=float(z),output_format="class")
            names=set(t.keys())
            if kin is None:
                check({"k (h/Mpc)","t_cdm","t_ncdm[0]","d_cdm","d_ncdm[0]","d_b"}.issubset(names),
                      "Direct CLASS Newtonian dTk/vTk required columns missing, STOP")
                kin=np.asarray(t["k (h/Mpc)"],float)
                check(kin.ndim==1 and np.all(np.diff(kin)>0)
                      and kin[0]<=min(KVEL[0],LONGK[0])
                      and kin[-1]>=KVEL[-1],
                      "CLASS k grid does NOT cover frozen KVEL and LONGK: extrapolation forbidden")
            check(np.array_equal(kin,np.asarray(t["k (h/Mpc)"],float)),
                  "Different CLASS transfer k grid across redshift")
            tcache[z]=t
        def interp(z,key,kh):
            return np.interp(kh,kin,np.asarray(tcache[z][key],float))
        Hc=float(cosmo.Hubble(ZC))
        k_mpc=KVEL*h
        tilt=np.sqrt(primordial_delta2(k_mpc,cro.A_S,base.NS))
        t0=tcache[ZC]
        th=np.asarray(interp(ZC,"t_ncdm[0]",KVEL)-interp(ZC,"t_cdm",KVEL))
        veld=-C*th/k_mpc*tilt
        sigma_direct=math.sqrt(var_to_cutoff(veld)/3.)
        check(sigma_direct>0 and np.isfinite(sigma_direct),
              "E12 direct theta-R16 sigma nonfinite/zero")
        proxy_by_step={}
        for hv in HSTEP:
            zp=zm[(hv,1)];zn=zm[(hv,-1)]
            dnu=(interp(zp,"d_ncdm[0]",KVEL)-interp(zn,"d_ncdm[0]",KVEL))/(zp-zn)
            dc=(interp(zp,"d_cdm",KVEL)-interp(zn,"d_cdm",KVEL))/(zp-zn)
            theta_proxy=Hc*(dnu-dc)
            velproxy=-C*theta_proxy/k_mpc*tilt
            sigma_proxy=math.sqrt(var_to_cutoff(velproxy)/3.)
            old_sigma=old["z_records"][[q["z"] for q in old["z_records"]].index(ZC)]["wind_by_original_relative_z_fd_step"][str(hv)]["v_parallel_plus_one_LOS_sigma_kms"]
            # signed dense operator v at each sample retains every frozen K
            per_long=[]
            LKM=LONGK*h
            dcl=interp(zp,"d_cdm",LONGK);dclm=interp(zn,"d_cdm",LONGK)
            dnl=interp(zp,"d_ncdm[0]",LONGK);dnlm=interp(zn,"d_ncdm[0]",LONGK)
            thetaL=interp(ZC,"t_ncdm[0]",LONGK)-interp(ZC,"t_cdm",LONGK)
            proxyL=Hc*((dnl-dnlm)-(dcl-dclm))/(zp-zn)
            db=interp(ZC,"d_b",LONGK)
            dc0=interp(ZC,"d_cdm",LONGK)
            w_b=float(cro.OMEGA_B)
            w_c=float(cro.OMEGA_CDM)
            cb=(w_b*db+w_c*dc0)/(w_b+w_c)
            metric_key="phi_prime" if "phi_prime" in names else None
            pmetric=interp(ZC,metric_key,LONGK) if metric_key else None
            theta_cdm=interp(ZC,"t_cdm",LONGK)
            CDMder=Hc*(dcl-dclm)/(zp-zn)
            cdm_res=(theta_cdm-CDMder-3*pmetric) if metric_key else None
            for idx,kval in enumerate(LONGK):
                sig=primordial_delta2(LKM[idx],cro.A_S,base.NS)
                imunit=cross_i_mu_real(cb[idx],thetaL[idx],LKM[idx],
                                       sigma_direct,cro.A_S,base.NS)
                res_norm=(float(cdm_res[idx]/max(abs(theta_cdm[idx]),abs(CDMder[idx]),abs(3*pmetric[idx]),1e-16))
                          if metric_key else None)
                per_long.append({
                    "K_comoving_h_per_Mpc":float(kval),"K_Mpc_inv":float(LKM[idx]),
                    "primordial_curvature_Delta2":float(sig),
                    "delta_cb_Newtonian_per_primordial_curvature":float(cb[idx]),
                    "theta_relative_direct_vTk_per_primordial_curvature_Mpc_inv":float(thetaL[idx]),
                    "theta_relative_proxy_H_dDeltaRel_dz_Mpc_inv":float(proxyL[idx]),
                    "theta_direct_over_proxy_signed":(float(thetaL[idx]/proxyL[idx])
                            if abs(proxyL[idx])>1e-15 else None),
                    "v_LOS_direct_transfer_per_primordial_curvature_at_mu_pos1_kms":
                        float(-C*thetaL[idx]/LKM[idx]),
                    "P_delta_cb_vrel_LOS_over_i_mu_Mpc3_kms":
                        float(imunit*sigma_direct),
                    "P_delta_cb_rdirect_LOS_over_i_mu_Mpc3":imunit,
                    "CDM_theta_vs_HdDeltaPlus3phiPrime_relative_closure":res_norm,
                    "CDM_metric_phi_prime_present":bool(metric_key)})
            proxy_by_step[str(hv)]={
                "center_z":ZC,
                "relative_density_derivative_step":hv,
                "original_E9_vLOS_proxy_kms":float(old_sigma),
                "recomputed_density_proxy_vLOS_kms":sigma_proxy,
                "replay_relative_difference":float(abs(sigma_proxy-old_sigma)/old_sigma),
                "vs_direct_sigma_ratio":float(sigma_direct/sigma_proxy),
                "v_direct_vs_proxy_signed_cosine_R16_weighted":float(np.sum(
                    veld*velproxy* (tophat(KVEL*R16)**2)*
                    np.gradient(np.log(KVEL)))/math.sqrt(
                    np.sum(veld*veld*(tophat(KVEL*R16)**2)*np.gradient(np.log(KVEL)))*
                    np.sum(velproxy*velproxy*(tophat(KVEL*R16)**2)*np.gradient(np.log(KVEL))))),
                "long_K":per_long}
        def gap(v,w):return float(abs(v-w)/max(abs(v),abs(w),1e-12))
        warn=[]
        for hv in HSTEP:
            v=proxy_by_step[str(hv)]
            if v["replay_relative_difference"]>.005:
                warn.append("Recomputed original E9 density-derived sigma does not reproduce frozen report for "+str(hv))
        if gap(proxy_by_step["0.004"]["recomputed_density_proxy_vLOS_kms"],
               proxy_by_step["0.002"]["recomputed_density_proxy_vLOS_kms"])>.01:
            warn.append("E9 density-proxy relative-z derivative finite-step instability > 1pct")
        if metric_key:
            vals=[abs(x["CDM_theta_vs_HdDeltaPlus3phiPrime_relative_closure"])
                  for step in proxy_by_step.values() for x in step["long_K"]]
            if max(vals)>.01: warn.append("CDM metric-continuity transfer closure > 1pct")
        result={
            "status":"E12_FROZEN_DIRECT_CLASS_VTK_LONG_SOURCE_AUDIT_COMPLETE_NOT_EBOSS_BISPECTRUM",
            "state":state,
            "CLASS_commit":CLASS_SHA,
            "E12_prospective_protocol_git_blob":PROTO_BLOB,
            "original_E8_4000q_CSV_SHA256":ORIG_E8,
            "original_E9_three_CLASS_JSON_SHA256":ORIG_E9,
            "original_E9_per_state_JSON_SHA256":a,
            "original_E11_Gamma_JSON_SHA256":ORIG_E11,
            "original_E9_PSD_SHA256":orig_psd_hash,
            "reconstructed_E12_PSD_SHA256":sha(psd.read_bytes()),
            "reconstructed_PSD_byte_exact_original_E9":psd_match,
            "parameters_same_as_E9_except_psd_output_path":True,
            "CLASS_output_direct_theta":{"gauge":"newtonian","output":"mPk,dTk,vTk",
                "transfer_keys":sorted(tcache[ZC].keys()),
                "raw_velocity_columns":["t_ncdm[0]","t_cdm"],
                "metric_phi_prime_present":bool(metric_key),
                "H_CLASS_Mpc_inverse":Hc,"h":h},
            "original_KVEL_comoving_h_per_Mpc":KVEL.tolist(),
            "center_theta_relative_direct_vTk_KVEL_Mpc_inv":th.tolist(),
            "center_direct_minus_proxy_density_KVEL_available_via_perK_steps":False,
            "direct_R16_sigma_LOS_kms":sigma_direct,
            "original_E9_R16_sigma_LOS_density_derived_model_not_identical_to_direct":True,
            "density_proxy_two_fixed_steps_and_four_long_K":proxy_by_step,
            "technical_QA_warnings_not_physical_rejection":warn,
            "interpretation":"Direct theta_ncdm-theta_cdm is CLASS velocity-divergence transfer. H d(dnu-dcdm)/dz is original E9 density proxy; ncdm pressure and metric terms may cause a physical difference. Direct normalized r requires direct_R16 sigma, E11 Gamma is based on E9 proxy sigma.",
            "P_delta_cb_rdirect_sign_and_units":"P_<delta_cb,v_LOS/sigma> = +i (Khat dot LOS) * real_table [Mpc^3] for v_LOS=-i (Khat dot LOS)c*theta_rel/K R, and Newtonian cb density transfer. Real table is not galaxy bispectrum and not gauge-independent density.",
            "E11_Gamma_multiplied_into_real_P_Dr":False,
            "full_Einstein_Vlasov_retarded_halo_or_real_eBOSS_triple_window":False,
            "observed_galaxy_random_or_odd_read":False,
            "new_catalogue_mock_download_or_science_cut_seed":False}
        save_new(out/("e12_direct_vTk_long_cross_"+state+".json"),result)
        print("E12_CLASS_DIRECT_VTK_STATE",state,
              "SIGMA_LOS_DIRECT_KMS",sigma_direct,
              "SIGMA_LOS_OLD_PROXY_KMS",proxy_by_step["0.004"]["original_E9_vLOS_proxy_kms"],
              "PROXY_OVER_DIRECT",proxy_by_step["0.004"]["recomputed_density_proxy_vLOS_kms"]/sigma_direct,
              "WARNINGS",len(warn),flush=True)
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()
    return 0

def aggregate(out):
    p,c,e9,files=original_inputs()
    algebraic_test()
    states={}
    for state in STATES:
        path=out/("e12_direct_vTk_long_cross_"+state+".json")
        raw=path.read_bytes(); d=json.loads(raw)
        check(d["state"]==state and d["CLASS_commit"]==CLASS_SHA and
              d["E12_prospective_protocol_git_blob"]==PROTO_BLOB and
              d["observed_galaxy_random_or_odd_read"] is False,
              "E12 source-only parent missing or not original")
        states[state]={"sha256":sha(raw),"record":d}
    comp={}
    for state in STATES:
        d=states[state]["record"]
        proxy=d["density_proxy_two_fixed_steps_and_four_long_K"]["0.004"]["recomputed_density_proxy_vLOS_kms"]
        direct=d["direct_R16_sigma_LOS_kms"]
        comp[state]={"direct_sigma_LOS_kms":direct,
                     "frozen_E9_proxy_sigma_LOS_kms":d["density_proxy_two_fixed_steps_and_four_long_K"]["0.004"]["original_E9_vLOS_proxy_kms"],
                     "proxy_to_direct":proxy/direct,
                     "fractional_direct_minus_proxy_over_proxy":(direct-proxy)/proxy,
                     "n_tech_warnings":len(d["technical_QA_warnings_not_physical_rejection"])}
    joint={"status":"E12_DIRECT_CLASS_VTK_VS_DENSITY_PROXY_LONG_CROSS_ALL_THREE_F_STATES_COMPLETE_CONDITIONAL_SOURCE_ONLY",
        "original_E8_CSV_SHA256":ORIG_E8,"original_E9_full_SHA256":ORIG_E9,
        "original_E11_full_SHA256":ORIG_E11,
        "preregistered_E12_protocol_git_blob":PROTO_BLOB,
        "per_state_output_SHA256":{s:states[s]["sha256"] for s in STATES},
        "matched_original_model_direct_vs_E9_density_proxy":comp,
        "original_fixed_K_LONG_comoving_h_Mpc":LONGK.tolist(),
        "all_per_state_full_records_included":False,
        "long_signed_Newtonian_cb_density_to_LOS_relative_wind_cross_tables_available_each_state":True,
        "CDM_continuity_metric_if_available_required_before_interpretation":True,
        "strict_E11_direct_rank_compatibility_not_assumed":True,
        "actual_full_eBOSS_bispectrum_or_covariance_calculated":False,
        "observed_odd_data_sealed":True,"main_untouched":True}
    save_new(out/"e12_all_states_summary_source_only.json",joint)
    print("E12_ALL_F_STATES_DIRECT_VTK_VS_DENSITY_PROXY_COMPLETED",
          json.dumps(comp,sort_keys=True),flush=True)
    return 0

def main():
    a=argparse.ArgumentParser(description=__doc__)
    g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--state",choices=STATES)
    g.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUTPUT)
    args=a.parse_args()
    if args.self_test:
        original_inputs();algebraic_test();return 0
    if args.state:return state_calc(args.state,args.output_dir)
    return aggregate(args.output_dir)

if __name__=="__main__":
    raise SystemExit(main())
