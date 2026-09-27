#!/usr/bin/env python3
"""E12 original three-state CLASS CDM-density x R16 ν-CDM relative-wind long modes.

Source inputs ONLY: original byte-pinned E8 CSV, E9 complete CLASS JSON and E11
Gaussian-rank conditional Fourier Γ response. For original long K fixed before
new CLASS, compute dimensionless per-logK signed long delta-CDM/LOS relative
wind cross spectrum, then formal squeezed normalized Gamma*C. This does not
predict observed eBOSS bispectrum, 24D odd xi, or real halo velocity.
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
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"code"))
import build_eboss_dr16_a03_e8_frozen_resonant_occupation as e8

PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e12_squeezed_long_CLASS_transfer_prereg_2026-09-27.json"
PROTOCOL_BLOB="57e979b72baaded262e24a653a619ff6a4026286"
E9_SOURCE=ROOT/"scripts/build_eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase.py"
E9_SOURCE_BLOB="700a80d9ea9bc45454925b9a06fec2106203c506"
BASE_SOURCE=ROOT/"code/wake_two_tracer_fisher.py"
BASE_BLOB="473e6d9b5941d8f67d1ecc1f5f170124d9ef8ea2"
ORIG_E8="bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
ORIG_E9="450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890"
ORIG_E11="fffdd6fac67256631e5981182e1b8ba6f48eb05ad14ddbec48e6d8eaa98b45eb"
KLONG=np.array([.003,.005,.01],float)
KVEL=np.geomspace(1e-4,.15,180)
STATES=("FD","plus","minus")
Z=.95
R=16.
K_CUT=.10

def need(ok,msg):
    if not ok:raise ValueError(msg)

def canonical(obj):
    return (json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()

def new_only(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        need(path.read_bytes()==raw,"E12 report already exists with different source bytes")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e12_",delete=False) as f:
        temp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(temp,path)
    finally:temp.unlink(missing_ok=True)

def pin_inputs(csv,e9_report,e11_report):
    for path,blob in ((PROTOCOL,PROTOCOL_BLOB),(E9_SOURCE,E9_SOURCE_BLOB),
                      (BASE_SOURCE,BASE_BLOB)):
        need(e8.git_blob(path.read_bytes())==blob,
             "E12 preregistration or ORIGINAL velocity/CLASS pipeline Git blob changed: "+str(path))
    e8.source_gate()
    p=json.loads(PROTOCOL.read_bytes())
    f=p["fixed_source_k_and_z"]
    need(f["long_K_h_per_Mpc"]==KLONG.tolist() and
         f["short_k_h_per_Mpc"]==.05 and
         f["center_z"]==Z and
         f["redshift_derivative_relative_z_steps"]==[.004,.002] and
         f["original_sigma_R16_truncation_K_upper_h_per_Mpc"]==K_CUT and
         f["original_R_filter_comoving_Mpc_per_h"]==R and
         p["physical_limits"]["observed_galaxies_or_randoms_or_odd_read"] is False and
         p["physical_limits"]["no_new_catalogue_or_mock_download"] is True,
         "Precomputed E12 K/z/mu and observed data guard changed")
    need(e8.sha(csv.read_bytes())==ORIG_E8 and
         e8.sha(e9_report.read_bytes())==ORIG_E9 and
         e8.sha(e11_report.read_bytes())==ORIG_E11,
         "Original E8/E9/E11 report byte SHA mismatch")
    prior=json.loads(e9_report.read_bytes())
    last=json.loads(e11_report.read_bytes())
    need(prior["CLASS_commit_expected"]==e8.cro.CLASS_COMMIT and
         prior["observed_odd_data_vector_read"] is False and
         last["zero_intrinsic_unconditional_odd_mean_under_symmetric_model"] is True and
         last["A04_valid_independent_bispectrum_covariance_available"] is False,
         "Original E9/E11 user seal or physical estimand changed")
    src=np.genfromtxt(csv,delimiter=",",names=True)
    need(src.dtype.names==("q_dimensionless","F0_CLASS_normalized",
         "Fplus_CLASS_normalized","Fminus_CLASS_normalized") and src.shape==(4000,),
         "Original archived F0/F± E8 4000-q CSV source changed")
    return p,prior,last,src

def e9_reference(e9,state,step):
    matches=[j for j in e9["all_own_and_common_FD_wind_cases"] if
             j["z"]==Z and j["wind_velocity_derivative_relative_z_step"]==step]
    need(len(matches)==1,"Original E9 z=.95 matching fixed velocity derivative step missing")
    return matches[0]["state"][state]

def class_state(out,state,csv,e9report,e11report):
    need(state in STATES,"Only predeclared original FD, plus, minus allowed")
    p,prev,last,src=pin_inputs(csv,e9report,e11report)
    from classy import Class
    import wake_two_tracer_fisher as base
    need(np.array_equal(base.KVEL,KVEL) and base.R_FILTER==R and
         base.NS==.9649 and base.KPIV_MPC==.05 and
         math.isclose(base.C_KMS,e8.C_KMS,rel_tol=1e-15),
         "Original three-state relative velocity implementation/units changed")
    out.mkdir(parents=True,exist_ok=True)
    q=src["q_dimensionless"]
    vals=np.asarray(src[{"FD":"F0_CLASS_normalized",
         "plus":"Fplus_CLASS_normalized","minus":"Fminus_CLASS_normalized"}[state]],float)
    need(np.isfinite(q).all() and np.isfinite(vals).all() and np.min(vals)>0,
         "Original matched E8 source positivity/numerical invalid")
    psd=out/("e12_archived_E8_state_"+state+".dat")
    with tempfile.NamedTemporaryFile(dir=out,prefix=".e12_psd_",mode="w",
           encoding="utf-8",delete=False) as fh:
        t=Path(fh.name)
        np.savetxt(fh,np.column_stack([q,vals]),fmt="%.14e")
        fh.flush();os.fsync(fh.fileno())
    try:
        if psd.exists():need(psd.read_bytes()==t.read_bytes(),
                             "Regenerated original E8 PSD state file changed")
        else:os.link(t,psd)
    finally:t.unlink(missing_ok=True)
    base.ZBINS=np.array([.90,1.],float)
    params=base.class_params(psd,e8.MASS_EV)
    params["output"]="mPk,dTk,vTk"
    need(params["gauge"]=="newtonian" and
         params["z_max_pk"]>=1.05 and
         params["P_k_max_h/Mpc"]>=KVEL[-1] and
         params["use_ncdm_psd_files"]==1 and
         params["m_ncdm"]==e8.MASS_EV,
         "Original CLASS backend/source/settings mismatch")
    c=Class();c.set(params);c.compute()
    try:
        h=float(c.h())
        need(math.isclose(h,e8.cro.H0/100.,rel_tol=5e-14),
             "E12 CLASS H0 not original E8/E9 universe")
        center=c.get_transfer(z=Z,output_format="class")
        kclass=np.asarray(center["k (h/Mpc)"],float)
        need(kclass.ndim==1 and np.all(np.diff(kclass)>0)
             and kclass[0]<=KVEL[0] and kclass[-1]>=KVEL[-1],
             "Original CLASS density transfer does not cover original sigma and long K without extrapolation")
        def vec(t,key,kk):
            kin=np.asarray(t["k (h/Mpc)"],float)
            need(kin[0]<=np.min(kk) and kin[-1]>=np.max(kk),
                 "CLASS transfer edge outside original, fail-closed (no extrapolation)")
            x=np.asarray(t[key],float)
            need(np.isfinite(x).all(),"CLASS source transfer nonfinite: "+key)
            return np.interp(kk,kin,x)
        Hc=float(c.Hubble(Z))
        amp=np.sqrt(e8.cro.A_S*(KLONG*h/base.KPIV_MPC)**(base.NS-1.))
        delta_c=vec(center,"d_cdm",KLONG)
        D_c=delta_c*amp
        stage={}
        for ds in (.004,.002):
            dz=ds*(1.+Z)
            tm=c.get_transfer(z=Z-dz,output_format="class")
            tp=c.get_transfer(z=Z+dz,output_format="class")
            def derivative(kh):
                dm=vec(tm,"d_ncdm[0]",kh)-vec(tm,"d_cdm",kh)
                dp=vec(tp,"d_ncdm[0]",kh)-vec(tp,"d_cdm",kh)
                return (dp-dm)/(2.*dz)
            ddV=derivative(KLONG)
            ddVgrid=derivative(KVEL)
            signedV=-e8.C_KMS*Hc/(KLONG*h)*ddV*amp
            grid_amp=np.sqrt(e8.cro.A_S*(KVEL*h/base.KPIV_MPC)**(base.NS-1.))
            gridV=-e8.C_KMS*Hc/(KVEL*h)*ddVgrid*grid_amp
            power=gridV**2*base.tophat(KVEL*R)**2
            intvar=base.cumulative_trapz_logk(power,KVEL)
            sigmaLOS=float(np.sqrt(float(np.interp(K_CUT,KVEL,intvar))/3.))
            ref=e9_reference(prev,state,ds)["own_CLASS_v_R16_1sigma_los_kms"]
            replay=abs(sigmaLOS-ref)/max(abs(sigmaLOS),abs(ref),1e-12)
            need(replay < p["preregistered_technical_QA"][
                 "E9_sigma_replay_vs_original_max_relative"],
                 "Original R16 CLASS sigma_R16 failed its separately preregistered E9 replay")
            filt=base.tophat(KLONG*R)
            cross=signedV*D_c*filt/sigmaLOS
            need(np.isfinite(cross).all() and np.isfinite(signedV).all(),
                 "Nonfinite longitudinal v-CDM dimensionless cross power")
            stage[str(ds)]={
              "original_E9_LOS_sigma_kms":ref,
              "recomputed_LOS_sigma_kms":sigmaLOS,
              "sigma_relative_replay_gap":replay,
              "signed_linear_Vrel_long_mode_kms_per_sqrt_DeltaR2":signedV.tolist(),
              "delta_c_transfer":delta_c.tolist(),
              "delta_c_amplitude_times_sqrt_DeltaR2":D_c.tolist(),
              "original_R16_tophat_window_W_KR":filt.tolist(),
              "signed_real_Delta2_delta_c_relative_LOS_rank_per_i_muLong":cross.tolist(),
              "velocity_h_phys_Hubble_CLASS_inverse_Mpc":Hc
            }
        a=np.asarray(stage["0.004"]["signed_real_Delta2_delta_c_relative_LOS_rank_per_i_muLong"])
        b=np.asarray(stage["0.002"]["signed_real_Delta2_delta_c_relative_LOS_rank_per_i_muLong"])
        gap=np.abs(a-b)/np.maximum.reduce([np.abs(a),np.abs(b),np.full_like(a,1e-15)])
        warnings=[int(i) for i in np.flatnonzero(gap>p["preregistered_technical_QA"][
                  "original_long_v_density_delta_stencil_0p004_vs_0p002_rel_warn"])]
        output={
          "status":"E12_PINNED_LONG_CDM_VREL_CLASS_STATE_COMPLETED_SOURCE_ONLY_NOT_EBOSS_TRACER",
          "state":state,"original_E8_CSV_sha256":ORIG_E8,
          "original_E9_JSON_sha256":ORIG_E9,
          "original_E11_Gamma_JSON_sha256":ORIG_E11,
          "original_CLASS_commit_expected":e8.cro.CLASS_COMMIT,
          "CLASS_params":params,
          "state_archived_E8_sourced_PSD_sha256":e8.sha(psd.read_bytes()),
          "z":Z,"short_k_com_h_per_Mpc":.05,
          "long_K_com_h_per_Mpc":KLONG.tolist(),
          "Delta_R2_as_input_A_s_and_original_ns":(amp*amp).tolist(),
          "long_mode_cross_spectrum_definition":"For delta(x)=int exp(iKx)delta(K), v_rel(K)=-i*Khat*V_signed(K)*zeta_unit(K); LOS Gaussian rank r=vrel_LOS,R16/sigma_LOS. Cross <delta_c(K)r(-K)> dimensionless Delta2 = +i*(Khat dot LOS)*[D_c*V_signed*W(KR16)/sigma_LOS]. Table records signed real coefficient in brackets. Not a measured/normalized galaxy bispectrum.",
          "recomputed_CLASS_source_each_ncdm_state":True,
          "stencil_results":stage,
          "stencil_relative_gaps":gap.tolist(),
          "technical_QA_warn_long_K_indices":warnings,
          "no_actual_LRG_ELG_bias_or_short_Pcb":True,
          "observed_odd_or_survey_data_read":False
        }
        raw=canonical(output);new_only(out/("e12_long_density_velocity_"+state+".json"),raw)
        print("E12_CLASS_LONG_CDM_WIND_STATE_COMPLETE",state,
              "SHA256",hashlib.sha256(raw).hexdigest(),
              "GAPS",gap.tolist(),"WARN",warnings,flush=True)
    finally:
        c.struct_cleanup();c.empty()

def aggregate(out,e8csv,e9report,e11report):
    proto,prev,gamma,src=pin_inputs(e8csv,e9report,e11report)
    table={}
    state_sha={}
    for state in STATES:
        path=out/("e12_long_density_velocity_"+state+".json")
        raw=path.read_bytes();j=json.loads(raw)
        need(j["state"]==state and
             j["original_CLASS_commit_expected"]==e8.cro.CLASS_COMMIT and
             j["original_E8_CSV_sha256"]==ORIG_E8 and
             j["original_E11_Gamma_JSON_sha256"]==ORIG_E11 and
             j["observed_odd_or_survey_data_read"] is False and
             j["technical_QA_warn_long_K_indices"]==[] and
             j["long_K_com_h_per_Mpc"]==KLONG.tolist(),
             "One of the three original CLASS long-transfer states did not pass fixed source or QA gate")
        table[state]=j
        state_sha[state]=hashlib.sha256(raw).hexdigest()
    original_g=gamma["all_angular_and_rank_quadrature_results"]["mu_512_GH_96"]
    output_rows=[]
    for ix,K in enumerate(KLONG):
        cell={"K_LONG_com_h_per_Mpc":float(K),
              "squeezed_K_over_short_k":float(K/.05),
              "per_unit_delta_b_Pcb_short_dimensionless_long_delta_c_response":{}}
        for state in STATES:
            c=table[state]["stencil_results"]["0.004"][
               "signed_real_Delta2_delta_c_relative_LOS_rank_per_i_muLong"][ix]
            ga={ell:original_g[state][ell][
                "Gaussian_squeezed_response_E_rank_times_T"] for ell in ("1","3")}
            cell["per_unit_delta_b_Pcb_short_dimensionless_long_delta_c_response"][state]={
              "signed_real_cross_Delta2_delta_c_r_per_i_muLong":c,
              "Gamma_ell":ga,
              "signed_real_Bhat_ell_per_i_muLong":{
                  ell:float(ga[ell]*c) for ell in ("1","3")}}
        aa=cell["per_unit_delta_b_Pcb_short_dimensionless_long_delta_c_response"]
        cell["plus_minus_difference"]={
          ell:aa["plus"]["signed_real_Bhat_ell_per_i_muLong"][ell]-
              aa["minus"]["signed_real_Bhat_ell_per_i_muLong"][ell]
              for ell in ("1","3")}
        cell["common_FD_long_transfer_control"]={
          ell:(aa["plus"]["Gamma_ell"][ell]-aa["minus"]["Gamma_ell"][ell])*
              aa["FD"]["signed_real_cross_Delta2_delta_c_r_per_i_muLong"]
              for ell in ("1","3")}
        cell["own_LONG_CLASS_transfer_correction"]={
          ell:cell["plus_minus_difference"][ell]-
              cell["common_FD_long_transfer_control"][ell]
              for ell in ("1","3")}
        for ell in ("1","3"):
            need(math.isclose(
                cell["plus_minus_difference"][ell],
                cell["common_FD_long_transfer_control"][ell]+
                cell["own_LONG_CLASS_transfer_correction"][ell],
                abs_tol=1e-16,rel_tol=1e-12),
                "Frozen E11 gamma/own CLASS long transfer algebraic difference failed")
        output_rows.append(cell)
    output={
      "date":"2026-09-27",
      "status":"E12_GENUINELY_SQUEEZED_ORIGINAL_CLASS_LONG_DENSITY_RELATIVE_WIND_SOURCE_KERNEL_COMPLETED_EBOSS_BISPECTRUM_INFERENCE_STOP",
      "preregistered_E12_protocol_git_blob":PROTOCOL_BLOB,
      "original_E8_E9_E11_source_SHA256":{"E8":ORIG_E8,"E9":ORIG_E9,"E11":ORIG_E11},
      "new_per_state_CLASS_long_density_wind_report_sha256":state_sha,
      "CLASS_commit":e8.cro.CLASS_COMMIT,
      "physical_fourier_convention":"delta(x)=int exp(iKx) delta(K); r(K)=-i muLong V_signed W/sigma. Delta2_{delta_c r}(K)=+i muLong C(K); leading squeezed Bhat(K)=C(K) Gamma_short. Values are signed real coefficients BEFORE i muLong, 2pi²/K³, delta_b, Pcb_short, long galaxy bias/window. No spurious nonzero 2point mean.",
      "original_same_common_r_Gaussian_ensemble_source_only":True,
      "long_K_predeclared_and_short_k":output_rows,
      "original_E10_unconditional_intrinsic_mean_zero":True,
      "actual_eBOSS_tracer_bias_and_3point_pair_triplet_window_available":False,
      "full_squeezed_geometry_halo_nonlocal_response_available":False,
      "original_v1_SGC_reverse_and_random_48k_stop_open":True,
      "three_point_valid_eBOSS_inferential_covariance_available":False,
      "new_catalogue_or_mock_download":False,
      "observed_galaxies_randoms_or_odd_read":False,
      "main_mutation":False
    }
    raw=canonical(output)
    new_only(out/"e12_CLASS_long_density_wind_squeezed_gamma_source_only.json",raw)
    print("EBOSS_A03_E12_ORIGINAL_CLASS_LONG_SQUEEZED_3POINT_KERNEL_COMPLETE",
          "SHA256",hashlib.sha256(raw).hexdigest(),
          "NO_EBOSS_XI NO_OBSERVED_ODD NO_BISPEC_DETECTION",flush=True)
    for item in output_rows:
        print("E12_LONG_K_DELTA_BHAT",item["K_LONG_com_h_per_Mpc"],
             json.dumps(item["plus_minus_difference"],sort_keys=True),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    m=p.add_mutually_exclusive_group(required=True)
    m.add_argument("--state",choices=STATES)
    m.add_argument("--aggregate",action="store_true")
    p.add_argument("--e8-csv",type=Path,required=True)
    p.add_argument("--e9-json",type=Path,required=True)
    p.add_argument("--e11-json",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    if a.aggregate:aggregate(a.outdir,a.e8_csv,a.e9_json,a.e11_json)
    else:class_state(a.outdir,a.state,a.e8_csv,a.e9_json,a.e11_json)
if __name__=="__main__":
    main()
