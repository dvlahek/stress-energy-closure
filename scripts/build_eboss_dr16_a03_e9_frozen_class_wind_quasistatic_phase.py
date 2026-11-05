#!/usr/bin/env python3
"""A03E9: three pinned CLASS states -> R16 velocity -> conditional Eq20 phase.

The frozen original F0/F+/F- sources are NOT refitted. CLASS density transfer
derivatives yield a model-defined, filtered one-sigma longitudinal wind under
a coherent-rank approximation, not the observed halo velocity realization.
This is a quasi-static exact-occupancy phase and its derivative, not a retarded
Einstein-Vlasov halo forecast, galaxy xi, survey window or inference.
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
import build_eboss_dr16_a03_e8b_conditional_static_halo_green_kernel as e8b

PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase_prereg_2026-09-27.json"
PROTOCOL_BLOB="a31249a6c429fb694aab525d539afd1392adf974"
E8C_ERRATUM=ROOT/"source_data/eboss_dr16_a03_e8c_exact_eq20_phase_and_bibliography_erratum_2026-09-27.json"
ERRATUM_BLOB="549dccc2b923800f06d0081741ea7c80f3181e4b"
BASE=ROOT/"code/wake_two_tracer_fisher.py"
BASE_BLOB="473e6d9b5941d8f67d1ecc1f5f170124d9ef8ea2"
STATES=("FD","plus","minus")
ZS=(.90,.94,.945,.95,.955,.96,1.)
KVEL=np.geomspace(1e-4,.15,180,dtype="f8")
SIGMA_CUTOFF=.10
R_FILTER=16.
OUT=ROOT/"eboss_workspace/a03_physics_source/e9_class_quasistatic_phase"

def require(cond,msg):
    if not cond:raise ValueError(msg)

def scope_gate():
    require(e8.git_blob(PROTOCOL.read_bytes())==PROTOCOL_BLOB,
            "Pre-CLASS E9 preregistration Git blob changed")
    require(e8.git_blob(BASE.read_bytes())==BASE_BLOB,
            "Original velocity-transfer implementation changed")
    require(e8.git_blob(E8C_ERRATUM.read_bytes())==ERRATUM_BLOB,
            "Post-original E8c exact FD/DOI erratum changed")
    e8.source_gate()
    pb=e8b.protocol()
    p=json.loads(PROTOCOL.read_bytes())
    q=p["fixed_numerical_grid"]
    require(p["frozen_physical_parameters"]["CLASS_commit"]==e8.cro.CLASS_COMMIT and
            p["frozen_physical_parameters"]["mass_eV"]==e8.MASS_EV and
            p["frozen_physical_parameters"]["frac_cap"]==e8.FRAC and
            p["frozen_physical_parameters"]["source_only_CSV_SHA256"]==
            "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0" and
            p["fixed_numerical_grid"]["compute_phase_z_samples"]==list(ZS) and
            p["fixed_numerical_grid"]["k_com_h_per_Mpc_phase"]==[.05] and
            p["fixed_numerical_grid"]["v_truncated_cumulative_kupper_h_per_Mpc"]==SIGMA_CUTOFF and
            p["fixed_numerical_grid"]["original_velocity_filter_R_comoving_Mpc_per_h"]==R_FILTER and
            q["velocity_transfer_derivative_relative_z_steps"]==[.004,.002] and
            q["phase_derivative_z_steps"]==[.01,.005] and
            p["limits"]["observed_odd_vector_read"] is False and
            p["limits"]["retarded_nonlinear_halo_wake_not_computed"] is True,
            "E9 science/numerical/observed STOP scope changed")
    return p,pb

def atomic_new(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,
                "Existing source-only output differs; do not overwrite")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e9_",delete=False) as fh:
        tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
    try:
        require(not path.exists(),"Concurrent E9 output appeared; refuse overwrite")
        os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)

def canonical(data):
    return (json.dumps(data,indent=2,allow_nan=False)+"\n").encode()

def state_name_to_f(src,state):
    return src[{"FD":"f0","plus":"fp","minus":"fm"}[state]]

def exact_occupation(src,state,q):
    """First exact Eq20 FD and frozen E8b interpolation-ratio convention."""
    require(math.isfinite(q) and src["q"][0]<=q<=src["q"][-1],
            "Resonance outside exact frozen q grid")
    fd=1./(1.+math.exp(q))
    if state=="FD":return fd
    f0=float(np.interp(q,src["q"],src["f0"]))
    fs=float(np.interp(q,src["q"],state_name_to_f(src,state)))
    require(f0>0 and fs>0 and fs/f0<=1.+e8.FRAC+1e-10 and
            fs/f0>=1.-e8.FRAC-1e-10,
            "Original E8b source ratio violates frozen 30% cap")
    return fd*fs/f0

def alpha_given_v(src,state,z,v_parallel_kms,kcom_h,physical):
    """Signed conditional static Eq20/33 instantaneous phase at one k_com."""
    require(state in STATES and z>=0 and kcom_h>0 and
            math.isfinite(v_parallel_kms) and math.isfinite(z),
            "Invalid phase source/state input")
    if v_parallel_kms==0.:return 0.
    r=physical["physical_reference"]
    mkg=e8.MASS_EV*r["rest_mass_conversion_eVc2_to_kg"]
    kcom=kcom_h*r["h"]/r["Mpc_m"]
    aa=1./(1.+z)
    vel=v_parallel_kms*1000.
    tnu0=e8.cro.T_NCDM*e8.cro.TCMB_K*e8.cro.KB_EV_K
    q=e8.MASS_EV*abs(v_parallel_kms)/(e8.C_KMS*tnu0*(1.+z))
    f=exact_occupation(src,state,q)
    return float(2.*r["Newton_G_m3_kg_s2"]*mkg**4*
                 aa*aa*vel*f/(r["hbar_Js"]**3*kcom*kcom))

def source_only_self_test():
    proto,pb=scope_gate()
    src=e8.kinetic_source()
    for z in (.9,.95,1.):
        for v in (0.,200.,-200.):
            aa=[alpha_given_v(src,s,z,v,.05,pb) for s in STATES]
            rev=[alpha_given_v(src,s,z,-v,.05,pb) for s in STATES]
            require(np.allclose(aa,-np.asarray(rev),atol=0,rtol=2e-13),
                    "Frozen exact FD phase loses velocity reversal")
            require(abs((aa[1]+aa[2])/2.-aa[0])<
                    max(abs(aa[0]),1e-300)*3e-13,
                    "Original exact Fplus/Fminus matched-source mean no longer FD")
    aa=alpha_given_v(src,"FD",.95,200.,.05,pb)
    require(math.isclose(aa,.0003087916839711642,rel_tol=5e-13),
            "E8c exact FD known comoving 0.05 phase changed")
    # Verify analytic derivative at fixed v; no CLASS transfer needed.
    z=.95;h=.005;v=200.;q=e8.MASS_EV*v/(e8.C_KMS*
         (e8.cro.T_NCDM*e8.cro.TCMB_K*e8.cro.KB_EV_K)*(1.+z))
    fd=1./(1.+math.exp(q))
    analyt=aa*(2.-q*(1.-fd))
    num=-(1.+z)*(alpha_given_v(src,"FD",z+h,v,.05,pb)-
         alpha_given_v(src,"FD",z-h,v,.05,pb))/(2*h)
    require(abs(analyt-num)/abs(analyt)<1.3e-5,
            "Fixed-wind dphase/dln a mathematical control failed")
    print("EBOSS_A03_E9_EXACT_SOURCE_ONLY_PHASE_FROZEN_FD_PARITY_AND_DLOGA_CONTROL_OK",
          "FD_ALPHA",aa,"FD_KINEMATIC_DLOGA",analyt,flush=True)
    return src,proto,pb

def run_state(outdir,state):
    require(state in STATES,"Only frozen FD/plus/minus states accepted")
    src,p,pb=source_only_self_test()
    from classy import Class
    import wake_two_tracer_fisher as base
    require(np.array_equal(KVEL,base.KVEL) and R_FILTER==base.R_FILTER and
            base.NS==.9649 and math.isclose(base.C_KMS,e8.C_KMS) and
            max(ZS)+.004*(1.+max(ZS))<1.05,
            "Original CLASS velocity grid/precision/velocity source mismatch")
    outdir.mkdir(parents=True,exist_ok=True)
    q=src["q"]
    distribution=state_name_to_f(src,state)
    psd=outdir/("frozen_class_state_"+state+".dat")
    with tempfile.NamedTemporaryFile(dir=outdir,prefix=".psd_",mode="w",
            delete=False,encoding="utf-8") as fh:
        temp=Path(fh.name)
        np.savetxt(fh,np.column_stack([q,distribution]),fmt="%.14e")
        fh.flush();os.fsync(fh.fileno())
    try:
        if psd.exists():require(psd.read_bytes()==temp.read_bytes(),
                                "Original source PSD worker bytes changed")
        else:os.link(temp,psd)
    finally:temp.unlink(missing_ok=True)
    base.ZBINS=np.array([.90,1.],dtype="f8")
    params=base.class_params(psd,e8.MASS_EV)
    params["output"]="mPk,dTk,vTk"
    require(params["gauge"]=="newtonian" and
            params["z_max_pk"]>=1.05 and
            params["P_k_max_h/Mpc"]>=KVEL[-1] and
            params["use_ncdm_psd_files"]==1,
            "Unsafe original pinned CLASS backend parameters")
    c=Class();c.set(params);c.compute()
    try:
        actual_h=float(c.h())
        require(math.isclose(actual_h,pb["physical_reference"]["h"],rel_tol=3e-6),
                "CLASS h changed from frozen E8b physical unit convention")
        cache={}
        def t_at(z):
            z=float(z)
            if z not in cache:
                t=c.get_transfer(z=z,output_format="class")
                kin=np.asarray(t["k (h/Mpc)"],dtype="f8")
                dnu=np.asarray(t["d_ncdm[0]"],dtype="f8")
                dcdm=np.asarray(t["d_cdm"],dtype="f8")
                require(kin.ndim==1 and kin.size>20 and
                        np.all(np.diff(kin)>0) and
                        np.isfinite(dnu).all() and np.isfinite(dcdm).all() and
                        kin[0]<=KVEL[0] and kin[-1]>=KVEL[-1],
                        "CLASS density transfer does not cover preregistered k grid; NO EXTRAPOLATION")
                cache[z]=np.interp(KVEL,kin,dnu-dcdm)
            return cache[z]
        z_records=[]
        print("E9_CLASS_START_STATE",state,"PINNED_CLASS",e8.cro.CLASS_COMMIT,
              "H",actual_h,flush=True)
        for z in ZS:
            stage={}
            Hc=float(c.Hubble(z))
            require(Hc>0 and math.isfinite(Hc),"Nonphysical CLASS H(z)")
            for delta in (.004,.002):
                dz=delta*(1.+z)
                ddelta=(t_at(z+dz)-t_at(z-dz))/(2.*dz)
                km=KVEL*actual_h
                tilt=np.sqrt(e8.cro.A_S*
                             (km/base.KPIV_MPC)**(base.NS-1.))
                vel=-Hc*ddelta/km*tilt*base.C_KMS
                require(np.isfinite(vel).all() and np.max(np.abs(vel))>0,
                        "Nonfinite or zero neutrino-CDM transfer-derived physical velocity proxy")
                weighted=vel**2*base.tophat(KVEL*R_FILTER)**2
                cum=base.cumulative_trapz_logk(weighted,KVEL)
                var=float(np.interp(SIGMA_CUTOFF,KVEL,cum))
                require(var>0 and math.isfinite(var),
                        "Nonpositive original R16-filtered transfer spectrum")
                stage[str(delta)]={"v_parallel_plus_one_LOS_sigma_kms":
                                    float(math.sqrt(var/3.)),
                                    "v_3D_R16_filtered_k_le_0p10_kms":
                                    float(math.sqrt(var)),
                                    "v_signed_transfer_at_kcom_0p05_kms":
                                    float(np.interp(.05,KVEL,vel))}
                if delta==.004 and z==.95:
                    ref=base.delta_v_kms(c,z,KVEL)
                    require(np.max(np.abs(np.abs(vel)-ref))<
                            5e-12*max(1.,float(np.max(np.abs(ref)))),
                            "Original CLASS density-derived relative velocity formula parity failed")
            rec={"z":float(z),"H_CLASS_inverse_Mpc":Hc,
                 "H_physical_SI_inverse_s":float(Hc*pb["physical_reference"]["c_m_per_s"]/
                                               pb["physical_reference"]["Mpc_m"]),
                 "wind_by_original_relative_z_fd_step":stage}
            z_records.append(rec)
            print("E9_CLASS_STATE_Z",state,format(z,".3f"),
                  "R16_1SIGMA_LOS_KMS",stage["0.004"]["v_parallel_plus_one_LOS_sigma_kms"],
                  "HALF_STEP",stage["0.002"]["v_parallel_plus_one_LOS_sigma_kms"],
                  flush=True)
        report={"status":"CLASS_PINNED_STATE_R16_WIND_ENGINEERING_COMPLETE_NOT_HALO_REALIZATION",
            "state":state,"CLASS_commit_expected":e8.cro.CLASS_COMMIT,
            "CLASS_params_from_original_wake_two_tracer_fisher":params,
            "frozen_source_input_CSV_sha256":
            p["frozen_physical_parameters"]["source_only_CSV_SHA256"],
            "frozen_state_PSD_SHA256":e8.sha(psd.read_bytes()),
            "transfer_expression":"-c*H_CLASS/k_Mpc * d(d_ncdm[0]-d_cdm)/dz * sqrt(A_s*(k_Mpc/0.05)^(n_s-1))",
            "original_KVEL_h_Mpc":KVEL.tolist(),
            "original_R_comoving_Mpc_per_h":R_FILTER,
            "truncation_k_upper_h_Mpc":SIGMA_CUTOFF,
            "one_sigma_positive_los_rank_not_measured_halo_wind":True,
            "z_records":z_records,
            "CLASS_native_get_transfer_calls_distinct_z":len(cache),
            "observed_odd_data_vector_read":False,
            "CLASS_state_galaxy_xi_not_computed":True}
        blob=canonical(report)
        atomic_new(outdir/("class_R16_wind_"+state+".json"),blob)
        print("EBOSS_A03_E9_FROZEN_CLASS_STATE_COMPLETE",
              state,"JSON_SHA256",e8.sha(blob),
              "TRANSFER_Z_CALLS",len(cache),flush=True)
    finally:
        c.struct_cleanup();c.empty()
    return 0

def phase_bundle(src,pb,wind,state,z,kcom=.05,windex=".004"):
    return alpha_given_v(src,state,z,
              wind[state][z]["wind_by_original_relative_z_fd_step"][windex]
              ["v_parallel_plus_one_LOS_sigma_kms"],kcom,pb)

def aggregated_results(outdir):
    src,p,pb=source_only_self_test()
    wind={}
    state_source_digests={}
    for state in STATES:
        path=outdir/("class_R16_wind_"+state+".json")
        raw=path.read_bytes()
        j=json.loads(raw)
        require(j["state"]==state and j["CLASS_commit_expected"]==e8.cro.CLASS_COMMIT and
                j["status"]=="CLASS_PINNED_STATE_R16_WIND_ENGINEERING_COMPLETE_NOT_HALO_REALIZATION" and
                j["frozen_source_input_CSV_sha256"]==
                p["frozen_physical_parameters"]["source_only_CSV_SHA256"] and
                j["one_sigma_positive_los_rank_not_measured_halo_wind"] is True and
                j["observed_odd_data_vector_read"] is False and
                [x["z"] for x in j["z_records"]]==list(ZS),
                "One of all THREE original E9 CLASS states is missing or non-frozen")
        wind[state]={x["z"]:x for x in j["z_records"]}
        state_source_digests[state]=e8.sha(raw)
    hzs=[.01,.005]
    cases={}
    derivative={}
    qa=[]
    for dzstep in (.004,.002):
        label=str(dzstep)
        for z in ZS:
            result={"z":z,"k_com_h_per_Mpc":.05,
                "wind_velocity_derivative_relative_z_step":dzstep,
                "state":{}}
            for state in STATES:
                rec=wind[state][z]
                v=rec["wind_by_original_relative_z_fd_step"][label][
                    "v_parallel_plus_one_LOS_sigma_kms"]
                vfd=wind["FD"][z]["wind_by_original_relative_z_fd_step"][label][
                    "v_parallel_plus_one_LOS_sigma_kms"]
                own=alpha_given_v(src,state,z,v,.05,pb)
                fixed_wind=alpha_given_v(src,state,z,vfd,.05,pb)
                result["state"][state]={
                    "own_CLASS_v_R16_1sigma_los_kms":v,
                    "counterfactual_common_FD_wind_kms":vfd,
                    "alpha_with_own_CLASS_velocity":own,
                    "alpha_with_common_FD_velocity":fixed_wind,
                    "H_physical_inverse_s":
                    rec["H_physical_SI_inverse_s"]}
            dp=result["state"]["plus"]
            dm=result["state"]["minus"]
            delta=dp["alpha_with_own_CLASS_velocity"]-dm["alpha_with_own_CLASS_velocity"]
            occ=(dp["alpha_with_common_FD_velocity"]-
                 dm["alpha_with_common_FD_velocity"])
            windshift=delta-occ
            result["delta_alpha_plus_minus_with_own_CLASS_winds"]=delta
            result["delta_alpha_occupation_at_common_FD_wind"]=occ
            result["delta_alpha_due_to_class_wind_difference_at_fixed_F_states"]=windshift
            require(math.isclose(occ+windshift,delta,abs_tol=1e-14,rel_tol=1e-12),
                    "Fixed-F frozen occupation/wind algebraic split failed")
            cases[(dzstep,z)]=result
        def val(z,name):
            return cases[(dzstep,z)]["state"][name]["alpha_with_own_CLASS_velocity"]
        def delta(z):
            return cases[(dzstep,z)]["delta_alpha_plus_minus_with_own_CLASS_winds"]
        def occ(z):
            return cases[(dzstep,z)]["delta_alpha_occupation_at_common_FD_wind"]
        def windshift(z):
            return cases[(dzstep,z)]["delta_alpha_due_to_class_wind_difference_at_fixed_F_states"]
        for hstep in hzs:
            z=.95
            deriv={
                state:float(-(1.+z)*(val(z+hstep,state)-val(z-hstep,state))/(2.*hstep))
                for state in STATES}
            deld=float(-(1.+z)*(delta(z+hstep)-delta(z-hstep))/(2.*hstep))
            occd=float(-(1.+z)*(occ(z+hstep)-occ(z-hstep))/(2.*hstep))
            windd=float(-(1.+z)*(windshift(z+hstep)-windshift(z-hstep))/(2.*hstep))
            require(math.isclose(deld,deriv["plus"]-deriv["minus"],
                                 abs_tol=2e-15,rel_tol=1e-10) and
                    math.isclose(deld,occd+windd,abs_tol=2e-15,rel_tol=1e-10),
                    "CLASS transfer quasi-static phase derivative/source+wind closure failed")
            derivkey=f"vel_{dzstep}_phase_z_{hstep}"
            derivative[derivkey]={
                "center_z":z,"k_com_h_per_Mpc":.05,
                "velocity_density_transfer_relative_z_step":dzstep,
                "phase_derivative_z_step":hstep,
                "dalpha_dln_a_by_state":deriv,
                "dalpha_plus_minus_dln_a":deld,
                "dalpha_occupation_split_dln_a":occd,
                "dalpha_class_wind_split_dln_a":windd,
                "dalpha_dt_SI_s_inverse_by_state":{
                    state:deriv[state]*wind[state][z]["H_physical_SI_inverse_s"]
                    for state in STATES},
                "normalized_to_eBOSS_xi_or_A":False}
    def err(a,b):
        return float(abs(a-b)/max(abs(a),abs(b),1e-12))
    for state in STATES:
        base=wind[state][.95]["wind_by_original_relative_z_fd_step"]
        er=err(base["0.004"]["v_parallel_plus_one_LOS_sigma_kms"],
               base["0.002"]["v_parallel_plus_one_LOS_sigma_kms"])
        qa.append({"check":state+"_velocity_fd_step","relative_gap":er,
                   "warn_over_5pct":er>.05})
        for velstep in (.004,.002):
            a=derivative[f"vel_{velstep}_phase_z_0.01"]["dalpha_dln_a_by_state"][state]
            b=derivative[f"vel_{velstep}_phase_z_0.005"]["dalpha_dln_a_by_state"][state]
            er=err(a,b)
            qa.append({"check":state+f"_phase_fd_step_with_vstep_{velstep}",
                       "relative_gap":er,"warn_over_5pct":er>.05})
    output={"status":"E9_CLASS_WIND_QUASISTATIC_PHASE_DERIVATIVE_READY_PHYSICAL_EBOSS_WINDOW_AND_GALAXY_AMPLITUDE_STOP",
        "preregistered_E9_protocol_git_blob":PROTOCOL_BLOB,
        "source_only_E8_original_CSV_sha256":
        p["frozen_physical_parameters"]["source_only_CSV_SHA256"],
        "CLASS_commit_expected":e8.cro.CLASS_COMMIT,
        "CLASS_state_JSON_sha256":state_source_digests,
        "physical_definition":"Exact Eq20 occupation, original CLASS transfer-derived R16 filtered plus-one-sigma LOS v; source states independently evolved; coherent Gaussian rank=+1 fixed over z; k=.05 h/Mpc COMOVING",
        "redshift_nodes_original_and_derivative_stencil":list(ZS),
        "R16_3d_truncated_k_to_0p10_not_full_halo_wind":True,
        "all_own_and_common_FD_wind_cases":[cases[key] for key in sorted(cases)],
        "derivative_center_0p95_two_phase_and_wind_steps":derivative,
        "prospectively_locked_engineering_QA":qa,
        "engineering_QA_warnings_not_survey_acceptance":[x for x in qa if x["warn_over_5pct"]],
        "original_Eq20_exact_FD_vs_printed_mu_approx_normalization_unresolved":True,
        "quasistatic_coherent_rank_model_not_retarded_Einstein_Vlasov":True,
        "CLASS_halo_potential_or_galaxy_tracer_calibration_computed":False,
        "A03_physical_pair_z_empirical_window_certified":False,
        "A04_eBOSS_covariance_available":False,
        "absolute_LRG_ELG_xi1_xi3_or_signal_to_noise_computed":False,
        "observed_galaxy_rows_read":False,"observed_random_rows_read":False,
        "observed_odd_data_vector_read":False,
        "new_catalogue_or_mock_download":False,
        "new_science_seeds_or_cuts":False}
    raw=canonical(output)
    atomic_new(outdir/"quasistatic_alpha_dln_a_source_only_plus_CLASS.json",raw)
    print("EBOSS_A03_E9_FROZEN_CLASS_QUASISTATIC_ALPHA_DLOGA_COMPLETE",
          "OUTPUT_SHA256",e8.sha(raw),
          "QA_WARN_COUNT",len(output["engineering_QA_warnings_not_survey_acceptance"]),
          "NO_GALAXY_XI NO_OBSERVED_ODD",flush=True)
    center=cases[(.004,.95)]
    print("E9_CENTER_OWN_WIND",json.dumps({
          s:center["state"][s]["own_CLASS_v_R16_1sigma_los_kms"]
          for s in STATES},sort_keys=True),flush=True)
    print("E9_CENTER_DALPHA_DLOGA",
          json.dumps(derivative["vel_0.004_phase_z_0.005"],sort_keys=True),flush=True)
    return 0

def main():
    a=argparse.ArgumentParser(description=__doc__)
    group=a.add_mutually_exclusive_group(required=True)
    group.add_argument("--self-test",action="store_true")
    group.add_argument("--state",choices=STATES)
    group.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    args=a.parse_args()
    if args.self_test:
        source_only_self_test()
        return 0
    if args.state:return run_state(args.output_dir,args.state)
    return aggregated_results(args.output_dir)

if __name__=="__main__":
    raise SystemExit(main())
