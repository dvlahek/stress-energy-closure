#!/usr/bin/env python3
"""K01E1R preregistered matched-phase single-mode self-gravity test."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PR=ROOT/"source_data/k01e1r_single_mode_self_gravity_revised_prereg_2026-11-06.json"
PR_BLOB="0edfd8c42f8080f72c7b2bbd9d422149411c0ff1"
E0=ROOT/"source_data/k01e0_e9_swift_cosmology_closure_result_freeze_2026-11-06.json"
E0_BLOB="7e8115d1c8779bd8e0dc1d82cef6a3c67ce1765d"
E0C=ROOT/"source_data/k01e0c_canonical_z1_pair_result_freeze_2026-11-06.json"
E0C_BLOB="80ddf7519217f49d4bf2b5ec0c4b74ad98249971"
PAIR=ROOT/"k01e0c_z1_pair/canonical_z1_moment_null_pair_4000q.csv"
PAIR_SHA="cc795bb304f8a4e2eff40f44e5332e5241cf34b018f948082a84884154ff3262"
WORK=ROOT/"k01e1r_self_gravity_work"
OUT=ROOT/"k01e1r_self_gravity_result.json"

MPC_CGS=3.0856775814913673e24
MSUN_CGS=1.98847e33
UNIT_MASS_CGS=1e10*MSUN_CGS
UNIT_VEL_CGS=1e5
UNIT_TIME_CGS=MPC_CGS/UNIT_VEL_CGS
C_KMS=299792.458
KB_EV_K=8.617333262e-5

STATES=("FD","Fplus","Fminus")
SIGNS=(+1,-1)
RESOLUTIONS=("coarse","fine")

def need(ok,msg):
    if not ok: raise RuntimeError("K01E1R_STOP: "+msg)

def blob(path):
    raw=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def loadp():
    need(blob(PR)==PR_BLOB,"prereg blob")
    need(blob(E0)==E0_BLOB,"E0 freeze blob")
    need(blob(E0C)==E0C_BLOB,"E0C freeze blob")
    p=json.loads(PR.read_text())
    need(p["stage"]=="K01E1R_SINGLE_MODE_MATCHED_PHASE_SELF_GRAVITY_REVISED_PREREG","stage")
    need(sha256(PAIR)==PAIR_SHA,"canonical pair CSV hash")
    return p

def lattice(nside,box,offset):
    x=(np.arange(nside,dtype=np.float64)+offset)*(box/nside)
    X,Y,Z=np.meshgrid(x,x,x,indexing="ij")
    return np.column_stack([X.ravel(),Y.ravel(),Z.ravel()])

def cdf_sample(q,f,u):
    pdf=q*q*f
    c=np.zeros_like(q)
    c[1:]=np.cumsum(0.5*(pdf[1:]+pdf[:-1])*(q[1:]-q[:-1]))
    need(c[-1]>0 and np.isfinite(c[-1]),"CDF norm")
    c/=c[-1]
    keep=np.concatenate(([True],np.diff(c)>0))
    return np.interp(u,c[keep],q[keep])

def mode_xy(pos,k,harm=1):
    ph=harm*k*pos[:,0]
    return float(2*np.mean(np.cos(ph))),float(2*np.mean(np.sin(ph)))

def relgap(a,b):
    return abs(a-b)/max(abs(a),abs(b),1e-300)

def build_phase_space(p,res):
    cfg=p["resolutions"][res]
    nside=int(cfg["nside_neutrino"])
    N=nside**3
    need(N%2==0,"even N")
    box=float(p["mode"]["box_Mpc"])
    pos_nu=lattice(nside,box,float(p["phase_space_design"]["neutrino_lattice_offset_in_cell_fraction"]))
    pos_cold=lattice(int(cfg["nside_cold"]),box,float(p["phase_space_design"]["cold_lattice_offset_in_cell_fraction"]))
    need(len(pos_nu)==N and len(pos_cold)==N,"equal species counts")

    seed=int(p["phase_space_design"]["neutrino_direction_seed"])
    rng=np.random.default_rng(seed)
    nb=N//2
    d=rng.normal(size=(nb,3))
    d/=np.linalg.norm(d,axis=1)[:,None]
    dirs=np.empty((N,3),dtype=np.float64)
    dirs[0::2]=d
    dirs[1::2]=-d
    perm=rng.permutation(N)
    dirs=dirs[perm]

    data=np.loadtxt(PAIR,delimiter=",",skiprows=1)
    qgrid=data[:,0]
    fs={"FD":data[:,1],"Fplus":data[:,2],"Fminus":data[:,3]}
    u=(np.arange(nb,dtype=np.float64)+0.5)/nb

    a0=float(p["cosmology"]["a_begin"])
    Tnu0=float(p["cosmology"]["T_nu_0_K"])*KB_EV_K
    mnu=float(p["cosmology"]["M_nu_eV"])
    pref=C_KMS*Tnu0/(a0*mnu)
    y=a0*mnu/Tnu0

    state={}
    for name in STATES:
        qb=cdf_sample(qgrid,fs[name],u)
        qfull=np.empty(N,dtype=np.float64)
        qfull[0::2]=qb
        qfull[1::2]=qb
        qfull=qfull[perm]
        vel=(pref*qfull)[:,None]*dirs
        eps=np.sqrt(qfull*qfull+y*y)
        state[name]={
          "q":qfull,
          "vel":vel,
          "energy_proxy":float(np.mean(eps)),
          "pressure_proxy":float(np.mean(qfull*qfull/(3*eps))),
          "mean_velocity":np.mean(vel,axis=0).tolist(),
        }
    return pos_cold,pos_nu,state

def write_ic(path,p,res,state,sign,pos_cold,pos_nu,vel_nu):
    cfg=p["resolutions"][res]
    box=float(p["mode"]["box_Mpc"])
    A=float(p["mode"]["finite_displacement_amplitude_A"])
    k=2*math.pi/box
    cold=pos_cold.copy()
    cold[:,0]=(cold[:,0]-sign*(A/k)*np.sin(k*cold[:,0]))%box
    N=len(cold)
    ids_dm=np.arange(1,N+1,dtype=np.uint64)
    ids_nu=np.arange(N+1,2*N+1,dtype=np.uint64)
    vdm=np.zeros((N,3),dtype=np.float32)
    mdm=np.full(N,float(cfg["cold_macro_mass_internal"]),dtype=np.float32)
    mnu=np.full(N,float(cfg["neutrino_macro_mass_internal"]),dtype=np.float32)

    with h5py.File(path,"w") as f:
        h=f.create_group("Header")
        counts=np.zeros(7,dtype=np.uint64);counts[1]=N;counts[6]=N
        h.attrs["BoxSize"]=np.array([box]*3,dtype=np.float64)
        h.attrs["Dimension"]=np.int32(3)
        h.attrs["NumPart_Total"]=counts
        h.attrs["NumPart_Total_HighWord"]=np.zeros(7,dtype=np.uint64)
        h.attrs["NumPart_ThisFile"]=counts
        h.attrs["NumFilesPerSnapshot"]=np.int32(1)
        h.attrs["MassTable"]=np.zeros(7,dtype=np.float64)
        h.attrs["Flag_Entropy_ICs"]=np.int32(0)
        h.attrs["Time"]=np.float64(p["cosmology"]["a_begin"])
        h.attrs["Redshift"]=np.float64(1.0/p["cosmology"]["a_begin"]-1.0)
        u=f.create_group("Units")
        u.attrs["Unit length in cgs (U_L)"]=np.float64(MPC_CGS)
        u.attrs["Unit mass in cgs (U_M)"]=np.float64(UNIT_MASS_CGS)
        u.attrs["Unit time in cgs (U_t)"]=np.float64(UNIT_TIME_CGS)
        u.attrs["Unit current in cgs (U_I)"]=np.float64(1.0)
        u.attrs["Unit temperature in cgs (U_T)"]=np.float64(1.0)
        g=f.create_group("PartType1")
        g.create_dataset("Coordinates",data=cold,dtype="f8")
        g.create_dataset("Velocities",data=vdm,dtype="f4")
        g.create_dataset("Masses",data=mdm,dtype="f4")
        g.create_dataset("ParticleIDs",data=ids_dm,dtype="u8")
        g=f.create_group("PartType6")
        g.create_dataset("Coordinates",data=pos_nu,dtype="f8")
        g.create_dataset("Velocities",data=np.asarray(vel_nu,dtype=np.float32),dtype="f4")
        g.create_dataset("Masses",data=mnu,dtype="f4")
        g.create_dataset("ParticleIDs",data=ids_nu,dtype="u8")

    return {
      "cold_C1":mode_xy(cold,k,1)[0],
      "cold_S1":mode_xy(cold,k,1)[1],
      "nu_C1":mode_xy(pos_nu,k,1)[0],
      "nu_S1":mode_xy(pos_nu,k,1)[1],
    }

def yml(p,res,state,sign):
    c=p["cosmology"];g=p["resolutions"][res]
    a1=float(c["a_end"])
    return f"""InternalUnitSystem:
  UnitMass_in_cgs:     {UNIT_MASS_CGS:.16e}
  UnitLength_in_cgs:   {MPC_CGS:.16e}
  UnitVelocity_in_cgs: {UNIT_VEL_CGS:.16e}
  UnitCurrent_in_cgs:  1.0
  UnitTemp_in_cgs:     1.0

Cosmology:
  Omega_cdm:     {c["Omega_cdm"]:.17g}
  Omega_b:       {c["Omega_b"]:.17g}
  Omega_lambda:  {c["Omega_lambda"]:.17g}
  h:             {c["h"]:.17g}
  a_begin:       {c["a_begin"]:.17g}
  a_end:         {c["a_end"]:.17g}
  N_nu:          1
  N_ur:          {c["N_ur"]:.17g}
  T_nu_0:        {c["T_nu_0_K"]:.17g}
  M_nu_eV:       {c["M_nu_eV"]:.17g}
  deg_nu:        {c["deg_nu"]:.17g}

TimeIntegration:
  dt_min: 1.0e-12
  dt_max: {p["numerical_engine"]["dt_max"]:.17g}

Scheduler:
  max_top_level_cells: 16
  tasks_per_cell: {p["numerical_engine"]["tasks_per_cell"]}
  cell_split_size: 100

Snapshots:
  scale_factor_first: {a1:.17g}
  delta_time: 2.0
  basename: e1_snap
  compression: 0

Statistics:
  scale_factor_first: {a1:.17g}
  delta_time: 2.0

InitialConditions:
  file_name: science_ic.hdf5
  periodic: 1

Neutrino:
  use_delta_f: 0
  use_delta_f_mesh_only: 0
  generate_ics: 0
  use_linear_response: 0
  use_model_none: 0

Gravity:
  mesh_side_length: {g["gravity_mesh_side"]}
  MAC: geometric
  eta: 0.025
  theta_cr: 0.5
  rebuild_frequency: 0.0
  comoving_DM_softening: {g["comoving_softening_Mpc"]:.17g}
  max_physical_DM_softening: {g["comoving_softening_Mpc"]:.17g}
  comoving_nu_softening: {g["comoving_softening_Mpc"]:.17g}
  max_physical_nu_softening: {g["comoving_softening_Mpc"]:.17g}

MetaData:
  run_name: K01E1R-{res}-{state}-{"plus" if sign>0 else "minus"}
"""

def prepare():
    p=loadp()
    WORK.mkdir(parents=True,exist_ok=True)
    qa={"pair_csv_sha256":sha256(PAIR),"resolutions":{}}
    A=float(p["mode"]["finite_displacement_amplitude_A"])
    k=2*math.pi/float(p["mode"]["box_Mpc"])

    for res in RESOLUTIONS:
        pos_cold,pos_nu,state=build_phase_space(p,res)
        qres={"states":{},"initial_modes":{}}
        Eplus=state["Fplus"]["energy_proxy"];Eminus=state["Fminus"]["energy_proxy"]
        Pplus=state["Fplus"]["pressure_proxy"];Pminus=state["Fminus"]["pressure_proxy"]
        egap=abs(Eplus-Eminus)/max(.5*(abs(Eplus)+abs(Eminus)),1e-300)
        pgap=abs(Pplus-Pminus)/max(.5*(abs(Pplus)+abs(Pminus)),1e-300)
        qres["sampled_energy_pair_gap"]=egap
        qres["sampled_pressure_pair_gap"]=pgap
        pre=p["pre_run_particle_QA"]
        if res=="fine":
            need(egap<=pre["fine_sampled_energy_pair_gap_at_most"],f"{res} energy gap {egap}")
            need(pgap<=pre["fine_sampled_pressure_pair_gap_at_most"],f"{res} pressure gap {pgap}")
        else:
            need(egap<=pre["coarse_sampled_energy_pair_gap_at_most"],f"{res} energy gap {egap}")
            need(pgap<=pre["coarse_sampled_pressure_pair_gap_at_most"],f"{res} pressure gap {pgap}")

        for state_name in STATES:
            qres["states"][state_name]={
              "energy_proxy":state[state_name]["energy_proxy"],
              "pressure_proxy":state[state_name]["pressure_proxy"],
              "mean_velocity_km_s":state[state_name]["mean_velocity"],
              "q_mean":float(np.mean(state[state_name]["q"])),
            }
            for sign in SIGNS:
                tag="plus" if sign>0 else "minus"
                run=WORK/res/state_name/tag
                run.mkdir(parents=True,exist_ok=True)
                modes=write_ic(run/"science_ic.hdf5",p,res,state_name,sign,pos_cold,pos_nu,state[state_name]["vel"])
                (run/"science.yml").write_text(yml(p,res,state_name,sign))
                qres["initial_modes"][f"{state_name}_{tag}"]=modes
                need(abs(modes["nu_C1"])<=pre["initial_neutrino_fundamental_cos_and_sin_abs_at_most"],"initial nu C1")
                need(abs(modes["nu_S1"])<=pre["initial_neutrino_fundamental_cos_and_sin_abs_at_most"],"initial nu S1")
                lo,hi=pre["cold_initial_cosine_amplitude_abs_range"]
                need(lo<=abs(modes["cold_C1"])<=hi,f"cold C1 {modes['cold_C1']}")
                need(np.sign(modes["cold_C1"])==sign,"cold mode sign")
                need(abs(modes["cold_S1"])<=pre["cold_initial_sine_amplitude_abs_at_most"],"cold S1")
        qa["resolutions"][res]=qres

    (WORK/"pre_run_qa.json").write_text(json.dumps(qa,indent=2,allow_nan=False)+"\n")
    print("K01E1R_PRE_RUN_QA_PASS",flush=True)
    for res in RESOLUTIONS:
        print("K01E1R_SAMPLE_QA",res,
              qa["resolutions"][res]["sampled_energy_pair_gap"],
              qa["resolutions"][res]["sampled_pressure_pair_gap"],flush=True)

def final_snapshot(run):
    fs=sorted(run.glob("e1_snap_*.hdf5"))
    need(fs,f"no snapshot {run}")
    vals=[]
    for path in fs:
        with h5py.File(path,"r") as f:
            if "Cosmology" in f and "Scale-factor" in f["Cosmology"].attrs:
                a=float(np.asarray(f["Cosmology"].attrs["Scale-factor"]).ravel()[0])
            else:
                a=float(np.asarray(f["Header"].attrs["Time"]).ravel()[0])
        vals.append((a,path))
    return max(vals,key=lambda x:x[0])

def read_modes(path,k):
    with h5py.File(path,"r") as f:
        dm=np.asarray(f["PartType1/Coordinates"][:],dtype=np.float64)
        nu=np.asarray(f["PartType6/Coordinates"][:],dtype=np.float64)
        idm=np.asarray(f["PartType1/ParticleIDs"][:],dtype=np.uint64)
        inu=np.asarray(f["PartType6/ParticleIDs"][:],dtype=np.uint64)
        mdm=np.asarray(f["PartType1/Masses"][:],dtype=np.float64)
        mnu=np.asarray(f["PartType6/Masses"][:],dtype=np.float64)
    need(len(np.unique(idm))==len(idm) and len(np.unique(inu))==len(inu),"IDs")
    need(np.all(mdm>0) and np.all(mnu>0),"masses")
    return {
      "cold_C1":mode_xy(dm,k,1)[0],"cold_S1":mode_xy(dm,k,1)[1],
      "cold_C2":mode_xy(dm,k,2)[0],"cold_S2":mode_xy(dm,k,2)[1],
      "nu_C1":mode_xy(nu,k,1)[0],"nu_S1":mode_xy(nu,k,1)[1],
      "nu_C2":mode_xy(nu,k,2)[0],"nu_S2":mode_xy(nu,k,2)[1],
      "mean_cold_mass":float(np.mean(mdm)),"mean_nu_mass":float(np.mean(mnu))
    }

def analyze():
    p=loadp()
    A=float(p["mode"]["finite_displacement_amplitude_A"])
    k=2*math.pi/float(p["mode"]["box_Mpc"])
    raw={}
    responses={}

    for res in RESOLUTIONS:
        raw[res]={};responses[res]={}
        for state in STATES:
            raw[res][state]={}
            for sign in SIGNS:
                tag="plus" if sign>0 else "minus"
                a,snap=final_snapshot(WORK/res/state/tag)
                need(abs(a-float(p["cosmology"]["a_end"]))<2e-6,f"final a {a}")
                raw[res][state][tag]=read_modes(snap,k)
            xp=raw[res][state]["plus"];xm=raw[res][state]["minus"]
            responses[res][state]={
              "Rcos":(xp["nu_C1"]-xm["nu_C1"])/(2*A),
              "Rsin":(xp["nu_S1"]-xm["nu_S1"])/(2*A),
              "R2cos":(xp["nu_C2"]-xm["nu_C2"])/(2*A),
              "cold_Rcos":(xp["cold_C1"]-xm["cold_C1"])/(2*A),
              "cold_Rsin":(xp["cold_S1"]-xm["cold_S1"])/(2*A),
            }

    diag={}
    for res in RESOLUTIONS:
        r=responses[res]
        H=r["Fminus"]["Rcos"]-r["Fplus"]["Rcos"]
        Hsin=r["Fminus"]["Rsin"]-r["Fplus"]["Rsin"]
        H2=r["Fminus"]["R2cos"]-r["Fplus"]["R2cos"]
        floor=max(abs(Hsin),abs(H2),1e-15)
        diag[res]={
          "hidden_contrast_H":H,
          "hidden_relative_to_abs_FD":H/max(abs(r["FD"]["Rcos"]),1e-300),
          "hidden_sine_null":Hsin,
          "hidden_2k_odd_null":H2,
          "null_floor":floor,
          "hidden_to_null_floor":abs(H)/floor,
          "pair_mean_minus_FD":0.5*(r["Fplus"]["Rcos"]+r["Fminus"]["Rcos"])-r["FD"]["Rcos"],
          "cold_Rcos_cross_state_range":max(v["cold_Rcos"] for v in r.values())-min(v["cold_Rcos"] for v in r.values()),
        }

    Hc=diag["coarse"]["hidden_contrast_H"];Hf=diag["fine"]["hidden_contrast_H"]
    hidden_gap=relgap(Hc,Hf)
    state_gaps={s:relgap(responses["coarse"][s]["Rcos"],responses["fine"][s]["Rcos"]) for s in STATES}
    g=p["preregistered_science_gate"]
    same_sign=Hc!=0 and Hf!=0 and np.sign(Hc)==np.sign(Hf)

    gate=bool(
      responses["coarse"]["FD"]["Rcos"]>0
      and responses["fine"]["FD"]["Rcos"]>0
      and same_sign
      and hidden_gap<=g["hidden_contrast_coarse_fine_relative_gap_at_most"]
      and max(state_gaps.values())<=g["each_state_Rcos_coarse_fine_relative_gap_at_most"]
      and diag["fine"]["hidden_to_null_floor"]>=g["fine_hidden_contrast_to_null_floor_at_least"]
      and diag["coarse"]["hidden_to_null_floor"]>=g["coarse_hidden_contrast_to_null_floor_at_least"]
    )

    result={
      "stage":p["stage"],
      "status":"PASS_NUMERICALLY_RESOLVED_MATCHED_PHASE_SELF_GRAVITY_HIDDEN_RESPONSE" if gate else "FAIL_K01E1R_PREREGISTERED_SELF_GRAVITY_GATE",
      "science_gate_pass":gate,
      "raw_final_modes":raw,
      "responses":responses,
      "diagnostics":diag,
      "resolution_checks":{
        "hidden_contrast_same_nonzero_sign":bool(same_sign),
        "hidden_contrast_coarse_fine_relative_gap":hidden_gap,
        "state_Rcos_coarse_fine_relative_gaps":state_gaps,
      },
      "interpretation":p["interpretation_boundary"],
      "stopping_rule":p["stopping_rule"],
      "guards":p["guards"],
      "next_stage_if_pass":p["next_if_pass"]
    }
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for res in RESOLUTIONS:
        print("K01E1R_RESPONSE",res,json.dumps(responses[res],sort_keys=True),flush=True)
        print("K01E1R_HIDDEN",res,json.dumps(diag[res],sort_keys=True),flush=True)
    print("K01E1R_HIDDEN_RESOLUTION_GAP",hidden_gap,flush=True)
    print("K01E1R_STATE_RESOLUTION_GAPS",state_gaps,flush=True)
    print("K01E1R_SCIENCE_GATE_PASS",gate,flush=True)
    print("K01E1R_STATUS",result["status"],flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("mode",choices=("prepare","analyze"))
    args=ap.parse_args()
    if args.mode=="prepare":prepare()
    else:analyze()

if __name__=="__main__":main()
