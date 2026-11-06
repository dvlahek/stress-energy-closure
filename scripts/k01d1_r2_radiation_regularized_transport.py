#!/usr/bin/env python3
"""K01D1-R2 radiation-regularized SWIFT transport-only isolation.

R1 was blocked before particle dynamics because its artificial pure-de-Sitter
background made SWIFT's age integral from a=0 non-integrable. R2 changes only
that engine-isolation background to Omega_r=1e-4, Omega_Lambda=0.9999.
All custom-F particle design and transport gates are inherited unchanged.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PREREG=ROOT/"source_data/k01d1_r2_radiation_regularized_transport_prereg_2026-11-06.json"
PREREG_BLOB="a8eeff05350859b0680a09d499472dcf96857a3f"
WORK=ROOT/"k01d1_r2_swift_work"
OUT=ROOT/"k01d1_r2_transport_result.json"

sys.path.insert(0,str(ROOT/"scripts"))
import k01d1_r1_swift_transport_only_engine_isolation as r1


def need(ok,msg):
    if not ok:
        raise RuntimeError("K01D1_R2_STOP: "+msg)


def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def load_prereg():
    raw=PREREG.read_bytes()
    need(git_blob(raw)==PREREG_BLOB,"preregistration changed")
    p=json.loads(raw)
    need(p["stage"]=="K01D1_R2_RADIATION_REGULARIZED_TRANSPORT_ONLY_PREREG","stage")
    c=p["artificial_transport_background"]
    need(c["Omega_cdm"]==0 and c["Omega_b"]==0,"matter background drift")
    need(c["Omega_r"]==1e-4 and c["Omega_lambda"]==0.9999,"radiation regularization drift")
    need(c["N_nu"]==0,"N_nu drift")
    need(p["engine_modes"]["self_gravity"] is False,"self-gravity forbidden")
    return p


def inherited_p(p):
    """Build R1-compatible parameter object with only background regularized."""
    return {
      "stage":p["stage"],
      "artificial_transport_background":{
        "h":p["artificial_transport_background"]["h"],
        "Omega_cdm":0.0,
        "Omega_b":0.0,
        "Omega_lambda":0.9999,
        "N_nu":0,
        "N_ur":0.0,
        "a_begin":p["artificial_transport_background"]["a_begin"],
        "a_end":p["artificial_transport_background"]["a_end"],
      },
      "custom_particle_kinematics":{
        "physical_T_nu_0_K_used_to_build_IC_velocities":1.951757805,
        "physical_M_nu_eV_used_to_build_IC_velocities":0.06,
        "macro_particle_mass_internal_units":1.0,
      },
      "particle_design":{
        "particles_per_state":1024,
        "box_comoving_Mpc":512.0,
        "seed_positions_directions":2026110601,
      },
      "preregistered_acceptance":p["preregistered_acceptance"],
      "interpretation_boundary":p["interpretation_boundary"],
      "guards":p["guards"],
      "stopping_and_next":p["stopping_and_next"],
      "swift":p["swift"],
    }


def yaml_text(p,state):
    c=p["artificial_transport_background"]
    a0=float(c["a_begin"]); a1=float(c["a_end"])
    return f"""InternalUnitSystem:
  UnitMass_in_cgs:     {r1.UNIT_MASS_CGS:.16e}
  UnitLength_in_cgs:   {r1.MPC_CGS:.16e}
  UnitVelocity_in_cgs: {r1.UNIT_VEL_CGS:.16e}
  UnitCurrent_in_cgs:  1.0000000000000000e+00
  UnitTemp_in_cgs:     1.0000000000000000e+00

Cosmology:
  Omega_cdm:     0.0
  Omega_b:       0.0
  Omega_r:       1.0e-4
  Omega_lambda:  0.9999
  h:             {c["h"]:.16g}
  a_begin:       {a0:.16g}
  a_end:         {a1:.16g}
  N_nu:          0

TimeIntegration:
  dt_min:    1.0000e-12
  dt_max:    1.0000e-03

Scheduler:
  max_top_level_cells: 16
  cell_split_size: 100

Snapshots:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0
  basename: r2_snap
  compression: 0

Statistics:
  scale_factor_first: {a1:.16g}
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

ConstantPotential:
  g_cgs: [0.0, 0.0, 0.0]

Gravity:
  mesh_side_length: 32
  MAC: geometric
  eta: 0.025
  theta_cr: 0.5
  rebuild_frequency: 0.0
  comoving_DM_softening: 0.04
  max_physical_DM_softening: 0.04
  comoving_nu_softening: 0.04
  max_physical_nu_softening: 0.04

MetaData:
  run_name: K01D1-R2-{state}
"""


def prepare(p):
    ip=inherited_p(p)
    src=r1.load_e8()
    qgrid=np.asarray(src["q"],dtype=np.float64)
    fs={"FD":np.asarray(src["f0"]),"Fplus":np.asarray(src["fp"]),"Fminus":np.asarray(src["fm"])}
    n=1024
    pos,dirs,ids=r1.deterministic_geometry(ip)
    a0=float(ip["artificial_transport_background"]["a_begin"])
    T0=1.951757805*r1.KB_EV_K
    pref=r1.C_KMS*T0/(a0*0.06)

    summary={}
    for name in ("FD","Fplus","Fminus"):
        q=r1.q_samples(qgrid,fs[name],n)
        vel=(pref*q)[:,None]*dirs
        mass=np.ones(n,dtype=np.float64)
        run=WORK/name
        run.mkdir(parents=True,exist_ok=True)
        r1.write_ic(run/"science_ic.hdf5",ip,pos,vel,mass,ids)
        (run/"science.yml").write_text(yaml_text(ip,name))
        np.savez_compressed(run/"input_reference.npz",
            ids=ids,positions=pos,directions=dirs,q=q,velocities=vel,masses=mass)
        summary[name]={
          "q_mean":float(q.mean()),
          "mean_input_speed_km_s":float(np.linalg.norm(vel,axis=1).mean()),
          "max_input_speed_km_s":float(np.linalg.norm(vel,axis=1).max()),
        }
        print("K01D1_R2_INPUT",name,summary[name],flush=True)
    (WORK/"input_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("K01D1_R2_INPUTS_READY",flush=True)


def analyze(p):
    # Reuse the already preregistered R1 analysis logic, only redirecting
    # work/output paths and final snapshot basename.
    ip=inherited_p(p)
    old_work,old_out,old_find=r1.WORK,r1.OUT,r1.find_final_snapshot

    def find_final_snapshot_r2(run):
        files=sorted(run.glob("r2_snap_*.hdf5"))
        r1.need(files,f"no snapshot in {run}")
        rec=[]
        import h5py
        for path in files:
            with h5py.File(path,"r") as f:
                if "Cosmology" in f and "Scale-factor" in f["Cosmology"].attrs:
                    a=float(np.asarray(f["Cosmology"].attrs["Scale-factor"]).ravel()[0])
                elif "Scale-factor" in f["Header"].attrs:
                    a=float(np.asarray(f["Header"].attrs["Scale-factor"]).ravel()[0])
                else:
                    a=float(np.asarray(f["Header"].attrs["Time"]).ravel()[0])
            rec.append((a,path))
        return max(rec,key=lambda z:z[0])

    try:
        r1.WORK=WORK
        r1.OUT=OUT
        r1.find_final_snapshot=find_final_snapshot_r2
        r1.analyze(ip)
    finally:
        r1.WORK, r1.OUT, r1.find_final_snapshot = old_work,old_out,old_find

    data=json.loads(OUT.read_text())
    data["stage"]=p["stage"]
    data["status"]="PASS_SWIFT_CUSTOM_F_RADIATION_REGULARIZED_TRANSPORT_ONLY"
    data["artificial_transport_background"]=p["artificial_transport_background"]
    data["interpretation"]=p["interpretation_boundary"]
    data["next_stage"]=p["stopping_and_next"]["if_pass"]
    OUT.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    print("K01D1_R2_SWIFT_CUSTOM_F_TRANSPORT_ONLY_PASS",flush=True)
    print("K01D1_R2_RESULT",OUT,flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("mode",choices=("prepare","analyze"))
    args=ap.parse_args()
    p=load_prereg()
    WORK.mkdir(parents=True,exist_ok=True)
    if args.mode=="prepare": prepare(p)
    else: analyze(p)

if __name__=="__main__":
    main()
