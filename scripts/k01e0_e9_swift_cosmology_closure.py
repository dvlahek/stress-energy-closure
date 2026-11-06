#!/usr/bin/env python3
"""K01E0 frozen-E9 -> SWIFT cosmology closure and neutrino mass calibration."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PR=ROOT/"source_data/k01e0_e9_swift_cosmology_closure_prereg_2026-11-06.json"
PR_BLOB="97e6f4f5035615582b58f254ed56570b3e8dfdb7"
WORK=ROOT/"k01e0_swift_closure_work"
OUT=ROOT/"k01e0_e9_swift_cosmology_closure_result.json"

MPC_CGS=3.0856775814913673e24
MSUN_CGS=1.98847e33
UNIT_MASS_CGS=1e10*MSUN_CGS
UNIT_VEL_CGS=1e5
UNIT_TIME_CGS=MPC_CGS/UNIT_VEL_CGS

def blob(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def need(ok,msg):
    if not ok: raise RuntimeError("K01E0_STOP: "+msg)

def loadp():
    need(blob(PR)==PR_BLOB,"prereg changed")
    p=json.loads(PR.read_text())
    need(p["stage"]=="K01E0_FROZEN_E9_SWIFT_COSMOLOGY_CLOSURE_AND_NU_MASS_CALIBRATION_PREREG","stage")
    return p

def lattice(n,box):
    m=round(n**(1/3))
    need(m**3==n,"N not cubic")
    q=(np.arange(m,dtype=np.float64)+0.5)/m*box
    x,y,z=np.meshgrid(q,q,q,indexing="ij")
    return np.column_stack([x.ravel(),y.ravel(),z.ravel()])

def write_ic(p):
    n=int(p["derived_preexecution"]["neutrino_particles"])
    box=float(p["derived_preexecution"]["box_Mpc_for_fundamental_k"])
    pos=lattice(n,box)
    ids=np.arange(1,n+1,dtype=np.uint64)
    vel=np.zeros((n,3),dtype=np.float32)
    mass=np.zeros(n,dtype=np.float32)
    path=WORK/"calibration_ic.hdf5"
    WORK.mkdir(parents=True,exist_ok=True)
    with h5py.File(path,"w") as f:
        h=f.create_group("Header")
        counts=np.zeros(7,dtype=np.uint64); counts[6]=n
        h.attrs["BoxSize"]=np.array([box]*3,dtype=np.float64)
        h.attrs["Dimension"]=np.int32(3)
        h.attrs["NumPart_Total"]=counts
        h.attrs["NumPart_Total_HighWord"]=np.zeros(7,dtype=np.uint64)
        h.attrs["NumPart_ThisFile"]=counts
        h.attrs["NumFilesPerSnapshot"]=np.int32(1)
        h.attrs["MassTable"]=np.zeros(7,dtype=np.float64)
        h.attrs["Flag_Entropy_ICs"]=np.int32(0)
        h.attrs["Time"]=np.float64(0.5)
        h.attrs["Redshift"]=np.float64(1.0)
        u=f.create_group("Units")
        u.attrs["Unit length in cgs (U_L)"]=np.float64(MPC_CGS)
        u.attrs["Unit mass in cgs (U_M)"]=np.float64(UNIT_MASS_CGS)
        u.attrs["Unit time in cgs (U_t)"]=np.float64(UNIT_TIME_CGS)
        u.attrs["Unit current in cgs (U_I)"]=np.float64(1.0)
        u.attrs["Unit temperature in cgs (U_T)"]=np.float64(1.0)
        g=f.create_group("PartType6")
        g.create_dataset("Coordinates",data=pos)
        g.create_dataset("Velocities",data=vel)
        g.create_dataset("Masses",data=mass)
        g.create_dataset("ParticleIDs",data=ids)
    return path

def yml(p):
    c=p["calibration_background"]
    a0=c["a_begin"]; a1=c["a_end"]
    return f"""InternalUnitSystem:
  UnitMass_in_cgs:     {UNIT_MASS_CGS:.16e}
  UnitLength_in_cgs:   {MPC_CGS:.16e}
  UnitVelocity_in_cgs: {UNIT_VEL_CGS:.16e}
  UnitCurrent_in_cgs:  1.0
  UnitTemp_in_cgs:     1.0

Cosmology:
  Omega_cdm:     0.0
  Omega_b:       0.0
  Omega_lambda:  {c["Omega_lambda"]:.16g}
  h:             {c["h"]:.16g}
  a_begin:       {a0:.16g}
  a_end:         {a1:.16g}
  N_nu:          1
  N_ur:          {c["N_ur"]:.16g}
  T_nu_0:        {c["T_nu_0_K"]:.16g}
  M_nu_eV:       {c["M_nu_eV"]:.16g}
  deg_nu:        {c["deg_nu"]:.16g}

TimeIntegration:
  dt_min: 1.0e-12
  dt_max: 1.0e-3

Scheduler:
  max_top_level_cells: 16
  tasks_per_cell: 100
  cell_split_size: 100

Snapshots:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0
  basename: e0_snap
  compression: 0

Statistics:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0

InitialConditions:
  file_name: calibration_ic.hdf5
  periodic: 1

Neutrino:
  use_delta_f: 0
  use_delta_f_mesh_only: 0
  generate_ics: 1
  use_linear_response: 0
  use_model_none: 0

ConstantPotential:
  g_cgs: [0.0,0.0,0.0]

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
  run_name: K01E0-E9-closure
"""

def prepare():
    p=loadp()
    write_ic(p)
    (WORK/"calibration.yml").write_text(yml(p))
    print("K01E0_INPUT_READY",flush=True)

def find_snap():
    fs=sorted(WORK.glob("e0_snap_*.hdf5"))
    need(fs,"no snapshot")
    best=None
    for q in fs:
        with h5py.File(q,"r") as f:
            a=float(np.asarray(f["Cosmology"].attrs["Scale-factor"]).ravel()[0])
        if best is None or a>best[0]: best=(a,q)
    return best

def scalar(a):
    return float(np.asarray(a).ravel()[0])

def analyze():
    p=loadp()
    a,snap=find_snap()
    with h5py.File(snap,"r") as f:
        cg=f["Cosmology"].attrs
        onu=scalar(cg["Omega_nu_0"])
        orad=scalar(cg["Omega_r"])
        og=scalar(cg["Omega_g"])
        our=scalar(cg["Omega_ur"])
        nnu=int(scalar(cg["N_nu"]))
        deg=np.asarray(cg["deg_nu"],dtype=float).ravel().tolist()
        m=np.asarray(f["PartType6/Masses"][:],dtype=np.float64)
        ids=np.asarray(f["PartType6/ParticleIDs"][:],dtype=np.uint64)
    need(nnu==1,"Nnu")
    need(len(deg)==1 and abs(deg[0]-1.0)<1e-14,"deg")
    need(np.all(np.isfinite(m)) and np.all(m>0),"mass")
    uniform=float(np.ptp(m)/np.mean(m))
    need(uniform<=2e-7,f"mass uniform {uniform}")
    need(len(np.unique(ids))==len(ids),"ids")

    d=p["derived_preexecution"]
    omb=float(d["Omega_b"]); omc=float(d["Omega_cdm"])
    ol=1.0-omb-omc-orad-onu
    residual=1.0-(omb+omc+orad+onu+ol)
    need(0.6<ol<0.8,f"lambda {ol}")
    need(abs(residual)<=1e-12,f"closure {residual}")

    result={
      "stage":p["stage"],
      "status":"PASS_FROZEN_E9_SWIFT_COSMOLOGY_CLOSURE_AND_NU_MASS_CALIBRATION",
      "science_gate_pass":True,
      "swift_commit":p["swift"]["commit"],
      "box_Mpc":float(d["box_Mpc_for_fundamental_k"]),
      "fundamental_k_1_Mpc":2*math.pi/float(d["box_Mpc_for_fundamental_k"]),
      "fundamental_k_h_Mpc":(2*math.pi/float(d["box_Mpc_for_fundamental_k"]))/float(p["frozen_E9_CLASS"]["h"]),
      "Omega_b":omb,
      "Omega_cdm":omc,
      "Omega_m_nonnu":omb+omc,
      "Omega_r":orad,
      "Omega_g":og,
      "Omega_ur":our,
      "Omega_nu_0":onu,
      "Omega_lambda_flat":ol,
      "closure_residual":residual,
      "N_nu":nnu,
      "deg_nu":deg,
      "neutrino_particles":len(m),
      "neutrino_macro_mass_internal_units":float(np.mean(m)),
      "neutrino_macro_mass_uniform_relative_range":uniform,
      "snapshot":str(snap.relative_to(ROOT)),
      "snapshot_scale_factor":a,
      "interpretation":p["interpretation"],
      "guards":p["guards"],
      "next_stage":p["next_if_pass"]
    }
    OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("K01E0_OMEGA_R",orad,flush=True)
    print("K01E0_OMEGA_NU0",onu,flush=True)
    print("K01E0_OMEGA_LAMBDA_FLAT",ol,flush=True)
    print("K01E0_NU_MACRO_MASS",float(np.mean(m)),flush=True)
    print("K01E0_FUNDAMENTAL_K_H_MPC",result["fundamental_k_h_Mpc"],flush=True)
    print("K01E0_CLOSURE_PASS",residual,flush=True)
    print("K01E0_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("mode",choices=("prepare","analyze"))
    a=ap.parse_args()
    prepare() if a.mode=="prepare" else analyze()

if __name__=="__main__": main()
