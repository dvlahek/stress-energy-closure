#!/usr/bin/env python3
"""K01E-D0 audit persistence of E8 n,rho,P matching away from z=1100."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
sys.path.insert(0,str(ROOT/"code"))
import build_eboss_dr16_a03_e8_frozen_resonant_occupation as e8
import class_response_optimize as cro

DECL=ROOT/"source_data/k01e_d0_e8_stress_energy_persistence_audit_2026-11-06.json"
OUT=ROOT/"k01e_d0_e8_stress_energy_persistence_result.json"

d=json.loads(DECL.read_text())
e8.source_gate()
src=e8.kinetic_source()
q=np.asarray(src["q"],dtype=float)
states={
 "FD":np.asarray(src["f0"],dtype=float),
 "Fplus":np.asarray(src["fp"],dtype=float),
 "Fminus":np.asarray(src["fm"],dtype=float),
}

def relsym(a,b):
    return np.abs(a-b)/np.maximum(0.5*(np.abs(a)+np.abs(b)),1e-300)

res={}
for z in d["inputs"]["audit_redshifts"]:
    _,weights,_,_,_,_,_=cro.kinetic_objects(q,e8.MASS_EV,float(z))
    m={name:np.asarray(cro.moments(f,q,weights),dtype=float) for name,f in states.items()}
    gpm=relsym(m["Fplus"],m["Fminus"])
    gp0=relsym(m["Fplus"],m["FD"])
    gm0=relsym(m["Fminus"],m["FD"])
    res[str(z)]={
      "FD_n_rho_P":m["FD"].tolist(),
      "Fplus_n_rho_P":m["Fplus"].tolist(),
      "Fminus_n_rho_P":m["Fminus"].tolist(),
      "Fplus_vs_Fminus_symmetric_relative_gaps":gpm.tolist(),
      "max_Fplus_vs_Fminus_gap":float(np.max(gpm)),
      "Fplus_vs_FD_symmetric_relative_gaps":gp0.tolist(),
      "Fminus_vs_FD_symmetric_relative_gaps":gm0.tolist(),
    }
    print("K01E_D0_Z",z,"FPM_GAPS",gpm.tolist(),"MAX",float(np.max(gpm)),flush=True)

z1=res["1.0"]["max_Fplus_vs_Fminus_gap"]
classification=(
  "PERSISTS_TO_Z1_AT_1E12_LEVEL"
  if z1<1e-12 else
  "DOES_NOT_PERSIST_TO_Z1_REQUIRE_PROSPECTIVE_Z1_NULL_PAIR"
)
out={
 "stage":d["stage"],
 "status":"DIAGNOSTIC_COMPLETE",
 "classification":classification,
 "results":res,
 "interpretation":d["interpretation"],
 "guards":d["guards"],
}
OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
print("K01E_D0_CLASSIFICATION",classification,flush=True)
