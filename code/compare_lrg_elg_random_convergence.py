#!/usr/bin/env python3
"""Compare exact LRGxELG odd-multipole vectors across random-catalog densities."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

def load(path):
    a=np.genfromtxt(Path(path)/"lrg_elg_exact_odd_multipoles.csv",delimiter=",",names=True)
    return (np.asarray(a["s_Mpc_over_h"],float),
            np.asarray(a["xi1_LRG_to_ELG"],float),
            np.asarray(a["xi3_LRG_to_ELG"],float))

def metrics(ref, cur):
    d=cur-ref
    rms=lambda x: float(np.sqrt(np.mean(np.asarray(x,float)**2)))
    return {
        "max_abs_change": float(np.max(np.abs(d))),
        "rms_change": rms(d),
        "relative_rms_change_to_reference": float(rms(d)/max(rms(ref),1e-300))
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--runs",nargs="+",required=True,help="Ordered output directories, lowest to highest random density.")
    ap.add_argument("--out",default="lrg_elg_random_convergence.json")
    args=ap.parse_args()
    vals=[]
    for d in args.runs:
        s,x1,x3=load(d); vals.append((d,s,x1,x3))
    sref=vals[-1][1]; x1ref=vals[-1][2]; x3ref=vals[-1][3]
    rows=[]
    for d,s,x1,x3 in vals:
        if not np.allclose(s,sref,rtol=0,atol=1e-12): raise RuntimeError(f"Separation grid differs for {d}")
        rows.append({"run":d,"dipole_vs_highest_random":metrics(x1ref,x1),
                     "octupole_vs_highest_random":metrics(x3ref,x3)})
    out={"reference_run":vals[-1][0],"runs":rows,
         "guardrail":"Random convergence is established only if changes become negligible compared with the final statistical uncertainty from mocks."}
    Path(args.out).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__": main()
