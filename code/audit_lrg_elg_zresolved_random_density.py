#!/usr/bin/env python3
"""Audit random-catalog density using paired z-resolved r1/r4 EZmocks."""
from __future__ import annotations
import argparse, glob, json
from pathlib import Path
import numpy as np

FILENAME="lrg_elg_exact_zresolved_odd_multipoles.csv"

def load(path):
    p=Path(path)
    if p.is_dir(): p=p/FILENAME
    a=np.genfromtxt(p,delimiter=",",names=True)
    keys=np.column_stack([a["zlo"],a["zhi"],a["s_Mpc_over_h"]]).astype(float)
    x1=np.asarray(a["xi1_LRG_to_ELG"],float)
    x3=np.asarray(a["xi3_LRG_to_ELG"],float)
    return keys,x1,x3

def discover(pattern):
    out={}
    for p in sorted(glob.glob(pattern)):
        pp=Path(p)
        d=pp if pp.is_dir() else pp.parent
        name=d.name
        try: idx=int(name.split("_")[-1])
        except Exception: continue
        out[idx]=d
    return out

def stats(A,B):
    D=B-A
    meanB=B.mean(axis=0)
    survey=B-meanB
    varD=np.var(D,axis=0,ddof=1)
    varB=np.var(B,axis=0,ddof=1)
    good=varB>0
    ratios=np.sqrt(varD[good]/varB[good])
    return {
        "global_rms_r4_minus_r1":float(np.sqrt(np.mean(D**2))),
        "global_rms_r4_mock_scatter":float(np.sqrt(np.mean(survey**2))),
        "global_rms_ratio":float(np.sqrt(np.mean(D**2))/np.sqrt(np.mean(survey**2))),
        "trace_covariance_ratio_delta_over_r4":float(np.sum(varD)/np.sum(varB)),
        "component_sigma_ratio_median":float(np.median(ratios)),
        "component_sigma_ratio_max":float(np.max(ratios)),
        "mean_shift_rms":float(np.sqrt(np.mean(D.mean(axis=0)**2))),
        "per_component_sigma_ratio":ratios.tolist(),
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--r1-glob",default="mocks_zresolved_desi_pair_r1/exact_z_*")
    ap.add_argument("--r4-glob",default="mocks_zresolved_desi_pair_r4/exact_z_*")
    ap.add_argument("--data-r1")
    ap.add_argument("--data-r4")
    ap.add_argument("--out",default="lrg_elg_zresolved_random_density_audit.json")
    args=ap.parse_args()

    R1,R4=discover(args.r1_glob),discover(args.r4_glob)
    ids=sorted(set(R1)&set(R4))
    if len(ids)<3: raise RuntimeError("Need at least three paired mocks")

    k0=None; a1=[]; a3=[]; b1=[]; b3=[]
    for i in ids:
        k,x1,x3=load(R1[i]); kb,y1,y3=load(R4[i])
        if not np.allclose(k,kb,rtol=0,atol=1e-12): raise RuntimeError(f"grid mismatch mock {i}")
        if k0 is None: k0=k
        elif not np.allclose(k,k0,rtol=0,atol=1e-12): raise RuntimeError(f"grid mismatch mock {i}")
        a1.append(x1);a3.append(x3);b1.append(y1);b3.append(y3)
    a1=np.asarray(a1);a3=np.asarray(a3);b1=np.asarray(b1);b3=np.asarray(b3)

    result={
        "scope":"Paired random-density audit for DESI-fiducial z-resolved estimator",
        "n_pairs":len(ids),
        "mock_ids":ids,
        "dipole":stats(a1,b1),
        "octupole":stats(a3,b3),
        "joint":stats(np.concatenate([a1,a3],axis=1),np.concatenate([b1,b3],axis=1)),
        "interpretation_guardrail":"This quantifies extra variance from using one instead of four random realizations. It does not replace the final survey covariance."
    }

    if args.data_r1 and args.data_r4:
        kd,d11,d13=load(args.data_r1); k4,d41,d43=load(args.data_r4)
        if not np.allclose(kd,k4,rtol=0,atol=1e-12) or not np.allclose(kd,k0,rtol=0,atol=1e-12):
            raise RuntimeError("data/mock grids differ")
        for name,dr1,dr4,mr4 in [
            ("dipole",d11,d41,b1),("octupole",d13,d43,b3),
            ("joint",np.r_[d11,d13],np.r_[d41,d43],np.concatenate([b1,b3],axis=1))
        ]:
            sd=np.std(mr4,axis=0,ddof=1)
            delta=dr4-dr1
            good=sd>0
            result.setdefault("data_r1_to_r4",{})[name]={
                "rms_shift":float(np.sqrt(np.mean(delta**2))),
                "median_abs_shift_in_component_sigma":float(np.median(np.abs(delta[good])/sd[good])),
                "max_abs_shift_in_component_sigma":float(np.max(np.abs(delta[good])/sd[good])),
            }

    Path(args.out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
