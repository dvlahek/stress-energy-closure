#!/usr/bin/env python3
"""Empirical finite-mock audit for the frozen z-resolved wake matched filter.

For each mock i:
  1. remove mock i from the covariance/mean estimate,
  2. project the same fixed wake template against the same fixed nuisance shape,
  3. score held-out mock i with the resulting Hartlap precision.

The observed-data score uses all mocks. The resulting empirical p+1 is a
finite-ensemble calibration of the fixed matched-filter statistic. With only
40 mocks its minimum attainable p is 1/41 = 0.02439, so it cannot by itself
establish a two-sided significance above about 2.25 sigma.

Also reports the Sellentin-Heavens modified-t likelihood ratio using the
unscaled sample-covariance inverse:
  -2 log Lambda = N log[(1+q0/(N-1))/(1+q1/(N-1))].
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

FILENAME = "lrg_elg_exact_zresolved_odd_multipoles.csv"


def load_measurement(path):
    p=Path(path)
    if p.is_dir(): p=p/FILENAME
    a=np.genfromtxt(p,delimiter=",",names=True)
    keys=np.column_stack([a["zlo"],a["zhi"],a["s_Mpc_over_h"]]).astype(float)
    x=np.asarray(a["xi1_LRG_to_ELG"],float)
    return p,keys,x


def load_template(path,col):
    a=np.genfromtxt(path,delimiter=",",names=True)
    keys=np.column_stack([a["zlo"],a["zhi"],a["s_Mpc_over_h"]]).astype(float)
    return keys,np.asarray(a[col],float)


def hartlap(n,p):
    return float((n-p-2)/(n-1)) if n>p+2 else None


def projected_stat(y,C,n,signal,nuisance):
    inv=np.linalg.pinv(C,rcond=1e-12)
    a=hartlap(n,len(y))
    if a is None:
        raise RuntimeError(f"not enough mocks: n={n}, p={len(y)}")
    P=a*inv

    nuisance=np.asarray(nuisance,float)
    if nuisance.size==0:
        nuisance=np.empty((len(y),0))
    if nuisance.ndim==1:
        nuisance=nuisance[:,None]

    if nuisance.shape[1]:
        G=nuisance.T@P@nuisance
        proj=np.linalg.pinv(G,rcond=1e-12)@(nuisance.T@P@signal)
        sperp=signal-nuisance@proj
        bn=np.linalg.pinv(G,rcond=1e-12)@(nuisance.T@P@y)
        r0=y-nuisance@bn
    else:
        sperp=signal.copy()
        r0=y.copy()

    den=float(sperp@P@sperp)
    z=float((sperp@P@y)/np.sqrt(den))
    A=float((sperp@P@y)/den)
    sigma=float(1.0/np.sqrt(den))

    q0_h=float(r0@P@r0)
    q0_raw=q0_h/a
    # Adding signal; because Hartlap is a scalar, GLS coefficients are
    # identical for raw and Hartlap-scaled precision.
    X=np.column_stack([nuisance,signal])
    invF=np.linalg.pinv(X.T@P@X,rcond=1e-12)
    b=invF@(X.T@P@y)
    r1=y-X@b
    q1_h=float(r1@P@r1)
    q1_raw=q1_h/a

    sh=float(n*np.log((1.0+q0_raw/(n-1.0))/(1.0+q1_raw/(n-1.0))))
    return {
        "z_signed":z,
        "abs_z":abs(z),
        "amplitude":A,
        "sigma_hartlap":sigma,
        "delta_chi2_hartlap":float(q0_h-q1_h),
        "q0_raw":q0_raw,
        "q1_raw":q1_raw,
        "minus2_log_lambda_sellentin_heavens":sh,
        "sqrt_minus2_log_lambda_SH_diagnostic":float(np.sqrt(max(sh,0.0))),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measurement",required=True)
    ap.add_argument("--mock-glob",required=True)
    ap.add_argument("--signal-template",required=True)
    ap.add_argument("--signal-col",default="wake_shape")
    ap.add_argument("--nuisance-template")
    ap.add_argument("--nuisance-cols",default="")
    ap.add_argument("--out",required=True)
    args=ap.parse_args()

    _,kd,data=load_measurement(args.measurement)
    ks,signal=load_template(args.signal_template,args.signal_col)
    if not np.allclose(kd,ks,rtol=0,atol=1e-12):
        raise RuntimeError("data/signal grids differ")

    names=[x.strip() for x in args.nuisance_cols.split(",") if x.strip()]
    nt=args.nuisance_template or args.signal_template
    nuisance=[]
    for name in names:
        kn,v=load_template(nt,name)
        if not np.allclose(kd,kn,rtol=0,atol=1e-12):
            raise RuntimeError(f"data/nuisance grid differs for {name}")
        nuisance.append(v)
    Nmat=np.column_stack(nuisance) if nuisance else np.empty((len(data),0))

    paths=sorted(glob.glob(args.mock_glob))
    if len(paths)<4:
        raise RuntimeError("need at least four mocks")
    X=[]
    for path in paths:
        _,k,x=load_measurement(path)
        if not np.allclose(kd,k,rtol=0,atol=1e-12):
            raise RuntimeError(f"mock grid differs: {path}")
        X.append(x)
    X=np.asarray(X,float)
    n=len(X)

    mean=X.mean(axis=0)
    C=np.cov(X,rowvar=False,ddof=1)
    data_score=projected_stat(data-mean,C,n,signal,Nmat)

    loo=[]
    for i in range(n):
        train=np.delete(X,i,axis=0)
        mt=train.mean(axis=0)
        Ct=np.cov(train,rowvar=False,ddof=1)
        sc=projected_stat(X[i]-mt,Ct,len(train),signal,Nmat)
        loo.append({
            "index":i,
            "path":paths[i],
            **sc,
        })

    absz=np.asarray([r["abs_z"] for r in loo],float)
    sh=np.asarray([r["minus2_log_lambda_sellentin_heavens"] for r in loo],float)
    nge=int(np.sum(absz>=data_score["abs_z"]))
    nge_sh=int(np.sum(sh>=data_score["minus2_log_lambda_sellentin_heavens"]))
    pz=float((nge+1)/(n+1))
    psh=float((nge_sh+1)/(n+1))

    result={
        "scope":"LOO empirical calibration of frozen z-resolved wake matched filter",
        "n_mocks":n,
        "vector_dimension":int(len(data)),
        "nuisance_columns":names,
        "data_score":data_score,
        "empirical_abs_z":{
            "n_ge_data":nge,
            "p_plus_one":pz,
            "two_sided_gaussian_equivalent_sigma":float(norm.isf(pz/2.0)),
            "loo_median_abs_z":float(np.median(absz)),
            "loo_max_abs_z":float(np.max(absz)),
        },
        "empirical_sellentin_heavens_lr":{
            "n_ge_data":nge_sh,
            "p_plus_one":psh,
            "two_sided_gaussian_equivalent_sigma":float(norm.isf(psh/2.0)),
            "loo_median_minus2loglambda":float(np.median(sh)),
            "loo_max_minus2loglambda":float(np.max(sh)),
        },
        "resolution_note":(
            f"With {n} mocks the minimum p+1 is 1/{n+1}={1/(n+1):.6g}; "
            "the empirical calibration therefore cannot resolve a smaller tail probability."
        ),
        "analysis_scope":(
            "Each held-out score uses N-1 mocks, while the data score uses all N mocks. "
            "The empirical tail probability uses the finite-sample +1 correction."
        ),
        "loo_scores":loo,
    }
    Path(args.out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
