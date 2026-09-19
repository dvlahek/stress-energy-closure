#!/usr/bin/env python3
"""Matched-filter tests for the exact 18D z-resolved LRG x ELG dipole."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm


def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def hartlap(n,p):
    return float((n-p-2)/(n-1)) if n>p+2 else None


def precision(C,n):
    a=hartlap(n,C.shape[0])
    if a is None:
        raise RuntimeError("Not enough mocks for Hartlap-corrected precision")
    return a*np.linalg.pinv(C,rcond=1e-12),a


def gls(y,C,n,signal,nuisance,names):
    P,a=precision(C,n)
    signal=np.asarray(signal,float)
    nuisance=np.asarray(nuisance,float)
    if nuisance.size==0: nuisance=np.empty((len(y),0))
    if nuisance.ndim==1: nuisance=nuisance[:,None]

    def fitX(X):
        if X.shape[1]==0:
            return np.empty(0),np.empty((0,0)),float(y@P@y)
        F=X.T@P@X
        Fi=np.linalg.pinv(F,rcond=1e-12)
        b=Fi@(X.T@P@y)
        r=y-X@b
        return b,Fi,float(r@P@r)

    _,_,q0=fitX(nuisance)
    X=np.column_stack([nuisance,signal])
    b,V,q1=fitX(X)

    if nuisance.shape[1]:
        proj=np.linalg.pinv(nuisance.T@P@nuisance,rcond=1e-12)@(nuisance.T@P@signal)
        sperp=signal-nuisance@proj
    else:
        sperp=signal.copy()
    den=float(sperp@P@sperp)
    z=float((sperp@P@y)/np.sqrt(den)) if den>0 else None

    sn2=float(signal@P@signal)
    retain=float(np.sqrt(max(den,0)/sn2)) if sn2>0 else None
    cos={}
    for i,name in enumerate(names):
        v=nuisance[:,i]
        vv=float(v@P@v)
        cos[name]=float((signal@P@v)/np.sqrt(sn2*vv)) if sn2>0 and vv>0 else None

    params={}
    for i,name in enumerate(names+["wake_shape"]):
        e=float(np.sqrt(max(V[i,i],0)))
        params[name]={"amplitude":float(b[i]),"sigma":e,"z_signed":float(b[i]/e) if e>0 else None}
    dq=max(0.0,float(q0-q1))
    pt=float(chi2.sf(dq,1))
    return {
        "parameters":params,
        "chi2_nuisance_only":q0,
        "chi2_nuisance_plus_signal":q1,
        "delta_chi2_signal":dq,
        "pvalue_signal_two_sided":pt,
        "gaussian_equivalent_sigma_two_sided":float(norm.isf(pt/2)) if pt>0 else float("inf"),
        "matched_filter_z_signed":z,
        "template_geometry":{
            "wake_vs_nuisance_metric_cosines":cos,
            "wake_metric_norm_retained_after_nuisance_projection":retain,
        },
        "hartlap_factor":a,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measurement",required=True)
    ap.add_argument("--covariance",required=True)
    ap.add_argument("--template",required=True)
    ap.add_argument("--nuisance-cols",default="")
    ap.add_argument("--subtract-mock-mean",action="store_true")
    ap.add_argument("--outdir",required=True)
    args=ap.parse_args()

    m=np.genfromtxt(args.measurement,delimiter=",",names=True)
    t=np.genfromtxt(args.template,delimiter=",",names=True)
    B=np.load(args.covariance,allow_pickle=False)

    keys_m=np.column_stack([m["zlo"],m["zhi"],m["s_Mpc_over_h"]]).astype(float)
    keys_t=np.column_stack([t["zlo"],t["zhi"],t["s_Mpc_over_h"]]).astype(float)
    keys_c=np.column_stack([B["zlo"],B["zhi"],B["separation"]]).astype(float)
    if not np.allclose(keys_m,keys_t,rtol=0,atol=1e-12) or not np.allclose(keys_m,keys_c,rtol=0,atol=1e-12):
        raise RuntimeError("measurement/template/covariance z-s grids differ")

    y=np.asarray(m["xi1_LRG_to_ELG"],float)
    mean=np.asarray(B["mean_xi1"],float)
    if args.subtract_mock_mean: y=y-mean
    C=np.asarray(B["cov_xi1"],float)
    n=int(np.asarray(B["n_mocks"]).item())
    wake=np.asarray(t["wake_shape"],float)
    names=[x.strip() for x in args.nuisance_cols.split(",") if x.strip()]
    N=np.column_stack([np.asarray(t[x],float) for x in names]) if names else np.empty((len(y),0))

    P,a=precision(C,n)
    q=float(y@P@y)
    p=float(chi2.sf(q,len(y)))
    result={
        "scope":"Exact z-resolved DESI DR1 LRGxELG dipole matched-filter test",
        "n_mocks":n,
        "vector_dimension":int(len(y)),
        "mock_mean_subtracted":bool(args.subtract_mock_mean),
        "hartlap_factor":a,
        "null_test":{
            "chi2":q,
            "dof":int(len(y)),
            "pvalue":p,
            "gaussian_equivalent_signed_upper_tail_z":float(norm.isf(p)) if 0<p<1 else None,
        },
        "nuisance_columns":names,
        "template_sha256":sha(args.template),
        "fit":gls(y,C,n,wake,N,names),
        "guardrail":(
            "The 0.80-0.90, 0.90-1.00 and 1.00-1.10 bins are frozen a priori. "
            "Do not alter redshift or separation cuts after inspecting this result. "
            "The standard odd basis remains diagnostic until tracer coefficients and survey window are physically linked."
        ),
    }
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    (out/"matched_filter_zresolved_summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
