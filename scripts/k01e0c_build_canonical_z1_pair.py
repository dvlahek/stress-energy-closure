#!/usr/bin/env python3
"""K01E0C: canonical closest z=1 moment-null pair for K01E."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
import class_response_optimize as cro

PR=ROOT/"source_data/k01e0c_canonical_z1_moment_null_pair_prereg_2026-11-06.json"
PR_BLOB="fc34036f53e70022dd75d3f710f1d0d476bd2426"
SRC=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
OUTDIR=ROOT/"k01e0c_z1_pair"
CSV=OUTDIR/"canonical_z1_moment_null_pair_4000q.csv"
REPORT=OUTDIR/"k01e0c_canonical_z1_moment_null_pair_qa.json"

def blob(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def need(ok,msg):
    if not ok: raise RuntimeError("K01E0C_STOP: "+msg)

def qweights(q):
    w=np.empty_like(q)
    w[0]=0.5*(q[1]-q[0]); w[-1]=0.5*(q[-1]-q[-2])
    w[1:-1]=0.5*(q[2:]-q[:-2])
    return w

def moment_weights(q,z,mass=.06):
    a=1.0/(1.0+z)
    tnu0=cro.T_NCDM*cro.TCMB_K*cro.KB_EV_K
    y=a*mass/tnu0
    eps=np.sqrt(q*q+y*y)
    return np.vstack([q*q,q*q*eps,q**4/eps]), y

def moments(f,q,Q,W):
    return np.array([np.sum(Q*w*f) for w in W],dtype=float)

def relsym(a,b):
    return np.abs(a-b)/np.maximum(0.5*(np.abs(a)+np.abs(b)),1e-300)

def main():
    p=json.loads(PR.read_text())
    need(blob(PR)==PR_BLOB,"prereg blob")
    need(p["stage"]=="K01E0C_CANONICAL_CLOSEST_Z1_MOMENT_NULL_PAIR_PREREG","stage")

    x=np.loadtxt(SRC,delimiter=",",skiprows=1)
    need(x.shape==(4000,4),"source shape")
    q,f0,fp0,fm0=x.T
    need(np.all(np.diff(q)>0) and q[0]==0 and q[-1]==20,"q grid")
    need(np.all(f0>0) and np.all(fp0>0) and np.all(fm0>0),"original positivity")

    Q=qweights(q)
    sqrtQ=np.sqrt(Q)
    r0=(fp0-fm0)/(0.6*f0)
    need(abs(np.max(np.abs(r0))-1)<5e-12,"original fractional direction cap")

    W1,y1=moment_weights(q,1.0)
    # Constraint A r=0 with A_ij=Q_j W_i,j F0_j.
    # In u=sqrt(Q)r coordinates, row matrix B=A Q^{-1/2}.
    B=W1*f0[None,:]*sqrtQ[None,:]

    # Normalize rows for conditioning; row space is unchanged.
    rn=np.linalg.norm(B,axis=1)
    need(np.all(rn>0),"constraint row norm")
    Bs=B/rn[:,None]

    U,S,Vh=np.linalg.svd(Bs,full_matrices=False)
    rank=int(np.sum(S > max(Bs.shape)*S[0]*np.finfo(float).eps))
    need(rank==3,f"constraint rank {rank}")

    u0=sqrtQ*r0
    Vr=Vh[:rank,:]
    u1=u0 - Vr.T@(Vr@u0)
    r1=u1/sqrtQ

    max_before=float(np.max(np.abs(r1)))
    need(np.isfinite(max_before) and max_before>0,"projected direction")
    r1=r1/max_before

    frac=float(p["immutable_inputs"]["pointwise_fractional_cap"])
    fp=f0*(1.0+frac*r1)
    fm=f0*(1.0-frac*r1)

    need(np.all(fp>0) and np.all(fm>0),"projected positivity")
    maxfrac=max(float(np.max(np.abs(fp-f0)/f0)),float(np.max(np.abs(fm-f0)/f0)))
    sym=float(np.max(np.abs(fp+fm-2*f0)/np.maximum(2*np.abs(f0),1e-300)))

    m0=moments(f0,q,Q,W1)
    mp=moments(fp,q,Q,W1)
    mm=moments(fm,q,Q,W1)
    gaps=relsym(mp,mm)

    dot=float(np.sum(Q*r0*r1))
    n0=float(np.sum(Q*r0*r0))
    n1=float(np.sum(Q*r1*r1))
    cosine=dot/np.sqrt(n0*n1)
    l2change=np.sqrt(float(np.sum(Q*(r1-r0)**2))/n0)

    W095,y095=moment_weights(q,.95)
    mp095=moments(fp,q,Q,W095)
    mm095=moments(fm,q,Q,W095)
    gaps095=relsym(mp095,mm095)

    a=p["acceptance"]
    need(rank==int(a["constraint_rank_exactly"]),"constraint rank acceptance")
    need(abs(maxfrac-frac)<=float(a["max_fractional_distortion_equals_0_30_abs_tolerance"]),"30% cap")
    need(sym<=float(a["pair_symmetry_max_relative_error_at_most"]),"symmetry")
    need(float(np.max(gaps))<=float(a["z1_max_symmetric_relative_n_rho_P_gap_at_most"]),f"z1 moment gap {gaps}")

    OUTDIR.mkdir(exist_ok=True)
    np.savetxt(
        CSV,
        np.column_stack([q,f0,fp,fm,r0,r1]),
        delimiter=",",
        header="q_dimensionless,F0_CLASS_normalized,Fplus_z1_null,Fminus_z1_null,original_E8_fractional_direction,z1_projected_fractional_direction",
        comments="",
        fmt="%.17e",
    )

    report={
      "stage":p["stage"],
      "status":"PASS_CANONICAL_CLOSEST_Z1_MOMENT_NULL_PAIR",
      "source_E8_git_blob":blob(SRC),
      "source_E8_sha256":sha256(SRC),
      "projection":{
        "constraint_singular_values_row_normalized":S.tolist(),
        "constraint_rank":rank,
        "projection_max_abs_before_renormalization":max_before,
        "weighted_fractional_direction_cosine_with_original":cosine,
        "weighted_fractional_direction_relative_L2_change":l2change,
      },
      "z1":{
        "mass_over_Tnu":y1,
        "FD_n_rho_P":m0.tolist(),
        "Fplus_n_rho_P":mp.tolist(),
        "Fminus_n_rho_P":mm.tolist(),
        "Fplus_vs_Fminus_symmetric_relative_gaps":gaps.tolist(),
        "max_gap":float(np.max(gaps)),
      },
      "z0p95_diagnostic":{
        "mass_over_Tnu":y095,
        "Fplus_n_rho_P":mp095.tolist(),
        "Fminus_n_rho_P":mm095.tolist(),
        "Fplus_vs_Fminus_symmetric_relative_gaps":gaps095.tolist(),
        "max_gap":float(np.max(gaps095)),
      },
      "positivity_and_symmetry":{
        "min_Fplus":float(np.min(fp)),
        "min_Fminus":float(np.min(fm)),
        "max_fractional_distortion":maxfrac,
        "max_pair_symmetry_relative_error":sym,
      },
      "output_csv":str(CSV.relative_to(ROOT)),
      "output_csv_sha256":None,
      "interpretation":p["interpretation"],
      "guards":p["guards"],
      "next_stage":p["next_if_pass"]
    }
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    report["output_csv_sha256"]=sha256(CSV)
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")

    print("K01E0C_Z1_MAX_MOMENT_GAP",float(np.max(gaps)),flush=True)
    print("K01E0C_DIRECTION_COSINE",cosine,flush=True)
    print("K01E0C_DIRECTION_REL_L2_CHANGE",l2change,flush=True)
    print("K01E0C_Z095_DIAGNOSTIC_MAX_GAP",float(np.max(gaps095)),flush=True)
    print("K01E0C_CSV_SHA256",report["output_csv_sha256"],flush=True)
    print("K01E0C_PASS_CANONICAL_CLOSEST_Z1_MOMENT_NULL_PAIR",flush=True)

if __name__=="__main__": main()
