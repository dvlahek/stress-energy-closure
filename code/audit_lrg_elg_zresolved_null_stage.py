#!/usr/bin/env python3
"""Stage-1 audit: inspect covariance and nuisance-only fit before wake significance."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import chi2

def hartlap(n,p):
    return float((n-p-2)/(n-1)) if n>p+2 else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measurement",required=True)
    ap.add_argument("--covariance",required=True)
    ap.add_argument("--nuisance-template",required=True)
    ap.add_argument("--nuisance-col",default="standard_linked_shape")
    ap.add_argument("--subtract-mock-mean",action="store_true")
    ap.add_argument("--out",required=True)
    args=ap.parse_args()

    m=np.genfromtxt(args.measurement,delimiter=",",names=True)
    t=np.genfromtxt(args.nuisance_template,delimiter=",",names=True)
    B=np.load(args.covariance,allow_pickle=False)
    km=np.column_stack([m["zlo"],m["zhi"],m["s_Mpc_over_h"]]).astype(float)
    kt=np.column_stack([t["zlo"],t["zhi"],t["s_Mpc_over_h"]]).astype(float)
    kc=np.column_stack([B["zlo"],B["zhi"],B["separation"]]).astype(float)
    if not np.allclose(km,kt,rtol=0,atol=1e-12) or not np.allclose(km,kc,rtol=0,atol=1e-12):
        raise RuntimeError("measurement/template/covariance grids differ")

    y=np.asarray(m["xi1_LRG_to_ELG"],float)
    if args.subtract_mock_mean:
        y=y-np.asarray(B["mean_xi1"],float)
    C=np.asarray(B["cov_xi1"],float)
    n=int(np.asarray(B["n_mocks"]).item()); p=len(y)
    a=hartlap(n,p)
    if a is None: raise RuntimeError("not enough mocks for Hartlap precision")
    P=a*np.linalg.pinv(C,rcond=1e-12)
    qzero=float(y@P@y)
    v=np.asarray(t[args.nuisance_col],float)
    den=float(v@P@v)
    amp=float((v@P@y)/den)
    sig=float(1/np.sqrt(den))
    r=y-amp*v
    qn=float(r@P@r)
    eig=np.linalg.eigvalsh(C)
    out={
      "scope":"Stage-1 covariance and nuisance-only audit; inspect before wake matched filter",
      "n_mocks":n,"dimension":p,"hartlap_factor":a,
      "covariance":{"condition_number":float(np.linalg.cond(C)),"min_eigenvalue":float(eig[0]),"max_eigenvalue":float(eig[-1]),"positive_definite":bool(np.all(eig>0))},
      "zero_null":{"chi2":qzero,"dof":p,"pvalue":float(chi2.sf(qzero,p))},
      "nuisance_only":{"amplitude":amp,"sigma":sig,"z_signed":amp/sig,"chi2":qn,"dof_residual":p-1,"pvalue_residual":float(chi2.sf(qn,p-1))},
      "mock_mean_subtracted":bool(args.subtract_mock_mean),
      "analysis_scope":"We evaluate the zero-null and nuisance-only residuals using the covariance adopted for the matched-filter analysis."
    }
    Path(args.out).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))

if __name__=="__main__":
    main()
