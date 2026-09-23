#!/usr/bin/env python3
"""Fit the genuine DESI DR1 LRGxELG dipole with the frozen hidden-state wake shape."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import chi2

def gls(y,cov,templates,names):
    ci=np.linalg.pinv(cov,rcond=1e-10); X=np.column_stack(templates)
    F=X.T@ci@X; Fi=np.linalg.pinv(F,rcond=1e-10)
    b=Fi@(X.T@ci@y); e=np.sqrt(np.maximum(np.diag(Fi),0))
    r=y-X@b; dof=max(len(y)-np.linalg.matrix_rank(X),0)
    o={n:{'amplitude':float(v),'sigma':float(s),'z':float(v/s) if s>0 else None}
       for n,v,s in zip(names,b,e)}
    o.update(chi2=float(r@ci@r),dof=int(dof),pvalue=float(chi2.sf(float(r@ci@r),dof)) if dof>0 else None)
    return o

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--measurement-dir',required=True); ap.add_argument('--template',required=True)
    ap.add_argument('--outdir',required=True)
    a=ap.parse_args(); md=Path(a.measurement_dir); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    d=np.genfromtxt(md/'lrg_elg_odd_multipoles.csv',delimiter=',',names=True)
    C=np.loadtxt(md/'lrg_elg_dipole_cov.csv',delimiter=',')
    t=np.genfromtxt(a.template,delimiter=',',names=True)
    s=np.asarray(d['s_Mpc_over_h']); y=np.asarray(d['xi1_lrg_to_elg'])
    wake=np.interp(s,t['s_Mpc_over_h'],t['wake_shape']); dop=np.interp(s,t['s_Mpc_over_h'],t['doppler_shape'])
    simple=gls(y,C,[dop,wake],['doppler','wake'])
    conservative=gls(y,C,[dop,60/s,(60/s)**2,wake],['doppler','broadband_inv_s','broadband_inv_s2','wake'])
    x3=np.asarray(d['xi3_lrg_to_elg']); C3=np.loadtxt(md/'lrg_elg_octupole_cov.csv',delimiter=',')
    ci3=np.linalg.pinv(C3,rcond=1e-10); q=float(x3@ci3@x3)
    S=dict(scope='genuine DESI DR1 LRGxELG matched-filter shape test; not an absolute wake-amplitude measurement',
           fit_doppler_plus_wake=simple,fit_conservative_odd_plus_wake=conservative,
           octupole_null={'chi2':q,'dof':len(x3),'pvalue':float(chi2.sf(q,len(x3)))},
           guardrail='Inspect the conservative wake coefficient only after pair-seed stability and mock/window closure.')
    (out/'fit_summary.json').write_text(json.dumps(S,indent=2)+'\n'); print(json.dumps(S,indent=2))
if __name__=='__main__': main()
