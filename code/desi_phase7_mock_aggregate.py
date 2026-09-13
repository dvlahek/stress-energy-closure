#!/usr/bin/env python3
"""Aggregate cut-sky Phase-7 mock vectors into a survey-calibrated likelihood.

Primary covariance is Oracle Approximating Shrinkage (OAS), appropriate for the
small Abacus ensemble relative to the 18-component data vector.  The raw sample
covariance and its Hartlap factor are reported as a conservative diagnostic.
The real-data vector is the predeclared full-sample five-tracer null-corrected
vector stored in source_data.
"""
from __future__ import annotations
import argparse, glob, json
from pathlib import Path
import numpy as np

import desi_dr1_phase7_lss as p
import desi_dr1_phase7_zresolved as zr


def oas_cov(X):
    X=np.asarray(X,float); X=X-X.mean(axis=0,keepdims=True)
    n,pdim=X.shape
    S=(X.T@X)/float(n)
    mu=float(np.trace(S))/pdim
    alpha=float(np.mean(S*S))
    den=(n+1.0)*(alpha-mu*mu/pdim)
    shrink=1.0 if den<=0 else min((alpha+mu*mu)/den,1.0)
    C=(1.0-shrink)*S+shrink*mu*np.eye(pdim)
    return C,float(shrink)


def fit_bundle(y,cov,s,zmeta,wake,dop):
    fit2=p.gls_fit(y,cov,[dop,wake],['doppler','wake_proxy'])
    ts,names=zr.conservative_templates(s,zmeta,dop,wake)
    names=[('wake_proxy' if n=='wake_matched_filter' else n) for n in names]
    fitc=p.gls_fit(y,cov,ts,names)
    return fit2,fitc,ts,names


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mock-root',required=True); ap.add_argument('--real-vector',required=True); ap.add_argument('--outdir',required=True)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    vf=sorted(glob.glob(str(Path(args.mock_root)/'**/mock_*_vector.csv'),recursive=True))
    wf=sorted(glob.glob(str(Path(args.mock_root)/'**/mock_*_window.csv'),recursive=True))
    if len(vf)<20 or len(vf)!=len(wf): raise RuntimeError(f'need >=20 matched mocks; vectors={len(vf)} windows={len(wf)}')
    X=[]; W=[]; D=[]; meta_ref=None; s_ref=None
    for vpath,wpath in zip(vf,wf):
        v=np.genfromtxt(vpath,delimiter=',',names=True)
        w=np.genfromtxt(wpath,delimiter=',',names=True)
        X.append(np.asarray(v['xi1_proxy_odd'],float)); W.append(np.asarray(w['wake_forward_response'],float)); D.append(np.asarray(w['doppler_forward_response'],float))
        s=np.asarray(v['s_Mpc_over_h'],float)
        zlo=np.asarray(v['zlo'],float); zhi=np.asarray(v['zhi'],float); ze=np.asarray(v['z_effective'],float)
        if s_ref is None:
            s_ref=np.unique(s)
            meta_ref=[]
            for lo,hi in zip(np.unique(zlo),np.unique(zhi)):
                q=(np.isclose(zlo,lo)&np.isclose(zhi,hi)); meta_ref.append({'zlo':float(lo),'zhi':float(hi),'z_effective':float(np.mean(ze[q]))})
    X=np.asarray(X); W=np.asarray(W); D=np.asarray(D); n,pdim=X.shape
    mean=X.mean(axis=0); Csample=np.cov(X,rowvar=False,ddof=1); Coas,shrink=oas_cov(X)
    rank=int(np.linalg.matrix_rank(Csample)); cond=float(np.linalg.cond(Coas))
    hartlap=float((n-pdim-2)/(n-1)) if n>pdim+2 else 0.0
    Wm=W.mean(axis=0); Dm=D.mean(axis=0)

    r=np.genfromtxt(args.real_vector,delimiter=',',names=True)
    yreal=np.asarray(r['xi1_null_corrected'],float)
    if len(yreal)!=pdim: raise RuntimeError('real/mock vector dimension mismatch')
    # Mock mean absorbs residual survey-window/estimator bias not removed by the real-data permutation control.
    y=yreal-mean
    fit2,fitc,ts,names=fit_bundle(y,Coas,s_ref,meta_ref,Wm,Dm)
    if hartlap>0:
        Ch=Csample/hartlap
        fit2_h,fitc_h,_,_=fit_bundle(y,Ch,s_ref,meta_ref,Wm,Dm)
    else:
        fit2_h=fitc_h=None

    # Null distribution of the conservative wake coefficient.
    aw=[]; sig=[]
    for row in X:
        f=p.gls_fit(row-mean,Coas,ts,names)
        aw.append(float(f['wake_proxy']['amplitude'])); sig.append(float(f['wake_proxy']['sigma']))
    aw=np.asarray(aw); sig=np.asarray(sig)
    areal=float(fitc['wake_proxy']['amplitude'])
    pemp=float((1+np.sum(np.abs(aw)>=abs(areal)))/(n+1))

    # End-to-end unit injection using each realization's own forward window, fit with the mean survey template.
    rec=[]
    for row,wi in zip(X,W):
        f=p.gls_fit((row-mean)+wi,Coas,ts,names)
        rec.append(float(f['wake_proxy']['amplitude']))
    rec=np.asarray(rec); sigma_ref=float(fitc['wake_proxy']['sigma'])
    coverage=float(np.mean(np.abs(rec-1.0)<=sigma_ref))

    np.savetxt(out/'abacus_mock_covariance_sample.csv',Csample,delimiter=',')
    np.savetxt(out/'abacus_mock_covariance_oas.csv',Coas,delimiter=',')
    np.savetxt(out/'abacus_mock_vectors.csv',X,delimiter=',')
    np.savetxt(out/'abacus_window_templates.csv',np.column_stack([Wm,Dm,W.std(axis=0,ddof=1),D.std(axis=0,ddof=1)]),delimiter=',',header='wake_mean,doppler_mean,wake_std,doppler_std',comments='')
    np.savetxt(out/'abacus_wake_null_and_injection.csv',np.column_stack([aw,rec]),delimiter=',',header='null_wake_amplitude,unit_injection_recovered_amplitude',comments='')
    summary={
      'scope':'AbacusSummit cut-sky covariance and pair-window calibration for the five-tracer luminosity-proxy estimator.',
      'mock_count':int(n),'vector_dimension':int(pdim),'sample_covariance_rank':rank,'oas_shrinkage':shrink,'oas_condition_number':cond,'hartlap_factor_for_raw_sample_covariance':hartlap,
      'mock_mean_vector':mean.tolist(),'mean_wake_forward_response':Wm.tolist(),'mean_doppler_forward_response':Dm.tolist(),
      'real_data_fit_oas':{'minimal':fit2,'conservative':fitc,'empirical_two_sided_wake_pvalue':pemp},
      'real_data_fit_hartlap_sample_covariance':{'minimal':fit2_h,'conservative':fitc_h},
      'mock_null_wake_amplitudes':aw.tolist(),
      'unit_injection_recovery':{'mean':float(rec.mean()),'std':float(rec.std(ddof=1)),'median':float(np.median(rec)),'nominal_one_sigma_coverage_fraction':coverage,'fit_sigma_reference':sigma_ref},
      'interpretation_guardrail':'The amplitude is calibrated to the predeclared luminosity-rank proxy pair model. A physical halo-mass wake amplitude requires a separate proxy-to-halo calibration.'
    }
    (out/'summary_abacus_window_covariance.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_MOCK_AGGREGATE',json.dumps(summary,indent=2))

if __name__=='__main__': main()
