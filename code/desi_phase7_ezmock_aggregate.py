#!/usr/bin/env python3
"""Aggregate the Phase-7 EZmock random-rank placebo covariance ensemble."""
from __future__ import annotations
import argparse, glob, json, re
from pathlib import Path
import numpy as np

import desi_dr1_phase7_lss as p
from desi_phase7_mock_aggregate import oas_cov, fit_bundle


def mock_id(path):
    m=re.search(r'mock_(\d+)_vector\.csv$',str(path))
    if not m: raise ValueError(f'cannot parse mock id from {path}')
    return int(m.group(1))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mock-root',required=True)
    ap.add_argument('--real-vector',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--min-mock-id',type=int,default=3,
                    help='minimum mock id included; default 3 excludes legacy smoke mocks 1-2')
    ap.add_argument('--require-shared-random',action='store_true',
                    help='require summary provenance for shared-angular/per-realization-z random mode')
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    all_v=glob.glob(str(Path(args.mock_root)/'**/mock_*_vector.csv'),recursive=True)
    vf=sorted([x for x in all_v if mock_id(x)>=args.min_mock_id],key=mock_id)
    pairs=[]
    for vpath in vf:
        mid=mock_id(vpath); base=Path(vpath).parent
        wpath=base/f'mock_{mid:02d}_window.csv'; spath=base/f'mock_{mid:02d}_summary.json'
        if not wpath.is_file():
            raise RuntimeError(f'missing window for mock {mid}: {wpath}')
        if args.require_shared_random:
            if not spath.is_file(): raise RuntimeError(f'missing summary for mock {mid}: {spath}')
            sm=json.loads(spath.read_text())
            mode=str(sm.get('random_geometry_mode',''))
            if 'shared angular cache' not in mode:
                raise RuntimeError(f'mock {mid} is not shared-random homogeneous: {mode!r}')
            if int(sm.get('random_count',0)) < int(1.8*sm.get('data_count',0)):
                raise RuntimeError(f'mock {mid} does not have approximately 2x random density')
        pairs.append((vpath,str(wpath)))

    if len(pairs)<40:
        raise RuntimeError(f'need >=40 homogeneous EZmocks with id>={args.min_mock_id}; found={len(pairs)}')

    X=[]; W=[]; D=[]; s_ref=None; meta_ref=None; used_ids=[]
    for vpath,wpath in pairs:
        mid=mock_id(vpath); used_ids.append(mid)
        v=np.genfromtxt(vpath,delimiter=',',names=True); w=np.genfromtxt(wpath,delimiter=',',names=True)
        X.append(np.asarray(v['xi1_proxy_odd'],float)); W.append(np.asarray(w['wake_forward_response'],float)); D.append(np.asarray(w['doppler_forward_response'],float))
        if s_ref is None:
            s=np.asarray(v['s_Mpc_over_h'],float); s_ref=np.unique(s); zlo=np.asarray(v['zlo']); zhi=np.asarray(v['zhi']); ze=np.asarray(v['z_effective'])
            meta_ref=[]
            for lo,hi in zip(np.unique(zlo),np.unique(zhi)):
                q=np.isclose(zlo,lo)&np.isclose(zhi,hi); meta_ref.append({'zlo':float(lo),'zhi':float(hi),'z_effective':float(np.mean(ze[q]))})
    X=np.asarray(X); W=np.asarray(W); D=np.asarray(D); n,pdim=X.shape
    mean=X.mean(axis=0); Csample=np.cov(X,rowvar=False,ddof=1); Coas,shrink=oas_cov(X)
    rank=int(np.linalg.matrix_rank(Csample)); cond=float(np.linalg.cond(Coas)); hartlap=float((n-pdim-2)/(n-1)) if n>pdim+2 else 0.0
    Wm=W.mean(axis=0); Dm=D.mean(axis=0)
    r=np.genfromtxt(args.real_vector,delimiter=',',names=True); y=np.asarray(r['xi1_null_corrected'],float)-mean
    fit2,fitc,ts,names=fit_bundle(y,Coas,s_ref,meta_ref,Wm,Dm)
    if hartlap>0:
        fit2_h,fitc_h,_,_=fit_bundle(y,Csample/hartlap,s_ref,meta_ref,Wm,Dm)
    else: fit2_h=fitc_h=None
    aw=[]; rec=[]
    for row,wi in zip(X,W):
        f=p.gls_fit(row-mean,Coas,ts,names); aw.append(float(f['wake_proxy']['amplitude']))
        fi=p.gls_fit((row-mean)+wi,Coas,ts,names); rec.append(float(fi['wake_proxy']['amplitude']))
    aw=np.asarray(aw); rec=np.asarray(rec); areal=float(fitc['wake_proxy']['amplitude'])
    pemp=float((1+np.sum(np.abs(aw)>=abs(areal)))/(n+1)); sigma=float(fitc['wake_proxy']['sigma']); coverage=float(np.mean(np.abs(rec-1.0)<=sigma))
    np.savetxt(out/'ezmock_placebo_covariance_sample.csv',Csample,delimiter=','); np.savetxt(out/'ezmock_placebo_covariance_oas.csv',Coas,delimiter=',')
    np.savetxt(out/'ezmock_placebo_vectors.csv',X,delimiter=','); np.savetxt(out/'ezmock_placebo_window_templates.csv',np.column_stack([Wm,Dm,W.std(axis=0,ddof=1),D.std(axis=0,ddof=1)]),delimiter=',',header='wake_mean,doppler_mean,wake_std,doppler_std',comments='')
    summary={
      'scope':'Custom cut-sky EZmock five-tracer equal-count random-rank placebo covariance/systematics control.',
      'proxy_kind':'RAN_NUM_0_1 equal-count rank within narrow-z and Galactic-cap cells',
      'mock_count':int(n),'mock_ids':used_ids,'minimum_mock_id':int(args.min_mock_id),
      'random_geometry_mode':'shared angular cache with per-realization narrow-z/cap redshift resampling' if args.require_shared_random else 'not enforced',
      'vector_dimension':int(pdim),'sample_covariance_rank':rank,'oas_shrinkage':shrink,'oas_condition_number':cond,'hartlap_factor':hartlap,
      'real_data_sensitivity_fit_oas':{'minimal':fit2,'conservative':fitc,'empirical_two_sided_wake_pvalue_against_placebo_null':pemp},
      'real_data_sensitivity_fit_hartlap_sample_covariance':{'minimal':fit2_h,'conservative':fitc_h},
      'placebo_null_wake_amplitudes':aw.tolist(),
      'unit_injection_recovery':{'mean':float(rec.mean()),'std':float(rec.std(ddof=1)),'median':float(np.median(rec)),'nominal_one_sigma_coverage_fraction':coverage,'fit_sigma_reference':sigma},
      'absolute_likelihood_claim':False,
      'guardrail':'The released DR1 EZmock BGS files used here do not contain R_MAG_APP/R_MAG_ABS. This homogeneous shared-angular ensemble tests geometry, pair compression, covariance conditioning and false-positive behaviour with equal-count random ranks. Random redshifts are redrawn from each mock data realization in the same narrow-z/cap cells. It is not a luminosity-matched covariance. Abacus is the physical luminosity-ranked validation.'
    }
    (out/'summary_ezmock_placebo_covariance.json').write_text(json.dumps(summary,indent=2)+'\n'); print('PHASE7_EZMOCK_PLACEBO_AGGREGATE',json.dumps(summary,indent=2))

if __name__=='__main__': main()
