#!/usr/bin/env python3
"""Redshift-resolved DESI DR1 BGS odd-sector deployment for Phase 7.

This is the second-stage real-survey validation. It keeps the public DESI DR1
clustering catalogs, DESI weights, survey randoms and the redshift-matched
50/50 luminosity split from the pilot, but removes two limitations exposed by
the first run:

1. the odd data vector is retained in broad redshift slices instead of being
   collapsed to one effective redshift, so the distinct z evolution of the
   hidden-wake and relativistic/Doppler templates can contribute to the fit;
2. the covariance uses many more delete-one sky regions than data-vector
   elements, and the permutation test uses an empirical Mahalanobis tail
   probability rather than a chi-square approximation from a tiny null set.

The wake coefficient remains a matched-filter coefficient, not an absolute
DESI wake constraint. A publication-level absolute constraint still requires
window-convolved theory and split-specific mock covariance/validation.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import chi2

import desi_dr1_phase7_lss as p


def prepare(cat, idx):
    return {
        'x':p.distance_xyz(cat['ra'][idx],cat['dec'][idx],cat['z'][idx]),
        'w':np.asarray(cat['w'][idx],float),
        'z':np.asarray(cat['z'][idx],float),
        'ra':np.asarray(cat['ra'][idx],float),
        'region':np.asarray(cat['region'][idx]),
    }


def estimator_prepared(D,labD,R,labR,edges,theta_min_deg,maskD=None,maskR=None):
    if maskD is None: maskD=np.ones(len(labD),bool)
    if maskR is None: maskR=np.ones(len(labR),bool)
    xd=D['x'][maskD]; wd=D['w'][maskD]; dl=labD[maskD]
    xr=R['x'][maskR]; wr=R['w'][maskR]; rl=labR[maskR]
    DH,DL=dl==1,dl==0; RH,RL=rl==1,rl==0
    if min(DH.sum(),DL.sum(),RH.sum(),RL.sum())<20:
        raise RuntimeError('too few objects in redshift/cap estimator subset')
    def pc(xa,wa,xb,wb): return p.pair_moments(xa,wa,xb,wb,edges,theta_min_deg)
    dd0,dd1=pc(xd[DH],wd[DH],xd[DL],wd[DL])
    dr0,dr1=pc(xd[DH],wd[DH],xr[RL],wr[RL])
    rd0,rd1=pc(xr[RH],wr[RH],xd[DL],wd[DL])
    rr0,rr1=pc(xr[RH],wr[RH],xr[RL],wr[RL])
    WDH,WDL,WRH,WRL=wd[DH].sum(),wd[DL].sum(),wr[RH].sum(),wr[RL].sum()
    dd0/=WDH*WDL; dd1/=WDH*WDL; dr0/=WDH*WRL; dr1/=WDH*WRL
    rd0/=WRH*WDL; rd1/=WRH*WDL; rr0/=WRH*WRL; rr1/=WRH*WRL
    den=np.maximum(rr0,1e-300)
    return (dd0-dr0-rd0+rr0)/den,(dd1-dr1-rd1+rr1)/den,{
        'N_bright':int(DH.sum()),'N_faint':int(DL.sum()),
        'R_bright':int(RH.sum()),'R_faint':int(RL.sum())}


def z_vector(D,labD,R,labR,edges,zedges,theta_min_deg,baseD=None,baseR=None):
    if baseD is None: baseD=np.ones(len(labD),bool)
    if baseR is None: baseR=np.ones(len(labR),bool)
    xi0=[]; xi1=[]; meta=[]
    for lo,hi in zip(zedges[:-1],zedges[1:]):
        md=baseD&(D['z']>=lo)&(D['z']<hi)
        mr=baseR&(R['z']>=lo)&(R['z']<hi)
        a,b,c=estimator_prepared(D,labD,R,labR,edges,theta_min_deg,md,mr)
        ze=float(np.average(D['z'][md],weights=D['w'][md]))
        xi0.extend(a.tolist()); xi1.extend(b.tolist())
        meta.append({'zlo':float(lo),'zhi':float(hi),'z_effective':ze,'counts':c})
    return np.asarray(xi0),np.asarray(xi1),meta


def template_vector(path,s,zmeta):
    wake=[]; dop=[]
    for m in zmeta:
        w,d=p.load_templates(path,s,float(m['z_effective']))
        wake.extend(w.tolist()); dop.extend(d.tolist())
    return np.asarray(wake),np.asarray(dop)


def metric_cosine(a,b,cinv):
    aa=float(a@cinv@a); bb=float(b@cinv@b); ab=float(a@cinv@b)
    return ab/np.sqrt(max(aa*bb,1e-300))


def block_template(v,iz,nz,nb):
    out=np.zeros(nz*nb,float); out[iz*nb:(iz+1)*nb]=v
    return out


def conservative_templates(s,zmeta,dop,wake):
    nz=len(zmeta); nb=len(s); ts=[]; names=[]
    b1=(s/60.0)**-1; b2=(s/60.0)**-2
    for iz in range(nz):
        ds=dop[iz*nb:(iz+1)*nb]
        ts.extend([block_template(ds,iz,nz,nb),block_template(b1,iz,nz,nb),block_template(b2,iz,nz,nb)])
        names.extend([f'doppler_z{iz}',f'broadband_inv_s_z{iz}',f'broadband_inv_s2_z{iz}'])
    ts.append(wake); names.append('wake_matched_filter')
    return ts,names


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',nargs='+',required=True); ap.add_argument('--random',nargs='+',required=True)
    ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40)
    ap.add_argument('--dz-split',type=float,default=0.02)
    ap.add_argument('--analysis-z-edges',default='0.10,0.20,0.30,0.40')
    ap.add_argument('--max-per-class',type=int,default=12000); ap.add_argument('--max-random-per-class',type=int,default=24000)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--jackknife',type=int,default=30)
    ap.add_argument('--permutations',type=int,default=48); ap.add_argument('--seed',type=int,default=20260913)
    ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    edges=np.array([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:])
    zedges=np.array([float(x) for x in args.analysis_z_edges.split(',')])
    if zedges[0]<args.zmin-1e-12 or zedges[-1]>args.zmax+1e-12: raise ValueError('analysis z edges outside selected range')

    catD=p.read_catalog(args.data,'data'); catR=p.read_catalog(args.random,'random')
    iD,lD,sD,thresholds=p.make_data_split(catD,args.zmin,args.zmax,args.dz_split)
    iD,lD,sD=p.stratified_subsample(iD,lD,sD,args.max_per_class,rng)
    iR,lR,sR=p.split_randoms(catR,args.zmin,args.zmax,args.dz_split,thresholds,args.max_random_per_class,rng)
    D=prepare(catD,iD); R=prepare(catR,iR)

    xi0,xi1,zmeta=z_vector(D,lD,R,lR,edges,zedges,args.theta_min_deg)
    pdim=len(xi1); nz=len(zmeta); nb=len(s)

    regD,nD=p.sky_jackknife_regions(D['ra'],D['region'],args.jackknife)
    regR,nR=p.sky_jackknife_regions(R['ra'],R['region'],args.jackknife)
    nj=min(nD,nR)
    if nj<=pdim+2:
        raise RuntimeError(f'jackknife regions ({nj}) must exceed vector dimension ({pdim}) by >2')
    jk=[]
    for j in range(nj):
        _,v,_=z_vector(D,lD,R,lR,edges,zedges,args.theta_min_deg,regD!=j,regR!=j)
        jk.append(v); print('JK_Z',j,json.dumps(v.tolist()))
    jk=np.asarray(jk); jm=jk.mean(axis=0)
    cov_raw=(nj-1)/nj*((jk-jm).T@(jk-jm))
    ev=np.linalg.eigvalsh(cov_raw); maxev=float(max(ev[-1],1e-300))
    ridge=max(maxev*1e-7,float(np.median(np.diag(cov_raw)))*1e-6,1e-14)
    cov=cov_raw+np.eye(pdim)*ridge
    rank_raw=int(np.linalg.matrix_rank(cov_raw,tol=maxev*1e-10))
    cond=float(np.linalg.cond(cov))
    cinv=np.linalg.pinv(cov,rcond=1e-10)

    null=[]
    for ip in range(args.permutations):
        lp=lD.copy()
        for st in np.unique(sD):
            q=np.where(sD==st)[0]; lp[q]=rng.permutation(lp[q])
        _,v,_=z_vector(D,lp,R,lR,edges,zedges,args.theta_min_deg)
        null.append(v); print('PERM_Z',ip,json.dumps(v.tolist()))
    null=np.asarray(null); null_mean=null.mean(axis=0); y=xi1-null_mean
    tdata=float(y@cinv@y)
    tperm=np.array([float((v-null_mean)@cinv@(v-null_mean)) for v in null])
    pemp=float((1+np.sum(tperm>=tdata))/(len(tperm)+1))

    wake,dop=template_vector(args.templates,s,zmeta)
    raw_fit=p.gls_fit(xi1,cov,[dop,wake],['doppler','wake_matched_filter'])
    fit2=p.gls_fit(y,cov,[dop,wake],['doppler','wake_matched_filter'])
    ctemps,cnames=conservative_templates(s,zmeta,dop,wake)
    fitc=p.gls_fit(y,cov,ctemps,cnames)
    td_corr=float(metric_cosine(wake,dop,cinv))

    perm_wake=[]
    for v in null:
        f=p.gls_fit(v-null_mean,cov,[dop,wake],['doppler','wake_matched_filter'])
        perm_wake.append(float(f['wake_matched_filter']['amplitude']))
    perm_wake=np.asarray(perm_wake)
    aw=float(fit2['wake_matched_filter']['amplitude'])
    wake_emp_two=float((1+np.sum(np.abs(perm_wake)>=abs(aw)))/(len(perm_wake)+1))

    caps={}
    for cap in ['NGC','SGC']:
        md=D['region']==cap; mr=R['region']==cap
        try:
            a,b,m=z_vector(D,lD,R,lR,edges,zedges,args.theta_min_deg,md,mr)
            caps[cap]={'xi1_odd':b.tolist(),'xi0_cross':a.tolist(),'z_bins':m}
        except Exception as e:
            caps[cap]={'error':str(e)}

    rows=[]; sig=np.sqrt(np.maximum(np.diag(cov),0))
    k=0
    for iz,m in enumerate(zmeta):
        for si in s:
            rows.append((m['zlo'],m['zhi'],m['z_effective'],si,xi0[k],xi1[k],null_mean[k],y[k],sig[k])); k+=1
    np.savetxt(out/'data_vector_zresolved.csv',np.asarray(rows),delimiter=',',
               header='zlo,zhi,z_effective,s_Mpc_over_h,xi0_cross,xi1_odd,permutation_null_mean,xi1_null_corrected,jackknife_sigma',comments='')
    np.savetxt(out/'jackknife_vectors_zresolved.csv',jk,delimiter=',')
    np.savetxt(out/'jackknife_covariance_zresolved.csv',cov,delimiter=',')
    np.savetxt(out/'permutation_vectors_zresolved.csv',null,delimiter=',')
    np.savetxt(out/'permutation_mahalanobis.csv',tperm,delimiter=',')

    summary={
      'scope':'DESI DR1 BGS redshift-resolved real-survey odd-sector deployment; matched-filter validation, not an absolute DESI wake constraint.',
      'catalogs_data':[Path(x).name for x in args.data],'catalogs_random':[Path(x).name for x in args.random],
      'analysis_z_edges':zedges.tolist(),'z_bins':zmeta,'separation_edges_Mpc_over_h':edges.tolist(),
      'luminosity_split':'weighted median dereddened r magnitude independently in dz=%.3g bins and survey cap'%args.dz_split,
      'selected_counts':{'data':int(len(lD)),'random':int(len(lR)),'data_bright':int((lD==1).sum()),'data_faint':int((lD==0).sum())},
      'theta_min_deg':args.theta_min_deg,'seed':args.seed,'vector_dimension':pdim,
      'jackknife':{'regions':int(nj),'raw_covariance_rank':rank_raw,'ridge_added':ridge,'regularized_condition_number':cond},
      'permutation_control':{'count':int(len(null)),'mahalanobis_data':tdata,'empirical_pvalue':pemp,'wake_amplitude_empirical_two_sided_pvalue':wake_emp_two},
      'template_metric_wake_doppler_cosine':td_corr,
      'fit_raw_doppler_plus_wake':raw_fit,
      'fit_null_corrected_doppler_plus_wake':fit2,
      'fit_null_corrected_conservative_per_z_odd_plus_wake':fitc,
      'cap_controls':caps,
      'xi1_odd':xi1.tolist(),'permutation_null_mean':null_mean.tolist(),'xi1_null_corrected':y.tolist(),
      'split_thresholds':thresholds,
      'wake_amplitude_scope':'One matched-filter coefficient for the production hidden-state wake shape with preserved relative redshift evolution. It is not yet an absolute DESI wake amplitude.',
      'required_upgrade_for_absolute_constraint':'Window-convolved theory template plus split-specific mock covariance and end-to-end validation, preferably Abacus for the luminosity split and EZmock/window controls.'
    }
    (out/'summary_zresolved.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_ZRESOLVED',json.dumps(summary,indent=2))

if __name__=='__main__': main()
