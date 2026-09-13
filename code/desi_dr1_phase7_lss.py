#!/usr/bin/env python3
"""Nature-Astronomy-grade Phase-7 DESI DR1 odd-sector deployment pilot.

This analysis intentionally uses the public clustering-ready BGS catalogs,
DESI clustering weights and survey randoms instead of the Gfinder-only smoke
test. A 50/50 luminosity-ranked split is made independently in narrow
redshift bins, which keeps the two tracer n(z) distributions closely matched
without requiring an external halo-mass catalog. The oriented bright-faint
cross-correlation dipole is estimated with a weighted Landy-Szalay estimator.

The primary covariance is delete-one sky jackknife. A stratified label-
permutation ensemble is also written as an estimator/null control. The final
fit compares a forecast-derived wake matched-filter shape against a free
standard odd/Doppler nuisance; a more conservative fit adds two smooth odd
broadband modes. This is a real-survey deployment/validation analysis, not a
DESI collaboration measurement and not a detection claim.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree
from scipy.stats import chi2

COSMO=FlatLambdaCDM(H0=67.4,Om0=0.315,Tcmb0=2.7255)
H=0.674


def col(d,*names,default=None):
    avail={n.upper():n for n in d.names}
    for n in names:
        if n.upper() in avail: return np.asarray(d[avail[n.upper()]])
    if default is not None: return default
    raise KeyError(f'missing columns {names}; available={list(d.names)}')


def read_catalog(paths,kind):
    chunks=[]
    for p in paths:
        region='NGC' if 'NGC' in Path(p).name.upper() else 'SGC' if 'SGC' in Path(p).name.upper() else 'UNK'
        with fits.open(p,memmap=True) as hdul:
            d=hdul[1].data
            ra=col(d,'RA').astype(float); dec=col(d,'DEC').astype(float); z=col(d,'Z').astype(float)
            w=col(d,'WEIGHT',default=np.ones(len(d))).astype(float)
            wf=col(d,'WEIGHT_FKP',default=np.ones(len(d))).astype(float)
            if kind=='data':
                fr=col(d,'FLUX_R_DERED','flux_r_dered').astype(float)
            else:
                fr=np.full(len(d),np.nan)
            chunks.append({'ra':ra,'dec':dec,'z':z,'w':w*wf,'fluxr':fr,
                           'region':np.full(len(d),region,dtype='U3')})
    return {k:np.concatenate([x[k] for x in chunks]) for k in chunks[0]}


def distance_xyz(ra,dec,z):
    r=np.asarray(COSMO.comoving_distance(z).value)*H
    rr=np.deg2rad(ra); dd=np.deg2rad(dec); cd=np.cos(dd)
    return np.column_stack((r*cd*np.cos(rr),r*cd*np.sin(rr),r*np.sin(dd)))


def weighted_median(x,w):
    q=np.argsort(x); x=np.asarray(x)[q]; w=np.asarray(w)[q]
    cs=np.cumsum(w); return float(x[np.searchsorted(cs,0.5*cs[-1],side='left')])


def make_data_split(cat,zmin,zmax,dz):
    good=np.isfinite(cat['ra'])&np.isfinite(cat['dec'])&np.isfinite(cat['z'])&np.isfinite(cat['w'])&(cat['w']>0)&np.isfinite(cat['fluxr'])&(cat['fluxr']>0)&(cat['z']>=zmin)&(cat['z']<zmax)
    idx=np.where(good)[0]; z=cat['z'][idx]; w=cat['w'][idx]
    mag=22.5-2.5*np.log10(cat['fluxr'][idx])
    labels=np.full(len(idx),-1,np.int8); strata=np.full(len(idx),-1,np.int32); thresholds={}; ib=0
    edges=np.arange(zmin,zmax+0.5*dz,dz)
    for reg in sorted(set(cat['region'][idx])):
        for lo,hi in zip(edges[:-1],edges[1:]):
            q=np.where((cat['region'][idx]==reg)&(z>=lo)&(z<hi))[0]
            if len(q)<4: continue
            med=weighted_median(mag[q],w[q]); labels[q]=(mag[q]<=med).astype(np.int8); strata[q]=ib
            thresholds[f'{reg}:{lo:.3f}-{hi:.3f}']={'median_rmag':med,'n':int(len(q)),'bright_weight_fraction':float(w[q][labels[q]==1].sum()/w[q].sum())}
            ib+=1
    keep=strata>=0
    return idx[keep],labels[keep],strata[keep],thresholds


def stratified_subsample(indices,labels,strata,max_per_class,rng):
    out=[]
    for lab in [0,1]:
        q=np.where(labels==lab)[0]
        if len(q)<=max_per_class: out.extend(q.tolist()); continue
        uniq=np.unique(strata[q]); alloc=[]
        for s in uniq:
            qq=q[strata[q]==s]; alloc.append((s,qq,len(qq)))
        total=sum(n for _,_,n in alloc); picked=[]
        for s,qq,n in alloc:
            k=max(1,int(round(max_per_class*n/total))); k=min(k,len(qq)); picked.extend(rng.choice(qq,size=k,replace=False).tolist())
        if len(picked)>max_per_class: picked=rng.choice(np.asarray(picked),size=max_per_class,replace=False).tolist()
        out.extend(picked)
    out=np.asarray(out,int)
    return indices[out],labels[out],strata[out]


def split_randoms(cat,zmin,zmax,dz,data_thresholds,max_per_class,rng):
    good=np.isfinite(cat['ra'])&np.isfinite(cat['dec'])&np.isfinite(cat['z'])&np.isfinite(cat['w'])&(cat['w']>0)&(cat['z']>=zmin)&(cat['z']<zmax)
    idx=np.where(good)[0]; labels=np.full(len(idx),-1,np.int8); strata=np.full(len(idx),-1,np.int32)
    edges=np.arange(zmin,zmax+0.5*dz,dz); ib=0
    for reg in sorted(set(cat['region'][idx])):
        for lo,hi in zip(edges[:-1],edges[1:]):
            q=np.where((cat['region'][idx]==reg)&(cat['z'][idx]>=lo)&(cat['z'][idx]<hi))[0]
            if len(q)==0: continue
            key=f'{reg}:{lo:.3f}-{hi:.3f}'; frac=data_thresholds.get(key,{}).get('bright_weight_fraction',0.5)
            k=int(round(frac*len(q))); pick=np.zeros(len(q),np.int8)
            if k>0: pick[rng.choice(len(q),size=min(k,len(q)),replace=False)]=1
            labels[q]=pick; strata[q]=ib; ib+=1
    keep=strata>=0
    idx,labels,strata=idx[keep],labels[keep],strata[keep]
    return stratified_subsample(idx,labels,strata,max_per_class,rng)


def angular_keep(xa,xb,row,col,theta_min_deg):
    if theta_min_deg<=0: return np.ones(len(row),bool)
    ua=xa[row]/np.linalg.norm(xa[row],axis=1)[:,None]
    ub=xb[col]/np.linalg.norm(xb[col],axis=1)[:,None]
    cost=np.sum(ua*ub,axis=1); th=np.rad2deg(np.arccos(np.clip(cost,-1,1)))
    return th>=theta_min_deg


def pair_moments(xa,wa,xb,wb,edges,theta_min_deg=0.05):
    ta=cKDTree(xa); tb=cKDTree(xb)
    coo=ta.sparse_distance_matrix(tb,float(edges[-1]),output_type='coo_matrix')
    nb=len(edges)-1; sw=np.zeros(nb); swmu=np.zeros(nb)
    if coo.nnz==0: return sw,swmu
    row=np.asarray(coo.row); col=np.asarray(coo.col); sep=np.asarray(coo.data)
    keep=(sep>=edges[0])&(sep<edges[-1]); row=row[keep]; col=col[keep]; sep=sep[keep]
    if len(sep)==0: return sw,swmu
    ak=angular_keep(xa,xb,row,col,theta_min_deg); row=row[ak]; col=col[ak]; sep=sep[ak]
    dv=xb[col]-xa[row]; mid=xb[col]+xa[row]
    mu=np.sum(dv*mid,axis=1)/(np.maximum(sep,1e-12)*np.maximum(np.linalg.norm(mid,axis=1),1e-12))
    pw=wa[row]*wb[col]
    b=np.searchsorted(edges,sep,side='right')-1
    for j in range(nb):
        q=b==j
        if q.any(): sw[j]=pw[q].sum(); swmu[j]=np.sum(pw[q]*(3.0*mu[q]))
    return sw,swmu


def estimator(catD,idxD,labD,catR,idxR,labR,edges,theta_min_deg,maskD=None,maskR=None):
    if maskD is None: maskD=np.ones(len(idxD),bool)
    if maskR is None: maskR=np.ones(len(idxR),bool)
    dsel=idxD[maskD]; dl=labD[maskD]; rsel=idxR[maskR]; rl=labR[maskR]
    xd=distance_xyz(catD['ra'][dsel],catD['dec'][dsel],catD['z'][dsel]); wd=catD['w'][dsel]
    xr=distance_xyz(catR['ra'][rsel],catR['dec'][rsel],catR['z'][rsel]); wr=catR['w'][rsel]
    DH,DL=dl==1,dl==0; RH,RL=rl==1,rl==0
    if min(DH.sum(),DL.sum(),RH.sum(),RL.sum())<20: raise RuntimeError('too few objects in estimator subset')
    def pc(xa,wa,xb,wb): return pair_moments(xa,wa,xb,wb,edges,theta_min_deg)
    dd0,dd1=pc(xd[DH],wd[DH],xd[DL],wd[DL]); dr0,dr1=pc(xd[DH],wd[DH],xr[RL],wr[RL])
    rd0,rd1=pc(xr[RH],wr[RH],xd[DL],wd[DL]); rr0,rr1=pc(xr[RH],wr[RH],xr[RL],wr[RL])
    WDH,WDL,WRH,WRL=wd[DH].sum(),wd[DL].sum(),wr[RH].sum(),wr[RL].sum()
    dd0/=WDH*WDL; dd1/=WDH*WDL; dr0/=WDH*WRL; dr1/=WDH*WRL
    rd0/=WRH*WDL; rd1/=WRH*WDL; rr0/=WRH*WRL; rr1/=WRH*WRL
    den=np.maximum(rr0,1e-300)
    xi0=(dd0-dr0-rd0+rr0)/den; xi1=(dd1-dr1-rd1+rr1)/den
    return xi0,xi1,{'N_bright':int(DH.sum()),'N_faint':int(DL.sum()),'R_bright':int(RH.sum()),'R_faint':int(RL.sum()),'RR0':rr0.tolist()}


def sky_jackknife_regions(ra,region,njack):
    out=np.full(len(ra),-1,int); regs=sorted(set(region)); per=max(2,njack//max(len(regs),1)); off=0
    for reg in regs:
        q=np.where(region==reg)[0]
        if len(q)==0: continue
        order=q[np.argsort(ra[q])]; chunks=np.array_split(order,per)
        for c in chunks:
            out[c]=off; off+=1
    return out,off


def load_templates(path,s,z_eff):
    a=np.genfromtxt(path,delimiter=',',names=True); zvals=np.unique(a['z']); W=[]; D=[]
    for z in zvals:
        q=a['z']==z
        W.append(np.interp(s,a['s_Mpc_over_h'][q],a['wake_shape'][q]))
        D.append(np.interp(s,a['s_Mpc_over_h'][q],a['doppler_shape'][q]))
    W=np.asarray(W); D=np.asarray(D)
    if z_eff<=zvals[0]: return W[0],D[0]
    if z_eff>=zvals[-1]: return W[-1],D[-1]
    j=np.searchsorted(zvals,z_eff)-1; t=(z_eff-zvals[j])/(zvals[j+1]-zvals[j])
    return (1-t)*W[j]+t*W[j+1],(1-t)*D[j]+t*D[j+1]


def gls_fit(y,cov,templates,names):
    Cinv=np.linalg.pinv(cov,rcond=1e-8); X=np.column_stack(templates)
    F=X.T@Cinv@X; Finv=np.linalg.pinv(F,rcond=1e-10); beta=Finv@(X.T@Cinv@y)
    resid=y-X@beta; ch=float(resid@Cinv@resid); dof=max(len(y)-np.linalg.matrix_rank(X),0)
    errs=np.sqrt(np.maximum(np.diag(Finv),0)); out={n:{'amplitude':float(b),'sigma':float(e),'z':float(b/e) if e>0 else None} for n,b,e in zip(names,beta,errs)}
    den=np.sqrt(np.maximum(np.outer(np.diag(Finv),np.diag(Finv)),1e-300))
    out.update({'chi2':ch,'dof':int(dof),'pvalue':float(chi2.sf(ch,dof)) if dof>0 else None,'parameter_correlation':(Finv/den).tolist()})
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',nargs='+',required=True); ap.add_argument('--random',nargs='+',required=True); ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40); ap.add_argument('--dz-split',type=float,default=0.02)
    ap.add_argument('--max-per-class',type=int,default=30000); ap.add_argument('--max-random-per-class',type=int,default=60000)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--jackknife',type=int,default=8); ap.add_argument('--permutations',type=int,default=16)
    ap.add_argument('--seed',type=int,default=20260913); ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(args.seed)
    edges=np.array([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:])
    D=read_catalog(args.data,'data'); R=read_catalog(args.random,'random')
    iD,lD,sD,thresholds=make_data_split(D,args.zmin,args.zmax,args.dz_split)
    iD,lD,sD=stratified_subsample(iD,lD,sD,args.max_per_class,rng)
    iR,lR,sR=split_randoms(R,args.zmin,args.zmax,args.dz_split,thresholds,args.max_random_per_class,rng)
    xi0,xi1,counts=estimator(D,iD,lD,R,iR,lR,edges,args.theta_min_deg)
    z_eff=float(np.average(D['z'][iD],weights=D['w'][iD]))

    regD,nD=sky_jackknife_regions(D['ra'][iD],D['region'][iD],args.jackknife)
    regR,nR=sky_jackknife_regions(R['ra'][iR],R['region'][iR],args.jackknife)
    nj=min(nD,nR); jk=[]
    for j in range(nj):
        _,v,_=estimator(D,iD,lD,R,iR,lR,edges,args.theta_min_deg,regD!=j,regR!=j); jk.append(v); print('JK',j,json.dumps(v.tolist()))
    jk=np.asarray(jk); jm=np.mean(jk,axis=0); cov=(nj-1)/nj*((jk-jm).T@(jk-jm)) if nj>1 else np.diag(np.full(len(s),1e6))
    floor=max(float(np.median(np.diag(cov)))*1e-6,1e-14); cov=cov+np.eye(len(s))*floor

    null=[]
    for ip in range(args.permutations):
        lp=lD.copy()
        for st in np.unique(sD):
            q=np.where(sD==st)[0]; lp[q]=rng.permutation(lp[q])
        _,v,_=estimator(D,iD,lp,R,iR,lR,edges,args.theta_min_deg); null.append(v); print('PERM',ip,json.dumps(v.tolist()))
    null=np.asarray(null); null_mean=np.mean(null,axis=0); null_cov=np.cov(null,rowvar=False,ddof=1) if len(null)>1 else np.diag(np.ones(len(s)))
    delta=xi1-null_mean; ni=np.linalg.pinv(null_cov,rcond=1e-8); null_chi=float(delta@ni@delta); null_p=float(chi2.sf(null_chi,len(s)))

    wake,dop=load_templates(args.templates,s,z_eff)
    fit2=gls_fit(xi1,cov,[dop,wake],['doppler','wake_matched_filter'])
    b1=(s/60.0)**-1; b2=(s/60.0)**-2
    fit4=gls_fit(xi1,cov,[dop,b1,b2,wake],['doppler','broadband_inv_s','broadband_inv_s2','wake_matched_filter'])
    cinv=np.linalg.pinv(cov,rcond=1e-8); chi0=float(xi1@cinv@xi1); p0=float(chi2.sf(chi0,len(s)))

    np.savetxt(out/'data_vector.csv',np.column_stack([s,xi0,xi1,np.sqrt(np.maximum(np.diag(cov),0)),null_mean]),delimiter=',',header='s_Mpc_over_h,xi0_cross,xi1_odd,jackknife_sigma,permutation_null_mean',comments='')
    np.savetxt(out/'jackknife_vectors.csv',jk,delimiter=','); np.savetxt(out/'permutation_vectors.csv',null,delimiter=','); np.savetxt(out/'jackknife_covariance.csv',cov,delimiter=',')
    summary={
      'scope':'DESI DR1 BGS clustering-ready real-survey odd-sector deployment; not a DESI collaboration likelihood or detection claim.',
      'catalogs_data':[Path(x).name for x in args.data],'catalogs_random':[Path(x).name for x in args.random],
      'z_range':[args.zmin,args.zmax],'z_effective':z_eff,'luminosity_split':'weighted median dereddened r magnitude independently in dz bins and survey cap',
      'dz_split':args.dz_split,'theta_min_deg':args.theta_min_deg,'separation_edges_Mpc_over_h':edges.tolist(),
      'counts':counts,'split_thresholds':thresholds,'seed':args.seed,
      'xi1_odd':xi1.tolist(),'xi0_cross':xi0.tolist(),'jackknife_regions':int(nj),'jackknife_covariance':cov.tolist(),
      'zero_vector_chi2':chi0,'zero_vector_dof':len(s),'zero_vector_pvalue':p0,
      'permutation_null_count':args.permutations,'permutation_null_mean':null_mean.tolist(),'permutation_null_chi2':null_chi,'permutation_null_dof':len(s),'permutation_null_pvalue':null_p,
      'fit_doppler_plus_wake':fit2,'fit_conservative_odd_plus_wake':fit4,
      'wake_amplitude_scope':'Matched-filter coefficient for the production hidden-state wake shape. It is not yet an absolute DESI wake amplitude because the full survey window is not convolved into the theory template.',
      'required_upgrade_for_absolute_constraint':'Window-convolved theory template and split-specific mock covariance/validation (Abacus luminosity split; EZmock survey-window control).'
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print('PHASE7_DESI',json.dumps(summary,indent=2))

if __name__=='__main__': main()
