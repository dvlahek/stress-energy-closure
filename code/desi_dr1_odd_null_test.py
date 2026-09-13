#!/usr/bin/env python3
"""End-to-end parity-odd null test on the public DESI DR1 Gfinder group VAC.

This is deliberately a pipeline demonstration, not a neutrino-wake likelihood.
The public group catalogue supplies sky positions, redshifts and a halo-mass
proxy (GRP_LOGM).  We split at log10(h M/Msun)=13.75, match the high- and
low-mass samples in narrow redshift bins, and measure an oriented cross-pair
P1(mu) moment in separation bins.  A permutation null is constructed by
shuffling the mass labels *within the same redshift bins*, preserving the sky
and redshift selection.

A statistically ordinary result demonstrates that the parity-odd estimator
family can be run on public DESI-linked data without generating a spurious
large odd signal.  It is not a substitute for a DESI collaboration LSS
likelihood, random-catalog correction, survey window treatment or a wake
constraint.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree
from scipy.stats import chi2


def xyz_from_radec_z(ra,dec,z):
    cosmo=FlatLambdaCDM(H0=67.4,Om0=0.315,Tcmb0=2.7255)
    h=0.674
    r=np.asarray(cosmo.comoving_distance(z).value)*h
    rr=np.deg2rad(ra); dd=np.deg2rad(dec); cd=np.cos(dd)
    return np.column_stack((r*cd*np.cos(rr),r*cd*np.sin(rr),r*np.sin(dd)))


def matched_indices(z,mass,split,zmin,zmax,dz,max_per_class,rng):
    edges=np.arange(zmin,zmax+0.5*dz,dz)
    chosen_hi=[]; chosen_lo=[]; bin_ids_hi=[]; bin_ids_lo=[]
    counts=[]
    for ib,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        zh=np.where((z>=lo)&(z<hi)&(mass>=split))[0]
        zl=np.where((z>=lo)&(z<hi)&(mass<split))[0]
        n=min(len(zh),len(zl)); counts.append((ib,lo,hi,len(zh),len(zl),n))
    ntot=sum(x[-1] for x in counts)
    frac=min(1.0,max_per_class/max(ntot,1))
    for ib,lo,hi,nh,nl,n in counts:
        k=int(np.floor(n*frac))
        if n>0 and k==0 and len(chosen_hi)<max_per_class: k=1
        if k<=0: continue
        zh=np.where((z>=lo)&(z<hi)&(mass>=split))[0]
        zl=np.where((z>=lo)&(z<hi)&(mass<split))[0]
        ih=rng.choice(zh,size=min(k,len(zh)),replace=False); il=rng.choice(zl,size=min(k,len(zl)),replace=False)
        kk=min(len(ih),len(il)); ih=ih[:kk]; il=il[:kk]
        chosen_hi.extend(ih.tolist()); chosen_lo.extend(il.tolist())
        bin_ids_hi.extend([ib]*kk); bin_ids_lo.extend([ib]*kk)
    return np.array(chosen_hi,int),np.array(chosen_lo,int),np.array(bin_ids_hi,int),np.array(bin_ids_lo,int),counts


def odd_pair_stat(xh,xl,sep_edges):
    th=cKDTree(xh); tl=cKDTree(xl)
    coo=th.sparse_distance_matrix(tl,float(sep_edges[-1]),output_type='coo_matrix')
    if coo.nnz==0: return np.full(len(sep_edges)-1,np.nan),np.zeros(len(sep_edges)-1,int)
    row=np.asarray(coo.row); col=np.asarray(coo.col); sep=np.asarray(coo.data)
    keep=(sep>=sep_edges[0])&(sep<sep_edges[-1]); row=row[keep]; col=col[keep]; sep=sep[keep]
    dv=xh[row]-xl[col]; mid=xh[row]+xl[col]
    midn=np.linalg.norm(mid,axis=1); mu=np.sum(dv*mid,axis=1)/(np.maximum(sep,1e-12)*np.maximum(midn,1e-12))
    vals=[]; counts=[]
    for lo,hi in zip(sep_edges[:-1],sep_edges[1:]):
        q=(sep>=lo)&(sep<hi); counts.append(int(q.sum()))
        vals.append(float(3.0*np.mean(mu[q])) if q.any() else np.nan)
    return np.array(vals,float),np.array(counts,int)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--fits',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--split',type=float,default=13.75); ap.add_argument('--zmin',type=float,default=0.05); ap.add_argument('--zmax',type=float,default=0.40)
    ap.add_argument('--dz-match',type=float,default=0.02); ap.add_argument('--max-per-class',type=int,default=12000)
    ap.add_argument('--permutations',type=int,default=40); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(args.seed)
    with fits.open(args.fits,memmap=True) as hdul:
        d=hdul[1].data; ra=np.asarray(d['GRP_RA'],float); dec=np.asarray(d['GRP_DEC'],float); z=np.asarray(d['GRP_Z'],float); mass=np.asarray(d['GRP_LOGM'],float)
    good=np.isfinite(ra)&np.isfinite(dec)&np.isfinite(z)&np.isfinite(mass)&(z>=args.zmin)&(z<args.zmax)&(mass>10)&(mass<16)
    ra,dec,z,mass=ra[good],dec[good],z[good],mass[good]
    ih,il,bh,bl,raw_counts=matched_indices(z,mass,args.split,args.zmin,args.zmax,args.dz_match,args.max_per_class,rng)
    if len(ih)<1000: raise RuntimeError(f'too few matched high-mass objects: {len(ih)}')
    # Build a combined, exactly redshift-stratified sample.
    pos_h=xyz_from_radec_z(ra[ih],dec[ih],z[ih]); pos_l=xyz_from_radec_z(ra[il],dec[il],z[il])
    sep_edges=np.array([20.,40.,60.,80.,100.])
    data_stat,pair_counts=odd_pair_stat(pos_h,pos_l,sep_edges)
    if not np.all(np.isfinite(data_stat)): raise RuntimeError('empty pair bin in data statistic')
    pos=np.vstack([pos_h,pos_l]); bins=np.concatenate([bh,bl]); n=len(ih)
    null=[]
    unique=np.unique(bins)
    for ip in range(args.permutations):
        hi_mask=np.zeros(2*n,dtype=bool)
        for b in unique:
            idx=np.where(bins==b)[0]; k=len(idx)//2
            pick=rng.choice(idx,size=k,replace=False); hi_mask[pick]=True
        ph=pos[hi_mask]; pl=pos[~hi_mask]
        # exact equality can differ by one only for pathological strata; matched construction makes it equal.
        stat,_=odd_pair_stat(ph,pl,sep_edges); null.append(stat)
        print('PERM',ip,json.dumps(stat.tolist()))
    null=np.asarray(null,float); null_mean=np.nanmean(null,axis=0); cov=np.cov(null,rowvar=False,ddof=1)
    std=np.sqrt(np.maximum(np.diag(cov),1e-300)); zscore=(data_stat-null_mean)/std
    delta=data_stat-null_mean; inv=np.linalg.pinv(cov,rcond=1e-10); chisq=float(delta@inv@delta); dof=len(data_stat); pval=float(chi2.sf(chisq,dof))
    summary={'catalog':'DESI DR1 Gfinder group VAC: DESIDR9.y1.v1_group.fits','statistic':'3 <mu> oriented high-minus-low cross-pair P1 moment',
             'interpretation_scope':'End-to-end parity-odd real-data null/smoke test only; not a DESI LSS wake likelihood.',
             'mass_split_log10_hM_Msun':args.split,'z_range':[args.zmin,args.zmax],'redshift_matching_width':args.dz_match,
             'matched_objects_per_class':int(n),'separation_edges_Mpc_over_h':sep_edges.tolist(),'data_dipole':data_stat.tolist(),
             'pair_counts':pair_counts.tolist(),'null_permutations':args.permutations,'null_mean':null_mean.tolist(),'null_std':std.tolist(),
             'per_bin_zscore':zscore.tolist(),'global_chi2':chisq,'global_dof':dof,'global_chi2_pvalue':pval,
             'seed':args.seed,'raw_selected_groups':int(len(z))}
    np.savetxt(out/'null_permutations.csv',null,delimiter=',',header=','.join([f'{sep_edges[i]:g}-{sep_edges[i+1]:g}' for i in range(len(sep_edges)-1)]),comments='')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print('DESI_NULL',json.dumps(summary,indent=2))

if __name__=='__main__': main()
