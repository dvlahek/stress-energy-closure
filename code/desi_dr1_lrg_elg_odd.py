#!/usr/bin/env python3
"""Genuine DESI DR1 LRG x ELG odd-multipole estimator.

Measures oriented LRG->ELG dipole and octupole in the genuine overlap
0.8 < z < 1.1 using independent clustering-ready data/random catalogs.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree

COSMO=FlatLambdaCDM(H0=67.4,Om0=0.315,Tcmb0=2.7255)
H=0.674

def col(d,*names,default=None):
    a={n.upper():n for n in d.names}
    for n in names:
        if n.upper() in a: return np.asarray(d[a[n.upper()]])
    if default is not None: return default
    raise KeyError(names)

def read(paths):
    chunks=[]
    for fn in paths:
        reg='NGC' if 'NGC' in Path(fn).name.upper() else 'SGC'
        with fits.open(fn,memmap=True) as h:
            d=h[1].data
            ra=col(d,'RA').astype(float); dec=col(d,'DEC').astype(float)
            z=col(d,'Z').astype(float)
            w=col(d,'WEIGHT',default=np.ones(len(d))).astype(float)
            wf=col(d,'WEIGHT_FKP',default=np.ones(len(d))).astype(float)
            chunks.append(dict(ra=ra,dec=dec,z=z,w=w*wf,
                               region=np.full(len(d),reg,dtype='U3')))
    return {k:np.concatenate([x[k] for x in chunks]) for k in chunks[0]}

def sel(c,z0,z1):
    g=np.isfinite(c['z'])&np.isfinite(c['w'])&(c['w']>0)&(c['z']>=z0)&(c['z']<z1)
    return {k:v[g] for k,v in c.items()}

def sub(c,n,rng):
    if n<=0 or len(c['z'])<=n:return c
    I=[]
    for reg in np.unique(c['region']):
        q=np.where(c['region']==reg)[0]
        k=max(1,min(len(q),int(round(n*len(q)/len(c['z'])))))
        I.extend(rng.choice(q,k,replace=False).tolist())
    I=np.asarray(sorted(I))
    return {k:v[I] for k,v in c.items()}

def xyz(c):
    r=np.asarray(COSMO.comoving_distance(c['z']).value)*H
    ra=np.deg2rad(c['ra']); de=np.deg2rad(c['dec']); cd=np.cos(de)
    x=np.c_[r*cd*np.cos(ra),r*cd*np.sin(ra),r*np.sin(de)]
    return x,x/np.maximum(np.linalg.norm(x,axis=1)[:,None],1e-300)

def jk_edges(A,B,nj):
    out={}; caps=sorted(set(A['region'])|set(B['region']))
    per=max(2,nj//len(caps))
    for reg in caps:
        ra=np.r_[A['ra'][A['region']==reg],B['ra'][B['region']==reg]]
        e=np.quantile(ra,np.linspace(0,1,per+1)); e[0]-=1e-6; e[-1]+=1e-6
        out[reg]=e
    return out

def jk_assign(c,E):
    a=np.full(len(c['z']),-1,np.int16); off=0
    for reg,e in E.items():
        q=np.where(c['region']==reg)[0]
        a[q]=off+np.clip(np.searchsorted(e,c['ra'][q],side='right')-1,0,len(e)-2)
        off+=len(e)-1
    return a,off

def prep(c,jk):
    x,u=xyz(c)
    return dict(x=x,u=u,w=c['w'].astype(float),z=c['z'].astype(float),jk=jk)

def pairs(A,B,edges,K,rng,theta,label):
    tr=cKDTree(B['x']); r0,r1=edges[0],edges[-1]
    aa=[];bb=[];bi=[];mu=[];imp=[];zp=[]
    for st in range(0,len(A['x']),512):
        ns=tr.query_ball_point(A['x'][st:st+512],r1,workers=-1)
        for loc,j0 in enumerate(ns):
            i=st+loc
            if not j0: continue
            j=np.asarray(j0,int); dv=B['x'][j]-A['x'][i]; rr=np.linalg.norm(dv,axis=1)
            g=(rr>=r0)&(rr<r1)
            if theta>0:
                th=np.rad2deg(np.arccos(np.clip(np.sum(A['u'][i]*B['u'][j],axis=1),-1,1)))
                g&=th>=theta
            j=j[g]; rr=rr[g]; dv=dv[g]
            if not len(j): continue
            nall=len(j)
            if nall>K:
                t=rng.choice(nall,K,replace=False); j=j[t];rr=rr[t];dv=dv[t]
            n=len(j); mid=B['x'][j]+A['x'][i]
            mm=np.sum(dv*mid,axis=1)/(np.maximum(rr,1e-12)*np.maximum(np.linalg.norm(mid,axis=1),1e-12))
            aa.append(np.full(n,i,np.int32));bb.append(j.astype(np.int32))
            bi.append((np.searchsorted(edges,rr,side='right')-1).astype(np.int16))
            mu.append(mm.astype(np.float32));imp.append(np.full(n,nall/n,np.float32))
            zp.append((0.5*(A['z'][i]+B['z'][j])).astype(np.float32))
        if st and st%(512*40)==0: print(label,st,'/',len(A['x']),flush=True)
    return dict(a=np.concatenate(aa),b=np.concatenate(bb),bin=np.concatenate(bi),
                mu=np.concatenate(mu).astype(float),imp=np.concatenate(imp).astype(float),
                zpair=np.concatenate(zp).astype(float))

def W(X,drop):
    g=np.ones(len(X['w']),bool) if drop is None else X['jk']!=drop
    return float(X['w'][g].sum())

def ph(t,A,B,nb,ell,drop):
    g=np.ones(len(t['a']),bool)
    if drop is not None:g&=(A['jk'][t['a']]!=drop)&(B['jk'][t['b']]!=drop)
    a=t['a'][g];b=t['b'][g];ib=t['bin'][g];m=t['mu'][g]
    pw=A['w'][a]*B['w'][b]*t['imp'][g]
    ang=np.ones_like(m) if ell==0 else 3*m if ell==1 else 3.5*(5*m**3-3*m)
    return np.bincount(ib,weights=pw*ang,minlength=nb)

def est(BL,ell,drop=None):
    A,B,RA,RB=BL['A'],BL['B'],BL['RA'],BL['RB']; nb=BL['nb']
    DD=ph(BL['DD'],A,B,nb,ell,drop)/(W(A,drop)*W(B,drop))
    DR=ph(BL['DR'],A,RB,nb,ell,drop)/(W(A,drop)*W(RB,drop))
    RD=ph(BL['RD'],RA,B,nb,ell,drop)/(W(RA,drop)*W(B,drop))
    RR0=ph(BL['RR'],RA,RB,nb,0,drop)/(W(RA,drop)*W(RB,drop))
    RR=ph(BL['RR'],RA,RB,nb,ell,drop)/(W(RA,drop)*W(RB,drop))
    return (DD-DR-RD+RR)/np.maximum(RR0,1e-300)

def cov(j):
    j=np.asarray(j);n=len(j);m=j.mean(0);C=(n-1)/n*((j-m).T@(j-m))
    ev=np.linalg.eigvalsh(C);mx=max(float(ev[-1]),1e-300)
    r=max(mx*1e-7,float(np.median(np.diag(C)))*1e-6,1e-14)
    return C+np.eye(C.shape[0])*r,r

def main():
    q=argparse.ArgumentParser()
    for x in ['lrg-data','elg-data','lrg-random','elg-random']:
        q.add_argument('--'+x,nargs='+',required=True)
    q.add_argument('--outdir',required=True);q.add_argument('--zmin',type=float,default=.8)
    q.add_argument('--zmax',type=float,default=1.1);q.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    q.add_argument('--neighbors-per-anchor',type=int,default=48);q.add_argument('--theta-min-deg',type=float,default=.05)
    q.add_argument('--jackknife',type=int,default=30);q.add_argument('--random-factor',type=float,default=2.)
    q.add_argument('--max-data-per-tracer',type=int,default=0);q.add_argument('--seed',type=int,default=20260918)
    a=q.parse_args(); rng=np.random.default_rng(a.seed); out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    ed=np.asarray([float(x) for x in a.sep_edges.split(',')]);s=.5*(ed[:-1]+ed[1:]);nb=len(s)
    L=sub(sel(read(a.lrg_data),a.zmin,a.zmax),a.max_data_per_tracer,rng)
    E=sub(sel(read(a.elg_data),a.zmin,a.zmax),a.max_data_per_tracer,rng)
    LR=sub(sel(read(a.lrg_random),a.zmin,a.zmax),int(a.random_factor*len(L['z'])),rng)
    ER=sub(sel(read(a.elg_random),a.zmin,a.zmax),int(a.random_factor*len(E['z'])),rng)
    print('counts',len(L['z']),len(E['z']),len(LR['z']),len(ER['z']),flush=True)
    je=jk_edges(L,E,a.jackknife);jl,nj=jk_assign(L,je);je2,_=jk_assign(E,je);jlr,_=jk_assign(LR,je);jer,_=jk_assign(ER,je)
    A=prep(L,jl);B=prep(E,je2);RA=prep(LR,jlr);RB=prep(ER,jer)
    BL=dict(A=A,B=B,RA=RA,RB=RB,nb=nb)
    BL['DD']=pairs(A,B,ed,a.neighbors_per_anchor,rng,a.theta_min_deg,'DD')
    BL['DR']=pairs(A,RB,ed,a.neighbors_per_anchor,rng,a.theta_min_deg,'DR')
    BL['RD']=pairs(RA,B,ed,a.neighbors_per_anchor,rng,a.theta_min_deg,'RD')
    BL['RR']=pairs(RA,RB,ed,a.neighbors_per_anchor,rng,a.theta_min_deg,'RR')
    x1=est(BL,1);x3=est(BL,3);j1=[];j3=[]
    for j in range(nj):
        j1.append(est(BL,1,j));j3.append(est(BL,3,j));print('JK',j,flush=True)
    C1,r1=cov(j1);C3,r3=cov(j3)
    pw=A['w'][BL['DD']['a']]*B['w'][BL['DD']['b']]*BL['DD']['imp']
    ze=float(np.sum(pw*BL['DD']['zpair'])/np.sum(pw))
    np.savetxt(out/'lrg_elg_odd_multipoles.csv',np.c_[s,x1,np.sqrt(np.diag(C1)),x3,np.sqrt(np.diag(C3))],
               delimiter=',',header='s_Mpc_over_h,xi1_lrg_to_elg,sigma_xi1,xi3_lrg_to_elg,sigma_xi3',comments='')
    np.savetxt(out/'lrg_elg_dipole_cov.csv',C1,delimiter=',');np.savetxt(out/'lrg_elg_octupole_cov.csv',C3,delimiter=',')
    S=dict(scope='genuine DESI DR1 LRG->ELG cross-population odd multipoles',z_range=[a.zmin,a.zmax],
           z_effective_pair_weighted=ze,counts={'LRG':len(L['z']),'ELG':len(E['z']),'LRG_random':len(LR['z']),'ELG_random':len(ER['z'])},
           jackknife_regions=nj,seed=a.seed,dipole_ridge=r1,octupole_ridge=r3,
           guardrail='data-vector measurement only; wake interpretation requires the separate frozen shape fit')
    (out/'summary_lrg_elg.json').write_text(json.dumps(S,indent=2)+'\\n');print(json.dumps(S,indent=2))

if __name__=='__main__':main()
