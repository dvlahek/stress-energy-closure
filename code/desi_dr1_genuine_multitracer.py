#!/usr/bin/env python3
"""DESI DR1 genuine cross-population odd-dipole test.

Primary: LRG x ELG_LOPnotqso over 0.8 <= z < 1.1.
Secondary: LRG x BGS_BRIGHT-21.5 across the z=0.4 boundary.

The tracer identity itself supplies the antisymmetric mark.  No luminosity
subdivision is used.  Geometry and pair sampling are frozen before the odd
statistic is evaluated.  This is an exploratory public-data validation, not a
DESI collaboration likelihood or a detection claim.
"""
from __future__ import annotations
import argparse, json, math
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
        if n.upper() in avail:
            return np.asarray(d[avail[n.upper()]])
    if default is not None:
        return default
    raise KeyError(f"missing columns {names}; available={list(d.names)}")


def read_catalog(paths,zmin,zmax):
    chunks=[]
    for path in paths:
        region='NGC' if 'NGC' in Path(path).name.upper() else 'SGC'
        with fits.open(path,memmap=True) as h:
            d=h[1].data
            ra=col(d,'RA').astype(float)
            dec=col(d,'DEC').astype(float)
            z=col(d,'Z').astype(float)
            w=col(d,'WEIGHT',default=np.ones(len(d))).astype(float)
            wf=col(d,'WEIGHT_FKP',default=np.ones(len(d))).astype(float)
            good=(np.isfinite(ra)&np.isfinite(dec)&np.isfinite(z)&np.isfinite(w)&
                  np.isfinite(wf)&(w>0)&(wf>0)&(z>=zmin)&(z<zmax))
            chunks.append({
                'ra':ra[good],'dec':dec[good],'z':z[good],
                'w':(w[good]*wf[good]),
                'region':np.full(int(good.sum()),region,dtype='U3')
            })
    if not chunks:
        raise RuntimeError("no catalog chunks")
    return {k:np.concatenate([x[k] for x in chunks]) for k in chunks[0]}


def distance_xyz(ra,dec,z):
    r=np.asarray(COSMO.comoving_distance(z).value)*H
    rr=np.deg2rad(ra); dd=np.deg2rad(dec); cd=np.cos(dd)
    return np.column_stack((r*cd*np.cos(rr),r*cd*np.sin(rr),r*np.sin(dd)))


def stratified_random_subset(R,D,zmin,zmax,dz,factor,rng):
    """Subsample official randoms to a fixed multiple of data density per cap/z stratum."""
    out=[]
    edges=np.arange(zmin,zmax+0.5*dz,dz)
    for cap in ('NGC','SGC'):
        for lo,hi in zip(edges[:-1],edges[1:]):
            qd=np.where((D['region']==cap)&(D['z']>=lo)&(D['z']<hi))[0]
            qr=np.where((R['region']==cap)&(R['z']>=lo)&(R['z']<hi))[0]
            if len(qd)==0 or len(qr)==0:
                continue
            target=min(len(qr),max(100,int(math.ceil(factor*len(qd)))))
            if target<len(qr):
                qr=rng.choice(qr,size=target,replace=False)
            out.append(np.asarray(qr,int))
    if not out:
        raise RuntimeError("no randoms survive stratified subset")
    idx=np.concatenate(out)
    return {k:np.asarray(v)[idx] for k,v in R.items()}


def common_jackknife(A,B,RA,RB,njack):
    """Common RA-quantile sky regions from pooled data geometry, then apply to randoms."""
    regs={}
    per=max(2,njack//2)
    rid=0
    bounds={}
    for cap in ('NGC','SGC'):
        vals=np.concatenate([
            A['ra'][A['region']==cap],
            B['ra'][B['region']==cap]
        ])
        if len(vals)==0:
            continue
        q=np.linspace(0,1,per+1)
        edges=np.quantile(vals,q)
        edges[0]-=1e-8; edges[-1]+=1e-8
        # Remove accidental duplicate quantiles but keep at least two regions if possible.
        edges=np.unique(edges)
        if len(edges)<3:
            edges=np.linspace(float(vals.min())-1e-8,float(vals.max())+1e-8,3)
        bounds[cap]=edges
        for _ in range(len(edges)-1):
            regs[rid]=cap
            rid+=1

    def assign(X):
        out=np.full(len(X['ra']),-1,int)
        off=0
        for cap in ('NGC','SGC'):
            if cap not in bounds:
                continue
            edges=bounds[cap]
            q=np.where(X['region']==cap)[0]
            b=np.searchsorted(edges,X['ra'][q],side='right')-1
            b=np.clip(b,0,len(edges)-2)
            out[q]=off+b
            off+=len(edges)-1
        return out

    return assign(A),assign(B),assign(RA),assign(RB),rid,bounds


def prepare(X,jk):
    x=distance_xyz(X['ra'],X['dec'],X['z'])
    u=x/np.maximum(np.linalg.norm(x,axis=1)[:,None],1e-300)
    return {
        'x':x,'u':u,'w':np.asarray(X['w'],float),'z':np.asarray(X['z'],float),
        'ra':np.asarray(X['ra'],float),'region':np.asarray(X['region']),
        'jk':np.asarray(jk,int)
    }


def sample_pairs(A,B,edges,kmax,rng,theta_min_deg,chunk=512,label='PAIR'):
    tree=cKDTree(B['x']); rmax=float(edges[-1]); rmin=float(edges[0])
    aa=[]; bb=[]; bins=[]; mus=[]; imps=[]
    n=len(A['x'])
    for start in range(0,n,chunk):
        stop=min(start+chunk,n)
        neigh=tree.query_ball_point(A['x'][start:stop],rmax,workers=-1,return_sorted=False)
        for loc,js0 in enumerate(neigh):
            i=start+loc
            if not js0:
                continue
            js=np.asarray(js0,int)
            dv=B['x'][js]-A['x'][i]
            rr=np.linalg.norm(dv,axis=1)
            keep=(rr>=rmin)&(rr<rmax)
            if theta_min_deg>0:
                c=np.sum(A['u'][i]*B['u'][js],axis=1)
                th=np.rad2deg(np.arccos(np.clip(c,-1,1)))
                keep&=(th>=theta_min_deg)
            js=js[keep]; rr=rr[keep]; dv=dv[keep]
            if len(js)==0:
                continue
            nall=len(js)
            if nall>kmax:
                take=rng.choice(nall,size=kmax,replace=False)
                js=js[take]; rr=rr[take]; dv=dv[take]
            ns=len(js); imp=float(nall/ns)
            mid=B['x'][js]+A['x'][i]
            mu=np.sum(dv*mid,axis=1)/(np.maximum(rr,1e-12)*np.maximum(np.linalg.norm(mid,axis=1),1e-12))
            ib=np.searchsorted(edges,rr,side='right')-1
            aa.append(np.full(ns,i,np.int32)); bb.append(js.astype(np.int32))
            bins.append(ib.astype(np.int16)); mus.append(mu.astype(np.float32))
            imps.append(np.full(ns,imp,np.float32))
        if start and start%(chunk*20)==0:
            print(label,'anchors',start,'/',n)
    if not aa:
        raise RuntimeError(f"no sampled pairs for {label}")
    out={
        'a':np.concatenate(aa),'b':np.concatenate(bb),'bin':np.concatenate(bins),
        'mu':np.concatenate(mus).astype(float),'imp':np.concatenate(imps).astype(float)
    }
    print(label,'sampled_pairs',len(out['a']),'estimated_pairs',float(out['imp'].sum()))
    return out


def pair_hist(tab,A,B,nb,ell,drop=None):
    keep=np.ones(len(tab['a']),bool)
    if drop is not None:
        keep&=(A['jk'][tab['a']]!=drop)&(B['jk'][tab['b']]!=drop)
    a=tab['a'][keep]; b=tab['b'][keep]; ib=tab['bin'][keep]
    pw=A['w'][a]*B['w'][b]*tab['imp'][keep]
    mu=tab['mu'][keep]
    if ell==0:
        basis=np.ones_like(mu)
    elif ell==1:
        basis=3.0*mu
    elif ell==3:
        basis=3.5*(5.0*mu**3-3.0*mu)
    else:
        raise ValueError("ell must be 0,1,3")
    return np.bincount(ib,weights=pw*basis,minlength=nb).astype(float)


def wsum(X,drop=None):
    keep=np.ones(len(X['w']),bool) if drop is None else X['jk']!=drop
    return float(X['w'][keep].sum())


def estimate(block,ell,drop=None):
    A,B,RA,RB=block['A'],block['B'],block['RA'],block['RB']
    nb=block['nb']
    dd=pair_hist(block['DD'],A,B,nb,ell,drop)/(max(wsum(A,drop)*wsum(B,drop),1e-300))
    dr=pair_hist(block['DR'],A,RB,nb,ell,drop)/(max(wsum(A,drop)*wsum(RB,drop),1e-300))
    rd=pair_hist(block['RD'],RA,B,nb,ell,drop)/(max(wsum(RA,drop)*wsum(B,drop),1e-300))
    rr=pair_hist(block['RR'],RA,RB,nb,ell,drop)/(max(wsum(RA,drop)*wsum(RB,drop),1e-300))
    rr0=pair_hist(block['RR'],RA,RB,nb,0,drop)/(max(wsum(RA,drop)*wsum(RB,drop),1e-300))
    return (dd-dr-rd+rr)/np.maximum(rr0,1e-300)


def load_template(path,s,z):
    a=np.genfromtxt(path,delimiter=',',names=True)
    zvals=np.unique(a['z'])
    W=[]; D=[]
    for zv in zvals:
        q=a['z']==zv
        W.append(np.interp(s,a['s_Mpc_over_h'][q],a['wake_shape'][q]))
        D.append(np.interp(s,a['s_Mpc_over_h'][q],a['doppler_shape'][q]))
    W=np.asarray(W); D=np.asarray(D)
    if z<=zvals[0]:
        return W[0],D[0]
    if z>=zvals[-1]:
        return W[-1],D[-1]
    j=np.searchsorted(zvals,z)-1
    t=(z-zvals[j])/(zvals[j+1]-zvals[j])
    return (1-t)*W[j]+t*W[j+1],(1-t)*D[j]+t*D[j+1]


def gls_fit(y,cov,templates,names):
    ci=np.linalg.pinv(cov,rcond=1e-10)
    X=np.column_stack(templates)
    F=X.T@ci@X
    Fi=np.linalg.pinv(F,rcond=1e-10)
    beta=Fi@(X.T@ci@y)
    err=np.sqrt(np.maximum(np.diag(Fi),0))
    resid=y-X@beta
    rank=int(np.linalg.matrix_rank(X))
    dof=max(len(y)-rank,0)
    out={n:{'amplitude':float(b),'sigma':float(e),'z':float(b/e) if e>0 else None}
         for n,b,e in zip(names,beta,err)}
    out.update({'chi2':float(resid@ci@resid),'dof':dof,
                'pvalue':float(chi2.sf(float(resid@ci@resid),dof)) if dof>0 else None})
    return out


def block_template(v,iz,nz,nb):
    out=np.zeros(nz*nb,float)
    out[iz*nb:(iz+1)*nb]=v
    return out


def run_primary(args,rng):
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')])
    s=0.5*(edges[:-1]+edges[1:]); nb=len(s)
    zedges=np.asarray([0.80,0.95,1.10])

    DA=read_catalog(args.lrg_data,0.80,1.10)
    DB=read_catalog(args.elg_data,0.80,1.10)
    RA0=read_catalog(args.lrg_random,0.80,1.10)
    RB0=read_catalog(args.elg_random,0.80,1.10)
    RA0=stratified_random_subset(RA0,DA,0.80,1.10,args.dz_random,args.random_factor,rng)
    RB0=stratified_random_subset(RB0,DB,0.80,1.10,args.dz_random,args.random_factor,rng)

    jkA,jkB,jkRA,jkRB,nj,bounds=common_jackknife(DA,DB,RA0,RB0,args.jackknife)
    if nj<=len(s)*2+2:
        raise RuntimeError(f"insufficient common jackknife regions: {nj}")

    blocks=[]; zmeta=[]
    for iz,(lo,hi) in enumerate(zip(zedges[:-1],zedges[1:])):
        def sub(X,jk):
            q=(X['z']>=lo)&(X['z']<hi)
            return {k:np.asarray(v)[q] for k,v in X.items()},np.asarray(jk)[q]
        a,ja=sub(DA,jkA); b,jb=sub(DB,jkB); ra,jra=sub(RA0,jkRA); rb,jrb=sub(RB0,jkRB)
        if min(len(a['z']),len(b['z']),len(ra['z']),len(rb['z']))<100:
            raise RuntimeError(f"too few objects in primary z block {lo}-{hi}")
        A=prepare(a,ja); B=prepare(b,jb); RA=prepare(ra,jra); RB=prepare(rb,jrb)
        DD=sample_pairs(A,B,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label=f'LRGxELG_DD_z{iz}')
        DR=sample_pairs(A,RB,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label=f'LRGxELG_DR_z{iz}')
        RD=sample_pairs(RA,B,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label=f'LRGxELG_RD_z{iz}')
        RR=sample_pairs(RA,RB,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label=f'LRGxELG_RR_z{iz}')
        blocks.append({'A':A,'B':B,'RA':RA,'RB':RB,'DD':DD,'DR':DR,'RD':RD,'RR':RR,'nb':nb})
        za=float(np.average(A['z'],weights=A['w'])); zb=float(np.average(B['z'],weights=B['w']))
        zmeta.append({'zlo':float(lo),'zhi':float(hi),'z_effective':0.5*(za+zb),
                      'N_LRG':len(A['z']),'N_ELG':len(B['z']),
                      'N_RLRG':len(RA['z']),'N_RELG':len(RB['z'])})

    def vec(ell,drop=None):
        return np.concatenate([estimate(b,ell,drop) for b in blocks])

    dip=vec(1); octu=vec(3)
    jk_d=[]; jk_o=[]
    for j in range(nj):
        jk_d.append(vec(1,j)); jk_o.append(vec(3,j))
        print('JK_PRIMARY',j)
    jk_d=np.asarray(jk_d); jk_o=np.asarray(jk_o)

    def covariance(jk):
        jm=jk.mean(axis=0)
        raw=(nj-1)/nj*((jk-jm).T@(jk-jm))
        ev=np.linalg.eigvalsh(raw); mx=float(max(ev[-1],1e-300))
        ridge=max(mx*1e-7,float(np.median(np.diag(raw)))*1e-6,1e-14)
        cov=raw+np.eye(raw.shape[0])*ridge
        return raw,cov,ridge,int(np.linalg.matrix_rank(raw,tol=mx*1e-10)),float(np.linalg.cond(cov))

    raw,cov,ridge,rank,cond=covariance(jk_d)
    rawo,covo,ridgeo,ranko,condo=covariance(jk_o)
    ci=np.linalg.pinv(cov,rcond=1e-10)
    cio=np.linalg.pinv(covo,rcond=1e-10)

    wake=[]; dop=[]; nz=len(blocks)
    for m in zmeta:
        w,d=load_template(args.templates,s,m['z_effective'])
        wake.extend(w.tolist()); dop.extend(d.tolist())
    wake=np.asarray(wake); dop=np.asarray(dop)

    nuisance=[]; names=[]
    b1=(s/60.0)**-1; b2=(s/60.0)**-2
    for iz in range(nz):
        ds=dop[iz*nb:(iz+1)*nb]
        nuisance.extend([block_template(ds,iz,nz,nb),
                         block_template(b1,iz,nz,nb),
                         block_template(b2,iz,nz,nb)])
        names.extend([f'doppler_z{iz}',f'broadband_inv_s_z{iz}',f'broadband_inv_s2_z{iz}'])
    fit=gls_fit(dip,cov,nuisance+[wake],names+['wake_matched_filter'])

    t0=float(dip@ci@dip)
    p0=float(chi2.sf(t0,len(dip)))
    to=float(octu@cio@octu)
    po=float(chi2.sf(to,len(octu)))

    out=Path(args.outdir)/'primary_lrg_elg'
    out.mkdir(parents=True,exist_ok=True)
    rows=[]; sig=np.sqrt(np.maximum(np.diag(cov),0)); sigo=np.sqrt(np.maximum(np.diag(covo),0))
    k=0
    for m in zmeta:
        for si in s:
            rows.append((m['zlo'],m['zhi'],m['z_effective'],si,dip[k],sig[k],octu[k],sigo[k])); k+=1
    np.savetxt(out/'data_vector.csv',np.asarray(rows),delimiter=',',
               header='zlo,zhi,z_effective,s_Mpc_over_h,xi1_LRG_to_ELG,jk_sigma_dipole,xi3_control,jk_sigma_octupole',comments='')
    np.savetxt(out/'jackknife_covariance_dipole.csv',cov,delimiter=',')
    np.savetxt(out/'jackknife_covariance_octupole.csv',covo,delimiter=',')
    np.savetxt(out/'jackknife_vectors_dipole.csv',jk_d,delimiter=',')
    np.savetxt(out/'jackknife_vectors_octupole.csv',jk_o,delimiter=',')

    summary={
        'status':'EXPLORATORY_GENUINE_CROSS_POPULATION',
        'primary_pair':'LRG -> ELG_LOPnotqso',
        'orientation':'fixed by tracer identity before inspecting odd statistic',
        'release':'DESI DR1 LSS v1.5 clustering catalogs',
        'z_bins':zmeta,
        'separation_edges_Mpc_over_h':edges.tolist(),
        'theta_min_deg':args.theta_min_deg,
        'neighbors_per_anchor':args.neighbors_per_anchor,
        'random_factor':args.random_factor,
        'seed':args.seed,
        'vector_dimension':len(dip),
        'jackknife':{'regions':nj,'raw_rank':rank,'ridge':ridge,'condition':cond},
        'zero_vector':{'chi2':t0,'dof':len(dip),'pvalue_asymptotic':p0},
        'conservative_fit':fit,
        'odd_octupole_control':{'chi2':to,'dof':len(octu),'pvalue_asymptotic':po,
                                'jackknife_raw_rank':ranko,'condition':condo},
        'template_scope':'Frozen hidden-state shape with free amplitude and no BGS-specific absolute calibration.',
        'interpretation_guardrail':'Exploratory cross-population validation only. A publication-level coefficient requires joint tracer mocks/window closure.'
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('GENUINE_PRIMARY',json.dumps(summary,indent=2))
    return summary


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--lrg-data',nargs='+',required=True)
    ap.add_argument('--lrg-random',nargs='+',required=True)
    ap.add_argument('--elg-data',nargs='+',required=True)
    ap.add_argument('--elg-random',nargs='+',required=True)
    ap.add_argument('--templates',required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    ap.add_argument('--dz-random',type=float,default=0.02)
    ap.add_argument('--random-factor',type=float,default=2.0)
    ap.add_argument('--neighbors-per-anchor',type=int,default=48)
    ap.add_argument('--theta-min-deg',type=float,default=0.05)
    ap.add_argument('--jackknife',type=int,default=30)
    ap.add_argument('--seed',type=int,default=20260918)
    args=ap.parse_args()
    Path(args.outdir).mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    run_primary(args,rng)


if __name__=='__main__':
    main()
