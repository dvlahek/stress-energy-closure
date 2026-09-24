#!/usr/bin/env python3
"""Full-sample five-tracer DESI DR1 Phase-7 odd-sector estimator.

This is a sensitivity upgrade of the redshift-resolved Phase-7 baseline.  The
tracer definition is fixed before looking at the odd signal: within narrow
redshift slices and separately in NGC/SGC, all BGS galaxies are ranked by
 dereddened r-band luminosity and assigned to five weighted quantile tracers.
The five tracers are compressed into one antisymmetric luminosity-rank mark,
so every selected galaxy contributes while the data-vector dimension remains
small enough for a stable jackknife covariance.

Pair counting uses an unbiased Monte-Carlo compression of the full pair set.
For every anchor object at most K eligible neighbours are sampled uniformly;
each sampled pair carries the inverse sampling probability.  Geometry is
frozen before permutations, so null tests cannot tune the pair selection.
This run remains a matched-filter validation.  The luminosity rank is a
mass/bias proxy, not an absolute halo-mass estimate.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import chi2

import desi_dr1_phase7_lss as p
import desi_dr1_phase7_zresolved as zr


def weighted_quantile_labels(mag, weight, ntracer):
    order=np.argsort(mag)
    lab=np.empty(len(mag),np.int16)
    cw=np.cumsum(weight[order]); total=float(cw[-1])
    # Small magnitude = brighter = larger proxy mark.
    q=np.minimum((ntracer*(cw-0.5*weight[order])/max(total,1e-300)).astype(int),ntracer-1)
    lab[order]=(ntracer-1-q).astype(np.int16)
    return lab


def make_data_proxy(cat,zmin,zmax,dz,ntracer):
    good=(np.isfinite(cat['ra'])&np.isfinite(cat['dec'])&np.isfinite(cat['z'])&
          np.isfinite(cat['w'])&(cat['w']>0)&np.isfinite(cat['fluxr'])&(cat['fluxr']>0)&
          (cat['z']>=zmin)&(cat['z']<zmax))
    idx=np.where(good)[0]
    mag=22.5-2.5*np.log10(cat['fluxr'][idx]); z=cat['z'][idx]; w=cat['w'][idx]
    labels=np.full(len(idx),-1,np.int16); strata=np.full(len(idx),-1,np.int32)
    info={}; sid=0; edges=np.arange(zmin,zmax+0.5*dz,dz)
    marks=np.linspace(-1.0,1.0,ntracer)
    for cap in sorted(set(cat['region'][idx])):
        for lo,hi in zip(edges[:-1],edges[1:]):
            q=np.where((cat['region'][idx]==cap)&(z>=lo)&(z<hi))[0]
            if len(q)<max(20,2*ntracer): continue
            l=weighted_quantile_labels(mag[q],w[q],ntracer)
            labels[q]=l; strata[q]=sid
            fr=[float(w[q][l==j].sum()/w[q].sum()) for j in range(ntracer)]
            med=[float(np.median(mag[q][l==j])) if np.any(l==j) else None for j in range(ntracer)]
            info[str(sid)]={'cap':str(cap),'zlo':float(lo),'zhi':float(hi),'n':int(len(q)),
                            'weighted_fractions':fr,'median_rmag_by_tracer':med}
            sid+=1
    keep=strata>=0; idx=idx[keep]; labels=labels[keep]; strata=strata[keep]
    mark=marks[labels]
    return idx,labels,mark,strata,info


def make_random_proxy(cat,zmin,zmax,dz,ntracer,data_info,random_factor,rng):
    good=(np.isfinite(cat['ra'])&np.isfinite(cat['dec'])&np.isfinite(cat['z'])&
          np.isfinite(cat['w'])&(cat['w']>0)&(cat['z']>=zmin)&(cat['z']<zmax))
    raw=np.where(good)[0]; edges=np.arange(zmin,zmax+0.5*dz,dz)
    marks=np.linspace(-1.0,1.0,ntracer)
    out_i=[]; out_l=[]; out_s=[]
    for sid_s,meta in data_info.items():
        sid=int(sid_s); cap=meta['cap']; lo=float(meta['zlo']); hi=float(meta['zhi'])
        q=raw[(cat['region'][raw]==cap)&(cat['z'][raw]>=lo)&(cat['z'][raw]<hi)]
        if len(q)==0: continue
        target=min(len(q),max(ntracer*10,int(math.ceil(random_factor*meta['n']))))
        pick=rng.choice(q,size=target,replace=False) if target<len(q) else q
        probs=np.asarray(meta['weighted_fractions'],float); probs/=probs.sum()
        lab=rng.choice(np.arange(ntracer),size=len(pick),p=probs)
        out_i.append(np.asarray(pick,int)); out_l.append(np.asarray(lab,np.int16)); out_s.append(np.full(len(pick),sid,np.int32))
    idx=np.concatenate(out_i); labels=np.concatenate(out_l); strata=np.concatenate(out_s)
    return idx,labels,marks[labels],strata


def prepare(cat,idx,mark,strata,reg):
    x=p.distance_xyz(cat['ra'][idx],cat['dec'][idx],cat['z'][idx])
    return {'x':x,'u':x/np.maximum(np.linalg.norm(x,axis=1)[:,None],1e-300),
            'w':np.asarray(cat['w'][idx],float),'z':np.asarray(cat['z'][idx],float),
            'ra':np.asarray(cat['ra'][idx],float),'region_name':np.asarray(cat['region'][idx]),
            'mark':np.asarray(mark,float),'strata':np.asarray(strata,int),'jk':np.asarray(reg,int)}


def sample_pairs(A,B,edges,kmax,rng,theta_min_deg,auto=False,chunk=512,label='PAIR'):
    tree=cKDTree(B['x']); rmax=float(edges[-1]); rmin=float(edges[0])
    aa=[]; bb=[]; bins=[]; mus=[]; imps=[]
    n=len(A['x'])
    for start in range(0,n,chunk):
        stop=min(start+chunk,n)
        neigh=tree.query_ball_point(A['x'][start:stop],rmax,workers=-1,return_sorted=False)
        for loc,js0 in enumerate(neigh):
            i=start+loc
            if not js0: continue
            js=np.asarray(js0,int)
            if auto: js=js[js>i]
            if len(js)==0: continue
            dv=B['x'][js]-A['x'][i]
            rr=np.linalg.norm(dv,axis=1)
            keep=(rr>=rmin)&(rr<rmax)
            if theta_min_deg>0:
                c=np.sum(A['u'][i]*B['u'][js],axis=1)
                th=np.rad2deg(np.arccos(np.clip(c,-1,1)))
                keep&=(th>=theta_min_deg)
            js=js[keep]; rr=rr[keep]; dv=dv[keep]
            if len(js)==0: continue
            nall=len(js)
            if nall>kmax:
                take=rng.choice(nall,size=kmax,replace=False)
                js=js[take]; rr=rr[take]; dv=dv[take]
            ns=len(js); imp=float(nall/ns)
            mid=B['x'][js]+A['x'][i]
            mu=np.sum(dv*mid,axis=1)/(np.maximum(rr,1e-12)*np.maximum(np.linalg.norm(mid,axis=1),1e-12))
            b=np.searchsorted(edges,rr,side='right')-1
            aa.append(np.full(ns,i,np.int32)); bb.append(js.astype(np.int32)); bins.append(b.astype(np.int16)); mus.append(mu.astype(np.float32)); imps.append(np.full(ns,imp,np.float32))
        if start and start%(chunk*20)==0: print(label,'anchors',start,'/',n)
    if not aa: raise RuntimeError(f'no sampled pairs for {label}')
    out={'a':np.concatenate(aa),'b':np.concatenate(bb),'bin':np.concatenate(bins),
         'mu':np.concatenate(mus).astype(float),'imp':np.concatenate(imps).astype(float)}
    print(label,'sampled_pairs',len(out['a']),'estimated_pairs',float(out['imp'].sum()))
    return out


def hist_pair(tab,A,B,markA,markB,nb,drop=None):
    keep=np.ones(len(tab['a']),bool)
    if drop is not None:
        keep&=(A['jk'][tab['a']]!=drop)&(B['jk'][tab['b']]!=drop)
    a=tab['a'][keep]; b=tab['b'][keep]; ib=tab['bin'][keep]
    pw=A['w'][a]*B['w'][b]*tab['imp'][keep]
    h0=np.bincount(ib,weights=pw,minlength=nb).astype(float)
    odd=pw*(markB[b]-markA[a])*(3.0*tab['mu'][keep])
    h1=np.bincount(ib,weights=odd,minlength=nb).astype(float)
    return h0,h1


def object_norm(X,drop=None):
    keep=np.ones(len(X['w']),bool) if drop is None else X['jk']!=drop
    w=X['w'][keep]; W=float(w.sum())
    return W,0.5*max(W*W-float(np.sum(w*w)),1e-300)


def estimate_block(block,markD=None,drop=None):
    D,R=block['D'],block['R']; markD=D['mark'] if markD is None else markD; markR=R['mark']
    nb=block['nb']; dd0,dd1=hist_pair(block['DD'],D,D,markD,markD,nb,drop)
    dr0,dr1=hist_pair(block['DR'],D,R,markD,markR,nb,drop)
    rr0,rr1=hist_pair(block['RR'],R,R,markR,markR,nb,drop)
    WD,NDD=object_norm(D,drop); WR,NRR=object_norm(R,drop); NDR=max(WD*WR,1e-300)
    dd0/=NDD; dd1/=NDD; dr0/=NDR; dr1/=NDR; rr0/=NRR; rr1/=NRR
    den=np.maximum(rr0,1e-300)
    return (dd0-2.0*dr0+rr0)/den,(dd1-2.0*dr1+rr1)/den


def vector(blocks,marks=None,drop=None):
    x0=[]; x1=[]
    for iz,b in enumerate(blocks):
        md=None if marks is None else marks[iz]
        a,c=estimate_block(b,md,drop)
        x0.extend(a.tolist()); x1.extend(c.tolist())
    return np.asarray(x0),np.asarray(x1)


def permuted_marks(blocks,rng):
    out=[]
    for b in blocks:
        D=b['D']; m=D['mark'].copy()
        for st in np.unique(D['strata']):
            q=np.where(D['strata']==st)[0]; m[q]=rng.permutation(m[q])
        out.append(m)
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',nargs='+',required=True); ap.add_argument('--random',nargs='+',required=True)
    ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40)
    ap.add_argument('--dz-proxy',type=float,default=0.02); ap.add_argument('--ntracer',type=int,default=5)
    ap.add_argument('--analysis-z-edges',default='0.10,0.20,0.30,0.40')
    ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    ap.add_argument('--random-factor',type=float,default=2.0); ap.add_argument('--neighbors-per-anchor',type=int,default=48)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--jackknife',type=int,default=30)
    ap.add_argument('--permutations',type=int,default=32); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:]); nb=len(s)
    zedges=np.asarray([float(x) for x in args.analysis_z_edges.split(',')]); nz=len(zedges)-1
    if args.ntracer<3: raise ValueError('ntracer must be >=3')

    catD=p.read_catalog(args.data,'data'); catR=p.read_catalog(args.random,'random')
    iD,lD,mD,sD,proxy_info=make_data_proxy(catD,args.zmin,args.zmax,args.dz_proxy,args.ntracer)
    iR,lR,mR,sR=make_random_proxy(catR,args.zmin,args.zmax,args.dz_proxy,args.ntracer,proxy_info,args.random_factor,rng)
    # Full data sample is retained. Randoms are a predeclared density-controlled subset.
    regD,njd=p.sky_jackknife_regions(catD['ra'][iD],catD['region'][iD],args.jackknife)
    regR,njr=p.sky_jackknife_regions(catR['ra'][iR],catR['region'][iR],args.jackknife)
    nj=min(njd,njr)
    if nj<=nz*nb+2: raise RuntimeError('insufficient jackknife regions')

    blocks=[]; zmeta=[]
    for iz,(lo,hi) in enumerate(zip(zedges[:-1],zedges[1:])):
        qd=np.where((catD['z'][iD]>=lo)&(catD['z'][iD]<hi))[0]
        qr=np.where((catR['z'][iR]>=lo)&(catR['z'][iR]<hi))[0]
        D=prepare(catD,iD[qd],mD[qd],sD[qd],regD[qd]); R=prepare(catR,iR[qr],mR[qr],sR[qr],regR[qr])
        if min(len(D['w']),len(R['w']))<100: raise RuntimeError(f'too few objects in z bin {lo}-{hi}')
        DD=sample_pairs(D,D,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'DD_z{iz}')
        DR=sample_pairs(D,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=False,label=f'DR_z{iz}')
        RR=sample_pairs(R,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'RR_z{iz}')
        blocks.append({'D':D,'R':R,'DD':DD,'DR':DR,'RR':RR,'nb':nb})
        ze=float(np.average(D['z'],weights=D['w']))
        zmeta.append({'zlo':float(lo),'zhi':float(hi),'z_effective':ze,'N_data':int(len(D['w'])),'N_random':int(len(R['w']))})

    xi0,xi1=vector(blocks); pdim=len(xi1)
    jk=[]
    for j in range(nj):
        _,v=vector(blocks,drop=j); jk.append(v); print('JK_MT',j,json.dumps(v.tolist()))
    jk=np.asarray(jk); jm=jk.mean(axis=0); cov_raw=(nj-1)/nj*((jk-jm).T@(jk-jm))
    ev=np.linalg.eigvalsh(cov_raw); maxev=float(max(ev[-1],1e-300)); ridge=max(maxev*1e-7,float(np.median(np.diag(cov_raw)))*1e-6,1e-14)
    cov=cov_raw+np.eye(pdim)*ridge; cinv=np.linalg.pinv(cov,rcond=1e-10)
    rank=int(np.linalg.matrix_rank(cov_raw,tol=maxev*1e-10)); cond=float(np.linalg.cond(cov))

    null=[]
    for ip in range(args.permutations):
        pm=permuted_marks(blocks,rng); _,v=vector(blocks,pm); null.append(v); print('PERM_MT',ip,json.dumps(v.tolist()))
    null=np.asarray(null); null_mean=null.mean(axis=0); y=xi1-null_mean
    tdata=float(y@cinv@y); tperm=np.asarray([float((v-null_mean)@cinv@(v-null_mean)) for v in null])
    pemp=float((1+np.sum(tperm>=tdata))/(len(tperm)+1))

    wake,dop=zr.template_vector(args.templates,s,zmeta)
    fit2=p.gls_fit(y,cov,[dop,wake],['doppler','wake_matched_filter'])
    ctemps,cnames=zr.conservative_templates(s,zmeta,dop,wake)
    fitc=p.gls_fit(y,cov,ctemps,cnames)
    cosine=float(zr.metric_cosine(wake,dop,cinv))
    perm_aw=[]
    for v in null:
        f=p.gls_fit(v-null_mean,cov,ctemps,cnames); perm_aw.append(float(f['wake_matched_filter']['amplitude']))
    perm_aw=np.asarray(perm_aw); aw=float(fitc['wake_matched_filter']['amplitude'])
    p_aw=float((1+np.sum(np.abs(perm_aw)>=abs(aw)))/(len(perm_aw)+1))

    sig=np.sqrt(np.maximum(np.diag(cov),0)); rows=[]; k=0
    for m in zmeta:
        for si in s:
            rows.append((m['zlo'],m['zhi'],m['z_effective'],si,xi0[k],xi1[k],null_mean[k],y[k],sig[k])); k+=1
    np.savetxt(out/'data_vector_multitracer.csv',np.asarray(rows),delimiter=',',
               header='zlo,zhi,z_effective,s_Mpc_over_h,xi0_proxy,xi1_proxy_odd,permutation_null_mean,xi1_null_corrected,jackknife_sigma',comments='')
    np.savetxt(out/'jackknife_covariance_multitracer.csv',cov,delimiter=',')
    np.savetxt(out/'permutation_vectors_multitracer.csv',null,delimiter=',')

    summary={
      'scope':'Full-data DESI DR1 BGS five-tracer luminosity-rank proxy odd-sector matched-filter test; not an absolute halo-mass or DESI wake constraint.',
      'proxy_definition':'Weighted dereddened-r luminosity quantiles in narrow redshift bins and survey cap; brighter tracers receive larger proxy mark.',
      'predeclared_ntracer':int(args.ntracer),'uses_all_selected_data_galaxies':True,
      'data_count':int(len(iD)),'random_count':int(len(iR)),'random_factor_target':args.random_factor,
      'analysis_z_edges':zedges.tolist(),'z_bins':zmeta,'separation_edges_Mpc_over_h':edges.tolist(),
      'pair_sampling':{'neighbors_per_anchor':int(args.neighbors_per_anchor),'rule':'uniform eligible-neighbour sampling per anchor with inverse-probability pair weight; geometry frozen before null permutations'},
      'vector_dimension':pdim,'jackknife':{'regions':int(nj),'raw_covariance_rank':rank,'regularized_condition_number':cond,'ridge_added':ridge},
      'permutation_control':{'count':int(len(null)),'mahalanobis_data':tdata,'empirical_pvalue':pemp,'conservative_wake_empirical_two_sided_pvalue':p_aw},
      'template_metric_wake_doppler_cosine':cosine,
      'fit_null_corrected_doppler_plus_wake':fit2,
      'fit_null_corrected_conservative_per_z_odd_plus_wake':fitc,
      'proxy_strata':proxy_info,
      'analysis_scope':'We use the specified tracer ranks, redshift and separation bins and pair-sampling rule. Independent pair-sampling seeds and survey mocks provide validation.'
    }
    (out/'summary_multitracer.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_MULTITRACER_FULLSAMPLE',json.dumps(summary,indent=2))

if __name__=='__main__': main()
