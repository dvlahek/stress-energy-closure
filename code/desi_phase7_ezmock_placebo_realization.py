#!/usr/bin/env python3
"""One EZmock realization for the Phase-7 five-tracer placebo covariance control.

The public DR1 EZmock BGS products currently exposed in the release do not
contain R_MAG_APP/R_MAG_ABS, despite the current data-model documentation.
For the EZmock-only covariance/systematics cross-check we therefore construct
five equal-count tracer marks from the released uniform RAN_NUM_0_1 field,
ranked independently in the same narrow-redshift and NGC/SGC cells used by the
real-data luminosity proxy.  This deliberately removes luminosity-dependent
bias and is *not* a replacement for the Abacus luminosity-ranked validation.
It tests survey geometry, pair compression, estimator window response and null
false-positive behaviour for the identical marked-dipole machinery.

For local ensemble runs, a compact shared angular-random cache can be used.
Only RA/DEC and random weights are shared.  Redshifts are reassigned for every
mock from that realization's own data in the same narrow-z and Galactic-cap
cells before the usual density-controlled random selection.  This avoids
re-downloading gigabyte random FITS files without freezing the radial
selection to one realization.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from astropy.io import fits

import desi_dr1_phase7_lss as p
import desi_dr1_phase7_multitracer_fullsample as mt


def _col(d, *names, default=None):
    avail={n.upper():n for n in d.names}
    for n in names:
        if n.upper() in avail:
            return np.asarray(d[avail[n.upper()]])
    if default is not None:
        return default
    raise KeyError(f'missing {names}; available={list(d.names)}')


def read_mock(paths, kind):
    chunks=[]
    for path in paths:
        name=Path(path).name.upper()
        cap='NGC' if 'NGC' in name else 'SGC' if 'SGC' in name else 'UNK'
        with fits.open(path,memmap=True) as hdul:
            hdul.verify('exception')
            d=hdul[1].data
            ra=_col(d,'RA').astype(float); dec=_col(d,'DEC').astype(float); z=_col(d,'Z','RSDZ').astype(float)
            w=_col(d,'WEIGHT',default=np.ones(len(d))).astype(float)
            wf=_col(d,'WEIGHT_FKP',default=np.ones(len(d))).astype(float)
            if kind=='data':
                u=_col(d,'RAN_NUM_0_1').astype(float)
                if not np.all(np.isfinite(u)):
                    raise RuntimeError(f'{path}: non-finite RAN_NUM_0_1')
                # make_data_proxy only requires a monotonic ranking variable.
                # The released uniform field gives an explicit equal-count placebo rank.
                flux=u
            else:
                flux=np.full(len(d),np.nan)
            chunks.append({'ra':ra,'dec':dec,'z':z,'w':w*wf,'fluxr':flux,
                           'region':np.full(len(d),cap,dtype='U3')})
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}


def read_shared_random_cache(path):
    with np.load(path, allow_pickle=False) as a:
        ra=np.asarray(a['ra'],float); dec=np.asarray(a['dec'],float); w=np.asarray(a['w'],float)
        cap=np.asarray(a['cap'],np.int8)
    if not (len(ra)==len(dec)==len(w)==len(cap)):
        raise RuntimeError('shared random cache arrays have inconsistent lengths')
    region=np.where(cap==0,'NGC','SGC').astype('U3')
    return {'ra':ra,'dec':dec,'z':np.full(len(ra),np.nan),'w':w,
            'fluxr':np.full(len(ra),np.nan),'region':region}


def remap_shared_random_redshifts(catR,catD,proxy_info,random_factor,ntracer,rng):
    """Assign per-realization z values to a shared angular random cache.

    Each narrow-z/cap stratum receives exactly the random density requested by
    the estimator.  Angular positions are used without replacement within a
    cap; redshifts are sampled from the current mock data in the same stratum.
    """
    selected=[]; zdraw=[]
    by_cap={}
    for cap in ('NGC','SGC'):
        pool=np.where(catR['region']==cap)[0]
        pool=np.asarray(pool,int).copy(); rng.shuffle(pool); by_cap[cap]=[pool,0]

    metas=sorted(proxy_info.items(),key=lambda kv:(kv[1]['cap'],float(kv[1]['zlo'])))
    for sid_s,meta in metas:
        cap=str(meta['cap']); lo=float(meta['zlo']); hi=float(meta['zhi'])
        target=max(ntracer*10,int(math.ceil(random_factor*int(meta['n']))))
        pool,cursor=by_cap[cap]
        if cursor+target>len(pool):
            raise RuntimeError(
                f'shared random cache too small for {cap}: need at least {cursor+target}, have {len(pool)}; '
                'rebuild cache with larger per-cap counts')
        pick=pool[cursor:cursor+target]; by_cap[cap][1]=cursor+target
        q=np.where((catD['region']==cap)&np.isfinite(catD['z'])&(catD['z']>=lo)&(catD['z']<hi))[0]
        if len(q)==0:
            raise RuntimeError(f'no data redshifts available for shared-random remap {cap} {lo}-{hi}')
        selected.append(pick)
        zdraw.append(rng.choice(catD['z'][q],size=target,replace=True))

    idx=np.concatenate(selected); z=np.concatenate(zdraw)
    out={k:np.asarray(v[idx]).copy() for k,v in catR.items()}
    out['z']=np.asarray(z,float)
    return out


def forward_block(block,template_path,z_eff):
    D,R=block['D'],block['R']; tab=block['DD']; nb=block['nb']
    a=tab['a']; b=tab['b']; ib=tab['bin']; mu=tab['mu']
    sep=np.linalg.norm(D['x'][b]-D['x'][a],axis=1)
    wake,dop=p.load_templates(template_path,sep,float(z_eff))
    dm=D['mark'][b]-D['mark'][a]
    pw=D['w'][a]*D['w'][b]*tab['imp']
    common=pw*(3.0*mu*mu*dm*dm)
    hw=np.bincount(ib,weights=common*wake,minlength=nb).astype(float)
    hd=np.bincount(ib,weights=common*dop,minlength=nb).astype(float)
    _,NDD=mt.object_norm(D)
    rr0,_=mt.hist_pair(block['RR'],R,R,R['mark'],R['mark'],nb)
    _,NRR=mt.object_norm(R); rr0=rr0/NRR
    den=np.maximum(rr0,1e-300)
    return (hw/NDD)/den,(hd/NDD)/den


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',nargs='+',required=True)
    rg=ap.add_mutually_exclusive_group(required=True)
    rg.add_argument('--random',nargs='+')
    rg.add_argument('--shared-random-cache')
    ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True); ap.add_argument('--mock-id',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40)
    ap.add_argument('--dz-proxy',type=float,default=0.02); ap.add_argument('--ntracer',type=int,default=5)
    ap.add_argument('--analysis-z-edges',default='0.10,0.20,0.30,0.40')
    ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    ap.add_argument('--random-factor',type=float,default=1.0); ap.add_argument('--neighbors-per-anchor',type=int,default=48)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    mid=int(args.mock_id); rng=np.random.default_rng(args.seed+1009*mid)
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:]); nb=len(s)
    zedges=np.asarray([float(x) for x in args.analysis_z_edges.split(',')])

    catD=read_mock(args.data,'data')
    iD,lD,mD,sD,proxy_info=mt.make_data_proxy(catD,args.zmin,args.zmax,args.dz_proxy,args.ntracer)

    if args.shared_random_cache:
        catR=read_shared_random_cache(args.shared_random_cache)
        catR=remap_shared_random_redshifts(catR,catD,proxy_info,args.random_factor,args.ntracer,rng)
        random_geometry_mode='shared angular cache with per-realization narrow-z/cap redshift resampling from current mock data'
    else:
        catR=read_mock(args.random,'random')
        random_geometry_mode='realization-specific released EZmock clustering random catalogs'

    iR,lR,mR,sR=mt.make_random_proxy(catR,args.zmin,args.zmax,args.dz_proxy,args.ntracer,proxy_info,args.random_factor,rng)
    regD=np.zeros(len(iD),int); regR=np.zeros(len(iR),int)

    blocks=[]; zmeta=[]
    for iz,(lo,hi) in enumerate(zip(zedges[:-1],zedges[1:])):
        qd=np.where((catD['z'][iD]>=lo)&(catD['z'][iD]<hi))[0]
        qr=np.where((catR['z'][iR]>=lo)&(catR['z'][iR]<hi))[0]
        D=mt.prepare(catD,iD[qd],mD[qd],sD[qd],regD[qd]); R=mt.prepare(catR,iR[qr],mR[qr],sR[qr],regR[qr])
        if min(len(D['w']),len(R['w']))<100: raise RuntimeError(f'too few mock objects in z={lo}-{hi}')
        DD=mt.sample_pairs(D,D,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'EZ{mid}_DD_z{iz}')
        DR=mt.sample_pairs(D,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=False,label=f'EZ{mid}_DR_z{iz}')
        RR=mt.sample_pairs(R,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'EZ{mid}_RR_z{iz}')
        block={'D':D,'R':R,'DD':DD,'DR':DR,'RR':RR,'nb':nb}; blocks.append(block)
        ze=float(np.average(D['z'],weights=D['w']))
        zmeta.append({'zlo':float(lo),'zhi':float(hi),'z_effective':ze,'N_data':int(len(D['w'])),'N_random':int(len(R['w']))})

    xi0,xi1=mt.vector(blocks)
    wake=[]; dop=[]
    for b,m in zip(blocks,zmeta):
        w,d=forward_block(b,args.templates,m['z_effective']); wake.extend(w.tolist()); dop.extend(d.tolist())
    wake=np.asarray(wake); dop=np.asarray(dop)

    stem=f'mock_{mid:02d}'
    rows=[]; k=0
    for m in zmeta:
        for si in s:
            rows.append((m['zlo'],m['zhi'],m['z_effective'],si,xi0[k],xi1[k])); k+=1
    np.savetxt(out/f'{stem}_vector.csv',np.asarray(rows),delimiter=',',header='zlo,zhi,z_effective,s_Mpc_over_h,xi0_proxy,xi1_proxy_odd',comments='')
    np.savetxt(out/f'{stem}_window.csv',np.column_stack([wake,dop]),delimiter=',',header='wake_forward_response,doppler_forward_response',comments='')
    summary={
      'mock_id':mid,
      'scope':'cut-sky BGS EZmock five-tracer equal-count random-rank placebo covariance/systematics validation',
      'rank_field':'RAN_NUM_0_1','rank_semantics':'uniform released mock field; re-ranked within the same narrow-z and Galactic-cap cells as the real luminosity proxy',
      'data_count':int(len(iD)),'random_count':int(len(iR)),'ntracer':int(args.ntracer),'z_bins':zmeta,
      'separation_edges_Mpc_over_h':edges.tolist(),'neighbors_per_anchor':int(args.neighbors_per_anchor),'theta_min_deg':args.theta_min_deg,
      'seed':int(args.seed+1009*mid),'random_geometry_mode':random_geometry_mode,
      'window_model':'pair-level response to xi_odd^{ij}=(m_j-m_i) mu T(s,z), evaluated through the identical sampled DD geometry and RR normalization',
      'wake_forward_response':wake.tolist(),'doppler_forward_response':dop.tolist(),'xi1_proxy_odd':xi1.tolist(),
      'guardrail':'EZmock public DR1 files used here lack released luminosity columns. This is an equal-count random-rank placebo covariance/systematics control, not a luminosity-matched covariance and not a halo-mass wake constraint. Shared-random mode reuses only angular survey geometry and redraws redshifts per realization. Abacus remains the physical luminosity-ranked mock validation.'
    }
    (out/f'{stem}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_EZMOCK_PLACEBO_REALIZATION',json.dumps(summary))

if __name__=='__main__': main()
