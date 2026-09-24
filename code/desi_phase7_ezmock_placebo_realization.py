#!/usr/bin/env python3
"""One EZmock realization for the Phase-7 five-tracer placebo covariance control.

The public DR1 EZmock BGS products currently exposed in the release do not
contain R_MAG_APP/R_MAG_ABS. For the EZmock-only covariance/systematics
cross-check we therefore construct five equal-count tracer marks from the
released uniform RAN_NUM_0_1 field, ranked independently in the same
narrow-redshift and NGC/SGC cells used by the real-data luminosity proxy. This
deliberately removes luminosity-dependent bias and is not a replacement for
the Abacus luminosity-ranked validation.

Each realization must use its own released NGC/SGC clustering random catalogs.
Those files are part of the realization-specific survey-selection and
fiber-assignment product. Reusing a random catalog from another mock is not a
valid homogeneous covariance construction and is deliberately unsupported.
"""
from __future__ import annotations
import argparse, json
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
                flux=u
            else:
                flux=np.full(len(d),np.nan)
            chunks.append({'ra':ra,'dec':dec,'z':z,'w':w*wf,'fluxr':flux,
                           'region':np.full(len(d),cap,dtype='U3')})
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}


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
    ap.add_argument('--random',nargs='+',required=True)
    ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True); ap.add_argument('--mock-id',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40)
    ap.add_argument('--dz-proxy',type=float,default=0.02); ap.add_argument('--ntracer',type=int,default=5)
    ap.add_argument('--analysis-z-edges',default='0.10,0.20,0.30,0.40')
    ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    ap.add_argument('--random-factor',type=float,default=2.0); ap.add_argument('--neighbors-per-anchor',type=int,default=48)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    mid=int(args.mock_id); rng=np.random.default_rng(args.seed+1009*mid)
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:]); nb=len(s)
    zedges=np.asarray([float(x) for x in args.analysis_z_edges.split(',')])

    catD=read_mock(args.data,'data')
    catR=read_mock(args.random,'random')
    iD,lD,mD,sD,proxy_info=mt.make_data_proxy(catD,args.zmin,args.zmax,args.dz_proxy,args.ntracer)
    iR,lR,mR,sR=mt.make_random_proxy(catR,args.zmin,args.zmax,args.dz_proxy,args.ntracer,proxy_info,args.random_factor,rng)
    random_geometry_mode='realization-specific released EZmock clustering random catalogs'
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
    max_xi0=float(np.max(np.abs(xi0))); max_xi1=float(np.max(np.abs(xi1)))
    if max_xi0>0.30 or max_xi1>0.20:
        raise RuntimeError(
            f'EZmock estimator sanity gate failed: max|xi0|={max_xi0:.6g}, '
            f'max|xi1|={max_xi1:.6g}; expected placebo-scale values. '
            'Check realization-specific random survey selection/normalization before continuing.')

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
      'sanity':{'max_abs_xi0':max_xi0,'max_abs_xi1':max_xi1,'gate_max_abs_xi0':0.30,'gate_max_abs_xi1':0.20},
      'window_model':'pair-level response to xi_odd^{ij}=(m_j-m_i) mu T(s,z), evaluated through the identical sampled DD geometry and RR normalization',
      'wake_forward_response':wake.tolist(),'doppler_forward_response':dop.tolist(),'xi1_proxy_odd':xi1.tolist(),
      'analysis_scope':'The available EZmock BGS catalogues do not provide the luminosity field required by the five-tracer measurement. We assign equal-count random ranks for a survey-geometry and covariance control, using the random catalogues supplied for each realization.'
    }
    (out/f'{stem}_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_EZMOCK_PLACEBO_REALIZATION',json.dumps(summary))

if __name__=='__main__': main()
