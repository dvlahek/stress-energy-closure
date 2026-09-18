#!/usr/bin/env python3
"""Secondary DESI DR1 genuine cross-population boundary diagnostic.

LRG 0.40<=z<0.50 -> BGS_BRIGHT-21.5 0.30<=z<0.40.
The two selections touch at z=0.4, so this is a boundary diagnostic only.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.stats import chi2

import desi_dr1_genuine_multitracer as gm


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--lrg-data',nargs='+',required=True)
    ap.add_argument('--lrg-random',nargs='+',required=True)
    ap.add_argument('--bgs-data',nargs='+',required=True)
    ap.add_argument('--bgs-random',nargs='+',required=True)
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

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')])
    s=0.5*(edges[:-1]+edges[1:]); nb=len(s)

    # Orientation fixed before inspection: A=LRG, B=BGS.
    DA=gm.read_catalog(args.lrg_data,0.40,0.50)
    DB=gm.read_catalog(args.bgs_data,0.30,0.40)
    RA0=gm.read_catalog(args.lrg_random,0.40,0.50)
    RB0=gm.read_catalog(args.bgs_random,0.30,0.40)
    RA0=gm.stratified_random_subset(RA0,DA,0.40,0.50,args.dz_random,args.random_factor,rng)
    RB0=gm.stratified_random_subset(RB0,DB,0.30,0.40,args.dz_random,args.random_factor,rng)

    jkA,jkB,jkRA,jkRB,nj,bounds=gm.common_jackknife(DA,DB,RA0,RB0,args.jackknife)
    if nj<=nb+2:
        raise RuntimeError(f'insufficient common jackknife regions: {nj}')

    A=gm.prepare(DA,jkA); B=gm.prepare(DB,jkB)
    RA=gm.prepare(RA0,jkRA); RB=gm.prepare(RB0,jkRB)

    DD=gm.sample_pairs(A,B,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label='LRGxBGS_DD')
    DR=gm.sample_pairs(A,RB,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label='LRGxBGS_DR')
    RD=gm.sample_pairs(RA,B,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label='LRGxBGS_RD')
    RR=gm.sample_pairs(RA,RB,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,label='LRGxBGS_RR')
    block={'A':A,'B':B,'RA':RA,'RB':RB,'DD':DD,'DR':DR,'RD':RD,'RR':RR,'nb':nb}

    dip=gm.estimate(block,1)
    octu=gm.estimate(block,3)
    jk_d=[]; jk_o=[]
    for j in range(nj):
        jk_d.append(gm.estimate(block,1,j))
        jk_o.append(gm.estimate(block,3,j))
        print('JK_SECONDARY',j)
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

    za=float(np.average(A['z'],weights=A['w']))
    zb=float(np.average(B['z'],weights=B['w']))
    zeff=0.5*(za+zb)
    wake,dop=gm.load_template(args.templates,s,zeff)
    b1=(s/60.0)**-1; b2=(s/60.0)**-2
    fit=gm.gls_fit(
        dip,cov,
        [dop,b1,b2,wake],
        ['doppler','broadband_inv_s','broadband_inv_s2','wake_matched_filter']
    )

    t0=float(dip@ci@dip); p0=float(chi2.sf(t0,len(dip)))
    to=float(octu@cio@octu); po=float(chi2.sf(to,len(octu)))

    sig=np.sqrt(np.maximum(np.diag(cov),0))
    sigo=np.sqrt(np.maximum(np.diag(covo),0))
    np.savetxt(out/'data_vector.csv',np.column_stack([s,dip,sig,octu,sigo]),delimiter=',',
               header='s_Mpc_over_h,xi1_LRG_to_BGS,jk_sigma_dipole,xi3_control,jk_sigma_octupole',comments='')
    np.savetxt(out/'jackknife_covariance_dipole.csv',cov,delimiter=',')
    np.savetxt(out/'jackknife_covariance_octupole.csv',covo,delimiter=',')
    np.savetxt(out/'jackknife_vectors_dipole.csv',jk_d,delimiter=',')
    np.savetxt(out/'jackknife_vectors_octupole.csv',jk_o,delimiter=',')

    summary={
        'status':'SECONDARY_BOUNDARY_DIAGNOSTIC',
        'pair':'LRG -> BGS_BRIGHT-21.5',
        'release':'DESI DR1 LSS v1.5 clustering catalogs',
        'selection':{'LRG_z':[0.40,0.50],'BGS_z':[0.30,0.40]},
        'z_effective_midpoint_of_weighted_means':zeff,
        'counts':{'LRG':len(A['z']),'BGS':len(B['z']),'R_LRG':len(RA['z']),'R_BGS':len(RB['z'])},
        'separation_edges_Mpc_over_h':edges.tolist(),
        'theta_min_deg':args.theta_min_deg,
        'neighbors_per_anchor':args.neighbors_per_anchor,
        'random_factor':args.random_factor,
        'seed':args.seed,
        'jackknife':{'regions':nj,'raw_rank':rank,'ridge':ridge,'condition':cond},
        'zero_vector':{'chi2':t0,'dof':len(dip),'pvalue_asymptotic':p0},
        'conservative_fit':fit,
        'odd_octupole_control':{'chi2':to,'dof':len(octu),'pvalue_asymptotic':po,
                                'jackknife_raw_rank':ranko,'condition':condo},
        'interpretation_guardrail':'Boundary diagnostic only because BGS and LRG clustering selections meet at z=0.4 rather than sharing a broad volume.'
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('GENUINE_SECONDARY',json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
