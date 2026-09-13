#!/usr/bin/env python3
"""DESI DR1 Phase-7 five-tracer halo-mass proxy using the public Gfinder VAC.

The primary LSS BGS TARGETID is decoded to the Legacy Surveys imaging tuple
(RELEASE, BRICKID, OBJID), which is matched exactly to the Gfinder galaxy VAC.
The published galaxy-group relation and group GRP_LOGM then attach a halo-mass
proxy to each matched BGS object.  Five weighted mass quantiles are fixed a
priori in the same narrow-z and NGC/SGC strata used by the luminosity analysis.

This is a physical-proxy robustness test, not a re-optimization of the odd
signal.  The Gfinder selection function is not identical to the primary BGS
clustering selection, so the result is reported alongside (not in place of)
the luminosity-proxy mock-calibrated analysis.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from astropy.io import fits
from desitarget.targets import decode_targetid

import desi_dr1_phase7_lss as p
import desi_dr1_phase7_multitracer_fullsample as mt
import desi_dr1_phase7_zresolved as zr


def _col(d,*names,default=None):
    avail={n.upper():n for n in d.names}
    for n in names:
        if n.upper() in avail: return np.asarray(d[avail[n.upper()]])
    if default is not None: return default
    raise KeyError(f'missing {names}; available={list(d.names)}')


def read_lss(paths,kind):
    chunks=[]
    for path in paths:
        cap='NGC' if 'NGC' in Path(path).name.upper() else 'SGC'
        with fits.open(path,memmap=True) as hdul:
            d=hdul[1].data
            ra=_col(d,'RA').astype(float); dec=_col(d,'DEC').astype(float); z=_col(d,'Z').astype(float)
            w=_col(d,'WEIGHT',default=np.ones(len(d))).astype(float)
            wf=_col(d,'WEIGHT_FKP',default=np.ones(len(d))).astype(float)
            tid=_col(d,'TARGETID',default=np.full(len(d),-1,dtype=np.int64)).astype(np.int64)
            chunks.append({'ra':ra,'dec':dec,'z':z,'w':w*wf,'targetid':tid,
                           'fluxr':np.full(len(d),np.nan), 'group_mass':np.full(len(d),np.nan),
                           'region':np.full(len(d),cap,dtype='U3')})
    return {k:np.concatenate([c[k] for c in chunks]) for k in chunks[0]}


def imaging_key(release,brickid,objid):
    release=np.asarray(release,np.int64); brickid=np.asarray(brickid,np.int64); objid=np.asarray(objid,np.int64)
    if np.any(release<0)|np.any(brickid<0)|np.any(objid<0):
        good=(release>=0)&(brickid>=0)&(objid>=0)
    else: good=np.ones(len(release),bool)
    if np.any(brickid[good]>=1_000_000) or np.any(objid[good]>=10_000_000) or np.any(release[good]>=100_000):
        raise RuntimeError('imaging-key bounds exceeded')
    key=np.full(len(release),-1,np.int64)
    key[good]=((release[good]*1_000_000+brickid[good])*10_000_000+objid[good])
    return key


def _match_sorted(wanted_values,wanted_indices,values):
    pos=np.searchsorted(wanted_values,values)
    good=(pos<len(wanted_values))
    outpos=np.zeros(len(values),int); outpos[good]=pos[good]
    good &= wanted_values[outpos]==values
    return good,wanted_indices[outpos]


def attach_group_mass(targetid,galaxy_path,relation_path,group_path,chunk=1_000_000):
    objid,brickid,release,mock,sky,gaiadr=decode_targetid(targetid)
    bkey=imaging_key(release,brickid,objid)
    valid=bkey>=0; order=np.argsort(bkey[valid]); base_idx=np.where(valid)[0][order]; skey=bkey[base_idx]
    igal=np.full(len(targetid),-1,np.int64)

    with fits.open(galaxy_path,memmap=True) as hdul:
        d=hdul[1].data
        for st in range(0,len(d),chunk):
            q=d[st:st+chunk]
            gkey=imaging_key(_col(q,'RELEASE'),_col(q,'BRICKID'),_col(q,'OBJID'))
            ok=gkey>=0
            good,bidx=_match_sorted(skey,base_idx,gkey[ok])
            if np.any(good): igal[bidx[good]]=_col(q,'IGAL')[ok][good].astype(np.int64)
    matched_igal=np.where(igal>=0)[0]
    if len(matched_igal)==0: raise RuntimeError('no BGS/Gfinder imaging matches')

    want_igal=igal[matched_igal]; oo=np.argsort(want_igal); sigal=want_igal[oo]; bidx_igal=matched_igal[oo]
    igrp=np.full(len(targetid),-1,np.int64)
    with fits.open(relation_path,memmap=True) as hdul:
        d=hdul[1].data
        for st in range(0,len(d),chunk):
            q=d[st:st+chunk]; vals=_col(q,'IGAL').astype(np.int64)
            good,bidx=_match_sorted(sigal,bidx_igal,vals)
            if np.any(good): igrp[bidx[good]]=_col(q,'IGRP')[good].astype(np.int64)

    matched_grp=np.where(igrp>=0)[0]
    want_grp=igrp[matched_grp]; oo=np.argsort(want_grp); sgrp=want_grp[oo]; bidx_grp=matched_grp[oo]
    mass=np.full(len(targetid),np.nan,float)
    with fits.open(group_path,memmap=True) as hdul:
        d=hdul[1].data
        for st in range(0,len(d),chunk):
            q=d[st:st+chunk]; vals=_col(q,'IGRP').astype(np.int64)
            good,bidx=_match_sorted(sgrp,bidx_grp,vals)
            if np.any(good): mass[bidx[good]]=_col(q,'GRP_LOGM')[good].astype(float)
    stats={'target_count':int(len(targetid)),'imaging_match_count':int(np.sum(igal>=0)),
           'group_relation_count':int(np.sum(igrp>=0)),'mass_match_count':int(np.sum(np.isfinite(mass)))}
    return mass,stats


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--data',nargs='+',required=True); ap.add_argument('--random',nargs='+',required=True)
    ap.add_argument('--gfinder-galaxy',required=True); ap.add_argument('--gfinder-relation',required=True); ap.add_argument('--gfinder-group',required=True)
    ap.add_argument('--templates',required=True); ap.add_argument('--outdir',required=True)
    ap.add_argument('--zmin',type=float,default=0.10); ap.add_argument('--zmax',type=float,default=0.40); ap.add_argument('--dz-proxy',type=float,default=0.02)
    ap.add_argument('--analysis-z-edges',default='0.10,0.20,0.30,0.40'); ap.add_argument('--sep-edges',default='20,40,60,80,100,120,140')
    ap.add_argument('--ntracer',type=int,default=5); ap.add_argument('--random-factor',type=float,default=2.0); ap.add_argument('--neighbors-per-anchor',type=int,default=48)
    ap.add_argument('--theta-min-deg',type=float,default=0.05); ap.add_argument('--jackknife',type=int,default=30); ap.add_argument('--permutations',type=int,default=32); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(args.seed)
    edges=np.asarray([float(x) for x in args.sep_edges.split(',')]); s=0.5*(edges[:-1]+edges[1:]); nb=len(s)
    zedges=np.asarray([float(x) for x in args.analysis_z_edges.split(',')]); nz=len(zedges)-1

    catD=read_lss(args.data,'data'); catR=read_lss(args.random,'random')
    mass,match_stats=attach_group_mass(catD['targetid'],args.gfinder_galaxy,args.gfinder_relation,args.gfinder_group)
    catD['group_mass']=mass
    # Synthetic monotonic flux: mt.make_data_proxy converts it to a magnitude and only uses the rank.
    good=np.isfinite(mass); catD['fluxr'][good]=np.power(10.0,0.4*(mass[good]-10.0))
    iD,lD,mD,sD,proxy_info=mt.make_data_proxy(catD,args.zmin,args.zmax,args.dz_proxy,args.ntracer)
    iR,lR,mR,sR=mt.make_random_proxy(catR,args.zmin,args.zmax,args.dz_proxy,args.ntracer,proxy_info,args.random_factor,rng)
    regD,njd=p.sky_jackknife_regions(catD['ra'][iD],catD['region'][iD],args.jackknife)
    regR,njr=p.sky_jackknife_regions(catR['ra'][iR],catR['region'][iR],args.jackknife); nj=min(njd,njr)
    if nj<=nz*nb+2: raise RuntimeError('insufficient jackknife regions')

    blocks=[]; zmeta=[]
    for iz,(lo,hi) in enumerate(zip(zedges[:-1],zedges[1:])):
        qd=np.where((catD['z'][iD]>=lo)&(catD['z'][iD]<hi))[0]; qr=np.where((catR['z'][iR]>=lo)&(catR['z'][iR]<hi))[0]
        D=mt.prepare(catD,iD[qd],mD[qd],sD[qd],regD[qd]); R=mt.prepare(catR,iR[qr],mR[qr],sR[qr],regR[qr])
        DD=mt.sample_pairs(D,D,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'GM_DD_z{iz}')
        DR=mt.sample_pairs(D,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=False,label=f'GM_DR_z{iz}')
        RR=mt.sample_pairs(R,R,edges,args.neighbors_per_anchor,rng,args.theta_min_deg,auto=True,label=f'GM_RR_z{iz}')
        blocks.append({'D':D,'R':R,'DD':DD,'DR':DR,'RR':RR,'nb':nb})
        zmeta.append({'zlo':float(lo),'zhi':float(hi),'z_effective':float(np.average(D['z'],weights=D['w'])),'N_data':int(len(D['w'])),'N_random':int(len(R['w']))})

    xi0,xi1=mt.vector(blocks); pdim=len(xi1)
    jk=[]
    for j in range(nj): _,v=mt.vector(blocks,drop=j); jk.append(v)
    jk=np.asarray(jk); jm=jk.mean(axis=0); cov_raw=(nj-1)/nj*((jk-jm).T@(jk-jm))
    ev=np.linalg.eigvalsh(cov_raw); maxev=float(max(ev[-1],1e-300)); ridge=max(maxev*1e-7,float(np.median(np.diag(cov_raw)))*1e-6,1e-14)
    cov=cov_raw+np.eye(pdim)*ridge; cinv=np.linalg.pinv(cov,rcond=1e-10)

    null=[]
    for ip in range(args.permutations):
        pm=mt.permuted_marks(blocks,rng); _,v=mt.vector(blocks,pm); null.append(v)
    null=np.asarray(null); null_mean=null.mean(axis=0); y=xi1-null_mean
    tdata=float(y@cinv@y); tperm=np.asarray([float((v-null_mean)@cinv@(v-null_mean)) for v in null]); pemp=float((1+np.sum(tperm>=tdata))/(len(null)+1))
    wake,dop=zr.template_vector(args.templates,s,zmeta)
    fit2=p.gls_fit(y,cov,[dop,wake],['doppler','wake_matched_filter'])
    ctemps,cnames=zr.conservative_templates(s,zmeta,dop,wake); fitc=p.gls_fit(y,cov,ctemps,cnames)
    perm_aw=[]
    for v in null:
        f=p.gls_fit(v-null_mean,cov,ctemps,cnames); perm_aw.append(float(f['wake_matched_filter']['amplitude']))
    aw=float(fitc['wake_matched_filter']['amplitude']); p_aw=float((1+np.sum(np.abs(perm_aw)>=abs(aw)))/(len(null)+1))

    med=[]
    for j in range(args.ntracer):
        q=lD==j; med.append({'tracer':j,'N':int(np.sum(q)),'median_GRP_LOGM':float(np.median(mass[iD][q])),'p16_GRP_LOGM':float(np.percentile(mass[iD][q],16)),'p84_GRP_LOGM':float(np.percentile(mass[iD][q],84))})
    np.savetxt(out/'data_vector_gfinder_massproxy.csv',np.column_stack([xi0,xi1,null_mean,y,np.sqrt(np.maximum(np.diag(cov),0))]),delimiter=',',header='xi0_proxy,xi1_proxy_odd,permutation_null_mean,xi1_null_corrected,jackknife_sigma',comments='')
    summary={'scope':'DESI DR1 BGS five-tracer Gfinder GRP_LOGM physical-proxy odd-sector robustness test.',
             'predeclared_split':'five weighted GRP_LOGM quantiles independently in dz=0.02 and NGC/SGC','match_stats':match_stats,'selected_data_count':int(len(iD)),'selected_random_count':int(len(iR)),
             'mass_tracer_summary':med,'z_bins':zmeta,'jackknife':{'regions':int(nj),'rank':int(np.linalg.matrix_rank(cov_raw)),'condition_number':float(np.linalg.cond(cov))},
             'permutation':{'count':int(len(null)),'global_empirical_pvalue':pemp,'wake_empirical_two_sided_pvalue':p_aw},
             'fit_minimal':fit2,'fit_conservative':fitc,
             'guardrail':'No mass threshold was selected from the odd data. Gfinder membership introduces an additional selection function, so this is a physical-proxy robustness test rather than the primary DESI likelihood.'}
    (out/'summary_gfinder_massproxy.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_GFINDER_MASSPROXY',json.dumps(summary,indent=2))

if __name__=='__main__': main()
