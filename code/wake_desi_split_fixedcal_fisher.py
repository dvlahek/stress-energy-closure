#!/usr/bin/env python3
"""Compare DESI-BGS halo-mass splits with one fixed physical wake calibration.

Previous threshold comparisons recalibrated the absolute Okoli wake amplitude
separately for each tracer split.  That is appropriate for matching each split
to a literature scaling diagnostic, but it is not appropriate for deciding
which split is optimal: the physical wake normalization cannot depend on how
the same parent galaxy catalogue is partitioned.

Here the per-redshift absolute calibration is determined ONCE from the
log10(Msplit/[h^-1 Msun])=13.75 reference sample and then frozen.  All other
splits change only tracer counts, biases and covariance.  Hidden kinetic
states are re-optimized after full odd-nuisance projection and top candidates
are validated with full CLASS F_+ and F_- transfer functions.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_hod_masssplit_fisher as refmod
from wake_desi_split_compare_fisher import SURVEYS

PROFILES={"13.75":list(refmod.w.SURVEY), **SURVEYS}


def make_calibration(q,f0,state0,mass,refsurvey):
    """Per-z physical normalization fixed from the 13.75 reference split."""
    cal={}
    for s in refsurvey:
        z=float(s['z']); info=state0['z'][z]
        bl,bf=float(s['b_bright']),float(s['b_faint'])
        nl,nf=float(s['n_bright_h3_Mpc3']),float(s['n_faint_h3_Mpc3'])
        V=float(s['volume_hminus3_Gpc3'])
        rows=w.mode_noise_weight(info,V,bl,bf,nl,nf)
        raw_base=0.0; cache={}
        for ik,im,wt in rows:
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,f0,f0,state0,state0,state0,mass,z,ik)[0]
            raw_base += wt*cache[ik]
        target=w.desired_fd_sn2(mass,z,s['N_faint']+s['N_bright'],bl-bf)
        cal[z]=target/max(raw_base,1e-300)
    return cal


def linear_fisher_fixed(q,f0,shapes,state0,mass,frac,survey,cal,nuisance_level='full'):
    nd=shapes.shape[0]; nz=len(survey); nn=nz*(1 if nuisance_level=='amp' else 4)
    F=np.zeros((nd,nd)); G=np.zeros((nd,nn)); N=np.zeros((nn,nn))
    for iz,s in enumerate(survey):
        z=float(s['z']); info=state0['z'][z]
        bl,bf=float(s['b_bright']),float(s['b_faint'])
        nl,nf=float(s['n_bright_h3_Mpc3']),float(s['n_faint_h3_Mpc3'])
        V=float(s['volume_hminus3_Gpc3'])
        rows=w.mode_noise_weight(info,V,bl,bf,nl,nf); cache={}
        for ik,im,w0 in rows:
            if ik not in cache:
                cache[ik]=base.linear_basis_moments(q,f0,shapes,state0,mass,z,ik,frac)
            m00,m0,mm=cache[ik]; wt=w0*cal[z]
            ns=base.nuisance_shapes(iz,base.KOBS[ik],base.MU[im],nz,nuisance_level)
            F += wt*mm
            G += wt*np.outer(m0,ns)
            N += wt*m00*np.outer(ns,ns)
    return F-G@np.linalg.pinv(N,rcond=1e-11)@G.T


def nonlinear_fixed(q,f0,fp,fm,state0,statep,statem,mass,survey,cal,nuisance_level='full'):
    nz=len(survey); nn=nz*(1 if nuisance_level=='amp' else 4)
    Fdd=0.0; g=np.zeros(nn); N=np.zeros((nn,nn)); raw=0.0; perz=[]
    for iz,s in enumerate(survey):
        z=float(s['z']); info=state0['z'][z]
        bl,bf=float(s['b_bright']),float(s['b_faint'])
        nl,nf=float(s['n_bright_h3_Mpc3']),float(s['n_faint_h3_Mpc3'])
        V=float(s['volume_hminus3_Gpc3'])
        rows=w.mode_noise_weight(info,V,bl,bf,nl,nf); cache={}; zun=0.0
        for ik,im,w0 in rows:
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,fp,fm,state0,statep,statem,mass,z,ik)
            m00,mdd,m0d=cache[ik]; wt=w0*cal[z]
            ns=base.nuisance_shapes(iz,base.KOBS[ik],base.MU[im],nz,nuisance_level)
            Fdd += wt*mdd; raw += wt*mdd; zun += wt*mdd
            g += wt*m0d*ns
            N += wt*m00*np.outer(ns,ns)
        perz.append({'z':z,'unprojected_hidden_SN':math.sqrt(max(zun,0.0))})
    proj=max(Fdd-float(g@np.linalg.pinv(N,rcond=1e-11)@g),0.0)
    return math.sqrt(max(raw,0.0)),math.sqrt(proj),perz


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--split',required=True,choices=sorted(PROFILES))
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--z-match',type=float,default=1100.0)
    ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--samples',type=int,default=50000)
    ap.add_argument('--validate',type=int,default=5)
    ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    profile=PROFILES[args.split]
    base.ZBINS=np.array([x['z'] for x in profile],float)
    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0)
    state0=base.build_state(p0,args.mass)
    # prepare_survey computes shell volume and number density for each profile
    target_survey=list(profile); w.SURVEY=target_survey; target_survey=w.prepare_survey(state0)
    w.SURVEY=list(PROFILES['13.75']); refsurvey=w.prepare_survey(state0)
    cal=make_calibration(q,f0,state0,args.mass,refsurvey)

    P=linear_fisher_fixed(q,f0,shapes,state0,args.mass,args.frac,target_survey,cal,'full')
    selected=base.candidate_pool(P,shapes,f0,args.frac,args.samples,args.seed,args.validate)
    vals=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f'cand_{rank:02d}_plus.dat'; pm=out/f'cand_{rank:02d}_minus.dat'
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)
        sn_un,sn_full,pz=nonlinear_fixed(q,f0,fp,fm,state0,statep,statem,args.mass,target_survey,cal,'full')
        _,sn_amp,_=nonlinear_fixed(q,f0,fp,fm,state0,statep,statem,args.mass,target_survey,cal,'amp')
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={'rank':rank,'predicted_full_projected_SN':math.sqrt(max(pred2,0.0)),
             'validated_unprojected_hidden_SN':sn_un,
             'validated_after_wake_amplitude_projection_SN':sn_amp,
             'validated_after_full_odd_nuisance_projection_SN':sn_full,
             'max_relative_moment_mismatch':float(mismatch.max()),'per_redshift':pz,
             'coefficients':coeff.tolist()}
        vals.append(row); print('FIXEDCAL',args.split,json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    sn,row,fp,fm,shape=best
    np.savetxt(out/'best_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',
               header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')
    summary={'mass_eV':args.mass,'z_match':args.z_match,'halo_mass_split_log10':float(args.split),
             'calibration_reference_split_log10':13.75,
             'calibration_by_z':{str(k):v for k,v in cal.items()},
             'calibration_rule':'Absolute per-z wake normalization fixed once from 13.75 and held unchanged across all tracer partitions.',
             'survey_bins':target_survey,'best':row,'all_validations':vals,
             'interpretation':'This is the appropriate threshold-comparison control. Split-dependent tracer biases/counts enter signal and covariance, but the physical wake amplitude normalization is not recalibrated.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('BEST_FIXEDCAL',args.split,json.dumps(row))

if __name__=='__main__': main()
