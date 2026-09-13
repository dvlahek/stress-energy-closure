#!/usr/bin/env python3
"""Uncalibrated relative ranking of DESI-BGS halo-mass splits.

This diagnostic removes the per-redshift Okoli Eq.35 calibration entirely.
Every split is evaluated with the same FD neutrino wake kernel, the same CLASS
transfer functions, and its own tracer covariance from the reconstructed HOD
counts/biases.  Absolute units are intentionally arbitrary.  Only ratios among
splits are interpreted.

Two rankings are reported:
  1) covariance_only: integral of the common signal prefactor before the wake
     velocity/occupation kernel, and
  2) raw_FD_wake: the same covariance multiplied by the identical FD wake
     stochastic kernel m00(k,z).
If both prefer thresholds above 10^13.75, the preference is not introduced by
our literature normalization step.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from classy import Class
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
from wake_desi_split_compare_fisher import SURVEYS

MASS=0.06
Z_MATCH=1100.0
AREA_DEG2=14000.0
OUT=Path('wake_split_raw_rank_output')


def prepare(profile):
    p=Class()
    p.set({"H0":cro.H0,"omega_b":cro.OMEGA_B,"omega_cdm":cro.OMEGA_CDM,
           "A_s":cro.A_S,"n_s":base.NS,"tau_reio":cro.TAU_REIO,
           "N_ur":cro.N_UR,"N_ncdm":0})
    p.compute()
    try:
        out=[]
        for s in profile:
            r=dict(s)
            V=w.shell_volume_from_class(p,r['zlo'],r['zhi'])
            r['volume_hminus3_Gpc3']=float(V)
            r['n_faint_h3_Mpc3']=float(r['N_faint']/(V*1e9))
            r['n_bright_h3_Mpc3']=float(r['N_bright']/(V*1e9))
            r['delta_b']=float(r['b_bright']-r['b_faint'])
            out.append(r)
        return out
    finally:
        p.struct_cleanup(); p.empty()


def rank_one(q,f0,state0,profile):
    survey=prepare(profile)
    cov_only=0.0
    raw_fd=0.0
    perz=[]
    for s in survey:
        z=float(s['z']); info=state0['z'][z]
        rows=w.mode_noise_weight(info,float(s['volume_hminus3_Gpc3']),
                                 float(s['b_bright']),float(s['b_faint']),
                                 float(s['n_bright_h3_Mpc3']),float(s['n_faint_h3_Mpc3']))
        zcov=0.0; zraw=0.0; cache={}
        for ik,im,wt in rows:
            zcov += wt
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,f0,f0,state0,state0,state0,MASS,z,ik)[0]
            zraw += wt*cache[ik]
        cov_only += zcov; raw_fd += zraw
        perz.append({'z':z,'covariance_only':zcov,'raw_FD_wake':zraw})
    return {'covariance_only':cov_only,'raw_FD_wake':raw_fd,'per_redshift':perz,
            'survey':survey}


def main():
    OUT.mkdir(exist_ok=True)
    # Add the already validated 13.75 profile explicitly for a common control.
    from wake_desi_hod_masssplit_fisher import w as wh
    profile_1375=list(wh.SURVEY)
    profiles={'13.75':profile_1375,**SURVEYS}
    base.ZBINS=np.array([0.075,0.125,0.175,0.225,0.275,0.325,0.375],float)
    q=np.linspace(0.0,20.0,4000)
    f0,*_=cro.kinetic_objects(q,MASS,Z_MATCH)
    fd=OUT/'fd_reference.dat'; base.write_psd(fd,q,f0)
    state0=base.build_state(fd,MASS)
    results={k:rank_one(q,f0,state0,v) for k,v in profiles.items()}
    ref=results['13.75']
    for k,v in results.items():
        v['covariance_ratio_to_13p75']=v['covariance_only']/ref['covariance_only']
        v['raw_FD_ratio_to_13p75']=v['raw_FD_wake']/ref['raw_FD_wake']
        print('RAW_SPLIT',k,json.dumps({x:v[x] for x in ['covariance_only','raw_FD_wake','covariance_ratio_to_13p75','raw_FD_ratio_to_13p75']}))
    best_cov=max(results,key=lambda k:results[k]['covariance_only'])
    best_raw=max(results,key=lambda k:results[k]['raw_FD_wake'])
    summary={'mass_eV':MASS,'z_match':Z_MATCH,'normalization':'none',
             'profiles_source':'validated BGS HOD threshold profiles, run 34754436382',
             'best_covariance_only':best_cov,'best_raw_FD_wake':best_raw,
             'results':results,
             'interpretation':'Only relative split ranking is meaningful. No Okoli Eq.35 or Ge et al. absolute S/N normalization is applied.'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__': main()
