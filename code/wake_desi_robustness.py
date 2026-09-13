#!/usr/bin/env python3
"""Final robustness sweeps for the DESI-BGS hidden-state wake forecast.

Uses the saturated nine-tracer halo-mass partition, the full N x N covariance,
and the fixed physical wake-normalization convention validated in the production
forecast. Three one-factor sweeps are supported:

  mass : m_nu = 0.05--0.10 eV, with the standard-FD physical calibration
         recomputed for each mass at the baseline kmax=0.10 h/Mpc.
  cap  : pointwise F_+-F0 and F_--F0 cap = 10, 20, 30 percent at m=0.06 eV.
  kmax : observed kmax varied while the absolute wake calibration is first fixed
         at kmax=0.10 h/Mpc and then held unchanged as high-k modes are removed.

Every point re-optimizes the matched-moment hidden direction after the full
odd-nuisance projection and validates the best candidates with nonlinear CLASS
F_+ and F_- transfer functions.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_multitracer_fisher as mt
import wake_desi_multitracer_fullgrid as fg
from wake_desi_split_fixedcal_fisher import make_calibration

BASE_MASS=0.06
BASE_FRAC=0.30
BASE_KMIN=0.003
BASE_KMAX=0.10
BASE_NK=48


def set_kgrid(kmax: float):
    # Preserve roughly the baseline logarithmic mode density as kmax changes.
    n=max(16,int(round(BASE_NK*np.log(kmax/BASE_KMIN)/np.log(BASE_KMAX/BASE_KMIN))))
    base.KOBS=np.geomspace(BASE_KMIN,float(kmax),n)
    return n


def survey_geometry():
    refsurvey=fg.reference_geometry()
    prepared=fg.cumulative_profiles(refsurvey)
    survey=mt.disjoint_tracers(tuple(fg.BOUNDARIES),prepared)
    return refsurvey,survey


def physical_calibration(q,f0,mass,refsurvey,out:Path):
    # Always determine the physical normalization on the baseline k-range.
    set_kgrid(BASE_KMAX)
    p0=out/'calibration_fd_reference.dat'
    base.write_psd(p0,q,f0)
    state0=base.build_state(p0,mass)
    w.SURVEY=list(refsurvey)
    ref=w.prepare_survey(state0)
    return make_calibration(q,f0,state0,mass,ref)


def evaluate(mode:str,value:float,out:Path,samples:int,validate:int,seed:int):
    mass=BASE_MASS; frac=BASE_FRAC; kmax=BASE_KMAX
    if mode=='mass': mass=float(value)
    elif mode=='cap': frac=float(value)
    elif mode=='kmax': kmax=float(value)
    else: raise ValueError(mode)

    refsurvey,survey=survey_geometry()
    base.ZBINS=np.array(sorted(survey),float)
    q=np.linspace(0.,20.,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,mass,1100.)

    # Physical normalization is defined before any k truncation.
    cal=physical_calibration(q,f0,mass,refsurvey,out)

    nk=set_kgrid(kmax)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0)
    state0=base.build_state(p0,mass)

    P=mt.linear_projected(q,f0,shapes,state0,mass,frac,survey,cal,'full')
    selected=base.candidate_pool(P,shapes,f0,frac,samples,seed,validate)
    vals=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+frac*shape; fm=f0-frac*shape
        pp=out/f'cand_{rank:02d}_plus.dat'; pm=out/f'cand_{rank:02d}_minus.dat'
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,mass); statem=base.build_state(pm,mass)
        sn_un,sn_full,pz=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,survey,cal,'full')
        _,sn_amp,_=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,survey,cal,'amp')
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={'rank':rank,'predicted_full_projected_SN':math.sqrt(max(pred2,0.0)),
             'validated_unprojected_hidden_SN':sn_un,
             'validated_after_wake_amplitude_projection_SN':sn_amp,
             'validated_after_full_odd_nuisance_projection_SN':sn_full,
             'max_relative_moment_mismatch':float(mismatch.max()),
             'per_redshift':pz,'coefficients':coeff.tolist()}
        vals.append(row); print('ROBUST',mode,value,json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    sn,row,fp,fm,shape=best
    np.savetxt(out/'best_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',
               header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')
    summary={'sweep':mode,'value':value,'mass_eV':mass,'pointwise_cap':frac,
             'kmin_h_Mpc':BASE_KMIN,'kmax_h_Mpc':kmax,'nk':nk,'ntracer':9,
             'partition_boundaries':fg.BOUNDARIES,
             'calibration_kmax_h_Mpc':BASE_KMAX,
             'calibration_rule':('Physical wake normalization is always fixed on the baseline '
                                 '0.003<=k<=0.10 h/Mpc range. For kmax sweeps the same calibration '
                                 'is held fixed after truncating observed modes.'),
             'best':row,'all_validations':vals}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('BEST_ROBUST',mode,value,json.dumps(row))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=['mass','cap','kmax'],required=True)
    ap.add_argument('--value',type=float,required=True)
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--samples',type=int,default=30000)
    ap.add_argument('--validate',type=int,default=3)
    ap.add_argument('--seed',type=int,default=20260913)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    evaluate(a.mode,a.value,out,a.samples,a.validate,a.seed)

if __name__=='__main__': main()
