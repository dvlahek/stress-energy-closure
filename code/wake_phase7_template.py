#!/usr/bin/env python3
"""Build a forecast-derived configuration-space template for the Phase-7 DESI test.

The hidden-state direction is the exact best zero-scatter nine-tracer direction
from the final publication validation run (2026-09-13, seed 20260913).  The
state is reconstructed from the seven coefficients and the production
kinetic basis, so no exploratory optimization is repeated here.

The output is intended as a *matched-filter shape* for an end-to-end survey
control.  The standard odd/Doppler template is included as an independent
free nuisance.  Absolute survey-likelihood normalization and a full DESI
window convolution are deliberately kept separate from this shape test.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as bonvin
import wake_desi_hod_masssplit_fisher as refmod
from wake_desi_split_fixedcal_fisher import make_calibration

# Exact coefficient vector recorded in the final sigma_log10M=0 publication run.
COEFF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)
REFERENCE_SN = 0.8605446904469461
REFERENCE_SEED = 20260913
REFERENCE_CLASS = cro.CLASS_COMMIT


def signed_hidden_response(q, f0, fp, fm, state0, statep, statem, mass, z):
    vals=[]
    for ik in range(len(base.KOBS)):
        m00,mdd,m0d = base.stochastic_moments(q,f0,fp,fm,state0,statep,statem,mass,z,ik)
        sign = 1.0 if m0d >= 0.0 else -1.0
        vals.append(sign*math.sqrt(max(mdd,0.0)))
    return np.asarray(vals,float)


def hankel_dipole(k, pk, amp, s):
    # Overall angular and Fourier-convention constants are irrelevant for the
    # matched-filter shape because amplitudes are fit freely in the data step.
    kk=k[:,None]; ss=s[None,:]
    integrand=(k*k*pk*amp)[:,None]*spherical_jn(1,kk*ss)
    return np.trapezoid(integrand,k,axis=0)/(2.0*np.pi**2)


def norm_shape(x):
    x=np.asarray(x,float)
    scale=float(np.max(np.abs(x)))
    if not np.isfinite(scale) or scale<=0: raise RuntimeError('degenerate template')
    return x/scale, scale


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--z-match',type=float,default=1100.0)
    ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--s-min',type=float,default=20.0)
    ap.add_argument('--s-max',type=float,default=140.0)
    ap.add_argument('--s-step',type=float,default=2.0)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    zgrid=np.array([0.125,0.175,0.225,0.275,0.325,0.375],float)
    base.ZBINS=zgrid
    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    if shapes.shape[0] != len(COEFF):
        raise RuntimeError(f'production basis changed: {shapes.shape[0]} != {len(COEFF)}')
    raw=COEFF@shapes
    maxrel=float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    shape=raw/maxrel
    fp=f0+args.frac*shape; fm=f0-args.frac*shape
    mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
    mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)

    p0=out/'fd_reference.dat'; pp=out/'hidden_plus.dat'; pm=out/'hidden_minus.dat'
    base.write_psd(p0,q,f0); base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
    state0=base.build_state(p0,args.mass); statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)

    # The publication calibration is fixed by the 13.75 halo-mass reference.
    bonvin.SURVEY=list(refmod.w.SURVEY)
    refsurvey=bonvin.prepare_survey(state0)
    cal=make_calibration(q,f0,state0,args.mass,refsurvey)

    s=np.arange(args.s_min,args.s_max+0.5*args.s_step,args.s_step)
    rows=[]
    for z in zgrid:
        info=state0['z'][float(z)]
        k=np.asarray(base.KOBS,float); pk=np.asarray(info['pk'],float)
        hid=signed_hidden_response(q,f0,fp,fm,state0,statep,statem,args.mass,float(z))
        # sqrt(cal) is the production signal-amplitude calibration. We retain it
        # for relative z/k weighting but still fit a free data amplitude.
        cz=float(cal[min(cal,key=lambda zz:abs(float(zz)-float(z)))])
        wake_amp=np.sqrt(cz)*hid
        # Standard relativistic odd sector has the characteristic H/k scaling;
        # only the shape is used here because its amplitude is a nuisance.
        H=float(info['H'])
        dop_amp=(H/np.maximum(k,1e-8))
        wake=hankel_dipole(k,pk,wake_amp,s)
        dop=hankel_dipole(k,pk,dop_amp,s)
        wake_n,wscale=norm_shape(wake); dop_n,dscale=norm_shape(dop)
        for si,wn,dn in zip(s,wake_n,dop_n):
            rows.append((z,si,wn,dn))

    arr=np.array(rows,float)
    np.savetxt(out/'phase7_templates.csv',arr,delimiter=',',
               header='z,s_Mpc_over_h,wake_shape,doppler_shape',comments='')
    summary={
        'scope':'forecast-derived matched-filter shapes for the DESI DR1 Phase-7 real-survey control',
        'absolute_likelihood_claim':False,
        'mass_eV':args.mass,'z_match':args.z_match,'pointwise_cap':args.frac,
        'reference_projected_SN':REFERENCE_SN,'reference_seed':REFERENCE_SEED,
        'class_commit':REFERENCE_CLASS,'coefficients':COEFF.tolist(),
        'max_relative_moment_mismatch':float(mismatch.max()),
        'z_grid':zgrid.tolist(),'s_grid_Mpc_over_h':s.tolist(),
        'wake_shape_definition':'Hankel j1 transform of P_cb(k,z) times signed hidden-state stochastic response, retaining sqrt(publication calibration) relative weighting; free amplitude in data fit.',
        'doppler_shape_definition':'Hankel j1 transform of P_cb(k,z) H(z)/k; free amplitude nuisance.',
        'normalization':'Each template is independently normalized to max(abs(template))=1 at each z.',
    }
    (out/'template_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PHASE7_TEMPLATE',json.dumps(summary,indent=2))

if __name__=='__main__': main()
