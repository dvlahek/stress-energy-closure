#!/usr/bin/env python3
"""Shape-only wake/Doppler templates for genuine DESI DR1 cross-population tests.

The hidden-state direction is frozen to the same seven-coefficient direction used
in the manuscript.  Unlike the BGS forecast calibration, this script does not
import a survey-specific absolute normalization.  The output is therefore a
shape template for a free matched-filter amplitude.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.special import spherical_jn

import class_response_optimize as cro
import wake_two_tracer_fisher as base

COEFF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def signed_hidden_response(q, f0, fp, fm, state0, statep, statem, mass, z):
    vals=[]
    for ik in range(len(base.KOBS)):
        m00,mdd,m0d = base.stochastic_moments(
            q,f0,fp,fm,state0,statep,statem,mass,float(z),ik
        )
        sign = 1.0 if m0d >= 0.0 else -1.0
        vals.append(sign*math.sqrt(max(mdd,0.0)))
    return np.asarray(vals,float)


def hankel_dipole(k, pk, amp, s):
    kk=k[:,None]; ss=s[None,:]
    integrand=(k*k*pk*amp)[:,None]*spherical_jn(1,kk*ss)
    return np.trapezoid(integrand,k,axis=0)/(2.0*np.pi**2)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--mass",type=float,default=0.06)
    ap.add_argument("--z-match",type=float,default=1100.0)
    ap.add_argument("--frac",type=float,default=0.30)
    ap.add_argument("--z-grid",default="0.375,0.425,0.825,0.875,0.925,0.975,1.025,1.075")
    ap.add_argument("--s-min",type=float,default=20.0)
    ap.add_argument("--s-max",type=float,default=140.0)
    ap.add_argument("--s-step",type=float,default=2.0)
    args=ap.parse_args()

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    zgrid=np.asarray([float(x) for x in args.z_grid.split(",")],float)
    base.ZBINS=zgrid.copy()

    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    if shapes.shape[0] != len(COEFF):
        raise RuntimeError(f"production basis changed: {shapes.shape[0]} != {len(COEFF)}")

    raw=COEFF@shapes
    maxrel=float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    shape=raw/maxrel
    fp=f0+args.frac*shape
    fm=f0-args.frac*shape
    mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
    mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)

    p0=out/"fd_reference.dat"; pp=out/"hidden_plus.dat"; pm=out/"hidden_minus.dat"
    base.write_psd(p0,q,f0); base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
    state0=base.build_state(p0,args.mass)
    statep=base.build_state(pp,args.mass)
    statem=base.build_state(pm,args.mass)

    s=np.arange(args.s_min,args.s_max+0.5*args.s_step,args.s_step)
    curves=[]
    for z in zgrid:
        info=state0["z"][float(z)]
        k=np.asarray(base.KOBS,float)
        pk=np.asarray(info["pk"],float)
        hid=signed_hidden_response(q,f0,fp,fm,state0,statep,statem,args.mass,float(z))
        wake=hankel_dipole(k,pk,hid,s)
        dop=hankel_dipole(k,pk,float(info["H"])/np.maximum(k,1e-8),s)
        curves.append((float(z),wake,dop))

    # One global normalization per template preserves relative z evolution.
    ws=max(float(np.max(np.abs(w))) for _,w,_ in curves)
    ds=max(float(np.max(np.abs(d))) for _,_,d in curves)
    if not np.isfinite(ws+ds) or ws<=0 or ds<=0:
        raise RuntimeError("degenerate template normalization")

    rows=[]
    for z,w,d in curves:
        for si,wv,dv in zip(s,w/ws,d/ds):
            rows.append((z,si,wv,dv))
    arr=np.asarray(rows,float)
    np.savetxt(
        out/"genuine_multitracer_templates.csv",arr,delimiter=",",
        header="z,s_Mpc_over_h,wake_shape,doppler_shape",comments=""
    )
    summary={
        "scope":"Shape-only templates for DESI DR1 genuine cross-population odd-dipole test.",
        "absolute_likelihood_claim":False,
        "mass_eV":args.mass,
        "z_match":args.z_match,
        "pointwise_cap":args.frac,
        "coefficients":COEFF.tolist(),
        "max_relative_moment_mismatch":float(mismatch.max()),
        "z_grid":zgrid.tolist(),
        "s_grid_Mpc_over_h":s.tolist(),
        "wake_definition":"Hankel j1 transform of P_cb times signed frozen hidden-state response; free global amplitude in the data fit.",
        "doppler_definition":"Hankel j1 transform of P_cb H(z)/k; nuisance amplitude in each analysis redshift block.",
        "normalization":"One global max-absolute normalization per template over the full z-s grid; no BGS-specific absolute calibration.",
        "global_wake_scale":ws,
        "global_doppler_scale":ds,
    }
    (out/"template_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("GENUINE_MULTITRACER_TEMPLATE",json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
