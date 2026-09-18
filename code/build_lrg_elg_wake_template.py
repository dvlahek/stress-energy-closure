#!/usr/bin/env python3
"""Build a shape-only hidden-state wake template at the measured LRGxELG effective redshift.

Reuses the frozen publication hidden-state direction from code/wake_phase7_template.py
but does NOT import the low-z BGS HOD amplitude calibration. For one broad z-bin,
the fitted wake amplitude absorbs the overall LRG-ELG normalization.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_phase7_template as wpt

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',required=True); ap.add_argument('--z',type=float,required=True)
    ap.add_argument('--mass',type=float,default=0.06); ap.add_argument('--z-match',type=float,default=1100.)
    ap.add_argument('--frac',type=float,default=.30); ap.add_argument('--s-min',type=float,default=20.)
    ap.add_argument('--s-max',type=float,default=140.); ap.add_argument('--s-step',type=float,default=20.)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    base.ZBINS=np.array([float(a.z)])
    q=np.linspace(0,20,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,a.mass,a.z_match)
    if shapes.shape[0]!=len(wpt.COEFF): raise RuntimeError('production kinetic basis changed')
    raw=wpt.COEFF@shapes
    shape=raw/float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    fp=f0+a.frac*shape; fm=f0-a.frac*shape

    p0=out/'fd_reference.dat'; pp=out/'hidden_plus.dat'; pm=out/'hidden_minus.dat'
    base.write_psd(p0,q,f0); base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
    st0=base.build_state(p0,a.mass); stp=base.build_state(pp,a.mass); stm=base.build_state(pm,a.mass)

    z=float(a.z); info=st0['z'][z]; k=np.asarray(base.KOBS,float); pk=np.asarray(info['pk'],float)
    hid=wpt.signed_hidden_response(q,f0,fp,fm,st0,stp,stm,a.mass,z)
    H=float(info['H']); s=np.arange(a.s_min,a.s_max+.5*a.s_step,a.s_step)
    wake=wpt.hankel_dipole(k,pk,hid,s); dop=wpt.hankel_dipole(k,pk,H/np.maximum(k,1e-8),s)
    wake/=max(float(np.max(np.abs(wake))),1e-300); dop/=max(float(np.max(np.abs(dop))),1e-300)
    np.savetxt(out/'lrg_elg_shape_template.csv',np.c_[s,wake,dop],delimiter=',',
               header='s_Mpc_over_h,wake_shape,doppler_shape',comments='')
    mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
    mis=np.abs(mp-mm)/np.maximum(.5*(np.abs(mp)+np.abs(mm)),1e-300)
    S=dict(scope='shape-only hidden-state wake matched-filter template',absolute_wake_prediction=False,
           z_effective=z,mass_eV=a.mass,z_match=a.z_match,pointwise_cap=a.frac,
           max_relative_moment_mismatch=float(mis.max()),
           caveat='No low-z BGS HOD amplitude calibration is used; fit only a free template coefficient.')
    (out/'template_summary.json').write_text(json.dumps(S,indent=2)+'\n')
    print(json.dumps(S,indent=2))
if __name__=='__main__': main()
