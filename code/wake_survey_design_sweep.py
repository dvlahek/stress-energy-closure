#!/usr/bin/env python3
"""Multi-tracer design sweep for the 0.06 eV hidden-state wake signal.

For each (b_l,b_f,n_l,n_f) design, re-optimize the exact matched-moment null
space using the conservative per-redshift wake-amplitude plus odd-broadband
projection. Rank designs at V=1 (Gpc/h)^3, then validate the most promising
ones with full CLASS F+/F- transfer functions. Report required survey volume
and galaxy count for 1,2,3 sigma hidden-state discrimination.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_calibration_hierarchy as hier

BIAS_PAIRS=[(2.0,1.0),(2.3,1.0),(2.5,1.0),(2.8,1.0),(3.0,1.0),(2.5,1.2),(3.0,1.2)]
DENSITY_PAIRS=[(0.003,0.010),(0.005,0.010),(0.005,0.020),(0.010,0.020),(0.020,0.020)]


def set_design(bl,bf,nl,nf):
    base.BL=float(bl); base.BF=float(bf); base.NL=float(nl); base.NF=float(nf)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--outdir',default='wake_survey_design_output')
    ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--volume',type=float,default=1.0)
    ap.add_argument('--samples',type=int,default=25000)
    ap.add_argument('--validate',type=int,default=6)
    ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,1100.0)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0)
    state0=base.build_state(p0,args.mass)

    rows=[]
    idx=0
    for bl,bf in BIAS_PAIRS:
        for nl,nf in DENSITY_PAIRS:
            set_design(bl,bf,nl,nf)
            Bp,Br,perz=hier.calibrated_linear_fisher(q,f0,shapes,state0,args.mass,args.frac,args.volume,'perbin')
            sel=hier.candidate_pool(Bp,shapes,f0,args.samples,args.seed+idx,1)[0]
            pred=math.sqrt(max(sel[0],0.0))
            rows.append({
                'index':idx,'b_l':bl,'b_f':bf,'delta_b':bl-bf,'n_l':nl,'n_f':nf,
                'predicted_projected_SN_V1':pred,
                'predicted_V_1sigma':1.0/max(pred*pred,1e-30),
                'predicted_V_3sigma':9.0/max(pred*pred,1e-30),
                '_candidate':sel,
            })
            idx+=1
    rows.sort(key=lambda r:r['predicted_projected_SN_V1'],reverse=True)

    validations=[]
    for rank,r in enumerate(rows[:args.validate]):
        set_design(r['b_l'],r['b_f'],r['n_l'],r['n_f'])
        pred2,coeff,norm,shape=r['_candidate']
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f'val_{rank:02d}_plus.dat'; pm=out/f'val_{rank:02d}_minus.dat'
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)
        un,proj,perz=hier.nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,args.mass,args.volume,'perbin')
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        v1=args.volume/max(proj*proj,1e-30); v2=4*v1; v3=9*v1
        nsum=r['n_l']+r['n_f']
        vr={
            'rank':rank,'index':r['index'],'b_l':r['b_l'],'b_f':r['b_f'],'delta_b':r['delta_b'],
            'n_l':r['n_l'],'n_f':r['n_f'],
            'predicted_projected_SN_V1':r['predicted_projected_SN_V1'],
            'validated_unprojected_SN_V1':un,
            'validated_projected_SN_V1':proj,
            'volume_1sigma_hminus3_Gpc3':v1,
            'volume_2sigma_hminus3_Gpc3':v2,
            'volume_3sigma_hminus3_Gpc3':v3,
            'galaxies_1sigma':nsum*v1*1e9,
            'galaxies_2sigma':nsum*v2*1e9,
            'galaxies_3sigma':nsum*v3*1e9,
            'max_relative_moment_mismatch':float(mismatch.max()),
            'per_redshift':perz,
        }
        validations.append(vr); print(json.dumps(vr))
        np.savetxt(out/f'best_design_rank{rank:02d}_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',
                   header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')

    public_rows=[]
    for r in rows:
        rr={k:v for k,v in r.items() if k!='_candidate'}
        public_rows.append(rr)
    summary={
        'class_commit':cro.CLASS_COMMIT,'mass_eV':args.mass,'fractional_distortion_cap':args.frac,
        'baseline_volume_hminus3_Gpc3':args.volume,
        'projection':'independent wake amplitude plus k^-1, k, k^2 odd broadband in each redshift bin',
        'linear_design_grid':public_rows,
        'nonlinear_validations':validations,
        'interpretation':'Literature-calibrated Okoli-type two-tracer survey-design sweep. Absolute FD normalization follows Eq.35 scaling; hidden-state shape and F+/F- transfer response are recomputed. Top designs are full-CLASS validated.'
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')

if __name__=='__main__': main()
