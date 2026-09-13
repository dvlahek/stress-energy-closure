#!/usr/bin/env python3
"""Halo-mass proxy scatter test for the final nine-tracer wake forecast.

The validated cumulative HOD profiles provide n(>M) and b*n(>M) on the mass
threshold grid.  For each redshift bin we reconstruct smooth cumulative number
and bias-weighted number densities with monotone PCHIP interpolation, derive the
differential distributions, and convolve them with a Gaussian scatter in
log10 halo mass.  Observed tracer counts and biases are then recomputed before
the full N x N covariance Fisher calculation.

This is a scatter in the mass-observable relation, not an arbitrary relabeling
of tracers.  The physical wake normalization remains fixed to the same 13.75
two-tracer reference used in the production forecast.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import ndtr
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_multitracer_fisher as mt
import wake_desi_multitracer_fullgrid as fg
from wake_desi_split_fixedcal_fisher import make_calibration

BOUND=np.array(fg.BOUNDARIES,float)
BASELINE_SN=0.8605446904469461


def row_at(profile,z):
    return min(profile,key=lambda r:abs(float(r['z'])-float(z)))


def reconstruct_scattered_tracers(prepared, refsurvey, sigma):
    if sigma<=0:
        return mt.disjoint_tracers(tuple(BOUND.tolist()),prepared), {'mass_grid_method':'exact_disjoint_bins','sigma_dex':0.0}
    vmap={round(float(r['z']),3):float(r['volume_hminus3_Gpc3']) for r in refsurvey}
    zs=[float(r['z']) for r in prepared[float(BOUND[0])]]
    out={}; diagnostics={}
    edges=np.concatenate(([-np.inf],BOUND,[np.inf]))
    for z in zs:
        rows=[row_at(prepared[float(t)],z) for t in BOUND]
        n_hi=np.array([float(r['n_bright_h3_Mpc3']) for r in rows])
        bn_hi=np.array([float(r['n_bright_h3_Mpc3'])*float(r['b_bright']) for r in rows])
        parent_n=np.median([float(r['n_faint_h3_Mpc3'])+float(r['n_bright_h3_Mpc3']) for r in rows])
        parent_bn=np.median([float(r['n_faint_h3_Mpc3'])*float(r['b_faint'])+float(r['n_bright_h3_Mpc3'])*float(r['b_bright']) for r in rows])
        # Endpoints are deliberately far enough from the observed threshold range
        # that <=0.3 dex scatter has negligible sensitivity to their exact values.
        x=np.concatenate(([11.50],BOUND,[15.50]))
        sn=np.concatenate(([parent_n],n_hi,[0.0]))
        sb=np.concatenate(([parent_bn],bn_hi,[0.0]))
        # enforce non-increasing cumulative functions against roundoff
        sn=np.minimum.accumulate(sn)
        sb=np.minimum.accumulate(sb)
        fn=PchipInterpolator(x,sn,extrapolate=False); fb=PchipInterpolator(x,sb,extrapolate=False)
        m=np.linspace(11.50,15.50,2401)
        dn=np.clip(-fn.derivative()(m),0.0,None)
        dbn=np.clip(-fb.derivative()(m),0.0,None)
        # renormalize derivative representations exactly to the parent integrals
        intn=np.trapezoid(dn,m); intb=np.trapezoid(dbn,m)
        dn*=parent_n/max(intn,1e-300); dbn*=parent_bn/max(intb,1e-300)
        V=vmap[round(z,3)]; tracers=[]; nsum=0.0; bnsum=0.0
        for j,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
            plo=np.zeros_like(m) if np.isneginf(lo) else ndtr((lo-m)/sigma)
            phi=np.ones_like(m) if np.isposinf(hi) else ndtr((hi-m)/sigma)
            prob=np.clip(phi-plo,0.0,1.0)
            nj=float(np.trapezoid(dn*prob,m)); bnj=float(np.trapezoid(dbn*prob,m))
            if nj<=0: raise RuntimeError(f'non-positive scattered tracer at z={z}, bin={j}')
            bj=bnj/nj; nsum+=nj; bnsum+=bnj
            label=(f'<{BOUND[0]:.2f}' if j==0 else f'>{BOUND[-1]:.2f}' if j==len(BOUND) else f'{BOUND[j-1]:.2f}-{BOUND[j]:.2f}')
            tracers.append({'label':label,'N':nj*V*1e9,'n':nj,'b':bj})
        out[z]=tracers
        diagnostics[str(z)]={'parent_n':parent_n,'reconstructed_n':nsum,'parent_bias':parent_bn/parent_n,
                             'reconstructed_bias':bnsum/nsum,'number_conservation_relative':abs(nsum-parent_n)/parent_n,
                             'bias_weight_conservation_relative':abs(bnsum-parent_bn)/max(abs(parent_bn),1e-300)}
        if diagnostics[str(z)]['number_conservation_relative']>2e-4 or diagnostics[str(z)]['bias_weight_conservation_relative']>2e-4:
            raise RuntimeError('scatter reconstruction conservation failed: '+json.dumps(diagnostics[str(z)]))
    return out, {'sigma_dex':sigma,'mass_grid_method':'PCHIP cumulative HOD reconstruction + Gaussian log10M convolution','per_redshift':diagnostics}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--sigma-dex',type=float,required=True)
    ap.add_argument('--outdir',required=True); ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--z-match',type=float,default=1100.0); ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--samples',type=int,default=40000); ap.add_argument('--validate',type=int,default=3)
    ap.add_argument('--seed',type=int,default=20260913); args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    refsurvey=fg.reference_geometry(); prepared=fg.cumulative_profiles(refsurvey)
    survey,scatter_diag=reconstruct_scattered_tracers(prepared,refsurvey,args.sigma_dex)
    base.ZBINS=np.array(sorted(survey),float)
    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0); state0=base.build_state(p0,args.mass)
    cal=make_calibration(q,f0,state0,args.mass,refsurvey)
    P=mt.linear_projected(q,f0,shapes,state0,args.mass,args.frac,survey,cal,'full')
    selected=base.candidate_pool(P,shapes,f0,args.frac,args.samples,args.seed,args.validate)
    vals=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f'cand_{rank:02d}_plus.dat'; pm=out/f'cand_{rank:02d}_minus.dat'
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)
        sn_un,sn_full,pz=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,args.mass,survey,cal,'full')
        _,sn_amp,_=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,args.mass,survey,cal,'amp')
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={'rank':rank,'predicted_full_projected_SN':math.sqrt(max(pred2,0.0)),
             'validated_unprojected_hidden_SN':sn_un,'validated_after_wake_amplitude_projection_SN':sn_amp,
             'validated_after_full_odd_nuisance_projection_SN':sn_full,
             'max_relative_moment_mismatch':float(mismatch.max()),'per_redshift':pz,'coefficients':coeff.tolist()}
        vals.append(row); print('SCATTER_VALIDATE',args.sigma_dex,json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    sn,row,fp,fm,shape=best
    np.savetxt(out/'best_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')
    summary={'sigma_log10M_dex':args.sigma_dex,'mass_eV':args.mass,'pointwise_cap':args.frac,'ntracer':9,
             'baseline_zero_scatter_SN_reference':BASELINE_SN,'best':row,'all_validations':vals,
             'SN_ratio_to_zero_scatter':sn/BASELINE_SN,'information_fraction_to_zero_scatter':(sn/BASELINE_SN)**2,
             'scatter_reconstruction':scatter_diag,'tracers_by_z':{str(z):survey[z] for z in sorted(survey)},
             'calibration_rule':'Same fixed physical wake normalization as the production nine-tracer forecast; scatter changes only observed tracer assignment, counts and effective biases.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print('BEST_SCATTER',json.dumps(summary,indent=2))

if __name__=='__main__': main()
