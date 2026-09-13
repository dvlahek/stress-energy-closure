#!/usr/bin/env python3
"""Covariance-consistent multi-tracer DESI-BGS hidden-state wake Fisher.

The parent BGS sample is partitioned into 2--5 disjoint halo-mass populations
using cumulative HOD profiles at log10(M/[h^-1 Msun]) = 13.75, 14.00, 14.25,
14.50.  For each requested tracer count N we scan all boundary combinations,
construct disjoint counts and number-weighted ST99 biases, and use the full
N x N Gaussian covariance per (k,mu,z) mode.

The parity-odd response matrix is
    D_ij = i P_cb mu^2 (b_i-b_j) Y_F,
so the geometric information factor is
    Tr(C^-1 D0 C^-1 D0),
with D0 omitting the stochastic wake scalar Y_F.  For N=2 this is verified
numerically to reproduce exactly the two-tracer weight used in the existing
forecast, preventing pairwise double counting.

Absolute wake calibration is determined ONCE from the 13.75 two-tracer
reference and held fixed for every partition.  Hidden kinetic directions are
optimized after full odd-nuisance projection, then the top candidates for the
best partition are validated with full CLASS F_+ and F_- transfer functions.
"""
from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_hod_masssplit_fisher as refmod
from wake_desi_split_compare_fisher import SURVEYS
from wake_desi_split_fixedcal_fisher import make_calibration

BOUNDARIES=[13.75,14.00,14.25,14.50]
CUM={13.75:list(refmod.w.SURVEY),14.00:SURVEYS['14.00'],14.25:SURVEYS['14.25'],14.50:SURVEYS['14.50']}


def prepare_profile(profile):
    w.SURVEY=list(profile)
    return w.prepare_survey(None)


def row_at(profile,z):
    return min(profile,key=lambda r:abs(float(r['z'])-float(z)))


def disjoint_tracers(boundaries, prepared):
    """Return {z:[{N,n,b,label},...]} from cumulative low/high profiles."""
    zs=[float(x['z']) for x in prepared[boundaries[0]]]
    out={}
    for z in zs:
        rows=[row_at(prepared[t],z) for t in boundaries]
        tr=[]
        # lowest cumulative population
        r=rows[0]
        tr.append({'label':f'<{boundaries[0]:.2f}','N':float(r['N_faint']),
                   'n':float(r['n_faint_h3_Mpc3']),'b':float(r['b_faint'])})
        # differences of cumulative low populations
        for i in range(len(boundaries)-1):
            ra,rb=rows[i],rows[i+1]
            Na,Nb=float(ra['N_faint']),float(rb['N_faint'])
            nn=Nb-Na
            if nn<=0: raise RuntimeError('non-positive disjoint tracer count')
            bb=(Nb*float(rb['b_faint'])-Na*float(ra['b_faint']))/nn
            V=float(rb['volume_hminus3_Gpc3'])
            tr.append({'label':f'{boundaries[i]:.2f}-{boundaries[i+1]:.2f}',
                       'N':nn,'n':nn/(V*1e9),'b':bb})
        # highest tail
        r=rows[-1]
        tr.append({'label':f'>{boundaries[-1]:.2f}','N':float(r['N_bright']),
                   'n':float(r['n_bright_h3_Mpc3']),'b':float(r['b_bright'])})
        out[z]=tr
    return out


def geom_factor(info,ik,mu,tracers):
    P=float(info['pk'][ik]); fg=float(info['f'][ik])
    b=np.array([t['b'] for t in tracers],float)
    n=np.array([t['n'] for t in tracers],float)
    a=b+fg*mu**2
    C=P*np.outer(a,a)+np.diag(1.0/n)
    Cinv=np.linalg.inv(C)
    Q=b[:,None]-b[None,:]
    D=1j*P*mu**2*Q
    val=float(np.trace(Cinv@D@Cinv@D).real)
    return max(val,0.0)


def two_tracer_equivalence(info,ik,mu,tracers):
    if len(tracers)!=2: return 0.0
    P=float(info['pk'][ik]); fg=float(info['f'][ik])
    b0,b1=tracers[0]['b'],tracers[1]['b']; n0,n1=tracers[0]['n'],tracers[1]['n']
    a0=b0+fg*mu**2; a1=b1+fg*mu**2
    C00=a0*a0*P+1/n0; C11=a1*a1*P+1/n1; C01=a0*a1*P
    det=C00*C11-C01*C01
    old=P*P*(b1-b0)**2*mu**4/(0.5*det)
    new=geom_factor(info,ik,mu,tracers)
    return abs(new-old)/max(abs(old),1e-300)


def mode_rows(info,V,tracers):
    dk=np.gradient(base.KOBS); rows=[]
    for ik,k in enumerate(base.KOBS):
        for im,mu in enumerate(base.MU):
            geom=geom_factor(info,ik,float(mu),tracers)
            mode=V*1e9/(4*np.pi**2)*k*k*dk[ik]*base.MUW[im]
            rows.append((ik,im,mode*geom))
    return rows


def linear_projected(q,f0,shapes,state0,mass,frac,survey_by_z,cal,nuisance='full'):
    nd=shapes.shape[0]; zs=sorted(survey_by_z); nz=len(zs)
    nn=nz*(1 if nuisance=='amp' else 4)
    F=np.zeros((nd,nd)); G=np.zeros((nd,nn)); N=np.zeros((nn,nn))
    for iz,z in enumerate(zs):
        info=state0['z'][z]; trs=survey_by_z[z]
        # shell volume is the same for every cumulative profile; recover from total N/n
        V=sum(t['N'] for t in trs)/sum(t['n'] for t in trs)/1e9
        rows=mode_rows(info,V,trs); cache={}
        for ik,im,w0 in rows:
            if ik not in cache:
                cache[ik]=base.linear_basis_moments(q,f0,shapes,state0,mass,z,ik,frac)
            m00,m0,mm=cache[ik]; wt=w0*cal[z]
            ns=base.nuisance_shapes(iz,base.KOBS[ik],base.MU[im],nz,nuisance)
            F+=wt*mm; G+=wt*np.outer(m0,ns); N+=wt*m00*np.outer(ns,ns)
    return F-G@np.linalg.pinv(N,rcond=1e-11)@G.T


def nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,survey_by_z,cal,nuisance='full'):
    zs=sorted(survey_by_z); nz=len(zs); nn=nz*(1 if nuisance=='amp' else 4)
    Fdd=0.; raw=0.; g=np.zeros(nn); N=np.zeros((nn,nn)); perz=[]
    for iz,z in enumerate(zs):
        info=state0['z'][z]; trs=survey_by_z[z]
        V=sum(t['N'] for t in trs)/sum(t['n'] for t in trs)/1e9
        rows=mode_rows(info,V,trs); cache={}; zraw=0.
        for ik,im,w0 in rows:
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,fp,fm,state0,statep,statem,mass,z,ik)
            m00,mdd,m0d=cache[ik]; wt=w0*cal[z]
            ns=base.nuisance_shapes(iz,base.KOBS[ik],base.MU[im],nz,nuisance)
            Fdd+=wt*mdd; raw+=wt*mdd; zraw+=wt*mdd
            g+=wt*m0d*ns; N+=wt*m00*np.outer(ns,ns)
        perz.append({'z':z,'unprojected_hidden_SN':math.sqrt(max(zraw,0.0))})
    proj=max(Fdd-float(g@np.linalg.pinv(N,rcond=1e-11)@g),0.0)
    return math.sqrt(max(raw,0.0)),math.sqrt(proj),perz


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ntracer',type=int,required=True,choices=[2,3,4,5])
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--z-match',type=float,default=1100.)
    ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--samples',type=int,default=50000)
    ap.add_argument('--validate',type=int,default=5)
    ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    prepared={t:prepare_profile(CUM[t]) for t in BOUNDARIES}
    zs=np.array([float(r['z']) for r in prepared[13.75]],float); base.ZBINS=zs
    q=np.linspace(0.,20.,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0); state0=base.build_state(p0,args.mass)
    # physical normalization is fixed once from the 13.75 two-tracer reference
    cal=make_calibration(q,f0,state0,args.mass,prepared[13.75])

    # exact N=2 matrix reduction check at representative and extreme modes
    ref2=disjoint_tracers((13.75,),prepared)
    checks=[]
    for z in [float(zs[0]),float(zs[len(zs)//2]),float(zs[-1])]:
        for ik in [0,len(base.KOBS)//2,len(base.KOBS)-1]:
            for im in [0,len(base.MU)//2,len(base.MU)-1]:
                checks.append(two_tracer_equivalence(state0['z'][z],ik,float(base.MU[im]),ref2[z]))
    maxeq=max(checks)
    if maxeq>1e-10: raise RuntimeError(f'N=2 covariance reduction failed: {maxeq}')

    configs=list(itertools.combinations(BOUNDARIES,args.ntracer-1))
    scan=[]; config_cache={}
    for boundaries in configs:
        survey=disjoint_tracers(boundaries,prepared); config_cache[boundaries]=survey
        P=linear_projected(q,f0,shapes,state0,args.mass,args.frac,survey,cal,'full')
        top=base.candidate_pool(P,shapes,f0,args.frac,args.samples,args.seed,1)[0]
        pred=math.sqrt(max(top[0],0.0))
        row={'boundaries':list(boundaries),'predicted_best_projected_SN':pred}
        scan.append(row); print('CONFIG_SCAN',json.dumps(row))
    scan.sort(key=lambda r:r['predicted_best_projected_SN'],reverse=True)
    bestb=tuple(scan[0]['boundaries']); survey=config_cache[bestb]
    P=linear_projected(q,f0,shapes,state0,args.mass,args.frac,survey,cal,'full')
    selected=base.candidate_pool(P,shapes,f0,args.frac,args.samples,args.seed,args.validate)
    vals=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f'cand_{rank:02d}_plus.dat'; pm=out/f'cand_{rank:02d}_minus.dat'
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)
        sn_un,sn_full,pz=nonlinear_projected(q,f0,fp,fm,state0,statep,statem,args.mass,survey,cal,'full')
        _,sn_amp,_=nonlinear_projected(q,f0,fp,fm,state0,statep,statem,args.mass,survey,cal,'amp')
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={'rank':rank,'predicted_full_projected_SN':math.sqrt(max(pred2,0.0)),
             'validated_unprojected_hidden_SN':sn_un,
             'validated_after_wake_amplitude_projection_SN':sn_amp,
             'validated_after_full_odd_nuisance_projection_SN':sn_full,
             'max_relative_moment_mismatch':float(mismatch.max()),
             'per_redshift':pz,'coefficients':coeff.tolist()}
        vals.append(row); print('MULTITRACER',args.ntracer,json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    sn,row,fp,fm,shape=best
    np.savetxt(out/'best_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',
               header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')
    summary={'mass_eV':args.mass,'ntracer':args.ntracer,'boundary_grid':BOUNDARIES,
             'best_boundaries':list(bestb),'configuration_scan':scan,
             'N2_reduction_max_relative_error':maxeq,
             'calibration_reference_split_log10':13.75,
             'calibration_rule':'One fixed physical per-z normalization from the 13.75 two-tracer reference, held fixed for all N and partitions.',
             'tracers_by_z':{str(z):survey[z] for z in sorted(survey)},
             'best':row,'all_validations':vals,
             'interpretation':'Full N x N covariance Fisher. Cross-spectra are not summed independently, so shared cosmic variance is retained and pairwise double counting is avoided.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('BEST_MULTITRACER',args.ntracer,json.dumps(row))

if __name__=='__main__': main()
