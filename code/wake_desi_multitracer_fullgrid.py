#!/usr/bin/env python3
"""Full-grid 5--9 tracer DESI-BGS hidden-state wake forecast.

Extends the covariance-consistent multi-tracer calculation to all validated HOD
mass thresholds 12.75, 13.00, ..., 14.50.  For each requested N, every choice
of N-1 boundaries is ranked using the already stable leading hidden-state
direction.  The best partition is then fully re-optimized in the matched-moment
null space and validated with nonlinear CLASS F_+ and F_- transfers.

The absolute wake normalization is fixed once from the 13.75 two-tracer
reference and never recalibrated between partitions.  This run therefore tests
only information recovery from finer tracer resolution.
"""
from __future__ import annotations
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_hod_masssplit_fisher as refmod
import wake_desi_multitracer_fisher as mt
from wake_desi_split_fixedcal_fisher import make_calibration
from bgs_full_mass_grid import GRID

BOUNDARIES=sorted(GRID)
CREF=np.array([0.10329109892348824,-0.060769580477018574,0.1691538145152833,
               0.4893798069085714,0.6861292155954458,0.4790989989937372,
               0.13123736993120932],float)


def reference_geometry():
    w.SURVEY=list(refmod.w.SURVEY)
    return w.prepare_survey(None)


def cumulative_profiles(refsurvey):
    """Convert validated cumulative n,b grid to mt.disjoint_tracers format."""
    vmap={round(float(r['z']),3):float(r['volume_hminus3_Gpc3']) for r in refsurvey}
    out={}
    for t,arr in GRID.items():
        prof=[]
        for z,nlo,blo,nhi,bhi in arr:
            V=vmap[round(float(z),3)]
            prof.append({'z':float(z),'zlo':float(z)-0.025,'zhi':float(z)+0.025,
                         'N_faint':float(nlo*V*1e9),'N_bright':float(nhi*V*1e9),
                         'b_faint':float(blo),'b_bright':float(bhi),
                         'n_faint_h3_Mpc3':float(nlo),'n_bright_h3_Mpc3':float(nhi),
                         'volume_hminus3_Gpc3':V,'delta_b':float(bhi-blo)})
        out[float(t)]=prof
    return out


def ref_score(P,shapes,f0,frac):
    c=CREF/np.linalg.norm(CREF)
    raw=c@shapes
    norm=1.0/float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    return math.sqrt(max(float(norm*norm*(c@P@c)),0.0))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--ntracer',type=int,required=True,choices=[5,6,7,8,9])
    ap.add_argument('--outdir',required=True)
    ap.add_argument('--mass',type=float,default=0.06)
    ap.add_argument('--z-match',type=float,default=1100.)
    ap.add_argument('--frac',type=float,default=0.30)
    ap.add_argument('--samples',type=int,default=50000)
    ap.add_argument('--validate',type=int,default=5)
    ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    refsurvey=reference_geometry(); prepared=cumulative_profiles(refsurvey)
    zs=np.array([float(r['z']) for r in refsurvey],float); base.ZBINS=zs
    q=np.linspace(0.,20.,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0); state0=base.build_state(p0,args.mass)
    cal=make_calibration(q,f0,state0,args.mass,refsurvey)

    configs=list(itertools.combinations(BOUNDARIES,args.ntracer-1))
    scan=[]; cache={}
    for boundaries in configs:
        survey=mt.disjoint_tracers(boundaries,prepared); cache[boundaries]=survey
        P=mt.linear_projected(q,f0,shapes,state0,args.mass,args.frac,survey,cal,'full')
        score=ref_score(P,shapes,f0,args.frac)
        scan.append({'boundaries':list(boundaries),'reference_direction_projected_SN':score})
        print('FULLGRID_SCAN',args.ntracer,json.dumps(scan[-1]))
    scan.sort(key=lambda r:r['reference_direction_projected_SN'],reverse=True)
    bestb=tuple(scan[0]['boundaries']); survey=cache[bestb]

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
             'validated_unprojected_hidden_SN':sn_un,
             'validated_after_wake_amplitude_projection_SN':sn_amp,
             'validated_after_full_odd_nuisance_projection_SN':sn_full,
             'max_relative_moment_mismatch':float(mismatch.max()),
             'per_redshift':pz,'coefficients':coeff.tolist()}
        vals.append(row); print('FULLGRID_VALIDATE',args.ntracer,json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    sn,row,fp,fm,shape=best
    np.savetxt(out/'best_pair.csv',np.column_stack([q,f0,fp,fm,shape]),delimiter=',',
               header='q,f_FD,f_plus,f_minus,normalized_delta_shape',comments='')
    summary={'mass_eV':args.mass,'ntracer':args.ntracer,'boundary_grid':BOUNDARIES,
             'number_of_partitions_scanned':len(configs),'best_boundaries':list(bestb),
             'configuration_scan':scan,'calibration_reference_split_log10':13.75,
             'calibration_rule':'One fixed physical per-z normalization from the 13.75 two-tracer reference.',
             'best':row,'all_validations':vals,
             'interpretation':'Full N x N covariance with all validated HOD thresholds. The scan tests saturation of information recovery as tracer mass resolution is increased.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('BEST_FULLGRID',args.ntracer,json.dumps(row))

if __name__=='__main__': main()
