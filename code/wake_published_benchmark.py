#!/usr/bin/env python3
"""Published wake benchmark comparison.

This script is a validation table, not a new forecast.  It checks two external
anchors used by the production wake pipeline:

1. Okoli et al. (2017) Eq. 35.  The code's absolute standard-wake calibration
   must reproduce the published galaxy-count scaling exactly.  This row is
   explicitly labelled a calibration identity, not an independent validation.
2. Ge, Pasquini & Tan (2024/2025) DESI-BGS halo-mass split.  Their published
   likelihood rises from -2 ln L = 2.6 at log10 Msplit=12.75 to a maximum 3.9
   at 13.75.  Here the same physical FD wake normalization is held fixed across
   split thresholds.  We compare the *shape* of our standard-FD two-tracer
   information curve after normalizing the 13.75 point to 3.9.  No hidden-state
   optimization enters this test.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as bonvin
import wake_desi_multitracer_fisher as mt
import wake_desi_multitracer_fullgrid as fg

PUBLISHED_GE={12.75:2.6,13.75:3.9}


def fd_raw_information(q,f0,state0,mass,survey):
    total=0.0
    for z in sorted(survey):
        info=state0['z'][z]; trs=survey[z]
        V=sum(t['N'] for t in trs)/sum(t['n'] for t in trs)/1e9
        rows=mt.mode_rows(info,V,trs); cache={}
        for ik,im,w in rows:
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,f0,f0,state0,state0,state0,mass,z,ik)[0]
            total += w*cache[ik]
    return float(total)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--outdir',default='wake_published_benchmark')
    ap.add_argument('--mass',type=float,default=0.06); ap.add_argument('--z-match',type=float,default=1100.0)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    refsurvey=fg.reference_geometry(); prepared=fg.cumulative_profiles(refsurvey)
    zs=np.array([float(r['z']) for r in refsurvey],float); base.ZBINS=zs
    q=np.linspace(0.0,20.0,4000)
    f0,_,_,_,_,_,_=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/'fd_reference.dat'; base.write_psd(p0,q,f0); state0=base.build_state(p0,args.mass)

    raw={}
    for t in sorted(prepared):
        survey=mt.disjoint_tracers((float(t),),prepared)
        raw[float(t)]=fd_raw_information(q,f0,state0,args.mass,survey)
    ref=raw[13.75]
    scaled={t:3.9*v/ref for t,v in raw.items()}
    best=max(scaled,key=scaled.get)

    # Okoli Eq.35 identity at a few published-style points.
    okoli=[]
    for mass,z,N,db in [(0.05,0.0,1.7e7,1.0),(0.07,0.0,2.0e6,1.0),(0.10,0.2,1.0e7,1.0)]:
        nth=1.7e7*(mass/0.05)**(-6.0)*(28.5**z)/(db**2)
        published_formula_sn=3.0*math.sqrt(N/nth)
        code_sn=math.sqrt(bonvin.desired_fd_sn2(mass,z,N,db))
        okoli.append({'mass_eV':mass,'z':z,'N_gal':N,'delta_b':db,
                      'published_Eq35_SN':published_formula_sn,'code_calibration_SN':code_sn,
                      'relative_difference':abs(code_sn-published_formula_sn)/max(abs(published_formula_sn),1e-300)})

    rows=[]
    for t in sorted(scaled):
        pub=PUBLISHED_GE.get(t,None)
        val=scaled[t]
        rows.append({'log10_Msplit_hinvMsun':t,'published_minus2logL':pub,
                     'our_shape_matched_minus2logL':val,
                     'relative_difference_if_published':None if pub is None else abs(val-pub)/pub})
    ge_1275=next(r for r in rows if r['log10_Msplit_hinvMsun']==12.75)
    summary={'mass_eV_used_for_shape_control':args.mass,
             'z_match_used_for_class_fd_state':args.z_match,
             'ge_published_anchor':{'12.75':2.6,'13.75':3.9,'published_best_threshold':13.75},
             'our_best_threshold':float(best),'our_scaled_curve':rows,
             'relative_difference_at_12.75':ge_1275['relative_difference_if_published'],
             'okoli_eq35_calibration_checks':okoli,
             'interpretation':'Ge comparison tests threshold-shape reproduction with one fixed physical normalization. The curve is normalized once at 13.75, so only the remaining threshold dependence is a validation. Okoli rows are calibration identities and are labelled as such.'}
    with open(out/'benchmark_table.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['log10_Msplit_hinvMsun','published_minus2logL','our_shape_matched_minus2logL','relative_difference_if_published']); w.writeheader(); w.writerows(rows)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    if abs(best-13.75)>1e-12:
        raise RuntimeError('Published optimum not reproduced: '+json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
