#!/usr/bin/env python3
"""Survey-design map for the DESI-BGS hidden-state wake observable.

Uses the validated five-tracer halo-mass partition and the stable leading
stress-energy-matched kinetic direction.  The physical wake normalization is
fixed to the same 13.75 reference calibration used in the production Fisher.
We then vary only survey area/volume and tracer number density while keeping
biases and the redshift distribution shape fixed.

Area factor multiplies N at fixed n (larger volume).  Density factor multiplies
both N and n at fixed volume (deeper spectroscopy).  No recalibration is
performed, so this is a genuine survey-design scaling test.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_multitracer_fisher as mt
from wake_desi_split_fixedcal_fisher import make_calibration

OUT=Path('wake_survey_design_output')
CREF=np.array([0.10329109892348824,-0.060769580477018574,0.1691538145152833,
               0.4893798069085714,0.6861292155954458,0.4790989989937372,
               0.13123736993120932],float)
AREA_FACTORS=[1.0,1.25,1.5,2.0,41253.0/14000.0]
DENSITY_FACTORS=[1.0,1.25,1.5,2.0,3.0,5.0,10.0]


def scaled_survey(survey,area,density):
    out={}
    for z,trs in survey.items():
        out[z]=[]
        for t in trs:
            q=dict(t)
            q['N']=float(t['N']*area*density)
            q['n']=float(t['n']*density)
            out[z].append(q)
    return out


def main():
    OUT.mkdir(exist_ok=True)
    mass=0.06; zmatch=1100.; frac=0.30
    prepared={t:mt.prepare_profile(mt.CUM[t]) for t in mt.BOUNDARIES}
    survey=mt.disjoint_tracers(tuple(mt.BOUNDARIES),prepared)
    zs=np.array(sorted(survey),float); base.ZBINS=zs
    q=np.linspace(0.,20.,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,mass,zmatch)
    c=CREF/np.linalg.norm(CREF); raw=c@shapes
    norm=1.0/float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    shape=raw*norm; fp=f0+frac*shape; fm=f0-frac*shape
    p0=OUT/'fd_reference.dat'; pp=OUT/'f_plus.dat'; pm=OUT/'f_minus.dat'
    base.write_psd(p0,q,f0); base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
    state0=base.build_state(p0,mass); statep=base.build_state(pp,mass); statem=base.build_state(pm,mass)
    w.SURVEY=list(mt.CUM[13.75]); refsurvey=w.prepare_survey(state0)
    cal=make_calibration(q,f0,state0,mass,refsurvey)

    rows=[]
    for A in AREA_FACTORS:
        for D in DENSITY_FACTORS:
            s=scaled_survey(survey,A,D)
            un,full,_=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,s,cal,'full')
            _,amp,_=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,s,cal,'amp')
            r={'area_factor':A,'density_factor':D,'unprojected_SN':un,
               'amplitude_projected_SN':amp,'full_projected_SN':full}
            rows.append(r); print('DESIGN',json.dumps(r))

    # Bisection: density factor needed for S/N=1 at current DESI area.
    lo,hi=1.0,20.0
    for _ in range(24):
        mid=0.5*(lo+hi)
        s=scaled_survey(survey,1.0,mid)
        _,sn,_=mt.nonlinear_projected(q,f0,fp,fm,state0,statep,statem,mass,s,cal,'full')
        if sn>=1.0: hi=mid
        else: lo=mid
    density_1sigma=hi
    # Area scaling is exactly sqrt(A) at fixed n. Infer exact requirements from baseline.
    baseline=next(r for r in rows if r['area_factor']==1.0 and r['density_factor']==1.0)
    area_1sigma=(1.0/baseline['full_projected_SN'])**2
    area_2sigma=(2.0/baseline['full_projected_SN'])**2
    summary={'mass_eV':mass,'partition_boundaries':mt.BOUNDARIES,
             'baseline_area_deg2':14000.0,'baseline':baseline,'grid':rows,
             'density_factor_for_1sigma_at_14000deg2':density_1sigma,
             'area_factor_for_1sigma_at_fixed_density':area_1sigma,
             'area_deg2_for_1sigma_at_fixed_density':14000.0*area_1sigma,
             'area_factor_for_2sigma_at_fixed_density':area_2sigma,
             'full_sky_SN_at_fixed_density':baseline['full_projected_SN']*math.sqrt(41253.0/14000.0),
             'calibration_rule':'Physical wake normalization fixed once to 13.75 reference. Area and density changes do not recalibrate the signal.'}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('SUMMARY',json.dumps(summary))

if __name__=='__main__': main()
