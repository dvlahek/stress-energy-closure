#!/usr/bin/env python3
"""Independent K01C central-orbit mathematics QA.

Reimplements the turning-point/angular integral without importing the K01C
runner, and compares the point-mass limit to the exact Rutherford/Kepler
scattering angle. Also verifies continuity of the frozen truncated-NFW
potential at r200.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as e27
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as e28

OUT=ROOT/"k01c_independent_math_qa.json"

def need(ok,msg):
    if not ok:
        raise RuntimeError("K01C_QA_STOP: "+msg)

def point_phi(r,M):
    return -e27.G*M/r

def orbit_coschi_point(v,b,M,nt=160):
    h=b*v
    lo=max(b*1e-14,1e-18)
    hi=b
    for _ in range(100):
        mid=.5*(lo+hi)
        R=v*v-2*point_phi(mid,M)-h*h/(mid*mid)
        if R>0:
            hi=mid
        else:
            lo=mid
    rmin=.5*(lo+hi)

    x,w=leggauss(nt)
    t=.5*(x+1)
    wt=.5*w
    u=1-t*t
    r=rmin/u
    rad=v*v-2*point_phi(r,M)-h*h/(r*r)
    need(float(np.min(rad))>-1e-7*max(v*v,1.0),"point radicand")
    rad=np.maximum(rad,np.finfo(float).tiny)
    integrand=2*t*(h/rmin)/np.sqrt(rad)
    phi0=float(np.sum(wt*integrand))
    chi=2*phi0-math.pi
    return math.cos(chi),chi,rmin

def exact_point(v,b,M):
    A=e27.G*M/(b*v*v)
    chi=2*math.atan(A)
    return math.cos(chi),chi

def nfw_phi(r,M,c,rs):
    r=np.asarray(r,dtype=float)
    r200=c*rs
    fc=math.log1p(c)-c/(1+c)
    out=np.empty_like(r)
    inside=r<r200
    if np.any(inside):
        x=r[inside]/rs
        out[inside]=-(e27.G*M/(rs*fc))*(np.log1p(x)/x-1/(1+c))
    if np.any(~inside):
        out[~inside]=-e27.G*M/r[~inside]
    return out

_,base,profiles=e28.sources()
_,_,_,halos=base
prof=profiles[0]
M=float(halos[0][0])
c=float(prof["c200c_population_median_DM14"])
rs=float(prof["r_s_comoving_Mpc_DM14"])/(1+e27.ZOBS)
r200=c*rs

speeds=[50.,100.,200.,500.,1000.,1800.]
bfacs=[.03,.1,.3,.7,1.,2.,4.]
rows=[]
max_cos_gap=0.0
max_chi_gap=0.0
for v in speeds:
    for f in bfacs:
        b=f*r200
        cn,chin,rmin=orbit_coschi_point(v,b,M)
        ce,chie=exact_point(v,b,M)
        cg=abs(cn-ce)
        ag=abs(chin-chie)
        max_cos_gap=max(max_cos_gap,cg)
        max_chi_gap=max(max_chi_gap,ag)
        rows.append({
            "v_km_s":v,"b_over_r200":f,
            "numeric_coschi":cn,"exact_coschi":ce,
            "abs_cos_gap":cg,"abs_chi_gap_rad":ag,
            "rmin_over_r200":rmin/r200,
        })

# Potential continuity at the truncation.
eps=1e-9
vals=nfw_phi(np.array([r200*(1-eps),r200,r200*(1+eps)]),M,c,rs)
target=-e27.G*M/r200
cont=max(abs(float(v)-target) for v in vals)/abs(target)

# Small-r finite NFW potential check against analytic x->0 limit.
phi0_expected=-(e27.G*M/(rs*(math.log1p(c)-c/(1+c))))*(1-1/(1+c))
phi_small=float(nfw_phi(np.array([rs*1e-10]),M,c,rs)[0])
center_rel=abs(phi_small-phi0_expected)/abs(phi0_expected)

need(max_cos_gap<2e-9,"point-mass coschi mismatch")
need(max_chi_gap<2e-9,"point-mass angle mismatch")
need(cont<2e-9,"NFW potential discontinuity")
need(center_rel<2e-9,"NFW central potential limit mismatch")

result={
    "stage":"K01C_INDEPENDENT_CENTRAL_ORBIT_MATH_QA",
    "status":"PASS_INDEPENDENT_POINT_MASS_AND_NFW_POTENTIAL_QA",
    "max_point_mass_abs_coschi_gap":max_cos_gap,
    "max_point_mass_abs_chi_gap_rad":max_chi_gap,
    "nfw_r200_potential_continuity_relative_gap":cont,
    "nfw_center_limit_relative_gap":center_rel,
    "cases":rows,
    "guards":{
        "does_not_import_K01C_runner":True,
        "observed_data_used":False,
        "survey_signal_to_noise_computed":False
    }
}
OUT.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
print("K01C_QA_MAX_POINT_COS_GAP",max_cos_gap,flush=True)
print("K01C_QA_MAX_POINT_CHI_GAP",max_chi_gap,flush=True)
print("K01C_QA_NFW_R200_CONTINUITY",cont,flush=True)
print("K01C_QA_NFW_CENTER_LIMIT",center_rel,flush=True)
print("K01C_INDEPENDENT_CENTRAL_ORBIT_MATH_QA_PASS",flush=True)
