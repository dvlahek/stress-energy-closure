#!/usr/bin/env python3
"""K01C-D1 arbitrary-precision orbit diagnostic.

Diagnostic only. It does not modify or reclassify the parent K01C PASS or the
formal K01C-R1 FAIL. It checks the central-orbit integral at 80-digit precision
and compares representative NFW orbits with the existing double-precision
Gauss-Legendre implementation.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as e27
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as e28

DECL=ROOT/"source_data/k01c_d1_arbitrary_precision_orbit_diagnostic_2026-11-05.json"
OUT=ROOT/"k01c_d1_arbitrary_precision_orbit_diagnostic_result.json"

mp.mp.dps=80

def need(ok,msg):
    if not ok:
        raise RuntimeError("K01C_D1_STOP: "+msg)

def setup():
    _,base,profiles=e28.sources()
    _,_,_,halos=base
    p=profiles[0]
    M=float(halos[0][0])
    c=float(p["c200c_population_median_DM14"])
    rs=float(p["r_s_comoving_Mpc_DM14"])/(1+e27.ZOBS)
    r200=c*rs
    return M,c,rs,r200

def mp_point_phi(r,M):
    return -mp.mpf(str(e27.G))*M/r

def mp_nfw_mfrac(r,c,rs):
    x=r/rs
    fc=mp.log(1+c)-c/(1+c)
    if x<c:
        return (mp.log(1+x)-x/(1+x))/fc
    return mp.mpf(1)

def mp_nfw_phi(r,M,c,rs):
    G=mp.mpf(str(e27.G))
    r200=c*rs
    fc=mp.log(1+c)-c/(1+c)
    if r<r200:
        x=r/rs
        return -(G*M/(rs*fc))*(mp.log(1+x)/x-1/(1+c))
    return -G*M/r

def bisect_root(v,b,M,potential,c=None,rs=None):
    h=b*v
    lo=b*mp.mpf("1e-30")
    hi=b
    def R(r):
        phi=potential(r,M) if c is None else potential(r,M,c,rs)
        return v*v-2*phi-h*h/(r*r)
    need(R(lo)<0 and R(hi)>0,"root bracket")
    for _ in range(260):
        mid=(lo+hi)/2
        if R(mid)>0:
            hi=mid
        else:
            lo=mid
    return (lo+hi)/2

def mp_orbit_chi(v,b,M,potential,c=None,rs=None):
    v=mp.mpf(str(v)); b=mp.mpf(str(b)); M=mp.mpf(str(M))
    if c is not None:
        c=mp.mpf(str(c)); rs=mp.mpf(str(rs))
    h=b*v
    rmin=bisect_root(v,b,M,potential,c,rs)

    def phi(r):
        return potential(r,M) if c is None else potential(r,M,c,rs)

    def mfrac(r):
        if c is None:
            return mp.mpf(1)
        return mp_nfw_mfrac(r,c,rs)

    G=mp.mpf(str(e27.G))
    # R'(rmin) = -2 Phi'(rmin)+2h^2/rmin^3,
    # with Phi'=G M(<r)/r^2.
    Rp=-2*G*M*mfrac(rmin)/(rmin*rmin)+2*h*h/(rmin**3)
    lim=2*h/(rmin*mp.sqrt(Rp*rmin))

    def integrand(t):
        if t == 0:
            return lim
        u=1-t*t
        r=rmin/u
        rad=v*v-2*phi(r)-h*h/(r*r)
        if rad<0 and abs(rad)<mp.mpf("1e-60")*max(v*v,1):
            rad=mp.mpf(0)
        return 2*t*(h/rmin)/mp.sqrt(rad)

    phi0=mp.quad(integrand,[0,mp.mpf("0.05"),mp.mpf("0.2"),mp.mpf("0.5"),mp.mpf("0.8"),mp.mpf("0.95"),1])
    chi=2*phi0-mp.pi
    return chi,rmin

def double_nfw_chi(v,b,M,c,rs,nt):
    h=b*v
    def phi(arr):
        arr=np.asarray(arr,dtype=float)
        r200=c*rs
        fc=math.log1p(c)-c/(1+c)
        out=np.empty_like(arr)
        inside=arr<r200
        if np.any(inside):
            x=arr[inside]/rs
            out[inside]=-(e27.G*M/(rs*fc))*(np.log1p(x)/x-1/(1+c))
        if np.any(~inside):
            out[~inside]=-e27.G*M/arr[~inside]
        return out

    lo=max(b*1e-14,c*rs*1e-16)
    hi=b
    for _ in range(100):
        mid=.5*(lo+hi)
        R=v*v-2*float(phi(np.array([mid]))[0])-h*h/(mid*mid)
        if R>0:
            hi=mid
        else:
            lo=mid
    rmin=.5*(lo+hi)

    x,w=leggauss(nt)
    t=.5*(x+1); wt=.5*w
    u=1-t*t
    r=rmin/u
    rad=v*v-2*phi(r)-h*h/(r*r)
    rad=np.maximum(rad,np.finfo(float).tiny)
    phi0=float(np.sum(wt*(2*t*(h/rmin)/np.sqrt(rad))))
    return 2*phi0-math.pi

decl=json.loads(DECL.read_text())
M,c,rs,r200=setup()

point=[]
max_point_chi_gap=mp.mpf(0)
max_point_cos_gap=mp.mpf(0)
for case in decl["fixed_point_mass_cases"]:
    v=case["v_km_s"]; b=case["b_over_r200"]*r200
    chi,_=mp_orbit_chi(v,b,M,mp_point_phi)
    exact=2*mp.atan(mp.mpf(str(e27.G))*mp.mpf(str(M))/(mp.mpf(str(b))*mp.mpf(str(v))**2))
    ag=abs(chi-exact)
    cg=abs(mp.cos(chi)-mp.cos(exact))
    max_point_chi_gap=max(max_point_chi_gap,ag)
    max_point_cos_gap=max(max_point_cos_gap,cg)
    point.append({
        "v_km_s":v,
        "b_over_r200":case["b_over_r200"],
        "mp_chi_rad":mp.nstr(chi,40),
        "exact_chi_rad":mp.nstr(exact,40),
        "abs_chi_gap":mp.nstr(ag,20),
        "abs_coschi_gap":mp.nstr(cg,20),
    })

nfw=[]
max_64=0.0
max_160=0.0
for case in decl["fixed_nfw_spot_cases"]:
    v=case["v_km_s"]; b=case["b_over_r200"]*r200
    chi_mp,_=mp_orbit_chi(v,b,M,mp_nfw_phi,c,rs)
    ref=float(chi_mp)
    chi64=double_nfw_chi(v,b,M,c,rs,64)
    chi160=double_nfw_chi(v,b,M,c,rs,160)
    g64=abs(chi64-ref)
    g160=abs(chi160-ref)
    max_64=max(max_64,g64)
    max_160=max(max_160,g160)
    nfw.append({
        "v_km_s":v,
        "b_over_r200":case["b_over_r200"],
        "mp_chi_rad":mp.nstr(chi_mp,40),
        "double64_chi_rad":chi64,
        "double160_chi_rad":chi160,
        "abs_gap_64_rad":g64,
        "abs_gap_160_rad":g160,
    })

result={
    "stage":decl["stage"],
    "status":"DIAGNOSTIC_COMPLETE_NO_SCIENCE_RECLASSIFICATION",
    "point_mass":{
        "max_abs_chi_gap":mp.nstr(max_point_chi_gap,20),
        "max_abs_coschi_gap":mp.nstr(max_point_cos_gap,20),
        "cases":point,
    },
    "nfw_spot_checks":{
        "max_abs_gap_64_rad":max_64,
        "max_abs_gap_160_rad":max_160,
        "cases":nfw,
    },
    "interpretation_guard":decl["guard"],
}
OUT.write_text(json.dumps(result,indent=2)+"\n")

print("K01C_D1_POINT_MAX_CHI_GAP",mp.nstr(max_point_chi_gap,20),flush=True)
print("K01C_D1_POINT_MAX_COS_GAP",mp.nstr(max_point_cos_gap,20),flush=True)
print("K01C_D1_NFW_MAX_GAP_64",max_64,flush=True)
print("K01C_D1_NFW_MAX_GAP_160",max_160,flush=True)
print("K01C_D1_DIAGNOSTIC_COMPLETE",flush=True)
