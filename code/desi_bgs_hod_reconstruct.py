#!/usr/bin/env python3
"""Sweep DESI-BGS halo-mass splits for neutrino-wake multi-tracer leverage.

Parent selection is fixed to the independently validated BGS-BRIGHT proxy:
r<19.5 and representative observer-frame g-r=1.0.  This selection gave
852.6 targets/deg^2 and dN/dz peaking at z=0.15, close to final DESI BGS.
We now vary only the halo-mass split.  The Ge, Pasquini & Tan choice
log10(M_split/[h^-1 Msun])=13.75 is one member of the sweep, but the ranking
here is based on the wake signal scaling rather than their standard focusing
likelihood.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from astropy.cosmology import FlatLambdaCDM
from hmf import MassFunction
from halomod.bias import ST99
from hodpy.hod_bgs import HOD_BGS
from hodpy.k_correction import GAMA_KCorrection

H0=67.36; H=H0/100.0
OMEGA_B=0.02237/H**2; OMEGA_CDM=0.1200/H**2
OMEGA_M=OMEGA_B+OMEGA_CDM+0.00064/H**2
SIGMA8=0.8111; NS=0.9649
COSMO=FlatLambdaCDM(H0=H0,Om0=OMEGA_M,Ob0=OMEGA_B,Tcmb0=2.7255)
LOGM_MIN=11.0; LOGM_MAX=16.0
LOGM_SPLITS=[12.75,13.0,13.25,13.5,13.75,14.0,14.25,14.5]
AREA_DEG2=14000.0; ZGRID=np.arange(0.05,0.401,0.025)
R_LIMIT=19.5; COLOUR=1.0; MASS_NU=0.06

class HODCosmologyWrapper:
    h0=H; OmegaM=OMEGA_M
    def comoving_distance(self,z): return np.asarray(COSMO.comoving_distance(z).value,float)*H

def integ(logm,y,lo,hi):
    m=(logm>=lo)&(logm<=hi)
    return float(simpson(y[m],x=np.log(10.0)*logm[m]))

def background(z):
    mf=MassFunction(Mmin=LOGM_MIN,Mmax=LOGM_MAX,dlog10m=0.01,z=float(z),
                    hmf_model='SMT',cosmo_model=COSMO,sigma_8=SIGMA8,n=NS)
    logm=np.log10(np.asarray(mf.m,float)); dnd=np.asarray(mf.dndlnm,float)
    nu=np.asarray(mf.nu,float)
    bh=np.asarray(ST99(nu=nu,delta_c=float(mf.delta_c),m=np.asarray(mf.m),
                       cosmo=COSMO,z=float(z)).bias(),float)
    dV=float(COSMO.differential_comoving_volume(z).value)*H**3*(np.pi/180.)**2
    return logm,dnd,bh,dV

def profile_for_split(hod,split,backs):
    rows=[]
    for z,(logm,dnd,bh,dV) in zip(ZGRID,backs):
        zz=np.array([z],float)
        mag=float(np.asarray(hod.kcorr.absolute_magnitude(np.array([R_LIMIT]),zz,
                                                          np.array([COLOUR])))[0])
        occ=np.asarray(hod.number_galaxies_mean(logm,np.full_like(logm,mag),
                                                 np.full_like(logm,z)),float)
        w=dnd*occ
        nlo=integ(logm,w,LOGM_MIN,split); nhi=integ(logm,w,split,LOGM_MAX)
        blo=integ(logm,w*bh,LOGM_MIN,split)/max(nlo,1e-300)
        bhi=integ(logm,w*bh,split,LOGM_MAX)/max(nhi,1e-300)
        rows.append({'z':float(z),'n_low':nlo,'n_high':nhi,'b_low':blo,'b_high':bhi,
                     'delta_b':bhi-blo,'dNdz_low':nlo*dV,'dNdz_high':nhi*dV,
                     'dNdz_total':(nlo+nhi)*dV})
    return rows

def summary(rows,split):
    z=np.array([r['z'] for r in rows]); dz=float(z[1]-z[0])
    dV=np.array([float(COSMO.differential_comoving_volume(x).value)*H**3*(np.pi/180.)**2 for x in z])
    nlo=np.array([r['n_low'] for r in rows]); nhi=np.array([r['n_high'] for r in rows])
    dndz=np.array([r['dNdz_total'] for r in rows]); db=np.array([r['delta_b'] for r in rows])
    low=AREA_DEG2*float(simpson(nlo*dV,x=z)); high=AREA_DEG2*float(simpson(nhi*dV,x=z))
    # Differential Okoli Eq.35 calibration diagnostic.  This is only for ranking
    # split thresholds before the full CLASS hidden-state Fisher.
    dN=AREA_DEG2*dndz*dz
    nth=1.7e7*(MASS_NU/0.05)**(-6.0)*(28.5**z)/(db**2)
    fd_sn=float(np.sqrt(np.sum(9.0*dN/nth)))
    return {'log10_Msplit':split,'Msplit_hminus1_Msun':10**split,
            'counts':{'low':low,'high':high,'total':low+high,'high_fraction':high/(low+high)},
            'peak_z':float(z[np.argmax(dndz)]),'peak_dNdz_deg2':float(np.max(dndz)),
            'bias':{'low_min':float(min(r['b_low'] for r in rows)),
                    'low_max':float(max(r['b_low'] for r in rows)),
                    'high_min':float(min(r['b_high'] for r in rows)),
                    'high_max':float(max(r['b_high'] for r in rows)),
                    'delta_b_min':float(db.min()),'delta_b_max':float(db.max())},
            'standard_FD_wake_SN_ranking_diagnostic':fd_sn,'rows':rows}

def main():
    out=Path('desi_bgs_hod_output'); out.mkdir(parents=True,exist_ok=True)
    hc=HODCosmologyWrapper(); hod=HOD_BGS(mass_function=None,cosmology=hc,kcorr=GAMA_KCorrection(hc))
    backs=[background(float(z)) for z in ZGRID]
    allsum={}
    for split in LOGM_SPLITS:
        rows=profile_for_split(hod,split,backs); s=summary(rows,split)
        key=f'{split:.2f}'; allsum[key]=s
        np.savetxt(out/f'bgs_masssplit_{key}.csv',
                   np.array([[r[k] for k in ['z','n_low','n_high','b_low','b_high','delta_b','dNdz_low','dNdz_high','dNdz_total']] for r in rows]),
                   delimiter=',',header='z,n_low,n_high,b_low,b_high,delta_b,dNdz_low,dNdz_high,dNdz_total',comments='')
        print('SPLIT',json.dumps({k:v for k,v in s.items() if k!='rows'}))
    best=max(allsum.values(),key=lambda s:s['standard_FD_wake_SN_ranking_diagnostic'])
    result={'selection':{'r_limit':R_LIMIT,'observer_g_minus_r':COLOUR,'area_deg2':AREA_DEG2,
                         'validated_parent_density_deg2':852.6064,'parent_peak_z':0.15},
            'mass_eV_for_ranking':MASS_NU,'splits':allsum,
            'best_by_standard_wake_scaling':{k:v for k,v in best.items() if k!='rows'}}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print('BEST_SPLIT',json.dumps(result['best_by_standard_wake_scaling']))
if __name__=='__main__': main()
