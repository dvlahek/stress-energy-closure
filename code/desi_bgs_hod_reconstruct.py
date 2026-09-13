#!/usr/bin/env python3
"""Diagnose the DESI-BGS HOD selection used by Ge, Pasquini & Tan (2024).

Ge et al. use the public MXXL/hodpy lightcone HOD, SMT HMF and ST99 bias.
Our first reproduction used hodpy.kcorr.magnitude_faint(), which deliberately
chooses the faintest absolute-magnitude limit among extreme galaxy colours.
That is a completeness envelope, not the same as applying an observer-frame
r<19.5 flux cut to a representative colour population.

This diagnostic keeps r=19.5 fixed and evaluates several observer-frame g-r
colours, plus the hodpy completeness-envelope prescription. We compare the
resulting n(z), halo-mass-split biases and counts to the published BGS profile.
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

H0 = 67.36
H = H0/100.0
OMEGA_B = 0.02237/H**2
OMEGA_CDM = 0.1200/H**2
OMEGA_M = OMEGA_B + OMEGA_CDM + 0.00064/H**2
SIGMA8 = 0.8111
NS = 0.9649
COSMO = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Ob0=OMEGA_B, Tcmb0=2.7255)

LOGM_MIN = 11.0
LOGM_MAX = 16.0
LOGM_SPLIT = 13.75
AREA_DEG2 = 14000.0
ZGRID = np.arange(0.05, 0.401, 0.025)
R_LIMIT = 19.5
COLOURS = [0.0, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0]


class HODCosmologyWrapper:
    h0 = H
    OmegaM = OMEGA_M
    def comoving_distance(self, redshift):
        return np.asarray(COSMO.comoving_distance(redshift).value, float) * H


def integrate_logm(logm, y, lo, hi):
    m = (logm >= lo) & (logm <= hi)
    return float(simpson(y[m], x=np.log(10.0)*logm[m]))


def mass_background(z):
    mf = MassFunction(Mmin=LOGM_MIN, Mmax=LOGM_MAX, dlog10m=0.01, z=float(z),
                      hmf_model='SMT', cosmo_model=COSMO, sigma_8=SIGMA8, n=NS)
    logm = np.log10(np.asarray(mf.m, float))
    dndlnm = np.asarray(mf.dndlnm, float)
    nu = np.asarray(mf.nu, float)
    bh = np.asarray(ST99(nu=nu, delta_c=float(mf.delta_c), m=np.asarray(mf.m),
                         cosmo=COSMO, z=float(z)).bias(), float)
    dV = float(COSMO.differential_comoving_volume(z).value)*H**3*(np.pi/180.)**2
    return logm, dndlnm, bh, dV


def abs_limit(hod, z, prescription):
    zz=np.array([z],float)
    if prescription == 'envelope':
        return float(np.asarray(hod.kcorr.magnitude_faint(zz, R_LIMIT))[0])
    colour=float(prescription)
    return float(np.asarray(hod.kcorr.absolute_magnitude(np.array([R_LIMIT]), zz,
                                                         np.array([colour])))[0])


def compute_profile(hod, prescription, backgrounds):
    rows=[]
    for z,(logm,dndlnm,bh,dV) in zip(ZGRID,backgrounds):
        mag=abs_limit(hod,float(z),prescription)
        ngal=np.asarray(hod.number_galaxies_mean(logm,np.full_like(logm,mag),
                                                  np.full_like(logm,z)),float)
        w=dndlnm*ngal
        nlo=integrate_logm(logm,w,LOGM_MIN,LOGM_SPLIT)
        nhi=integrate_logm(logm,w,LOGM_SPLIT,LOGM_MAX)
        blo=integrate_logm(logm,w*bh,LOGM_MIN,LOGM_SPLIT)/max(nlo,1e-300)
        bhi=integrate_logm(logm,w*bh,LOGM_SPLIT,LOGM_MAX)/max(nhi,1e-300)
        rows.append({'z':float(z),'absolute_magnitude_limit':mag,
                     'n_low_h3_Mpc3':nlo,'n_high_h3_Mpc3':nhi,
                     'b_low':blo,'b_high':bhi,'delta_b':bhi-blo,
                     'dN_dz_deg2_low':nlo*dV,'dN_dz_deg2_high':nhi*dV,
                     'dN_dz_deg2_total':(nlo+nhi)*dV})
    return rows


def summarize(rows, label):
    z=np.array([r['z'] for r in rows])
    dV=np.array([float(COSMO.differential_comoving_volume(x).value)*H**3*(np.pi/180.)**2 for x in z])
    nlo=np.array([r['n_low_h3_Mpc3'] for r in rows]); nhi=np.array([r['n_high_h3_Mpc3'] for r in rows])
    dndz=np.array([r['dN_dz_deg2_total'] for r in rows])
    low=AREA_DEG2*float(simpson(nlo*dV,x=z)); high=AREA_DEG2*float(simpson(nhi*dV,x=z))
    return {'selection':label,
            'integrated_counts_0p05_0p4':{'low':low,'high':high,'total':low+high,
                                         'per_deg2_total':(low+high)/AREA_DEG2},
            'peak_dNdz_total_z':float(z[np.argmax(dndz)]),
            'peak_dNdz_total_deg2':float(np.max(dndz)),
            'bias_ranges':{'low_min':float(min(r['b_low'] for r in rows)),
                           'low_max':float(max(r['b_low'] for r in rows)),
                           'high_min':float(min(r['b_high'] for r in rows)),
                           'high_max':float(max(r['b_high'] for r in rows)),
                           'delta_b_min':float(min(r['delta_b'] for r in rows)),
                           'delta_b_max':float(max(r['delta_b'] for r in rows))},
            'rows':rows}


def main():
    out=Path('desi_bgs_hod_output'); out.mkdir(parents=True,exist_ok=True)
    hcosmo=HODCosmologyWrapper(); kcorr=GAMA_KCorrection(hcosmo)
    hod=HOD_BGS(mass_function=None,cosmology=hcosmo,kcorr=kcorr)
    backgrounds=[mass_background(float(z)) for z in ZGRID]

    prescriptions=[('envelope','envelope')]+[(f'colour_{c:.1f}',c) for c in COLOURS]
    summaries={}
    for label,prescription in prescriptions:
        rows=compute_profile(hod,prescription,backgrounds)
        s=summarize(rows,label); summaries[label]=s
        np.savetxt(out/f'bgs_profiles_{label}.csv',
                   np.array([[r[k] for k in ['z','n_low_h3_Mpc3','n_high_h3_Mpc3','b_low','b_high','delta_b','dN_dz_deg2_low','dN_dz_deg2_high','dN_dz_deg2_total']] for r in rows]),
                   delimiter=',',header='z,n_low,n_high,b_low,b_high,delta_b,dNdz_deg2_low,dNdz_deg2_high,dNdz_deg2_total',comments='')
        print('COLOUR_SWEEP',json.dumps({k:v for k,v in s.items() if k!='rows'}))

    result={'r_limit':R_LIMIT,
            'paper_targets':{'count_order':1.0e7,'peak_z_about':0.15,
                             'high_bias_gt':2.4,'delta_b_order':1.0},
            'cosmology':{'H0':H0,'Omega_m':OMEGA_M,'sigma8':SIGMA8,'ns':NS},
            'prescriptions':summaries}
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__': main()
