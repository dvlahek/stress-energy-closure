#!/usr/bin/env python3
"""Reconstruct the DESI-BGS halo-mass split used in Ge, Pasquini & Tan (2024).

The calculation follows their Eqs. (4.1)-(4.2):
  n_g(z) = int dlnM [dn/dlnM] <N(M)>
  b_g(z) = n_g^-1 int dlnM [dn/dlnM] <N(M)> b_h(M,z)

This version performs a compact apparent-magnitude sweep to identify the
observable BGS selection corresponding to Fig. 1 of Ge, Pasquini & Tan. The
published paper describes about 10 million bright galaxies over 14,000 deg^2,
with dN/dz peaking near z~0.15. DESI BGS Bright is the r<19.5 sample.

Inputs:
  * BGS HOD, target LF, slide factors and k-correction from amjsmith/hodpy
  * SMT halo mass function from hmf
  * ST99 peak-background-split halo bias from halomod
  * Planck-2018-like flat LCDM parameters
  * M_min=1e11, M_max=1e16 h^-1 Msun
  * benchmark split M_split=1e13.75 h^-1 Msun
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
MAG_LIMITS = [19.3, 19.5, 19.7, 20.0, 20.2]


class HODCosmologyWrapper:
    """Only the Mpc/h comoving-distance method required by hodpy k-correction."""
    h0 = H
    OmegaM = OMEGA_M
    def comoving_distance(self, redshift):
        return np.asarray(COSMO.comoving_distance(redshift).value, float) * H


def integrate_logm(logm, y, lo, hi):
    m = (logm >= lo) & (logm <= hi)
    return float(simpson(y[m], x=np.log(10.0)*logm[m]))


def compute_mass_background(z):
    mf = MassFunction(
        Mmin=LOGM_MIN,
        Mmax=LOGM_MAX,
        dlog10m=0.01,
        z=float(z),
        hmf_model='SMT',
        cosmo_model=COSMO,
        sigma_8=SIGMA8,
        n=NS,
    )
    logm = np.log10(np.asarray(mf.m, float))
    dndlnm = np.asarray(mf.dndlnm, float)
    nu = np.asarray(mf.nu, float)
    bh = np.asarray(ST99(nu=nu, delta_c=float(mf.delta_c), m=np.asarray(mf.m), cosmo=COSMO, z=float(z)).bias(), float)
    dV_dz_sr_mpc3 = float(COSMO.differential_comoving_volume(z).value)
    dV_dz_deg2_hm3 = dV_dz_sr_mpc3 * H**3 * (np.pi/180.0)**2
    return logm, dndlnm, bh, dV_dz_deg2_hm3


def compute_profile(hod, mag_limit):
    rows = []
    for z in ZGRID:
        logm, dndlnm, bh, dV = compute_mass_background(z)
        zz = np.array([z], dtype=float)
        maglim_abs = float(np.asarray(hod.kcorr.magnitude_faint(zz, mag_limit))[0])
        ngal_halo = np.asarray(
            hod.number_galaxies_mean(logm, np.full_like(logm, maglim_abs), np.full_like(logm, z)),
            float,
        )
        weight = dndlnm * ngal_halo
        nlo = integrate_logm(logm, weight, LOGM_MIN, LOGM_SPLIT)
        nhi = integrate_logm(logm, weight, LOGM_SPLIT, LOGM_MAX)
        blo = integrate_logm(logm, weight*bh, LOGM_MIN, LOGM_SPLIT)/max(nlo,1e-300)
        bhi = integrate_logm(logm, weight*bh, LOGM_SPLIT, LOGM_MAX)/max(nhi,1e-300)
        ntot = nlo+nhi
        btot = (nlo*blo+nhi*bhi)/max(ntot,1e-300)
        rows.append({
            'z': float(z),
            'absolute_magnitude_limit': maglim_abs,
            'n_low_h3_Mpc3': nlo,
            'n_high_h3_Mpc3': nhi,
            'n_total_h3_Mpc3': ntot,
            'b_low': blo,
            'b_high': bhi,
            'b_total': btot,
            'delta_b': bhi-blo,
            'dN_dz_deg2_low': nlo*dV,
            'dN_dz_deg2_high': nhi*dV,
            'dN_dz_deg2_total': ntot*dV,
        })
    return rows


def summarize(rows, mag_limit):
    z = np.array([r['z'] for r in rows])
    dV = np.array([float(COSMO.differential_comoving_volume(x).value)*H**3*(np.pi/180.)**2 for x in z])
    nlo = np.array([r['n_low_h3_Mpc3'] for r in rows])
    nhi = np.array([r['n_high_h3_Mpc3'] for r in rows])
    dndz = np.array([r['dN_dz_deg2_total'] for r in rows])
    total_low = AREA_DEG2*float(simpson(nlo*dV, x=z))
    total_high = AREA_DEG2*float(simpson(nhi*dV, x=z))
    return {
        'apparent_magnitude_limit': mag_limit,
        'integrated_counts_0p05_0p4': {
            'low': total_low,
            'high': total_high,
            'total': total_low+total_high,
            'per_deg2_total': (total_low+total_high)/AREA_DEG2,
        },
        'peak_dNdz_total_z': float(z[np.argmax(dndz)]),
        'peak_dNdz_total_deg2': float(np.max(dndz)),
        'bias_ranges': {
            'low_min': float(min(r['b_low'] for r in rows)),
            'low_max': float(max(r['b_low'] for r in rows)),
            'high_min': float(min(r['b_high'] for r in rows)),
            'high_max': float(max(r['b_high'] for r in rows)),
            'delta_b_min': float(min(r['delta_b'] for r in rows)),
            'delta_b_max': float(max(r['delta_b'] for r in rows)),
        },
        'rows': rows,
    }


def main():
    out = Path('desi_bgs_hod_output')
    out.mkdir(parents=True, exist_ok=True)

    hcosmo = HODCosmologyWrapper()
    kcorr = GAMA_KCorrection(hcosmo)
    hod = HOD_BGS(mass_function=None, cosmology=hcosmo, kcorr=kcorr)

    all_summaries = {}
    for mag_limit in MAG_LIMITS:
        rows = compute_profile(hod, mag_limit)
        summary = summarize(rows, mag_limit)
        key = f'r{mag_limit:.1f}'
        all_summaries[key] = summary
        np.savetxt(
            out/f'bgs_profiles_{key}.csv',
            np.array([[r[k] for k in ['z','n_low_h3_Mpc3','n_high_h3_Mpc3','b_low','b_high','delta_b','dN_dz_deg2_low','dN_dz_deg2_high','dN_dz_deg2_total']] for r in rows]),
            delimiter=',',
            header='z,n_low,n_high,b_low,b_high,delta_b,dNdz_deg2_low,dNdz_deg2_high,dNdz_deg2_total',
            comments='',
        )
        compact = {k:v for k,v in summary.items() if k != 'rows'}
        print('MAG_SWEEP', json.dumps(compact))

    result = {
        'paper_benchmark': {
            'M_min_hminus1_Msun': 1e11,
            'M_split_hminus1_Msun': 10**LOGM_SPLIT,
            'M_max_hminus1_Msun': 1e16,
            'area_deg2': AREA_DEG2,
            'published_peak_redshift_about': 0.15,
            'published_bgs_count_order': 1.0e7,
            'published_high_bias_greater_than': 2.4,
            'published_delta_b_order': 1.0,
        },
        'cosmology': {'H0':H0,'Omega_m':OMEGA_M,'Omega_b':OMEGA_B,'sigma8':SIGMA8,'ns':NS},
        'magnitude_sweep': all_summaries,
    }
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print('SUMMARY')
    print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk!='rows'} for k,v in all_summaries.items()}, indent=2))

if __name__=='__main__':
    main()
