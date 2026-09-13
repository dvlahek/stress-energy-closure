#!/usr/bin/env python3
"""Reconstruct the DESI-BGS halo-mass split used in Ge, Pasquini & Tan (2024).

The calculation follows their Eqs. (4.1)-(4.2):
  n_g(z) = int dlnM [dn/dlnM] <N(M)>
  b_g(z) = n_g^-1 int dlnM [dn/dlnM] <N(M)> b_h(M,z)

Inputs follow the paper as closely as possible with public software:
  * BGS HOD from amjsmith/hodpy (MXXL lightcone HOD)
  * SMT halo mass function from hmf
  * ST99 peak-background-split halo bias from halomod
  * Planck-2018-like flat LCDM parameters
  * M_min=1e11, M_max=1e16 h^-1 Msun
  * benchmark split M_split=1e13.75 h^-1 Msun

The observable BGS absolute-magnitude threshold is obtained from the public
HOD_BGS k-correction and its default apparent-magnitude limit, preserving the
lightcone catalogue model used by the cited paper.
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

H0 = 67.36
OMEGA_B = 0.02237/(H0/100.)**2
OMEGA_CDM = 0.1200/(H0/100.)**2
OMEGA_M = OMEGA_B + OMEGA_CDM + 0.00064/(H0/100.)**2
SIGMA8 = 0.8111
NS = 0.9649
COSMO = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Ob0=OMEGA_B, Tcmb0=2.7255)

LOGM_MIN = 11.0
LOGM_MAX = 16.0
LOGM_SPLIT = 13.75
AREA_DEG2 = 14000.0
ZGRID = np.arange(0.05, 0.401, 0.025)


def integrate_logm(logm, y, lo, hi):
    m = (logm >= lo) & (logm <= hi)
    return float(simpson(y[m], x=np.log(10.0)*logm[m]))


def main():
    out = Path('desi_bgs_hod_output')
    out.mkdir(parents=True, exist_ok=True)

    hod = HOD_BGS()
    rows = []

    for z in ZGRID:
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

        # Public lightcone HOD applies its apparent-magnitude threshold through
        # the GAMA k-correction to obtain the redshift-dependent absolute cut.
        zz = np.array([z], dtype=float)
        maglim = float(np.asarray(hod.kcorr.magnitude_faint(zz, hod.mag_faint))[0])
        ngal_halo = np.asarray(hod.number_galaxies_mean(logm, np.full_like(logm, maglim), np.full_like(logm, z)), float)

        weight = dndlnm * ngal_halo
        nlo = integrate_logm(logm, weight, LOGM_MIN, LOGM_SPLIT)
        nhi = integrate_logm(logm, weight, LOGM_SPLIT, LOGM_MAX)
        blo = integrate_logm(logm, weight*bh, LOGM_MIN, LOGM_SPLIT)/max(nlo,1e-300)
        bhi = integrate_logm(logm, weight*bh, LOGM_SPLIT, LOGM_MAX)/max(nhi,1e-300)
        ntot = nlo+nhi
        btot = (nlo*blo+nhi*bhi)/max(ntot,1e-300)

        # Number per deg^2 per dz. Astropy gives Mpc^3/sr. HMF densities are
        # h^3/Mpc^3. Convert volume to (Mpc/h)^3 before multiplying.
        dV_dz_sr_mpc3 = float(COSMO.differential_comoving_volume(z).value)
        h = H0/100.0
        dV_dz_deg2_hm3 = dV_dz_sr_mpc3 * h**3 * (np.pi/180.0)**2

        row = {
            'z': float(z),
            'absolute_magnitude_limit': maglim,
            'n_low_h3_Mpc3': nlo,
            'n_high_h3_Mpc3': nhi,
            'n_total_h3_Mpc3': ntot,
            'b_low': blo,
            'b_high': bhi,
            'b_total': btot,
            'delta_b': bhi-blo,
            'dN_dz_deg2_low': nlo*dV_dz_deg2_hm3,
            'dN_dz_deg2_high': nhi*dV_dz_deg2_hm3,
            'dN_dz_deg2_total': ntot*dV_dz_deg2_hm3,
        }
        rows.append(row)
        print(json.dumps(row))

    z = np.array([r['z'] for r in rows])
    nlo = np.array([r['n_low_h3_Mpc3'] for r in rows])
    nhi = np.array([r['n_high_h3_Mpc3'] for r in rows])
    dV = np.array([float(COSMO.differential_comoving_volume(x).value)*(H0/100.)**3*(np.pi/180.)**2 for x in z])
    total_low = AREA_DEG2*float(simpson(nlo*dV, x=z))
    total_high = AREA_DEG2*float(simpson(nhi*dV, x=z))

    summary = {
        'paper_benchmark': {
            'M_min_hminus1_Msun': 1e11,
            'M_split_hminus1_Msun': 10**LOGM_SPLIT,
            'M_max_hminus1_Msun': 1e16,
            'area_deg2': AREA_DEG2,
        },
        'cosmology': {'H0':H0,'Omega_m':OMEGA_M,'Omega_b':OMEGA_B,'sigma8':SIGMA8,'ns':NS},
        'hod_apparent_magnitude_limit': float(hod.mag_faint),
        'integrated_counts_0p05_0p4': {
            'low': total_low,
            'high': total_high,
            'total': total_low+total_high,
        },
        'peak_dNdz_total_z': float(z[np.argmax([r['dN_dz_deg2_total'] for r in rows])]),
        'bias_ranges': {
            'low_min': float(min(r['b_low'] for r in rows)),
            'low_max': float(max(r['b_low'] for r in rows)),
            'high_min': float(min(r['b_high'] for r in rows)),
            'high_max': float(max(r['b_high'] for r in rows)),
            'delta_b_min': float(min(r['delta_b'] for r in rows)),
            'delta_b_max': float(max(r['delta_b'] for r in rows)),
        },
        'paper_qualitative_targets': {
            'peak_redshift_about': 0.15,
            'low_bias_between': [1.0,2.0],
            'high_bias_greater_than': 2.4,
            'delta_b_about': 1.0,
        },
        'rows': rows,
    }
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    np.savetxt(out/'bgs_profiles.csv', np.array([[r[k] for k in ['z','n_low_h3_Mpc3','n_high_h3_Mpc3','b_low','b_high','delta_b','dN_dz_deg2_low','dN_dz_deg2_high']] for r in rows]), delimiter=',', header='z,n_low,n_high,b_low,b_high,delta_b,dNdz_deg2_low,dNdz_deg2_high', comments='')
    print('SUMMARY')
    print(json.dumps(summary, indent=2))

if __name__=='__main__':
    main()
