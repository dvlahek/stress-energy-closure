#!/usr/bin/env python3
"""Unit-corrected wrapper for the survey-aware CMB likelihood forecast.

The CLASS temperature/polarization spectra used by the forecast are dimensionless
D_ell spectra in DeltaT/T units.  Map depths are specified in microkelvin-arcmin,
so the white-noise amplitude must first be converted to DeltaT/T before the beam
and D_ell factors are applied.

This wrapper leaves the validated CLASS spectra and Fisher implementation unchanged
and replaces only the map-noise conversion.  It is intended to re-summarize the
already-computed survey spectra from workflow run 34677322082; no CLASS rerun is
required for this correction.
"""
from __future__ import annotations

import math

import numpy as np

import class_cmb_survey_likelihood_forecast as forecast
from class_response_optimize import TCMB_K


def noise_cl_dimensionless(ell: np.ndarray, exp: dict, pol: bool) -> np.ndarray:
    """Beam-deconvolved white map noise in dimensionless C_ell units.

    Input depths are microkelvin-arcmin.  CLASS scalar TT/TE/EE outputs used by
    this forecast are in dimensionless temperature units, so convert

        microK -> K -> DeltaT/T_CMB

    before the usual arcmin-to-radian and Gaussian-beam factors.
    """
    inv = np.zeros_like(ell, dtype=float)
    for ch in exp["channels"]:
        depth_uK_arcmin = float(
            ch["noise_P_uK_arcmin"] if pol else ch["noise_T_uK_arcmin"]
        )
        fwhm_arcmin = float(ch["beam_arcmin"])

        delta_dimensionless_rad = (
            depth_uK_arcmin
            * 1.0e-6
            / TCMB_K
            * math.pi
            / (180.0 * 60.0)
        )
        sigma_b = (
            fwhm_arcmin
            * math.pi
            / (180.0 * 60.0)
            / math.sqrt(8.0 * math.log(2.0))
        )
        ncl = (
            delta_dimensionless_rad**2
            * np.exp(ell * (ell + 1.0) * sigma_b * sigma_b)
        )
        inv += np.where(ncl > 0, 1.0 / ncl, 0.0)

    return np.where(inv > 0, 1.0 / inv, np.inf)


# The original Fisher code resolves this global at call time.  Only the unit
# conversion is replaced; all covariance, derivatives, cuts and marginalization
# remain exactly as in class_cmb_survey_likelihood_forecast.py.
forecast.noise_cl = noise_cl_dimensionless


if __name__ == "__main__":
    forecast.main()
