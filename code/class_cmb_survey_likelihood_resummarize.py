#!/usr/bin/env python3
"""Corrected survey-aware CMB likelihood summary using existing CLASS spectra.

The original survey summary mixed CLASS dimensionless D_ell spectra with detector
noise expressed in thermodynamic microkelvin units.  This wrapper converts map
depths in uK-arcmin to dimensionless DeltaT/Tcmb noise before reusing the validated
likelihood/Fisher machinery.  No CLASS spectra are recomputed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import class_cmb_survey_likelihood_forecast as base
from class_response_optimize import TCMB_K


def dimensionless_noise_cl(ell: np.ndarray, exp: dict, pol: bool) -> np.ndarray:
    inv = np.zeros_like(ell, dtype=float)
    for ch in exp["channels"]:
        depth_uK_arcmin = float(
            ch["noise_P_uK_arcmin"] if pol else ch["noise_T_uK_arcmin"]
        )
        fwhm = float(ch["beam_arcmin"])
        # CLASS default CMB spectra are dimensionless D_ell = l(l+1)C_l/(2pi).
        # Convert thermodynamic map depth from uK-arcmin to dimensionless
        # DeltaT/Tcmb-radian before forming N_l.
        delta = (
            depth_uK_arcmin
            * 1.0e-6
            / TCMB_K
            * math.pi
            / (180.0 * 60.0)
        )
        sigma_b = (
            fwhm
            * math.pi
            / (180.0 * 60.0)
            / math.sqrt(8.0 * math.log(2.0))
        )
        ncl = delta * delta * np.exp(ell * (ell + 1.0) * sigma_b * sigma_b)
        inv += np.where(ncl > 0, 1.0 / ncl, 0.0)
    return np.where(inv > 0, 1.0 / inv, np.inf)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    base.noise_cl = dimensionless_noise_cl
    base.summarize(args.results, args.manifest, args.output)

    out = json.loads(args.output.read_text())
    out["unit_correction"] = {
        "status": "corrected",
        "class_spectrum_units": "dimensionless D_ell = l(l+1) C_l /(2 pi)",
        "detector_noise_input_units": "thermodynamic uK-arcmin",
        "conversion": "depth * 1e-6 / Tcmb * pi/(180*60)",
        "source_spectra_reused": True,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
