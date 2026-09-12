#!/usr/bin/env python3
"""Final survey-aware forecast for the fixed-omega0 + fixed-early-Neff winner.

This wrapper applies two required corrections to the original survey forecast:
1. detector map depths in microkelvin-arcmin are converted to dimensionless
   DeltaT/T noise before combination with CLASS spectra;
2. the CLASS inputs use the final control's N_ur from control_parameters.json
   instead of the historical reference value N_ur=2.0328.

All covariance, derivative, multipole-cut and Fisher-projection logic is reused
unchanged from class_cmb_survey_likelihood_forecast.py.
"""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import numpy as np

import class_cmb_survey_likelihood_forecast as forecast
import class_cmb_survey_likelihood_units_fixed as units_fixed  # patches forecast.noise_cl
from class_response_optimize import H0, OMEGA_B, A_S, T_NCDM

N_S = forecast.N_S
TAU_REIO = forecast.TAU_REIO


def write_ini_final(path: Path, psd: Path, root: Path, params: dict,
                    overrides: dict, lmax: int) -> None:
    vals = {
        "H0": H0,
        "omega_b": OMEGA_B,
        "omega_cdm": float(params["omega_cdm"]),
        "A_s": A_S,
        "n_s": N_S,
        "tau_reio": TAU_REIO,
    }
    vals.update(overrides)
    n_ur = float(params["N_ur"])
    text = f"""# Final fixed-omega0 + fixed-early-Neff survey forecast input.
output = tCl,pCl,lCl,mPk
modes = s
lensing = yes
H0 = {vals['H0']:.15g}
omega_b = {vals['omega_b']:.15g}
omega_cdm = {vals['omega_cdm']:.15g}
A_s = {vals['A_s']:.15e}
n_s = {vals['n_s']:.15g}
tau_reio = {vals['tau_reio']:.15g}
N_ur = {n_ur:.15g}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd.resolve()}
m_ncdm = {float(params['mass_eV']):.15g}
T_ncdm = {T_NCDM}
deg_ncdm = {float(params['deg_ncdm']):.15g}
l_max_scalars = {lmax}
P_k_max_h/Mpc = 0.30
z_pk = 0
root = {root.resolve()}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


_original_prepare = forecast.prepare


def prepare_final(source: Path, work: Path, lmax: int) -> None:
    params_path = forecast.unique_find(source, "control_parameters.json")
    params = json.loads(params_path.read_text())
    if params.get("control") != "fixed-relic-omega0-fixed-early-neff":
        raise RuntimeError(f"Unexpected source control: {params.get('control')}")
    if abs(float(params["mass_eV"]) - 0.60) > 1e-12:
        raise RuntimeError(f"Expected final m=0.60 eV winner, got {params['mass_eV']}")
    if abs(float(params["controlled_early_Neff"]) - 3.046) > 1e-10:
        raise RuntimeError("Final source does not preserve early Neff=3.046")
    if "N_ur" not in params:
        raise RuntimeError("Final source control parameters do not contain N_ur")

    _original_prepare(source, work, lmax)
    manifest_path = work / "forecast_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["final_control_validation"] = {
        "control": params["control"],
        "mass_eV": float(params["mass_eV"]),
        "N_ur": float(params["N_ur"]),
        "deg_ncdm": float(params["deg_ncdm"]),
        "controlled_omega_ncdm_z0": float(params["controlled_omega_ncdm_z0"]),
        "controlled_early_Neff": float(params["controlled_early_Neff"]),
        "map_noise_units": "microK-arcmin converted to dimensionless DeltaT/T before D_ell conversion",
    }
    manifest["forecast_scope"] += (
        "; final hidden-state background keeps present-day relic omega_ncdm and "
        "early total neutrino Neff fixed"
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["final_control_validation"], indent=2))


forecast.write_ini = write_ini_final
forecast.prepare = prepare_final


if __name__ == "__main__":
    forecast.main()
