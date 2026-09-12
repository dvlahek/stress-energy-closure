#!/usr/bin/env python3
"""Survey-aware Gaussian CMB likelihood forecast for the hidden-state winner.

This is a validation forecast, not an official Planck/SO/CMB-S4 likelihood.
It upgrades the ideal full-sky cosmic-variance S/N by including:

* finite sky fraction;
* beam-deconvolved white detector noise;
* experiment-specific multipole cuts;
* the joint TT/TE/EE Gaussian band-power covariance;
* linearized marginalization over the six standard background/primordial
  parameters H0, omega_b, omega_cdm, ln(A_s), n_s and tau_reio.

The nuisance derivatives are generated with CLASS around the midpoint hidden
state.  A Planck-like low-E prior sigma(tau)=0.007 is included in all three
survey configurations.  Foreground nuisance templates, non-Gaussian lensing
covariance, reconstruction noise, calibration systematics and the exact Planck
Plik nuisance model are deliberately not included.  The output must therefore
be described as a survey-aware Fisher/Gaussian-likelihood forecast rather than
as a reproduction of an official experimental likelihood.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np

from class_response_optimize import (
    H0, OMEGA_B, A_S, N_UR, T_NCDM, write_psd,
)

N_S = 0.9649
TAU_REIO = 0.054

# Central finite-difference steps.  They are deliberately small enough for a
# local Fisher projection but large compared with CLASS numerical noise.
DERIVATIVE_STEPS = {
    "H0": 0.50,
    "omega_b": 1.5e-4,
    "omega_cdm": 8.0e-4,
    "lnAs": 1.0e-2,
    "n_s": 5.0e-3,
    "tau_reio": 4.0e-3,
}

# Reference survey configurations.  Planck multipole ranges follow the 2018
# Plik high-l likelihood (TT to 2508, TE/EE to 1996).  The effective f_sky=0.50
# is deliberately conservative relative to the frequency-dependent Plik masks.
# The Planck 143/217 GHz map sensitivities correspond to approximately
# 33/47 uK-arcmin in T and 70/105 uK-arcmin in polarization, with 7.3/5.0 arcmin
# beams.  SO uses the published LAT combined 93+145 GHz baseline depth of
# 6 uK-arcmin in T and sqrt(2) higher polarization depth over f_sky~0.4; an
# effective 1.4 arcmin beam is conservative for the combination.  CMB-S4 is
# retained as the commonly used reference-design benchmark: 1 uK-arcmin T,
# 1.4 uK-arcmin P, 1.4 arcmin beam, f_sky=0.4.
EXPERIMENTS = {
    "Planck2018_like": {
        "fsky": 0.50,
        "lmin": 30,
        "lmax_TT": 2500,
        "lmax_P": 1996,
        "channels": [
            {"name": "143", "beam_arcmin": 7.3, "noise_T_uK_arcmin": 33.0,
             "noise_P_uK_arcmin": 70.0},
            {"name": "217", "beam_arcmin": 5.0, "noise_T_uK_arcmin": 47.0,
             "noise_P_uK_arcmin": 105.0},
        ],
        "tau_prior_sigma": 0.007,
        "description": "Planck-2018-like high-l TT/TE/EE Gaussian forecast",
    },
    "SimonsObservatory_LAT_baseline": {
        "fsky": 0.40,
        "lmin": 30,
        "lmax_TT": 2500,
        "lmax_P": 2500,
        "channels": [
            {"name": "93+145 effective", "beam_arcmin": 1.4,
             "noise_T_uK_arcmin": 6.0,
             "noise_P_uK_arcmin": 6.0 * math.sqrt(2.0)},
        ],
        "tau_prior_sigma": 0.007,
        "description": "SO-LAT baseline-depth reference forecast",
    },
    "CMB_S4_reference_design": {
        "fsky": 0.40,
        "lmin": 30,
        "lmax_TT": 2500,
        "lmax_P": 2500,
        "channels": [
            {"name": "effective", "beam_arcmin": 1.4,
             "noise_T_uK_arcmin": 1.0,
             "noise_P_uK_arcmin": 1.4},
        ],
        "tau_prior_sigma": 0.007,
        "description": "CMB-S4 archived reference-design benchmark",
    },
}


def unique_find(root: Path, name: str) -> Path:
    hits = list(root.rglob(name))
    if len(hits) != 1:
        raise RuntimeError(f"Expected exactly one {name} below {root}, found {hits}")
    return hits[0]


def read_table(path: Path):
    header = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                break
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr[None, :]
    labels = {}
    import re
    for line in header:
        for num, name in re.findall(r"(\d+):([^\s]+)", line):
            labels[name.strip().lower()] = int(num) - 1
    return arr, labels


def col(labels: dict, *names: str) -> int:
    for n in names:
        if n.lower() in labels:
            return labels[n.lower()]
    raise KeyError(f"Missing columns {names}; labels={labels}")


def load_ttee(path: Path):
    a, labels = read_table(path)
    return {
        "ell": a[:, 0].astype(int),
        "TT": a[:, col(labels, "tt")],
        "EE": a[:, col(labels, "ee")],
        "TE": a[:, col(labels, "te")],
    }


def write_ini(path: Path, psd: Path, root: Path, params: dict,
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
    text = f"""# Survey-aware hidden-state likelihood forecast input.
output = tCl,pCl,lCl,mPk
modes = s
lensing = yes
H0 = {vals['H0']:.15g}
omega_b = {vals['omega_b']:.15g}
omega_cdm = {vals['omega_cdm']:.15g}
A_s = {vals['A_s']:.15e}
n_s = {vals['n_s']:.15g}
tau_reio = {vals['tau_reio']:.15g}
N_ur = {N_UR}
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


def prepare(source: Path, work: Path, lmax: int) -> None:
    work.mkdir(parents=True, exist_ok=True)
    ini_dir = work / "ini"
    dat_dir = work / "data"
    out_dir = work / "outputs"
    for d in (ini_dir, dat_dir, out_dir):
        d.mkdir(parents=True, exist_ok=True)

    plus_src = unique_find(source, "winner_plus.dat")
    minus_src = unique_find(source, "winner_minus.dat")
    params_path = unique_find(source, "control_parameters.json")
    params = json.loads(params_path.read_text())

    plus = np.loadtxt(plus_src)
    minus = np.loadtxt(minus_src)
    if plus.shape != minus.shape or not np.allclose(plus[:, 0], minus[:, 0], rtol=0, atol=1e-13):
        raise RuntimeError("Winner PSD grids do not agree")
    if plus.shape[1] < 2:
        raise RuntimeError("Winner PSD files require at least two columns")

    pplus = dat_dir / "winner_plus.dat"
    pminus = dat_dir / "winner_minus.dat"
    pmid = dat_dir / "winner_midpoint.dat"
    shutil.copy2(plus_src, pplus)
    shutil.copy2(minus_src, pminus)
    q = plus[:, 0]
    fmid = 0.5 * (plus[:, 1] + minus[:, 1])
    if np.any(fmid <= 0):
        raise RuntimeError("Midpoint PSD is not positive")
    write_psd(pmid, q, fmid)

    jobs = {}
    for name, psd in [("winner_plus", pplus), ("winner_minus", pminus), ("fiducial", pmid)]:
        ini = ini_dir / f"{name}.ini"
        write_ini(ini, psd, out_dir / f"{name}_", params, {}, lmax)
        jobs[name] = {"ini": str(ini), "kind": name}

    for par, step in DERIVATIVE_STEPS.items():
        for sign, fac in (("plus", +1.0), ("minus", -1.0)):
            ov = {}
            if par == "lnAs":
                ov["A_s"] = A_S * math.exp(fac * step)
            elif par == "H0":
                ov["H0"] = H0 + fac * step
            elif par == "omega_b":
                ov["omega_b"] = OMEGA_B + fac * step
            elif par == "omega_cdm":
                ov["omega_cdm"] = float(params["omega_cdm"]) + fac * step
            elif par == "n_s":
                ov["n_s"] = N_S + fac * step
            elif par == "tau_reio":
                ov["tau_reio"] = TAU_REIO + fac * step
            else:
                raise KeyError(par)
            name = f"d_{par}_{sign}"
            ini = ini_dir / f"{name}.ini"
            write_ini(ini, pmid, out_dir / f"{name}_", params, ov, lmax)
            jobs[name] = {"ini": str(ini), "kind": "derivative", "parameter": par,
                          "sign": sign, "step": step}

    manifest = {
        "control_parameters": params,
        "lmax": lmax,
        "derivative_steps": DERIVATIVE_STEPS,
        "experiments": EXPERIMENTS,
        "jobs": jobs,
        "forecast_scope": (
            "Gaussian TT/TE/EE band-power likelihood with survey f_sky, white noise, "
            "Gaussian beams and linearized six-parameter LambdaCDM marginalization"
        ),
        "limitations": [
            "not the official Planck Plik likelihood",
            "no foreground nuisance templates or frequency cross-spectrum nuisance model",
            "no lensing-reconstruction likelihood",
            "no non-Gaussian lensing covariance",
            "no calibration or beam-uncertainty nuisance parameters",
            "CMB-S4 is used only as an archived reference-design benchmark",
        ],
    }
    (work / "forecast_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def noise_cl(ell: np.ndarray, exp: dict, pol: bool) -> np.ndarray:
    inv = np.zeros_like(ell, dtype=float)
    for ch in exp["channels"]:
        depth = float(ch["noise_P_uK_arcmin"] if pol else ch["noise_T_uK_arcmin"])
        fwhm = float(ch["beam_arcmin"])
        delta = depth * math.pi / (180.0 * 60.0)
        sigma_b = fwhm * math.pi / (180.0 * 60.0) / math.sqrt(8.0 * math.log(2.0))
        ncl = delta * delta * np.exp(ell * (ell + 1.0) * sigma_b * sigma_b)
        inv += np.where(ncl > 0, 1.0 / ncl, 0.0)
    return np.where(inv > 0, 1.0 / inv, np.inf)


def cl_to_dl_noise(ell: np.ndarray, ncl: np.ndarray) -> np.ndarray:
    return ell * (ell + 1.0) * ncl / (2.0 * math.pi)


def align_and_load(results: Path, name: str):
    p = unique_find(results, f"{name}_00_cl_lensed.dat")
    return load_ttee(p)


def fisher_forecast(results: Path, manifest: dict, exp: dict) -> dict:
    fid = align_and_load(results, "fiducial")
    wp = align_and_load(results, "winner_plus")
    wm = align_and_load(results, "winner_minus")
    ell = fid["ell"]
    for x in (wp, wm):
        if not np.array_equal(ell, x["ell"]):
            raise RuntimeError("ell grids do not match")

    deriv = {}
    for par, step in DERIVATIVE_STEPS.items():
        ap = align_and_load(results, f"d_{par}_plus")
        am = align_and_load(results, f"d_{par}_minus")
        if not np.array_equal(ell, ap["ell"]) or not np.array_equal(ell, am["ell"]):
            raise RuntimeError(f"ell grid mismatch for {par}")
        deriv[par] = {
            k: (ap[k] - am[k]) / (2.0 * step) for k in ("TT", "EE", "TE")
        }

    dsignal = {k: wp[k] - wm[k] for k in ("TT", "EE", "TE")}
    nT = cl_to_dl_noise(ell, noise_cl(ell.astype(float), exp, pol=False))
    nP = cl_to_dl_noise(ell, noise_cl(ell.astype(float), exp, pol=True))
    pars = list(DERIVATIVE_STEPS)
    np_ = len(pars)
    F = np.zeros((np_, np_))
    b = np.zeros(np_)
    sn2 = 0.0
    sn2_tt = 0.0
    sn2_ee = 0.0
    used_joint = 0

    fsky = float(exp["fsky"])
    lmin = int(exp["lmin"])
    lmaxT = int(exp["lmax_TT"])
    lmaxP = int(exp["lmax_P"])

    for i, l in enumerate(ell):
        if l < lmin or l > max(lmaxT, lmaxP):
            continue
        fac = (2.0 * l + 1.0) * fsky
        T = float(fid["TT"][i])
        E = float(fid["EE"][i])
        X = float(fid["TE"][i])
        Tt = T + float(nT[i])
        Et = E + float(nP[i])

        if l <= lmaxT and Tt > 0:
            varT = 2.0 * Tt * Tt / fac
            sn2_tt += float(dsignal["TT"][i]) ** 2 / varT
        if l <= lmaxP and Et > 0:
            varE = 2.0 * Et * Et / fac
            sn2_ee += float(dsignal["EE"][i]) ** 2 / varE

        if l <= lmaxP and l <= lmaxT and Tt > 0 and Et > 0:
            cov = np.array([
                [2*Tt*Tt, 2*X*X, 2*Tt*X],
                [2*X*X, 2*Et*Et, 2*Et*X],
                [2*Tt*X, 2*Et*X, X*X + Tt*Et],
            ]) / fac
            s = np.array([dsignal["TT"][i], dsignal["EE"][i], dsignal["TE"][i]])
            D = np.array([[deriv[p]["TT"][i], deriv[p]["EE"][i], deriv[p]["TE"][i]]
                          for p in pars]).T
        elif l <= lmaxT and Tt > 0:
            cov = np.array([[2*Tt*Tt / fac]])
            s = np.array([dsignal["TT"][i]])
            D = np.array([[deriv[p]["TT"][i] for p in pars]])
        else:
            continue

        inv = np.linalg.pinv(cov, rcond=1e-12)
        term = float(s @ inv @ s)
        if not np.isfinite(term) or term < -1e-10:
            continue
        sn2 += max(term, 0.0)
        b += D.T @ inv @ s
        F += D.T @ inv @ D
        used_joint += 1

    # A common low-E prior is used to prevent the high-l As-tau degeneracy from
    # becoming an artificial singular direction in the linearized projection.
    tau_sigma = float(exp.get("tau_prior_sigma", 0.0) or 0.0)
    if tau_sigma > 0:
        j = pars.index("tau_reio")
        F[j, j] += 1.0 / (tau_sigma * tau_sigma)

    Finv = np.linalg.pinv(F, rcond=1e-10)
    absorbed = float(b @ Finv @ b)
    marg2 = max(sn2 - absorbed, 0.0)
    shift = Finv @ b
    evals = np.linalg.eigvalsh(0.5 * (F + F.T))

    return {
        "fixed_cosmology": {
            "joint_TT_TE_EE_SN": math.sqrt(max(sn2, 0.0)),
            "TT_only_SN": math.sqrt(max(sn2_tt, 0.0)),
            "EE_only_SN": math.sqrt(max(sn2_ee, 0.0)),
            "delta_chi2": sn2,
        },
        "linearized_LCDM_marginalized": {
            "joint_TT_TE_EE_SN": math.sqrt(marg2),
            "delta_chi2": marg2,
            "fraction_delta_chi2_retained": (marg2 / sn2 if sn2 > 0 else None),
            "fraction_signal_norm_absorbed_by_LCDM": (absorbed / sn2 if sn2 > 0 else None),
            "best_fit_parameter_shifts": {p: float(v) for p, v in zip(pars, shift)},
        },
        "fisher_diagnostics": {
            "parameters": pars,
            "eigenvalues": [float(x) for x in evals],
            "multipoles_used": used_joint,
        },
        "experiment": exp,
    }


def summarize(results: Path, manifest_path: Path, out: Path) -> None:
    manifest = json.loads(manifest_path.read_text())
    forecasts = {}
    for name, exp in EXPERIMENTS.items():
        forecasts[name] = fisher_forecast(results, manifest, exp)
    final = {
        "forecast_type": "survey-aware Gaussian TT/TE/EE likelihood + linearized LCDM Fisher marginalization",
        "source_control_parameters": manifest["control_parameters"],
        "forecasts": forecasts,
        "limitations": manifest["limitations"],
        "interpretation": (
            "Fixed-cosmology S/N includes finite sky, beam and white detector noise. "
            "Marginalized S/N additionally projects the hidden-state spectral difference "
            "against local variations of the six standard LambdaCDM parameters."
        ),
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(final, indent=2) + "\n")
    print(json.dumps(final, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["prepare", "summarize"], required=True)
    ap.add_argument("--source", type=Path)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path)
    ap.add_argument("--manifest", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--lmax", type=int, default=2500)
    args = ap.parse_args()

    if args.mode == "prepare":
        if args.source is None:
            ap.error("--source is required for prepare")
        prepare(args.source, args.work, args.lmax)
    else:
        if args.results is None:
            ap.error("--results is required for summarize")
        manifest = args.manifest or (args.work / "forecast_manifest.json")
        output = args.output or (args.work / "cmb_survey_likelihood_forecast.json")
        summarize(args.results, manifest, output)


if __name__ == "__main__":
    main()
