#!/usr/bin/env python3
"""Exact bounded Planck-2018 Plik-lite fit for the frozen hidden-state family.

This is a robustness control for the earlier local linearized real-data fit.
For alpha in {-1,0,+1}, the hidden kinetic distribution is frozen and CLASS is
recomputed at every nuisance/cosmological parameter evaluation.  The fit uses
the public Planck 2018 Plik-lite high-l TTTEEE data and covariance, a Gaussian
tau prior, and a Gaussian Planck calibration prior.

The purpose is not to claim a realistic 0.60-eV neutrino model.  The tested
state is the controlled fixed-present-day-relic-abundance benchmark used in the
manuscript validation.  The result therefore quantifies robustness of the
Planck likelihood comparison within that benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import scipy.linalg
from scipy.optimize import minimize

TCMB_K = 2.7255
UK2_PER_DIMENSIONLESS = (TCMB_K * 1.0e6) ** 2
TAU0 = 0.054
TAU_SIGMA = 0.007
CAL_SIGMA = 0.0025
N_UR = 2.0328
T_NCDM = 0.71611

# Broad but finite bounds keep the exact nonlinear fit inside a region where the
# controlled cosmological benchmark remains interpretable and CLASS is stable.
BOUNDS = {
    "H0": (55.0, 80.0),
    "omega_b": (0.0200, 0.0250),
    "omega_cdm": (0.095, 0.145),
    "ln10As": (2.80, 3.30),
    "n_s": (0.90, 1.03),
    "tau_reio": (0.020, 0.100),
    "A_planck": (0.985, 1.015),
}
NAMES = list(BOUNDS)
DEFAULT_X0 = np.array([67.36, 0.02237, 0.1200, 3.044, 0.9649, TAU0, 1.0], dtype=float)


def unique_find(root: Path, name: str) -> Path:
    hits = list(root.rglob(name))
    if len(hits) != 1:
        raise RuntimeError(f"Expected exactly one {name} below {root}, found {hits}")
    return hits[0]


def load_planck_lite(root: Path):
    module_path = root / "planck_lite_py.py"
    spec = importlib.util.spec_from_file_location("planck_lite_py_external", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {module_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.PlanckLitePy(
        data_directory=str(root / "data"), year=2018,
        spectra="TTTEEE", use_low_ell_bins=False,
    )


def read_class_lensed(path: Path) -> dict[str, np.ndarray]:
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
    import re
    labels = {}
    for line in header:
        for num, name in re.findall(r"(\d+):([^\s]+)", line):
            labels[name.strip().lower()] = int(num) - 1
    def col(name: str) -> int:
        if name not in labels:
            raise RuntimeError(f"Missing {name} in {path}; labels={labels}")
        return labels[name]
    return {
        "ell": arr[:, 0].astype(int),
        "TT": arr[:, col("tt")] * UK2_PER_DIMENSIONLESS,
        "TE": arr[:, col("te")] * UK2_PER_DIMENSIONLESS,
        "EE": arr[:, col("ee")] * UK2_PER_DIMENSIONLESS,
    }


def valid_bin_indices(like, ell_max: int):
    keep = []
    counts = {"TT": 0, "TE": 0, "EE": 0}
    for i in range(like.nbintt):
        if int(like.blmax_TT[i] + like.plmin_TT) <= ell_max:
            keep.append(i); counts["TT"] += 1
    te0 = like.nbintt
    for i in range(like.nbinte):
        if int(like.blmax[i] + like.plmin) <= ell_max:
            keep.append(te0 + i); counts["TE"] += 1
    ee0 = like.nbintt + like.nbinte
    for i in range(like.nbinee):
        if int(like.blmax[i] + like.plmin) <= ell_max:
            keep.append(ee0 + i); counts["EE"] += 1
    return keep, counts


def binned_model(like, sp: dict[str, np.ndarray], keep: list[int]) -> np.ndarray:
    ell = sp["ell"]
    ellmin, ellmax = int(ell[0]), int(ell[-1])
    ls = ell.astype(float)
    fac = ls * (ls + 1.0) / (2.0 * math.pi)
    fac = np.where(fac > 0, fac, np.inf)
    Cltt, Clte, Clee = sp["TT"] / fac, sp["TE"] / fac, sp["EE"] / fac
    full = np.full(like.nbin_tot, np.nan, dtype=float)
    for i in range(like.nbintt):
        lo = int(like.blmin_TT[i] + like.plmin_TT)
        hi = int(like.blmax_TT[i] + like.plmin_TT)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w_TT[like.blmin_TT[i]:like.blmax_TT[i] + 1]
        full[i] = float(np.sum(Cltt[sl] * w))
    off = like.nbintt
    for i in range(like.nbinte):
        lo = int(like.blmin[i] + like.plmin)
        hi = int(like.blmax[i] + like.plmin)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w[like.blmin[i]:like.blmax[i] + 1]
        full[off + i] = float(np.sum(Clte[sl] * w))
    off = like.nbintt + like.nbinte
    for i in range(like.nbinee):
        lo = int(like.blmin[i] + like.plmin)
        hi = int(like.blmax[i] + like.plmin)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w[like.blmin[i]:like.blmax[i] + 1]
        full[off + i] = float(np.sum(Clee[sl] * w))
    vec = full[np.asarray(keep, dtype=int)]
    if not np.all(np.isfinite(vec)):
        raise RuntimeError("Non-finite binned model")
    return vec


def make_psd(source: Path, alpha: float, out: Path) -> dict:
    plus = np.loadtxt(unique_find(source, "winner_plus.dat"))
    minus = np.loadtxt(unique_find(source, "winner_minus.dat"))
    if plus.shape != minus.shape or not np.allclose(plus[:, 0], minus[:, 0], rtol=0, atol=1e-13):
        raise RuntimeError("winner PSD grids differ")
    mid = 0.5 * (plus[:, 1] + minus[:, 1])
    delta = 0.5 * (plus[:, 1] - minus[:, 1])
    f = mid + alpha * delta
    if np.any(f <= 0):
        raise RuntimeError("Requested alpha leaves positive PSD family")
    np.savetxt(out, np.column_stack([plus[:, 0], f]), fmt="%.17e")
    return {"min_psd": float(np.min(f)), "max_relative_to_mid": float(np.max(np.abs(f-mid)/mid))}


def write_ini(path: Path, root: Path, psd: Path, params: dict, x: np.ndarray, lmax: int):
    H0, ob, oc, ln10As, ns, tau, _ = map(float, x)
    As = math.exp(ln10As) * 1.0e-10
    text = f"""output = tCl,pCl,lCl
modes = s
lensing = yes
H0 = {H0:.15g}
omega_b = {ob:.15g}
omega_cdm = {oc:.15g}
A_s = {As:.15e}
n_s = {ns:.15g}
tau_reio = {tau:.15g}
N_ur = {N_UR}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd.resolve()}
m_ncdm = {float(params['mass_eV']):.15g}
T_ncdm = {T_NCDM}
deg_ncdm = {float(params['deg_ncdm']):.15g}
l_max_scalars = {lmax}
root = {root.resolve()}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


def initial_from_json(path: Path | None) -> np.ndarray:
    if path is None:
        return DEFAULT_X0.copy()
    obj = json.loads(path.read_text())
    bp = obj.get("best_parameters", obj)
    vals = []
    for n, default in zip(NAMES, DEFAULT_X0):
        vals.append(float(bp.get(n, default)))
    x = np.array(vals, dtype=float)
    for i, n in enumerate(NAMES):
        lo, hi = BOUNDS[n]
        x[i] = np.clip(x[i], lo, hi)
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--class-exe", type=Path, required=True)
    ap.add_argument("--precision-pre", type=Path, required=True)
    ap.add_argument("--planck-lite-root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--initial-json", type=Path)
    ap.add_argument("--maxfev", type=int, default=140)
    ap.add_argument("--lmax", type=int, default=2500)
    args = ap.parse_args()
    if args.alpha not in (-1.0, 0.0, 1.0):
        raise SystemExit("This robustness control is restricted to alpha=-1,0,+1")
    args.outdir.mkdir(parents=True, exist_ok=True)
    runs = args.outdir / "class_runs"
    runs.mkdir(exist_ok=True)
    psd = args.outdir / f"alpha_{args.alpha:+.0f}.dat"
    psd_info = make_psd(args.source, args.alpha, psd)
    params = json.loads(unique_find(args.source, "control_parameters.json").read_text())

    like = load_planck_lite(args.planck_lite_root)
    keep, counts = valid_bin_indices(like, args.lmax)
    idx = np.asarray(keep, dtype=int)
    cov_full = scipy.linalg.inv(like.fisher)
    cov = cov_full[np.ix_(idx, idx)]
    precision = scipy.linalg.cho_solve(scipy.linalg.cho_factor(cov), np.eye(cov.shape[0]))
    data = np.asarray(like.X_data, dtype=float)[idx]

    cache = {}
    history = []

    def objective(x):
        x = np.asarray(x, dtype=float)
        key = tuple(round(float(v), 10) for v in x)
        if key in cache:
            return cache[key][0]
        token = hashlib.sha1(repr(key).encode()).hexdigest()[:12]
        root = runs / f"r_{token}_"
        ini = runs / f"r_{token}.ini"
        write_ini(ini, root, psd, params, x, args.lmax)
        log = runs / f"r_{token}.log"
        with log.open("w") as fh:
            cp = subprocess.run([str(args.class_exe.resolve()), str(ini.resolve()), str(args.precision_pre.resolve())], stdout=fh, stderr=subprocess.STDOUT)
        cl = runs / f"r_{token}_00_cl_lensed.dat"
        if cp.returncode != 0 or not cl.exists():
            val = 1.0e30
            cache[key] = (val, None)
            history.append({"x": x.tolist(), "chi2": val, "status": "CLASS_FAIL"})
            return val
        sp = read_class_lensed(cl)
        model = binned_model(like, sp, keep)
        cal = float(x[6])
        model = model / (cal * cal)
        r = data - model
        chi2_data = float(r @ precision @ r)
        chi2_prior = ((float(x[5]) - TAU0) / TAU_SIGMA) ** 2 + ((cal - 1.0) / CAL_SIGMA) ** 2
        val = chi2_data + chi2_prior
        cache[key] = (val, {"chi2_data": chi2_data, "chi2_prior": chi2_prior})
        history.append({"x": x.tolist(), "chi2": val, "chi2_data": chi2_data, "chi2_prior": chi2_prior, "status": "ok"})
        print(json.dumps({"alpha": args.alpha, "eval": len(history), "chi2": val, "x": x.tolist()}), flush=True)
        return val

    x0 = initial_from_json(args.initial_json)
    scipy_bounds = [BOUNDS[n] for n in NAMES]
    f0 = objective(x0)
    res = minimize(objective, x0, method="Powell", bounds=scipy_bounds,
                   options={"maxfev": args.maxfev, "xtol": 2e-4, "ftol": 2e-5, "disp": True})
    xbest = np.asarray(res.x, dtype=float)
    fbest = float(objective(xbest))
    at_bounds = {}
    for i, n in enumerate(NAMES):
        lo, hi = BOUNDS[n]
        scale = max(abs(hi-lo), 1.0)
        at_bounds[n] = bool(abs(xbest[i]-lo) < 1e-4*scale or abs(xbest[i]-hi) < 1e-4*scale)
    out = {
        "experiment": "Planck 2018 exact bounded nonlinear hidden-mode profile",
        "alpha": args.alpha,
        "chi2_start": f0,
        "chi2_best": fbest,
        "best_parameters": {n: float(v) for n, v in zip(NAMES, xbest)},
        "bounds": {n: list(map(float, BOUNDS[n])) for n in NAMES},
        "at_bounds": at_bounds,
        "optimizer": {"method": "Powell", "success": bool(res.success), "message": str(res.message), "nfev": int(res.nfev), "nit": int(getattr(res, "nit", -1)), "maxfev": args.maxfev},
        "data_selection": {"retained_bins": counts, "n_bins": len(keep), "ell_max": args.lmax},
        "hidden_family": {"mass_eV": float(params["mass_eV"]), "deg_ncdm": float(params["deg_ncdm"]), "control": "fixed present-day relic omega_ncdm", "alpha": args.alpha, **psd_info},
        "priors": {"tau": [TAU0, TAU_SIGMA], "A_planck": [1.0, CAL_SIGMA]},
        "scope": "Exact CLASS recomputation at each bounded optimization point. Controlled massive-relic benchmark; not a realistic standard-neutrino fit.",
    }
    (args.outdir / "result.json").write_text(json.dumps(out, indent=2) + "\n")
    (args.outdir / "history.json").write_text(json.dumps(history, indent=2) + "\n")
    print(json.dumps(out, indent=2))

if __name__ == "__main__":
    main()
