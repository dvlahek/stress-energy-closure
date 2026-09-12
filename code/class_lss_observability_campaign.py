#!/usr/bin/env python3
"""Idealized 3D large-scale-structure observability campaign for the final hidden mode.

This campaign asks a deliberately narrow question before any survey-specific
galaxy/weak-lensing forecast: can three-dimensional matter clustering recover
kinetic information that is hidden from the matched instantaneous moments?

The physical class is frozen to the validated final control:
  * m_ncdm = 0.60 eV;
  * matched (n, rho, P) at z_match = 1100;
  * fixed omega_ncdm(z=0);
  * fixed early total N_eff = 3.046;
  * smooth positive distributions with a +/-30% pointwise cap.

The code performs three diagnostics:
  A. evaluate the already validated frozen CMB winner in P_m(k,z);
  B. optimize the null-space direction directly for *LCDM-marginalized*
     tomographic 3D matter-power distinguishability;
  C. validate the optimized pair with fresh CLASS runs.

The covariance is the ideal Gaussian matter-field cosmic-variance covariance.
It is not a realistic galaxy-survey likelihood.  Tomographic bins are assigned
an explicitly stated equal share of a total effective volume; outputs are
reported for several total volumes and k_max values.  This makes the result a
clean observability screen and a decision gate for a later realistic
galaxy/weak-lensing forecast.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
from pathlib import Path

import numpy as np

from class_response_optimize import (
    H0, OMEGA_B, A_S, T_NCDM,
    kinetic_objects, pair_stats, write_psd,
)
import class_scalar_fixed_omega0_neff_run as fixed
from class_cmb_survey_likelihood_forecast import (
    N_S, TAU_REIO, DERIVATIVE_STEPS,
)

REDSHIFTS = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)
KMAX_VALUES = (0.10, 0.20, 0.30)
TOTAL_VOLUMES_GPCH3 = (10.0, 30.0, 50.0, 100.0)
PROBE_FRAC = 0.08
FINAL_FRAC = 0.30
RANDOM_SAMPLES = 100000
RANDOM_SEED = 314159


def ztag(z: float) -> str:
    return f"z{int(round(100.0*z)):03d}"


def unique_find(root: Path, name: str) -> Path:
    hits = list(root.rglob(name))
    if len(hits) != 1:
        raise RuntimeError(f"Expected exactly one {name} below {root}, found {hits}")
    return hits[0]


def read_psd(path: Path):
    a = np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 2:
        raise RuntimeError(f"Unexpected PSD format: {path}")
    return a[:, 0], a[:, 1]


def write_mpk_ini(path: Path, psd: Path, root: Path, params: dict,
                  z: float, pkmax: float, overrides: dict | None = None) -> None:
    overrides = {} if overrides is None else dict(overrides)
    vals = {
        "H0": H0,
        "omega_b": OMEGA_B,
        "omega_cdm": float(params["omega_cdm"]),
        "A_s": A_S,
        "n_s": N_S,
        "tau_reio": TAU_REIO,
    }
    vals.update(overrides)
    text = f"""# Idealized LSS hidden-state observability input.
output = mPk
modes = s
H0 = {vals['H0']:.15g}
omega_b = {vals['omega_b']:.15g}
omega_cdm = {vals['omega_cdm']:.15g}
A_s = {vals['A_s']:.15e}
n_s = {vals['n_s']:.15g}
tau_reio = {vals['tau_reio']:.15g}
N_ur = {float(params['N_ur']):.15g}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd.resolve()}
m_ncdm = {float(params['mass_eV']):.15g}
T_ncdm = {T_NCDM}
deg_ncdm = {float(params['deg_ncdm']):.15g}
P_k_max_h/Mpc = {pkmax:.15g}
z_pk = {z:.15g}
root = {root.resolve()}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


def parameter_override(par: str, sign: float, params: dict) -> dict:
    step = float(DERIVATIVE_STEPS[par])
    if par == "lnAs":
        return {"A_s": A_S * math.exp(sign * step)}
    if par == "H0":
        return {"H0": H0 + sign * step}
    if par == "omega_b":
        return {"omega_b": OMEGA_B + sign * step}
    if par == "omega_cdm":
        return {"omega_cdm": float(params["omega_cdm"]) + sign * step}
    if par == "n_s":
        return {"n_s": N_S + sign * step}
    if par == "tau_reio":
        return {"tau_reio": TAU_REIO + sign * step}
    raise KeyError(par)


def write_case_inis(work: Path, case: str, psd: Path, params: dict,
                    pkmax: float, overrides: dict | None = None) -> None:
    for z in REDSHIFTS:
        tag = ztag(z)
        write_mpk_ini(
            work / "ini" / f"{case}_{tag}.ini",
            psd,
            work / "outputs" / f"{case}_{tag}_",
            params,
            z,
            pkmax,
            overrides,
        )


def prepare(source: Path, work: Path, pkmax: float) -> None:
    for d in (work, work / "data", work / "ini", work / "outputs"):
        d.mkdir(parents=True, exist_ok=True)

    plus_src = unique_find(source, "winner_plus.dat")
    minus_src = unique_find(source, "winner_minus.dat")
    params = json.loads(unique_find(source, "control_parameters.json").read_text())

    if params.get("control") != fixed.CONTROL:
        raise RuntimeError(f"Unexpected source control: {params.get('control')}")
    if abs(float(params["mass_eV"]) - 0.60) > 1e-12:
        raise RuntimeError(f"Expected m=0.60 eV source, got {params['mass_eV']}")
    if abs(float(params["controlled_early_Neff"]) - 3.046) > 1e-10:
        raise RuntimeError("Source does not preserve early N_eff=3.046")

    qp, fp = read_psd(plus_src)
    qm, fm = read_psd(minus_src)
    if qp.shape != qm.shape or not np.allclose(qp, qm, rtol=0, atol=1e-13):
        raise RuntimeError("Winner PSD q grids differ")
    q = qp
    fmid = 0.5 * (fp + fm)
    if np.any(fmid <= 0):
        raise RuntimeError("Non-positive winner midpoint")

    ffd, weights, basis, N, shapes, y, M = kinetic_objects(
        q, float(params["mass_eV"]), float(params["z_match"])
    )
    rel_mid = float(np.max(np.abs(fmid - ffd) / np.maximum(ffd, 1e-300)))
    if rel_mid > 1e-7:
        raise RuntimeError(
            f"Source midpoint is not the expected FD baseline: max relative mismatch={rel_mid}"
        )

    np.save(work / "q.npy", q)
    np.save(work / "fmid.npy", fmid)
    np.save(work / "weights.npy", weights)
    np.save(work / "null_shapes.npy", shapes)

    frozen_plus = work / "data" / "frozen_plus.dat"
    frozen_minus = work / "data" / "frozen_minus.dat"
    midpoint = work / "data" / "fiducial.dat"
    shutil.copy2(plus_src, frozen_plus)
    shutil.copy2(minus_src, frozen_minus)
    write_psd(midpoint, q, fmid)
    write_case_inis(work, "frozen_plus", frozen_plus, params, pkmax)
    write_case_inis(work, "frozen_minus", frozen_minus, params, pkmax)
    write_case_inis(work, "fiducial", midpoint, params, pkmax)

    probe_cases = []
    for j, shape in enumerate(shapes):
        pp = fmid + PROBE_FRAC * shape
        pm = fmid - PROBE_FRAC * shape
        if np.min(pp) <= 0 or np.min(pm) <= 0:
            raise RuntimeError(f"Probe positivity failed for null direction {j}")
        for sign, arr in (("plus", pp), ("minus", pm)):
            case = f"probe_{j:02d}_{sign}"
            p = work / "data" / f"{case}.dat"
            write_psd(p, q, arr)
            write_case_inis(work, case, p, params, pkmax)
            probe_cases.append(case)

    nuisance_cases = ["fiducial"]
    for par in DERIVATIVE_STEPS:
        for sign_name, sign in (("plus", +1.0), ("minus", -1.0)):
            case = f"d_{par}_{sign_name}"
            write_case_inis(
                work, case, midpoint, params, pkmax,
                parameter_override(par, sign, params),
            )
            nuisance_cases.append(case)

    manifest = {
        "analysis": "idealized 3D LSS observability screen",
        "control_parameters": params,
        "mass_eV": float(params["mass_eV"]),
        "z_match": float(params["z_match"]),
        "redshifts": list(REDSHIFTS),
        "kmax_values_h_Mpc": list(KMAX_VALUES),
        "total_effective_volumes_Gpch3": list(TOTAL_VOLUMES_GPCH3),
        "tomography_volume_rule": (
            "total effective volume divided equally among the listed disjoint redshift bins"
        ),
        "probe_fractional_distortion": PROBE_FRAC,
        "final_fractional_distortion_cap": FINAL_FRAC,
        "null_dimension": int(shapes.shape[0]),
        "basis_dimension": int(basis.shape[0]),
        "moment_matrix_rank": int(np.linalg.matrix_rank(M)),
        "mass_over_Tnu_at_match": float(y),
        "midpoint_max_relative_difference_from_FD": rel_mid,
        "probe_cases": probe_cases,
        "nuisance_cases": nuisance_cases,
        "derivative_steps": DERIVATIVE_STEPS,
        "random_search_samples": RANDOM_SAMPLES,
        "random_seed": RANDOM_SEED,
        "scope": (
            "Ideal Gaussian matter-field cosmic variance only. This is a pre-survey "
            "observability gate, not a galaxy-clustering or weak-lensing likelihood."
        ),
        "limitations": [
            "linear matter power spectrum only",
            "no galaxy bias or redshift-space distortions",
            "no shot noise or selection function",
            "no weak-lensing shape noise or kernels",
            "no nonlinear covariance or baryonic uncertainty",
            "equal effective volume weighting across redshift bins is a diagnostic convention",
        ],
    }
    (work / "lss_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (work / "control_parameters.json").write_text(json.dumps(params, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def load_pk(results: Path, case: str, z: float):
    tag = ztag(z)
    p = unique_find(results, f"{case}_{tag}_00_pk.dat")
    a = np.loadtxt(p)
    if a.ndim == 1:
        a = a[None, :]
    if a.shape[1] < 2:
        raise RuntimeError(f"Unexpected P(k) table: {p}")
    k, pk = a[:, 0], a[:, 1]
    ok = np.isfinite(k) & np.isfinite(pk) & (k > 0) & (pk > 0)
    k, pk = k[ok], pk[ok]
    order = np.argsort(k)
    return k[order], pk[order]


def interp_log(kref: np.ndarray, k: np.ndarray, p: np.ndarray) -> np.ndarray:
    if np.any(p <= 0):
        raise RuntimeError("P(k) must be positive for log interpolation")
    if kref.min() < k.min() or kref.max() > k.max():
        mask = (kref >= k.min()) & (kref <= k.max())
        out = np.full_like(kref, np.nan, dtype=float)
        out[mask] = np.exp(np.interp(np.log(kref[mask]), np.log(k), np.log(p)))
        return out
    return np.exp(np.interp(np.log(kref), np.log(k), np.log(p)))


def common_pair(results: Path, case_plus: str, case_minus: str, z: float):
    kp, pp = load_pk(results, case_plus, z)
    km, pm = load_pk(results, case_minus, z)
    lo, hi = max(kp.min(), km.min()), min(kp.max(), km.max())
    mask = (kp >= lo) & (kp <= hi)
    k = kp[mask]
    if len(k) < 8:
        raise RuntimeError(f"Insufficient common P(k) grid at z={z}")
    pm_i = interp_log(k, km, pm)
    return k, pp[mask], pm_i


def trap_matrix(k: np.ndarray, X: np.ndarray) -> np.ndarray:
    nvec = X.shape[1]
    G = np.zeros((nvec, nvec))
    for i in range(nvec):
        for j in range(i, nvec):
            v = float(np.trapezoid(X[:, i] * X[:, j], k))
            G[i, j] = G[j, i] = v
    return G


def response_arrays(results: Path, manifest: dict, kmax: float):
    nd = int(manifest["null_dimension"])
    pars = list(DERIVATIVE_STEPS)
    nz = len(REDSHIFTS)
    Vbin = 1.0e9 / nz

    G = np.zeros((nd, nd))
    C = np.zeros((nd, len(pars)))
    F = np.zeros((len(pars), len(pars)))

    frozen_sn2 = 0.0
    per_z = {}

    for z in REDSHIFTS:
        k0, p0 = load_pk(results, "fiducial", z)
        sel0 = (k0 > 0) & (k0 <= kmax)
        k = k0[sel0]
        p0 = p0[sel0]
        if len(k) < 8:
            raise RuntimeError(f"Too few k points at z={z}, kmax={kmax}")

        B = np.zeros((len(k), nd))
        for j in range(nd):
            kp, pp = load_pk(results, f"probe_{j:02d}_plus", z)
            km, pm = load_pk(results, f"probe_{j:02d}_minus", z)
            pp_i = interp_log(k, kp, pp)
            pm_i = interp_log(k, km, pm)
            B[:, j] = pp_i - pm_i

        D = np.zeros((len(k), len(pars)))
        for a, par in enumerate(pars):
            kp, pp = load_pk(results, f"d_{par}_plus", z)
            km, pm = load_pk(results, f"d_{par}_minus", z)
            pp_i = interp_log(k, kp, pp)
            pm_i = interp_log(k, km, pm)
            step = float(DERIVATIVE_STEPS[par])
            D[:, a] = (pp_i - pm_i) / (2.0 * step)

        pref = math.sqrt(Vbin / (4.0 * math.pi**2))
        W = pref * k[:, None] / p0[:, None]
        G += trap_matrix(k, W * B)
        F += trap_matrix(k, W * D)
        for j in range(nd):
            for a in range(len(pars)):
                C[j, a] += float(np.trapezoid(
                    (pref * k * B[:, j] / p0) *
                    (pref * k * D[:, a] / p0), k
                ))

        kf, pfp, pfm = common_pair(results, "frozen_plus", "frozen_minus", z)
        sf = (kf > 0) & (kf <= kmax)
        kf, pfp, pfm = kf[sf], pfp[sf], pfm[sf]
        pbar = 0.5 * (pfp + pfm)
        sn2z = Vbin / (4.0 * math.pi**2) * float(np.trapezoid(
            kf**2 * ((pfp - pfm) / pbar)**2, kf
        ))
        frozen_sn2 += sn2z
        per_z[str(z)] = {
            "frozen_SN_per_sqrt_total_1_Gpch3_equal_bin_weight": float(math.sqrt(max(sn2z, 0.0))),
            "max_abs_relative_difference_percent": float(
                100.0 * np.max(np.abs(pfp - pfm) / pbar)
            ),
        }

    evals = np.linalg.eigvalsh(0.5 * (F + F.T))
    rank = int(np.linalg.matrix_rank(F, tol=max(np.max(np.abs(evals)), 1.0) * 1e-12))
    Fpinv = np.linalg.pinv(F, rcond=1e-10)
    raw_proj = G - C @ Fpinv @ C.T
    Gproj = 0.5 * (raw_proj + raw_proj.T)

    return {
        "G_fixed": G,
        "G_projected": Gproj,
        "cross": C,
        "F_nuisance": F,
        "fisher_eigenvalues": evals,
        "fisher_rank": rank,
        "parameters": pars,
        "frozen_sn2_per_total_1_Gpch3": frozen_sn2,
        "per_z": per_z,
    }


def top_eigvec(G: np.ndarray) -> np.ndarray:
    w, v = np.linalg.eigh(0.5 * (G + G.T))
    return v[:, int(np.argmax(w))]


def optimize(work: Path, results: Path, output: Path, pkmax: float,
             random_samples: int, seed: int) -> None:
    manifest = json.loads((work / "lss_manifest.json").read_text())
    params = manifest["control_parameters"]
    q = np.load(work / "q.npy")
    fmid = np.load(work / "fmid.npy")
    weights = np.load(work / "weights.npy")
    shapes = np.load(work / "null_shapes.npy")
    nd = shapes.shape[0]

    info = response_arrays(results, manifest, pkmax)
    Gfix = info["G_fixed"]
    Gproj = info["G_projected"]

    seeds = [top_eigvec(Gproj), top_eigvec(Gfix)]
    for j in range(nd):
        e = np.zeros(nd)
        e[j] = 1.0
        seeds.extend([e.copy(), -e.copy()])

    rng = np.random.default_rng(seed)
    batch = rng.normal(size=(random_samples, nd))
    batch /= np.linalg.norm(batch, axis=1, keepdims=True)

    best = None

    def evaluate(c: np.ndarray):
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(fmid, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            return None
        norm = 1.0 / maxrel
        scale = (FINAL_FRAC / PROBE_FRAC) * norm
        fixed_sn = float(scale * math.sqrt(max(float(c @ Gfix @ c), 0.0)))
        marg_sn = float(scale * math.sqrt(max(float(c @ Gproj @ c), 0.0)))
        return {
            "coeff": c.copy(),
            "norm": norm,
            "scale_from_probe_response": scale,
            "fixed_SN_per_sqrt_total_1_Gpch3": fixed_sn,
            "marginalized_SN_per_sqrt_total_1_Gpch3": marg_sn,
            "shape": raw * norm,
        }

    for c in list(seeds) + list(batch):
        r = evaluate(c)
        if r is None:
            continue
        if best is None or r["marginalized_SN_per_sqrt_total_1_Gpch3"] > best["marginalized_SN_per_sqrt_total_1_Gpch3"]:
            best = r

    if best is None:
        raise RuntimeError("No valid LSS hidden direction found")

    shape = best["shape"]
    fp = fmid + FINAL_FRAC * shape
    fm = fmid - FINAL_FRAC * shape
    if fp.min() <= 0 or fm.min() <= 0:
        raise RuntimeError("Optimized LSS pair failed positivity")

    plus = work / "data" / "optimized_plus.dat"
    minus = work / "data" / "optimized_minus.dat"
    write_psd(plus, q, fp)
    write_psd(minus, q, fm)
    write_case_inis(work, "optimized_plus", plus, params, max(KMAX_VALUES))
    write_case_inis(work, "optimized_minus", minus, params, max(KMAX_VALUES))

    stats = pair_stats(fmid, fp, fm, q, weights)
    pred = {}
    frozen = {}
    for V in TOTAL_VOLUMES_GPCH3:
        pred[str(V)] = {
            "fixed_SN": best["fixed_SN_per_sqrt_total_1_Gpch3"] * math.sqrt(V),
            "LCDM_marginalized_SN": best["marginalized_SN_per_sqrt_total_1_Gpch3"] * math.sqrt(V),
        }
        frozen[str(V)] = math.sqrt(max(info["frozen_sn2_per_total_1_Gpch3"], 0.0) * V)

    out = {
        "analysis": "LSS null-space optimization with ideal 3D Gaussian covariance",
        "optimization_kmax_h_Mpc": pkmax,
        "optimization_target": "LCDM-marginalized tomographic matter-power S/N",
        "coefficients": best["coeff"].tolist(),
        "normalization_to_unit_pointwise_shape": float(best["norm"]),
        "final_fractional_distortion_cap": FINAL_FRAC,
        "pair_stats": stats,
        "predicted_optimized": pred,
        "frozen_current_winner_SN": frozen,
        "fisher_parameters": info["parameters"],
        "fisher_rank": info["fisher_rank"],
        "fisher_eigenvalues": info["fisher_eigenvalues"].tolist(),
        "per_redshift_frozen_diagnostics": info["per_z"],
        "scope": manifest["scope"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def actual_signal_and_projection(results: Path, case_plus: str, case_minus: str,
                                 manifest: dict, kmax: float):
    pars = list(DERIVATIVE_STEPS)
    nz = len(REDSHIFTS)
    Vbin = 1.0e9 / nz

    F = np.zeros((len(pars), len(pars)))
    b = np.zeros(len(pars))
    sn2 = 0.0
    max_rel = 0.0
    per_z = {}

    for z in REDSHIFTS:
        k0, p0 = load_pk(results, "fiducial", z)
        sel0 = (k0 > 0) & (k0 <= kmax)
        k = k0[sel0]
        p0 = p0[sel0]

        kp, pp = load_pk(results, case_plus, z)
        km, pm = load_pk(results, case_minus, z)
        pp_i = interp_log(k, kp, pp)
        pm_i = interp_log(k, km, pm)
        h = pp_i - pm_i
        pbar = 0.5 * (pp_i + pm_i)
        max_rel = max(max_rel, float(np.max(np.abs(h) / pbar)))

        D = np.zeros((len(k), len(pars)))
        for a, par in enumerate(pars):
            k1, p1 = load_pk(results, f"d_{par}_plus", z)
            k2, p2 = load_pk(results, f"d_{par}_minus", z)
            step = float(DERIVATIVE_STEPS[par])
            D[:, a] = (interp_log(k, k1, p1) - interp_log(k, k2, p2)) / (2.0 * step)

        pref2 = Vbin / (4.0 * math.pi**2)
        sn2z = pref2 * float(np.trapezoid(k**2 * (h / p0)**2, k))
        sn2 += sn2z
        per_z[str(z)] = {
            "SN_per_sqrt_total_1_Gpch3_equal_bin_weight": float(math.sqrt(max(sn2z, 0.0))),
            "max_abs_relative_difference_percent": float(100.0 * np.max(np.abs(h) / pbar)),
        }

        for a in range(len(pars)):
            da = D[:, a]
            b[a] += pref2 * float(np.trapezoid(k**2 * h * da / (p0**2), k))
            for c in range(a, len(pars)):
                dc = D[:, c]
                v = pref2 * float(np.trapezoid(k**2 * da * dc / (p0**2), k))
                F[a, c] += v
                if c != a:
                    F[c, a] += v

    Fpinv = np.linalg.pinv(F, rcond=1e-10)
    absorbed = float(b @ Fpinv @ b)
    marg2 = max(float(sn2 - absorbed), 0.0)
    shifts = Fpinv @ b
    return {
        "fixed_SN2_per_total_1_Gpch3": float(sn2),
        "marginalized_SN2_per_total_1_Gpch3": float(marg2),
        "fraction_delta_chi2_retained": float(marg2 / sn2 if sn2 > 0 else 0.0),
        "fraction_signal_norm_absorbed_by_LCDM": float(absorbed / sn2 if sn2 > 0 else 0.0),
        "max_abs_relative_difference_percent": float(100.0 * max_rel),
        "best_fit_parameter_shifts": {
            p: float(shifts[i]) for i, p in enumerate(pars)
        },
        "fisher_rank": int(np.linalg.matrix_rank(F, tol=max(float(np.max(np.abs(F))), 1.0) * 1e-12)),
        "per_redshift": per_z,
    }


def summarize(work: Path, screen_results: Path, optimized_results: Path,
              output: Path) -> None:
    manifest = json.loads((work / "lss_manifest.json").read_text())

    class CombinedPath:
        def __init__(self, roots):
            self.roots = roots
        def rglob(self, name):
            out = []
            for r in self.roots:
                out.extend(r.rglob(name))
            return out

    all_results = CombinedPath([screen_results, optimized_results])

    def load_combined(case, z):
        tag = ztag(z)
        hits = all_results.rglob(f"{case}_{tag}_00_pk.dat")
        if len(hits) != 1:
            raise RuntimeError(f"Expected one {case}_{tag}_00_pk.dat, found {hits}")
        a = np.loadtxt(hits[0])
        if a.ndim == 1:
            a = a[None, :]
        k, p = a[:, 0], a[:, 1]
        ok = np.isfinite(k) & np.isfinite(p) & (k > 0) & (p > 0)
        order = np.argsort(k[ok])
        return k[ok][order], p[ok][order]

    global load_pk
    original_load_pk = load_pk
    def _loader(_unused, case, z):
        return load_combined(case, z)
    load_pk = _loader
    try:
        out = {
            "analysis": "validated idealized 3D LSS observability campaign",
            "control_parameters": manifest["control_parameters"],
            "volume_convention": manifest["tomography_volume_rule"],
            "scope": manifest["scope"],
            "limitations": manifest["limitations"],
            "results_by_kmax": {},
        }
        for kmax in KMAX_VALUES:
            frozen = actual_signal_and_projection(
                Path("."), "frozen_plus", "frozen_minus", manifest, kmax
            )
            opt = actual_signal_and_projection(
                Path("."), "optimized_plus", "optimized_minus", manifest, kmax
            )
            entry = {
                "frozen_current_winner": frozen,
                "LSS_optimized_hidden_mode": opt,
                "volumes": {},
            }
            for V in TOTAL_VOLUMES_GPCH3:
                entry["volumes"][str(V)] = {
                    "frozen_fixed_SN": math.sqrt(max(frozen["fixed_SN2_per_total_1_Gpch3"] * V, 0.0)),
                    "frozen_LCDM_marginalized_SN": math.sqrt(max(frozen["marginalized_SN2_per_total_1_Gpch3"] * V, 0.0)),
                    "optimized_fixed_SN": math.sqrt(max(opt["fixed_SN2_per_total_1_Gpch3"] * V, 0.0)),
                    "optimized_LCDM_marginalized_SN": math.sqrt(max(opt["marginalized_SN2_per_total_1_Gpch3"] * V, 0.0)),
                }
            out["results_by_kmax"][str(kmax)] = entry
    finally:
        load_pk = original_load_pk

    key = out["results_by_kmax"]["0.3"]["volumes"]["50.0"]["optimized_LCDM_marginalized_SN"]
    out["decision_gate"] = {
        "reference": "optimized LCDM-marginalized S/N at kmax=0.3 h/Mpc and total effective volume 50 (Gpc/h)^3",
        "SN": float(key),
        "classification": (
            "promising_for_realistic_survey_forecast" if key >= 1.0
            else "weak_in_idealized_LSS_screen"
        ),
        "note": (
            "This classification only decides the next numerical step. It is not a detection claim."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["prepare", "optimize", "summarize"])
    ap.add_argument("--source", type=Path)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path)
    ap.add_argument("--optimized-results", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--pkmax", type=float, default=0.30)
    ap.add_argument("--random-samples", type=int, default=RANDOM_SAMPLES)
    ap.add_argument("--seed", type=int, default=RANDOM_SEED)
    args = ap.parse_args()

    if args.mode == "prepare":
        if args.source is None:
            ap.error("--source is required in prepare mode")
        prepare(args.source, args.work, args.pkmax)
    elif args.mode == "optimize":
        if args.results is None or args.output is None:
            ap.error("--results and --output are required in optimize mode")
        optimize(args.work, args.results, args.output, args.pkmax,
                 args.random_samples, args.seed)
    else:
        if args.results is None or args.optimized_results is None or args.output is None:
            ap.error("--results, --optimized-results and --output are required in summarize mode")
        summarize(args.work, args.results, args.optimized_results, args.output)


if __name__ == "__main__":
    main()
