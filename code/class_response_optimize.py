#!/usr/bin/env python3
"""Optimize stress-energy-matched relic distributions for CMB tensor B-mode response.

This is a diagnostic forecast, not a data fit. For a fixed relic mass, the script:
1. builds a smooth basis of nonthermal deformations around the CLASS Fermi-Dirac PSD;
2. projects that basis into the exact numerical null space of n, rho and P at z_match;
3. probes each null direction with small +/- deformations and CLASS tensor spectra;
4. uses the CLASS response Jacobian to optimize the null direction for an ideal
   full-sky cosmic-variance-limited B-mode S/N over a chosen ell range;
5. rescales the optimized direction to a fixed positivity/fractional-distortion cap;
6. writes a final +/- pair for independent CLASS confirmation.

The optimization target is deliberately an optimistic upper bound. Instrument noise,
foregrounds, lensing residuals and parameter degeneracies are not included.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import numpy as np

H0 = 67.36
OMEGA_B = 0.02237
OMEGA_CDM = 0.1200
A_S = float(np.exp(3.044) * 1e-10)
TAU_REIO = 0.054
N_UR = 2.0328
T_NCDM = 0.71611
TCMB_K = 2.7255
KB_EV_K = 8.617333262e-5
CLASS_COMMIT = "e85808324f51fc694d12e3ed7439552a3c3f9540"
CLASS_FD_NORM = 2.0 / (2.0 * np.pi) ** 3

CENTERS = np.array([0.45, 0.85, 1.35, 2.0, 2.8, 3.8, 5.0, 6.5, 8.5, 11.0])
WIDTHS = np.array([0.28, 0.32, 0.38, 0.46, 0.55, 0.65, 0.78, 0.95, 1.20, 1.55])


def trap(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def fd_occ(q):
    return 1.0 / (np.exp(q) + 1.0)


def f0_class(q):
    return CLASS_FD_NORM * fd_occ(q)


def kinetic_objects(q, mass_eV, z_match):
    f0 = f0_class(q)
    occ = fd_occ(q)
    a = 1.0 / (1.0 + z_match)
    tnu0_eV = T_NCDM * TCMB_K * KB_EV_K
    y = a * mass_eV / tnu0_eV
    eps = np.sqrt(q * q + y * y)
    weights = np.vstack([q**2, q**2 * eps, q**4 / eps])
    basis = np.array([
        f0 * (1.0 - occ) * np.exp(-0.5 * ((q - c) / s) ** 2)
        for c, s in zip(CENTERS, WIDTHS)
    ])
    M = np.array([[trap(w * b, q) for b in basis] for w in weights])
    _, sv, vh = np.linalg.svd(M, full_matrices=True)
    tol = max(M.shape) * np.max(sv) * np.finfo(float).eps
    rank = int(np.sum(sv > tol))
    N = vh[rank:].T
    if N.shape[1] == 0:
        raise RuntimeError("Null space is empty")

    shapes = []
    for j in range(N.shape[1]):
        d = N[:, j] @ basis
        rel = np.max(np.abs(d) / np.maximum(f0, 1e-300))
        if not np.isfinite(rel) or rel <= 0:
            raise RuntimeError("Invalid null direction normalization")
        shapes.append(d / rel)
    shapes = np.asarray(shapes)
    return f0, weights, basis, N, shapes, y, M


def moments(f, q, weights):
    return np.array([trap(w * f, q) for w in weights])


def pair_stats(f0, fp, fm, q, weights):
    mp, mm = moments(fp, q, weights), moments(fm, q, weights)
    denom = np.maximum(0.5 * (np.abs(mp) + np.abs(mm)), 1e-300)
    mismatch = np.abs(mp - mm) / denom
    return {
        "relative_n_mismatch": float(mismatch[0]),
        "relative_rho_mismatch": float(mismatch[1]),
        "relative_P_mismatch": float(mismatch[2]),
        "max_relative_moment_mismatch": float(mismatch.max()),
        "max_fractional_distortion_plus": float(np.max(np.abs(fp - f0) / np.maximum(f0, 1e-300))),
        "max_fractional_distortion_minus": float(np.max(np.abs(fm - f0) / np.maximum(f0, 1e-300))),
        "distribution_L2_separation_over_FD": float(np.sqrt(trap((fp - fm)**2, q) / trap(f0**2, q))),
    }


def write_ini(path, psd_path, root, mass_eV, lmax=350, r=0.01, deg=1.0):
    text = f"""# Response-optimization CLASS tensor forecast.
output = tCl,pCl
modes = t
tensor_method = exact
lensing = no
H0 = {H0}
omega_b = {OMEGA_B}
omega_cdm = {OMEGA_CDM}
A_s = {A_S:.12e}
tau_reio = {TAU_REIO}
N_ur = {N_UR}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd_path}
m_ncdm = {mass_eV}
T_ncdm = {T_NCDM}
deg_ncdm = {deg}
r = {r}
n_t = 0.0
k_pivot = 0.05
l_max_tensors = {lmax}
root = {root}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


def write_psd(path, q, f):
    if np.min(f) <= 0:
        raise RuntimeError(f"Non-positive PSD in {path}")
    np.savetxt(path, np.column_stack([q, f]), fmt="%.14e")


def read_cl(path):
    arr = np.loadtxt(path)
    if arr.ndim != 2 or arr.shape[1] < 5:
        raise RuntimeError(f"Unexpected CLASS Cl format: {path}")
    return arr[:, 0].astype(int), arr[:, 4]


def find_cl(outdir, prefix):
    cand = sorted(outdir.glob(f"{prefix}*cl*.dat"))
    cand = [p for p in cand if "lensed" not in p.name]
    if not cand:
        raise FileNotFoundError(f"No Cl file for prefix={prefix} in {outdir}")
    return cand[0]


def relative_pair(plus_path, minus_path):
    ep, bp = read_cl(plus_path)
    em, bm = read_cl(minus_path)
    if not np.array_equal(ep, em):
        raise RuntimeError("ell grids differ")
    den = 0.5 * (np.abs(bp) + np.abs(bm))
    mask = den > max(float(np.max(den)) * 1e-12, 1e-300)
    rel = np.zeros_like(den)
    rel[mask] = (bp[mask] - bm[mask]) / den[mask]
    return ep, bp, bm, rel, mask


def prepare_probes(outdir, mass_eV, z_match, probe_frac, lmax):
    outdir.mkdir(parents=True, exist_ok=True)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, N, shapes, y, M = kinetic_objects(q, mass_eV, z_match)

    np.savetxt(outdir / "q_grid.csv", q, delimiter=",", header="q", comments="")
    np.save(outdir / "null_shapes.npy", shapes)
    np.save(outdir / "q.npy", q)
    np.save(outdir / "f0.npy", f0)
    np.save(outdir / "weights.npy", weights)

    for j, shape in enumerate(shapes):
        fp = f0 + probe_frac * shape
        fm = f0 - probe_frac * shape
        if np.min(fp) <= 0 or np.min(fm) <= 0:
            raise RuntimeError(f"Probe positivity failed for direction {j}")
        pp = outdir / f"probe_{j:02d}_plus.dat"
        pm = outdir / f"probe_{j:02d}_minus.dat"
        write_psd(pp, q, fp)
        write_psd(pm, q, fm)
        write_ini(outdir / f"probe_{j:02d}_plus.ini", pp.resolve(), outdir / f"probe_{j:02d}_plus_", mass_eV, lmax=lmax)
        write_ini(outdir / f"probe_{j:02d}_minus.ini", pm.resolve(), outdir / f"probe_{j:02d}_minus_", mass_eV, lmax=lmax)

    meta = {
        "mass_eV": mass_eV,
        "z_match": z_match,
        "mass_over_Tnu_at_match": y,
        "probe_fractional_distortion": probe_frac,
        "null_dimension": int(shapes.shape[0]),
        "basis_dimension": int(basis.shape[0]),
        "moment_matrix_rank": int(np.linalg.matrix_rank(M)),
        "lmax": lmax,
        "class_commit": CLASS_COMMIT,
        "optimization_target": "ideal full-sky cosmic-variance-limited BB S/N",
    }
    (outdir / "probe_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


def predicted_metrics(rel, ell, lmin, lmax):
    sel = (ell >= lmin) & (ell <= lmax)
    rr = rel[sel]
    ee = ell[sel]
    sn = float(np.sqrt(np.sum(0.5 * (2.0 * ee + 1.0) * rr**2)))
    return sn


def optimize(outdir, final_frac, lmin, lmax, random_samples, seed):
    meta = json.loads((outdir / "probe_meta.json").read_text())
    probe_frac = float(meta["probe_fractional_distortion"])
    mass_eV = float(meta["mass_eV"])
    z_match = float(meta["z_match"])
    shapes = np.load(outdir / "null_shapes.npy")
    q = np.load(outdir / "q.npy")
    f0 = np.load(outdir / "f0.npy")
    weights = np.load(outdir / "weights.npy")
    nd = shapes.shape[0]

    rel_cols = []
    ell0 = None
    for j in range(nd):
        ell, _, _, rel, _ = relative_pair(find_cl(outdir, f"probe_{j:02d}_plus_"), find_cl(outdir, f"probe_{j:02d}_minus_"))
        if ell0 is None:
            ell0 = ell
        elif not np.array_equal(ell0, ell):
            raise RuntimeError("Probe ell grids differ")
        rel_cols.append(rel)
    R = np.column_stack(rel_cols)
    ell = ell0
    sel = (ell >= lmin) & (ell <= lmax)
    W = np.sqrt(0.5 * (2.0 * ell[sel] + 1.0))[:, None]
    WR = W * R[sel]

    # SVD seed direction, then deterministic random search with the exact pointwise
    # fractional-distortion normalization included in the predicted objective.
    _, _, vh = np.linalg.svd(WR, full_matrices=False)
    candidates = [vh[0]]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0; candidates.append(e)
        candidates.append(-e)

    rng = np.random.default_rng(seed)
    batch = rng.normal(size=(random_samples, nd))
    batch /= np.linalg.norm(batch, axis=1, keepdims=True)

    best = None
    best_sn = -np.inf

    def score(c):
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            return None
        norm = 1.0 / maxrel
        pred_rel = (final_frac / probe_frac) * norm * (R @ c)
        sn = predicted_metrics(pred_rel, ell, lmin, lmax)
        return sn, norm, pred_rel, raw * norm

    for c in candidates:
        s = score(c)
        if s and s[0] > best_sn:
            best_sn = s[0]; best = (c.copy(),) + s[1:]

    # Evaluate random candidates in chunks to keep memory small.
    for c in batch:
        s = score(c)
        if s and s[0] > best_sn:
            best_sn = s[0]; best = (c.copy(),) + s[1:]

    coeff, norm, pred_rel, shape = best
    fp = f0 + final_frac * shape
    fm = f0 - final_frac * shape
    if np.min(fp) <= 0 or np.min(fm) <= 0:
        raise RuntimeError("Optimized final pair failed positivity")

    stats = pair_stats(f0, fp, fm, q, weights)
    write_psd(outdir / "final_plus.dat", q, fp)
    write_psd(outdir / "final_minus.dat", q, fm)
    write_ini(outdir / "final_plus.ini", (outdir / "final_plus.dat").resolve(), outdir / "final_plus_", mass_eV, lmax=int(meta["lmax"]))
    write_ini(outdir / "final_minus.ini", (outdir / "final_minus.dat").resolve(), outdir / "final_minus_", mass_eV, lmax=int(meta["lmax"]))
    np.savetxt(outdir / "optimized_pair.csv", np.column_stack([q, f0, fp, fm, shape]), delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape", comments="")
    np.savetxt(outdir / "predicted_optimized_response.csv", np.column_stack([ell, pred_rel]), delimiter=",",
               header="ell,predicted_symmetric_relative_difference", comments="")

    summary = dict(meta)
    summary.update(stats)
    summary.update({
        "final_fractional_distortion_cap": final_frac,
        "optimization_ell_min": lmin,
        "optimization_ell_max": lmax,
        "random_search_samples": random_samples,
        "random_seed": seed,
        "predicted_ideal_fullsky_cv_SN": float(best_sn),
        "optimized_coefficients": coeff.tolist(),
        "pointwise_normalization_factor": float(norm),
    })
    (outdir / "optimization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def summarize_final(outdir, lmin, lmax, tag):
    plus = find_cl(outdir, "final_plus_")
    minus = find_cl(outdir, "final_minus_")
    ell, bp, bm, rel, mask = relative_pair(plus, minus)
    sel = (ell >= lmin) & (ell <= lmax) & mask
    if not np.any(sel):
        raise RuntimeError("No valid multipoles in requested range")
    idxs = np.where(sel)[0]
    im = idxs[np.argmax(np.abs(rel[sel]))]
    sn = predicted_metrics(rel, ell, lmin, lmax)

    def at(target):
        i = int(np.argmin(np.abs(ell - target)))
        return int(ell[i]), float(100.0 * rel[i])

    ell80, p80 = at(80)
    ell233, p233 = at(233)
    out_csv = outdir / f"final_bmode_comparison_{tag}.csv"
    np.savetxt(out_csv, np.column_stack([ell, bp, bm, 100.0 * rel]), delimiter=",",
               header="ell,Dl_BB_plus,Dl_BB_minus,symmetric_relative_difference_percent", comments="")

    base = json.loads((outdir / "optimization_summary.json").read_text())
    summary = dict(base)
    summary.update({
        "summary_tag": tag,
        "actual_ideal_fullsky_cv_SN": float(sn),
        "actual_max_abs_relative_difference_percent": float(100.0 * abs(rel[im])),
        "ell_at_actual_max": int(ell[im]),
        "relative_difference_percent_at_ell80": p80,
        "ell_used_for_ell80": ell80,
        "relative_difference_percent_at_ell233": p233,
        "ell_used_for_ell233": ell233,
        "plus_class_output": str(plus),
        "minus_class_output": str(minus),
    })
    (outdir / f"final_summary_{tag}.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def aggregate(root, lmin, lmax):
    rows = []
    for d in sorted(root.glob("mass_*")):
        p = d / "final_summary_default.json"
        if not p.exists():
            continue
        s = json.loads(p.read_text())
        rows.append({
            "directory": d.name,
            "mass_eV": s["mass_eV"],
            "mass_over_Tnu_at_match": s["mass_over_Tnu_at_match"],
            "moment_mismatch": s["max_relative_moment_mismatch"],
            "predicted_cv_SN": s["predicted_ideal_fullsky_cv_SN"],
            "actual_cv_SN": s["actual_ideal_fullsky_cv_SN"],
            "max_abs_percent": s["actual_max_abs_relative_difference_percent"],
            "ell_at_max": s["ell_at_actual_max"],
            "percent_at_ell80": s["relative_difference_percent_at_ell80"],
            "percent_at_ell233": s["relative_difference_percent_at_ell233"],
        })
    if not rows:
        raise RuntimeError("No default final summaries found")
    rows.sort(key=lambda r: r["mass_eV"])
    with (root / "mass_sweep_summary.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    best = max(rows, key=lambda r: r["actual_cv_SN"])
    payload = {
        "optimization_ell_min": lmin,
        "optimization_ell_max": lmax,
        "best_directory": best["directory"],
        "best_mass_eV": best["mass_eV"],
        "best_actual_cv_SN": best["actual_cv_SN"],
        "best_max_abs_percent": best["max_abs_percent"],
        "rows": rows,
    }
    (root / "best_candidate.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["prepare-probes", "optimize", "summarize-final", "aggregate"], required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mass-eV", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.02)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--lmin", type=int, default=20)
    ap.add_argument("--lmax", type=int, default=300)
    ap.add_argument("--class-lmax", type=int, default=350)
    ap.add_argument("--random-samples", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=271828)
    ap.add_argument("--tag", default="default")
    args = ap.parse_args()

    if args.mode == "prepare-probes":
        prepare_probes(args.outdir, args.mass_eV, args.z_match, args.probe_frac, args.class_lmax)
    elif args.mode == "optimize":
        optimize(args.outdir, args.final_frac, args.lmin, args.lmax, args.random_samples, args.seed)
    elif args.mode == "summarize-final":
        summarize_final(args.outdir, args.lmin, args.lmax, args.tag)
    else:
        aggregate(args.outdir, args.lmin, args.lmax)


if __name__ == "__main__":
    main()
