#!/usr/bin/env python3
"""Planck-anchored CLASS forecast for stress-energy-matched nonthermal relics.

This script does not fit observational data. It uses a Planck 2018 reference
cosmology and the standard CLASS Fermi-Dirac relic-neutrino spectrum as a
baseline, constructs two positive nonthermal distributions with identical
n, rho and P at z_match, writes CLASS ncdm PSD/INI files, and post-processes
tensor B-mode spectra produced by CLASS with tensor_method=exact.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

PLANCK = {
    "H0": 67.36,
    "omega_b": 0.02237,
    "omega_cdm": 0.1200,
    "A_s": float(np.exp(3.044) * 1e-10),
    "n_s": 0.9649,
    "tau_reio": 0.054,
    "N_ur": 2.0328,
    "N_ncdm": 1,
    "m_ncdm": 0.06,
    "T_ncdm": 0.71611,
}
TCMB_K = 2.7255
KB_EV_K = 8.617333262e-5
CLASS_COMMIT = "e85808324f51fc694d12e3ed7439552a3c3f9540"
CLASS_FD_NORM = 2.0 / (2.0 * np.pi) ** 3


def trap(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def fd_occupancy(q):
    return 1.0 / (np.exp(q) + 1.0)


def fermi_dirac(q):
    """CLASS default zero-chemical-potential particle+antiparticle PSD."""
    return CLASS_FD_NORM * fd_occupancy(q)


def angular_A(x):
    mu, wg = np.polynomial.legendre.leggauss(64)
    shape = wg * (1.0 - mu * mu) ** 2
    return np.cos(np.outer(np.atleast_1d(x), mu)) @ shape


def construct_pair(q, z_match=1100.0, max_fractional_distortion=0.30):
    occ = fd_occupancy(q)
    f0 = CLASS_FD_NORM * occ
    a = 1.0 / (1.0 + z_match)
    Tnu0_eV = PLANCK["T_ncdm"] * TCMB_K * KB_EV_K
    y = a * PLANCK["m_ncdm"] / Tnu0_eV
    eps = np.sqrt(q * q + y * y)

    centers = np.array([0.45, 0.85, 1.35, 2.0, 2.8, 3.8, 5.0, 6.5, 8.5, 11.0])
    widths = np.array([0.28, 0.32, 0.38, 0.46, 0.55, 0.65, 0.78, 0.95, 1.20, 1.55])
    basis = np.array([
        f0 * (1.0 - occ) * np.exp(-0.5 * ((q - c) / s) ** 2)
        for c, s in zip(centers, widths)
    ])

    weights = np.vstack([q ** 2, q ** 2 * eps, q ** 4 / eps])
    M = np.array([[trap(w * b, q) for b in basis] for w in weights])
    _, _, vh = np.linalg.svd(M, full_matrices=True)
    rank = np.linalg.matrix_rank(M)
    N = vh[rank:].T
    if N.shape[1] == 0:
        raise RuntimeError("No null space found for matched-moment construction")

    dbasis = np.gradient(basis, q, axis=1)
    v = q / eps
    w_resp = q ** 5 / eps
    delays = np.linspace(0.0, 30.0, 180)
    R = np.empty((len(delays), len(basis)))
    for i, t in enumerate(delays):
        Avec = angular_A(v * t)
        R[i] = [trap(w_resp * db * Avec, q) for db in dbasis]
    _, _, vhr = np.linalg.svd(R @ N, full_matrices=False)
    delta = (N @ vhr[0]) @ basis

    mask = np.abs(delta) > 1e-24
    alpha_pos = np.min(f0[mask] / np.abs(delta[mask]))
    alpha_frac = max_fractional_distortion / np.max(np.abs(delta[mask]) / f0[mask])
    alpha = min(0.98 * alpha_pos, alpha_frac)
    fp = f0 + alpha * delta
    fm = f0 - alpha * delta
    if fp.min() <= 0 or fm.min() <= 0:
        raise RuntimeError("Positivity failed")

    def moms(f):
        return np.array([trap(w * f, q) for w in weights])

    m0, mp, mm = moms(f0), moms(fp), moms(fm)
    rel_match = np.abs(mp - mm) / np.maximum(0.5 * (np.abs(mp) + np.abs(mm)), 1e-300)
    summary = {
        "z_match": z_match,
        "a_match": a,
        "m_ncdm_eV": PLANCK["m_ncdm"],
        "T_ncdm_over_Tcmb": PLANCK["T_ncdm"],
        "mass_over_Tnu_at_match": y,
        "class_fd_normalization": CLASS_FD_NORM,
        "max_relative_moment_mismatch": float(rel_match.max()),
        "relative_n_mismatch": float(rel_match[0]),
        "relative_rho_mismatch": float(rel_match[1]),
        "relative_P_mismatch": float(rel_match[2]),
        "max_fractional_distortion_plus": float(np.max(np.abs(fp - f0) / f0)),
        "max_fractional_distortion_minus": float(np.max(np.abs(fm - f0) / f0)),
        "distribution_L2_separation_over_FD": float(np.sqrt(trap((fp - fm) ** 2, q) / trap(f0 ** 2, q))),
        "n_baseline": float(m0[0]),
        "rho_baseline_dimensionless": float(m0[1]),
        "P3_baseline_dimensionless": float(m0[2]),
    }
    return f0, fp, fm, summary


def write_ini(path, psd_path, root, lmax=600, r=0.01):
    p = PLANCK
    text = f"""# Planck-anchored tensor forecast with custom ncdm PSD.
output = tCl,pCl
modes = t
tensor_method = exact
lensing = no
H0 = {p['H0']}
omega_b = {p['omega_b']}
omega_cdm = {p['omega_cdm']}
A_s = {p['A_s']:.12e}
tau_reio = {p['tau_reio']}
N_ur = {p['N_ur']}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd_path}
m_ncdm = {p['m_ncdm']}
T_ncdm = {p['T_ncdm']}
deg_ncdm = 1.0
r = {r}
n_t = 0.0
k_pivot = 0.05
l_max_tensors = {lmax}
root = {root}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


def prepare(outdir: Path, z_match: float, frac: float):
    outdir.mkdir(parents=True, exist_ok=True)
    q = np.linspace(0.0, 20.0, 4000)
    f0, fp, fm, summary = construct_pair(q, z_match, frac)
    np.savetxt(outdir / "psd_fd_reference.dat", np.column_stack([q, f0]), fmt="%.12e")
    np.savetxt(outdir / "psd_plus.dat", np.column_stack([q, fp]), fmt="%.12e")
    np.savetxt(outdir / "psd_minus.dat", np.column_stack([q, fm]), fmt="%.12e")
    np.savetxt(outdir / "matched_distribution_pair.csv", np.column_stack([q, f0, fp, fm]), delimiter=",", header="q,f_FD,f_plus,f_minus", comments="")
    write_ini(outdir / "class_plus.ini", (outdir / "psd_plus.dat").resolve(), outdir / "plus_")
    write_ini(outdir / "class_minus.ini", (outdir / "psd_minus.dat").resolve(), outdir / "minus_")
    summary["class_commit"] = CLASS_COMMIT
    summary["planck_reference"] = "Planck 2018 TT,TE,EE+lowE+lensing central values"
    (outdir / "pair_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


def discover_cl(outdir: Path, stem: str, explicit: Path | None):
    if explicit is not None:
        return explicit
    candidates = sorted(outdir.glob(f"{stem}*cl*.dat"))
    if not candidates:
        candidates = sorted(Path(".").glob(f"**/{stem}*cl*.dat"))
    candidates = [p for p in candidates if "lensed" not in p.name]
    if not candidates:
        produced = "\n".join(str(p) for p in sorted(outdir.glob("*")))
        raise FileNotFoundError(f"No CLASS Cl file found for {stem}. Produced files:\n{produced}")
    print(f"Using {stem} CLASS spectrum: {candidates[0]}")
    return candidates[0]


def read_class_cl(path: Path):
    arr = np.loadtxt(path)
    if arr.shape[1] < 5:
        raise RuntimeError(f"Expected CLASS cl.dat with at least 5 columns, got {arr.shape[1]} in {path}")
    # CLASS headers for CMB Cl output: l, TT, EE, TE, BB (additional columns may follow).
    return arr[:, 0].astype(int), arr[:, 4]


def summarize(outdir: Path, plus_cl: Path | None, minus_cl: Path | None):
    plus_cl = discover_cl(outdir, "plus", plus_cl)
    minus_cl = discover_cl(outdir, "minus", minus_cl)
    ellp, bbp = read_class_cl(plus_cl)
    ellm, bbm = read_class_cl(minus_cl)
    if not np.array_equal(ellp, ellm):
        raise RuntimeError("CLASS ell grids differ")
    ell = ellp
    denom = 0.5 * (np.abs(bbp) + np.abs(bbm))
    rel = np.zeros_like(denom)
    mask = denom > max(np.max(denom) * 1e-10, 1e-300)
    rel[mask] = 100.0 * (bbp[mask] - bbm[mask]) / denom[mask]
    absrel = np.abs(rel)
    imax = np.argmax(absrel)

    np.savetxt(outdir / "class_bmode_comparison.csv", np.column_stack([ell, bbp, bbm, rel]), delimiter=",",
               header="ell,Dl_BB_plus,Dl_BB_minus,symmetric_relative_difference_percent", comments="")
    summary = json.loads((outdir / "pair_summary.json").read_text())
    summary.update({
        "plus_class_output": str(plus_cl),
        "minus_class_output": str(minus_cl),
        "max_abs_BB_relative_difference_percent": float(absrel[imax]),
        "ell_at_max_abs_BB_difference": int(ell[imax]),
        "Dl_BB_plus_at_max": float(bbp[imax]),
        "Dl_BB_minus_at_max": float(bbm[imax]),
        "r_used_for_normalization": 0.01,
    })
    (outdir / "class_forecast_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7.0, 4.5))
        plt.semilogy(ell, np.abs(bbp), label=r"$F_+$")
        plt.semilogy(ell, np.abs(bbm), label=r"$F_-$", linestyle="--")
        plt.xlabel(r"$\ell$")
        plt.ylabel(r"$D_\ell^{BB}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "class_bmode_spectra.pdf")
        plt.close()
        plt.figure(figsize=(7.0, 4.5))
        plt.plot(ell[mask], rel[mask])
        plt.axhline(0.0, linewidth=0.8)
        plt.xlabel(r"$\ell$")
        plt.ylabel(r"$\Delta D_\ell^{BB}/D_\ell^{BB}$ (\%)")
        plt.tight_layout()
        plt.savefig(outdir / "class_bmode_relative_difference.pdf")
        plt.close()
    except Exception as exc:
        print(f"Plotting skipped: {exc}")
    print(json.dumps(summary, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["prepare", "summarize"], required=True)
    ap.add_argument("--outdir", default="observational_forecast_output")
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--max-fractional-distortion", type=float, default=0.30)
    ap.add_argument("--plus-cl", type=Path)
    ap.add_argument("--minus-cl", type=Path)
    args = ap.parse_args()
    out = Path(args.outdir)
    if args.mode == "prepare":
        prepare(out, args.z_match, args.max_fractional_distortion)
    else:
        summarize(out, args.plus_cl, args.minus_cl)


if __name__ == "__main__":
    main()
