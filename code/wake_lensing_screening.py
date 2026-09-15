#!/usr/bin/env python3
r"""Screen hidden-state sensitivity of a neutrino-wake weak-lensing dipole.

This is deliberately a *screening* calculation, not a final survey forecast.
It asks a narrow question before we build a full lensing pipeline:

    If a fiducial neutrino wake were detected in weak lensing with S/N = S0,
    how much of that significance would separate the stress-energy-matched
    hidden states F_+ and F_- used in the paper?

The kinetic dependence is taken from the same nonrelativistic Landau-pole
factor used by the wake calculation,

    Y_F(u,z) \propto |u| F(q_res),
    q_res = m |u| / [c T_nu,0 (1+z)].

All lensing geometry, halo weighting and survey noise multiply this kinetic
factor in the first screening approximation.  They therefore cancel in the
fractional pair separation

    R_pair = |A(F_+) - A(F_-)| / A(F_FD).

The implied hidden-state significance is S_hidden = S0 * R_pair.  The script
reports the fiducial wake-lensing significance required for a 1-sigma or
3-sigma hidden-state separation.  It evaluates two velocity averages:

  * coherent_mean: <|x| F(alpha |x|)>, appropriate to a coherently oriented
    stacked dipole amplitude;
  * rms_power: sqrt(<x^2 F(alpha |x|)^2>), a power-like cross-check.

A large and stable R_pair in both summaries is a reason to proceed to the full
weak-lensing forecast.  A small value means the channel is not worth adding to
the manuscript.  No manuscript claim should be based on this screening layer
alone.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import class_response_optimize as cro

C_KMS = 299792.458
TNU0_EV = cro.T_NCDM * cro.TCMB_K * cro.KB_EV_K

# Reference direction used as the stable leading hidden-state direction in the
# validated DESI-wake multi-tracer calculation.  A final pair CSV can be passed
# with --pair-csv to avoid relying on this reconstruction.
CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def reconstruct_reference_pair(mass_eV: float, z_match: float, frac: float):
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, _, _, shapes, _, _ = cro.kinetic_objects(q, mass_eV, z_match)
    if shapes.shape[0] != CREF.size:
        raise RuntimeError(
            f"Reference coefficient length {CREF.size} does not match null-space "
            f"dimension {shapes.shape[0]}. Use --pair-csv with the validated pair."
        )
    coeff = CREF / np.linalg.norm(CREF)
    raw = coeff @ shapes
    maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
    if not np.isfinite(maxrel) or maxrel <= 0:
        raise RuntimeError("Invalid reference hidden-state direction")
    shape = raw / maxrel
    fp = f0 + frac * shape
    fm = f0 - frac * shape
    if np.min(fp) <= 0 or np.min(fm) <= 0:
        raise RuntimeError("Reference pair violates positivity")
    stats = cro.pair_stats(f0, fp, fm, q, weights)
    return q, f0, fp, fm, shape, stats


def load_pair(path: Path):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    required = {"q", "f_FD", "f_plus", "f_minus"}
    if arr.dtype.names is None or not required.issubset(arr.dtype.names):
        raise RuntimeError(
            f"{path} must contain columns q,f_FD,f_plus,f_minus; "
            f"found {arr.dtype.names}"
        )
    q = np.asarray(arr["q"], float)
    f0 = np.asarray(arr["f_FD"], float)
    fp = np.asarray(arr["f_plus"], float)
    fm = np.asarray(arr["f_minus"], float)
    if np.any(np.diff(q) <= 0):
        raise RuntimeError("pair CSV q grid must be strictly increasing")
    if min(np.min(f0), np.min(fp), np.min(fm)) <= 0:
        raise RuntimeError("pair CSV contains non-positive distribution values")
    shape = (fp - fm) / 2.0
    return q, f0, fp, fm, shape


def gh_normal(n: int):
    x0, w0 = np.polynomial.hermite.hermgauss(n)
    return np.sqrt(2.0) * x0, w0 / np.sqrt(np.pi)


def resonant_q(abs_x, mass_eV: float, z: float, sigma_kms: float):
    return mass_eV * abs_x * sigma_kms / (C_KMS * TNU0_EV * (1.0 + z))


def weighted_quantile(values, weights, quantile: float):
    order = np.argsort(values)
    vv = np.asarray(values, float)[order]
    ww = np.asarray(weights, float)[order]
    cdf = np.cumsum(ww)
    cdf /= cdf[-1]
    return float(np.interp(quantile, cdf, vv))


def response_summaries(q, f, mass_eV: float, z: float, sigma_kms: float, x, w):
    qres = resonant_q(np.abs(x), mass_eV, z, sigma_kms)
    F = np.interp(qres, q, f, left=f[0], right=f[-1])
    coherent = float(np.sum(w * np.abs(x) * F))
    rms = float(np.sqrt(np.sum(w * x * x * F * F)))
    q50 = weighted_quantile(qres, w, 0.50)
    q90 = weighted_quantile(qres, w, 0.90)
    return coherent, rms, q50, q90


def row_for(q, f0, fp, fm, mass_eV, z, sigma_kms, x, w):
    a0, r0, q50, q90 = response_summaries(q, f0, mass_eV, z, sigma_kms, x, w)
    ap, rp, _, _ = response_summaries(q, fp, mass_eV, z, sigma_kms, x, w)
    am, rm, _, _ = response_summaries(q, fm, mass_eV, z, sigma_kms, x, w)
    if a0 <= 0 or r0 <= 0:
        raise RuntimeError("Non-positive FD screening amplitude")
    mean_pair = abs(ap - am) / a0
    rms_pair = abs(rp - rm) / r0
    robust = min(mean_pair, rms_pair)
    return {
        "z": float(z),
        "sigma_v_kms": float(sigma_kms),
        "qres_median": q50,
        "qres_p90": q90,
        "coherent_plus_over_FD": ap / a0,
        "coherent_minus_over_FD": am / a0,
        "coherent_pair_fraction": mean_pair,
        "rms_plus_over_FD": rp / r0,
        "rms_minus_over_FD": rm / r0,
        "rms_pair_fraction": rms_pair,
        "robust_pair_fraction": robust,
        "FD_SN_required_for_hidden_1sigma": math.inf if robust == 0 else 1.0 / robust,
        "FD_SN_required_for_hidden_3sigma": math.inf if robust == 0 else 3.0 / robust,
    }


def write_csv(path: Path, rows):
    keys = list(rows[0].keys())
    with path.open("w", encoding="utf-8") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            vals = []
            for k in keys:
                v = row[k]
                vals.append("inf" if isinstance(v, float) and math.isinf(v) else f"{v:.12g}")
            f.write(",".join(vals) + "\n")


def parse_grid(text: str):
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one number")
    return vals


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="weak_lensing_screening")
    ap.add_argument("--pair-csv", type=Path, default=None,
                    help="Validated best_pair.csv from a wake run. If omitted, reconstruct CREF.")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--z-grid", type=parse_grid,
                    default=parse_grid("0.05,0.10,0.20,0.30,0.40,0.60,0.80"))
    ap.add_argument("--sigma-grid", type=parse_grid,
                    default=parse_grid("100,150,200,250,300,350,500,700,1000"),
                    help="One-component neutrino-CDM relative-velocity dispersion [km/s].")
    ap.add_argument("--quadrature", type=int, default=96)
    ap.add_argument("--baseline-fd-sn", type=float, default=None,
                    help="Optional fiducial FD wake-lensing detection S/N used only to rescale R_pair.")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.pair_csv is None:
        q, f0, fp, fm, shape, pair_stats = reconstruct_reference_pair(
            args.mass, args.z_match, args.frac
        )
        pair_source = "reconstructed_CREF"
    else:
        q, f0, fp, fm, shape = load_pair(args.pair_csv)
        pair_stats = None
        pair_source = str(args.pair_csv)

    np.savetxt(
        out / "pair_used.csv",
        np.column_stack([q, f0, fp, fm, shape]),
        delimiter=",",
        header="q,f_FD,f_plus,f_minus,half_pair_difference",
        comments="",
    )

    x, w = gh_normal(args.quadrature)
    rows = []
    for z in args.z_grid:
        for sigma in args.sigma_grid:
            rows.append(row_for(q, f0, fp, fm, args.mass, z, sigma, x, w))
    write_csv(out / "screening_grid.csv", rows)

    best = max(rows, key=lambda r: r["robust_pair_fraction"])
    # A representative at z=0.3 is useful because the published wake-lensing
    # geometry peaks near this redshift for source galaxies around z~0.8.
    z03 = [r for r in rows if abs(r["z"] - 0.30) < 1e-12]
    best_z03 = max(z03, key=lambda r: r["robust_pair_fraction"]) if z03 else None

    summary = {
        "calculation": "weak-lensing kinetic-response screening",
        "status": "screening_only_not_a_survey_forecast",
        "pair_source": pair_source,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "fractional_deformation_cap": args.frac,
        "quadrature_order": args.quadrature,
        "kinetic_kernel": "|u| F(m|u|/[Tnu0(1+z)])",
        "interpretation": (
            "R_pair is the hidden-state pair separation as a fraction of the fiducial FD "
            "wake-lensing amplitude. A full lensing calculation is warranted only if this "
            "fraction is large and stable over a physically relevant velocity range."
        ),
        "best_grid_point": best,
        "best_z_0p3_point": best_z03,
        "pair_stats_if_reconstructed": pair_stats,
    }
    if args.baseline_fd_sn is not None:
        summary["baseline_fd_sn"] = args.baseline_fd_sn
        summary["best_implied_hidden_sn"] = args.baseline_fd_sn * best["robust_pair_fraction"]
        if best_z03 is not None:
            summary["z0p3_implied_hidden_sn"] = (
                args.baseline_fd_sn * best_z03["robust_pair_fraction"]
            )

    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("WEAK_LENSING_SCREENING")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote: {out / 'screening_grid.csv'}")
    print(f"Wrote: {out / 'summary.json'}")
    print(f"Wrote: {out / 'pair_used.csv'}")


if __name__ == "__main__":
    main()
