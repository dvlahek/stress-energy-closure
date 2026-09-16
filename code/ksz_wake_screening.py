#!/usr/bin/env python3
r"""Screen kSZ velocity-tagging of the hidden-state neutrino wake response.

This is deliberately a screening calculation, not a survey forecast and not a
manuscript result.  The goal is to decide if a kSZ-tagged observable is worth a
full covariance calculation.

The nonrelativistic wake factor used in the paper is odd in the line-of-sight
relative velocity,

    W_F(x) ~ x F(q_res),
    q_res = m |x| sigma_v / [c T_nu,0 (1+z)].

A kSZ temperature/velocity tracer contributes another factor proportional to
x (up to optical-depth and reconstruction factors that are common to F_FD,
F_+, and F_- at this screening stage).  The simplest velocity-tagged cross
response is therefore

    K_F ~ < x^2 F(q_res) >.

For comparison we also evaluate the untagged coherent wake amplitude

    U_F ~ < |x| F(q_res) >,

and a power-like tagged statistic

    P_F ~ sqrt(< x^4 F(q_res)^2 >).

The key outputs are the pair separations relative to the fiducial FD response,

    R_ksz = |K_+ - K_-| / K_FD,
    R_wake = |U_+ - U_-| / U_FD,
    R_power = |P_+ - P_-| / P_FD.

The ratio R_ksz/R_wake tells us if adding a velocity tag enhances or suppresses
hidden-state information at fixed underlying wake amplitude.  No absolute kSZ
S/N is inferred here.  A full calculation would need optical-depth scatter,
velocity-reconstruction noise, halo/galaxy selection, survey geometry, and the
joint covariance with the wake observable.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import class_response_optimize as cro

C_KMS = 299792.458
TNU0_EV = cro.T_NCDM * cro.TCMB_K * cro.KB_EV_K

# Stable leading hidden-state direction used in the validated wake forecast.
CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def reconstruct_pair(mass_eV: float, z_match: float, frac: float):
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, _, _, shapes, _, _ = cro.kinetic_objects(q, mass_eV, z_match)
    if shapes.shape[0] != CREF.size:
        raise RuntimeError(
            f"CREF length {CREF.size} != null-space dimension {shapes.shape[0]}"
        )
    coeff = CREF / np.linalg.norm(CREF)
    raw = coeff @ shapes
    maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
    if not np.isfinite(maxrel) or maxrel <= 0:
        raise RuntimeError("invalid hidden-state direction")
    shape = raw / maxrel
    fp = f0 + frac * shape
    fm = f0 - frac * shape
    if np.min(fp) <= 0 or np.min(fm) <= 0:
        raise RuntimeError("pair violates positivity")
    return q, f0, fp, fm, shape, cro.pair_stats(f0, fp, fm, q, weights)


def load_pair(path: Path):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    required = {"q", "f_FD", "f_plus", "f_minus"}
    if arr.dtype.names is None or not required.issubset(set(arr.dtype.names)):
        raise RuntimeError(
            f"{path} must contain q,f_FD,f_plus,f_minus; found {arr.dtype.names}"
        )
    q = np.asarray(arr["q"], float)
    f0 = np.asarray(arr["f_FD"], float)
    fp = np.asarray(arr["f_plus"], float)
    fm = np.asarray(arr["f_minus"], float)
    if np.any(np.diff(q) <= 0):
        raise RuntimeError("q grid must be strictly increasing")
    if min(np.min(f0), np.min(fp), np.min(fm)) <= 0:
        raise RuntimeError("pair contains non-positive values")
    return q, f0, fp, fm, 0.5 * (fp - fm), None


def gh_normal(order: int):
    x0, w0 = np.polynomial.hermite.hermgauss(order)
    return np.sqrt(2.0) * x0, w0 / np.sqrt(np.pi)


def parse_grid(text: str):
    vals = [float(v.strip()) for v in text.split(",") if v.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return vals


def q_resonant(abs_x, mass_eV: float, z: float, sigma_kms: float):
    return mass_eV * abs_x * sigma_kms / (C_KMS * TNU0_EV * (1.0 + z))


def responses(q, f, mass_eV: float, z: float, sigma_kms: float, x, w):
    qres = q_resonant(np.abs(x), mass_eV, z, sigma_kms)
    F = np.interp(qres, q, f, left=f[0], right=f[-1])

    # Odd wake amplitude after coherent orientation.
    wake = float(np.sum(w * np.abs(x) * F))

    # kSZ-like velocity tag adds one power of LOS velocity.
    ksz_tagged = float(np.sum(w * x * x * F))

    # Power-like tagged cross-check; useful because a real estimator may weight
    # high-velocity objects differently from the coherent mean.
    tagged_power = float(np.sqrt(np.sum(w * x**4 * F**2)))

    # Dimensionless resonant-momentum diagnostics.
    qmean = float(np.sum(w * qres))
    order = np.argsort(qres)
    qr = qres[order]
    ww = w[order]
    cdf = np.cumsum(ww) / np.sum(ww)
    q90 = float(np.interp(0.90, cdf, qr))
    return wake, ksz_tagged, tagged_power, qmean, q90


def row_for(q, f0, fp, fm, mass_eV, z, sigma_kms, x, w):
    u0, k0, p0, qmean, q90 = responses(q, f0, mass_eV, z, sigma_kms, x, w)
    up, kp, pp, _, _ = responses(q, fp, mass_eV, z, sigma_kms, x, w)
    um, km, pm, _, _ = responses(q, fm, mass_eV, z, sigma_kms, x, w)
    if min(u0, k0, p0) <= 0:
        raise RuntimeError("non-positive fiducial screening response")

    r_wake = abs(up - um) / u0
    r_ksz = abs(kp - km) / k0
    r_power = abs(pp - pm) / p0
    gain = r_ksz / r_wake if r_wake > 0 else math.inf
    robust = min(r_ksz, r_power)

    return {
        "z": float(z),
        "sigma_v_kms": float(sigma_kms),
        "qres_mean": qmean,
        "qres_p90": q90,
        "wake_plus_over_FD": up / u0,
        "wake_minus_over_FD": um / u0,
        "wake_pair_fraction": r_wake,
        "ksz_plus_over_FD": kp / k0,
        "ksz_minus_over_FD": km / k0,
        "ksz_pair_fraction": r_ksz,
        "tagged_power_pair_fraction": r_power,
        "ksz_information_gain_over_wake": gain,
        "robust_tagged_pair_fraction": robust,
        "fiducial_tagged_SN_required_for_hidden_1sigma": (
            math.inf if robust == 0 else 1.0 / robust
        ),
        "fiducial_tagged_SN_required_for_hidden_3sigma": (
            math.inf if robust == 0 else 3.0 / robust
        ),
    }


def write_csv(path: Path, rows):
    keys = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="ksz_wake_screening")
    ap.add_argument("--pair-csv", type=Path, default=None,
                    help="Validated wake best_pair.csv; otherwise reconstruct CREF.")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--z-grid", type=parse_grid,
                    default=parse_grid("0.05,0.10,0.20,0.30,0.40,0.60,0.80,1.0"))
    ap.add_argument("--sigma-grid", type=parse_grid,
                    default=parse_grid("100,150,200,250,300,350,500,700,1000"))
    ap.add_argument("--quadrature", type=int, default=128)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.pair_csv is None:
        q, f0, fp, fm, shape, pair_stats = reconstruct_pair(
            args.mass, args.z_match, args.frac
        )
        pair_source = "reconstructed_CREF"
    else:
        q, f0, fp, fm, shape, pair_stats = load_pair(args.pair_csv)
        pair_source = str(args.pair_csv)

    np.savetxt(
        out / "pair_used.csv",
        np.column_stack([q, f0, fp, fm, shape]),
        delimiter=",",
        header="q,f_FD,f_plus,f_minus,half_pair_difference",
        comments="",
    )

    x, w = gh_normal(args.quadrature)
    rows = [
        row_for(q, f0, fp, fm, args.mass, z, sigma, x, w)
        for z in args.z_grid
        for sigma in args.sigma_grid
    ]
    write_csv(out / "screening_grid.csv", rows)

    best_robust = max(rows, key=lambda r: r["robust_tagged_pair_fraction"])
    best_gain = max(rows, key=lambda r: r["ksz_information_gain_over_wake"])
    z03 = [r for r in rows if abs(r["z"] - 0.30) < 1e-12]
    best_z03 = max(z03, key=lambda r: r["robust_tagged_pair_fraction"]) if z03 else None

    summary = {
        "calculation": "kSZ-tagged neutrino-wake kinetic screening",
        "status": "screening_only_not_a_survey_forecast",
        "pair_source": pair_source,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "fractional_deformation_cap": args.frac,
        "quadrature_order": args.quadrature,
        "wake_kernel": "x F(m|x| sigma_v/[c Tnu0(1+z)])",
        "ksz_tagged_kernel": "x^2 F(m|x| sigma_v/[c Tnu0(1+z)])",
        "interpretation": (
            "The kSZ tag contributes one LOS-velocity factor to the odd wake. "
            "R_ksz/R_wake tests only fractional hidden-state information retention. "
            "No absolute kSZ-wake detection significance is inferred."
        ),
        "best_robust_tagged_point": best_robust,
        "best_information_gain_point": best_gain,
        "best_z_0p3_point": best_z03,
        "pair_stats_if_reconstructed": pair_stats,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("KSZ_WAKE_SCREENING")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out / 'summary.json'}")
    print(f"Wrote {out / 'screening_grid.csv'}")
    print(f"Wrote {out / 'pair_used.csv'}")


if __name__ == "__main__":
    main()
