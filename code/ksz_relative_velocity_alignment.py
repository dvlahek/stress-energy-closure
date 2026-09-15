#!/usr/bin/env python3
"""Screen if kSZ-traced baryon velocities align with the neutrino-CDM relative velocity.

This is a cheap prerequisite for any kSZ-tagged neutrino-wake forecast. The wake
is controlled by the neutrino-CDM relative velocity, while kSZ traces the baryon
electron bulk velocity. If the two large-scale velocity fields are poorly aligned,
a kSZ tag cannot recover the wake direction efficiently even if the kinetic
response itself has large hidden-state contrast.

For each redshift and smoothing scale R, the script computes signed linear-theory
velocity transfer amplitudes from CLASS using finite redshift derivatives of the
density transfer functions. It reports the filtered field correlation

    r(rel,b) = <v_rel v_b> / sqrt(<v_rel^2><v_b^2>),

with v_rel = v_ncdm - v_cdm, together with the analogous CDM correlation and an
"effective" alignment after multiplying by representative velocity-reconstruction
correlations r_rec. This is still a screening layer, not a survey forecast.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from classy import Class

import class_response_optimize as cro

C_KMS = 299792.458
NS = 0.9649
KPIV_MPC = 0.05
KGRID = np.geomspace(1.0e-4, 0.30, 360)  # h/Mpc
CLASS_COMMIT_EXPECTED = cro.CLASS_COMMIT


def parse_grid(text: str):
    vals = [float(v.strip()) for v in text.split(",") if v.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return vals


def params(mass_eV: float):
    return {
        "output": "mPk,dTk",
        "modes": "s",
        "gauge": "newtonian",
        "H0": cro.H0,
        "omega_b": cro.OMEGA_B,
        "omega_cdm": cro.OMEGA_CDM,
        "A_s": cro.A_S,
        "n_s": NS,
        "tau_reio": cro.TAU_REIO,
        "N_ur": cro.N_UR,
        "N_ncdm": 1,
        "m_ncdm": mass_eV,
        "T_ncdm": cro.T_NCDM,
        "deg_ncdm": 1.0,
        "P_k_max_h/Mpc": 0.40,
        "z_max_pk": 1.5,
    }


def transfer(cosmo, z: float, key: str):
    t = cosmo.get_transfer(z=float(z), output_format="class")
    kin = np.asarray(t["k (h/Mpc)"], float)
    val = np.asarray(t[key], float)
    return np.interp(KGRID, kin, val)


def signed_velocity_from_density(cosmo, z: float, key: str):
    """Return a signed velocity-amplitude proxy with the same convention for all species.

    From continuity, a*dot(delta)/k = -H d(delta)/dz / k because a(1+z)=1.
    The common normalization cancels in the field-correlation coefficient.
    """
    dz = 0.003 * (1.0 + z)
    zm = max(0.0, z - dz)
    zp = z + dz
    if zm == z:
        d0 = transfer(cosmo, z, key)
        dp = transfer(cosmo, zp, key)
        d_dz = (dp - d0) / (zp - z)
    else:
        dm = transfer(cosmo, zm, key)
        dp = transfer(cosmo, zp, key)
        d_dz = (dp - dm) / (zp - zm)

    h = float(cosmo.h())
    k_mpc = KGRID * h
    H_mpc = float(cosmo.Hubble(float(z)))
    primordial = np.sqrt(cro.A_S * (k_mpc / KPIV_MPC) ** (NS - 1.0))
    return -H_mpc * d_dz / np.maximum(k_mpc, 1e-30) * primordial * C_KMS


def tophat(x):
    x = np.asarray(x, float)
    out = np.ones_like(x)
    m = np.abs(x) > 1e-5
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / xm**3
    return out


def integrate_logk(y):
    lk = np.log(KGRID)
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, lk))
    return float(np.trapz(y, lk))


def filtered_stats(vrel, vtag, R):
    w = tophat(KGRID * R)
    aa = integrate_logk(vrel * vrel * w * w)
    bb = integrate_logk(vtag * vtag * w * w)
    ab = integrate_logk(vrel * vtag * w * w)
    aa = max(aa, 0.0)
    bb = max(bb, 0.0)
    den = math.sqrt(max(aa * bb, 1e-300))
    r = ab / den
    slope = ab / max(bb, 1e-300)
    return {
        "r": float(np.clip(r, -1.0, 1.0)),
        "sigma_rel_kms": math.sqrt(aa),
        "sigma_tag_kms": math.sqrt(bb),
        "sigma_rel_over_tag": math.sqrt(aa / max(bb, 1e-300)),
        "regression_rel_on_tag": slope,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="ksz_relative_velocity_alignment")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-grid", type=parse_grid,
                    default=parse_grid("0.05,0.10,0.20,0.30,0.40,0.60,0.80,1.0"))
    ap.add_argument("--R-grid", type=parse_grid,
                    default=parse_grid("5,8,12,16,25,40,60"),
                    help="Top-hat smoothing radii in Mpc/h.")
    ap.add_argument("--rrec-grid", type=parse_grid,
                    default=parse_grid("0.50,0.65,0.70,0.80,1.0"),
                    help="Representative correlation of reconstructed and true baryon velocity.")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    cosmo = Class()
    cosmo.set(params(args.mass))
    cosmo.compute()
    rows = []
    try:
        for z in args.z_grid:
            vb = signed_velocity_from_density(cosmo, z, "d_b")
            vc = signed_velocity_from_density(cosmo, z, "d_cdm")
            vn = signed_velocity_from_density(cosmo, z, "d_ncdm[0]")
            vrel = vn - vc

            for R in args.R_grid:
                sb = filtered_stats(vrel, vb, R)
                sc = filtered_stats(vrel, vc, R)
                row = {
                    "z": float(z),
                    "R_Mpc_h": float(R),
                    "r_rel_baryon": sb["r"],
                    "r_rel_cdm": sc["r"],
                    "abs_r_rel_baryon": abs(sb["r"]),
                    "sigma_rel_kms": sb["sigma_rel_kms"],
                    "sigma_baryon_kms": sb["sigma_tag_kms"],
                    "sigma_rel_over_baryon": sb["sigma_rel_over_tag"],
                    "regression_rel_on_baryon": sb["regression_rel_on_tag"],
                }
                for rrec in args.rrec_grid:
                    row[f"effective_abs_alignment_rrec_{rrec:.2f}"] = abs(sb["r"]) * rrec
                rows.append(row)
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()

    keys = list(rows[0].keys())
    with (out / "alignment_grid.csv").open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)

    best = max(rows, key=lambda r: r["abs_r_rel_baryon"])
    worst = min(rows, key=lambda r: r["abs_r_rel_baryon"])
    z03 = [r for r in rows if abs(r["z"] - 0.30) < 1e-12]
    best_z03 = max(z03, key=lambda r: r["abs_r_rel_baryon"]) if z03 else None

    arr = np.array([r["abs_r_rel_baryon"] for r in rows], float)
    summary = {
        "calculation": "kSZ-baryon / neutrino-CDM relative-velocity alignment screening",
        "status": "screening_only_not_a_survey_forecast",
        "mass_eV": args.mass,
        "class_commit_expected": CLASS_COMMIT_EXPECTED,
        "k_range_h_Mpc": [float(KGRID[0]), float(KGRID[-1])],
        "definition": "r=<v_rel v_b>/sqrt(<v_rel^2><v_b^2>) after top-hat smoothing",
        "best_alignment": best,
        "worst_alignment": worst,
        "best_z_0p3_alignment": best_z03,
        "median_abs_alignment": float(np.median(arr)),
        "min_abs_alignment": float(np.min(arr)),
        "max_abs_alignment": float(np.max(arr)),
        "decision_guide": {
            "strong": "|r| >= 0.7 over relevant R,z: kSZ velocity tag is physically well aligned; continue to a covariance forecast.",
            "intermediate": "0.4 <= |r| < 0.7: possible but reconstruction noise will matter strongly.",
            "weak": "|r| < 0.4 over relevant R,z: stop; kSZ traces the wrong velocity direction for this wake test."
        }
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("KSZ_RELATIVE_VELOCITY_ALIGNMENT")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out / 'summary.json'}")
    print(f"Wrote {out / 'alignment_grid.csv'}")


if __name__ == "__main__":
    main()
