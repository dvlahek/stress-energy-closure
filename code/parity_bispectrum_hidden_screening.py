#!/usr/bin/env python3
r"""Screen hidden-state sensitivity of a parity-odd neutrino relative-velocity channel.

This is a *screening* calculation, not a bispectrum forecast and not a manuscript
result.  The motivation is the parity-odd relative-velocity observable proposed
by Zhu & Castorina (Phys. Rev. D 101, 023525), whose signal is linear in the
CDM-neutrino relative-velocity transfer at leading order.

For the stress-energy-matched hidden states F_+ and F_- used in this repository,
we recompute the linear relative-velocity transfer with CLASS and ask how strongly
it changes at fixed neutrino mass.  We also report a simple density-weighted proxy
v_rel P_cb to check that the conclusion is not an artefact of looking at velocity
alone.  No absolute bispectrum S/N is inferred here.
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
CLASS_COMMIT = cro.CLASS_COMMIT

CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def parse_grid(text: str):
    vals = [float(v.strip()) for v in text.split(",") if v.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return vals


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
    stats = cro.pair_stats(f0, fp, fm, q, weights)
    return q, f0, fp, fm, shape, stats


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


def write_psd(path: Path, q, f):
    np.savetxt(path, np.column_stack([q, f]), fmt="%.14e")


def class_params(psd: Path, mass: float, zmax: float, kmax_h: float):
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
        "use_ncdm_psd_files": 1,
        "ncdm_psd_filenames": str(psd.resolve()),
        "m_ncdm": mass,
        "T_ncdm": cro.T_NCDM,
        "deg_ncdm": 1.0,
        "P_k_max_h/Mpc": max(0.35, 1.1 * kmax_h),
        "z_max_pk": max(1.0, zmax + 0.2),
    }


def interp_transfer(cosmo, z, key, kh):
    t = cosmo.get_transfer(z=float(z), output_format="class")
    kin = np.asarray(t["k (h/Mpc)"], float)
    val = np.asarray(t[key], float)
    return np.interp(kh, kin, val)


def delta_rel(cosmo, z, kh):
    return interp_transfer(cosmo, z, "d_ncdm[0]", kh) - interp_transfer(cosmo, z, "d_cdm", kh)


def relative_velocity_kms(cosmo, z, kh):
    dz = 0.004 * (1.0 + z)
    zm = max(0.0, z - dz)
    zp = z + dz
    if zm == 0.0 and z < dz:
        d0 = delta_rel(cosmo, z, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - d0) / (zp - z)
    else:
        dm = delta_rel(cosmo, zm, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - dm) / (zp - zm)
    h = float(cosmo.h())
    k_mpc = kh * h
    H_mpc = float(cosmo.Hubble(float(z)))
    primordial = np.sqrt(cro.A_S * (k_mpc / KPIV_MPC) ** (NS - 1.0))
    return -H_mpc * d_dz / np.maximum(k_mpc, 1e-30) * primordial * C_KMS


def pk_cb_h(cosmo, z, kh):
    h = float(cosmo.h())
    return np.array([cosmo.pk_cb_lin(float(k * h), float(z)) * h**3 for k in kh])


def build_state(psd: Path, mass: float, zgrid, kh):
    c = Class()
    c.set(class_params(psd, mass, max(zgrid), max(kh)))
    c.compute()
    try:
        out = {}
        for z in zgrid:
            out[float(z)] = {
                "vrel": relative_velocity_kms(c, float(z), kh),
                "pk": pk_cb_h(c, float(z), kh),
            }
    finally:
        c.struct_cleanup()
        c.empty()
    return out


def weighted_rms(frac, weights):
    weights = np.asarray(weights, float)
    weights = np.maximum(weights, 0.0)
    if not np.any(weights > 0):
        return float("nan")
    return float(np.sqrt(np.sum(weights * frac * frac) / np.sum(weights)))


def row_for_z(z, kh, s0, sp, sm):
    v0, vp, vm = s0["vrel"], sp["vrel"], sm["vrel"]
    p0, pp, pm = s0["pk"], sp["pk"], sm["pk"]

    vmax = max(float(np.max(np.abs(v0))), 1e-300)
    mask = np.abs(v0) > 1e-4 * vmax
    if np.sum(mask) < 5:
        raise RuntimeError(f"too few stable relative-velocity modes at z={z}")

    # Symmetric hidden-state response around FD.  If the parity-odd signal is
    # linear in v_rel, this is the leading fractional pair response.
    f_v = np.zeros_like(v0)
    f_v[mask] = np.abs(vp[mask] - vm[mask]) / (2.0 * np.abs(v0[mask]))

    # Simple density-weighted proxy, not the full bispectrum kernel.
    a0, ap, am = v0 * p0, vp * pp, vm * pm
    amax = max(float(np.max(np.abs(a0))), 1e-300)
    mask_a = np.abs(a0) > 1e-4 * amax
    f_a = np.zeros_like(a0)
    f_a[mask_a] = np.abs(ap[mask_a] - am[mask_a]) / (2.0 * np.abs(a0[mask_a]))

    lnk = np.log(kh)
    dlnk = np.gradient(lnk)
    mode_w = kh**3 * dlnk
    power_w = kh**3 * p0 * dlnk

    mv = f_v[mask]
    ma = f_a[mask_a]
    return {
        "z": float(z),
        "velocity_fraction_median": float(np.median(mv)),
        "velocity_fraction_p90": float(np.quantile(mv, 0.90)),
        "velocity_fraction_max": float(np.max(mv)),
        "velocity_fraction_rms_logk": weighted_rms(f_v[mask], dlnk[mask]),
        "velocity_fraction_rms_modecount": weighted_rms(f_v[mask], mode_w[mask]),
        "velocity_fraction_rms_powerweighted": weighted_rms(f_v[mask], power_w[mask]),
        "proxy_fraction_median": float(np.median(ma)),
        "proxy_fraction_p90": float(np.quantile(ma, 0.90)),
        "proxy_fraction_rms_modecount": weighted_rms(f_a[mask_a], mode_w[mask_a]),
        "pair_difference_sign_consistency": float(
            abs(np.mean(np.sign((vp - vm)[mask])))
        ),
        "vrel_FD_rms_kms": float(np.sqrt(np.mean(v0[mask] ** 2))),
    }


def write_csv(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="parity_bispectrum_screening")
    ap.add_argument("--pair-csv", type=Path, default=None,
                    help="Validated wake best_pair.csv; otherwise reconstruct CREF.")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--z-grid", type=parse_grid,
                    default=parse_grid("0.05,0.10,0.20,0.30,0.40,0.60,0.80,1.0,1.5,2.0"))
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=96)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.pair_csv is None:
        q, f0, fp, fm, shape, pair_stats = reconstruct_pair(args.mass, args.z_match, args.frac)
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
    p0 = out / "fd.dat"
    pp = out / "plus.dat"
    pm = out / "minus.dat"
    write_psd(p0, q, f0)
    write_psd(pp, q, fp)
    write_psd(pm, q, fm)

    kh = np.geomspace(args.kmin, args.kmax, args.nk)
    s0 = build_state(p0, args.mass, args.z_grid, kh)
    sp = build_state(pp, args.mass, args.z_grid, kh)
    sm = build_state(pm, args.mass, args.z_grid, kh)

    rows = [row_for_z(float(z), kh, s0[float(z)], sp[float(z)], sm[float(z)])
            for z in args.z_grid]
    write_csv(out / "screening_by_redshift.csv", rows)

    best = max(rows, key=lambda r: r["proxy_fraction_rms_modecount"])
    z03 = min(rows, key=lambda r: abs(r["z"] - 0.30))
    med_proxy = float(np.median([r["proxy_fraction_rms_modecount"] for r in rows]))
    med_vel = float(np.median([r["velocity_fraction_rms_modecount"] for r in rows]))

    summary = {
        "calculation": "parity-odd relative-velocity hidden-state screening",
        "status": "screening_only_not_a_bispectrum_forecast",
        "literature_target": "Zhu & Castorina, Phys. Rev. D 101, 023525 (2020)",
        "class_commit_expected": CLASS_COMMIT,
        "pair_source": pair_source,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "fractional_deformation_cap": args.frac,
        "k_range_h_Mpc": [args.kmin, args.kmax],
        "nk": args.nk,
        "definition_velocity_fraction": "|vrel_plus-vrel_minus|/(2|vrel_FD|)",
        "definition_proxy_fraction": "|vrel_plus Pcb_plus-vrel_minus Pcb_minus|/(2|vrel_FD Pcb_FD|)",
        "median_velocity_fraction_rms_modecount_across_z": med_vel,
        "median_proxy_fraction_rms_modecount_across_z": med_proxy,
        "best_redshift": best,
        "z_0p3": z03,
        "pair_stats_if_reconstructed": pair_stats,
        "decision_guide": {
            "strong": "mode-count RMS hidden fraction >=0.10 over useful z: build the full parity-odd bispectrum Fisher forecast.",
            "intermediate": "0.03-0.10: only continue if the fiducial literature forecast has very high S/N.",
            "weak": "<0.03: hidden-state direction is too suppressed in this channel."
        },
        "analysis_scope": "We evaluate a kernel-level bispectrum response; survey sensitivity is not estimated."
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("PARITY_BISPECTRUM_HIDDEN_SCREENING")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out / 'summary.json'}")
    print(f"Wrote {out / 'screening_by_redshift.csv'}")
    print(f"Wrote {out / 'pair_used.csv'}")


if __name__ == "__main__":
    main()
