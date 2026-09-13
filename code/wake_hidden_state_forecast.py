#!/usr/bin/env python3
"""Matched-state neutrino-wake detectability diagnostic.

This script optimizes the exact n-rho-P null space for the parity-odd neutrino-wake
response. For an isotropic collisionless distribution the Landau-pole contribution
samples the occupation at q_res = m v_parallel/[T_nu0 (1+z)]. We therefore optimize
matched F_+/F_- pairs for their RMS fractional wake separation over a low-z velocity
window and translate that separation into the galaxy-count scaling of Okoli et al.
(MNRAS 468, 2164, 2017, Eq. 35). This is a detectability diagnostic, not a modern
survey likelihood.
"""
from __future__ import annotations

import argparse, json
from pathlib import Path
import numpy as np
import class_response_optimize as cro

C_KMS = 299792.458
TNU0_EV = cro.T_NCDM * cro.TCMB_K * cro.KB_EV_K


def interp(x, xp, fp):
    return np.interp(x, xp, fp, left=fp[0], right=fp[-1])


def make_samples(mass, zgrid, vgrid):
    zz, vv = np.meshgrid(zgrid, vgrid, indexing="ij")
    qres = mass * (vv / C_KMS) / (TNU0_EV * (1.0 + zz))
    return zz.ravel(), vv.ravel(), qres.ravel()


def optimize_for_mass(mass, z_match, frac, samples, seed, zgrid, vgrid):
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, N, shapes, y, M = cro.kinetic_objects(q, mass, z_match)
    z_s, v_s, qres = make_samples(mass, zgrid, vgrid)
    fbase = interp(qres, q, f0)

    # Null directions are normalized so max |shape|/f0 = 1. The +/- pair at cap
    # frac has Delta F/F0 = 2 frac shape/f0.
    R = np.array([2.0 * frac * interp(qres, q, s) / np.maximum(fbase, 1e-300) for s in shapes])
    B = R @ R.T / R.shape[1]
    evals, evecs = np.linalg.eigh(B)
    candidates = [evecs[:, np.argmax(evals)]]
    nd = shapes.shape[0]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0; candidates.extend([e, -e])
    rng = np.random.default_rng(seed)
    rr = rng.normal(size=(samples, nd)); rr /= np.linalg.norm(rr, axis=1, keepdims=True)
    candidates.extend(rr)

    best = None
    for c in candidates:
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0: continue
        shape = raw / maxrel
        rel = 2.0 * frac * interp(qres, q, shape) / np.maximum(fbase, 1e-300)
        score = float(np.sqrt(np.mean(rel**2)))
        if best is None or score > best[0]:
            best = (score, shape.copy(), rel.copy(), c.copy()/maxrel)

    score, shape, rel, coeff = best
    fp = f0 + frac * shape
    fm = f0 - frac * shape
    mp, mm = cro.moments(fp, q, weights), cro.moments(fm, q, weights)
    mismatch = np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)

    by_z = []
    for z in zgrid:
        sel = np.isclose(z_s, z)
        r = rel[sel]
        # Okoli+17 Eq.35: galaxy count for >3 sigma baseline wake with N_nu=1, Delta b=1.
        nbase = 1.7e7 * (mass/0.05)**(-6.0) * 28.5**float(z)
        reff = float(np.sqrt(np.mean(r**2)))
        nhidden = float(nbase / max(reff**2, 1e-30))
        by_z.append({
            "z": float(z),
            "wake_fractional_RMS": reff,
            "wake_fractional_median_abs": float(np.median(np.abs(r))),
            "wake_fractional_max_abs": float(np.max(np.abs(r))),
            "Okoli17_baseline_Ngal_3sigma": float(nbase),
            "scaled_hidden_state_Ngal_3sigma": nhidden,
        })

    return {
        "mass_eV": mass,
        "mass_over_Tnu_at_match": float(y),
        "max_relative_moment_mismatch": float(mismatch.max()),
        "global_wake_fractional_RMS": float(np.sqrt(np.mean(rel**2))),
        "global_wake_fractional_median_abs": float(np.median(np.abs(rel))),
        "global_wake_fractional_max_abs": float(np.max(np.abs(rel))),
        "qres_min": float(qres.min()), "qres_max": float(qres.max()),
        "by_redshift": by_z,
        "q": q.tolist(), "f0": f0.tolist(), "fp": fp.tolist(), "fm": fm.tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="wake_hidden_state_output")
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    masses = [0.06,0.10,0.20,0.60]
    zgrid = np.array([0.0,0.25,0.50])
    vgrid = np.linspace(100.0,500.0,81)
    results=[]
    for i,m in enumerate(masses):
        r=optimize_for_mass(m,args.z_match,args.frac,args.samples,args.seed+i,zgrid,vgrid)
        q=np.asarray(r.pop("q")); f0=np.asarray(r.pop("f0")); fp=np.asarray(r.pop("fp")); fm=np.asarray(r.pop("fm"))
        np.savetxt(out/f"pair_m{m:.2f}.csv",np.column_stack([q,f0,fp,fm]),delimiter=",",header="q,f0,fplus,fminus",comments="")
        results.append(r)
        print(json.dumps(r,indent=2))
    summary={
        "class_commit": cro.CLASS_COMMIT,
        "Tnu0_eV": TNU0_EV,
        "fractional_distortion_cap": args.frac,
        "velocity_window_km_s": [100.0,500.0],
        "redshift_grid": zgrid.tolist(),
        "mass_results": results,
        "interpretation": "Parity-odd wake response diagnostic using the exact resonant occupation and Okoli et al. 2017 Eq.35 galaxy-count scaling. Not a survey likelihood and not a replacement for a full relative-velocity covariance calculation."
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")

if __name__ == "__main__": main()
