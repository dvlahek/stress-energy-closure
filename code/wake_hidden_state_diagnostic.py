#!/usr/bin/env python3
"""Parity-odd neutrino-wake diagnostic for stress-energy-matched relic states.

For an isotropic collisionless distribution F(p) viewed in the rest frame of a
halo with neutrino-CDM relative velocity u, the Landau-pole contribution to the
steady linear Vlasov wake has an odd density response proportional to
u_parallel * F(m |u_parallel|), up to common known factors. The wake therefore
samples the distribution at a resonant momentum instead of only low moments.
This script quantifies the response-shape separation for a matched pair.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

T_NCDM = 0.71611
TCMB_K = 2.7255
KB_EV_K = 8.617333262e-5
C_KMS = 299792.458


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pair", required=True)
    ap.add_argument("--mass", type=float, default=0.60)
    ap.add_argument("--qmin", type=float, default=0.5)
    ap.add_argument("--qmax", type=float, default=6.0)
    ap.add_argument("--out", default="wake_hidden_state_summary.json")
    args = ap.parse_args()

    data = np.genfromtxt(args.pair, delimiter=",", names=True)
    q = data["q"]
    fp = data["f_plus"]
    fm = data["f_minus"]
    mask = (q >= args.qmin) & (q <= args.qmax)
    if not np.any(mask):
        raise RuntimeError("Empty resonant-momentum window")

    rel = 2.0 * (fp - fm) / np.maximum(fp + fm, 1e-300)
    rr = rel[mask]
    qq = q[mask]
    ii = int(np.argmax(np.abs(rr)))
    Tnu0_eV = T_NCDM * TCMB_K * KB_EV_K
    vel = qq * Tnu0_eV / args.mass * C_KMS

    summary = {
        "mass_eV": args.mass,
        "q_window": [args.qmin, args.qmax],
        "velocity_window_km_s": [float(vel.min()), float(vel.max())],
        "rms_symmetric_relative_wake_difference": float(np.sqrt(np.mean(rr**2))),
        "mean_abs_symmetric_relative_wake_difference": float(np.mean(np.abs(rr))),
        "median_abs_symmetric_relative_wake_difference": float(np.median(np.abs(rr))),
        "max_abs_symmetric_relative_wake_difference": float(np.max(np.abs(rr))),
        "q_at_max": float(qq[ii]),
        "velocity_km_s_at_max": float(vel[ii]),
        "signed_relative_difference_at_max": float(rr[ii]),
        "interpretation": "Odd Landau-pole wake response shape only, not a survey S/N forecast."
    }
    Path(args.out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
