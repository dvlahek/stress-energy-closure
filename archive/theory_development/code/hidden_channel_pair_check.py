#!/usr/bin/env python3
r"""Checkpointed low-memory validation for one hidden-state pair.

This is the decisive, OOM-safe follow-up to hidden_channel_forensics.py.  Each
invocation evaluates only one redshift, one hidden direction and one precision
setting.  CLASS itself runs in short-lived subprocesses through
hidden_channel_forensics_lowmem.isolated_build_state, so a failed run does not
lose previous checkpoints and the long operator-probe campaign is not repeated.

Typical use:
  python code/hidden_channel_pair_check.py --direction cref --z 0.3
  python code/hidden_channel_pair_check.py --direction selected --z 0.3 \
      --pair-csv hidden_channel_operator_diagnostics/best_proxy_pair.csv
  python code/hidden_channel_pair_check.py --direction cref --z 0.3 --high-precision

The reported linear quantity is the same half-pair fractional density-proxy
response used in the forensics scripts.  The wake quantity uses the same
half-pair convention, so their ratio is directly comparable.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import class_response_optimize as cro
import hidden_channel_operator_diagnostics as hd
import hidden_channel_forensics as hf
from hidden_channel_forensics_lowmem import isolated_build_state

CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def reconstruct_cref(mass: float, z_match: float, frac: float):
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, _, _, shapes, _, _ = cro.kinetic_objects(q, mass, z_match)
    c = CREF / np.linalg.norm(CREF)
    rel_shapes = shapes / np.maximum(f0[None, :], 1e-300)
    nrm = hd.normalization_for_coeff(c, rel_shapes)
    shape = nrm * (c @ shapes)
    fp = f0 + frac * shape
    fm = f0 - frac * shape
    if min(float(np.min(fp)), float(np.min(fm))) <= 0:
        raise RuntimeError("CREF pair violates positivity")
    return q, f0, fp, fm, weights


def load_selected(path: Path, mass: float, z_match: float):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    required = {"q", "f_FD", "f_plus", "f_minus"}
    if arr.dtype.names is None or not required.issubset(set(arr.dtype.names)):
        raise RuntimeError(f"{path} must contain q,f_FD,f_plus,f_minus")
    q = np.asarray(arr["q"], float)
    f0 = np.asarray(arr["f_FD"], float)
    fp = np.asarray(arr["f_plus"], float)
    fm = np.asarray(arr["f_minus"], float)
    # Moment weights on exactly the same q grid.
    _, weights, _, _, _, _, _ = cro.kinetic_objects(q, mass, z_match)
    return q, f0, fp, fm, weights


def write_psd(path: Path, q, f):
    hd.write_psd(path, q, f)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="hidden_channel_pair_checks")
    ap.add_argument("--direction", choices=("cref", "selected"), required=True)
    ap.add_argument("--pair-csv", type=Path, default=Path("hidden_channel_operator_diagnostics/best_proxy_pair.csv"))
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30,
                    help="CREF cap. Selected pair is read directly from --pair-csv.")
    ap.add_argument("--z", type=float, default=0.30)
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=96)
    ap.add_argument("--high-precision", action="store_true")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.direction == "cref":
        q, f0, fp, fm, weights = reconstruct_cref(args.mass, args.z_match, args.frac)
    else:
        if not args.pair_csv.exists():
            raise FileNotFoundError(
                f"Selected-pair file not found: {args.pair_csv}. "
                "Use the best_proxy_pair.csv produced by hidden_channel_operator_diagnostics.py"
            )
        q, f0, fp, fm, weights = load_selected(args.pair_csv, args.mass, args.z_match)

    tag = f"{args.direction}_z{args.z:.3f}_{'hp' if args.high_precision else 'std'}"
    p0 = out / f"{tag}_fd.dat"
    pp = out / f"{tag}_plus.dat"
    pm = out / f"{tag}_minus.dat"
    write_psd(p0, q, f0)
    write_psd(pp, q, fp)
    write_psd(pm, q, fm)

    kh = np.geomspace(args.kmin, args.kmax, args.nk)
    zgrid = [float(args.z)]

    print(f"STATE 1/3 FD  direction={args.direction} z={args.z} hp={args.high_precision}", flush=True)
    s0, theta0, keys = isolated_build_state(p0, args.mass, zgrid, kh, high_precision=args.high_precision)
    print("STATE 2/3 PLUS", flush=True)
    sp, thetap, _ = isolated_build_state(pp, args.mass, zgrid, kh, high_precision=args.high_precision)
    print("STATE 3/3 MINUS", flush=True)
    sm, thetam, _ = isolated_build_state(pm, args.mass, zgrid, kh, high_precision=args.high_precision)

    density_rms = hf.actual_rms(zgrid, kh, s0, sp, sm, "density_proxy")
    theta_ok = bool(theta0 and thetap and thetam)
    theta_rms = hf.actual_rms(zgrid, kh, s0, sp, sm, "theta_proxy") if theta_ok else None

    wake = hd.wake_half_pair_fraction(
        q, f0, fp, fm, args.mass, float(args.z), 200.0
    )
    stats = cro.pair_stats(f0, fp, fm, q, weights)

    summary = {
        "calculation": "checkpointed one-pair hidden-channel validation",
        "direction": args.direction,
        "precision": "high" if args.high_precision else "standard",
        "mass_eV": args.mass,
        "z": args.z,
        "k_range_h_Mpc": [args.kmin, args.kmax],
        "nk": args.nk,
        "density_proxy_half_pair_rms": density_rms,
        "theta_proxy_half_pair_rms": theta_rms,
        "theta_transfer_available": theta_ok,
        "wake_half_pair_fraction_sigma200": wake,
        "wake_to_density_proxy_ratio": wake / max(density_rms, 1e-300),
        "max_relative_moment_mismatch": stats["max_relative_moment_mismatch"],
        "transfer_keys": keys,
        "guardrail": "Single-redshift direct pair validation; no survey S/N is inferred."
    }

    js = out / f"{tag}.json"
    js.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("HIDDEN_CHANNEL_PAIR_CHECK")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {js}")


if __name__ == "__main__":
    main()
