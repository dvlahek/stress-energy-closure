#!/usr/bin/env python3
r"""Direct-theta validation of hidden kinetic information retention.

This is the velocity-transfer counterpart of hidden_channel_pair_check.py.  It
requests CLASS vTk output explicitly and compares the same source-matched pair
in three transfer-level quantities at one redshift:

1. the existing density-derived relative-velocity proxy,
2. the direct CLASS relative velocity-divergence transfer theta_ncdm-theta_cdm,
3. the density-weighted direct-theta proxy (theta_rel P_cb).

The local resonant wake is evaluated for the same pair using the same half-pair
normalization.  These are response-level diagnostics, not survey S/N forecasts
and not a kSZ forecast.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

import class_response_optimize as cro
import hidden_channel_operator_diagnostics as hd

WORKER = Path(__file__).with_name("class_state_worker_velocity.py")

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
    _, weights, _, _, _, _, _ = cro.kinetic_objects(q, mass, z_match)
    return q, f0, fp, fm, weights


def write_psd(path: Path, q, f):
    hd.write_psd(path, q, f)


def build_state_isolated(psd: Path, mass: float, zgrid, kh):
    workdir = psd.parent
    workdir.mkdir(parents=True, exist_ok=True)

    kh_fd, kh_name = tempfile.mkstemp(prefix="theta_kh_", suffix=".npy", dir=workdir)
    os.close(kh_fd)
    out_fd, out_name = tempfile.mkstemp(prefix="theta_state_", suffix=".npz", dir=workdir)
    os.close(out_fd)
    kh_path = Path(kh_name)
    out_path = Path(out_name)

    try:
        np.save(kh_path, np.asarray(kh, dtype=float))
        cmd = [
            sys.executable,
            str(WORKER),
            "--psd", str(psd.resolve()),
            "--mass", repr(float(mass)),
            "--z-grid", ",".join(repr(float(z)) for z in zgrid),
            "--kh-npy", str(kh_path.resolve()),
            "--out-npz", str(out_path.resolve()),
        ]
        subprocess.run(cmd, check=True)

        with np.load(out_path, allow_pickle=False) as d:
            zs = np.asarray(d["z_grid"], dtype=float)
            state = {}
            for i, z in enumerate(zs):
                theta = np.asarray(d[f"theta_rel_{i}"], dtype=float)
                state[float(z)] = {
                    "pk": np.asarray(d[f"pk_{i}"], dtype=float).copy(),
                    "vrel_density": np.asarray(d[f"vrel_density_{i}"], dtype=float).copy(),
                    "theta_rel": None if theta.size == 0 else theta.copy(),
                }
            return {
                "state": state,
                "theta_available": bool(int(d["theta_available"])),
                "transfer_keys": json.loads(str(d["transfer_keys_json"].item())),
                "cdm_velocity_key": str(d["cdm_velocity_key"].item()),
                "ncdm_velocity_key": str(d["ncdm_velocity_key"].item()),
            }
    finally:
        for p in (kh_path, out_path):
            try:
                p.unlink()
            except FileNotFoundError:
                pass


def stable_mask(x0, threshold=1e-4):
    xmax = max(float(np.max(np.abs(x0))), 1e-300)
    return np.abs(x0) > threshold * xmax


def half_pair_rms(x0, xp, xm, kh, threshold=1e-4):
    mask = stable_mask(x0, threshold)
    if np.sum(mask) < 5:
        return float("nan")
    frac = (xp[mask] - xm[mask]) / (2.0 * x0[mask])
    dlnk = np.gradient(np.log(kh))
    w = kh[mask] ** 3 * dlnk[mask]
    return float(np.sqrt(np.sum(w * frac * frac) / np.sum(w)))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="hidden_channel_theta_checks")
    ap.add_argument("--direction", choices=("cref", "selected"), required=True)
    ap.add_argument("--pair-csv", type=Path,
                    default=Path("hidden_channel_operator_diagnostics/best_proxy_pair.csv"))
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--z", type=float, default=0.30)
    ap.add_argument("--sigma-kms", type=float, default=200.0)
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=48)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    if args.direction == "cref":
        q, f0, fp, fm, weights = reconstruct_cref(args.mass, args.z_match, args.frac)
    else:
        if not args.pair_csv.exists():
            raise FileNotFoundError(
                f"Selected-pair file not found: {args.pair_csv}. "
                "Use hidden_channel_operator_diagnostics/best_proxy_pair.csv"
            )
        q, f0, fp, fm, weights = load_selected(args.pair_csv, args.mass, args.z_match)

    tag = f"{args.direction}_z{args.z:.3f}_theta_std"
    p0 = out / f"{tag}_fd.dat"
    pp = out / f"{tag}_plus.dat"
    pm = out / f"{tag}_minus.dat"
    write_psd(p0, q, f0)
    write_psd(pp, q, fp)
    write_psd(pm, q, fm)

    kh = np.geomspace(args.kmin, args.kmax, args.nk)
    zgrid = [float(args.z)]

    print(f"STATE 1/3 FD direction={args.direction} z={args.z}", flush=True)
    s0 = build_state_isolated(p0, args.mass, zgrid, kh)
    print("STATE 2/3 PLUS", flush=True)
    sp = build_state_isolated(pp, args.mass, zgrid, kh)
    print("STATE 3/3 MINUS", flush=True)
    sm = build_state_isolated(pm, args.mass, zgrid, kh)

    z = float(args.z)
    r0 = s0["state"][z]
    rp = sp["state"][z]
    rm = sm["state"][z]

    density_proxy = half_pair_rms(
        r0["vrel_density"] * r0["pk"],
        rp["vrel_density"] * rp["pk"],
        rm["vrel_density"] * rm["pk"],
        kh,
    )

    theta_ok = bool(s0["theta_available"] and sp["theta_available"] and sm["theta_available"])
    theta_only = None
    theta_pk = None
    if theta_ok:
        theta_only = half_pair_rms(r0["theta_rel"], rp["theta_rel"], rm["theta_rel"], kh)
        theta_pk = half_pair_rms(
            r0["theta_rel"] * r0["pk"],
            rp["theta_rel"] * rp["pk"],
            rm["theta_rel"] * rm["pk"],
            kh,
        )

    wake = hd.wake_half_pair_fraction(
        q, f0, fp, fm, args.mass, z, args.sigma_kms
    )
    stats = cro.pair_stats(f0, fp, fm, q, weights)

    summary = {
        "calculation": "direct CLASS velocity-transfer hidden-channel validation",
        "direction": args.direction,
        "precision": "standard",
        "mass_eV": args.mass,
        "z": z,
        "k_range_h_Mpc": [args.kmin, args.kmax],
        "nk": args.nk,
        "class_output_requested": "mPk,dTk,vTk",
        "theta_transfer_available": theta_ok,
        "cdm_velocity_key": s0["cdm_velocity_key"],
        "ncdm_velocity_key": s0["ncdm_velocity_key"],
        "density_derived_velocity_Pcb_half_pair_rms": density_proxy,
        "direct_theta_half_pair_rms": theta_only,
        "direct_theta_Pcb_half_pair_rms": theta_pk,
        "wake_half_pair_fraction_sigma_kms": wake,
        "wake_to_density_derived_proxy_ratio": wake / max(density_proxy, 1e-300),
        "wake_to_direct_theta_ratio": None if theta_only is None else wake / max(theta_only, 1e-300),
        "wake_to_direct_theta_Pcb_ratio": None if theta_pk is None else wake / max(theta_pk, 1e-300),
        "max_relative_moment_mismatch": stats["max_relative_moment_mismatch"],
        "transfer_keys": s0["transfer_keys"],
        "guardrail": "Transfer-level single-redshift control. Direct theta is not a kSZ or RSD survey forecast."
    }

    js = out / f"{tag}.json"
    js.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("HIDDEN_CHANNEL_THETA_PAIR_CHECK")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {js}")


if __name__ == "__main__":
    main()
