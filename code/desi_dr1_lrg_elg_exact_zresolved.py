#!/usr/bin/env python3
"""Exact redshift-resolved DESI DR1 LRG x ELG odd multipoles.

Frozen production split:
  0.80 <= z < 0.90
  0.90 <= z < 1.00
  1.00 <= z < 1.10

Uses the same pycorr/Corrfunc cross Landy-Szalay estimator, midpoint LOS,
weights, angular cut, separation bins and mu grid as desi_dr1_lrg_elg_exact.py.
Catalogs are read once over the full redshift range and then sliced in memory.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

import desi_dr1_lrg_elg_exact as exact


def subset(cat, lo, hi):
    m = (cat[2] >= lo) & (cat[2] < hi)
    return tuple(np.asarray(x[m]) for x in cat)


def parse_edges(text):
    e = np.asarray([float(x) for x in text.split(",")], dtype="f8")
    if len(e) < 2 or np.any(np.diff(e) <= 0):
        raise ValueError("Redshift edges must be strictly increasing")
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lrg-data", nargs="+", required=True)
    ap.add_argument("--elg-data", nargs="+", required=True)
    ap.add_argument("--lrg-random", nargs="+", required=True)
    ap.add_argument("--elg-random", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--z-edges", default="0.80,0.90,1.00,1.10")
    ap.add_argument("--sep-edges", default="20,40,60,80,100,120,140")
    ap.add_argument("--mu-bins", type=int, default=240)
    ap.add_argument("--theta-min-deg", type=float, default=0.05)
    ap.add_argument("--nthreads", type=int, default=max(1, os.cpu_count() or 1))
    ap.add_argument("--skip-reverse", action="store_true")
    args = ap.parse_args()

    zedges = parse_edges(args.z_edges)
    sedges = np.asarray([float(x) for x in args.sep_edges.split(",")], dtype="f8")
    if np.any(np.diff(sedges) <= 0):
        raise ValueError("Separation edges must be strictly increasing")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    zmin, zmax = float(zedges[0]), float(zedges[-1])
    print("Reading full-range catalogs once", flush=True)
    lrg = exact.read_catalog(args.lrg_data, zmin, zmax)
    elg = exact.read_catalog(args.elg_data, zmin, zmax)
    lrg_r = exact.read_catalog(args.lrg_random, zmin, zmax)
    elg_r = exact.read_catalog(args.elg_random, zmin, zmax)

    rows = []
    summaries = []
    reverse_rows = []

    for iz, (lo, hi) in enumerate(zip(zedges[:-1], zedges[1:])):
        lo, hi = float(lo), float(hi)
        d1, d2 = subset(lrg, lo, hi), subset(elg, lo, hi)
        r1, r2 = subset(lrg_r, lo, hi), subset(elg_r, lo, hi)
        counts = {
            "LRG": int(len(d1[2])),
            "ELG": int(len(d2[2])),
            "LRG_random": int(len(r1[2])),
            "ELG_random": int(len(r2[2])),
        }
        if min(counts.values()) == 0:
            raise RuntimeError(f"Empty catalog slice in z={lo}-{hi}: {counts}")

        zeff = exact.weighted_zeff(d1, d2)
        print(f"Exact z-bin {iz}: {lo:.2f}-{hi:.2f}; zeff_proxy={zeff:.6f}; counts={counts}", flush=True)
        corr, sep, xi1, xi3 = exact.run_corr(
            d1, d2, r1, r2, sedges, args.mu_bins, args.nthreads, args.theta_min_deg
        )
        exact.finite_or_raise(f"zbin{iz} LRG->ELG", sep, xi1, xi3)
        corr.save(str(out / f"pycorr_lrg_to_elg_z{iz}.npy"))

        reverse = None
        if not args.skip_reverse:
            corr_rev, sep_rev, x1r, x3r = exact.run_corr(
                d2, d1, r2, r1, sedges, args.mu_bins, args.nthreads, args.theta_min_deg
            )
            exact.finite_or_raise(f"zbin{iz} ELG->LRG", sep_rev, x1r, x3r)
            corr_rev.save(str(out / f"pycorr_elg_to_lrg_z{iz}.npy"))
            if not np.allclose(sep, sep_rev, rtol=0.0, atol=1e-12):
                raise RuntimeError(f"Forward/reverse separation mismatch in z-bin {iz}")
            eps = 1e-300
            reverse = {
                "dipole_max_abs_forward_plus_reverse": float(np.max(np.abs(xi1 + x1r))),
                "octupole_max_abs_forward_plus_reverse": float(np.max(np.abs(xi3 + x3r))),
                "dipole_relative_sign_residual": float(np.max(np.abs(xi1 + x1r)) / max(np.max(np.abs(xi1)), eps)),
                "octupole_relative_sign_residual": float(np.max(np.abs(xi3 + x3r)) / max(np.max(np.abs(xi3)), eps)),
            }
            for s, a, ar, b, br in zip(sep, xi1, x1r, xi3, x3r):
                reverse_rows.append((lo, hi, zeff, s, a, ar, b, br))

        for s, a, b in zip(sep, xi1, xi3):
            rows.append((lo, hi, zeff, s, a, b))
        summaries.append({
            "index": iz,
            "zlo": lo,
            "zhi": hi,
            "object_weighted_z_proxy": zeff,
            "counts": counts,
            "sign_reversal": reverse,
        })

    np.savetxt(
        out / "lrg_elg_exact_zresolved_odd_multipoles.csv",
        np.asarray(rows, float),
        delimiter=",",
        header="zlo,zhi,z_effective,s_Mpc_over_h,xi1_LRG_to_ELG,xi3_LRG_to_ELG",
        comments="",
    )
    if reverse_rows:
        np.savetxt(
            out / "lrg_elg_exact_zresolved_forward_reverse.csv",
            np.asarray(reverse_rows, float),
            delimiter=",",
            header="zlo,zhi,z_effective,s_Mpc_over_h,xi1_LRG_to_ELG,xi1_ELG_to_LRG,xi3_LRG_to_ELG,xi3_ELG_to_LRG",
            comments="",
        )

    summary = {
        "scope": "Exact-pair redshift-resolved genuine DESI DR1 LRGxELG odd-sector vector",
        "frozen_z_edges": zedges.tolist(),
        "vector_order": "z-bin major, separation-bin minor",
        "multipoles": [1, 3],
        "orientation": "LRG->ELG",
        "los": "midpoint",
        "separation_edges_Mpc_over_h": sedges.tolist(),
        "mu_bins": int(args.mu_bins),
        "theta_min_deg": float(args.theta_min_deg),
        "weights": "WEIGHT * WEIGHT_FKP for data and randoms",
        "bins": summaries,
        "guardrail": (
            "The three redshift bins are frozen before inspection of the z-resolved data. "
            "Use identically processed mocks for covariance; do not merge, move or select bins "
            "after inspecting which choice increases significance."
        ),
    }
    (out / "summary_exact_zresolved.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
