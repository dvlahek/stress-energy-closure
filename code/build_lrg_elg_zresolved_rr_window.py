#!/usr/bin/env python3
"""Build fine DESI LRG-random x ELG-random counts for correlation-function window convolution.

This is deliberately data-blind: it reads only the DESI random catalogs.  The
pair definition is identical to the frozen production estimator (LRG->ELG,
midpoint LOS, theta >= 0.05 deg, DESI fiducial distances, WEIGHT*WEIGHT_FKP).

Fine RR(s,mu) counts are later used to average a fixed theory through exactly
the same finite (s,mu) support as the measured Landy-Szalay correlation
function.  No observed correlation-function values enter this step.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from pycorr import TwoPointCounter

import desi_dr1_lrg_elg_exact as exact
import desi_dr1_lrg_elg_exact_zresolved as zresolved


def parse_edges(text):
    edges = np.asarray([float(x) for x in text.split(",")], dtype="f8")
    if len(edges) < 2 or np.any(np.diff(edges) <= 0):
        raise ValueError("Edges must be strictly increasing")
    return edges


def fine_edges(smin, smax, step):
    n = int(round((smax - smin) / step))
    if n <= 0 or not np.isclose(smin + n * step, smax, rtol=0.0, atol=1e-10):
        raise ValueError("fine-step must exactly divide the requested separation range")
    return np.linspace(smin, smax, n + 1, dtype="f8")


def asymmetry_summary(wcounts, sedges, output_edges):
    """Quantify raw RR mu-asymmetry without using any data vector."""
    out = []
    for lo, hi in zip(output_edges[:-1], output_edges[1:]):
        m = (sedges[:-1] >= lo - 1e-12) & (sedges[1:] <= hi + 1e-12)
        w = np.sum(wcounts[m], axis=0)
        wr = w[::-1]
        den = float(np.sum(np.abs(w)))
        l1 = float(np.sum(np.abs(w - wr)) / den) if den > 0 else None
        out.append({"slo": float(lo), "shi": float(hi), "mu_mirror_l1_fraction": l1})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lrg-random", nargs="+", required=True)
    ap.add_argument("--elg-random", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--z-edges", default="0.80,0.90,1.00,1.10")
    ap.add_argument("--output-sep-edges", default="20,40,60,80,100,120,140")
    ap.add_argument("--fine-step", type=float, default=1.0)
    ap.add_argument("--mu-bins", type=int, default=240)
    ap.add_argument("--theta-min-deg", type=float, default=0.05)
    ap.add_argument("--distance-cosmology", choices=("desi", "legacy_astropy"), default="desi")
    ap.add_argument("--nthreads", type=int, default=max(1, os.cpu_count() or 1))
    args = ap.parse_args()

    zedges = parse_edges(args.z_edges)
    outedges = parse_edges(args.output_sep_edges)
    sedges = fine_edges(float(outedges[0]), float(outedges[-1]), float(args.fine_step))
    muedges = np.linspace(-1.0 - 1e-7, 1.0 + 1e-7, int(args.mu_bins) + 1, dtype="f8")

    exact.set_distance_cosmology(args.distance_cosmology)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("Reading full-range random catalogs once", flush=True)
    lrg_r = exact.read_catalog(args.lrg_random, float(zedges[0]), float(zedges[-1]))
    elg_r = exact.read_catalog(args.elg_random, float(zedges[0]), float(zedges[-1]))

    selection = None
    if args.theta_min_deg > 0:
        selection = {"theta": (float(args.theta_min_deg), 180.0)}

    all_wcounts = []
    all_wnorm = []
    all_counts = []
    all_asymmetry = []

    for iz, (lo, hi) in enumerate(zip(zedges[:-1], zedges[1:])):
        lo, hi = float(lo), float(hi)
        r1 = zresolved.subset(lrg_r, lo, hi)
        r2 = zresolved.subset(elg_r, lo, hi)
        counts = {"LRG_random": int(len(r1[2])), "ELG_random": int(len(r2[2]))}
        if min(counts.values()) == 0:
            raise RuntimeError(f"Empty random slice in z={lo}-{hi}: {counts}")

        print(
            f"Fine RR z-bin {iz}: {lo:.2f}-{hi:.2f}; counts={counts}; "
            f"ns={len(sedges)-1}; nmu={len(muedges)-1}",
            flush=True,
        )

        rr = TwoPointCounter(
            mode="smu",
            edges=(sedges, muedges),
            positions1=exact.rdd(r1),
            positions2=exact.rdd(r2),
            weights1=exact.weights(r1),
            weights2=exact.weights(r2),
            position_type="rdd",
            los="midpoint",
            selection_attrs=selection,
            engine="corrfunc",
            nthreads=int(args.nthreads),
            compute_sepsavg=False,
        )
        w = np.asarray(rr.wcounts, dtype="f8")
        if w.shape != (len(sedges) - 1, len(muedges) - 1):
            raise RuntimeError(f"Unexpected RR shape {w.shape}")
        if not np.all(np.isfinite(w)):
            raise RuntimeError(f"Non-finite RR counts in z-bin {iz}")

        rr.save(str(out / f"pycorr_R1R2_fine_z{iz}.npy"))
        all_wcounts.append(w)
        all_wnorm.append(float(rr.wnorm))
        all_counts.append(counts)
        all_asymmetry.append(asymmetry_summary(w, sedges, outedges))

    W = np.asarray(all_wcounts, dtype="f8")
    np.savez_compressed(
        out / "lrg_elg_zresolved_rr_window_counts.npz",
        zlo=np.asarray(zedges[:-1], dtype="f8"),
        zhi=np.asarray(zedges[1:], dtype="f8"),
        sedges=sedges,
        muedges=muedges,
        wcounts=W,
        wnorm=np.asarray(all_wnorm, dtype="f8"),
        output_sedges=outedges,
    )

    summary = {
        "scope": "Fine random-pair window counts for frozen DESI DR1 LRG->ELG z-resolved estimator",
        "uses_observed_data_vector": False,
        "orientation": "LRG->ELG",
        "los": "midpoint",
        "distance_cosmology": args.distance_cosmology,
        "weights": "WEIGHT * WEIGHT_FKP",
        "z_edges": zedges.tolist(),
        "output_separation_edges_Mpc_over_h": outedges.tolist(),
        "fine_separation_step_Mpc_over_h": float(args.fine_step),
        "fine_separation_bins": int(len(sedges) - 1),
        "mu_bins": int(args.mu_bins),
        "theta_min_deg": float(args.theta_min_deg),
        "random_files": {
            "LRG": [str(Path(x)) for x in args.lrg_random],
            "ELG": [str(Path(x)) for x in args.elg_random],
        },
        "per_zbin_counts": all_counts,
        "rr_wnorm": all_wnorm,
        "mu_asymmetry_by_output_bin": all_asymmetry,
        "analysis_scope": (
            "We construct the random-pair response using the same separation bins, angular selection, "
            "pair orientation and catalogue weights as the measured estimator."
        ),
    }
    (out / "rr_window_counts_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
