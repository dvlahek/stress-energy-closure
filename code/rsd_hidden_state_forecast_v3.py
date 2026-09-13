#!/usr/bin/env python3
"""Nonlinear validation/optimization of the matched-state RSD handle.

This extends rsd_hidden_state_forecast.py by ranking many candidate directions
with the small-amplitude response Jacobian and then evaluating the best
candidates with full CLASS runs at the final 30 per cent pointwise cap.
The objective is the nuisance-projected Kaiser redshift-space S/N after
marginalizing H0, omega_b, omega_cdm, ln A_s, n_s and one independent galaxy
bias in each redshift bin.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import rsd_hidden_state_forecast as base
import class_response_optimize as cro


def candidate_pool(Bproj, shapes, f0, final_frac, samples, seed, keep):
    nd = shapes.shape[0]
    evals, evecs = np.linalg.eigh(Bproj)
    candidates = [evecs[:, j] for j in np.argsort(evals)[::-1]]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0
        candidates.extend([e, -e])
    rng = np.random.default_rng(seed)
    rr = rng.normal(size=(samples, nd))
    rr /= np.linalg.norm(rr, axis=1, keepdims=True)
    candidates.extend(rr)

    scored = []
    for c in candidates:
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            continue
        norm = 1.0 / maxrel
        amp = 2.0 * final_frac * norm
        s2 = float(amp * amp * (c @ Bproj @ c))
        scored.append((s2, c.copy(), norm, raw * norm))
    scored.sort(key=lambda x: x[0], reverse=True)

    selected = []
    for item in scored:
        c = item[1] / np.linalg.norm(item[1])
        if all(abs(float(c @ (q[1] / np.linalg.norm(q[1])))) < 0.9995 for q in selected):
            selected.append(item)
        if len(selected) >= keep:
            break
    return selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="rsd_hidden_state_output_v3")
    ap.add_argument("--mass", type=float, default=0.60)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.01)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--volume", type=float, default=50.0)
    ap.add_argument("--samples", type=int, default=30000)
    ap.add_argument("--keep", type=int, default=28)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    q = np.linspace(0.0, 20.0, 4000)
    f0, weights_mom, basis, N, shapes, y, M = cro.kinetic_objects(q, args.mass, args.z_match)
    f0_path = out / "fd_reference.dat"
    base.write_psd(f0_path, q, f0)

    ref = base.observable_from_psd(f0_path, args.mass)
    fisher_w = base.integration_weights(args.volume)
    nuisance, nuisance_names = base.nuisance_matrix(f0_path, args.mass, ref)
    qn = base.weighted_projector_basis(nuisance, fisher_w)

    response_cols = []
    for j, shape in enumerate(shapes):
        fp = f0 + args.probe_frac * shape
        fm = f0 - args.probe_frac * shape
        pp = out / f"probe_{j:02d}_plus.dat"
        pm = out / f"probe_{j:02d}_minus.dat"
        base.write_psd(pp, q, fp)
        base.write_psd(pm, q, fm)
        op = base.observable_from_psd(pp, args.mass)
        om = base.observable_from_psd(pm, args.mass)
        response_cols.append((op - om) / (2.0 * args.probe_frac))
    R = np.column_stack(response_cols)
    rw = np.sqrt(fisher_w)[:, None] * R
    rproj = rw - qn @ (qn.T @ rw)
    Bproj = rproj.T @ rproj

    selected = candidate_pool(Bproj, shapes, f0, args.final_frac, args.samples, args.seed, args.keep)
    rows = []
    best = None
    for rank, (pred2, coeff, norm, shape) in enumerate(selected):
        fp = f0 + args.final_frac * shape
        fm = f0 - args.final_frac * shape
        pp = out / f"cand_{rank:02d}_plus.dat"
        pm = out / f"cand_{rank:02d}_minus.dat"
        base.write_psd(pp, q, fp)
        base.write_psd(pm, q, fm)
        op = base.observable_from_psd(pp, args.mass)
        om = base.observable_from_psd(pm, args.mass)
        fixed, projected, retained = base.sn_metrics(op - om, fisher_w, qn)
        row = {
            "rank": rank,
            "predicted_projected_SN": float(np.sqrt(max(pred2, 0.0))),
            "validated_fixed_SN": fixed,
            "validated_projected_SN": projected,
            "retained_delta_chi2_fraction": retained,
            "normalization": float(norm),
            "coefficients": coeff.tolist(),
        }
        rows.append(row)
        if best is None or projected > best[0]:
            best = (projected, row, shape.copy(), fp.copy(), fm.copy())
        print(json.dumps(row))

    if best is None:
        raise RuntimeError("No candidate validated")
    projected, brow, shape, fp, fm = best
    base.write_psd(out / "best_plus.dat", q, fp)
    base.write_psd(out / "best_minus.dat", q, fm)
    np.savetxt(out / "best_pair.csv", np.column_stack([q, f0, fp, fm, shape]), delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape", comments="")

    mp = cro.moments(fp, q, weights_mom)
    mm = cro.moments(fm, q, weights_mom)
    denom = np.maximum(0.5 * (np.abs(mp) + np.abs(mm)), 1e-300)
    mismatch = np.abs(mp - mm) / denom

    summary = {
        "class_commit": cro.CLASS_COMMIT,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "mass_over_Tnu_at_match": float(y),
        "redshift_bins": base.Z_BINS.tolist(),
        "k_h_Mpc_min": float(base.K_H.min()),
        "k_h_Mpc_max": float(base.K_H.max()),
        "total_effective_volume_hminus3_Gpc3": args.volume,
        "nuisance_parameters": nuisance_names,
        "final_fractional_distortion_cap": args.final_frac,
        "max_relative_moment_mismatch": float(mismatch.max()),
        "nonlinear_candidates_validated": len(rows),
        "best": brow,
        "all_candidates": rows,
        "interpretation": "Ideal CV-limited linear Kaiser RSD diagnostic with full CLASS nonlinear-in-distribution validation at the final cap, not a survey likelihood."
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("BEST")
    print(json.dumps(summary["best"], indent=2))


if __name__ == "__main__":
    main()
