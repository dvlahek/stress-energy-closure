#!/usr/bin/env python3
"""Observable eigenspectrum of the seven-dimensional hidden kinetic subspace.

This is a post-processing diagnostic for the fixed-omega0 + fixed-early-Neff
m=0.60 eV control.  It reuses the already-computed symmetric CLASS probes.

The raw null directions are not an orthonormal physical basis, so diagonalizing
the observable Gram matrix alone would depend on their arbitrary individual
normalizations.  We therefore solve the generalized eigenproblem

    G v = lambda H v,

where G is the ideal full-sky response Gram matrix and H is a kinetic-state
metric.  H is the number-density-weighted RMS fractional-deformation metric

    ||delta f||_H^2 = int q^2 f0 (delta f/f0)^2 dq / int q^2 f0 dq.

For manuscript-facing comparison, every eigenmode is also rescaled to the same
30% maximum pointwise fractional-deformation cap used for the optimized pair.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import class_scalar_response_optimize as scalar


def trap(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def kinetic_metric(q: np.ndarray, f0: np.ndarray, shapes: np.ndarray) -> np.ndarray:
    denom = trap(q * q * f0, q)
    if not np.isfinite(denom) or denom <= 0:
        raise RuntimeError("Invalid kinetic-metric normalization")
    weight = q * q / np.maximum(f0, 1e-300)
    nd = shapes.shape[0]
    H = np.empty((nd, nd), dtype=float)
    for i in range(nd):
        for j in range(i, nd):
            val = trap(weight * shapes[i] * shapes[j], q) / denom
            H[i, j] = H[j, i] = val
    return 0.5 * (H + H.T)


def generalized_modes(G: np.ndarray, H: np.ndarray):
    G = 0.5 * (G + G.T)
    H = 0.5 * (H + H.T)
    hw, hu = np.linalg.eigh(H)
    htol = max(float(np.max(np.abs(hw))) * 1e-12, 1e-15)
    keep = hw > htol
    if int(np.sum(keep)) != H.shape[0]:
        raise RuntimeError(f"Kinetic metric is rank deficient: eigenvalues={hw.tolist()}")
    B = hu[:, keep] @ np.diag(1.0 / np.sqrt(hw[keep])) @ hu[:, keep].T
    A = 0.5 * (B @ G @ B + (B @ G @ B).T)
    w, u = np.linalg.eigh(A)
    order = np.argsort(w)[::-1]
    w = np.maximum(w[order], 0.0)
    u = u[:, order]
    coeff = B @ u
    # Numerical cleanup: enforce unit H norm column by column.
    for j in range(coeff.shape[1]):
        n2 = float(coeff[:, j] @ H @ coeff[:, j])
        coeff[:, j] /= np.sqrt(max(n2, 1e-300))
    return w, coeff


def sn_from_gram(c: np.ndarray, G: np.ndarray, scale: float = 1.0) -> float:
    return float(abs(scale) * np.sqrt(max(float(c @ G @ c), 0.0)))


def analyze_channel(name: str, G: np.ndarray, Gs: np.ndarray, Gp: np.ndarray,
                    H: np.ndarray, shapes: np.ndarray, f0: np.ndarray,
                    probe_frac: float, final_frac: float):
    eig, vec = generalized_modes(G, H)
    total = float(np.sum(eig))
    cumulative = 0.0
    modes = []
    for rank, (lam, c) in enumerate(zip(eig, vec.T), start=1):
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            raise RuntimeError(f"Invalid pointwise amplitude for {name} mode {rank}")
        cap_coeff_scale = final_frac / (probe_frac * maxrel)
        frac = float(lam / total) if total > 0 else 0.0
        cumulative += frac
        modes.append({
            "rank": rank,
            "generalized_eigenvalue_SN2_at_probe_scale": float(lam),
            "SN_at_probe_scale_H_unit": float(np.sqrt(lam)),
            "information_fraction": frac,
            "cumulative_information_fraction": float(cumulative),
            "coefficients_in_stored_null_basis": c.tolist(),
            "max_fractional_shape_per_H_unit": maxrel,
            "coefficient_scale_to_30pct_cap": float(cap_coeff_scale),
            "one_sided_H_rms_fraction_at_30pct_cap": float(final_frac / maxrel),
            "predicted_TTEE_SN_at_30pct_cap": sn_from_gram(c, Gs, cap_coeff_scale),
            "predicted_phiphi_SN_at_30pct_cap": sn_from_gram(c, Gp, cap_coeff_scale),
            "predicted_block_diagonal_TTEE_plus_phiphi_SN_at_30pct_cap":
                float(np.hypot(sn_from_gram(c, Gs, cap_coeff_scale),
                               sn_from_gram(c, Gp, cap_coeff_scale))),
        })
    tol = max(float(eig[0]) * 1e-10 if len(eig) else 0.0, 1e-18)
    return {
        "channel": name,
        "rank_above_1e-10_relative": int(np.sum(eig > tol)),
        "eigenvalues_SN2": eig.tolist(),
        "top_mode_information_fraction": modes[0]["information_fraction"] if modes else 0.0,
        "top_two_cumulative_information_fraction": modes[min(1, len(modes)-1)]["cumulative_information_fraction"] if modes else 0.0,
        "modes": modes,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--lmin", type=int, default=30)
    ap.add_argument("--lmax", type=int, default=600)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    meta_path = args.work / "scalar_probe_meta.json"
    if not meta_path.exists():
        raise FileNotFoundError(meta_path)
    meta = json.loads(meta_path.read_text())
    probe_frac = float(meta["probe_fractional_distortion"])
    shapes = np.load(args.work / "null_shapes.npy")
    q = np.load(args.work / "q.npy")
    f0 = np.load(args.work / "f0.npy")
    if shapes.shape[0] != 7:
        raise RuntimeError(f"Expected seven hidden directions, got {shapes.shape[0]}")

    Gs, Gp = scalar.build_gram(args.work, args.lmin, args.lmax)
    H = kinetic_metric(q, f0, shapes)
    Gsum = Gs + Gp

    heig = np.linalg.eigvalsh(H)
    out = {
        "analysis": "seven-dimensional hidden-state observable generalized eigenspectrum",
        "source_mass_eV": float(meta["mass_eV"]),
        "source_z_match": float(meta["z_match"]),
        "null_dimension": int(shapes.shape[0]),
        "probe_fractional_distortion": probe_frac,
        "comparison_pointwise_fractional_cap": float(args.final_frac),
        "ell_range_TTEE": [int(args.lmin), int(args.lmax)],
        "ell_range_phiphi": [8, int(args.lmax)],
        "kinetic_metric": (
            "number-density-weighted RMS fractional deformation: "
            "int q^2 f0 (delta f/f0)^2 dq / int q^2 f0 dq"
        ),
        "kinetic_metric_eigenvalues": heig.tolist(),
        "important_scope_note": (
            "The TTEE+phiphi spectrum treats the two response blocks as independent and is an "
            "optimistic diagnostic only; manuscript claims should use the separate TTEE and phiphi spectra."
        ),
        "TTEE": analyze_channel("TTEE", Gs, Gs, Gp, H, shapes, f0, probe_frac, args.final_frac),
        "phiphi": analyze_channel("phiphi", Gp, Gs, Gp, H, shapes, f0, probe_frac, args.final_frac),
        "block_diagonal_TTEE_plus_phiphi": analyze_channel(
            "block_diagonal_TTEE_plus_phiphi", Gsum, Gs, Gp, H, shapes, f0,
            probe_frac, args.final_frac),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")

    rows = []
    for ch in ("TTEE", "phiphi", "block_diagonal_TTEE_plus_phiphi"):
        for m in out[ch]["modes"]:
            rows.append([
                ch, m["rank"], m["generalized_eigenvalue_SN2_at_probe_scale"],
                m["information_fraction"], m["cumulative_information_fraction"],
                m["one_sided_H_rms_fraction_at_30pct_cap"],
                m["predicted_TTEE_SN_at_30pct_cap"],
                m["predicted_phiphi_SN_at_30pct_cap"],
                m["predicted_block_diagonal_TTEE_plus_phiphi_SN_at_30pct_cap"],
            ])
    csv_path = args.output.with_suffix(".csv")
    with csv_path.open("w", encoding="utf-8") as fh:
        fh.write("channel,rank,eigenvalue_SN2_at_probe_scale,information_fraction,cumulative_information_fraction,one_sided_H_rms_fraction_at_30pct_cap,predicted_TTEE_SN_at_30pct_cap,predicted_phiphi_SN_at_30pct_cap,predicted_blockdiag_SN_at_30pct_cap\n")
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
