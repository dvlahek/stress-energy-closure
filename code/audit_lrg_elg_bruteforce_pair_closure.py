#!/usr/bin/env python3
"""Independent brute-force closure for the exact DESI LRGxELG estimator.

This deliberately does NOT reuse pycorr pair counts. It:
  1. reads one small deterministic catalog subset,
  2. computes weighted D1D2, D1R2, R1D2 and R1R2 histograms in pure NumPy,
  3. forms cross Landy-Szalay manually,
  4. projects to ell=1,3 with the same bin-integrated Legendre convention,
  5. compares every stage against pycorr/Corrfunc on the identical subset.

The test is numerical/algorithmic only. It does not evaluate wake significance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import special

import desi_dr1_lrg_elg_exact as exact


def sample_catalog(cat, n, seed):
    nall = len(cat[2])
    n = min(int(n), nall)
    rng = np.random.default_rng(int(seed))
    idx = np.sort(rng.choice(nall, size=n, replace=False))
    return tuple(np.asarray(x[idx]) for x in cat), idx


def hash_catalog(cat):
    h = hashlib.sha256()
    for arr in cat:
        a = np.ascontiguousarray(np.asarray(arr, dtype="f8"))
        h.update(a.tobytes())
    return h.hexdigest()


def xyz_unit(cat):
    ra, dec, z, w = cat
    chi = exact.distance_mpc_over_h(z)
    rar = np.deg2rad(ra)
    decr = np.deg2rad(dec)
    c = np.cos(decr)
    unit = np.column_stack([c * np.cos(rar), c * np.sin(rar), np.sin(decr)])
    xyz = chi[:, None] * unit
    return xyz, unit, np.asarray(w, dtype="f8")


def brute_pair_counts(cat1, cat2, sedges, muedges, theta_min_deg, block=256):
    """Weighted cross-pair histogram using Corrfunc's midpoint convention.

    With x1 and x2 the two comoving position vectors:
      svec = x2 - x1
      lvec = (x1 + x2) / 2
      mu = dot(svec,lvec) / (|svec||lvec|)
    This matches pycorr's reference convention: the separation vector points
    from catalog 1 to catalog 2. The factor 1/2 cancels from mu.
    """
    xyz1, unit1, w1 = xyz_unit(cat1)
    xyz2, unit2, w2 = xyz_unit(cat2)

    ns = len(sedges) - 1
    nm = len(muedges) - 1
    hist = np.zeros((ns, nm), dtype="f8")
    r1sq = np.einsum("ij,ij->i", xyz1, xyz1)
    r2sq = np.einsum("ij,ij->i", xyz2, xyz2)
    cos_theta_max = np.cos(np.deg2rad(theta_min_deg))

    for i0 in range(0, len(xyz1), block):
        i1 = min(i0 + block, len(xyz1))
        x1 = xyz1[i0:i1]
        u1 = unit1[i0:i1]
        ww1 = w1[i0:i1]

        dot = x1 @ xyz2.T
        s2 = r1sq[i0:i1, None] + r2sq[None, :] - 2.0 * dot
        m2 = r1sq[i0:i1, None] + r2sq[None, :] + 2.0 * dot
        np.maximum(s2, 0.0, out=s2)
        np.maximum(m2, 0.0, out=m2)
        s = np.sqrt(s2)

        denom = np.sqrt(s2 * m2)
        mu = np.zeros_like(s)
        good_denom = denom > 0.0
        mu[good_denom] = (
            r2sq[None, :] - r1sq[i0:i1, None]
        )[good_denom] / denom[good_denom]

        cosang = u1 @ unit2.T
        si = np.searchsorted(sedges, s, side="right") - 1
        mi = np.searchsorted(muedges, mu, side="right") - 1

        valid = (
            good_denom
            & (si >= 0) & (si < ns)
            & (mi >= 0) & (mi < nm)
            & (cosang > -1.0)
            & (cosang <= cos_theta_max)
        )
        if not np.any(valid):
            continue

        flat = si[valid] * nm + mi[valid]
        pairw = (ww1[:, None] * w2[None, :])[valid]
        hist += np.bincount(flat, weights=pairw, minlength=ns * nm).reshape(ns, nm)

    norm = float(np.sum(w1, dtype="f8") * np.sum(w2, dtype="f8"))
    return hist, norm


def compare_array(manual, reference):
    manual = np.asarray(manual, dtype="f8")
    reference = np.asarray(reference, dtype="f8")
    diff = manual - reference
    scale = max(float(np.max(np.abs(reference))), 1e-300)
    denom_l1 = max(float(np.sum(np.abs(reference))), 1e-300)
    return {
        "max_abs": float(np.max(np.abs(diff))),
        "max_abs_over_reference_peak": float(np.max(np.abs(diff)) / scale),
        "l1_relative": float(np.sum(np.abs(diff)) / denom_l1),
        "rms_abs": float(np.sqrt(np.mean(diff * diff))),
    }


def manual_landy_szalay(counts, norms):
    n = {k: counts[k] / norms[k] for k in counts}
    out = np.full_like(n["R1R2"], np.nan, dtype="f8")
    good = counts["R1R2"] != 0.0
    out[good] = (
        n["D1D2"][good]
        - n["D1R2"][good]
        - n["R1D2"][good]
        + n["R1R2"][good]
    ) / n["R1R2"][good]
    return out


def manual_poles(xi, muedges, ells=(1, 3)):
    dmu = np.diff(muedges)
    poles = []
    for ell in ells:
        integ = special.legendre(ell).integ()(muedges)
        legendre = (2 * ell + 1) * np.diff(integ)
        poles.append(np.sum(xi * legendre[None, :], axis=-1) / np.sum(dmu))
    return np.asarray(poles, dtype="f8")


def scalar_rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lrg-data", nargs="+", required=True)
    ap.add_argument("--elg-data", nargs="+", required=True)
    ap.add_argument("--lrg-random", nargs="+", required=True)
    ap.add_argument("--elg-random", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--zlo", type=float, default=0.80)
    ap.add_argument("--zhi", type=float, default=0.90)
    ap.add_argument("--ndata", type=int, default=1200)
    ap.add_argument("--nrandom", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--sep-edges", default="20,40,60,80,100,120,140")
    ap.add_argument("--mu-bins", type=int, default=24)
    ap.add_argument("--theta-min-deg", type=float, default=0.05)
    ap.add_argument("--distance-cosmology", choices=("desi", "legacy_astropy"), default="desi")
    ap.add_argument("--nthreads", type=int, default=8)
    ap.add_argument("--block", type=int, default=256)
    ap.add_argument("--tol-pair-rel", type=float, default=5e-9)
    ap.add_argument("--tol-xi-abs", type=float, default=5e-9)
    ap.add_argument("--tol-pole-abs", type=float, default=5e-9)
    args = ap.parse_args()

    exact.set_distance_cosmology(args.distance_cosmology)
    sedges = np.asarray([float(x) for x in args.sep_edges.split(",")], dtype="f8")
    muedges = np.linspace(-1.0 - 1e-7, 1.0 + 1e-7, args.mu_bins + 1)

    print("Reading closure catalogs", flush=True)
    D1 = exact.read_catalog(args.lrg_data, args.zlo, args.zhi)
    D2 = exact.read_catalog(args.elg_data, args.zlo, args.zhi)
    R1 = exact.read_catalog(args.lrg_random, args.zlo, args.zhi)
    R2 = exact.read_catalog(args.elg_random, args.zlo, args.zhi)

    D1, iD1 = sample_catalog(D1, args.ndata, args.seed + 11)
    D2, iD2 = sample_catalog(D2, args.ndata, args.seed + 22)
    R1, iR1 = sample_catalog(R1, args.nrandom, args.seed + 33)
    R2, iR2 = sample_catalog(R2, args.nrandom, args.seed + 44)

    subsets = {"D1": D1, "D2": D2, "R1": R1, "R2": R2}
    print("Subset sizes", {k: len(v[2]) for k, v in subsets.items()}, flush=True)

    print("Running pycorr/Corrfunc on identical subset", flush=True)
    corr, sep_py, xi1_py, xi3_py = exact.run_corr(
        D1, D2, R1, R2, sedges, args.mu_bins, args.nthreads, args.theta_min_deg
    )

    pairs = {
        "D1D2": (D1, D2),
        "D1R2": (D1, R2),
        "R1D2": (R1, D2),
        "R1R2": (R1, R2),
    }
    manual_counts, manual_norms = {}, {}
    for name, (a, b) in pairs.items():
        print(f"Brute-force {name}", flush=True)
        manual_counts[name], manual_norms[name] = brute_pair_counts(
            a, b, sedges, muedges, args.theta_min_deg, block=args.block
        )

    pair_metrics = {}
    norm_metrics = {}
    for name in pairs:
        refcount = np.asarray(getattr(corr, name).wcounts, dtype="f8")
        refnorm = float(np.asarray(getattr(corr, name).wnorm))
        pair_metrics[name] = compare_array(manual_counts[name], refcount)
        norm_metrics[name] = {
            "manual": manual_norms[name],
            "pycorr": refnorm,
            "relative_difference": scalar_rel(manual_norms[name], refnorm),
        }

    xi_manual = manual_landy_szalay(manual_counts, manual_norms)
    xi_py = np.asarray(corr.corr, dtype="f8")
    common = np.isfinite(xi_manual) & np.isfinite(xi_py)
    if not np.any(common):
        raise RuntimeError("No finite common (s,mu) bins in brute-force closure")
    xi_metrics = compare_array(xi_manual[common], xi_py[common])

    full_rows = np.all(np.isfinite(xi_manual), axis=1) & np.all(np.isfinite(xi_py), axis=1)
    if not np.any(full_rows):
        raise RuntimeError(
            "No separation row has all mu bins populated; increase --nrandom or lower --mu-bins."
        )

    poles_manual = manual_poles(xi_manual[full_rows], muedges, ells=(1, 3))
    poles_py = np.vstack([xi1_py[full_rows], xi3_py[full_rows]])
    pole_metrics = {
        "ell1": compare_array(poles_manual[0], poles_py[0]),
        "ell3": compare_array(poles_manual[1], poles_py[1]),
    }

    sep_expected = 0.5 * (sedges[:-1] + sedges[1:])
    sep_metrics = compare_array(sep_expected[full_rows], np.asarray(sep_py)[full_rows])

    max_pair_rel = max(v["max_abs_over_reference_peak"] for v in pair_metrics.values())
    max_norm_rel = max(v["relative_difference"] for v in norm_metrics.values())
    max_pole_abs = max(v["max_abs"] for v in pole_metrics.values())
    passed = bool(
        max_pair_rel <= args.tol_pair_rel
        and max_norm_rel <= args.tol_pair_rel
        and xi_metrics["max_abs"] <= args.tol_xi_abs
        and max_pole_abs <= args.tol_pole_abs
    )

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    corr.save(str(out / "pycorr_subset.npy"))

    result = {
        "scope": "Independent brute-force pair-count closure of the DESI LRGxELG exact estimator",
        "status": "PASS" if passed else "FAIL",
        "purpose": "Numerical estimator validation only; no wake significance is evaluated.",
        "subset": {
            "region": "as specified by input files",
            "z_range": [args.zlo, args.zhi],
            "sizes": {k: int(len(v[2])) for k, v in subsets.items()},
            "seed": int(args.seed),
            "hashes_sha256": {k: hash_catalog(v) for k, v in subsets.items()},
            "sample_index_ranges": {
                "D1": [int(iD1.min()), int(iD1.max())],
                "D2": [int(iD2.min()), int(iD2.max())],
                "R1": [int(iR1.min()), int(iR1.max())],
                "R2": [int(iR2.min()), int(iR2.max())],
            },
        },
        "estimator": {
            "orientation": "LRG->ELG",
            "pair_vector_convention": "x2-x1 (catalog 1 -> catalog 2; matches pinned pycorr reference implementation)",
            "los": "midpoint",
            "landy_szalay": "(D1D2-D1R2-R1D2+R1R2)/R1R2 using separately normalized weighted cross-counts",
            "distance_cosmology": args.distance_cosmology,
            "separation_edges_Mpc_over_h": sedges.tolist(),
            "mu_edges": muedges.tolist(),
            "mu_bins": int(args.mu_bins),
            "theta_min_deg": float(args.theta_min_deg),
            "weights": "WEIGHT * WEIGHT_FKP",
        },
        "normalization_comparison": norm_metrics,
        "weighted_pair_count_comparison": pair_metrics,
        "smu_correlation_comparison": {
            **xi_metrics,
            "finite_common_bins": int(np.sum(common)),
            "total_bins": int(xi_py.size),
        },
        "multipole_comparison": {
            **pole_metrics,
            "full_finite_separation_rows": np.flatnonzero(full_rows).astype(int).tolist(),
            "separation_grid_comparison": sep_metrics,
        },
        "thresholds": {
            "pair_and_norm_relative": args.tol_pair_rel,
            "smu_max_abs": args.tol_xi_abs,
            "pole_max_abs": args.tol_pole_abs,
        },
        "guardrail": (
            "This closure must not be used to tune redshift bins, separation bins, nuisance parameters, "
            "or wake templates. A failure indicates an estimator implementation issue to diagnose before "
            "unblinding the 200-mock wake fit."
        ),
    }
    (out / "bruteforce_pair_closure_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
