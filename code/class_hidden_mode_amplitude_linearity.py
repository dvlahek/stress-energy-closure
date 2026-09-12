#!/usr/bin/env python3
"""High-precision amplitude-linearity control for the final hidden-mode winner.

The direction is frozen to the final m=0.60 eV fixed-omega0 + fixed-early-Neff
winner.  No re-optimization is performed.  The existing +/-30% pair defines
its midpoint and direction; +/-10% and +/-20% pairs are exact rescalings of
that same direction.  The summary compares both scalar spectra and phi-phi in
the same ideal covariance metrics used by the manuscript diagnostics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import class_scalar_response_optimize as scalar
import class_scalar_fixed_omega0_neff_run as fixed
from class_response_optimize import kinetic_objects, pair_stats, write_psd


def read_psd(path: Path):
    a = np.loadtxt(path)
    if a.ndim != 2 or a.shape[1] < 2:
        raise RuntimeError(f"Unexpected PSD format: {path}")
    return a[:, 0], a[:, 1]


def prepare(source: Path, work: Path, mass_eV: float, z_match: float,
            lmax: int, pkmax: float) -> None:
    work.mkdir(parents=True, exist_ok=True)
    qp, fp30 = read_psd(source / "winner_plus.dat")
    qm, fm30 = read_psd(source / "winner_minus.dat")
    if not np.array_equal(qp, qm):
        raise RuntimeError("Source winner q grids differ")
    q = qp
    f0 = 0.5 * (fp30 + fm30)
    d30 = 0.5 * (fp30 - fm30)
    if np.min(f0) <= 0:
        raise RuntimeError("Non-positive midpoint distribution")

    # Confirm that the source really is the intended final controlled winner.
    cp = json.loads((source / "control_parameters.json").read_text())
    if abs(float(cp["mass_eV"]) - mass_eV) > 1e-12:
        raise RuntimeError(f"Source mass mismatch: {cp['mass_eV']} vs {mass_eV}")
    if abs(float(cp["controlled_early_Neff"]) - 3.046) > 1e-10:
        raise RuntimeError("Source does not have fixed early Neff=3.046")

    params = fixed.control_parameters(mass_eV, z_match)
    fixed.install_writer(params)
    _, weights, _, _, _, _, _ = kinetic_objects(q, mass_eV, z_match)

    manifest = {
        "control": fixed.CONTROL,
        "mass_eV": mass_eV,
        "z_match": z_match,
        "lmax": lmax,
        "pkmax_h_Mpc": pkmax,
        "source_pointwise_cap": 0.30,
        "tested_pointwise_caps": [0.10, 0.20, 0.30],
        "direction_policy": "frozen final HP winner; exact amplitude rescaling; no re-optimization",
        "control_parameters": params,
        "pairs": {},
    }

    for frac, tag in [(0.10, "010"), (0.20, "020"), (0.30, "030")]:
        scale = frac / 0.30
        fp = f0 + scale * d30
        fm = f0 - scale * d30
        if np.min(fp) <= 0 or np.min(fm) <= 0:
            raise RuntimeError(f"Positivity failed for amplitude {frac}")
        stats = pair_stats(f0, fp, fm, q, weights)
        manifest["pairs"][tag] = {
            "fractional_cap": frac,
            "scale_from_30pct_direction": scale,
            "pair_stats": stats,
        }
        if tag == "030":
            # The 30% CLASS spectra are reused exactly from the validated source.
            continue
        pp = work / f"amp{tag}_plus.dat"
        pm = work / f"amp{tag}_minus.dat"
        write_psd(pp, q, fp)
        write_psd(pm, q, fm)
        scalar.write_scalar_ini(
            work / f"amp{tag}_plus.ini", pp, work / f"amp{tag}_plus_",
            mass_eV, lmax, pkmax,
        )
        scalar.write_scalar_ini(
            work / f"amp{tag}_minus.ini", pm, work / f"amp{tag}_minus_",
            mass_eV, lmax, pkmax,
        )

    np.savetxt(
        work / "frozen_direction.csv",
        np.column_stack([q, f0, d30, fp30, fm30]),
        delimiter=",",
        header="q,f_midpoint,delta_at_30pct,f_plus_30pct,f_minus_30pct",
        comments="",
    )
    (work / "amplitude_linearity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


def ttee_norm_and_residual(p, ref, ratio: float, lmin: int, lmax: int):
    sn2_target = 0.0
    sn2_res = 0.0
    used = 0
    for i, l in enumerate(ref["ell"]):
        if l < lmin or l > lmax:
            continue
        T = 0.5 * (ref["TTp"][i] + ref["TTm"][i])
        E = 0.5 * (ref["EEp"][i] + ref["EEm"][i])
        X = 0.5 * (ref["TEp"][i] + ref["TEm"][i])
        if T <= 0 or E <= 0:
            continue
        cov = np.array([
            [2*T*T, 2*X*X, 2*T*X],
            [2*X*X, 2*E*E, 2*E*X],
            [2*T*X, 2*E*X, X*X + T*E],
        ]) / (2*l + 1.0)
        inv = np.linalg.pinv(cov, rcond=1e-12)
        dr = np.array([
            ref["TTp"][i] - ref["TTm"][i],
            ref["EEp"][i] - ref["EEm"][i],
            ref["TEp"][i] - ref["TEm"][i],
        ])
        ds = np.array([
            p["TTp"][i] - p["TTm"][i],
            p["EEp"][i] - p["EEm"][i],
            p["TEp"][i] - p["TEm"][i],
        ])
        target = ratio * dr
        resid = ds - target
        vt = float(target @ inv @ target)
        vr = float(resid @ inv @ resid)
        if np.isfinite(vt) and np.isfinite(vr) and vt >= 0 and vr >= 0:
            sn2_target += vt
            sn2_res += vr
            used += 1
    target_sn = float(np.sqrt(sn2_target))
    residual_sn = float(np.sqrt(sn2_res))
    return {
        "expected_scaled_response_SN": target_sn,
        "nonlinear_residual_SN": residual_sn,
        "relative_CV_weighted_nonlinearity": (
            residual_sn / target_sn if target_sn > 0 else None
        ),
        "multipoles_used": used,
    }


def phiphi_norm_and_residual(p, ref, ratio: float, lmin: int, lmax: int):
    if not np.array_equal(p["ellu"], ref["ellu"]):
        raise RuntimeError("phi-phi ell grids differ")
    target2 = 0.0
    resid2 = 0.0
    used = 0
    for i, l in enumerate(ref["ellu"]):
        if l < lmin or l > lmax:
            continue
        P0 = 0.5 * (ref["PPp"][i] + ref["PPm"][i])
        if P0 == 0 or not np.isfinite(P0):
            continue
        dr = ref["PPp"][i] - ref["PPm"][i]
        ds = p["PPp"][i] - p["PPm"][i]
        target = ratio * dr
        resid = ds - target
        w = 0.5 * (2*l + 1.0) / (P0 * P0)
        target2 += w * target * target
        resid2 += w * resid * resid
        used += 1
    target_sn = float(np.sqrt(max(target2, 0.0)))
    residual_sn = float(np.sqrt(max(resid2, 0.0)))
    return {
        "expected_scaled_response_SN": target_sn,
        "nonlinear_residual_SN": residual_sn,
        "relative_CV_weighted_nonlinearity": (
            residual_sn / target_sn if target_sn > 0 else None
        ),
        "multipoles_used": used,
    }


def summarize(source: Path, work: Path, results: Path, output: Path,
              lmin: int, lmax: int, pkmax: float) -> None:
    manifest = json.loads((work / "amplitude_linearity_manifest.json").read_text())

    # The 30% reference spectra are the already-validated final HP spectra.
    ref_dir = source / "highprec"
    ref_pair = scalar.load_pair(ref_dir, "winner")
    metrics = {
        "030": scalar.actual_metrics(ref_dir, "winner", lmin, lmax, pkmax)
    }
    pairs = {"030": ref_pair}

    for tag in ("010", "020"):
        pairs[tag] = scalar.load_pair(results, f"amp{tag}")
        metrics[tag] = scalar.actual_metrics(results, f"amp{tag}", lmin, lmax, pkmax)

    checks = {}
    ref_sn_ttee = metrics["030"]["TTEE_TE_combined"]["ideal_fullsky_gaussian_cv_SN"]
    ref_sn_pp = metrics["030"]["phiphi"]["ideal_fullsky_cv_SN_auto_only"]
    for tag, ratio in [("010", 1.0/3.0), ("020", 2.0/3.0)]:
        mt = metrics[tag]["TTEE_TE_combined"]["ideal_fullsky_gaussian_cv_SN"]
        mp = metrics[tag]["phiphi"]["ideal_fullsky_cv_SN_auto_only"]
        checks[tag] = {
            "expected_amplitude_ratio_to_30pct": ratio,
            "TTEE_SN_ratio_to_30pct": float(mt / ref_sn_ttee),
            "TTEE_SN_ratio_fractional_error_from_linear": float(mt / ref_sn_ttee / ratio - 1.0),
            "phiphi_SN_ratio_to_30pct": float(mp / ref_sn_pp),
            "phiphi_SN_ratio_fractional_error_from_linear": float(mp / ref_sn_pp / ratio - 1.0),
            "TTEE_full_response_linearity": ttee_norm_and_residual(
                pairs[tag], ref_pair, ratio, lmin, lmax
            ),
            "phiphi_full_response_linearity": phiphi_norm_and_residual(
                pairs[tag], ref_pair, ratio, 8, lmax
            ),
        }

    # A compact manuscript-facing pass/fail: the complete CV-weighted response
    # must follow the frozen-direction linear prediction to better than 5%.
    max_rel = max(
        checks[tag][ch]["relative_CV_weighted_nonlinearity"]
        for tag in checks
        for ch in ("TTEE_full_response_linearity", "phiphi_full_response_linearity")
    )
    out = {
        "analysis": "frozen-direction high-precision amplitude-linearity control",
        "manifest": manifest,
        "metrics_by_pointwise_cap": {
            "0.10": metrics["010"],
            "0.20": metrics["020"],
            "0.30": metrics["030"],
        },
        "linearity_checks": checks,
        "maximum_relative_CV_weighted_nonlinearity": float(max_rel),
        "five_percent_linearity_test_pass": bool(max_rel <= 0.05),
        "scope": (
            "This tests amplitude scaling along the single frozen final winner direction; "
            "it does not re-optimize the hidden direction at each amplitude."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["prepare", "summarize"], required=True)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--mass-eV", type=float, default=0.60)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--lmin", type=int, default=30)
    ap.add_argument("--lmax", type=int, default=1200)
    ap.add_argument("--pkmax", type=float, default=0.30)
    args = ap.parse_args()

    if args.mode == "prepare":
        prepare(args.source, args.work, args.mass_eV, args.z_match, args.lmax, args.pkmax)
    else:
        if args.results is None or args.output is None:
            ap.error("--results and --output are required in summarize mode")
        summarize(args.source, args.work, args.results, args.output,
                  args.lmin, args.lmax, args.pkmax)


if __name__ == "__main__":
    main()
