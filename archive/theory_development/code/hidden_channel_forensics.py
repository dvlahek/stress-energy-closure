#!/usr/bin/env python3
r"""Cross-validated forensics for the tiny hidden-state linear-transfer response.

This script is intentionally diagnostic.  It addresses a specific failure mode
seen in ``hidden_channel_operator_diagnostics.py``: a small-probe response
operator predicted a ~3 percent optimized proxy response, while direct CLASS
validation at the final 30 percent cap remained ~6e-4 and was non-monotonic in
amplitude.

The calculation therefore asks if the apparent optimized mode is physical,
nonlinear, or a numerical/noise overfit.  It performs four independent checks:

1. Rebuild the seven-dimensional response operator at two probe amplitudes and
   compare eigenvalues, optimized directions and cross-predictions.
2. Measure the same hidden response from the direct CLASS velocity-divergence
   transfer, theta_ncdm-theta_cdm, avoiding the redshift finite-difference used
   to infer v_rel from density transfers.  Only fractional pair responses are
   compared, so any common theta normalization cancels.
3. Validate CREF and every optimized direction with full +/- CLASS runs at the
   final cap instead of extrapolating the small-probe Jacobian.
4. Repeat the decisive final-cap runs with the repository high-precision CLASS
   settings to test if the ~1e-3 response is above the numerical floor.

No survey S/N is inferred here.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from classy import Class

import class_response_optimize as cro
import hidden_channel_operator_diagnostics as hd

CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)

HIGH_PRECISION = {
    "tol_ncdm_bg": 1.0e-10,
    "tol_ncdm_synchronous": 1.0e-10,
    "tol_ncdm_newtonian": 1.0e-10,
    "tol_perturbations_integration": 1.0e-7,
    "perturbations_sampling_stepsize": 0.01,
    "l_max_ncdm": 50,
    "q_linstep": 0.10,
    "q_logstep_spline": 30.0,
    "q_logstep_trapzd": 0.25,
    "q_numstep_transition": 400,
    "l_logstep": 1.015,
    "l_linstep": 15,
}


def parse_grid(text: str):
    vals = [float(v.strip()) for v in text.split(",") if v.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return vals


def class_params(psd: Path, mass: float, zmax: float, kmax_h: float, high_precision=False):
    p = hd.class_params(psd, mass, zmax, kmax_h)
    if high_precision:
        p.update(HIGH_PRECISION)
    return p


def interp_from_transfer(t, key, kh):
    kin = np.asarray(t["k (h/Mpc)"], float)
    return np.interp(kh, kin, np.asarray(t[key], float))


def build_state(psd: Path, mass: float, zgrid, kh, high_precision=False):
    c = Class(); c.set(class_params(psd, mass, max(zgrid), max(kh), high_precision)); c.compute()
    try:
        out = {}
        theta_available = None
        transfer_keys = None
        for z in zgrid:
            z = float(z)
            # Density-derivative relative velocity used by the existing forecast.
            vrel = hd.relative_velocity_kms(c, z, kh, 0.004)
            pk = hd.pk_cb_h(c, z, kh)
            t = c.get_transfer(z=z, output_format="class")
            if transfer_keys is None:
                transfer_keys = sorted(t.keys())
                theta_available = ("t_ncdm[0]" in t and "t_cdm" in t)
            theta_rel = None
            if theta_available:
                theta_rel = interp_from_transfer(t, "t_ncdm[0]", kh) - interp_from_transfer(t, "t_cdm", kh)
            out[z] = {"vrel": vrel, "pk": pk, "theta_rel": theta_rel}
    finally:
        c.struct_cleanup(); c.empty()
    return out, bool(theta_available), transfer_keys


def stable_mask(x0, threshold=1e-4):
    xmax = max(float(np.max(np.abs(x0))), 1e-300)
    return np.abs(x0) > threshold * xmax


def operator_from_states(zgrid, kh, s0, plus_states, minus_states, probe_frac, channel):
    nd = len(plus_states)
    blocks, weights = [], []
    dlnk = np.gradient(np.log(kh))
    for z in zgrid:
        z = float(z)
        if channel == "density_proxy":
            x0 = s0[z]["vrel"] * s0[z]["pk"]
        elif channel == "theta_proxy":
            if s0[z]["theta_rel"] is None:
                raise RuntimeError("CLASS transfer output does not expose t_ncdm[0] and t_cdm")
            x0 = s0[z]["theta_rel"] * s0[z]["pk"]
        else:
            raise ValueError(channel)
        mask = stable_mask(x0)
        J = np.zeros((int(np.sum(mask)), nd), float)
        for j in range(nd):
            sp, sm = plus_states[j][z], minus_states[j][z]
            if channel == "density_proxy":
                xp = sp["vrel"] * sp["pk"]
                xm = sm["vrel"] * sm["pk"]
            else:
                xp = sp["theta_rel"] * sp["pk"]
                xm = sm["theta_rel"] * sm["pk"]
            J[:, j] = ((xp - xm) / (2.0 * probe_frac))[mask] / x0[mask]
        blocks.append(J)
        weights.append(kh[mask] ** 3 * dlnk[mask])
    return np.vstack(blocks), np.concatenate(weights)


def actual_rms(zgrid, kh, s0, sp, sm, channel):
    vals, ww = [], []
    dlnk = np.gradient(np.log(kh))
    for z in zgrid:
        z = float(z)
        if channel == "density_proxy":
            x0 = s0[z]["vrel"] * s0[z]["pk"]
            xp = sp[z]["vrel"] * sp[z]["pk"]
            xm = sm[z]["vrel"] * sm[z]["pk"]
        elif channel == "theta_proxy":
            if s0[z]["theta_rel"] is None:
                return float("nan")
            x0 = s0[z]["theta_rel"] * s0[z]["pk"]
            xp = sp[z]["theta_rel"] * sp[z]["pk"]
            xm = sm[z]["theta_rel"] * sm[z]["pk"]
        else:
            raise ValueError(channel)
        mask = stable_mask(x0)
        f = (xp[mask] - xm[mask]) / (2.0 * x0[mask])
        w = kh[mask] ** 3 * dlnk[mask]
        vals.append(f); ww.append(w)
    f = np.concatenate(vals); w = np.concatenate(ww)
    return float(np.sqrt(np.sum(w * f * f) / np.sum(w)))


def normalize_shape(coeff, shapes, f0):
    rel_shapes = shapes / np.maximum(f0[None, :], 1e-300)
    nrm = hd.normalization_for_coeff(coeff, rel_shapes)
    return nrm * (coeff @ shapes)


def shape_cosine(q, a, b):
    return hd.cosine_shapes(q, a, b)


def weighted_operator_cosine(J1, w1, J2, w2):
    if J1.shape != J2.shape or w1.shape != w2.shape:
        return float("nan")
    w = np.sqrt(np.maximum(0.5 * (w1 + w2), 0.0))[:, None]
    a = (w * J1).ravel(); b = (w * J2).ravel()
    return float(np.dot(a, b) / math.sqrt(max(np.dot(a, a) * np.dot(b, b), 1e-300)))


def eigengap(opt):
    ev = np.asarray(opt["quadratic_eigenvalues"], float)
    return float(ev[0] / max(ev[1], 1e-300)) if ev.size > 1 else float("inf")


def run_pair(out: Path, tag: str, q, f0, shape, frac, mass, zgrid, kh, high_precision=False):
    fp = f0 + frac * shape; fm = f0 - frac * shape
    if min(float(np.min(fp)), float(np.min(fm))) <= 0:
        raise RuntimeError(f"{tag}: pair violates positivity")
    pp = out / f"{tag}_plus.dat"; pm = out / f"{tag}_minus.dat"
    hd.write_psd(pp, q, fp); hd.write_psd(pm, q, fm)
    sp, ta, _ = build_state(pp, mass, zgrid, kh, high_precision=high_precision)
    sm, tb, _ = build_state(pm, mass, zgrid, kh, high_precision=high_precision)
    return fp, fm, sp, sm, bool(ta and tb)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="hidden_channel_forensics")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-fracs", type=parse_grid, default=parse_grid("0.03,0.10"))
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--z-grid", type=parse_grid, default=parse_grid("0.10,0.30,0.60,1.00,2.00"))
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=192)
    ap.add_argument("--samples", type=int, default=30000)
    ap.add_argument("--seed", type=int, default=20260916)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights_mom, _, _, shapes, _, _ = cro.kinetic_objects(q, args.mass, args.z_match)
    if shapes.shape[0] != CREF.size:
        raise RuntimeError(f"expected {CREF.size} null directions, got {shapes.shape[0]}")
    kh = np.geomspace(args.kmin, args.kmax, args.nk)

    p0 = out / "fd.dat"; hd.write_psd(p0, q, f0)
    s0, theta_available, transfer_keys = build_state(p0, args.mass, args.z_grid, kh, high_precision=False)
    if not theta_available:
        print("WARNING: direct theta transfer unavailable; theta diagnostics will be skipped", flush=True)

    operators = {}
    optimized = {}
    for pf in args.probe_fracs:
        plus_states, minus_states = [], []
        for j, shape in enumerate(shapes):
            fp = f0 + pf * shape; fm = f0 - pf * shape
            pp = out / f"probe_{pf:.3f}_{j:02d}_plus.dat"
            pm = out / f"probe_{pf:.3f}_{j:02d}_minus.dat"
            hd.write_psd(pp, q, fp); hd.write_psd(pm, q, fm)
            sp, _, _ = build_state(pp, args.mass, args.z_grid, kh, high_precision=False)
            sm, _, _ = build_state(pm, args.mass, args.z_grid, kh, high_precision=False)
            plus_states.append(sp); minus_states.append(sm)
            print(f"PROBE frac={pf:.3f} dir={j+1}/{shapes.shape[0]}", flush=True)

        Jd, wd = operator_from_states(args.z_grid, kh, s0, plus_states, minus_states, pf, "density_proxy")
        od = hd.optimize_direction(Jd, wd, shapes, f0, args.final_frac, args.samples, args.seed + int(1000*pf))
        operators[(pf, "density_proxy")] = (Jd, wd)
        optimized[(pf, "density_proxy")] = od
        if theta_available:
            Jt, wt = operator_from_states(args.z_grid, kh, s0, plus_states, minus_states, pf, "theta_proxy")
            ot = hd.optimize_direction(Jt, wt, shapes, f0, args.final_frac, args.samples, args.seed + 10000 + int(1000*pf))
            operators[(pf, "theta_proxy")] = (Jt, wt)
            optimized[(pf, "theta_proxy")] = ot

    # Direct final-cap validation of every optimized direction.
    validation = []
    for (pf, channel), opt in optimized.items():
        shape = np.asarray(opt["shape"], float)
        _, _, sp, sm, theta_ok = run_pair(
            out, f"opt_{channel}_{pf:.3f}_final", q, f0, shape, args.final_frac,
            args.mass, args.z_grid, kh, high_precision=False,
        )
        rd = actual_rms(args.z_grid, kh, s0, sp, sm, "density_proxy")
        rt = actual_rms(args.z_grid, kh, s0, sp, sm, "theta_proxy") if theta_ok and theta_available else float("nan")
        validation.append({
            "probe_fraction": float(pf),
            "optimized_channel": channel,
            "predicted_final_rms": float(opt["predicted_rms"]),
            "actual_final_density_proxy_rms": rd,
            "actual_final_theta_proxy_rms": rt,
            "actual_over_prediction_density": rd / max(float(opt["predicted_rms"]), 1e-300),
            "eigengap_lambda1_over_lambda2": eigengap(opt),
        })
        print(f"VALIDATE probe={pf:.3f} channel={channel}", flush=True)

    # CREF amplitude sweep: direct physical scaling, independent of operator extrapolation.
    cref = CREF / np.linalg.norm(CREF)
    cref_shape = normalize_shape(cref, shapes, f0)
    cref_sweep = []
    for frac in sorted(set([float(x) for x in args.probe_fracs] + [args.final_frac])):
        fp, fm, sp, sm, theta_ok = run_pair(
            out, f"cref_{frac:.3f}", q, f0, cref_shape, frac,
            args.mass, args.z_grid, kh, high_precision=False,
        )
        cref_sweep.append({
            "fractional_cap": frac,
            "density_proxy_rms": actual_rms(args.z_grid, kh, s0, sp, sm, "density_proxy"),
            "theta_proxy_rms": actual_rms(args.z_grid, kh, s0, sp, sm, "theta_proxy") if theta_ok and theta_available else float("nan"),
            "max_relative_moment_mismatch": cro.pair_stats(f0, fp, fm, q, weights_mom)["max_relative_moment_mismatch"],
        })

    # Operator and optimized-shape consistency across probe amplitudes.
    consistency = []
    pfs = list(args.probe_fracs)
    for i in range(len(pfs)):
        for j in range(i+1, len(pfs)):
            p1, p2 = pfs[i], pfs[j]
            for channel in ("density_proxy", "theta_proxy"):
                if (p1, channel) not in operators or (p2, channel) not in operators:
                    continue
                J1, w1 = operators[(p1, channel)]; J2, w2 = operators[(p2, channel)]
                s1 = np.asarray(optimized[(p1, channel)]["shape"], float)
                s2 = np.asarray(optimized[(p2, channel)]["shape"], float)
                consistency.append({
                    "channel": channel,
                    "probe_fraction_1": float(p1),
                    "probe_fraction_2": float(p2),
                    "weighted_operator_cosine": weighted_operator_cosine(J1, w1, J2, w2),
                    "optimized_shape_cosine": shape_cosine(q, s1, s2),
                    "predicted_rms_1": float(optimized[(p1, channel)]["predicted_rms"]),
                    "predicted_rms_2": float(optimized[(p2, channel)]["predicted_rms"]),
                })

    # High-precision decisive check at final cap for CREF and the density-proxy
    # direction from the largest probe amplitude (least vulnerable to a noise floor).
    pf_ref = max(args.probe_fracs)
    selected_shape = np.asarray(optimized[(pf_ref, "density_proxy")]["shape"], float)
    p0hp = out / "fd_high_precision.dat"; hd.write_psd(p0hp, q, f0)
    s0hp, theta_hp, _ = build_state(p0hp, args.mass, args.z_grid, kh, high_precision=True)
    high_precision = []
    for tag, shape in (("cref", cref_shape), ("selected_density_opt", selected_shape)):
        fp, fm, sp, sm, tok = run_pair(
            out, f"hp_{tag}", q, f0, shape, args.final_frac,
            args.mass, args.z_grid, kh, high_precision=True,
        )
        high_precision.append({
            "direction": tag,
            "density_proxy_rms": actual_rms(args.z_grid, kh, s0hp, sp, sm, "density_proxy"),
            "theta_proxy_rms": actual_rms(args.z_grid, kh, s0hp, sp, sm, "theta_proxy") if tok and theta_hp else float("nan"),
            "max_relative_moment_mismatch": cro.pair_stats(f0, fp, fm, q, weights_mom)["max_relative_moment_mismatch"],
        })
        print(f"HIGH_PRECISION {tag}", flush=True)

    # Same-direction wake comparison using validated final-cap responses.
    wake_cref = hd.wake_half_pair_fraction(q, f0, f0 + args.final_frac*cref_shape,
                                           f0 - args.final_frac*cref_shape,
                                           args.mass, 0.30, 200.0)
    wake_selected = hd.wake_half_pair_fraction(q, f0, f0 + args.final_frac*selected_shape,
                                               f0 - args.final_frac*selected_shape,
                                               args.mass, 0.30, 200.0)

    def dump_csv(path, rows):
        if not rows: return
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    dump_csv(out / "final_validation.csv", validation)
    dump_csv(out / "cref_amplitude_sweep.csv", cref_sweep)
    dump_csv(out / "operator_consistency.csv", consistency)
    dump_csv(out / "high_precision_validation.csv", high_precision)

    hp_map = {r["direction"]: r for r in high_precision}
    summary = {
        "calculation": "cross-validated hidden-channel operator forensics",
        "status": "diagnostic_not_survey_forecast",
        "mass_eV": args.mass,
        "class_commit_expected": cro.CLASS_COMMIT,
        "probe_fractions": args.probe_fracs,
        "final_fractional_cap": args.final_frac,
        "null_dimension": int(shapes.shape[0]),
        "theta_transfer_available": theta_available,
        "transfer_keys": transfer_keys,
        "operator_consistency": consistency,
        "final_validation": validation,
        "cref_amplitude_sweep": cref_sweep,
        "high_precision_validation": high_precision,
        "wake_half_pair_fraction": {
            "cref": wake_cref,
            "selected_density_opt": wake_selected,
        },
        "high_precision_wake_to_linear_ratio": {
            "cref_over_density_proxy": wake_cref / max(hp_map["cref"]["density_proxy_rms"], 1e-300),
            "selected_over_density_proxy": wake_selected / max(hp_map["selected_density_opt"]["density_proxy_rms"], 1e-300),
            "cref_over_theta_proxy": None if not theta_hp else wake_cref / max(hp_map["cref"]["theta_proxy_rms"], 1e-300),
            "selected_over_theta_proxy": None if not theta_hp else wake_selected / max(hp_map["selected_density_opt"]["theta_proxy_rms"], 1e-300),
        },
        "interpretation_guide": {
            "operator_is_physical": "response operators and optimized shapes remain stable across probe amplitudes, direct final-cap validation follows prediction, and standard/high-precision plus density-derivative/theta estimators agree",
            "operator_overfit_or_noise": "optimized prediction collapses out of sample, operator/shape changes with probe amplitude, or high-precision/theta diagnostics disagree",
            "robust_structural_suppression": "direct final-cap and high-precision responses remain below 0.01 for CREF and independently optimized directions, with stable direct-theta confirmation",
        },
        "guardrail": "All statements are restricted to the tested seven-dimensional smooth source-matched class."
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("HIDDEN_CHANNEL_FORENSICS")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out/'summary.json'}")


if __name__ == "__main__":
    main()
