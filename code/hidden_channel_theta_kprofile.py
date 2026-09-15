#!/usr/bin/env python3
r"""k-resolved decomposition of direct-theta hidden-state retention.

This is the follow-up to hidden_channel_theta_pair_check.py. It evaluates a
source-matched pair at one redshift but writes the fractional response as a
function of k instead of only an RMS summary. The goal is to distinguish the
pure velocity-divergence response from the density-weighted parity-odd proxy and
to verify exactly how P_cb changes the response.

For each k the code reports

  delta_theta      = (theta_+ - theta_-)/(2 theta_0)
  delta_P          = (P_+ - P_-)/(2 P_0)
  bar_theta        = (theta_+ + theta_-)/(2 theta_0)
  bar_P            = (P_+ + P_-)/(2 P_0)
  delta_thetaP     = (theta_+ P_+ - theta_- P_-)/(2 theta_0 P_0)

with the exact identity

  delta_thetaP = bar_theta * delta_P + bar_P * delta_theta.

The density-derived velocity-times-P_cb proxy is written on the same grid as an
independent control. The ``orthogonal`` direction is a deterministic,
response-independent coefficient-space control: the first canonical null-basis
axis is Gram-Schmidt orthogonalized against CREF, then pointwise normalized and
used with the same deformation cap. No CLASS response amplitude enters its
definition. No survey S/N is inferred.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import class_response_optimize as cro
import hidden_channel_operator_diagnostics as hd
import hidden_channel_theta_pair_check as ht


def safe_frac(num, den, threshold=1e-4):
    den = np.asarray(den, float)
    num = np.asarray(num, float)
    scale = max(float(np.max(np.abs(den))), 1e-300)
    good = np.abs(den) > threshold * scale
    out = np.full_like(den, np.nan, dtype=float)
    out[good] = num[good] / den[good]
    return out, good


def weighted_rms(frac, kh, mask):
    dlnk = np.gradient(np.log(kh))
    w = kh**3 * dlnk
    m = mask & np.isfinite(frac)
    if np.sum(m) < 2:
        return float("nan")
    return float(np.sqrt(np.sum(w[m] * frac[m]**2) / np.sum(w[m])))


def zero_crossings(x, mask):
    ids = np.flatnonzero(mask & np.isfinite(x))
    if ids.size < 2:
        return 0
    y = x[ids]
    s = np.sign(y)
    nz = s != 0
    s = s[nz]
    if s.size < 2:
        return 0
    return int(np.sum(s[1:] * s[:-1] < 0))


def reconstruct_orthogonal_control(mass: float, z_match: float, frac: float):
    """Build a deterministic response-independent null direction.

    Start from the first canonical coefficient-space basis vector, remove its
    projection onto normalized CREF, normalize the remaining coefficient vector,
    then impose the same pointwise relative-shape normalization and +/- cap as
    the reference pair. This uses no CLASS transfer response or optimizer output.
    """
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, _, _, shapes, _, _ = cro.kinetic_objects(q, mass, z_match)
    cref = ht.CREF / np.linalg.norm(ht.CREF)
    e0 = np.zeros_like(cref)
    e0[0] = 1.0
    coeff = e0 - float(np.dot(e0, cref)) * cref
    ncoeff = float(np.linalg.norm(coeff))
    if not np.isfinite(ncoeff) or ncoeff <= 0:
        raise RuntimeError("failed to construct coefficient-space orthogonal control")
    coeff /= ncoeff

    rel_shapes = shapes / np.maximum(f0[None, :], 1e-300)
    nrm = hd.normalization_for_coeff(coeff, rel_shapes)
    if nrm is None:
        raise RuntimeError("failed to pointwise-normalize orthogonal control")
    shape = nrm * (coeff @ shapes)
    fp = f0 + frac * shape
    fm = f0 - frac * shape
    if min(float(np.min(fp)), float(np.min(fm))) <= 0:
        raise RuntimeError("orthogonal control violates positivity")

    metadata = {
        "definition": "Gram-Schmidt first canonical null-basis axis against normalized CREF; no response optimization",
        "coefficients": coeff.tolist(),
        "coefficient_dot_normalized_cref": float(np.dot(coeff, cref)),
        "pointwise_normalization": float(nrm),
    }
    return q, f0, fp, fm, weights, metadata


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="hidden_channel_theta_kprofiles")
    ap.add_argument("--direction", choices=("cref", "orthogonal", "selected"), required=True)
    ap.add_argument("--pair-csv", type=Path,
                    default=Path("hidden_channel_operator_diagnostics/best_proxy_pair.csv"))
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--z", type=float, default=0.30)
    ap.add_argument("--sigma-kms", type=float, default=200.0)
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=96)
    ap.add_argument("--mask-threshold", type=float, default=1e-4)
    ap.add_argument("--precision", choices=("standard", "moderate"), default="standard")
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    direction_definition = None
    if args.direction == "cref":
        q, f0, fp, fm, weights = ht.reconstruct_cref(args.mass, args.z_match, args.frac)
        direction_definition = {
            "definition": "fixed reference/CREF null direction",
            "coefficients": (ht.CREF / np.linalg.norm(ht.CREF)).tolist(),
        }
    elif args.direction == "orthogonal":
        q, f0, fp, fm, weights, direction_definition = reconstruct_orthogonal_control(
            args.mass, args.z_match, args.frac
        )
    else:
        if not args.pair_csv.exists():
            raise FileNotFoundError(f"Selected-pair file not found: {args.pair_csv}")
        q, f0, fp, fm, weights = ht.load_selected(args.pair_csv, args.mass, args.z_match)
        direction_definition = {
            "definition": "response-selected development direction; precision-sensitive and not retained for quantitative claims",
            "pair_csv": str(args.pair_csv),
        }

    tag = f"{args.direction}_z{args.z:.3f}_nk{args.nk}_{args.precision}"
    p0 = out / f"{tag}_fd.dat"
    pp = out / f"{tag}_plus.dat"
    pm = out / f"{tag}_minus.dat"
    ht.write_psd(p0, q, f0)
    ht.write_psd(pp, q, fp)
    ht.write_psd(pm, q, fm)

    kh = np.geomspace(args.kmin, args.kmax, args.nk)
    zgrid = [float(args.z)]

    print(f"STATE 1/3 FD direction={args.direction} z={args.z} precision={args.precision}", flush=True)
    s0 = ht.build_state_isolated(p0, args.mass, zgrid, kh, precision=args.precision)
    print("STATE 2/3 PLUS", flush=True)
    sp = ht.build_state_isolated(pp, args.mass, zgrid, kh, precision=args.precision)
    print("STATE 3/3 MINUS", flush=True)
    sm = ht.build_state_isolated(pm, args.mass, zgrid, kh, precision=args.precision)

    if not (s0["theta_available"] and sp["theta_available"] and sm["theta_available"]):
        raise RuntimeError("direct theta transfer unavailable")

    z = float(args.z)
    r0 = s0["state"][z]
    rp = sp["state"][z]
    rm = sm["state"][z]

    th0, thp, thm = r0["theta_rel"], rp["theta_rel"], rm["theta_rel"]
    P0, Pp, Pm = r0["pk"], rp["pk"], rm["pk"]
    v0, vp, vm = r0["vrel_density"], rp["vrel_density"], rm["vrel_density"]

    delta_theta, m_theta = safe_frac(thp - thm, 2.0 * th0, args.mask_threshold)
    delta_P, m_P = safe_frac(Pp - Pm, 2.0 * P0, args.mask_threshold)
    bar_theta, _ = safe_frac(thp + thm, 2.0 * th0, args.mask_threshold)
    bar_P, _ = safe_frac(Pp + Pm, 2.0 * P0, args.mask_threshold)
    delta_thetaP, m_thetaP = safe_frac(thp * Pp - thm * Pm, 2.0 * th0 * P0, args.mask_threshold)
    delta_vP, m_vP = safe_frac(vp * Pp - vm * Pm, 2.0 * v0 * P0, args.mask_threshold)

    term_P = bar_theta * delta_P
    term_theta = bar_P * delta_theta
    reconstructed = term_P + term_theta
    decomp_resid = delta_thetaP - reconstructed

    common = m_theta & m_P & m_thetaP & m_vP & np.isfinite(reconstructed)

    rows = []
    for i, k in enumerate(kh):
        rows.append({
            "k_h_Mpc": float(k),
            "valid_common": int(bool(common[i])),
            "delta_theta": float(delta_theta[i]) if np.isfinite(delta_theta[i]) else "",
            "delta_Pcb": float(delta_P[i]) if np.isfinite(delta_P[i]) else "",
            "bar_theta": float(bar_theta[i]) if np.isfinite(bar_theta[i]) else "",
            "bar_Pcb": float(bar_P[i]) if np.isfinite(bar_P[i]) else "",
            "delta_theta_Pcb": float(delta_thetaP[i]) if np.isfinite(delta_thetaP[i]) else "",
            "theta_bar_times_delta_Pcb": float(term_P[i]) if np.isfinite(term_P[i]) else "",
            "Pcb_bar_times_delta_theta": float(term_theta[i]) if np.isfinite(term_theta[i]) else "",
            "decomposition_residual": float(decomp_resid[i]) if np.isfinite(decomp_resid[i]) else "",
            "delta_density_derived_velocity_Pcb": float(delta_vP[i]) if np.isfinite(delta_vP[i]) else "",
        })

    csv_path = out / f"{tag}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    wake = hd.wake_half_pair_fraction(q, f0, fp, fm, args.mass, z, args.sigma_kms)

    edges = np.geomspace(args.kmin, args.kmax, 4)
    bins = []
    for j in range(3):
        if j < 2:
            mb = common & (kh >= edges[j]) & (kh < edges[j + 1])
        else:
            mb = common & (kh >= edges[j]) & (kh <= edges[j + 1])
        bin_theta = weighted_rms(delta_theta, kh, mb)
        bin_thetaP = weighted_rms(delta_thetaP, kh, mb)
        bin_vP = weighted_rms(delta_vP, kh, mb)
        bins.append({
            "k_range_h_Mpc": [float(edges[j]), float(edges[j + 1])],
            "n_modes": int(np.sum(mb)),
            "theta_rms": bin_theta,
            "theta_Pcb_rms": bin_thetaP,
            "density_derived_velocity_Pcb_rms": bin_vP,
            "wake_to_theta_Pcb_contrast_same_wake_numerator": wake / max(bin_thetaP, 1e-300),
            "wake_to_density_derived_proxy_contrast_same_wake_numerator": wake / max(bin_vP, 1e-300),
        })

    stats = cro.pair_stats(f0, fp, fm, q, weights)

    max_profile = max(float(np.nanmax(np.abs(delta_thetaP[common]))), 1e-300)
    max_resid = float(np.nanmax(np.abs(decomp_resid[common]))) if np.any(common) else float("nan")
    rms_theta = weighted_rms(delta_theta, kh, common)
    rms_thetaP = weighted_rms(delta_thetaP, kh, common)
    rms_vP = weighted_rms(delta_vP, kh, common)

    summary = {
        "calculation": "k-resolved direct-theta retention decomposition",
        "direction": args.direction,
        "direction_definition": direction_definition,
        "precision": args.precision,
        "mass_eV": args.mass,
        "z": z,
        "k_range_h_Mpc": [args.kmin, args.kmax],
        "nk": args.nk,
        "mask_threshold": args.mask_threshold,
        "class_output_requested": "mPk,dTk,vTk",
        "cdm_velocity_key": s0["cdm_velocity_key"],
        "ncdm_velocity_key": s0["ncdm_velocity_key"],
        "valid_common_modes": int(np.sum(common)),
        "direct_theta_half_pair_rms": rms_theta,
        "direct_theta_Pcb_half_pair_rms": rms_thetaP,
        "density_derived_velocity_Pcb_half_pair_rms": rms_vP,
        "wake_half_pair_fraction_sigma_kms": wake,
        "wake_to_direct_theta_ratio": wake / max(rms_theta, 1e-300),
        "wake_to_direct_theta_Pcb_ratio": wake / max(rms_thetaP, 1e-300),
        "wake_to_density_derived_proxy_ratio": wake / max(rms_vP, 1e-300),
        "theta_zero_crossings": zero_crossings(delta_theta, common),
        "theta_Pcb_zero_crossings": zero_crossings(delta_thetaP, common),
        "density_derived_proxy_zero_crossings": zero_crossings(delta_vP, common),
        "max_abs_decomposition_residual": max_resid,
        "max_abs_decomposition_residual_relative_to_thetaP_peak": max_resid / max_profile,
        "log_k_thirds": bins,
        "max_relative_moment_mismatch": stats["max_relative_moment_mismatch"],
        "identity": "delta(theta P)=bar(theta) delta(P)+bar(P) delta(theta)",
        "interpretation_guardrail": (
            "Pure theta and theta*P_cb are different transfer-level quantities; their wake ratios "
            "need not agree. The wake fraction used here has no k dependence after common potential "
            "factors cancel, so bin-specific contrasts reuse the same wake numerator and are not "
            "independent wake observables. This diagnostic is not a kSZ or RSD survey forecast."
        ),
        "profile_csv": str(csv_path),
    }

    js = out / f"{tag}.json"
    js.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("HIDDEN_CHANNEL_THETA_KPROFILE")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {csv_path}")
    print(f"Wrote {js}")


if __name__ == "__main__":
    main()
