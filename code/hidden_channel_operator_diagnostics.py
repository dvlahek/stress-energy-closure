#!/usr/bin/env python3
r"""Deep diagnostics for hidden kinetic information in linear relative-velocity channels.

This script tests if the very small F_+ / F_- response found for the parity-odd
relative-velocity proxy is a property of one wake-optimized deformation or a
more structural property of the tested smooth stress-energy-null class.

It performs five checks at fixed relic mass:

1. Build the seven-dimensional smooth null space that exactly matches n, rho
   and P at z_match.
2. Probe every null direction with small symmetric CLASS runs and construct the
   linear response operators for v_(nu-cdm) and v_(nu-cdm) P_cb.
3. Optimize the hidden direction for each operator under the same 30 percent
   pointwise deformation cap used in the paper.
4. Validate the best proxy direction nonlinearly at several amplitudes and test
   finite-difference, k-grid and small-denominator-mask stability.
5. Compare the optimized linear-transfer response with the local resonant wake
   response for the same deformation.

The central question is therefore not just "does CREF give a small signal?" but
"can any smooth source-matched direction in the tested class evade that
suppression?"

This is a diagnostic calculation, not a survey Fisher forecast and not a
manuscript result by itself.
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
import parity_bispectrum_hidden_screening as pb

C_KMS = 299792.458
NS = 0.9649
KPIV_MPC = 0.05
TNU0_EV = cro.T_NCDM * cro.TCMB_K * cro.KB_EV_K

CREF = np.array([
    0.10329109892348824,
   -0.060769580477018574,
    0.16915381451528330,
    0.48937980690857140,
    0.68612921559544580,
    0.47909899899373720,
    0.13123736993120932,
], dtype=float)


def parse_grid(text: str):
    vals = [float(v.strip()) for v in text.split(",") if v.strip()]
    if not vals:
        raise argparse.ArgumentTypeError("grid must contain at least one value")
    return vals


def write_psd(path: Path, q, f):
    if np.min(f) <= 0:
        raise RuntimeError(f"non-positive PSD in {path}: {np.min(f)}")
    np.savetxt(path, np.column_stack([q, f]), fmt="%.14e")


def class_params(psd: Path, mass: float, zmax: float, kmax_h: float):
    return {
        "output": "mPk,dTk",
        "modes": "s",
        "gauge": "newtonian",
        "H0": cro.H0,
        "omega_b": cro.OMEGA_B,
        "omega_cdm": cro.OMEGA_CDM,
        "A_s": cro.A_S,
        "n_s": NS,
        "tau_reio": cro.TAU_REIO,
        "N_ur": cro.N_UR,
        "N_ncdm": 1,
        "use_ncdm_psd_files": 1,
        "ncdm_psd_filenames": str(psd.resolve()),
        "m_ncdm": mass,
        "T_ncdm": cro.T_NCDM,
        "deg_ncdm": 1.0,
        "P_k_max_h/Mpc": max(0.35, 1.1 * kmax_h),
        "z_max_pk": max(1.0, zmax + 0.25),
    }


def interp_transfer(cosmo, z, key, kh):
    t = cosmo.get_transfer(z=float(z), output_format="class")
    kin = np.asarray(t["k (h/Mpc)"], float)
    val = np.asarray(t[key], float)
    return np.interp(kh, kin, val)


def delta_rel(cosmo, z, kh):
    return (
        interp_transfer(cosmo, z, "d_ncdm[0]", kh)
        - interp_transfer(cosmo, z, "d_cdm", kh)
    )


def relative_velocity_kms(cosmo, z, kh, dz_factor):
    dz = dz_factor * (1.0 + z)
    zm = max(0.0, z - dz)
    zp = z + dz
    if zm == 0.0 and z < dz:
        d0 = delta_rel(cosmo, z, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - d0) / (zp - z)
    else:
        dm = delta_rel(cosmo, zm, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - dm) / (zp - zm)
    h = float(cosmo.h())
    k_mpc = kh * h
    H_mpc = float(cosmo.Hubble(float(z)))
    primordial = np.sqrt(cro.A_S * (k_mpc / KPIV_MPC) ** (NS - 1.0))
    return -H_mpc * d_dz / np.maximum(k_mpc, 1e-30) * primordial * C_KMS


def pk_cb_h(cosmo, z, kh):
    h = float(cosmo.h())
    return np.array([cosmo.pk_cb_lin(float(k * h), float(z)) * h**3 for k in kh])


def build_state(psd: Path, mass: float, zgrid, kh, dz_factors=(0.004,)):
    c = Class()
    c.set(class_params(psd, mass, max(zgrid), max(kh)))
    c.compute()
    try:
        out = {}
        for z in zgrid:
            z = float(z)
            row = {"pk": pk_cb_h(c, z, kh), "vrel": {}}
            for dzf in dz_factors:
                row["vrel"][float(dzf)] = relative_velocity_kms(c, z, kh, float(dzf))
            out[z] = row
    finally:
        c.struct_cleanup()
        c.empty()
    return out


def pair_stats(q, f0, fp, fm, weights):
    return cro.pair_stats(f0, fp, fm, q, weights)


def stable_mask(x0, threshold):
    xmax = max(float(np.max(np.abs(x0))), 1e-300)
    return np.abs(x0) > threshold * xmax


def response_operator(zgrid, kh, s0, plus_states, minus_states, probe_frac,
                      field="proxy", dz_factor=0.004, threshold=1e-4):
    nd = len(plus_states)
    blocks = []
    weights = []
    metadata = []
    dlnk = np.gradient(np.log(kh))
    for z in zgrid:
        z = float(z)
        v0 = s0[z]["vrel"][float(dz_factor)]
        p0 = s0[z]["pk"]
        if field == "velocity":
            x0 = v0
        elif field == "proxy":
            x0 = v0 * p0
        else:
            raise ValueError(field)
        mask = stable_mask(x0, threshold)
        if np.sum(mask) < 5:
            raise RuntimeError(f"too few stable modes for {field} at z={z}")
        J = np.zeros((int(np.sum(mask)), nd), float)
        for j in range(nd):
            vp = plus_states[j][z]["vrel"][float(dz_factor)]
            vm = minus_states[j][z]["vrel"][float(dz_factor)]
            if field == "velocity":
                xp, xm = vp, vm
            else:
                xp = vp * plus_states[j][z]["pk"]
                xm = vm * minus_states[j][z]["pk"]
            deriv = (xp - xm) / (2.0 * probe_frac)
            # Signed fractional derivative. The sign drops out of the quadratic
            # optimization but is retained for reconstruction and diagnostics.
            J[:, j] = deriv[mask] / x0[mask]
        w = kh[mask] ** 3 * dlnk[mask]
        blocks.append(J)
        weights.append(w)
        metadata.extend([(z, float(k)) for k in kh[mask]])
    return np.vstack(blocks), np.concatenate(weights), metadata


def normalization_for_coeff(c, rel_shapes):
    raw_rel = c @ rel_shapes
    maxrel = float(np.max(np.abs(raw_rel)))
    if not np.isfinite(maxrel) or maxrel <= 0:
        return None
    return 1.0 / maxrel


def quadratic_matrix(J, w):
    return J.T @ (w[:, None] * J), float(np.sum(w))


def optimize_direction(J, w, shapes, f0, final_frac, samples, seed):
    rel_shapes = shapes / np.maximum(f0[None, :], 1e-300)
    P, wsum = quadratic_matrix(J, w)
    evals, evecs = np.linalg.eigh(0.5 * (P + P.T))
    candidates = [evecs[:, j] for j in np.argsort(evals)[::-1]]
    nd = shapes.shape[0]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0
        candidates.extend([e, -e])
    candidates.append(CREF / np.linalg.norm(CREF))

    rng = np.random.default_rng(seed)
    rr = rng.normal(size=(samples, nd))
    rr /= np.linalg.norm(rr, axis=1, keepdims=True)

    best = None

    def evaluate(c):
        nonlocal best
        nrm = normalization_for_coeff(c, rel_shapes)
        if nrm is None:
            return
        qv = float(c @ P @ c)
        rms = final_frac * nrm * math.sqrt(max(qv, 0.0) / max(wsum, 1e-300))
        if best is None or rms > best[0]:
            best = (rms, c.copy(), nrm)

    for c in candidates:
        evaluate(c)
    # Random search in chunks so the exact pointwise cap can be evaluated
    # without constructing a very large samples x q array.
    chunk = 2000
    for start in range(0, samples, chunk):
        block = rr[start:start + chunk]
        rel = block @ rel_shapes
        maxrel = np.max(np.abs(rel), axis=1)
        good = np.isfinite(maxrel) & (maxrel > 0)
        if not np.any(good):
            continue
        bg = block[good]
        ng = 1.0 / maxrel[good]
        qv = np.einsum("bi,ij,bj->b", bg, P, bg)
        rms = final_frac * ng * np.sqrt(np.maximum(qv, 0.0) / max(wsum, 1e-300))
        idx = int(np.argmax(rms))
        if best is None or float(rms[idx]) > best[0]:
            best = (float(rms[idx]), bg[idx].copy(), float(ng[idx]))

    rms, coeff, nrm = best
    shape = nrm * (coeff @ shapes)
    return {
        "predicted_rms": rms,
        "coefficients": coeff.tolist(),
        "pointwise_normalization": nrm,
        "shape": shape,
        "quadratic_eigenvalues": np.sort(evals)[::-1].tolist(),
    }


def rms_for_coeff(J, w, c, shapes, f0, frac):
    rel_shapes = shapes / np.maximum(f0[None, :], 1e-300)
    nrm = normalization_for_coeff(c, rel_shapes)
    pred = frac * nrm * (J @ c)
    return float(np.sqrt(np.sum(w * pred * pred) / np.sum(w)))


def actual_fractional_rms(zgrid, kh, s0, sp, sm, field="proxy",
                          dz_factor=0.004, threshold=1e-4, stride=1):
    vals = []
    ww = []
    dlnk = np.gradient(np.log(kh))
    ids = np.arange(0, len(kh), stride)
    for z in zgrid:
        z = float(z)
        v0 = s0[z]["vrel"][float(dz_factor)]
        vp = sp[z]["vrel"][float(dz_factor)]
        vm = sm[z]["vrel"][float(dz_factor)]
        if field == "velocity":
            x0, xp, xm = v0, vp, vm
        else:
            x0 = v0 * s0[z]["pk"]
            xp = vp * sp[z]["pk"]
            xm = vm * sm[z]["pk"]
        mask = stable_mask(x0, threshold)
        keep = mask.copy()
        index_mask = np.zeros_like(mask, dtype=bool)
        index_mask[ids] = True
        keep &= index_mask
        frac = np.zeros_like(x0)
        frac[keep] = (xp[keep] - xm[keep]) / (2.0 * x0[keep])
        w = kh ** 3 * dlnk
        vals.append(frac[keep])
        ww.append(w[keep])
    v = np.concatenate(vals)
    w = np.concatenate(ww)
    return float(np.sqrt(np.sum(w * v * v) / np.sum(w)))


def wake_half_pair_fraction(q, f0, fp, fm, mass, z, sigma_kms, quadrature=128):
    x0, w0 = np.polynomial.hermite.hermgauss(quadrature)
    x = np.sqrt(2.0) * x0
    w = w0 / np.sqrt(np.pi)
    qres = mass * np.abs(x) * sigma_kms / (C_KMS * TNU0_EV * (1.0 + z))
    def amp(f):
        F = np.interp(qres, q, f, left=f[0], right=f[-1])
        return float(np.sum(w * np.abs(x) * F))
    a0, ap, am = amp(f0), amp(fp), amp(fm)
    return 0.5 * abs(ap - am) / max(abs(a0), 1e-300)


def cosine_shapes(q, a, b):
    aa = float(np.trapezoid(a * a, q))
    bb = float(np.trapezoid(b * b, q))
    ab = float(np.trapezoid(a * b, q))
    return ab / math.sqrt(max(aa * bb, 1e-300))


def save_row_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader(); writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="hidden_channel_operator_diagnostics")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.03)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--z-grid", type=parse_grid,
                    default=parse_grid("0.10,0.30,0.60,1.00,2.00"))
    ap.add_argument("--kmin", type=float, default=0.003)
    ap.add_argument("--kmax", type=float, default=0.20)
    ap.add_argument("--nk", type=int, default=192)
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=20260915)
    ap.add_argument("--validation-fracs", type=parse_grid,
                    default=parse_grid("0.10,0.20,0.30"))
    ap.add_argument("--dz-factors", type=parse_grid,
                    default=parse_grid("0.002,0.004,0.008"))
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights_mom, _, _, shapes, _, _ = cro.kinetic_objects(q, args.mass, args.z_match)
    nd = shapes.shape[0]
    if nd != CREF.size:
        raise RuntimeError(f"expected {CREF.size} null directions, got {nd}")
    kh = np.geomspace(args.kmin, args.kmax, args.nk)

    p0 = out / "fd.dat"; write_psd(p0, q, f0)
    s0 = build_state(p0, args.mass, args.z_grid, kh, dz_factors=args.dz_factors)

    plus_states, minus_states = [], []
    for j, shape in enumerate(shapes):
        fp = f0 + args.probe_frac * shape
        fm = f0 - args.probe_frac * shape
        if min(np.min(fp), np.min(fm)) <= 0:
            raise RuntimeError(f"probe direction {j} violates positivity")
        pp = out / f"probe_{j:02d}_plus.dat"
        pm = out / f"probe_{j:02d}_minus.dat"
        write_psd(pp, q, fp); write_psd(pm, q, fm)
        plus_states.append(build_state(pp, args.mass, args.z_grid, kh, dz_factors=(0.004,)))
        minus_states.append(build_state(pm, args.mass, args.z_grid, kh, dz_factors=(0.004,)))
        print(f"PROBE_DONE {j+1}/{nd}", flush=True)

    Jv, wv, _ = response_operator(args.z_grid, kh, s0, plus_states, minus_states,
                                   args.probe_frac, field="velocity")
    Ja, wa, _ = response_operator(args.z_grid, kh, s0, plus_states, minus_states,
                                   args.probe_frac, field="proxy")

    opt_v = optimize_direction(Jv, wv, shapes, f0, args.final_frac, args.samples, args.seed)
    opt_a = optimize_direction(Ja, wa, shapes, f0, args.final_frac, args.samples, args.seed + 1)
    cref = CREF / np.linalg.norm(CREF)
    cref_v = rms_for_coeff(Jv, wv, cref, shapes, f0, args.final_frac)
    cref_a = rms_for_coeff(Ja, wa, cref, shapes, f0, args.final_frac)

    # Nonlinear validation of the most favorable proxy direction.
    best_shape = np.asarray(opt_a.pop("shape"), float)
    cref_shape = normalization_for_coeff(cref, shapes / np.maximum(f0[None, :], 1e-300)) * (cref @ shapes)
    validation = []
    final_states = None
    for frac in args.validation_fracs:
        fp = f0 + frac * best_shape
        fm = f0 - frac * best_shape
        pp = out / f"best_proxy_plus_{frac:.3f}.dat"
        pm = out / f"best_proxy_minus_{frac:.3f}.dat"
        write_psd(pp, q, fp); write_psd(pm, q, fm)
        dzs = args.dz_factors if abs(frac - args.final_frac) < 1e-12 else (0.004,)
        sp = build_state(pp, args.mass, args.z_grid, kh, dz_factors=dzs)
        sm = build_state(pm, args.mass, args.z_grid, kh, dz_factors=dzs)
        rv = actual_fractional_rms(args.z_grid, kh, s0, sp, sm, "velocity")
        ra = actual_fractional_rms(args.z_grid, kh, s0, sp, sm, "proxy")
        pred_a = float(opt_a["predicted_rms"]) * frac / args.final_frac
        validation.append({
            "fractional_cap": float(frac),
            "actual_velocity_rms": rv,
            "actual_proxy_rms": ra,
            "linear_predicted_proxy_rms": pred_a,
            "actual_over_prediction_proxy": ra / max(pred_a, 1e-300),
            "max_relative_moment_mismatch": pair_stats(q, f0, fp, fm, weights_mom)["max_relative_moment_mismatch"],
        })
        if abs(frac - args.final_frac) < 1e-12:
            final_states = (fp, fm, sp, sm)
        print(f"VALIDATION_DONE frac={frac}", flush=True)

    if final_states is None:
        # Ensure the final cap is validated even if the user supplied a custom list.
        frac = args.final_frac
        fp = f0 + frac * best_shape; fm = f0 - frac * best_shape
        pp = out / "best_proxy_plus_final.dat"; pm = out / "best_proxy_minus_final.dat"
        write_psd(pp, q, fp); write_psd(pm, q, fm)
        sp = build_state(pp, args.mass, args.z_grid, kh, dz_factors=args.dz_factors)
        sm = build_state(pm, args.mass, args.z_grid, kh, dz_factors=args.dz_factors)
        final_states = (fp, fm, sp, sm)

    fp_best, fm_best, sp_best, sm_best = final_states

    # Numerical stability checks on the nonlinear best direction.
    dz_checks = []
    for dzf in args.dz_factors:
        rv = actual_fractional_rms(args.z_grid, kh, s0, sp_best, sm_best,
                                   "velocity", dz_factor=dzf)
        ra = actual_fractional_rms(args.z_grid, kh, s0, sp_best, sm_best,
                                   "proxy", dz_factor=dzf)
        dz_checks.append({"dz_factor": float(dzf), "velocity_rms": rv, "proxy_rms": ra})

    k_checks = []
    for stride in (1, 2, 4):
        k_checks.append({
            "effective_nk": int(math.ceil(args.nk / stride)),
            "stride": stride,
            "proxy_rms": actual_fractional_rms(args.z_grid, kh, s0, sp_best, sm_best,
                                                "proxy", stride=stride),
        })

    mask_checks = []
    for thr in (1e-3, 1e-4, 1e-5):
        mask_checks.append({
            "relative_denominator_threshold": thr,
            "proxy_rms": actual_fractional_rms(args.z_grid, kh, s0, sp_best, sm_best,
                                                "proxy", threshold=thr),
        })

    # Standardize the local wake comparison to HALF-pair fractional response,
    # matching the transfer convention used above.
    cref_fp = f0 + args.final_frac * cref_shape
    cref_fm = f0 - args.final_frac * cref_shape
    wake_cref = wake_half_pair_fraction(q, f0, cref_fp, cref_fm, args.mass, 0.30, 200.0)
    wake_best = wake_half_pair_fraction(q, f0, fp_best, fm_best, args.mass, 0.30, 200.0)

    np.savetxt(out / "best_proxy_pair.csv",
               np.column_stack([q, f0, fp_best, fm_best, best_shape]), delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape", comments="")
    save_row_csv(out / "nonlinear_validation.csv", validation)
    save_row_csv(out / "dz_stability.csv", dz_checks)
    save_row_csv(out / "kgrid_stability.csv", k_checks)
    save_row_csv(out / "mask_stability.csv", mask_checks)

    summary = {
        "calculation": "deep hidden-channel response-operator diagnostics",
        "status": "diagnostic_not_survey_forecast",
        "class_commit_expected": cro.CLASS_COMMIT,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "null_dimension": nd,
        "probe_fraction": args.probe_frac,
        "final_fractional_cap": args.final_frac,
        "k_range_h_Mpc": [args.kmin, args.kmax],
        "nk": args.nk,
        "z_grid": args.z_grid,
        "cref_predicted_velocity_rms": cref_v,
        "cref_predicted_proxy_rms": cref_a,
        "optimized_velocity": {k: v for k, v in opt_v.items() if k != "shape"},
        "optimized_proxy": opt_a,
        "optimized_proxy_gain_over_cref": float(opt_a["predicted_rms"] / max(cref_a, 1e-300)),
        "shape_cosine_optimized_proxy_vs_cref": cosine_shapes(q, best_shape, cref_shape),
        "nonlinear_validation": validation,
        "dz_stability": dz_checks,
        "kgrid_stability": k_checks,
        "mask_stability": mask_checks,
        "wake_half_pair_fraction_at_z0p3_sigma200": {
            "cref": wake_cref,
            "optimized_proxy_direction": wake_best,
        },
        "wake_to_linear_proxy_ratio_same_convention": {
            "cref": float(wake_cref / max(cref_a, 1e-300)),
            "optimized_proxy_direction": float(wake_best / max(opt_a["predicted_rms"], 1e-300)),
        },
        "decision_guide": {
            "structural_suppression": "optimized nonlinear proxy RMS remains <0.01 and numerical/linearity checks are stable",
            "direction_specific": "optimization raises proxy RMS by orders of magnitude or above 0.03",
            "numerical_warning": "result shifts materially with dz factor, k grid, mask threshold, or deformation amplitude",
        },
        "analysis_scope": "The optimization covers the tested ten-function smooth deformation class with three matched source moments.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print("HIDDEN_CHANNEL_OPERATOR_DIAGNOSTICS")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote {out / 'summary.json'}")


if __name__ == "__main__":
    main()
