#!/usr/bin/env python3
"""Ideal redshift-space forecast for stress-energy-matched relic distributions.

The calculation searches the exact n-rho-P null space already used by the
manuscript and maximizes the part of the linear Kaiser redshift-space response
that remains after projection over standard cosmological parameters and one
independent galaxy-bias amplitude in each redshift bin.

This is deliberately an optimistic diagnostic. It is not a survey likelihood.
A positive result identifies a parameter-orthogonal observational channel worth
promoting to a survey-specific forecast. A null result is equally useful because
it rules out linear RSD as an easy escape from the degeneracy found in P(k).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import class_response_optimize as cro

try:
    from classy import Class
except Exception as exc:
    raise RuntimeError("classy is required. Install the pinned CLASS checkout with `pip install ./class_public`.") from exc

NS = 0.9649
Z_BINS = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0])
K_H = np.geomspace(0.01, 0.30, 36)
MU, MU_W = np.polynomial.legendre.leggauss(12)
BIAS = 1.0 + 0.84 * Z_BINS


def write_psd(path, q, f):
    if np.min(f) <= 0:
        raise RuntimeError(f"Non-positive PSD in {path}")
    np.savetxt(path, np.column_stack([q, f]), fmt="%.14e")


def class_parameters(psd_path, mass_eV, overrides=None):
    p = {
        "output": "mPk",
        "modes": "s",
        "H0": cro.H0,
        "omega_b": cro.OMEGA_B,
        "omega_cdm": cro.OMEGA_CDM,
        "A_s": cro.A_S,
        "n_s": NS,
        "tau_reio": cro.TAU_REIO,
        "N_ur": cro.N_UR,
        "N_ncdm": 1,
        "use_ncdm_psd_files": 1,
        "ncdm_psd_filenames": str(Path(psd_path).resolve()),
        "m_ncdm": mass_eV,
        "T_ncdm": cro.T_NCDM,
        "deg_ncdm": 1.0,
        "P_k_max_h/Mpc": 0.6,
        "z_max_pk": 3.2,
    }
    if overrides:
        p.update(overrides)
    return p


def pk_cb_grid(cosmo, z, k_h):
    h = float(cosmo.h())
    return np.array([cosmo.pk_cb_lin(float(k * h), float(z)) for k in k_h])


def growth_grid(cosmo, z, k_h):
    p0 = pk_cb_grid(cosmo, z, k_h)
    dz = 0.01 * (1.0 + z)
    if z > dz:
        pm = pk_cb_grid(cosmo, z - dz, k_h)
        pp = pk_cb_grid(cosmo, z + dz, k_h)
        dlnp_dz = (np.log(pp) - np.log(pm)) / (2.0 * dz)
    else:
        pp = pk_cb_grid(cosmo, z + dz, k_h)
        dlnp_dz = (np.log(pp) - np.log(p0)) / dz
    f = -0.5 * (1.0 + z) * dlnp_dz
    return p0, f


def observable_from_psd(psd_path, mass_eV, overrides=None):
    cosmo = Class()
    cosmo.set(class_parameters(psd_path, mass_eV, overrides))
    cosmo.compute()
    blocks = []
    try:
        for z, b in zip(Z_BINS, BIAS):
            pk, f = growth_grid(cosmo, float(z), K_H)
            ps = pk[:, None] * (b + f[:, None] * MU[None, :] ** 2) ** 2
            blocks.append(np.log(ps).reshape(-1))
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()
    return np.concatenate(blocks)


def integration_weights(total_volume_h3_gpc3):
    dk = np.gradient(K_H)
    volume_per_bin = total_volume_h3_gpc3 * 1.0e9 / len(Z_BINS)
    one = (volume_per_bin / (8.0 * np.pi**2)) * (K_H[:, None] ** 2) * dk[:, None] * MU_W[None, :]
    return np.tile(one.reshape(-1), len(Z_BINS))


def bias_derivatives(reference_logp, psd_path, mass_eV):
    cosmo = Class()
    cosmo.set(class_parameters(psd_path, mass_eV))
    cosmo.compute()
    out = []
    block = len(K_H) * len(MU)
    try:
        for iz, (z, b) in enumerate(zip(Z_BINS, BIAS)):
            _, f = growth_grid(cosmo, float(z), K_H)
            d = 2.0 * b / (b + f[:, None] * MU[None, :] ** 2)
            full = np.zeros_like(reference_logp)
            full[iz * block:(iz + 1) * block] = np.broadcast_to(d, (len(K_H), len(MU))).reshape(-1)
            out.append(full)
    finally:
        cosmo.struct_cleanup()
        cosmo.empty()
    return out


def nuisance_matrix(psd_path, mass_eV, reference):
    specs = [
        ("H0", "H0", 0.25, cro.H0),
        ("omega_b", "omega_b", 1.0e-4, cro.OMEGA_B),
        ("omega_cdm", "omega_cdm", 4.0e-4, cro.OMEGA_CDM),
        ("lnA_s", "A_s", 0.01, cro.A_S),
        ("n_s", "n_s", 0.003, NS),
    ]
    cols = []
    names = []
    for label, key, step, center in specs:
        if label == "lnA_s":
            op = observable_from_psd(psd_path, mass_eV, {key: center * np.exp(step)})
            om = observable_from_psd(psd_path, mass_eV, {key: center * np.exp(-step)})
            deriv = (op - om) / (2.0 * step)
        else:
            op = observable_from_psd(psd_path, mass_eV, {key: center + step})
            om = observable_from_psd(psd_path, mass_eV, {key: center - step})
            deriv = (op - om) / (2.0 * step)
        cols.append(deriv)
        names.append(label)
    for iz, d in enumerate(bias_derivatives(reference, psd_path, mass_eV)):
        cols.append(d)
        names.append(f"lnb_z{Z_BINS[iz]:g}")
    return np.column_stack(cols), names


def weighted_projector_basis(nuisance, weights):
    nw = np.sqrt(weights)[:, None] * nuisance
    q, r = np.linalg.qr(nw, mode="reduced")
    diag = np.abs(np.diag(r))
    keep = diag > 1e-10 * max(1.0, float(np.max(diag)))
    return q[:, keep]


def sn_metrics(delta, weights, qn):
    dw = np.sqrt(weights) * delta
    fixed2 = float(dw @ dw)
    proj = dw - qn @ (qn.T @ dw)
    proj2 = float(proj @ proj)
    fixed = float(np.sqrt(max(fixed2, 0.0)))
    projected = float(np.sqrt(max(proj2, 0.0)))
    retained = float(proj2 / fixed2) if fixed2 > 0 else 0.0
    return fixed, projected, retained


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="rsd_hidden_state_output")
    ap.add_argument("--mass", type=float, default=0.60)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.01)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--volume", type=float, default=50.0, help="total effective volume in (h^-1 Gpc)^3")
    ap.add_argument("--samples", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    q = np.linspace(0.0, 20.0, 4000)
    f0, weights_mom, basis, N, shapes, y, M = cro.kinetic_objects(q, args.mass, args.z_match)
    f0_path = out / "fd_reference.dat"
    write_psd(f0_path, q, f0)

    ref = observable_from_psd(f0_path, args.mass)
    fisher_w = integration_weights(args.volume)
    nuisance, nuisance_names = nuisance_matrix(f0_path, args.mass, ref)
    qn = weighted_projector_basis(nuisance, fisher_w)

    response_cols = []
    for j, shape in enumerate(shapes):
        fp = f0 + args.probe_frac * shape
        fm = f0 - args.probe_frac * shape
        pp = out / f"probe_{j:02d}_plus.dat"
        pm = out / f"probe_{j:02d}_minus.dat"
        write_psd(pp, q, fp)
        write_psd(pm, q, fm)
        op = observable_from_psd(pp, args.mass)
        om = observable_from_psd(pm, args.mass)
        response_cols.append((op - om) / (2.0 * args.probe_frac))
    R = np.column_stack(response_cols)

    rw = np.sqrt(fisher_w)[:, None] * R
    rproj = rw - qn @ (qn.T @ rw)
    Bfixed = rw.T @ rw
    Bproj = rproj.T @ rproj

    rng = np.random.default_rng(args.seed)
    candidates = []
    evals, evecs = np.linalg.eigh(Bproj)
    candidates.append(evecs[:, np.argmax(evals)])
    for j in range(shapes.shape[0]):
        e = np.zeros(shapes.shape[0])
        e[j] = 1.0
        candidates.extend([e, -e])
    random = rng.normal(size=(args.samples, shapes.shape[0]))
    random /= np.linalg.norm(random, axis=1, keepdims=True)
    candidates.extend(random)

    best = None
    best_projected = -np.inf
    for c in candidates:
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            continue
        norm = 1.0 / maxrel
        amp = 2.0 * args.final_frac * norm
        proj2 = float(amp * amp * (c @ Bproj @ c))
        if proj2 > best_projected:
            fixed2 = float(amp * amp * (c @ Bfixed @ c))
            best_projected = proj2
            best = (c.copy(), norm, raw * norm, fixed2)

    if best is None:
        raise RuntimeError("No valid optimized direction found")

    coeff, norm, shape, pred_fixed2 = best
    fp = f0 + args.final_frac * shape
    fm = f0 - args.final_frac * shape
    pplus = out / "optimized_plus.dat"
    pminus = out / "optimized_minus.dat"
    write_psd(pplus, q, fp)
    write_psd(pminus, q, fm)

    op = observable_from_psd(pplus, args.mass)
    om = observable_from_psd(pminus, args.mass)
    delta = op - om
    fixed, projected, retained = sn_metrics(delta, fisher_w, qn)

    mp = cro.moments(fp, q, weights_mom)
    mm = cro.moments(fm, q, weights_mom)
    denom = np.maximum(0.5 * (np.abs(mp) + np.abs(mm)), 1e-300)
    mismatch = np.abs(mp - mm) / denom

    summary = {
        "class_commit": cro.CLASS_COMMIT,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "mass_over_Tnu_at_match": float(y),
        "redshift_bins": Z_BINS.tolist(),
        "k_h_Mpc_min": float(K_H.min()),
        "k_h_Mpc_max": float(K_H.max()),
        "n_k": int(len(K_H)),
        "n_mu": int(len(MU)),
        "total_effective_volume_hminus3_Gpc3": args.volume,
        "bias_model": "b(z)=1+0.84 z, with independent ln b nuisance in every redshift bin",
        "nuisance_parameters": nuisance_names,
        "probe_fractional_distortion": args.probe_frac,
        "final_fractional_distortion_cap": args.final_frac,
        "max_relative_moment_mismatch": float(mismatch.max()),
        "predicted_linear_fixed_SN": float(np.sqrt(max(pred_fixed2, 0.0))),
        "predicted_linear_projected_SN": float(np.sqrt(max(best_projected, 0.0))),
        "validated_fixed_SN": fixed,
        "validated_projected_SN": projected,
        "validated_retained_delta_chi2_fraction": retained,
        "optimized_coefficients": coeff.tolist(),
        "pointwise_normalization_factor": float(norm),
        "random_samples": args.samples,
        "random_seed": args.seed,
        "interpretation": "Ideal cosmic-variance-limited linear Kaiser RSD diagnostic, not a survey likelihood."
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    np.savetxt(out / "optimized_pair.csv", np.column_stack([q, f0, fp, fm, shape]), delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape", comments="")
    np.savetxt(out / "delta_log_Ps.csv", np.column_stack([np.arange(delta.size), delta]), delimiter=",",
               header="flattened_index,delta_log_Ps", comments="")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
