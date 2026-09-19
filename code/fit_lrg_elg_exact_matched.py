#!/usr/bin/env python3
"""Production matched-filter tests for the exact DESI DR1 LRG x ELG odd sector.

Consumes the exact data vector and a covariance bundle produced by
build_lrg_elg_mock_covariance.py. The wake fit is a fixed-shape linear
matched filter. Optional nuisance templates are profiled simultaneously.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_measurement(path):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    s = np.atleast_1d(np.asarray(arr["s_Mpc_over_h"], dtype="f8"))
    x1 = np.atleast_1d(np.asarray(arr["xi1_LRG_to_ELG"], dtype="f8"))
    x3 = np.atleast_1d(np.asarray(arr["xi3_LRG_to_ELG"], dtype="f8"))
    return s, x1, x3


def hartlap(nmock, ndim):
    if nmock is None or nmock <= ndim + 2:
        return 1.0, False
    return float((nmock - ndim - 2) / (nmock - 1)), True


def precision(cov, nmock=None):
    cov = np.asarray(cov, dtype="f8")
    alpha, applied = hartlap(nmock, cov.shape[0])
    return alpha * np.linalg.pinv(cov, rcond=1e-12), alpha, applied


def null_test(y, cov, nmock=None):
    psi, alpha, applied = precision(cov, nmock=nmock)
    q = float(y @ psi @ y)
    dof = int(len(y))
    p = float(chi2.sf(q, dof))
    z_raw = float(norm.isf(p)) if 0.0 < p < 1.0 else (float("inf") if p == 0.0 else float("-inf"))
    z_excess = max(0.0, z_raw)
    return {
        "chi2": q,
        "dof": dof,
        "pvalue": p,
        "gaussian_equivalent_excess_sigma_one_sided": z_excess,
        "hartlap_factor": alpha,
        "hartlap_applied": applied,
    }


def interp_columns(path, s, columns):
    tab = np.genfromtxt(path, delimiter=",", names=True)
    names = list(tab.dtype.names or [])
    if "s_Mpc_over_h" not in names:
        raise KeyError(f"{path} lacks s_Mpc_over_h")
    st = np.asarray(tab["s_Mpc_over_h"], dtype="f8")
    out = {}
    for name in columns:
        if name not in names:
            raise KeyError(f"{path} lacks requested template column {name}; available={names}")
        out[name] = np.interp(s, st, np.asarray(tab[name], dtype="f8"))
    return out


def gls_profile(y, cov, signal, nuisance, signal_name, nuisance_names, nmock=None):
    psi, alpha, applied = precision(cov, nmock=nmock)
    signal = np.asarray(signal, dtype="f8")
    nuisance = np.asarray(nuisance, dtype="f8")
    if nuisance.ndim == 1:
        nuisance = nuisance[:, None]
    if nuisance.size == 0:
        nuisance = np.empty((len(y), 0), dtype="f8")

    X0 = nuisance
    X1 = np.column_stack([nuisance, signal])

    def fit(X):
        if X.shape[1] == 0:
            resid = y.copy()
            return np.empty(0), np.empty((0, 0)), float(resid @ psi @ resid)
        fisher = X.T @ psi @ X
        fisher_inv = np.linalg.pinv(fisher, rcond=1e-12)
        beta = fisher_inv @ (X.T @ psi @ y)
        resid = y - X @ beta
        return beta, fisher_inv, float(resid @ psi @ resid)

    _, _, chi0 = fit(X0)
    b1, c1, chi1 = fit(X1)

    amp = float(b1[-1])
    sig = float(np.sqrt(max(c1[-1, -1], 0.0)))
    z_signed = float(amp / sig) if sig > 0.0 else None
    delta = max(0.0, float(chi0 - chi1))
    p_two_sided = float(chi2.sf(delta, 1))
    z_abs = float(norm.isf(p_two_sided / 2.0)) if p_two_sided > 0.0 else float("inf")

    if nuisance.shape[1]:
        fn = nuisance.T @ psi @ nuisance
        proj_coeff = np.linalg.pinv(fn, rcond=1e-12) @ (nuisance.T @ psi @ signal)
        signal_perp = signal - nuisance @ proj_coeff
    else:
        signal_perp = signal.copy()
    denom = float(signal_perp @ psi @ signal_perp)
    numerator = float(signal_perp @ psi @ y)
    matched_z = float(numerator / np.sqrt(denom)) if denom > 0.0 else None

    params = {}
    for i, name in enumerate(nuisance_names):
        err = float(np.sqrt(max(c1[i, i], 0.0)))
        params[name] = {
            "amplitude": float(b1[i]),
            "sigma": err,
            "z_signed": float(b1[i] / err) if err > 0.0 else None,
        }
    params[signal_name] = {
        "amplitude": amp,
        "sigma": sig,
        "z_signed": z_signed,
    }

    return {
        "parameters": params,
        "chi2_nuisance_only": chi0,
        "chi2_nuisance_plus_signal": chi1,
        "delta_chi2_signal": delta,
        "pvalue_signal_two_sided": p_two_sided,
        "gaussian_equivalent_sigma_two_sided": z_abs,
        "matched_filter_z_signed": matched_z,
        "hartlap_factor": alpha,
        "hartlap_applied": applied,
        "nuisance_rank": int(np.linalg.matrix_rank(nuisance)) if nuisance.shape[1] else 0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measurement", required=True,
                    help="Exact lrg_elg_exact_odd_multipoles.csv")
    ap.add_argument("--covariance", required=True,
                    help=".npz from build_lrg_elg_mock_covariance.py")
    ap.add_argument("--wake-template", required=True,
                    help="CSV containing s_Mpc_over_h and wake_shape")
    ap.add_argument("--nuisance-template",
                    help="Optional CSV containing nuisance template columns; defaults to wake-template")
    ap.add_argument("--nuisance-cols", default="doppler_shape",
                    help="Comma-separated dipole nuisance columns; use empty string for none")
    ap.add_argument("--wake-col", default="wake_shape")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--subtract-mock-mean", action="store_true",
                    help="Diagnostic only: subtract measured mock mean before tests")
    args = ap.parse_args()

    measurement = Path(args.measurement)
    covariance = Path(args.covariance)
    wake_template = Path(args.wake_template)
    nuisance_template = Path(args.nuisance_template) if args.nuisance_template else wake_template
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    s, x1, x3 = load_measurement(measurement)
    bundle = np.load(covariance, allow_pickle=False)
    scov = np.asarray(bundle["separation"], dtype="f8")
    if not np.allclose(s, scov, rtol=0.0, atol=1e-12):
        raise RuntimeError("Measurement and covariance separation grids differ")

    C1 = np.asarray(bundle["cov_xi1"], dtype="f8")
    C3 = np.asarray(bundle["cov_xi3"], dtype="f8")
    Cj = np.asarray(bundle["cov_joint"], dtype="f8")
    nmock = int(np.asarray(bundle["n_mocks"]).item()) if "n_mocks" in bundle.files else None
    nb = len(s)
    if Cj.shape != (2 * nb, 2 * nb):
        raise RuntimeError(f"Unexpected joint covariance shape {Cj.shape}")

    mock_mean1 = np.asarray(bundle["mean_xi1"], dtype="f8") if "mean_xi1" in bundle.files else np.zeros(nb)
    mock_mean3 = np.asarray(bundle["mean_xi3"], dtype="f8") if "mean_xi3" in bundle.files else np.zeros(nb)
    if args.subtract_mock_mean:
        x1_test = x1 - mock_mean1
        x3_test = x3 - mock_mean3
    else:
        x1_test = x1.copy()
        x3_test = x3.copy()
    joint = np.concatenate([x1_test, x3_test])

    wake = interp_columns(wake_template, s, [args.wake_col])[args.wake_col]
    nuisance_names = [x.strip() for x in args.nuisance_cols.split(",") if x.strip()]
    nuisance_dict = interp_columns(nuisance_template, s, nuisance_names) if nuisance_names else {}
    nuisance = (
        np.column_stack([nuisance_dict[name] for name in nuisance_names])
        if nuisance_names else np.empty((nb, 0), dtype="f8")
    )

    result = {
        "scope": "Exact DESI DR1 LRGxELG odd-sector null and fixed-shape matched-filter tests",
        "measurement": str(measurement),
        "covariance": str(covariance),
        "wake_template": str(wake_template),
        "wake_template_sha256": sha256_file(wake_template),
        "nuisance_template": str(nuisance_template),
        "nuisance_template_sha256": sha256_file(nuisance_template),
        "nuisance_columns": nuisance_names,
        "wake_column": args.wake_col,
        "n_mocks": nmock,
        "separation_Mpc_over_h": s.tolist(),
        "mock_mean_subtracted": bool(args.subtract_mock_mean),
        "null_tests": {
            "dipole": null_test(x1_test, C1, nmock=nmock),
            "octupole": null_test(x3_test, C3, nmock=nmock),
            "joint_dipole_octupole": null_test(joint, Cj, nmock=nmock),
        },
        "dipole_matched_filter": gls_profile(
            x1_test, C1, wake, nuisance, args.wake_col, nuisance_names, nmock=nmock
        ),
        "guardrails": [
            "The scale range is inherited from the exact measurement and is not optimized on the observed data.",
            "A higher matched-filter significance is interpretable only for templates fixed independently of this data vector.",
            "The final physical fit must replace provisional nuisance shapes with the frozen relativistic/Doppler, wide-angle, evolution and magnification model.",
            "Do not choose nuisance columns, redshift cuts or scale cuts after inspecting which choice maximizes significance.",
        ],
    }
    out = outdir / "matched_filter_summary.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
