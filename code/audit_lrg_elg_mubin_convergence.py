#!/usr/bin/env python3
"""Audit numerical convergence of the DESI z-resolved odd estimator in mu binning."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

FILENAME = "lrg_elg_exact_zresolved_odd_multipoles.csv"


def load(path):
    p = Path(path)
    if p.is_dir():
        p = p / FILENAME
    a = np.genfromtxt(p, delimiter=",", names=True)
    keys = np.column_stack([a["zlo"], a["zhi"], a["s_Mpc_over_h"]]).astype(float)
    return p, keys, np.asarray(a["xi1_LRG_to_ELG"], float), np.asarray(a["xi3_LRG_to_ELG"], float)


def metric(delta, C, n):
    p = len(delta)
    h = (n - p - 2) / (n - 1) if n > p + 2 else None
    if h is None:
        return None
    P = h * np.linalg.pinv(C, rcond=1e-12)
    return float(delta @ P @ delta)


def summarize(delta, sigma=None, C=None, n=None):
    out = {
        "rms_abs": float(np.sqrt(np.mean(delta * delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }
    if sigma is not None:
        z = np.abs(delta) / sigma
        out["median_abs_shift_in_mock_sigma"] = float(np.median(z))
        out["max_abs_shift_in_mock_sigma"] = float(np.max(z))
    if C is not None and n is not None:
        q = metric(delta, C, n)
        out["hartlap_metric_delta_chi2"] = q
        out["hartlap_metric_sqrt_delta_chi2"] = float(np.sqrt(max(q, 0.0))) if q is not None else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu120", required=True)
    ap.add_argument("--mu240", required=True)
    ap.add_argument("--mu480", required=True)
    ap.add_argument("--covariance")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    p120, k120, d120, o120 = load(args.mu120)
    p240, k240, d240, o240 = load(args.mu240)
    p480, k480, d480, o480 = load(args.mu480)
    if not (np.allclose(k120, k240, rtol=0, atol=1e-12) and np.allclose(k480, k240, rtol=0, atol=1e-12)):
        raise RuntimeError("120/240/480 measurement grids differ")

    cov = None
    if args.covariance:
        cov = np.load(args.covariance, allow_pickle=False)
        n = int(np.asarray(cov["n_mocks"]).item())
        C1 = np.asarray(cov["cov_xi1"], float)
        C3 = np.asarray(cov["cov_xi3"], float)
        s1 = np.sqrt(np.diag(C1))
        s3 = np.sqrt(np.diag(C3))
    else:
        n = None
        C1 = C3 = s1 = s3 = None

    result = {
        "scope": "DESI DR1 r4 full-sky mu-bin convergence; 240 is the frozen production baseline",
        "files": {"mu120": str(p120), "mu240": str(p240), "mu480": str(p480)},
        "dipole": {
            "mu120_minus_mu240": summarize(d120-d240, s1, C1, n),
            "mu480_minus_mu240": summarize(d480-d240, s1, C1, n),
            "mu480_minus_mu120": summarize(d480-d120, s1, C1, n),
        },
        "octupole": {
            "mu120_minus_mu240": summarize(o120-o240, s3, C3, n),
            "mu480_minus_mu240": summarize(o480-o240, s3, C3, n),
            "mu480_minus_mu120": summarize(o480-o120, s3, C3, n),
        },
        "analysis_scope": (
            "We compare 120, 240 and 480 angular bins at fixed catalogue selection and pair weighting. "
            "The numerical differences are evaluated against mock scatter."
        ),
    }
    Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
