#!/usr/bin/env python3
"""Audit finite-mock significance and data/mock processing consistency.

This is a diagnostic, not a final likelihood. It reports:
- centered vs zero-centered quadratic forms,
- raw vs Hartlap-corrected precision,
- the maximum possible fixed linear matched-filter |Z| within each vector space,
- r1 -> r4 random-density shift in mock-covariance units,
- leave-one-out empirical Mahalanobis ranks for the supplied mock ensemble.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm

FILENAME = "lrg_elg_exact_odd_multipoles.csv"


def load_vec(path):
    p = Path(path)
    if p.is_dir():
        p = p / FILENAME
    a = np.genfromtxt(p, delimiter=",", names=True)
    s = np.atleast_1d(np.asarray(a["s_Mpc_over_h"], dtype="f8"))
    x1 = np.atleast_1d(np.asarray(a["xi1_LRG_to_ELG"], dtype="f8"))
    x3 = np.atleast_1d(np.asarray(a["xi3_LRG_to_ELG"], dtype="f8"))
    return s, x1, x3


def alpha_hartlap(n, p):
    return float((n - p - 2) / (n - 1)) if n > p + 2 else None


def quad(y, C, nmock):
    inv = np.linalg.pinv(C, rcond=1e-12)
    qraw = float(y @ inv @ y)
    a = alpha_hartlap(nmock, len(y))
    qh = float(a * qraw) if a is not None else None
    def stats(q):
        if q is None:
            return None
        p = float(chi2.sf(q, len(y)))
        z = float(norm.isf(p)) if 0.0 < p < 1.0 else None
        return {
            "chi2": q,
            "dof": int(len(y)),
            "upper_tail_p": p,
            "gaussian_equivalent_signed_upper_tail_z": z,
            "max_abs_fixed_linear_matched_filter_z": float(np.sqrt(max(q, 0.0))),
        }
    return {
        "raw": stats(qraw),
        "hartlap": stats(qh),
        "hartlap_factor": a,
    }


def pack(x1, x3, which):
    if which == "dipole":
        return x1
    if which == "octupole":
        return x3
    return np.concatenate([x1, x3])


def cov_for(bundle, which):
    return np.asarray(bundle[
        {"dipole": "cov_xi1", "octupole": "cov_xi3", "joint": "cov_joint"}[which]
    ], dtype="f8")


def loo_rank(data_vec, mock_vectors, which):
    X = np.asarray(mock_vectors, dtype="f8")
    n = len(X)
    if n < 4:
        return None
    # Data score from all mocks.
    mean = X.mean(axis=0)
    C = np.cov(X, rowvar=False, ddof=1)
    yd = data_vec - mean
    ad = alpha_hartlap(n, len(yd))
    if ad is None:
        return None
    qd = float(ad * (yd @ np.linalg.pinv(C, rcond=1e-12) @ yd))

    qs = []
    for i in range(n):
        train = np.delete(X, i, axis=0)
        nt = len(train)
        a = alpha_hartlap(nt, X.shape[1])
        if a is None:
            return None
        m = train.mean(axis=0)
        c = np.cov(train, rowvar=False, ddof=1)
        y = X[i] - m
        q = float(a * (y @ np.linalg.pinv(c, rcond=1e-12) @ y))
        qs.append(q)
    qs = np.asarray(qs)
    nge = int(np.sum(qs >= qd))
    return {
        "data_q_hartlap": qd,
        "loo_mock_q_median": float(np.median(qs)),
        "loo_mock_q_min": float(np.min(qs)),
        "loo_mock_q_max": float(np.max(qs)),
        "n_mock_scores_ge_data": nge,
        "empirical_rank_p_plus_one": float((nge + 1) / (n + 1)),
        "note": "Finite-N diagnostic only; LOO mock scores use N-1 training mocks while the data score uses all N.",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-r4", required=True)
    ap.add_argument("--covariance", required=True)
    ap.add_argument("--mock-glob", required=True)
    ap.add_argument("--data-r1", help="Optional r1 data vector for random-density shift audit")
    ap.add_argument("--out", default="lrg_elg_significance_audit.json")
    args = ap.parse_args()

    s4, d41, d43 = load_vec(args.data_r4)
    b = np.load(args.covariance, allow_pickle=False)
    sc = np.asarray(b["separation"], dtype="f8")
    if not np.allclose(s4, sc, rtol=0.0, atol=1e-12):
        raise RuntimeError("Data and covariance separation grids differ")
    nmock = int(np.asarray(b["n_mocks"]).item())
    m1 = np.asarray(b["mean_xi1"], dtype="f8")
    m3 = np.asarray(b["mean_xi3"], dtype="f8")

    mock_paths = sorted(glob.glob(args.mock_glob))
    mock_vecs = {"dipole": [], "octupole": [], "joint": []}
    for p in mock_paths:
        s, x1, x3 = load_vec(p)
        if not np.allclose(s, s4, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"Mock grid differs: {p}")
        for w in mock_vecs:
            mock_vecs[w].append(pack(x1, x3, w))

    out = {
        "n_mocks_covariance": nmock,
        "n_mocks_glob": len(mock_paths),
        "tests": {},
    }
    for w in ("dipole", "octupole", "joint"):
        d = pack(d41, d43, w)
        mean = pack(m1, m3, w)
        C = cov_for(b, w)
        out["tests"][w] = {
            "data_minus_mock_mean": quad(d - mean, C, nmock),
            "data_vs_zero": quad(d, C, nmock),
            "leave_one_out_rank": loo_rank(d, mock_vecs[w], w) if len(mock_paths) >= 2 else None,
        }

    if args.data_r1:
        s1, d11, d13 = load_vec(args.data_r1)
        if not np.allclose(s1, s4, rtol=0.0, atol=1e-12):
            raise RuntimeError("r1 and r4 data grids differ")
        out["random_density_shift_r1_to_r4"] = {}
        for w in ("dipole", "octupole", "joint"):
            delta = pack(d41, d43, w) - pack(d11, d13, w)
            C = cov_for(b, w)
            out["random_density_shift_r1_to_r4"][w] = {
                "rms_shift": float(np.sqrt(np.mean(delta**2))),
                "survey_covariance_metric": quad(delta, C, nmock),
                "note": "This is a scale diagnostic, not an estimate of random-catalog covariance; paired r1/r4 mocks are needed for that.",
            }

    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
