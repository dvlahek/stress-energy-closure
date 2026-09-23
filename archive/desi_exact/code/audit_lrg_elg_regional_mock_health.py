#!/usr/bin/env python3
"""Technical health audit for one regional z-resolved EZmock realization.

This is a pipeline diagnostic only. It does not decide whether a statistically
extreme mock should be removed. A mock should remain in the ensemble unless an
objective processing/catalog failure is identified.
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

CSV = "lrg_elg_exact_zresolved_odd_multipoles.csv"
SUMMARY = "summary_exact_zresolved.json"


def load_vec(path: Path):
    a = np.genfromtxt(path, delimiter=",", names=True)
    return np.asarray(a["xi1_LRG_to_ELG"], float)


def safe_z(x, arr):
    arr = np.asarray(arr, float)
    s = arr.std(ddof=1)
    return float((x - arr.mean()) / s) if s > 0 else None


def percentile_rank(x, arr):
    arr = np.asarray(arr, float)
    return float(np.mean(arr <= x))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-root", required=True)
    ap.add_argument("--mock-id", type=int, required=True)
    ap.add_argument("--rerun-dir", help="Optional independent rerun directory with reverse closure enabled")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    root = Path(args.mock_root)
    dirs = [Path(p) for p in sorted(glob.glob(str(root / "exact_z_*")))]
    dirs = [d for d in dirs if (d / CSV).is_file() and (d / SUMMARY).is_file()]
    if len(dirs) < 4:
        raise RuntimeError("Need at least four complete mock directories")

    target_name = f"exact_z_{args.mock_id:04d}"
    target = root / target_name
    if not (target / CSV).is_file() or not (target / SUMMARY).is_file():
        raise RuntimeError(f"Target mock output missing: {target}")

    vectors = np.asarray([load_vec(d / CSV) for d in dirs])
    names = [d.name for d in dirs]
    i = names.index(target_name)
    x = vectors[i]
    train = np.delete(vectors, i, axis=0)
    mean = train.mean(axis=0)
    C = np.cov(train, rowvar=False, ddof=1)
    inv = np.linalg.pinv(C, rcond=1e-12)
    p = len(x)
    n = len(train)
    h = float((n - p - 2) / (n - 1)) if n > p + 2 else None
    dx = x - mean
    q_raw = float(dx @ inv @ dx)
    q_h = float(h * q_raw) if h is not None else None

    summaries = [json.loads((d / SUMMARY).read_text()) for d in dirs]
    ts = summaries[i]
    count_audit = []
    for iz in range(len(ts["bins"])):
        for key in ("LRG", "ELG", "LRG_random", "ELG_random"):
            vals = [s["bins"][iz]["counts"][key] for s in summaries]
            val = ts["bins"][iz]["counts"][key]
            count_audit.append({
                "zbin": iz,
                "catalog": key,
                "target": int(val),
                "ensemble_mean": float(np.mean(vals)),
                "ensemble_std": float(np.std(vals, ddof=1)),
                "zscore_including_target": safe_z(val, vals),
                "percentile_rank_including_target": percentile_rank(val, vals),
                "min": int(np.min(vals)),
                "max": int(np.max(vals)),
            })

    result = {
        "scope": "Technical health audit for one regional z-resolved EZmock",
        "mock_root": str(root),
        "mock_id": args.mock_id,
        "n_complete_mocks": len(dirs),
        "dipole_dimension": p,
        "target_dipole": {
            "rms": float(np.sqrt(np.mean(x*x))),
            "max_abs_component": float(np.max(np.abs(x))),
            "loo_zero_mean_mahalanobis_q_raw": q_raw,
            "loo_hartlap_factor": h,
            "loo_zero_mean_mahalanobis_q_hartlap": q_h,
        },
        "count_audit": count_audit,
        "stored_summary_sign_reversal": ts.get("bins", []),
        "rerun_comparison": None,
        "guardrail": (
            "Statistical extremeness alone is not a removal criterion. Exclude a mock only for an "
            "objective catalog, random, estimator, or closure failure defined independently of its matched-filter score."
        ),
    }

    if args.rerun_dir:
        rr = Path(args.rerun_dir)
        y = load_vec(rr / CSV)
        rs = json.loads((rr / SUMMARY).read_text())
        dif = y - x
        result["rerun_comparison"] = {
            "rerun_dir": str(rr),
            "max_abs_dipole_difference": float(np.max(np.abs(dif))),
            "rms_dipole_difference": float(np.sqrt(np.mean(dif*dif))),
            "rerun_bins": [
                {
                    "index": b["index"],
                    "counts": b["counts"],
                    "sign_reversal": b["sign_reversal"],
                } for b in rs["bins"]
            ],
        }

    Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
