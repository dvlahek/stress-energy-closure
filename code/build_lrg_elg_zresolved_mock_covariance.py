#!/usr/bin/env python3
"""Build covariance for exact z-resolved DESI DR1 LRG x ELG odd multipoles."""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np

FILENAME = "lrg_elg_exact_zresolved_odd_multipoles.csv"


def hartlap(n, p):
    return float((n-p-2)/(n-1)) if n > p + 2 else None


def diag(c):
    ev = np.linalg.eigvalsh(c)
    return {
        "dimension": int(c.shape[0]),
        "condition_number": float(np.linalg.cond(c)),
        "min_eigenvalue": float(ev[0]),
        "max_eigenvalue": float(ev[-1]),
        "positive_definite": bool(np.all(ev > 0)),
    }


def load(path):
    p = Path(path)
    if p.is_dir():
        p = p / FILENAME
    a = np.genfromtxt(p, delimiter=",", names=True)
    keys = np.column_stack([
        np.asarray(a["zlo"], float),
        np.asarray(a["zhi"], float),
        np.asarray(a["s_Mpc_over_h"], float),
    ])
    zef = np.asarray(a["z_effective"], float)
    x1 = np.asarray(a["xi1_LRG_to_ELG"], float)
    x3 = np.asarray(a["xi3_LRG_to_ELG"], float)
    if not np.all(np.isfinite(np.c_[keys, zef, x1, x3])):
        raise RuntimeError(f"Non-finite values in {p}")
    return p, keys, zef, x1, x3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", dest="pattern", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    runs = sorted(glob.glob(args.pattern))
    if len(runs) < 2:
        raise RuntimeError("Need at least two mock runs")

    keys0 = None
    X1, X3, zefs, sources = [], [], [], []
    for r in runs:
        p, keys, ze, x1, x3 = load(r)
        if keys0 is None:
            keys0 = keys
        elif not np.allclose(keys, keys0, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"z/s grid differs in {p}")
        X1.append(x1); X3.append(x3); zefs.append(ze); sources.append(str(p))

    X1 = np.asarray(X1, float)
    X3 = np.asarray(X3, float)
    Xj = np.concatenate([X1, X3], axis=1)
    n = X1.shape[0]
    p1 = X1.shape[1]
    pj = Xj.shape[1]

    C1 = np.cov(X1, rowvar=False, ddof=1)
    C3 = np.cov(X3, rowvar=False, ddof=1)
    Cj = np.cov(Xj, rowvar=False, ddof=1)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        zlo=keys0[:,0],
        zhi=keys0[:,1],
        separation=keys0[:,2],
        mean_z_effective=np.mean(np.asarray(zefs,float), axis=0),
        mean_xi1=X1.mean(axis=0),
        mean_xi3=X3.mean(axis=0),
        cov_xi1=C1,
        cov_xi3=C3,
        cov_joint=Cj,
        n_mocks=np.asarray(n, dtype="i8"),
        vector_order=np.asarray("xi1_zmajor_then_xi3_zmajor"),
    )

    summary = {
        "scope": "Mock covariance for exact z-resolved DESI DR1 LRGxELG odd sector",
        "n_mocks": int(n),
        "dipole_dimension": int(p1),
        "joint_dimension": int(pj),
        "frozen_grid": [
            {"zlo":float(a),"zhi":float(b),"s_Mpc_over_h":float(s)}
            for a,b,s in keys0
        ],
        "hartlap": {
            "xi1": hartlap(n,p1),
            "xi3": hartlap(n,p1),
            "joint": hartlap(n,pj),
        },
        "diagnostics": {
            "xi1": diag(C1),
            "xi3": diag(C3),
            "joint": diag(Cj),
        },
        "guardrail": (
            "With 40 mocks the 18D dipole covariance is preliminary; the 36D joint precision "
            "is expected to be extremely noisy and should not be used for final inference. "
            "Increase the mock ensemble before interpreting high-dimensional tails."
        ),
        "sources_count": len(sources),
    }
    out.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
