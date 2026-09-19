#!/usr/bin/env python3
"""Build covariance matrices for the exact DESI DR1 LRG x ELG odd-sector estimator.

Each input run must contain lrg_elg_exact_odd_multipoles.csv produced by
code/desi_dr1_lrg_elg_exact.py. The vector ordering is
[xi1(s_1..s_N), xi3(s_1..s_N)].
"""
from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np


FILENAME = "lrg_elg_exact_odd_multipoles.csv"


def load_run(path):
    path = Path(path)
    fn = path / FILENAME if path.is_dir() else path
    if not fn.exists():
        raise FileNotFoundError(f"{fn} not found")
    arr = np.genfromtxt(fn, delimiter=",", names=True)
    s = np.atleast_1d(np.asarray(arr["s_Mpc_over_h"], dtype="f8"))
    x1 = np.atleast_1d(np.asarray(arr["xi1_LRG_to_ELG"], dtype="f8"))
    x3 = np.atleast_1d(np.asarray(arr["xi3_LRG_to_ELG"], dtype="f8"))
    if not (np.all(np.isfinite(s)) and np.all(np.isfinite(x1)) and np.all(np.isfinite(x3))):
        raise RuntimeError(f"Non-finite values in {fn}")
    return fn, s, x1, x3


def hartlap(nmock, ndim):
    if nmock <= ndim + 2:
        return None
    return float((nmock - ndim - 2) / (nmock - 1))


def diagnostics(cov):
    eig = np.linalg.eigvalsh(cov)
    return {
        "dimension": int(cov.shape[0]),
        "condition_number": float(np.linalg.cond(cov)),
        "min_eigenvalue": float(eig[0]),
        "max_eigenvalue": float(eig[-1]),
        "positive_definite": bool(np.all(eig > 0.0)),
    }


def collect_runs(args):
    runs = []
    if args.runs:
        runs.extend(args.runs)
    if args.runs_file:
        for line in Path(args.runs_file).read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                runs.append(line)
    if args.glob_pattern:
        runs.extend(sorted(glob.glob(args.glob_pattern)))
    out = []
    seen = set()
    for run in runs:
        key = str(Path(run))
        if key not in seen:
            out.append(run)
            seen.add(key)
    if len(out) < 2:
        raise RuntimeError("Need at least two mock runs")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=[])
    ap.add_argument("--runs-file")
    ap.add_argument("--glob", dest="glob_pattern")
    ap.add_argument("--out", required=True, help="Output .npz covariance bundle")
    ap.add_argument("--summary", help="Optional JSON diagnostics path")
    args = ap.parse_args()

    runs = collect_runs(args)
    sref = None
    vectors = []
    sources = []
    for run in runs:
        fn, s, x1, x3 = load_run(run)
        if sref is None:
            sref = s
        elif not np.allclose(s, sref, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"Separation grid differs in {fn}")
        vectors.append(np.concatenate([x1, x3]))
        sources.append(str(fn))

    X = np.asarray(vectors, dtype="f8")
    if not np.all(np.isfinite(X)):
        raise RuntimeError("Mock matrix contains non-finite values")
    nmock, ndim = X.shape
    nb = len(sref)
    if ndim != 2 * nb:
        raise RuntimeError("Unexpected mock-vector dimension")

    mean = X.mean(axis=0)
    cov_joint = np.cov(X, rowvar=False, ddof=1)
    cov_xi1 = cov_joint[:nb, :nb]
    cov_xi3 = cov_joint[nb:, nb:]
    cov_xi13 = cov_joint[:nb, nb:]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        separation=np.asarray(sref, dtype="f8"),
        mean_joint=mean,
        mean_xi1=mean[:nb],
        mean_xi3=mean[nb:],
        cov_joint=cov_joint,
        cov_xi1=cov_xi1,
        cov_xi3=cov_xi3,
        cov_xi13=cov_xi13,
        n_mocks=np.asarray(nmock, dtype="i8"),
        vector_order=np.asarray("xi1_then_xi3"),
    )

    summary = {
        "scope": "Mock covariance for exact DESI DR1 LRGxELG odd-sector estimator",
        "n_mocks": int(nmock),
        "n_bins_per_multipole": int(nb),
        "vector_order": "xi1(s bins), then xi3(s bins)",
        "separation_Mpc_over_h": sref.tolist(),
        "hartlap": {
            "xi1": hartlap(nmock, nb),
            "xi3": hartlap(nmock, nb),
            "joint": hartlap(nmock, 2 * nb),
        },
        "diagnostics": {
            "xi1": diagnostics(cov_xi1),
            "xi3": diagnostics(cov_xi3),
            "joint": diagnostics(cov_joint),
        },
        "mean_mock_vector": {
            "xi1": mean[:nb].tolist(),
            "xi3": mean[nb:].tolist(),
        },
        "guardrail": (
            "Covariance is valid only for data vectors measured with the identical estimator, "
            "selection, weights, separation bins, mu bins, theta cut, and survey window."
        ),
        "sources_count": len(sources),
    }
    summary_path = Path(args.summary) if args.summary else out.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
