#!/usr/bin/env python3
"""Build a compact shared angular-random cache for local Phase-7 EZmock runs.

The public EZmock clustering random FITS files are very large (~GB per
realization), while the placebo estimator only needs a density-controlled
subset.  This utility samples only angular positions and random weights from
one validated NGC/SGC random pair.  Redshifts are deliberately not retained:
the realization script assigns redshifts afresh from each mock data
realization within the same narrow-z and Galactic-cap cells.

This cache is a transport/performance device for the EZmock placebo
covariance/systematics control.  It does not turn the placebo split into a
luminosity- or halo-mass-matched mock analysis.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from astropy.io import fits


def _col(d, *names, default=None):
    avail = {n.upper(): n for n in d.names}
    for n in names:
        if n.upper() in avail:
            return np.asarray(d[avail[n.upper()]])
    if default is not None:
        return default
    raise KeyError(f"missing {names}; available={list(d.names)}")


def sample_one(path: str, n_keep: int, cap_code: int, rng: np.random.Generator):
    with fits.open(path, memmap=True) as h:
        h.verify("exception")
        d = h[1].data
        n = len(d)
        if n_keep > n:
            raise ValueError(f"requested {n_keep} rows from {path}, only {n} available")
        idx = rng.choice(n, size=n_keep, replace=False)
        idx.sort()
        ra = _col(d, "RA")[idx].astype(np.float64)
        dec = _col(d, "DEC")[idx].astype(np.float64)
        w = _col(d, "WEIGHT", default=np.ones(n))[idx].astype(np.float64)
        wf = _col(d, "WEIGHT_FKP", default=np.ones(n))[idx].astype(np.float64)
        weight = w * wf
    good = np.isfinite(ra) & np.isfinite(dec) & np.isfinite(weight) & (weight > 0)
    return (
        ra[good],
        dec[good],
        weight[good],
        np.full(int(good.sum()), cap_code, dtype=np.int8),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ngc", required=True)
    ap.add_argument("--sgc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ngc-count", type=int, default=1_000_000)
    ap.add_argument("--sgc-count", type=int, default=500_000)
    ap.add_argument("--seed", type=int, default=20260914)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    parts = [
        sample_one(args.ngc, args.ngc_count, 0, rng),
        sample_one(args.sgc, args.sgc_count, 1, rng),
    ]
    ra = np.concatenate([x[0] for x in parts])
    dec = np.concatenate([x[1] for x in parts])
    w = np.concatenate([x[2] for x in parts])
    cap = np.concatenate([x[3] for x in parts])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        ra=ra,
        dec=dec,
        w=w,
        cap=cap,
        seed=np.asarray([args.seed], dtype=np.int64),
        ngc_source=np.asarray([str(Path(args.ngc))]),
        sgc_source=np.asarray([str(Path(args.sgc))]),
    )
    print(
        "PHASE7_EZMOCK_SHARED_RANDOM_CACHE",
        f"out={out}",
        f"rows={len(ra)}",
        f"NGC={(cap == 0).sum()}",
        f"SGC={(cap == 1).sum()}",
        f"bytes={out.stat().st_size}",
    )


if __name__ == "__main__":
    main()
