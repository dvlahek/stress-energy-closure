#!/usr/bin/env python3
"""
Noise sweep for the finite gravitational-response inversion.

This is a stability illustration, not part of the injectivity proof.
For each nonzero noise level, three deterministic noise realizations are used
for each of the two stress-energy-matched hidden states.

Outputs are written to response_identifiability_output/.
"""

from pathlib import Path
import json
import numpy as np
from scipy.optimize import nnls
import matplotlib.pyplot as plt

from injectivity_tomography import (
    PMAX, TAUMAX, N_P, N_TAU,
    construct_matched_pair, spline_forward_operator
)

NOISE_LEVELS = [0.0, 1e-5, 1e-4, 1e-3, 3e-3, 1e-2]
SEEDS = [9, 17, 31]
LAMBDAS = np.logspace(-9, -3, 13)


def unpack_pair(result):
    if isinstance(result, dict):
        return (
            result["F_plus"], result["F_minus"],
            result["K_plus"], result["K_minus"]
        )
    return result


def reconstruct_fast(M, B, L, K_true, F_true, noise_fraction, seed):
    scale = np.max(np.abs(K_true))
    rng = np.random.default_rng(seed)
    y = K_true + rng.normal(0.0, noise_fraction * scale, len(K_true))

    idx = np.arange(len(y))
    train = idx[idx % 4 != 0]
    valid = idx[idx % 4 == 0]

    best = None
    for lam in LAMBDAS:
        A = np.vstack([M[train] / scale, np.sqrt(lam) * L])
        b = np.r_[y[train] / scale, np.zeros(L.shape[0])]
        x, _ = nnls(A, b, maxiter=2000)
        val = np.linalg.norm(M[valid] @ x - y[valid]) / (
            np.sqrt(len(valid)) * scale
        )
        if best is None or val < best[0]:
            best = (val, lam)

    lam = best[1]
    A = np.vstack([M / scale, np.sqrt(lam) * L])
    b = np.r_[y / scale, np.zeros(L.shape[0])]
    x, _ = nnls(A, b, maxiter=2000)

    F_hat = x @ B
    K_hat = M @ x

    f_error = np.linalg.norm(F_hat - F_true) / np.linalg.norm(F_true)
    k_error = np.linalg.norm(K_hat - K_true) / np.linalg.norm(K_true)

    return float(f_error), float(k_error), float(lam)


def main():
    out = Path("response_identifiability_output")
    out.mkdir(exist_ok=True)

    p = np.linspace(0.0, PMAX, N_P)
    taus = np.linspace(0.0, TAUMAX, N_TAU)
    Fp, Fm, Kp, Km = unpack_pair(construct_matched_pair(p, taus))

    M, B = spline_forward_operator(p, taus, n_basis=24)
    n_coeff = M.shape[1]
    L = np.zeros((n_coeff - 2, n_coeff))
    for i in range(n_coeff - 2):
        L[i, i:i+3] = [1.0, -2.0, 1.0]

    rows = []
    for noise in NOISE_LEVELS:
        seeds = [SEEDS[0]] if noise == 0.0 else SEEDS
        for seed in seeds:
            for state, K, F in [("plus", Kp, Fp), ("minus", Km, Fm)]:
                f_error, k_error, lam = reconstruct_fast(
                    M, B, L, K, F, noise, seed
                )
                rows.append([
                    noise, seed, 1 if state == "plus" else -1,
                    f_error, k_error, lam
                ])

    arr = np.asarray(rows, dtype=float)
    np.savetxt(
        out / "noise_sweep.csv", arr, delimiter=",",
        header="noise_fraction,seed,state_sign,relative_F_error,relative_K_error,lambda",
        comments=""
    )

    summary = []
    for noise in NOISE_LEVELS:
        vals = arr[arr[:, 0] == noise, 3]
        summary.append({
            "noise_fraction": float(noise),
            "n_reconstructions": int(len(vals)),
            "median_relative_F_error": float(np.median(vals)),
            "min_relative_F_error": float(np.min(vals)),
            "max_relative_F_error": float(np.max(vals)),
            "q25_relative_F_error": float(np.quantile(vals, 0.25)),
            "q75_relative_F_error": float(np.quantile(vals, 0.75)),
        })

    with open(out / "noise_sweep_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    x = np.array([r["noise_fraction"] for r in summary])
    y = np.array([r["median_relative_F_error"] for r in summary])
    lo = y - np.array([r["min_relative_F_error"] for r in summary])
    hi = np.array([r["max_relative_F_error"] for r in summary]) - y

    xp = x.copy()
    xp[0] = 1e-6

    plt.figure(figsize=(7.0, 4.5))
    plt.errorbar(xp, y, yerr=np.vstack([lo, hi]), marker="o", capsize=3)
    plt.xscale("log")
    plt.xlabel("response-noise fraction")
    plt.ylabel("relative kinetic-profile reconstruction error")
    plt.title("Finite noisy inversion degrades despite exact identifiability")
    plt.tight_layout()
    plt.savefig(out / "noise_sweep.png", dpi=180)
    plt.close()

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
