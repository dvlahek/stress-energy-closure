#!/usr/bin/env python3
"""Numerical support fix for the frozen LSS observability campaign.

The original campaign correctly avoids extrapolating P(k): interp_log returns NaN
outside each CLASS table's support.  The Fisher integrator then passed those NaNs
into the nuisance matrix when tiny endpoint differences occurred between cases.
This wrapper changes no physics, parameters, grids, caps, derivative steps, volume
convention, or decision gates.  It only restricts each redshift integral to the
common finite support of the already-computed fiducial, probe, and nuisance tables.
"""
from __future__ import annotations

import math
import numpy as np

import class_lss_observability_campaign as base


def response_arrays_finite_support(results, manifest, kmax):
    nd = int(manifest["null_dimension"])
    pars = list(base.DERIVATIVE_STEPS)
    nz = len(base.REDSHIFTS)
    Vbin = 1.0e9 / nz

    G = np.zeros((nd, nd))
    C = np.zeros((nd, len(pars)))
    F = np.zeros((len(pars), len(pars)))
    frozen_sn2 = 0.0
    per_z = {}

    for z in base.REDSHIFTS:
        k0, p0 = base.load_pk(results, "fiducial", z)
        sel0 = (k0 > 0) & (k0 <= kmax)
        k = k0[sel0]
        p0 = p0[sel0]
        if len(k) < 8:
            raise RuntimeError(f"Too few k points at z={z}, kmax={kmax}")

        B = np.zeros((len(k), nd))
        for j in range(nd):
            kp, pp = base.load_pk(results, f"probe_{j:02d}_plus", z)
            km, pm = base.load_pk(results, f"probe_{j:02d}_minus", z)
            B[:, j] = base.interp_log(k, kp, pp) - base.interp_log(k, km, pm)

        D = np.zeros((len(k), len(pars)))
        for a, par in enumerate(pars):
            kp, pp = base.load_pk(results, f"d_{par}_plus", z)
            km, pm = base.load_pk(results, f"d_{par}_minus", z)
            step = float(base.DERIVATIVE_STEPS[par])
            D[:, a] = (base.interp_log(k, kp, pp) - base.interp_log(k, km, pm)) / (2.0 * step)

        # No extrapolation: integrate only where every response entering the
        # same Fisher block is finite. This removes endpoint NaNs caused by
        # slightly different CLASS k supports across otherwise valid runs.
        finite = np.isfinite(k) & np.isfinite(p0) & (p0 > 0)
        finite &= np.all(np.isfinite(B), axis=1)
        finite &= np.all(np.isfinite(D), axis=1)
        k = k[finite]
        p0 = p0[finite]
        B = B[finite]
        D = D[finite]
        if len(k) < 8:
            raise RuntimeError(f"Too few common finite k points at z={z}, kmax={kmax}")

        pref = math.sqrt(Vbin / (4.0 * math.pi**2))
        W = pref * k[:, None] / p0[:, None]
        G += base.trap_matrix(k, W * B)
        F += base.trap_matrix(k, W * D)
        for j in range(nd):
            for a in range(len(pars)):
                C[j, a] += float(np.trapezoid(
                    (pref * k * B[:, j] / p0) *
                    (pref * k * D[:, a] / p0), k
                ))

        kf, pfp, pfm = base.common_pair(results, "frozen_plus", "frozen_minus", z)
        sf = (kf > 0) & (kf <= kmax)
        kf, pfp, pfm = kf[sf], pfp[sf], pfm[sf]
        finite_f = np.isfinite(kf) & np.isfinite(pfp) & np.isfinite(pfm)
        kf, pfp, pfm = kf[finite_f], pfp[finite_f], pfm[finite_f]
        pbar = 0.5 * (pfp + pfm)
        good = pbar > 0
        kf, pfp, pfm, pbar = kf[good], pfp[good], pfm[good], pbar[good]
        sn2z = Vbin / (4.0 * math.pi**2) * float(np.trapezoid(
            kf**2 * ((pfp - pfm) / pbar)**2, kf
        ))
        frozen_sn2 += sn2z
        per_z[str(z)] = {
            "frozen_SN_per_sqrt_total_1_Gpch3_equal_bin_weight": float(math.sqrt(max(sn2z, 0.0))),
            "max_abs_relative_difference_percent": float(100.0 * np.max(np.abs(pfp - pfm) / pbar)),
            "common_finite_k_points": int(len(k)),
            "common_finite_k_min_h_Mpc": float(k.min()),
            "common_finite_k_max_h_Mpc": float(k.max()),
        }

    if not (np.all(np.isfinite(G)) and np.all(np.isfinite(C)) and np.all(np.isfinite(F))):
        raise RuntimeError("Non-finite Fisher block remains after common-support restriction")

    Fs = 0.5 * (F + F.T)
    evals = np.linalg.eigvalsh(Fs)
    rank = int(np.linalg.matrix_rank(F, tol=max(np.max(np.abs(evals)), 1.0) * 1e-12))
    Fpinv = np.linalg.pinv(F, rcond=1e-10)
    raw_proj = G - C @ Fpinv @ C.T
    Gproj = 0.5 * (raw_proj + raw_proj.T)

    return {
        "G_fixed": G,
        "G_projected": Gproj,
        "cross": C,
        "F_nuisance": F,
        "fisher_eigenvalues": evals,
        "fisher_rank": rank,
        "parameters": pars,
        "frozen_sn2_per_total_1_Gpch3": frozen_sn2,
        "per_z": per_z,
        "numerical_support_policy": "intersection of finite CLASS P(k) support across fiducial, probes, and nuisance derivatives; no extrapolation",
    }


base.response_arrays = response_arrays_finite_support

if __name__ == "__main__":
    base.main()
