#!/usr/bin/env python3
"""Adversarial post-processing audit of the final LSS Fisher projection.

No CLASS calculations are performed here. The script reuses the already validated
screen and optimized spectra and asks if the strong LambdaCDM absorption reported
by the final LSS campaign is sensitive to

* inclusion of the nearly-null tau_reio matter-power derivative;
* which nuisance directions are projected out; or
* the Moore-Penrose pseudoinverse cutoff.

All nuisance-subset comparisons use the same common finite k support defined by
ALL six derivative spectra. Thus any change is caused by the projection itself,
not by a changing integration domain.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import class_lss_observability_campaign as base


PARAMS = list(base.DERIVATIVE_STEPS)
VOLUME_GPCH3 = 50.0
KMAX_VALUES = (0.10, 0.20, 0.30)
RCONDS = (1.0e-8, 1.0e-10, 1.0e-12, 1.0e-14)


def signal_projection_blocks(screen: Path, signal_root: Path,
                             case_plus: str, case_minus: str,
                             kmax: float):
    nz = len(base.REDSHIFTS)
    vbin = 1.0e9 / nz  # (Mpc/h)^3 for total 1 (Gpc/h)^3
    npar = len(PARAMS)
    F = np.zeros((npar, npar))
    b = np.zeros(npar)
    sn2 = 0.0
    points = {}

    for z in base.REDSHIFTS:
        k0, p0 = base.load_pk(screen, "fiducial", z)
        sel = (k0 > 0) & (k0 <= kmax)
        k = k0[sel]
        p0 = p0[sel]

        kp, pp = base.load_pk(signal_root, case_plus, z)
        km, pm = base.load_pk(signal_root, case_minus, z)
        pp_i = base.interp_log(k, kp, pp)
        pm_i = base.interp_log(k, km, pm)
        h = pp_i - pm_i

        D = np.zeros((len(k), npar))
        for a, par in enumerate(PARAMS):
            k1, p1 = base.load_pk(screen, f"d_{par}_plus", z)
            k2, p2 = base.load_pk(screen, f"d_{par}_minus", z)
            step = float(base.DERIVATIVE_STEPS[par])
            D[:, a] = (
                base.interp_log(k, k1, p1) - base.interp_log(k, k2, p2)
            ) / (2.0 * step)

        # Freeze the support to the intersection required by ALL six derivatives.
        finite = np.isfinite(k) & np.isfinite(p0) & (p0 > 0)
        finite &= np.isfinite(pp_i) & np.isfinite(pm_i) & np.isfinite(h)
        finite &= np.all(np.isfinite(D), axis=1)
        k, p0, h, D = k[finite], p0[finite], h[finite], D[finite]
        if len(k) < 8:
            raise RuntimeError(f"Too few common finite points at z={z}, kmax={kmax}")

        pref2 = vbin / (4.0 * math.pi**2)
        sn2 += pref2 * float(np.trapezoid(k**2 * (h / p0)**2, k))
        for a in range(npar):
            da = D[:, a]
            b[a] += pref2 * float(np.trapezoid(k**2 * h * da / p0**2, k))
            for c in range(a, npar):
                dc = D[:, c]
                val = pref2 * float(np.trapezoid(k**2 * da * dc / p0**2, k))
                F[a, c] += val
                if c != a:
                    F[c, a] += val

        points[str(z)] = {
            "n": int(len(k)),
            "kmin_h_Mpc": float(k.min()),
            "kmax_h_Mpc": float(k.max()),
        }

    if not (np.isfinite(sn2) and np.all(np.isfinite(F)) and np.all(np.isfinite(b))):
        raise RuntimeError("Non-finite Fisher block")
    return sn2, F, b, points


def project(sn2: float, F: np.ndarray, b: np.ndarray,
            subset: list[str], rcond: float):
    if not subset:
        marg2 = sn2
        absorbed = 0.0
        eig = []
        rank = 0
    else:
        idx = [PARAMS.index(p) for p in subset]
        Fs = F[np.ix_(idx, idx)]
        bs = b[idx]
        Fpinv = np.linalg.pinv(Fs, rcond=rcond)
        absorbed = float(bs @ Fpinv @ bs)
        # Round-off can only make this tiny negative near exact projection.
        marg2 = max(float(sn2 - absorbed), 0.0)
        eig = np.linalg.eigvalsh(0.5 * (Fs + Fs.T)).tolist()
        if eig:
            tol = max(max(abs(x) for x in eig), 1.0) * rcond
            rank = int(sum(abs(x) > tol for x in eig))
        else:
            rank = 0
    return {
        "parameters": subset,
        "rcond": float(rcond),
        "fixed_SN_V50": float(math.sqrt(max(sn2 * VOLUME_GPCH3, 0.0))),
        "marginalized_SN_V50": float(math.sqrt(max(marg2 * VOLUME_GPCH3, 0.0))),
        "retained_delta_chi2_fraction": float(marg2 / sn2 if sn2 > 0 else 0.0),
        "absorbed_delta_chi2_fraction": float(absorbed / sn2 if sn2 > 0 else 0.0),
        "fisher_eigenvalues": eig,
        "effective_rank_at_rcond": rank,
    }


def audit_signal(screen: Path, signal_root: Path,
                 case_plus: str, case_minus: str):
    result = {}
    for kmax in KMAX_VALUES:
        sn2, F, b, points = signal_projection_blocks(
            screen, signal_root, case_plus, case_minus, kmax
        )
        subsets = {
            "fixed": [],
            "lnAs_only": ["lnAs"],
            "ns_only": ["n_s"],
            "H0_only": ["H0"],
            "omega_b_only": ["omega_b"],
            "omega_cdm_only": ["omega_cdm"],
            "tau_only": ["tau_reio"],
            "amplitude_tilt": ["lnAs", "n_s"],
            "physical5_no_tau": ["H0", "omega_b", "omega_cdm", "lnAs", "n_s"],
            "full6": PARAMS,
        }
        baseline = {name: project(sn2, F, b, pars, 1.0e-10)
                    for name, pars in subsets.items()}
        sensitivity = {
            f"rcond_{r:.0e}": {
                "physical5_no_tau": project(
                    sn2, F, b,
                    ["H0", "omega_b", "omega_cdm", "lnAs", "n_s"], r
                ),
                "full6": project(sn2, F, b, PARAMS, r),
            }
            for r in RCONDS
        }
        single_absorption = {
            p: project(sn2, F, b, [p], 1.0e-10)["absorbed_delta_chi2_fraction"]
            for p in PARAMS
        }
        result[str(kmax)] = {
            "common_support": points,
            "baseline_projections": baseline,
            "pseudoinverse_sensitivity": sensitivity,
            "single_parameter_absorbed_fraction": single_absorption,
        }
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", type=Path, required=True)
    ap.add_argument("--optimized", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    out = {
        "analysis": "adversarial LSS Fisher projection audit using existing CLASS spectra",
        "new_CLASS_runs": False,
        "volume_Gpch3": VOLUME_GPCH3,
        "support_policy": "same intersection of finite support across signal, fiducial and all six nuisance derivatives for every nuisance subset",
        "parameter_order": PARAMS,
        "frozen": audit_signal(
            args.screen, args.screen, "frozen_plus", "frozen_minus"
        ),
        "optimized": audit_signal(
            args.screen, args.optimized, "optimized_plus", "optimized_minus"
        ),
    }

    # Reference result must reproduce the published recovery summary.
    ref = out["optimized"]["0.3"]["baseline_projections"]["full6"]["marginalized_SN_V50"]
    out["reference_check"] = {
        "expected": 0.06521388853657625,
        "recomputed": float(ref),
        "absolute_difference": float(abs(ref - 0.06521388853657625)),
        "pass_1e-6": bool(abs(ref - 0.06521388853657625) < 1.0e-6),
    }
    if not out["reference_check"]["pass_1e-6"]:
        raise RuntimeError(f"Reference reproduction failed: {out['reference_check']}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
