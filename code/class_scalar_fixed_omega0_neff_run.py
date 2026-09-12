#!/usr/bin/env python3
"""Scalar/lensing control with fixed present-day relic density and early radiation.

For each mass, deg_ncdm is first chosen exactly as in the strict fixed-omega0
control so that omega_ncdm(z=0) equals the m=0.10 eV reference value.  The
ultra-relativistic component is then adjusted so that the early-time total
neutrino effective number is the same as in the reference setup.

The reference setup uses T_ncdm=0.71611, deg_ncdm=1 and N_ur=2.0328.  At this
temperature one thermal ncdm family contributes 1.0132 to early N_eff, hence
the reference total is 3.046.  With arbitrary degeneracy d we therefore use

    N_ur(d) = 3.046 - 1.0132 d.

This control is intended only for masses where the required N_ur is nonnegative.
The tested grid m=0.10,0.18,0.30,0.60 eV satisfies this condition.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import class_scalar_control_run as ctrl
import class_scalar_fixed_omega0_run as omega0
from class_response_optimize import N_UR

CONTROL = "fixed-relic-omega0-fixed-early-neff"
REF_MASS_EV = omega0.REF_MASS_EV
EARLY_NEFF_TARGET = 3.046
NCDM_NEFF_PER_DEG = EARLY_NEFF_TARGET - N_UR  # 1.0132 for T_ncdm=0.71611


def control_parameters(mass_eV: float, z_match: float) -> dict:
    p = omega0.control_parameters(mass_eV, z_match)
    deg = float(p["deg_ncdm"])
    n_ur = EARLY_NEFF_TARGET - NCDM_NEFF_PER_DEG * deg
    if n_ur < 0:
        raise RuntimeError(
            f"Cannot keep early N_eff={EARLY_NEFF_TARGET} with deg_ncdm={deg}: "
            f"required N_ur={n_ur} < 0"
        )
    p.update({
        "control": CONTROL,
        "N_ur": float(n_ur),
        "reference_N_ur": float(N_UR),
        "early_Neff_target": float(EARLY_NEFF_TARGET),
        "ncdm_early_Neff_per_deg": float(NCDM_NEFF_PER_DEG),
        "controlled_early_Neff": float(n_ur + NCDM_NEFF_PER_DEG * deg),
        "invariants": [
            "omega_ncdm at z=0",
            "early-time total neutrino N_eff for the thermal baseline normalization",
        ],
        "interpretation": (
            "Present-day relic abundance and the baseline early radiation density are fixed; "
            "the intermediate-time massive transition and momentum-dependent velocity response vary."
        ),
    })
    return p


def install_writer(params: dict) -> None:
    # Reuse the validated control writer, changing only the value of N_ur that
    # it resolves when constructing each CLASS input file.
    ctrl.N_UR = float(params["N_ur"])
    ctrl.install_writer(params)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=[
        "prepare-probes", "optimize", "summarize-candidates",
        "prepare-winner-highprec", "summarize-winner-highprec",
        "print-control-parameters",
    ])
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mass-eV", type=float, default=REF_MASS_EV)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.08)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--lmin", type=int, default=30)
    ap.add_argument("--lmax", type=int, default=600)
    ap.add_argument("--pkmax", type=float, default=0.30)
    ap.add_argument("--random-samples", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=314159)
    args = ap.parse_args()

    params = control_parameters(args.mass_eV, args.z_match)
    install_writer(params)

    if args.mode == "print-control-parameters":
        print(json.dumps(params, indent=2))
        return

    if args.mode == "prepare-probes":
        ctrl.base.prepare_probes(
            args.outdir, args.mass_eV, args.z_match,
            args.probe_frac, args.lmax, args.pkmax,
        )
        ctrl.stamp_metadata(args.outdir, params)
    elif args.mode == "optimize":
        ctrl.base.optimize(
            args.outdir, args.final_frac, args.lmin, args.lmax,
            args.random_samples, args.seed,
        )
        ctrl.stamp_metadata(args.outdir, params)
    elif args.mode == "summarize-candidates":
        ctrl.base.summarize_candidates(
            args.outdir, args.lmin, args.lmax, args.pkmax
        )
        ctrl.stamp_metadata(args.outdir, params)
    elif args.mode == "prepare-winner-highprec":
        ctrl.base.prepare_winner_highprec(
            args.outdir, args.mass_eV, args.lmax, args.pkmax,
        )
        ctrl.stamp_metadata(args.outdir, params)
    else:
        ctrl.base.summarize_winner_highprec(
            args.outdir, args.lmin, args.lmax, args.pkmax,
        )
        ctrl.stamp_metadata(args.outdir, params)


if __name__ == "__main__":
    main()
