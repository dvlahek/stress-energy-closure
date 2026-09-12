#!/usr/bin/env python3
"""Scalar/lensing control with fixed present-day relic density.

This is the strict abundance control for the hidden-state mass sweep.  For each
mass, deg_ncdm is chosen so that the present-day physical density of the ncdm
component equals the m=0.10 eV reference value.  omega_cdm is left at the
baseline value.  Therefore the present-day relic gravitational weight and the
total matter density are both fixed while m/T and the velocity map vary.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import class_scalar_control_run as ctrl
from class_response_optimize import OMEGA_CDM, CLASS_COMMIT

CONTROL = "fixed-relic-omega0"
REF_MASS_EV = ctrl.REF_MASS_EV


def control_parameters(mass_eV: float, z_match: float) -> dict:
    omega_ref = ctrl.natural_omega_ncdm(REF_MASS_EV)
    omega_nat = ctrl.natural_omega_ncdm(mass_eV)
    deg = omega_ref / omega_nat
    if deg <= 0:
        raise RuntimeError(f"Unphysical deg_ncdm={deg}")
    return {
        "control": CONTROL,
        "reference_mass_eV": REF_MASS_EV,
        "mass_eV": mass_eV,
        "z_match": z_match,
        "omega_cdm": float(OMEGA_CDM),
        "deg_ncdm": float(deg),
        "natural_omega_ncdm_deg1": float(omega_nat),
        "reference_natural_omega_ncdm_deg1": float(omega_ref),
        "controlled_omega_ncdm_z0": float(deg * omega_nat),
        "invariant": "omega_ncdm at z=0",
        "invariant_reference_value": float(omega_ref),
        "class_omega_anchor_mass_eV": ctrl.CLASS_OMEGA_ANCHOR_MASS_EV,
        "class_omega_anchor": ctrl.CLASS_OMEGA_ANCHOR,
        "class_commit": CLASS_COMMIT,
    }


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
    ctrl.install_writer(params)

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
