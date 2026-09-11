#!/usr/bin/env python3
"""Control-aware wrapper for scalar/lensing hidden-state response optimization.

Two controls are supported without changing the validated baseline implementation:

1. fixed-total-matter
   Keep deg_ncdm=1 and compensate the natural present-day ncdm density by
   omega_cdm so that omega_cdm + omega_ncdm equals the m=0.10 eV reference.

2. fixed-relic-rho-match
   Keep omega_cdm fixed and rescale deg_ncdm so that the Fermi-Dirac relic
   energy density at z_match equals the m=0.10 eV reference. This holds the
   gravitational weight of the hidden component fixed at the matching epoch
   while varying m/T and the momentum-dependent velocity map.

The natural ncdm density is normalized to the pinned CLASS default value
omega_ncdm=0.0006451439 for m=0.06 eV and T_ncdm=0.71611. Ratios are computed
from the same Fermi-Dirac phase-space integral used by the kinetic setup.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import class_scalar_response_optimize as base
from class_response_optimize import (
    H0, OMEGA_B, OMEGA_CDM, A_S, N_UR, T_NCDM, TCMB_K, KB_EV_K,
    CLASS_COMMIT, f0_class, trap,
)

N_S = 0.9649
TAU_REIO = 0.054
REF_MASS_EV = 0.10
CLASS_OMEGA_ANCHOR_MASS_EV = 0.06
CLASS_OMEGA_ANCHOR = 0.0006451439
CONTROL_CHOICES = ("fixed-total-matter", "fixed-relic-rho-match")


def rho_integral(mass_eV: float, z: float) -> float:
    q = np.linspace(0.0, 20.0, 20000)
    f0 = f0_class(q)
    a = 1.0 / (1.0 + z)
    tnu0_eV = T_NCDM * TCMB_K * KB_EV_K
    y = a * mass_eV / tnu0_eV
    eps = np.sqrt(q * q + y * y)
    return float(trap(q**2 * eps * f0, q))


def natural_omega_ncdm(mass_eV: float) -> float:
    """Present-day omega_ncdm for deg=1, calibrated to pinned CLASS."""
    num = rho_integral(mass_eV, 0.0)
    den = rho_integral(CLASS_OMEGA_ANCHOR_MASS_EV, 0.0)
    return CLASS_OMEGA_ANCHOR * num / den


def control_parameters(control: str, mass_eV: float, z_match: float) -> dict:
    if control not in CONTROL_CHOICES:
        raise ValueError(control)

    omega_ref = natural_omega_ncdm(REF_MASS_EV)
    omega_nat = natural_omega_ncdm(mass_eV)

    if control == "fixed-total-matter":
        deg = 1.0
        omega_cdm = OMEGA_CDM + omega_ref - omega_nat
        invariant = "omega_cdm + natural omega_ncdm at z=0"
        invariant_ref = OMEGA_CDM + omega_ref
    else:
        rho_ref = rho_integral(REF_MASS_EV, z_match)
        rho_mass = rho_integral(mass_eV, z_match)
        deg = rho_ref / rho_mass
        omega_cdm = OMEGA_CDM
        invariant = "FD relic rho at z_match"
        invariant_ref = rho_ref

    if omega_cdm <= 0 or deg <= 0:
        raise RuntimeError(
            f"Unphysical control parameters: omega_cdm={omega_cdm}, deg={deg}"
        )

    return {
        "control": control,
        "reference_mass_eV": REF_MASS_EV,
        "mass_eV": mass_eV,
        "z_match": z_match,
        "omega_cdm": float(omega_cdm),
        "deg_ncdm": float(deg),
        "natural_omega_ncdm_deg1": float(omega_nat),
        "reference_natural_omega_ncdm_deg1": float(omega_ref),
        "invariant": invariant,
        "invariant_reference_value": float(invariant_ref),
        "class_omega_anchor_mass_eV": CLASS_OMEGA_ANCHOR_MASS_EV,
        "class_omega_anchor": CLASS_OMEGA_ANCHOR,
        "class_commit": CLASS_COMMIT,
    }


def install_writer(params: dict) -> None:
    omega_cdm = float(params["omega_cdm"])
    deg_ncdm = float(params["deg_ncdm"])
    control = params["control"]

    def write_scalar_ini(path: Path, psd: Path, root: Path, mass_eV: float,
                         lmax: int, pkmax: float) -> None:
        text = f"""# Scalar/lensing hidden-state control: {control}.
output = tCl,pCl,lCl,mPk
modes = s
lensing = yes
H0 = {H0}
omega_b = {OMEGA_B}
omega_cdm = {omega_cdm:.15g}
A_s = {A_S:.12e}
n_s = {N_S}
tau_reio = {TAU_REIO}
N_ur = {N_UR}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd.resolve()}
m_ncdm = {mass_eV}
T_ncdm = {T_NCDM}
deg_ncdm = {deg_ncdm:.15g}
l_max_scalars = {lmax}
P_k_max_h/Mpc = {pkmax}
z_pk = 0
root = {root}
headers = yes
write warnings = yes
"""
        path.write_text(text, encoding="utf-8")

    # The validated implementation resolves this global at call time, so the
    # wrapper can change only the background normalization while reusing all
    # null-space construction, optimization and metric code unchanged.
    base.write_scalar_ini = write_scalar_ini


def stamp_metadata(outdir: Path, params: dict) -> None:
    p = outdir / "scalar_probe_meta.json"
    if p.exists():
        meta = json.loads(p.read_text())
        meta["control_parameters"] = params
        p.write_text(json.dumps(meta, indent=2) + "\n")
    (outdir / "control_parameters.json").write_text(
        json.dumps(params, indent=2) + "\n"
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=[
        "prepare-probes", "optimize", "summarize-candidates",
        "prepare-winner-highprec", "summarize-winner-highprec",
        "print-control-parameters",
    ])
    ap.add_argument("--control", required=True, choices=CONTROL_CHOICES)
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

    params = control_parameters(args.control, args.mass_eV, args.z_match)
    install_writer(params)

    if args.mode == "print-control-parameters":
        print(json.dumps(params, indent=2))
        return

    if args.mode == "prepare-probes":
        base.prepare_probes(
            args.outdir, args.mass_eV, args.z_match,
            args.probe_frac, args.lmax, args.pkmax,
        )
        stamp_metadata(args.outdir, params)
    elif args.mode == "optimize":
        base.optimize(
            args.outdir, args.final_frac, args.lmin, args.lmax,
            args.random_samples, args.seed,
        )
        stamp_metadata(args.outdir, params)
    elif args.mode == "summarize-candidates":
        base.summarize_candidates(args.outdir, args.lmin, args.lmax, args.pkmax)
        stamp_metadata(args.outdir, params)
    elif args.mode == "prepare-winner-highprec":
        base.prepare_winner_highprec(
            args.outdir, args.mass_eV, args.lmax, args.pkmax,
        )
        stamp_metadata(args.outdir, params)
    else:
        base.summarize_winner_highprec(
            args.outdir, args.lmin, args.lmax, args.pkmax,
        )
        stamp_metadata(args.outdir, params)


if __name__ == "__main__":
    main()
