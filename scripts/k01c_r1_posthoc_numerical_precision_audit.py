#!/usr/bin/env python3
"""K01C-R1 posthoc numerical precision audit.

The parent K01C result is preserved. This audit increases only numerical
quadrature precision for the same deterministic local scattering estimand and
adds an adaptive point-mass benchmark against the exact Kepler/Rutherford
deflection formula.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import k01c_deterministic_local_nfw_phase_space_scattering as k01c
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as e27
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as e28

PREREG = ROOT / "source_data/k01c_r1_posthoc_numerical_precision_audit_prereg_2026-11-05.json"
PREREG_BLOB = "6d44e3a9cd9c147e9089a30898aa97e2042d19b8"
OUT = ROOT / "k01c_r1_precision_audit_result.json"


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01C_R1_STOP: " + msg)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "preregistration changed")
    p = json.loads(raw)
    hp = p["high_precision_replay"]
    need(
        p["stage"] == "K01C_R1_POSTHOC_NUMERICAL_PRECISION_AUDIT_PREREG"
        and hp["q_legendre_nodes_per_interval"] == 12
        and hp["mu_legendre_nodes"] == 32
        and hp["b_legendre_nodes_0_to_r200"] == 64
        and hp["orbit_angle_legendre_nodes"] == 160,
        "precision-audit design drift",
    )
    return p


def adaptive_point_mass_coschi(v: float, b: float, M: float):
    """Independent numerical orbit-angle integral for a Kepler point mass."""
    A = e27.G * M / (v * v)
    rmin = math.sqrt(A * A + b * b) - A
    h = b * v

    def phi(r: float) -> float:
        return -e27.G * M / r

    def integrand(t: float) -> float:
        # u=rmin/r=1-t^2 smooths the square-root turning-point singularity.
        if t == 0.0:
            t = 1e-12
        u = 1.0 - t * t
        r = rmin / u
        rad = v * v - 2.0 * phi(r) - h * h / (r * r)
        if rad <= 0.0:
            # Only removable endpoint roundoff is allowed.
            rad = np.finfo(float).tiny
        return 2.0 * t * (h / rmin) / math.sqrt(rad)

    phi0, err = quad(
        integrand,
        0.0,
        1.0,
        epsabs=2e-13,
        epsrel=2e-13,
        limit=500,
        points=[0.0],
    )
    chi = 2.0 * phi0 - math.pi
    return math.cos(chi), chi, err


def exact_point_mass(v: float, b: float, M: float):
    chi = 2.0 * math.atan(e27.G * M / (b * v * v))
    return math.cos(chi), chi


def relgap(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main() -> None:
    p = load_prereg()

    # Frozen source/halo setup exactly as parent K01C.
    k01c.e8.source_gate()
    src = k01c.e8.kinetic_source()
    _, base, profiles = e28.sources()
    _, _, _, halos = base
    prof = profiles[0]

    M = float(halos[0][0])
    c = float(prof["c200c_population_median_DM14"])
    rs = float(prof["r_s_comoving_Mpc_DM14"]) / (1.0 + e27.ZOBS)
    r200 = c * rs

    need(
        math.isclose(M, float(prof["M200c_physical_Msun_DM14"]), rel_tol=1e-13),
        "frozen halo mismatch",
    )

    # Independent adaptive Kepler benchmark.
    qa = p["independent_point_mass_adaptive_QA"]
    cases = []
    max_cos_gap = 0.0
    max_chi_gap = 0.0
    max_quad_err = 0.0

    for v in qa["speeds_km_s"]:
        for bf in qa["b_over_r200"]:
            b = float(bf) * r200
            cn, an, err = adaptive_point_mass_coschi(float(v), b, M)
            ce, ae = exact_point_mass(float(v), b, M)
            cg = abs(cn - ce)
            ag = abs(an - ae)
            max_cos_gap = max(max_cos_gap, cg)
            max_chi_gap = max(max_chi_gap, ag)
            max_quad_err = max(max_quad_err, abs(err))
            cases.append({
                "v_km_s": float(v),
                "b_over_r200": float(bf),
                "numeric_coschi": cn,
                "exact_coschi": ce,
                "abs_coschi_gap": cg,
                "abs_chi_gap_rad": ag,
                "quad_reported_error": err,
            })

    point_pass = bool(
        max_cos_gap <= float(qa["max_abs_coschi_gap"])
        and max_chi_gap <= float(qa["max_abs_chi_gap_rad"])
    )

    # Higher-precision replay of the actual NFW integrated parent estimand.
    hp = p["high_precision_replay"]
    fine = k01c.evaluate_grid(
        src,
        M,
        c,
        rs,
        e27.VHALO,
        hp["q_intervals"],
        hp["q_legendre_nodes_per_interval"],
        hp["mu_legendre_nodes"],
        hp["b_legendre_nodes_0_to_r200"],
        hp["orbit_angle_legendre_nodes"],
    )

    parent = p["parent"]
    parent_vals = {
        "FD": float(parent["primary_FD_D"]),
        "Fplus": float(parent["primary_Fplus_D"]),
        "Fminus": float(parent["primary_Fminus_D"]),
    }
    key = "D_bmax_1r200"
    replay_vals = {s: float(fine[s][key]) for s in parent_vals}
    state_shifts = {
        s: relgap(replay_vals[s], parent_vals[s]) for s in parent_vals
    }

    C = float(fine["contrast"][key])
    parent_C = float(parent["primary_hidden_C"])
    C_shift = relgap(C, parent_C)

    keys = ("D_bmax_1r200", "D_bmax_2r200", "D_bmax_4r200")
    signs = [int(np.sign(fine["contrast"][kk])) for kk in keys]
    same_cutoff_sign = all(x != 0 for x in signs) and len(set(signs)) == 1
    same_parent_sign = C != 0 and parent_C != 0 and np.sign(C) == np.sign(parent_C)

    a = p["acceptance"]
    pass_gate = bool(
        point_pass
        and C_shift <= float(a["parent_primary_C_relative_shift_at_most"])
        and max(state_shifts.values()) <= float(a["each_parent_primary_state_D_relative_shift_at_most"])
        and same_parent_sign
        and same_cutoff_sign
    )

    result = {
        "stage": p["stage"],
        "status": (
            "PASS_K01C_HIGH_PRECISION_AND_ADAPTIVE_POINT_MASS_QA"
            if pass_gate
            else "FAIL_K01C_POSTHOC_PRECISION_AUDIT"
        ),
        "science_precision_gate_pass": pass_gate,
        "interpretation": p["interpretation"],
        "preregistration_git_blob": PREREG_BLOB,
        "adaptive_point_mass_QA": {
            "pass": point_pass,
            "max_abs_coschi_gap": max_cos_gap,
            "max_abs_chi_gap_rad": max_chi_gap,
            "max_quad_reported_error": max_quad_err,
            "cases": cases,
        },
        "high_precision_replay": fine,
        "parent_comparison": {
            "parent_primary_state_D": parent_vals,
            "replay_primary_state_D": replay_vals,
            "state_relative_shifts": state_shifts,
            "parent_primary_C": parent_C,
            "replay_primary_C": C,
            "primary_C_relative_shift": C_shift,
            "primary_C_sign_unchanged": bool(same_parent_sign),
            "replay_contrast_signs_bmax_1_2_4r200": signs,
            "contrast_sign_same_bmax_1_2_4r200": bool(same_cutoff_sign),
        },
        "guards": p["guards"],
    }

    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print("K01C_R1_PREREG_PIN_PASS", PREREG_BLOB, flush=True)
    print("K01C_R1_POINT_MAX_COS_GAP", max_cos_gap, flush=True)
    print("K01C_R1_POINT_MAX_CHI_GAP", max_chi_gap, flush=True)
    print("K01C_R1_POINT_QA_PASS", point_pass, flush=True)
    print("K01C_R1_REPLAY_FD", replay_vals["FD"], flush=True)
    print("K01C_R1_REPLAY_FPLUS", replay_vals["Fplus"], flush=True)
    print("K01C_R1_REPLAY_FMINUS", replay_vals["Fminus"], flush=True)
    print("K01C_R1_REPLAY_C", C, flush=True)
    print("K01C_R1_C_REL_SHIFT", C_shift, flush=True)
    print("K01C_R1_STATE_REL_SHIFTS", state_shifts, flush=True)
    print("K01C_R1_CUTOFF_SIGNS", signs, flush=True)
    print("K01C_R1_GATE_PASS", pass_gate, flush=True)
    print("K01C_R1_STATUS", result["status"], flush=True)


if __name__ == "__main__":
    main()
