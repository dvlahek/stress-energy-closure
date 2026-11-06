#!/usr/bin/env python3
"""K01B-R1: one preregistered higher-N numerical confirmation.

The original K01B result remains FAIL. R1 changes only particle count by the
predeclared factor of four and applies the same local dynamics and unchanged
resolution gates.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import k01b_remote_local_nfw_matched_particle_dynamics as k

PREREG = ROOT / "source_data/k01b_r1_higher_n_resolution_confirmation_prereg_2026-11-05.json"
PREREG_BLOB = "4bb4cefd78d261e498ab11e88f87623c26fa55b3"
OUT = ROOT / "k01b_r1_ci_result.json"


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01B_R1_STOP: " + msg)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "R1 preregistration changed")
    p = json.loads(raw)
    need(
        p["stage"] == "K01B_R1_HIGHER_N_NUMERICAL_RESOLUTION_CONFIRMATION_PREREG"
        and p["only_resolution_change"]["particles_per_state"] == 262144
        and p["only_resolution_change"]["nested_half_sample"] == 131072
        and p["only_resolution_change"]["factor_vs_parent"] == 4
        and p["immutable_physics"]["physical_periodic_box_Mpc"] == 4
        and p["immutable_physics"]["position_direction_seed"] == 2026110502
        and p["immutable_physics"]["coarse_steps"] == 128
        and p["immutable_physics"]["fine_steps"] == 256
        and p["immutable_physics"]["halo_relative_wind_km_s"] == 200,
        "R1 declared design drift",
    )
    return p


def main() -> None:
    p = load_prereg()

    # Reuse the same frozen source and halo gates as original K01B.
    k.e8.source_gate()
    src = k.e8.kinetic_source()
    _, base, profiles = k.e28.sources()
    _, _, backgrounds, halos = base
    profile = profiles[0]
    mass_msun = float(halos[0][0])
    c200 = float(profile["c200c_population_median_DM14"])
    rs_comoving = float(profile["r_s_comoving_Mpc_DM14"])
    rs_phys = rs_comoving / (1.0 + k.e27.ZOBS)
    duration = k.proper_interval_gyr(backgrounds["FD"])

    need(
        math.isclose(mass_msun, float(profile["M200c_physical_Msun_DM14"]), rel_tol=1e-13),
        "same frozen halo mismatch",
    )

    # Translate R1 preregistration into the exact original K01B phase-space helper schema.
    phase_prereg = {
        "particle_design": {
            "particles_per_state": int(p["only_resolution_change"]["particles_per_state"]),
            "position_direction_seed": int(p["immutable_physics"]["position_direction_seed"]),
            "physical_periodic_box_Mpc": float(p["immutable_physics"]["physical_periodic_box_Mpc"]),
        }
    }

    pos, dirs, qs = k.build_matched_phase_space(src, phase_prereg)
    need(len(pos) == 262144, "R1 particle count")

    box = float(p["immutable_physics"]["physical_periodic_box_Mpc"])
    ncoarse = int(p["immutable_physics"]["coarse_steps"])
    nfine = int(p["immutable_physics"]["fine_steps"])

    fine = k.run_wind(pos, dirs, qs, +1, duration, nfine, box, mass_msun, c200, rs_phys)
    coarse = k.run_wind(pos, dirs, qs, +1, duration, ncoarse, box, mass_msun, c200, rs_phys)
    reverse = k.run_wind(pos, dirs, qs, -1, duration, nfine, box, mass_msun, c200, rs_phys)
    free = k.run_wind(pos, dirs, qs, +1, duration, 8, box, 0.0, c200, rs_phys)

    free_max = max(abs(free[s]["D_full_km_s"]) for s in ("FD", "Fplus", "Fminus"))
    need(free_max == 0.0, "R1 no-halo control")

    C = fine["contrast"]["C_full_km_s"]
    Cc = coarse["contrast"]["C_full_km_s"]
    Ch = fine["contrast"]["C_half_km_s"]
    Cr = reverse["contrast"]["C_full_km_s"]

    tiny = 1e-300
    timestep_gap = abs(C - Cc) / max(abs(C), tiny)
    nested_gap = abs(C - Ch) / max(abs(C), tiny)
    floor = max(abs(C - Cc), abs(C - Ch))
    contrast_to_floor = abs(C) / max(floor, tiny)

    parity_parts = []
    for s in ("FD", "Fplus", "Fminus"):
        a = fine[s]["D_full_km_s"]
        b = reverse[s]["D_full_km_s"]
        parity_parts.append(abs(a + b) / max(abs(a), abs(b), tiny))
    parity_parts.append(abs(C + Cr) / max(abs(C), abs(Cr), tiny))
    parity = max(parity_parts)

    sign_consistent = (
        C != 0 and Cc != 0 and Ch != 0
        and np.sign(C) == np.sign(Cc) == np.sign(Ch)
    )

    g = p["unchanged_gate"]
    gate_pass = bool(
        sign_consistent
        and timestep_gap <= float(g["timestep_relative_gap_at_most"])
        and nested_gap <= float(g["nested_sample_relative_gap_at_most"])
        and contrast_to_floor >= float(g["contrast_to_numerical_floor_at_least"])
        and parity <= float(g["wind_reversal_relative_parity_residual_at_most"])
    )

    Dfd = fine["FD"]["D_full_km_s"]
    rel = C / abs(Dfd) if Dfd != 0 else None
    parent_C = float(p["parent"]["parent_fine_contrast_km_s"])
    cross_parent_relative_difference = abs(C - parent_C) / max(abs(C), abs(parent_C), tiny)

    result = {
        "stage": p["stage"],
        "status": (
            "PASS_HIGHER_N_CONFIRMATION_RESOLVES_LOCAL_HIDDEN_STATE_DYNAMICS"
            if gate_pass
            else "FAIL_HIGHER_N_CONFIRMATION_UNRESOLVED_STOP_LOCAL_PARTICLE_CHASING"
        ),
        "science_gate_pass": gate_pass,
        "original_K01B_status_remains": p["parent"]["result"],
        "interpretation": p["interpretation"],
        "stopping_rule": p["stopping_rule"],
        "preregistration_git_blob": PREREG_BLOB,
        "frozen_inputs": {
            "duration_Gyr": duration,
            "halo_mass_Msun": mass_msun,
            "c200c": c200,
            "r_s_physical_Mpc": rs_phys,
            "particles_per_state": int(len(pos)),
            "nested_half_sample": int(len(pos) // 2),
            "box_physical_Mpc": box,
            "coarse_steps": ncoarse,
            "fine_steps": nfine,
        },
        "fine_plus_wind": fine,
        "coarse_plus_wind": coarse,
        "fine_reversed_wind": reverse,
        "controls": {
            "no_halo_max_abs_D_km_s": free_max,
            "wind_reversal_component_relative_residuals": parity_parts,
            "wind_reversal_max_relative_residual": parity,
        },
        "hidden_contrast": {
            "C_fine_full_km_s": C,
            "C_coarse_full_km_s": Cc,
            "C_fine_half_km_s": Ch,
            "C_reversed_fine_full_km_s": Cr,
            "relative_to_abs_FD_drag": rel,
            "timestep_relative_gap": timestep_gap,
            "nested_sample_relative_gap": nested_gap,
            "numerical_floor_km_s": floor,
            "contrast_to_numerical_floor": contrast_to_floor,
            "sign_consistent": bool(sign_consistent),
            "parent_K01B_C_fine_km_s": parent_C,
            "cross_parent_relative_difference": cross_parent_relative_difference,
        },
        "guards": p["guards"],
    }

    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print("K01B_R1_PREREG_PIN_PASS", PREREG_BLOB, flush=True)
    print("K01B_R1_PARTICLES", len(pos), flush=True)
    print("K01B_R1_FD_D_KMS", Dfd, flush=True)
    print("K01B_R1_HIDDEN_C_FINE_KMS", C, flush=True)
    print("K01B_R1_HIDDEN_RELATIVE_TO_FD", rel, flush=True)
    print("K01B_R1_TIMESTEP_GAP", timestep_gap, flush=True)
    print("K01B_R1_NESTED_GAP", nested_gap, flush=True)
    print("K01B_R1_CONTRAST_TO_NUMERICAL_FLOOR", contrast_to_floor, flush=True)
    print("K01B_R1_WIND_REVERSAL_MAX_RESIDUAL", parity, flush=True)
    print("K01B_R1_CROSS_PARENT_RELATIVE_DIFFERENCE", cross_parent_relative_difference, flush=True)
    print("K01B_R1_SCIENCE_GATE_PASS", gate_pass, flush=True)
    print("K01B_R1_STATUS", result["status"], flush=True)
    print("K01B_R1_RESULT", OUT, flush=True)


if __name__ == "__main__":
    main()
