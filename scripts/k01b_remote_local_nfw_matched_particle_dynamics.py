#!/usr/bin/env python3
"""K01B: remote local matched-particle hidden-state dynamics smoke test.

Directly evolves matched neutrino test particles through the same frozen NFW
halo potential for the original E8 F0/F+/F- momentum distributions. This is a
local frozen-z collisionless dynamics check only. It is not a cosmological
N-body simulation, galaxy/HOD model, survey observable, or detection forecast.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "source_data/k01b_remote_local_nfw_matched_particle_dynamics_prereg_2026-11-05.json"
PREREG_BLOB = "495e98a885a1ff69462f6c2b1e19b4d9b15b39e6"
OUT = ROOT / "k01b_ci_result.json"

sys.path.insert(0, str(ROOT / "scripts"))
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as e27
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as e28

E8_PATH = ROOT / "scripts/build_eboss_dr16_a03_e8_frozen_resonant_occupation.py"
spec = importlib.util.spec_from_file_location("k01b_e8_frozen", E8_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("Cannot load frozen E8 source builder")
e8 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = e8
spec.loader.exec_module(e8)

MPC_KM = 3.0856775814913673e19
SEC_GYR = 365.25 * 86400.0 * 1.0e9
KMS_TO_MPC_PER_GYR = SEC_GYR / MPC_KM
MPC_OVER_C_GYR = MPC_KM / e27.C / SEC_GYR


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01B_IMPLEMENTATION_STOP: " + msg)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "preregistration changed after lock")
    p = json.loads(raw)
    need(
        p["stage"] == "K01B_REMOTE_LOCAL_NFW_MATCHED_PARTICLE_DYNAMICS_PREREG"
        and p["particle_design"]["particles_per_state"] == 65536
        and p["particle_design"]["nested_half_sample"] == 32768
        and p["particle_design"]["physical_periodic_box_Mpc"] == 4
        and p["dynamics"]["primary_steps"] == 256
        and p["dynamics"]["coarse_steps"] == 128
        and p["frozen_inputs"]["halo_relative_wind_km_s"] == 200
        and p["frozen_inputs"]["z_obs"] == 0.95,
        "preregistered numerical design changed",
    )
    return p


def inverse_cdf_samples(q: np.ndarray, f: np.ndarray, u: np.ndarray) -> np.ndarray:
    pdf = q * q * f
    cdf = np.zeros_like(q)
    cdf[1:] = np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * (q[1:] - q[:-1]))
    need(np.isfinite(cdf[-1]) and cdf[-1] > 0, "invalid momentum CDF")
    cdf /= cdf[-1]
    keep = np.concatenate(([True], np.diff(cdf) > 0))
    s = np.interp(u, cdf[keep], q[keep])
    need(np.all(np.isfinite(s)) and np.all(s >= 0), "invalid momentum quantiles")
    return s


def build_matched_phase_space(src: dict, prereg: dict):
    n = int(prereg["particle_design"]["particles_per_state"])
    need(n % 2 == 0, "antithetic design requires even particle count")
    nbase = n // 2
    seed = int(prereg["particle_design"]["position_direction_seed"])
    box = float(prereg["particle_design"]["physical_periodic_box_Mpc"])

    rng = np.random.default_rng(seed)

    # Randomly permute the stratified quantiles so the preregistered first-half
    # nested sample remains representative rather than containing only low-q.
    u0 = (np.arange(nbase, dtype=np.float64) + 0.5) / nbase
    u = u0[rng.permutation(nbase)]

    pos0 = rng.uniform(-0.5 * box, 0.5 * box, size=(nbase, 3))
    d0 = rng.normal(size=(nbase, 3))
    d0 /= np.linalg.norm(d0, axis=1)[:, None]

    pos = np.empty((n, 3), dtype=np.float64)
    dirs = np.empty((n, 3), dtype=np.float64)
    pos[0::2] = pos0
    pos[1::2] = -pos0
    dirs[0::2] = d0
    dirs[1::2] = -d0

    qgrid = np.asarray(src["q"], dtype=np.float64)
    states = {
        "FD": np.asarray(src["f0"], dtype=np.float64),
        "Fplus": np.asarray(src["fp"], dtype=np.float64),
        "Fminus": np.asarray(src["fm"], dtype=np.float64),
    }

    qs = {}
    for name, f in states.items():
        qb = inverse_cdf_samples(qgrid, f, u)
        qfull = np.empty(n, dtype=np.float64)
        qfull[0::2] = qb
        qfull[1::2] = qb
        qs[name] = qfull

    return pos, dirs, qs


def proper_interval_gyr(background) -> float:
    zz, hh = background
    z = np.linspace(e27.ZOBS, e27.ZINIT, 2049, dtype=np.float64)
    H = np.interp(z, zz, hh)
    need(np.all(H > 0), "invalid frozen CLASS H(z)")
    dt_mpc = float(np.trapezoid(1.0 / ((1.0 + z) * H), z))
    dt_gyr = dt_mpc * MPC_OVER_C_GYR
    need(0 < dt_gyr < 2.0, "unexpected z=1 to z=.95 proper-time interval")
    return dt_gyr


def thermal_velocity_kms(qs: np.ndarray) -> np.ndarray:
    p_ev = qs * e27.TNU0 * (1.0 + e27.ZOBS)
    return e27.C * p_ev / np.sqrt(p_ev * p_ev + e27.MNU * e27.MNU)


def wrap(pos: np.ndarray, box: float) -> None:
    pos[:] = (pos + 0.5 * box) % box - 0.5 * box


def nfw_accel_kms_per_gyr(
    pos: np.ndarray,
    mass_msun: float,
    c200: float,
    rs_phys_mpc: float,
) -> np.ndarray:
    if mass_msun == 0:
        return np.zeros_like(pos)

    r2 = np.einsum("ij,ij->i", pos, pos)
    r = np.sqrt(r2)
    r_safe = np.maximum(r, rs_phys_mpc * 1e-8)
    x = r_safe / rs_phys_mpc
    fc = math.log1p(c200) - c200 / (1.0 + c200)
    fx = np.log1p(x) - x / (1.0 + x)
    frac = np.where(x < c200, fx / fc, 1.0)
    menc = mass_msun * frac

    invr3 = 1.0 / (r_safe * r_safe * r_safe)
    a_native = -e27.G * menc[:, None] * pos * invr3[:, None]
    a_native[r == 0] = 0.0
    return a_native * KMS_TO_MPC_PER_GYR


def evolve(
    pos_initial: np.ndarray,
    vel_initial: np.ndarray,
    duration_gyr: float,
    steps: int,
    box: float,
    mass_msun: float,
    c200: float,
    rs_phys_mpc: float,
):
    pos = pos_initial.copy()
    vel = vel_initial.copy()
    dt = duration_gyr / steps

    for _ in range(steps):
        vel += 0.5 * dt * nfw_accel_kms_per_gyr(pos, mass_msun, c200, rs_phys_mpc)
        pos += dt * KMS_TO_MPC_PER_GYR * vel
        wrap(pos, box)
        vel += 0.5 * dt * nfw_accel_kms_per_gyr(pos, mass_msun, c200, rs_phys_mpc)

    return pos, vel


def run_wind(
    pos: np.ndarray,
    dirs: np.ndarray,
    qs: dict,
    wind_sign: int,
    duration_gyr: float,
    steps: int,
    box: float,
    mass_msun: float,
    c200: float,
    rs_phys_mpc: float,
):
    need(wind_sign in (-1, 1), "wind sign")
    V = e27.VHALO * wind_sign
    out = {}

    for name in ("FD", "Fplus", "Fminus"):
        speed = thermal_velocity_kms(qs[name])
        v0 = speed[:, None] * dirs
        # +wind_sign means halo moves +x, so the neutrino bath has mean -V in halo frame.
        v0[:, 0] -= V
        _, vf = evolve(pos, v0, duration_gyr, steps, box, mass_msun, c200, rs_phys_mpc)
        dvx = vf[:, 0] - v0[:, 0]
        out[name] = {
            "D_full_km_s": float(np.mean(dvx)),
            "D_half_km_s": float(np.mean(dvx[: len(dvx) // 2])),
            "mean_initial_vx_km_s": float(np.mean(v0[:, 0])),
            "mean_final_vx_km_s": float(np.mean(vf[:, 0])),
        }

    out["contrast"] = {
        "C_full_km_s": out["Fminus"]["D_full_km_s"] - out["Fplus"]["D_full_km_s"],
        "C_half_km_s": out["Fminus"]["D_half_km_s"] - out["Fplus"]["D_half_km_s"],
    }
    return out


def main() -> None:
    prereg = load_prereg()

    e8.source_gate()
    src = e8.kinetic_source()

    # Existing E28 source gate pins the original E27 halo benchmark and the
    # DM14/NFW low/high profile anchors. We use only the preregistered lower one.
    _, base, profiles = e28.sources()
    _, _, backgrounds, halos = base
    profile = profiles[0]
    mass_msun = float(halos[0][0])
    c200 = float(profile["c200c_population_median_DM14"])
    rs_comoving = float(profile["r_s_comoving_Mpc_DM14"])
    rs_phys = rs_comoving / (1.0 + e27.ZOBS)

    need(
        math.isclose(mass_msun, float(profile["M200c_physical_Msun_DM14"]), rel_tol=1e-13),
        "frozen E27 mass and E28 NFW profile mismatch",
    )

    duration = proper_interval_gyr(backgrounds["FD"])
    box = float(prereg["particle_design"]["physical_periodic_box_Mpc"])
    ncoarse = int(prereg["dynamics"]["coarse_steps"])
    nfine = int(prereg["dynamics"]["primary_steps"])

    pos, dirs, qs = build_matched_phase_space(src, prereg)

    fine = run_wind(
        pos, dirs, qs, +1, duration, nfine, box, mass_msun, c200, rs_phys
    )
    coarse = run_wind(
        pos, dirs, qs, +1, duration, ncoarse, box, mass_msun, c200, rs_phys
    )
    reverse = run_wind(
        pos, dirs, qs, -1, duration, nfine, box, mass_msun, c200, rs_phys
    )

    # No-halo control. Since acceleration is identically zero, velocity must
    # remain bitwise unchanged modulo ordinary floating-point copies.
    free = run_wind(pos, dirs, qs, +1, duration, 8, box, 0.0, c200, rs_phys)
    free_max = max(abs(free[s]["D_full_km_s"]) for s in ("FD", "Fplus", "Fminus"))
    need(free_max == 0.0, "no-halo control produced nonzero momentum transfer")

    C = fine["contrast"]["C_full_km_s"]
    Cc = coarse["contrast"]["C_full_km_s"]
    Ch = fine["contrast"]["C_half_km_s"]
    Cr = reverse["contrast"]["C_full_km_s"]

    tiny = 1e-300
    timestep_gap = abs(C - Cc) / max(abs(C), tiny)
    nested_gap = abs(C - Ch) / max(abs(C), tiny)
    numerical_floor = max(abs(C - Cc), abs(C - Ch))
    contrast_to_floor = abs(C) / max(numerical_floor, tiny)

    parity_components = []
    for s in ("FD", "Fplus", "Fminus"):
        a = fine[s]["D_full_km_s"]
        b = reverse[s]["D_full_km_s"]
        parity_components.append(abs(a + b) / max(abs(a), abs(b), tiny))
    parity_components.append(abs(C + Cr) / max(abs(C), abs(Cr), tiny))
    parity_residual = max(parity_components)

    sign_consistent = (
        C != 0.0
        and Cc != 0.0
        and Ch != 0.0
        and np.sign(C) == np.sign(Cc) == np.sign(Ch)
    )

    gate = prereg["preregistered_resolution_gate"]
    science_gate_pass = bool(
        sign_consistent
        and timestep_gap <= float(gate["timestep_relative_gap_at_most"])
        and nested_gap <= float(gate["nested_sample_relative_gap_at_most"])
        and contrast_to_floor >= float(gate["contrast_to_numerical_floor_at_least"])
        and parity_residual <= float(gate["wind_reversal_relative_parity_residual_at_most"])
    )

    Dfd = fine["FD"]["D_full_km_s"]
    relative_contrast = C / abs(Dfd) if Dfd != 0 else None

    result = {
        "stage": "K01B_REMOTE_LOCAL_NFW_MATCHED_PARTICLE_DYNAMICS",
        "status": (
            "PASS_NUMERICALLY_RESOLVED_LOCAL_HIDDEN_STATE_DYNAMICS"
            if science_gate_pass
            else "FAIL_LOCAL_HIDDEN_STATE_CONTRAST_NOT_RESOLVED_UNDER_PREREGISTERED_GATE"
        ),
        "science_gate_pass": science_gate_pass,
        "interpretation": (
            "Direct collisionless test-particle dynamics in one frozen static periodic NFW halo benchmark. "
            "A pass means only that F+/F- hidden-state dependence is numerically resolved in this local setup."
        ),
        "preregistration_git_blob": PREREG_BLOB,
        "frozen_inputs": {
            "z_obs": e27.ZOBS,
            "z_init_for_duration": e27.ZINIT,
            "duration_Gyr": duration,
            "halo_mass_Msun": mass_msun,
            "c200c": c200,
            "r_s_comoving_Mpc": rs_comoving,
            "r_s_physical_Mpc": rs_phys,
            "halo_wind_km_s": e27.VHALO,
            "box_physical_Mpc": box,
            "particles_per_state": int(len(pos)),
            "coarse_steps": ncoarse,
            "fine_steps": nfine,
        },
        "fine_plus_wind": fine,
        "coarse_plus_wind": coarse,
        "fine_reversed_wind": reverse,
        "controls": {
            "no_halo_max_abs_D_km_s": free_max,
            "wind_reversal_component_relative_residuals": parity_components,
            "wind_reversal_max_relative_residual": parity_residual,
        },
        "hidden_contrast": {
            "C_fine_full_km_s": C,
            "C_coarse_full_km_s": Cc,
            "C_fine_half_km_s": Ch,
            "C_reversed_fine_full_km_s": Cr,
            "relative_to_abs_FD_drag": relative_contrast,
            "timestep_relative_gap": timestep_gap,
            "nested_sample_relative_gap": nested_gap,
            "numerical_floor_km_s": numerical_floor,
            "contrast_to_numerical_floor": contrast_to_floor,
            "sign_consistent": bool(sign_consistent),
        },
        "guards": {
            "observed_data_used": False,
            "observed_odd_used": False,
            "E63_used_as_science_input": False,
            "galaxy_response_fitted": False,
            "survey_covariance_used": False,
            "pvalue_or_detection_sigma_computed": False,
            "cosmological_self_consistent_nbody": False,
            "neutrino_self_gravity": False,
            "delta_f_method": False,
        },
        "next_stage_if_pass": "K01C_REMOTE_COSMOLOGICAL_MATCHED_PHASE_ENGINE_FEASIBILITY",
        "next_stage_if_fail": "Diagnose only the preregistered local setup; do not tune to observed data.",
    }

    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print("K01B_PREREG_PIN_PASS", PREREG_BLOB, flush=True)
    print("K01B_DURATION_GYR", duration, flush=True)
    print("K01B_HALO", mass_msun, c200, rs_phys, flush=True)
    print("K01B_FD_D_KMS", Dfd, flush=True)
    print("K01B_HIDDEN_C_FINE_KMS", C, flush=True)
    print("K01B_HIDDEN_RELATIVE_TO_FD", relative_contrast, flush=True)
    print("K01B_TIMESTEP_GAP", timestep_gap, flush=True)
    print("K01B_NESTED_GAP", nested_gap, flush=True)
    print("K01B_CONTRAST_TO_NUMERICAL_FLOOR", contrast_to_floor, flush=True)
    print("K01B_WIND_REVERSAL_MAX_RESIDUAL", parity_residual, flush=True)
    print("K01B_SCIENCE_GATE_PASS", science_gate_pass, flush=True)
    print("K01B_STATUS", result["status"], flush=True)
    print("K01B_RESULT", OUT, flush=True)


if __name__ == "__main__":
    main()
