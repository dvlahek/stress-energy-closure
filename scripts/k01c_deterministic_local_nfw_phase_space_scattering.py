#!/usr/bin/env python3
"""K01C deterministic local NFW phase-space scattering quadrature.

This replaces the failed random-particle mean-drag estimator of K01B with a
fully deterministic scattering kernel. It keeps the same frozen E8 states and
the same lower-mass E28 NFW halo anchor, but it is a new local estimand:
an incoming-flux-weighted momentum-transfer cross-section proxy.

It is not a cosmological N-body calculation, total halo drag, galaxy response,
survey observable, or detection forecast.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "source_data/k01c_deterministic_local_nfw_phase_space_scattering_prereg_2026-11-05.json"
PREREG_BLOB = "c7a7dbdb1609ec3b94827018e6ec8b0e47b78a9d"
OUT = ROOT / "k01c_ci_result.json"

sys.path.insert(0, str(ROOT / "scripts"))
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as e27
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as e28

E8_PATH = ROOT / "scripts/build_eboss_dr16_a03_e8_frozen_resonant_occupation.py"
spec = importlib.util.spec_from_file_location("k01c_e8_frozen", E8_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("Cannot load frozen E8 source builder")
e8 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = e8
spec.loader.exec_module(e8)


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01C_STOP: " + msg)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "preregistration changed")
    p = json.loads(raw)
    n = p["deterministic_numerics"]
    need(
        p["stage"] == "K01C_DETERMINISTIC_LOCAL_NFW_PHASE_SPACE_SCATTERING_PREREG"
        and p["frozen_physics"]["evaluation_redshift"] == 0.95
        and p["frozen_physics"]["halo_relative_wind_km_s"] == 200
        and n["fine_q_legendre_nodes_per_interval"] == 12
        and n["coarse_q_legendre_nodes_per_interval"] == 8
        and n["fine_mu_legendre_nodes"] == 32
        and n["coarse_mu_legendre_nodes"] == 20
        and n["fine_b_legendre_nodes_0_to_r200"] == 48
        and n["coarse_b_legendre_nodes_0_to_r200"] == 32
        and n["fine_orbit_angle_legendre_nodes"] == 64
        and n["coarse_orbit_angle_legendre_nodes"] == 40
        and n["turning_point_solver"].startswith("80-step"),
        "preregistered design drift",
    )
    return p


def q_quadrature(intervals, nodes_per_interval: int):
    x, w = leggauss(nodes_per_interval)
    qs = []
    ws = []
    for a, b in intervals:
        qs.append(0.5 * (b - a) * x + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(qs), np.concatenate(ws)


def nfw_phi(r: np.ndarray, M: float, c: float, rs: float) -> np.ndarray:
    """Truncated NFW potential, zero at infinity, in (km/s)^2."""
    r = np.asarray(r, dtype=np.float64)
    r200 = c * rs
    fc = math.log1p(c) - c / (1.0 + c)

    out = np.empty_like(r)
    inside = r < r200

    if np.any(inside):
        x = r[inside] / rs
        ratio = np.log1p(x) / x
        out[inside] = -(e27.G * M / (rs * fc)) * (
            ratio - 1.0 / (1.0 + c)
        )
    if np.any(~inside):
        out[~inside] = -e27.G * M / r[~inside]
    return out


def scattering_J(
    speeds: np.ndarray,
    M: float,
    c: float,
    rs: float,
    nb: int,
    nt: int,
    chunk: int = 64,
):
    """Return J(v;bmax) for bmax=r200,2r200,4r200.

    J = integral b [cos(chi)-1] db. The inner 0< b <r200 scattering angle
    is computed from the exact central-orbit integral in the truncated NFW
    potential. Outside r200, the added impact-parameter annuli use the exact
    Kepler point-mass formula.
    """
    speeds = np.asarray(speeds, dtype=np.float64)
    need(np.all(np.isfinite(speeds)) and np.all(speeds > 0), "nonpositive speed")

    r200 = c * rs

    xb, wb = leggauss(nb)
    bnodes = 0.5 * r200 * (xb + 1.0)
    bweights = 0.5 * r200 * wb

    xt, wt = leggauss(nt)
    tnodes = 0.5 * (xt + 1.0)
    tweights = 0.5 * wt
    u_nodes = 1.0 - tnodes * tnodes

    J1 = np.empty_like(speeds)

    for start in range(0, len(speeds), chunk):
        stop = min(start + chunk, len(speeds))
        v = speeds[start:stop, None]
        b = bnodes[None, :]
        h = v * b

        lo = np.maximum(b * 1e-14, r200 * 1e-16)
        hi = np.broadcast_to(b, h.shape).copy()

        # Unique outer turning point lies in (0,b] for an attractive potential.
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            R = v * v - 2.0 * nfw_phi(mid, M, c, rs) - (h * h) / (mid * mid)
            choose_hi = R > 0.0
            hi = np.where(choose_hi, mid, hi)
            lo = np.where(choose_hi, lo, mid)

        rmin = 0.5 * (lo + hi)

        # Smooth the turning-point square-root singularity with
        # u=rmin/r=1-t^2.
        r = rmin[:, :, None] / u_nodes[None, None, :]
        phi = nfw_phi(r, M, c, rs)
        rad = (
            v[:, :, None] * v[:, :, None]
            - 2.0 * phi
            - (h[:, :, None] * h[:, :, None]) / (r * r)
        )
        # Only roundoff-level negatives are admissible.
        scale = np.maximum(v[:, :, None] * v[:, :, None], 1.0)
        min_scaled = float(np.min(rad / scale))
        need(min_scaled > -2e-10, f"orbit radicand negative: {min_scaled}")
        rad = np.maximum(rad, np.finfo(np.float64).tiny)

        integrand = (
            2.0
            * tnodes[None, None, :]
            * (h / rmin)[:, :, None]
            / np.sqrt(rad)
        )
        phi0 = np.sum(integrand * tweights[None, None, :], axis=2)
        chi = 2.0 * phi0 - math.pi
        coschi = np.cos(chi)

        J1[start:stop] = np.sum(
            bweights[None, :] * b * (coschi - 1.0), axis=1
        )

    # Exact Kepler additions outside r200. A=GM/v^2 has length units.
    A = e27.G * M / (speeds * speeds)

    def outer(b1, b2):
        return -(A * A) * np.log(
            (b2 * b2 + A * A) / (b1 * b1 + A * A)
        )

    J2 = J1 + outer(r200, 2.0 * r200)
    J4 = J2 + outer(2.0 * r200, 4.0 * r200)

    need(
        np.all(np.isfinite(J1))
        and np.all(np.isfinite(J2))
        and np.all(np.isfinite(J4))
        and np.all(J1 <= 0)
        and np.all(J2 <= J1)
        and np.all(J4 <= J2),
        "invalid scattering impact-parameter integral",
    )

    return np.stack([J1, J2, J4], axis=1)


def evaluate_grid(
    src: dict,
    M: float,
    c: float,
    rs: float,
    wind: float,
    intervals,
    nq: int,
    nmu: int,
    nb: int,
    nt: int,
):
    q, wq = q_quadrature(intervals, nq)
    mu, wmu_raw = leggauss(nmu)
    wmu = 0.5 * wmu_raw

    p_ev = q * e27.TNU0 * (1.0 + e27.ZOBS)
    vth = e27.C * p_ev / np.sqrt(p_ev * p_ev + e27.MNU * e27.MNU)

    # q x mu halo-frame relative velocity.
    wx = vth[:, None] * mu[None, :] - wind
    vrel = np.sqrt(
        vth[:, None] * vth[:, None]
        + wind * wind
        - 2.0 * vth[:, None] * wind * mu[None, :]
    )
    need(np.all(vrel > 0), "zero relative speed on quadrature node")
    mu_rel = wx / vrel

    flat_v = vrel.ravel()
    J = scattering_J(flat_v, M, c, rs, nb, nt).reshape(
        len(q), len(mu), 3
    )

    kernel = (
        2.0
        * math.pi
        * (vrel[:, :, None] ** 2)
        * mu_rel[:, :, None]
        * J
    )

    fq = {
        "FD": np.interp(q, src["q"], src["f0"]),
        "Fplus": np.interp(q, src["q"], src["fp"]),
        "Fminus": np.interp(q, src["q"], src["fm"]),
    }

    result = {}
    for name, f in fq.items():
        qw = wq * q * q * f
        denom = float(np.sum(qw))
        need(denom > 0, f"{name} q normalization")

        mu_integrated = np.sum(
            kernel * wmu[None, :, None], axis=1
        )
        D = np.sum(
            qw[:, None] * mu_integrated, axis=0
        ) / denom

        result[name] = {
            "D_bmax_1r200": float(D[0]),
            "D_bmax_2r200": float(D[1]),
            "D_bmax_4r200": float(D[2]),
            "q_number_norm": denom,
        }

    result["contrast"] = {}
    for key in ("D_bmax_1r200", "D_bmax_2r200", "D_bmax_4r200"):
        result["contrast"][key] = (
            result["Fminus"][key] - result["Fplus"][key]
        )

    return result


def relgap(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def main() -> None:
    p = load_prereg()
    e8.source_gate()
    src = e8.kinetic_source()

    _, base, profiles = e28.sources()
    _, _, _, halos = base
    profile = profiles[0]

    M = float(halos[0][0])
    c = float(profile["c200c_population_median_DM14"])
    rs_comoving = float(profile["r_s_comoving_Mpc_DM14"])
    rs = rs_comoving / (1.0 + e27.ZOBS)
    r200 = c * rs

    need(
        math.isclose(M, float(profile["M200c_physical_Msun_DM14"]), rel_tol=1e-13),
        "frozen halo mass/profile mismatch",
    )

    n = p["deterministic_numerics"]
    intervals = n["q_intervals"]

    fine = evaluate_grid(
        src, M, c, rs, e27.VHALO, intervals,
        n["fine_q_legendre_nodes_per_interval"],
        n["fine_mu_legendre_nodes"],
        n["fine_b_legendre_nodes_0_to_r200"],
        n["fine_orbit_angle_legendre_nodes"],
    )
    coarse = evaluate_grid(
        src, M, c, rs, e27.VHALO, intervals,
        n["coarse_q_legendre_nodes_per_interval"],
        n["coarse_mu_legendre_nodes"],
        n["coarse_b_legendre_nodes_0_to_r200"],
        n["coarse_orbit_angle_legendre_nodes"],
    )

    # Symmetry controls use the fine deterministic grid.
    zero = evaluate_grid(
        src, M, c, rs, 0.0, intervals,
        n["fine_q_legendre_nodes_per_interval"],
        n["fine_mu_legendre_nodes"],
        n["fine_b_legendre_nodes_0_to_r200"],
        n["fine_orbit_angle_legendre_nodes"],
    )
    reverse = evaluate_grid(
        src, M, c, rs, -e27.VHALO, intervals,
        n["fine_q_legendre_nodes_per_interval"],
        n["fine_mu_legendre_nodes"],
        n["fine_b_legendre_nodes_0_to_r200"],
        n["fine_orbit_angle_legendre_nodes"],
    )

    keys = ("D_bmax_1r200", "D_bmax_2r200", "D_bmax_4r200")
    primary = keys[0]

    state_gaps = {
        s: relgap(fine[s][primary], coarse[s][primary])
        for s in ("FD", "Fplus", "Fminus")
    }
    Cfine = fine["contrast"][primary]
    Ccoarse = coarse["contrast"][primary]
    contrast_gap = relgap(Cfine, Ccoarse)

    contrast_signs = [
        int(np.sign(fine["contrast"][k])) for k in keys
    ]
    same_cutoff_sign = (
        all(s != 0 for s in contrast_signs)
        and len(set(contrast_signs)) == 1
    )
    same_coarse_fine_sign = (
        Cfine != 0 and Ccoarse != 0 and np.sign(Cfine) == np.sign(Ccoarse)
    )

    fd_scale = max(abs(fine["FD"][primary]), 1e-300)
    zero_residual = max(
        abs(zero[s][primary]) for s in ("FD", "Fplus", "Fminus")
    ) / fd_scale

    parity_parts = []
    for s in ("FD", "Fplus", "Fminus"):
        parity_parts.append(
            abs(fine[s][primary] + reverse[s][primary])
            / max(abs(fine[s][primary]), abs(reverse[s][primary]), 1e-300)
        )
    parity_parts.append(
        abs(Cfine + reverse["contrast"][primary])
        / max(abs(Cfine), abs(reverse["contrast"][primary]), 1e-300)
    )
    parity_residual = max(parity_parts)

    gate = p["preregistered_gate"]
    gate_pass = bool(
        same_coarse_fine_sign
        and contrast_gap <= float(gate["primary_contrast_coarse_fine_relative_gap_at_most"])
        and max(state_gaps.values()) <= float(gate["each_state_primary_D_coarse_fine_relative_gap_at_most"])
        and same_cutoff_sign
        and zero_residual <= float(gate["zero_wind_absolute_residual_relative_to_primary_FD_at_most"])
        and parity_residual <= float(gate["wind_reversal_relative_parity_residual_at_most"])
    )

    result = {
        "stage": p["stage"],
        "status": (
            "PASS_NUMERICALLY_RESOLVED_DETERMINISTIC_LOCAL_SCATTERING_CONTRAST"
            if gate_pass
            else "FAIL_DETERMINISTIC_LOCAL_SCATTERING_GATE"
        ),
        "science_gate_pass": gate_pass,
        "interpretation": p["interpretation_boundary"],
        "preregistration_git_blob": PREREG_BLOB,
        "frozen_halo": {
            "M200c_physical_Msun": M,
            "c200c": c,
            "r_s_physical_Mpc": rs,
            "r200_physical_Mpc": r200,
            "wind_km_s": e27.VHALO,
        },
        "units_D": "(km/s)^2 Mpc^2",
        "fine": fine,
        "coarse": coarse,
        "zero_wind": zero,
        "reversed_wind": reverse,
        "numerical_diagnostics": {
            "primary_state_coarse_fine_relative_gaps": state_gaps,
            "primary_contrast_coarse_fine_relative_gap": contrast_gap,
            "primary_contrast_same_sign_coarse_fine": bool(same_coarse_fine_sign),
            "fine_contrast_signs_bmax_1_2_4r200": contrast_signs,
            "contrast_sign_same_bmax_1_2_4r200": bool(same_cutoff_sign),
            "zero_wind_residual_relative_to_primary_FD": zero_residual,
            "wind_reversal_component_relative_residuals": parity_parts,
            "wind_reversal_max_relative_residual": parity_residual,
        },
        "hidden_contrast": {
            "primary_C_bmax_1r200": Cfine,
            "primary_relative_to_abs_FD": Cfine / fd_scale,
            "C_bmax_2r200": fine["contrast"][keys[1]],
            "C_bmax_4r200": fine["contrast"][keys[2]],
        },
        "guards": p["guards"],
        "next_stage_if_pass": "K01D_COSMOLOGICAL_ENGINE_DESIGN_NO_OBSERVATIONS",
        "next_stage_if_fail": p["stopping_rule"],
    }

    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    print("K01C_PREREG_PIN_PASS", PREREG_BLOB, flush=True)
    print("K01C_HALO", M, c, rs, r200, flush=True)
    print("K01C_FD_PRIMARY_D", fine["FD"][primary], flush=True)
    print("K01C_FPLUS_PRIMARY_D", fine["Fplus"][primary], flush=True)
    print("K01C_FMINUS_PRIMARY_D", fine["Fminus"][primary], flush=True)
    print("K01C_PRIMARY_HIDDEN_C", Cfine, flush=True)
    print("K01C_PRIMARY_RELATIVE_TO_FD", Cfine / fd_scale, flush=True)
    print("K01C_CONTRAST_COARSE_FINE_GAP", contrast_gap, flush=True)
    print("K01C_STATE_GAPS", state_gaps, flush=True)
    print("K01C_CUTOFF_SIGNS", contrast_signs, flush=True)
    print("K01C_ZERO_WIND_REL_RESIDUAL", zero_residual, flush=True)
    print("K01C_WIND_REVERSAL_MAX_RESIDUAL", parity_residual, flush=True)
    print("K01C_SCIENCE_GATE_PASS", gate_pass, flush=True)
    print("K01C_STATUS", result["status"], flush=True)
    print("K01C_RESULT", OUT, flush=True)


if __name__ == "__main__":
    main()
