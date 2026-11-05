#!/usr/bin/env python3
"""K01A: CI-only feasibility of matched particle sampling for frozen F0/F+/F- states.

This is a numerical implementation check, not an N-body run, survey forecast,
observed-data analysis, or detection claim. It reuses the frozen E8 source
construction and verifies that the three positive isotropic distributions can
be mapped to matched equal-weight particle momentum magnitudes using identical
stratified quantiles.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
E8_PATH = ROOT / "scripts" / "build_eboss_dr16_a03_e8_frozen_resonant_occupation.py"
OUT = ROOT / "k01a_ci_result.json"

spec = importlib.util.spec_from_file_location("e8_frozen", E8_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError("Cannot load frozen E8 module")
e8 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = e8
spec.loader.exec_module(e8)

# Reuse the original source pins and frozen physical construction.
e8.source_gate()
src = e8.kinetic_source()

q = np.asarray(src["q"], dtype=np.float64)
states = {
    "FD": np.asarray(src["f0"], dtype=np.float64),
    "Fplus": np.asarray(src["fp"], dtype=np.float64),
    "Fminus": np.asarray(src["fm"], dtype=np.float64),
}

if q.shape != (4000,) or not np.all(np.diff(q) > 0):
    raise RuntimeError("Frozen q grid changed")
for name, f in states.items():
    if f.shape != q.shape or not np.all(np.isfinite(f)) or not np.all(f > 0):
        raise RuntimeError(f"{name}: invalid frozen distribution")

symmetry_gap = float(
    np.max(
        np.abs(states["Fplus"] + states["Fminus"] - 2.0 * states["FD"])
        / np.maximum(2.0 * np.abs(states["FD"]), 1e-300)
    )
)
if symmetry_gap >= 5e-12:
    raise RuntimeError(f"F+/F- symmetry gap too large: {symmetry_gap}")

mom0 = np.asarray(src["moments_0"], dtype=np.float64)
exact_moment_gaps = {}
for name, key in (("Fplus", "moments_plus"), ("Fminus", "moments_minus")):
    mom = np.asarray(src[key], dtype=np.float64)
    gap = np.abs(mom - mom0) / np.maximum(np.abs(mom0), 1e-300)
    exact_moment_gaps[name] = gap.tolist()
    if float(np.max(gap)) >= 1e-12:
        raise RuntimeError(f"{name}: exact source moments no longer match FD")

# Equal-weight isotropic particle sampling uses p(q) proportional to q^2 F(q).
# Identical stratified quantiles give a matched common-random-number ensemble.
N = 4096
u = (np.arange(N, dtype=np.float64) + 0.5) / N

def cumulative_trapezoid(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    c = np.zeros_like(x)
    c[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]))
    return c

def quantile_sample(f: np.ndarray) -> np.ndarray:
    pdf = q * q * f
    cdf = cumulative_trapezoid(pdf, q)
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise RuntimeError("Invalid momentum CDF")
    cdf /= cdf[-1]
    if np.any(np.diff(cdf) < 0):
        raise RuntimeError("Non-monotone momentum CDF")
    keep = np.concatenate(([True], np.diff(cdf) > 0))
    s = np.interp(u, cdf[keep], q[keep])
    if s.shape != (N,) or not np.all(np.isfinite(s)) or np.any(s < 0):
        raise RuntimeError("Invalid quantile sample")
    return s

samples = {name: quantile_sample(f) for name, f in states.items()}

# At z_match, equal-weight particle averages of epsilon and q^2/(3 epsilon)
# estimate rho/n and P/n. These are implementation QA only.
tnu0_ev = e8.cro.T_NCDM * e8.cro.TCMB_K * e8.cro.KB_EV_K
a_match = 1.0 / (1.0 + e8.Z_MATCH)
y = a_match * e8.MASS_EV / tnu0_ev

sample_stats = {}
for name, s in samples.items():
    eps = np.sqrt(s * s + y * y)
    sample_stats[name] = {
        "q_mean": float(np.mean(s)),
        "epsilon_mean_zmatch": float(np.mean(eps)),
        "pressure_proxy_mean_zmatch": float(np.mean(s * s / (3.0 * eps))),
        "q_min": float(np.min(s)),
        "q_max": float(np.max(s)),
    }

sample_gaps = {}
for name in ("Fplus", "Fminus"):
    sample_gaps[name] = {}
    for key in ("epsilon_mean_zmatch", "pressure_proxy_mean_zmatch"):
        a = sample_stats[name][key]
        b = sample_stats["FD"][key]
        sample_gaps[name][key] = float(abs(a - b) / max(abs(a), abs(b), 1e-300))

max_sample_gap = max(v for d in sample_gaps.values() for v in d.values())
if max_sample_gap >= 1e-2:
    raise RuntimeError(f"Matched quantile particle QA unexpectedly poor: {max_sample_gap}")

result = {
    "stage": "K01A_MATCHED_CUSTOM_F_PARTICLE_SAMPLING_CI",
    "status": "PASS_MATCHED_F0_FPLUS_FMINUS_PARTICLE_SAMPLING_FEASIBILITY",
    "interpretation": (
        "Implementation feasibility only: frozen positive isotropic F0/F+/F- states "
        "can be sampled as matched equal-weight particle momentum magnitudes using "
        "identical stratified quantiles. This is not an N-body, halo-response, "
        "survey, or detectability result."
    ),
    "guards": {
        "observed_data_used": False,
        "observed_odd_used": False,
        "nbody_run_performed": False,
        "survey_signal_to_noise_computed": False,
        "galaxy_response_fitted": False,
        "e63_used_as_science_input": False,
    },
    "frozen_source": {
        "q_nodes": int(len(q)),
        "mass_eV": float(e8.MASS_EV),
        "source_matching_z": float(e8.Z_MATCH),
        "maximum_fractional_distortion": float(e8.FRAC),
        "max_Fplus_Fminus_symmetry_relative_gap": symmetry_gap,
        "exact_moment_relative_gaps_vs_FD": exact_moment_gaps,
    },
    "matched_particle_design": {
        "N_test_particles_per_state": N,
        "same_stratified_quantiles": True,
        "same_particle_ids_compatible": True,
        "same_positions_compatible": True,
        "same_direction_realization_compatible": True,
        "momentum_magnitude_pdf": "p(q) proportional to q^2 F_state(q)",
    },
    "sample_stats": sample_stats,
    "sample_moment_proxy_relative_gaps_vs_FD": sample_gaps,
    "max_sample_proxy_gap": float(max_sample_gap),
    "next_stage": "K01B_REMOTE_SMALL_MATCHED_PHASE_DYNAMICS_SMOKE_TEST",
}

OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")

print("K01A_SOURCE_PINS_PASS", flush=True)
print("K01A_EXACT_MOMENT_MATCH_PASS", exact_moment_gaps, flush=True)
print("K01A_FPM_SYMMETRY_PASS", symmetry_gap, flush=True)
print("K01A_MATCHED_QUANTILE_SAMPLE_MAX_PROXY_GAP", max_sample_gap, flush=True)
print("K01A_MATCHED_CUSTOM_F_PARTICLE_SAMPLING_FEASIBILITY_PASS", flush=True)
print("K01A_RESULT", OUT, flush=True)
