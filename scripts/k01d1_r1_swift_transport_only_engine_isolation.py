#!/usr/bin/env python3
"""K01D1-R1 executable SWIFT transport-only engine isolation.

This is a child of the engine-blocked physical-background K01D1 attempt.
It intentionally sets N_nu=0 only to bypass SWIFT's massive-neutrino
background GSL integral while retaining PartType6 particles and SWIFT's
species-based relativistic neutrino drift. No self-gravity is enabled.

Modes:
  prepare  -> create matched FD/Fplus/Fminus PartType6 ICs and YAMLs
  analyze  -> validate completed SWIFT snapshots against preregistered gates
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PREREG = ROOT / "source_data/k01d1_r1_transport_only_engine_isolation_prereg_2026-11-06.json"
PREREG_BLOB = "2daff97168f843ce53bd766810cbad063d15ab85"
WORK = ROOT / "k01d1_r1_swift_work"
OUT = ROOT / "k01d1_r1_transport_only_result.json"
E8_PATH = ROOT / "scripts/build_eboss_dr16_a03_e8_frozen_resonant_occupation.py"

MPC_CGS = 3.0856775814913673e24
MSUN_CGS = 1.98847e33
UNIT_MASS_CGS = 1.0e10 * MSUN_CGS
UNIT_VEL_CGS = 1.0e5
UNIT_TIME_CGS = MPC_CGS / UNIT_VEL_CGS
C_KMS = 299792.458
KB_EV_K = 8.617333262e-5


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01D1_R1_STOP: " + msg)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def load_prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "preregistration changed")
    p = json.loads(raw)
    need(p["stage"] == "K01D1_R1_TRANSPORT_ONLY_ENGINE_ISOLATION_PREREG", "stage")
    need(p["artificial_transport_background"]["N_nu"] == 0, "N_nu isolation drift")
    need(p["engine_modes"]["self_gravity"] is False, "self-gravity forbidden")
    need(
        p["engine_modes"]["physical_external_acceleration"] == [0.0, 0.0, 0.0],
        "external force not zero",
    )
    need(p["particle_design"]["particles_per_state"] == 1024, "particle count drift")
    need(p["particle_design"]["box_comoving_Mpc"] == 512.0, "box drift")
    return p


def load_e8():
    spec = importlib.util.spec_from_file_location("k01d1_r1_e8", E8_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import frozen E8 builder")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    mod.source_gate()
    return mod.kinetic_source()


def make_header(f: h5py.File, p: dict, nnu: int) -> None:
    box = float(p["particle_design"]["box_comoving_Mpc"])
    a0 = float(p["artificial_transport_background"]["a_begin"])
    z0 = 1.0 / a0 - 1.0

    h = f.create_group("Header")
    counts = np.zeros(7, dtype=np.uint64)
    counts[6] = nnu
    h.attrs["BoxSize"] = np.array([box, box, box], dtype=np.float64)
    h.attrs["Dimension"] = np.int32(3)
    h.attrs["NumPart_Total"] = counts
    h.attrs["NumPart_Total_HighWord"] = np.zeros(7, dtype=np.uint64)
    h.attrs["NumPart_ThisFile"] = counts
    h.attrs["NumFilesPerSnapshot"] = np.int32(1)
    h.attrs["MassTable"] = np.zeros(7, dtype=np.float64)
    h.attrs["Flag_Entropy_ICs"] = np.int32(0)
    h.attrs["Time"] = np.float64(a0)
    h.attrs["Redshift"] = np.float64(z0)

    u = f.create_group("Units")
    u.attrs["Unit length in cgs (U_L)"] = np.float64(MPC_CGS)
    u.attrs["Unit mass in cgs (U_M)"] = np.float64(UNIT_MASS_CGS)
    u.attrs["Unit time in cgs (U_t)"] = np.float64(UNIT_TIME_CGS)
    u.attrs["Unit current in cgs (U_I)"] = np.float64(1.0)
    u.attrs["Unit temperature in cgs (U_T)"] = np.float64(1.0)


def write_ic(path: Path, p: dict, pos, vel, mass, ids) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        make_header(f, p, len(ids))
        g = f.create_group("PartType6")
        g.create_dataset("Coordinates", data=np.asarray(pos, dtype=np.float64))
        g.create_dataset("Velocities", data=np.asarray(vel, dtype=np.float32))
        g.create_dataset("Masses", data=np.asarray(mass, dtype=np.float32))
        g.create_dataset("ParticleIDs", data=np.asarray(ids, dtype=np.uint64))


def yaml_text(p: dict, state: str) -> str:
    c = p["artificial_transport_background"]
    a0 = float(c["a_begin"])
    a1 = float(c["a_end"])
    return f"""InternalUnitSystem:
  UnitMass_in_cgs:     {UNIT_MASS_CGS:.16e}
  UnitLength_in_cgs:   {MPC_CGS:.16e}
  UnitVelocity_in_cgs: {UNIT_VEL_CGS:.16e}
  UnitCurrent_in_cgs:  1.0000000000000000e+00
  UnitTemp_in_cgs:     1.0000000000000000e+00

Cosmology:
  Omega_cdm:     {c["Omega_cdm"]:.16g}
  Omega_b:       {c["Omega_b"]:.16g}
  Omega_lambda:  {c["Omega_lambda"]:.16g}
  h:             {c["h"]:.16g}
  a_begin:       {a0:.16g}
  a_end:         {a1:.16g}
  N_nu:          0
  N_ur:          0.0

TimeIntegration:
  dt_min:    1.0000e-12
  dt_max:    1.0000e-03

Scheduler:
  max_top_level_cells: 16
  cell_split_size: 100

Snapshots:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0
  basename: r1_snap
  compression: 0

Statistics:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0

InitialConditions:
  file_name: science_ic.hdf5
  periodic: 1

Neutrino:
  use_delta_f: 0
  use_delta_f_mesh_only: 0
  generate_ics: 0
  use_linear_response: 0
  use_model_none: 0

ConstantPotential:
  g_cgs: [0.0, 0.0, 0.0]

Gravity:
  mesh_side_length: 32
  MAC: geometric
  eta: 0.025
  theta_cr: 0.5
  rebuild_frequency: 0.0
  comoving_DM_softening: 0.04
  max_physical_DM_softening: 0.04
  comoving_nu_softening: 0.04
  max_physical_nu_softening: 0.04

MetaData:
  run_name: K01D1-R1-{state}
"""


def deterministic_geometry(p: dict):
    n = int(p["particle_design"]["particles_per_state"])
    box = float(p["particle_design"]["box_comoving_Mpc"])
    seed = int(p["particle_design"]["seed_positions_directions"])
    need(n % 2 == 0, "antithetic design needs even N")
    half = n // 2
    rng = np.random.default_rng(seed)

    dpos = rng.uniform(-0.20 * box, 0.20 * box, size=(half, 3))
    pos = np.empty((n, 3), dtype=np.float64)
    pos[0::2] = 0.5 * box + dpos
    pos[1::2] = 0.5 * box - dpos

    d = rng.normal(size=(half, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    dirs = np.empty((n, 3), dtype=np.float64)
    dirs[0::2] = d
    dirs[1::2] = -d

    ids = np.arange(1, n + 1, dtype=np.uint64)
    return pos, dirs, ids


def q_samples(qgrid: np.ndarray, f: np.ndarray, n: int) -> np.ndarray:
    half = n // 2
    pdf = qgrid * qgrid * f
    cdf = np.zeros_like(qgrid, dtype=np.float64)
    cdf[1:] = np.cumsum(
        0.5 * (pdf[1:] + pdf[:-1]) * (qgrid[1:] - qgrid[:-1])
    )
    need(np.isfinite(cdf[-1]) and cdf[-1] > 0, "invalid q CDF")
    cdf /= cdf[-1]
    keep = np.concatenate(([True], np.diff(cdf) > 0))
    u = (np.arange(half, dtype=np.float64) + 0.5) / half
    qb = np.interp(u, cdf[keep], qgrid[keep])
    q = np.empty(n, dtype=np.float64)
    q[0::2] = qb
    q[1::2] = qb
    return q


def average_ranks(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    order = np.argsort(x, kind="mergesort")
    s = x[order]
    ranks = np.empty(len(x), dtype=np.float64)
    i = 0
    while i < len(x):
        j = i + 1
        while j < len(x) and s[j] == s[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * ((i + 1) + j)
        i = j
    return ranks


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float:
    rx = average_ranks(x)
    ry = average_ranks(y)
    rx -= rx.mean()
    ry -= ry.mean()
    den = math.sqrt(float(np.dot(rx, rx) * np.dot(ry, ry)))
    need(den > 0, "degenerate Spearman ranks")
    return float(np.dot(rx, ry) / den)


def find_final_snapshot(run: Path):
    files = sorted(run.glob("r1_snap_*.hdf5"))
    need(files, f"no snapshot in {run}")
    rec = []
    for path in files:
        with h5py.File(path, "r") as f:
            if "Cosmology" in f and "Scale-factor" in f["Cosmology"].attrs:
                a = float(np.asarray(f["Cosmology"].attrs["Scale-factor"]).ravel()[0])
            elif "Scale-factor" in f["Header"].attrs:
                a = float(np.asarray(f["Header"].attrs["Scale-factor"]).ravel()[0])
            else:
                a = float(np.asarray(f["Header"].attrs["Time"]).ravel()[0])
        rec.append((a, path))
    return max(rec, key=lambda z: z[0])


def prepare(p: dict) -> None:
    src = load_e8()
    qgrid = np.asarray(src["q"], dtype=np.float64)
    fs = {
        "FD": np.asarray(src["f0"], dtype=np.float64),
        "Fplus": np.asarray(src["fp"], dtype=np.float64),
        "Fminus": np.asarray(src["fm"], dtype=np.float64),
    }

    n = int(p["particle_design"]["particles_per_state"])
    pos, dirs, ids = deterministic_geometry(p)
    a0 = float(p["artificial_transport_background"]["a_begin"])
    kin = p["custom_particle_kinematics"]
    T0_eV = float(kin["physical_T_nu_0_K_used_to_build_IC_velocities"]) * KB_EV_K
    mnu = float(kin["physical_M_nu_eV_used_to_build_IC_velocities"])
    pref = C_KMS * T0_eV / (a0 * mnu)
    mass_value = float(kin["macro_particle_mass_internal_units"])

    summary = {}
    for name in ("FD", "Fplus", "Fminus"):
        q = q_samples(qgrid, fs[name], n)
        vel = (pref * q)[:, None] * dirs
        mass = np.full(n, mass_value, dtype=np.float64)
        run = WORK / name
        run.mkdir(parents=True, exist_ok=True)
        write_ic(run / "science_ic.hdf5", p, pos, vel, mass, ids)
        (run / "science.yml").write_text(yaml_text(p, name))
        np.savez_compressed(
            run / "input_reference.npz",
            ids=ids,
            positions=pos,
            directions=dirs,
            q=q,
            velocities=vel,
            masses=mass,
        )
        summary[name] = {
            "q_mean": float(q.mean()),
            "q_min": float(q.min()),
            "q_max": float(q.max()),
            "mean_input_speed_km_s": float(np.linalg.norm(vel, axis=1).mean()),
            "max_input_speed_km_s": float(np.linalg.norm(vel, axis=1).max()),
        }
        print(
            "K01D1_R1_INPUT",
            name,
            "QMEAN", summary[name]["q_mean"],
            "VMAX_KMS", summary[name]["max_input_speed_km_s"],
            flush=True,
        )
    (WORK / "input_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("K01D1_R1_INPUTS_READY", flush=True)


def analyze(p: dict) -> None:
    a0 = float(p["artificial_transport_background"]["a_begin"])
    a1 = float(p["artificial_transport_background"]["a_end"])
    expected_scale = a0 / a1
    box = float(p["particle_design"]["box_comoving_Mpc"])
    gate = p["preregistered_acceptance"]
    tol = float(gate["expected_velocity_scaling_max_relative_error_at_most"])

    result_states = {}
    q_by = {}
    speed_by = {}
    disp_by = {}
    mass_means = []

    for name in ("FD", "Fplus", "Fminus"):
        run = WORK / name
        ref = np.load(run / "input_reference.npz")
        ids0 = np.asarray(ref["ids"], dtype=np.uint64)
        x0 = np.asarray(ref["positions"], dtype=np.float64)
        v0 = np.asarray(ref["velocities"], dtype=np.float64)
        q0 = np.asarray(ref["q"], dtype=np.float64)
        m0 = np.asarray(ref["masses"], dtype=np.float64)

        a_snap, snap = find_final_snapshot(run)
        with h5py.File(snap, "r") as f:
            g = f["PartType6"]
            ids = np.asarray(g["ParticleIDs"][:], dtype=np.uint64)
            x = np.asarray(g["Coordinates"][:], dtype=np.float64)
            v = np.asarray(g["Velocities"][:], dtype=np.float64)
            m = np.asarray(g["Masses"][:], dtype=np.float64)
            weights = np.asarray(g["Weights"][:], dtype=np.float64) if "Weights" in g else None
            sampled = np.asarray(g["SampledSpeeds"][:], dtype=np.float64) if "SampledSpeeds" in g else None

        o = np.argsort(ids)
        o0 = np.argsort(ids0)
        ids, x, v, m = ids[o], x[o], v[o], m[o]
        ids0, x0, v0, q0, m0 = (
            ids0[o0], x0[o0], v0[o0], q0[o0], m0[o0]
        )
        if weights is not None:
            weights = weights[o]
        if sampled is not None:
            sampled = sampled[o]

        need(np.array_equal(ids, ids0), f"{name} particle IDs changed")
        need(np.all(m > 0), f"{name} nonpositive masses")
        mass_rel = float(np.max(np.abs(m - m0) / np.maximum(np.abs(m0), 1e-300)))
        need(mass_rel <= 2e-6, f"{name} mass changed {mass_rel}")
        mass_means.append(float(m.mean()))

        vexp = v0 * expected_scale
        vrel = np.linalg.norm(v - vexp, axis=1) / np.maximum(
            np.linalg.norm(vexp, axis=1), 1e-30
        )
        max_vrel = float(vrel.max())
        rms_vrel = float(np.sqrt(np.mean(vrel * vrel)))
        need(max_vrel <= tol, f"{name} velocity scaling error {max_vrel}")

        if weights is not None:
            wgap = float(np.max(np.abs(weights - 1.0)))
            need(wgap <= 1e-12, f"{name} non-unity weights {wgap}")
        else:
            wgap = None
        if sampled is not None:
            smax = float(np.max(np.abs(sampled)))
            need(smax == 0.0, f"{name} nonzero SampledSpeeds {smax}")
        else:
            smax = None

        dx = (x - x0 + 0.5 * box) % box - 0.5 * box
        disp = np.linalg.norm(dx, axis=1)
        rho = spearman_rho(q0, disp)
        need(
            rho >= float(gate["per_state_q_vs_periodic_displacement_spearman_at_least"]),
            f"{name} q-displacement Spearman {rho}",
        )

        speed0 = np.linalg.norm(v0, axis=1)
        speed = np.linalg.norm(v, axis=1)
        max_travel = float(disp.max())
        need(max_travel < 0.5 * box, f"{name} periodic displacement ambiguous {max_travel}")

        result_states[name] = {
            "snapshot": str(snap.relative_to(ROOT)),
            "snapshot_scale_factor": a_snap,
            "max_velocity_vector_relative_error": max_vrel,
            "rms_velocity_vector_relative_error": rms_vrel,
            "q_displacement_spearman": rho,
            "mass_max_relative_change": mass_rel,
            "weights_max_abs_from_unity_if_present": wgap,
            "sampled_speeds_max_abs_if_present": smax,
            "mean_input_speed_km_s": float(speed0.mean()),
            "mean_output_speed_km_s": float(speed.mean()),
            "mean_displacement_Mpc": float(disp.mean()),
            "max_displacement_Mpc": max_travel,
        }
        q_by[name] = q0
        speed_by[name] = speed
        disp_by[name] = disp

    mass_cross = (max(mass_means) - min(mass_means)) / max(np.mean(mass_means), 1e-300)
    need(mass_cross <= 1e-12, f"cross-state macro mass mismatch {mass_cross}")

    pairs = (("FD", "Fplus"), ("FD", "Fminus"), ("Fplus", "Fminus"))
    ordering = {}
    ordering_all = True
    displacement_rms = {}
    distinct = True

    for a, b in pairs:
        dq = q_by[a] - q_by[b]
        ds = speed_by[a] - speed_by[b]
        mask = np.abs(dq) > 1e-13
        frac = 1.0 if not np.any(mask) else float(np.mean(np.sign(dq[mask]) == np.sign(ds[mask])))
        ordering[f"{a}_vs_{b}"] = {
            "non_tied_particles": int(mask.sum()),
            "fraction_matching_input_q_order": frac,
        }
        ordering_all &= frac == 1.0

        gap = disp_by[a] - disp_by[b]
        rms = float(np.sqrt(np.mean(gap * gap)))
        displacement_rms[f"{a}_vs_{b}"] = rms
        distinct &= rms > 1e-8

    need(ordering_all, "cross-state final speed ordering changed")
    need(distinct, "state-specific transport collapsed to identical displacements")

    result = {
        "stage": p["stage"],
        "status": "PASS_SWIFT_CUSTOM_F_TRANSPORT_ONLY_ENGINE_ISOLATION",
        "science_gate_pass": True,
        "swift_commit": p["swift"]["commit"],
        "artificial_transport_background": p["artificial_transport_background"],
        "custom_particle_kinematics": p["custom_particle_kinematics"],
        "expected_force_free_snapshot_velocity_scale": expected_scale,
        "states": result_states,
        "cross_state_macro_mass_relative_gap": mass_cross,
        "state_speed_ordering": ordering,
        "state_speed_ordering_all_match": bool(ordering_all),
        "pairwise_displacement_rms_difference_Mpc": displacement_rms,
        "state_specific_displacement_not_identical": bool(distinct),
        "interpretation": p["interpretation_boundary"],
        "guards": p["guards"],
        "next_stage": p["stopping_and_next"]["if_pass"],
    }
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    print("K01D1_R1_EXPECTED_VELOCITY_SCALE", expected_scale, flush=True)
    for name in ("FD", "Fplus", "Fminus"):
        s = result_states[name]
        print(
            "K01D1_R1_STATE", name,
            "VEL_MAX_REL", s["max_velocity_vector_relative_error"],
            "SPEARMAN", s["q_displacement_spearman"],
            "MEAN_DX_MPC", s["mean_displacement_Mpc"],
            "MAX_DX_MPC", s["max_displacement_Mpc"],
            flush=True,
        )
    print("K01D1_R1_STATE_ORDERING_PASS", ordering_all, flush=True)
    print("K01D1_R1_DISTINCT_TRANSPORT_PASS", distinct, displacement_rms, flush=True)
    print("K01D1_R1_SWIFT_CUSTOM_F_TRANSPORT_ONLY_ENGINE_ISOLATION_PASS", flush=True)
    print("K01D1_R1_RESULT", OUT, flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("prepare", "analyze"))
    args = ap.parse_args()
    p = load_prereg()
    WORK.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare":
        prepare(p)
    else:
        analyze(p)


if __name__ == "__main__":
    main()
