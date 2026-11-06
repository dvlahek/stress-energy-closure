#!/usr/bin/env python3
"""K01D1 executable SWIFT custom-F cosmological transport smoke test helper.

Subcommands:
  prepare-calibration  -> write placeholder PartType6 IC + calibration YAML
  extract-calibration  -> read SWIFT calibration snapshot and freeze macro mass
  prepare-science      -> write matched FD/Fplus/Fminus ICs + YAMLs
  analyze              -> validate the three completed science runs

The actual SWIFT executable is built and invoked by GitHub Actions. This helper
never enables self-gravity and never uses observed data.
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
PREREG = ROOT / "source_data/k01d1_swift_custom_f_transport_smoke_prereg_2026-11-06.json"
PREREG_BLOB = "4519708bde9cfab729651375e22e5176ee353609"
WORK = ROOT / "k01d1_swift_work"
CAL = WORK / "calibration"
SCI = WORK / "science"
CAL_MASS = WORK / "calibrated_macro_mass.json"
OUT = ROOT / "k01d1_swift_transport_result.json"

E8_PATH = ROOT / "scripts/build_eboss_dr16_a03_e8_frozen_resonant_occupation.py"

MPC_CGS = 3.0856775814913673e24
MSUN_CGS = 1.98847e33
UNIT_MASS_CGS = 1.0e10 * MSUN_CGS
UNIT_VEL_CGS = 1.0e5  # 1 km/s
UNIT_TIME_CGS = MPC_CGS / UNIT_VEL_CGS
C_KMS = 299792.458
KB_EV_K = 8.617333262e-5


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def need(ok: bool, msg: str) -> None:
    if not ok:
        raise RuntimeError("K01D1_STOP: " + msg)


def prereg() -> dict:
    raw = PREREG.read_bytes()
    need(git_blob(raw) == PREREG_BLOB, "preregistration changed after execution lock")
    p = json.loads(raw)
    need(
        p["stage"] == "K01D1_EXECUTABLE_SWIFT_CUSTOM_F_COSMOLOGICAL_TRANSPORT_SMOKE_PREREG",
        "stage",
    )
    need(
        p["pre_execution_geometry_amendment"]["status"]
        == "AMENDED_BEFORE_ANY_K01D1_EXECUTION",
        "missing pre-execution geometry amendment",
    )
    need(p["particle_design"]["box_comoving_Mpc"] == 512, "box amendment drift")
    return p


def load_e8():
    spec = importlib.util.spec_from_file_location("k01d1_e8", E8_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import frozen E8 builder")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    mod.source_gate()
    return mod, mod.kinetic_source()


def make_header(f: h5py.File, p: dict, nnu: int) -> None:
    box = float(p["particle_design"]["box_comoving_Mpc"])
    a0 = float(p["cosmology"]["a_begin"])
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


def yaml_text(p: dict, file_name: str, basename: str, generate_ics: int) -> str:
    c = p["cosmology"]
    box = p["particle_design"]["box_comoving_Mpc"]
    a0 = c["a_begin"]
    a1 = c["a_end"]
    # Snapshot exactly at the requested final scale factor. delta_time is
    # multiplicative in cosmological mode and chosen > run extent afterwards.
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
  N_nu:          {int(c["N_nu"])}
  N_ur:          {c["N_ur"]:.16g}
  T_nu_0:        {c["T_nu_0_K"]:.16g}
  M_nu_eV:       {c["M_nu_eV"]:.16g}
  deg_nu:        {c["deg_nu"]:.16g}

TimeIntegration:
  dt_min:    1.0000e-12
  dt_max:    1.0000e-03

Scheduler:
  max_top_level_cells: 16
  cell_split_size: 100

Snapshots:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0
  basename: {basename}
  compression: 0

Statistics:
  scale_factor_first: {a1:.16g}
  delta_time: 2.0

InitialConditions:
  file_name: {file_name}
  periodic: 1

Neutrino:
  use_delta_f: 0
  use_delta_f_mesh_only: 0
  generate_ics: {generate_ics}
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
"""


def deterministic_geometry(p: dict):
    n = int(p["particle_design"]["particles_per_state"])
    need(n % 2 == 0, "antithetic particle count must be even")
    box = float(p["particle_design"]["box_comoving_Mpc"])
    seed = int(p["particle_design"]["seed_positions_directions"])
    rng = np.random.default_rng(seed)
    half = n // 2

    # Center positions away from exact periodic boundaries. Antithetic geometry
    # around the box centre is shared by all states.
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


def q_samples(qgrid, f, n: int, rng: np.random.Generator):
    half = n // 2
    pdf = qgrid * qgrid * f
    cdf = np.zeros_like(qgrid, dtype=np.float64)
    cdf[1:] = np.cumsum(
        0.5 * (pdf[1:] + pdf[:-1]) * (qgrid[1:] - qgrid[:-1])
    )
    need(cdf[-1] > 0 and np.isfinite(cdf[-1]), "invalid q CDF")
    cdf /= cdf[-1]
    keep = np.concatenate(([True], np.diff(cdf) > 0))
    u = (np.arange(half, dtype=np.float64) + 0.5) / half
    # Same deterministic permutation for all states is generated outside this
    # function by resetting the same state-independent ordering.
    qb = np.interp(u, cdf[keep], qgrid[keep])
    q = np.empty(n, dtype=np.float64)
    q[0::2] = qb
    q[1::2] = qb
    return q


def final_snapshot(run_dir: Path, basename: str):
    files = sorted(run_dir.glob(basename + "_*.hdf5"))
    need(files, f"no SWIFT snapshots for {run_dir}")
    records = []
    for path in files:
        with h5py.File(path, "r") as f:
            if "Cosmology" in f and "Scale-factor" in f["Cosmology"].attrs:
                a = float(np.asarray(f["Cosmology"].attrs["Scale-factor"]).ravel()[0])
            elif "Scale-factor" in f["Header"].attrs:
                a = float(np.asarray(f["Header"].attrs["Scale-factor"]).ravel()[0])
            else:
                a = float(np.asarray(f["Header"].attrs["Time"]).ravel()[0])
        records.append((a, path))
    return max(records, key=lambda z: z[0])


def prepare_calibration(p: dict):
    CAL.mkdir(parents=True, exist_ok=True)
    n = int(p["particle_design"]["particles_per_state"])
    pos, _, ids = deterministic_geometry(p)
    vel = np.zeros((n, 3), dtype=np.float32)
    mass = np.zeros(n, dtype=np.float32)
    write_ic(CAL / "calibration_ic.hdf5", p, pos, vel, mass, ids)
    (CAL / "calibration.yml").write_text(
        yaml_text(p, "calibration_ic.hdf5", "cal_snap", 1)
    )
    print("K01D1_CALIBRATION_INPUT_READY", CAL, flush=True)


def extract_calibration(p: dict):
    a, snap = final_snapshot(CAL, "cal_snap")
    with h5py.File(snap, "r") as f:
        m = np.asarray(f["PartType6/Masses"][:], dtype=np.float64)
        ids = np.asarray(f["PartType6/ParticleIDs"][:], dtype=np.uint64)
    need(len(m) == p["particle_design"]["particles_per_state"], "calibration N")
    need(np.all(np.isfinite(m)) and np.all(m > 0), "calibration mass not positive")
    uniform_gap = float(np.ptp(m) / np.mean(m))
    need(uniform_gap <= 2e-7, f"calibration mass not uniform: {uniform_gap}")
    need(len(np.unique(ids)) == len(ids), "calibration IDs not unique")
    result = {
        "snapshot": str(snap.relative_to(ROOT)),
        "snapshot_scale_factor": a,
        "macro_mass_internal_units": float(np.mean(m)),
        "macro_mass_uniform_relative_range": uniform_gap,
    }
    CAL_MASS.write_text(json.dumps(result, indent=2) + "\n")
    print("K01D1_CALIBRATED_MACRO_MASS", result["macro_mass_internal_units"], flush=True)
    print("K01D1_CALIBRATION_MASS_PASS", flush=True)


def prepare_science(p: dict):
    need(CAL_MASS.is_file(), "calibrated macro mass missing")
    cal = json.loads(CAL_MASS.read_text())
    mass0 = float(cal["macro_mass_internal_units"])
    need(mass0 > 0, "bad calibrated macro mass")

    e8, src = load_e8()
    qgrid = np.asarray(src["q"], dtype=np.float64)
    fs = {
        "FD": np.asarray(src["f0"], dtype=np.float64),
        "Fplus": np.asarray(src["fp"], dtype=np.float64),
        "Fminus": np.asarray(src["fm"], dtype=np.float64),
    }

    n = int(p["particle_design"]["particles_per_state"])
    pos, dirs, ids = deterministic_geometry(p)
    a0 = float(p["cosmology"]["a_begin"])
    T0 = float(p["cosmology"]["T_nu_0_K"]) * KB_EV_K
    mnu = float(p["cosmology"]["M_nu_eV"])
    pref = C_KMS * T0 / (a0 * mnu)

    SCI.mkdir(parents=True, exist_ok=True)
    source_summary = {}

    for name in ("FD", "Fplus", "Fminus"):
        # The same stratified u ordering is used by construction.
        q = q_samples(qgrid, fs[name], n, np.random.default_rng(0))
        vel = (pref * q)[:, None] * dirs
        masses = np.full(n, mass0, dtype=np.float64)

        run = SCI / name
        run.mkdir(parents=True, exist_ok=True)
        write_ic(run / "science_ic.hdf5", p, pos, vel, masses, ids)
        (run / "science.yml").write_text(
            yaml_text(p, "science_ic.hdf5", "sci_snap", 0)
        )
        np.savez_compressed(
            run / "input_reference.npz",
            ids=ids,
            positions=pos,
            directions=dirs,
            q=q,
            velocities=vel,
            masses=masses,
        )
        source_summary[name] = {
            "q_min": float(q.min()),
            "q_max": float(q.max()),
            "q_mean": float(q.mean()),
            "speed_input_mean_km_s": float(np.linalg.norm(vel, axis=1).mean()),
            "speed_input_max_km_s": float(np.linalg.norm(vel, axis=1).max()),
        }
        print(
            "K01D1_SCIENCE_INPUT",
            name,
            "QMEAN", source_summary[name]["q_mean"],
            "VMAX", source_summary[name]["speed_input_max_km_s"],
            flush=True,
        )

    (WORK / "science_input_summary.json").write_text(
        json.dumps(source_summary, indent=2) + "\n"
    )
    print("K01D1_SCIENCE_INPUTS_READY", flush=True)


def average_ranks(values):
    x=np.asarray(values,dtype=np.float64)
    order=np.argsort(x,kind="mergesort")
    ranks=np.empty(len(x),dtype=np.float64)
    s=x[order]
    i=0
    while i<len(x):
        j=i+1
        while j<len(x) and s[j]==s[i]:
            j+=1
        # Average of 1-based ranks i+1 ... j.
        ranks[order[i:j]]=0.5*((i+1)+j)
        i=j
    return ranks


def spearman_rho(x,y):
    rx=average_ranks(x)
    ry=average_ranks(y)
    rx-=rx.mean(); ry-=ry.mean()
    den=math.sqrt(float(np.dot(rx,rx)*np.dot(ry,ry)))
    need(den>0,"degenerate Spearman ranks")
    return float(np.dot(rx,ry)/den)


def analyze(p: dict):
    a0 = float(p["cosmology"]["a_begin"])
    a1 = float(p["cosmology"]["a_end"])
    expected_scale = a0 / a1
    box = float(p["particle_design"]["box_comoving_Mpc"])
    tol = float(p["acceptance"]["expected_velocity_scaling_relative_error_at_most"])

    state = {}
    final_speed_by_state = {}
    input_q_by_state = {}
    displacement_by_state = {}

    for name in ("FD", "Fplus", "Fminus"):
        run = SCI / name
        ref = np.load(run / "input_reference.npz")
        ids0 = np.asarray(ref["ids"], dtype=np.uint64)
        x0 = np.asarray(ref["positions"], dtype=np.float64)
        v0 = np.asarray(ref["velocities"], dtype=np.float64)
        q0 = np.asarray(ref["q"], dtype=np.float64)
        m0 = np.asarray(ref["masses"], dtype=np.float64)

        a_snap, snap = final_snapshot(run, "sci_snap")
        with h5py.File(snap, "r") as f:
            ids = np.asarray(f["PartType6/ParticleIDs"][:], dtype=np.uint64)
            x = np.asarray(f["PartType6/Coordinates"][:], dtype=np.float64)
            v = np.asarray(f["PartType6/Velocities"][:], dtype=np.float64)
            m = np.asarray(f["PartType6/Masses"][:], dtype=np.float64)
            weights = np.asarray(f["PartType6/Weights"][:], dtype=np.float64) if "Weights" in f["PartType6"] else None
            sampled = np.asarray(f["PartType6/SampledSpeeds"][:], dtype=np.float64) if "SampledSpeeds" in f["PartType6"] else None

        order = np.argsort(ids)
        ids = ids[order]; x = x[order]; v = v[order]; m = m[order]
        order0 = np.argsort(ids0)
        ids0 = ids0[order0]; x0 = x0[order0]; v0 = v0[order0]
        q0 = q0[order0]; m0 = m0[order0]

        need(np.array_equal(ids, ids0), f"{name} IDs changed")
        mass_rel = float(np.max(np.abs(m - m0) / np.maximum(np.abs(m0), 1e-300)))
        need(mass_rel <= 2e-6, f"{name} mass changed {mass_rel}")

        v_expected = v0 * expected_scale
        vel_rel = np.linalg.norm(v - v_expected, axis=1) / np.maximum(
            np.linalg.norm(v_expected, axis=1), 1e-30
        )
        max_vel_rel = float(np.max(vel_rel))
        rms_vel_rel = float(np.sqrt(np.mean(vel_rel * vel_rel)))
        need(max_vel_rel <= tol, f"{name} force-free velocity scaling {max_vel_rel}")

        # Periodic minimum-image displacement is unambiguous because the
        # amended 512 Mpc box is much wider than twice the travel distance.
        dx = (x - x0 + 0.5 * box) % box - 0.5 * box
        disp = np.linalg.norm(dx, axis=1)
        rho = spearman_rho(q0, disp)
        need(np.isfinite(rho), f"{name} displacement Spearman nonfinite")
        need(
            rho >= float(p["acceptance"]["per_state_q_vs_displacement_spearman_at_least"]),
            f"{name} q-displacement Spearman {rho}",
        )

        speed0 = np.linalg.norm(v0, axis=1)
        speed = np.linalg.norm(v, axis=1)
        speed_scale = speed / np.maximum(speed0, 1e-30)
        scale_err = float(np.max(np.abs(speed_scale - expected_scale) / expected_scale))

        if weights is not None:
            need(
                float(np.max(np.abs(weights - 1.0))) <= 1e-12,
                f"{name} non-unity neutrino weights imply hidden delta-f behavior",
            )
        if sampled is not None:
            # Plain particle mode intentionally reports zero because it does
            # not reconstruct an FD initial momentum.
            need(
                float(np.max(np.abs(sampled))) == 0.0,
                f"{name} SampledSpeeds nonzero implies FD reconstruction",
            )

        state[name] = {
            "snapshot": str(snap.relative_to(ROOT)),
            "snapshot_scale_factor": a_snap,
            "max_velocity_vector_relative_error": max_vel_rel,
            "rms_velocity_vector_relative_error": rms_vel_rel,
            "max_speed_scaling_relative_error": scale_err,
            "q_displacement_spearman": rho,
            "mass_max_relative_change": mass_rel,
            "mean_displacement_Mpc": float(np.mean(disp)),
            "max_displacement_Mpc": float(np.max(disp)),
            "mean_input_speed_km_s": float(np.mean(speed0)),
            "mean_output_speed_km_s": float(np.mean(speed)),
            "weights_unity_if_present": None if weights is None else True,
            "sampled_speeds_zero_if_present": None if sampled is None else True,
        }
        final_speed_by_state[name] = speed
        input_q_by_state[name] = q0
        displacement_by_state[name] = disp

    # All states must use the same macro masses.
    means = []
    for name in ("FD", "Fplus", "Fminus"):
        ref = np.load(SCI / name / "input_reference.npz")
        means.append(float(np.mean(ref["masses"])))
    cross_mass_gap = (max(means) - min(means)) / max(np.mean(means), 1e-300)
    need(cross_mass_gap <= 1e-12, "science state macro masses differ")

    # Pairwise state speed ordering must preserve the ordering implied by q for
    # every particle whose input q values are not numerically tied.
    pairs = (("FD", "Fplus"), ("FD", "Fminus"), ("Fplus", "Fminus"))
    ordering = {}
    ordering_all = True
    for a, b in pairs:
        dq = input_q_by_state[a] - input_q_by_state[b]
        ds = final_speed_by_state[a] - final_speed_by_state[b]
        mask = np.abs(dq) > 1e-13
        good = np.sign(dq[mask]) == np.sign(ds[mask])
        frac = 1.0 if np.sum(mask) == 0 else float(np.mean(good))
        ordering[f"{a}_vs_{b}"] = {
            "non_tied_particles": int(np.sum(mask)),
            "fraction_matching": frac,
        }
        ordering_all &= bool(frac == 1.0)
    need(ordering_all, "final state speed ordering changed")

    # State-specific transport must not collapse to one identical trajectory.
    pair_disp_rms = {}
    distinct = True
    for a, b in pairs:
        gap = displacement_by_state[a] - displacement_by_state[b]
        rms = float(np.sqrt(np.mean(gap * gap)))
        pair_disp_rms[f"{a}_vs_{b}"] = rms
        distinct &= rms > 1e-8
    need(distinct, "state-specific displacement became numerically identical")

    result = {
        "stage": p["stage"],
        "status": "PASS_SWIFT_CUSTOM_F_COSMOLOGICAL_TRANSPORT_SMOKE",
        "science_gate_pass": True,
        "swift_commit": p["swift"]["commit"],
        "pre_execution_geometry_amendment": p["pre_execution_geometry_amendment"],
        "calibration": json.loads(CAL_MASS.read_text()),
        "expected_force_free_velocity_scale": expected_scale,
        "states": state,
        "cross_state_macro_mass_relative_gap": cross_mass_gap,
        "state_speed_ordering": ordering,
        "state_speed_ordering_all_match": ordering_all,
        "pairwise_displacement_rms_difference_Mpc": pair_disp_rms,
        "state_specific_displacement_not_identical": distinct,
        "interpretation": p["interpretation"],
        "guards": p["guards"],
        "next_stage": p["next_if_pass"],
    }
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")

    print("K01D1_EXPECTED_VELOCITY_SCALE", expected_scale, flush=True)
    for name in ("FD", "Fplus", "Fminus"):
        print(
            "K01D1_STATE",
            name,
            "VEL_MAX_REL", state[name]["max_velocity_vector_relative_error"],
            "SPEARMAN", state[name]["q_displacement_spearman"],
            "MEAN_DX", state[name]["mean_displacement_Mpc"],
            flush=True,
        )
    print("K01D1_STATE_ORDERING_PASS", ordering_all, flush=True)
    print("K01D1_DISTINCT_TRANSPORT_PASS", distinct, pair_disp_rms, flush=True)
    print("K01D1_SWIFT_CUSTOM_F_COSMOLOGICAL_TRANSPORT_SMOKE_PASS", flush=True)
    print("K01D1_RESULT", OUT, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "mode",
        choices=("prepare-calibration", "extract-calibration", "prepare-science", "analyze"),
    )
    args = ap.parse_args()
    p = prereg()
    WORK.mkdir(parents=True, exist_ok=True)
    if args.mode == "prepare-calibration":
        prepare_calibration(p)
    elif args.mode == "extract-calibration":
        extract_calibration(p)
    elif args.mode == "prepare-science":
        prepare_science(p)
    else:
        analyze(p)


if __name__ == "__main__":
    main()
