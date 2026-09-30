#!/usr/bin/env python3
"""
E62: Quijote z=1 periodic-box truth calibration of the E59/E61
density-to-velocity reconstruction.

This stage is deliberately independent of eBOSS observations. It uses one
complete periodic N-body halo catalogue with true halo peculiar velocities to
test the reconstruction method itself.

Frozen primary design:
  - Quijote fiducial realization 0, FoF snapshot 2 (nominal z=1);
  - exactly 165107 most massive FoF halos, unit number weights;
  - exact E59 top-hat-smoothed velocity kernel, R=16 Mpc/h;
  - primary cutoff 256 Mpc/h, check cutoff 192 Mpc/h;
  - 4096 deterministic tracer-halo probes;
  - compare reconstructed x/y/z components with true halo peculiar velocities;
  - all three axes must satisfy the preregistered truth and cutoff gates.

A PASS authorizes only a later end-to-end survey/lightcone truth-transfer test.
It does not calibrate the Einstein-Vlasov amplitude and does not authorize
opening observed eBOSS galaxies or the observed odd vector.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "source_data/e61_elg_fullpool_reconstruction_lock_compact_summary_2026-09-30.json"
OUT = ROOT / "source_data/e62_quijote_z1_truth_velocity_calibration_result.json"

SNAPNUM = 2
Z = 1.0
BOX = 1000.0
NTRACER = 165107
NPROBE = 4096
SEED = 202609620000

R_SMOOTH = 16.0
RMAX_PRIMARY = 256.0
RMAX_CHECK = 192.0

TRUTH_R_MIN = 0.70
TRUTH_SIGN_MIN = 0.70
CUTOFF_R_MIN = 0.90
CUTOFF_SIGN_MIN = 0.90


def need(c, msg):
    if not c:
        raise RuntimeError(msg)


def atomic(path: Path, obj):
    raw = (json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".e62_", delete=False) as f:
        q = Path(f.name)
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())
    try:
        os.replace(q, path)
    finally:
        q.unlink(missing_ok=True)


def arr_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def minimum_image(delta, box=BOX):
    return delta - box * np.rint(delta / box)


def top_hat_velocity_kernel(delta, R=R_SMOOTH):
    r2 = np.einsum("ij,ij->i", delta, delta)
    r = np.sqrt(r2)
    fac = np.empty_like(r)
    inside = r < R
    fac[inside] = 1.0 / (R ** 3)
    outside = ~inside
    fac[outside] = 1.0 / np.maximum(r[outside] ** 3, 1e-300)
    return delta * fac[:, None], r


def pearson(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    aa = a - np.mean(a)
    bb = b - np.mean(b)
    den = float(np.linalg.norm(aa) * np.linalg.norm(bb))
    return float(aa @ bb / den) if den > 0 else 0.0


def sign_agreement(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    nz = (a != 0) & (b != 0)
    need(np.count_nonzero(nz) >= 0.95 * len(a), "Too many exact-zero values")
    return float(np.mean(np.sign(a[nz]) == np.sign(b[nz])))


def reconstruct_periodic(source_pos, probe_pos, chunk=32):
    tree = cKDTree(source_pos, boxsize=BOX, leafsize=32)
    primary = np.zeros((len(probe_pos), 3), float)
    check = np.zeros((len(probe_pos), 3), float)

    for a0 in range(0, len(probe_pos), chunk):
        a1 = min(a0 + chunk, len(probe_pos))
        nbs = tree.query_ball_point(probe_pos[a0:a1], RMAX_PRIMARY, workers=2)
        for local, js in enumerate(nbs):
            if not len(js):
                continue
            j = np.asarray(js, dtype=np.int64)
            delta = minimum_image(source_pos[j] - probe_pos[a0 + local])
            K, r = top_hat_velocity_kernel(delta)
            primary[a0 + local] = K.sum(axis=0)
            m = r <= RMAX_CHECK
            if np.any(m):
                check[a0 + local] = K[m].sum(axis=0)

    # Overall positive normalization is irrelevant to Pearson/sign, but divide
    # by the fixed tracer count so stored values remain numerically tame.
    primary /= float(len(source_pos))
    check /= float(len(source_pos))
    return primary, check


def self_test():
    d = np.array([[490.0, 0.0, 0.0], [-490.0, 0.0, 0.0]])
    w = minimum_image(d, 1000.0)
    need(np.allclose(w[:, 0], [490.0, -490.0]), "minimum-image interior fail")

    d2 = np.array([[510.0, 0.0, 0.0], [-510.0, 0.0, 0.0]])
    w2 = minimum_image(d2, 1000.0)
    need(np.allclose(w2[:, 0], [-490.0, 490.0]), "minimum-image wrap fail")

    q = np.array([[0., 0., 0.], [8., 0., 0.], [16., 0., 0.], [32., 0., 0.]])
    K, r = top_hat_velocity_kernel(q)
    need(np.allclose(K[0], 0), "kernel r=0 fail")
    need(np.isclose(K[1, 0], 8.0 / R_SMOOTH**3), "kernel inside fail")
    need(np.isclose(K[2, 0], 16.0 / 16.0**3), "kernel boundary fail")
    need(np.isclose(K[3, 0], 32.0 / 32.0**3), "kernel outside fail")

    a = np.array([1., -2., 3., -4.])
    b = 3.0 * a
    need(abs(pearson(a, b) - 1.0) < 1e-14, "Pearson self-test fail")
    need(sign_agreement(a, b) == 1.0, "sign self-test fail")
    print("E62_SYNTHETIC_PERIODIC_TRUTH_CALIBRATION_SELF_TEST_PASS", flush=True)


def load_quijote_fof(catalog_dir: Path):
    try:
        import readfof
    except Exception as e:
        raise RuntimeError(
            "Missing readfof/Pylians. Install Pylians3 in the project venv before E62."
        ) from e

    need(catalog_dir.is_dir(), f"Missing Quijote FoF directory: {catalog_dir}")

    fof = readfof.FoF_catalog(
        str(catalog_dir), SNAPNUM,
        long_ids=False, swap=False, SFR=False, read_IDs=False
    )

    pos = np.asarray(fof.GroupPos, dtype=np.float64) / 1e3
    mass = np.asarray(fof.GroupMass, dtype=np.float64) * 1e10
    vel = np.asarray(fof.GroupVel, dtype=np.float64) * (1.0 + Z)
    npart = np.asarray(fof.GroupLen)

    need(pos.ndim == 2 and pos.shape[1] == 3, "Bad GroupPos shape")
    need(vel.shape == pos.shape, "Bad GroupVel shape")
    need(len(mass) == len(pos) == len(npart), "FoF array-length mismatch")
    need(np.isfinite(pos).all() and np.isfinite(vel).all() and np.isfinite(mass).all(),
         "Nonfinite FoF values")
    need(np.all((pos >= -1e-6) & (pos <= BOX + 1e-6)), "Position outside fixed 1 Gpc/h box")
    # cKDTree(boxsize=BOX) requires [0,BOX); wrap only possible roundoff-level edge values.
    pos = np.mod(pos, BOX)
    need(len(pos) >= NTRACER, f"FoF catalog has only {len(pos)} halos < frozen {NTRACER}")

    return pos, mass, vel, npart


def choose_tracer(pos, mass, vel, npart):
    # Stable descending mass rank, then original row index as deterministic tie break.
    idx = np.lexsort((np.arange(len(mass), dtype=np.int64), -mass))[:NTRACER]
    return (
        np.asarray(pos[idx], float),
        np.asarray(mass[idx], float),
        np.asarray(vel[idx], float),
        np.asarray(npart[idx]),
        np.asarray(idx, np.int64),
    )


def axis_metrics(rec, chk, truth):
    out = {}
    for j, axis in enumerate(("x", "y", "z")):
        rtruth = pearson(rec[:, j], truth[:, j])
        struth = sign_agreement(rec[:, j], truth[:, j])
        rcut = pearson(rec[:, j], chk[:, j])
        scut = sign_agreement(rec[:, j], chk[:, j])
        out[axis] = {
            "truth_pearson": rtruth,
            "truth_sign_agreement": struth,
            "cutoff_192_vs_256_pearson": rcut,
            "cutoff_192_vs_256_sign_agreement": scut,
            "reconstruction_rms": float(np.sqrt(np.mean(rec[:, j] ** 2))),
            "truth_velocity_rms_km_s": float(np.sqrt(np.mean(truth[:, j] ** 2))),
            "gate_pass": bool(
                rtruth >= TRUTH_R_MIN
                and struth >= TRUTH_SIGN_MIN
                and rcut >= CUTOFF_R_MIN
                and scut >= CUTOFF_SIGN_MIN
            ),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--catalog-dir", type=Path)
    ap.add_argument("--chunk", type=int, default=32)
    args = ap.parse_args()

    self_test()
    if args.self_test:
        return

    need(args.run, "Use --run for the real E62 truth calibration")
    need(args.catalog_dir is not None, "--catalog-dir is required")
    need(8 <= args.chunk <= 128, "Unsafe chunk")
    need(PARENT.is_file(), "Missing E61 compact parent")

    parent = json.loads(PARENT.read_text())
    need(parent["status"] == "PASS_ELG_FULLPOOL_RECONSTRUCTION_LOCK_QUANTIFIED",
         "E61 status changed")
    need(parent["decision"]["full_pool_reconstruction_operator_locked"] is True,
         "E61 full-pool lock changed")
    need(parent["observed_odd_used"] is False and parent["observed_galaxy_rows_used"] is False,
         "E61 observation guardrail changed")

    pos, mass, vel, npart = load_quijote_fof(args.catalog_dir)
    tpos, tmass, tvel, tnpart, tidx = choose_tracer(pos, mass, vel, npart)

    rng = np.random.default_rng(SEED)
    pidx = np.sort(rng.choice(NTRACER, size=NPROBE, replace=False))
    probes = tpos[pidx]
    truth = tvel[pidx]

    rec, chk = reconstruct_periodic(tpos, probes, chunk=args.chunk)
    metrics = axis_metrics(rec, chk, truth)
    passed = all(metrics[a]["gate_pass"] for a in ("x", "y", "z"))

    state = {
        "stage": "E62_QUIJOTE_Z1_TRUTH_VELOCITY_CALIBRATION",
        "date": "2026-09-30",
        "status": "PASS_TRUTH_CALIBRATION_QUANTIFIED" if passed else "FAIL_TRUTH_CALIBRATION_QUANTIFIED",
        "dataset": {
            "suite": "Quijote",
            "cosmology": "fiducial",
            "realization": 0,
            "halo_catalog": "FoF",
            "snapshot_number": SNAPNUM,
            "nominal_redshift": Z,
            "box_size_Mpc_over_h": BOX,
            "catalog_dir": str(args.catalog_dir),
            "catalog_total_halos": int(len(pos)),
        },
        "tracer": {
            "policy": "exactly 165107 most massive FoF halos",
            "n": NTRACER,
            "minimum_selected_mass_Msun_over_h": float(np.min(tmass)),
            "maximum_selected_mass_Msun_over_h": float(np.max(tmass)),
            "minimum_selected_particle_count": int(np.min(tnpart)),
            "selected_source_row_indices_SHA256": arr_sha(tidx),
            "selected_positions_SHA256": arr_sha(tpos),
            "selected_velocities_SHA256": arr_sha(tvel),
        },
        "probes": {
            "n": NPROBE,
            "seed": SEED,
            "selected_tracer_indices_SHA256": arr_sha(pidx),
            "positions_SHA256": arr_sha(probes),
            "truth_velocities_SHA256": arr_sha(truth),
        },
        "reconstruction": {
            "smoothing_R_Mpc_over_h": R_SMOOTH,
            "primary_cutoff_Mpc_over_h": RMAX_PRIMARY,
            "check_cutoff_Mpc_over_h": RMAX_CHECK,
            "periodic_minimum_image": True,
            "unit_number_weights": True,
            "same_kernel_as_E59": True,
        },
        "gates": {
            "each_axis_truth_pearson_min": TRUTH_R_MIN,
            "each_axis_truth_sign_agreement_min": TRUTH_SIGN_MIN,
            "each_axis_cutoff_pearson_min": CUTOFF_R_MIN,
            "each_axis_cutoff_sign_agreement_min": CUTOFF_SIGN_MIN,
        },
        "axes": metrics,
        "decision": {
            "truth_velocity_gate_pass": bool(passed),
            "interpretation": (
                "PERIODIC_NBODY_TRUTH_CALIBRATION_PASS; NEXT_END_TO_END_SURVEY_LIGHTCONE_TRUTH_TRANSFER"
                if passed else
                "PERIODIC_NBODY_TRUTH_CALIBRATION_FAIL; DO_NOT_USE_OBSERVED_VELOCITY_TAG"
            ),
            "important_limit": (
                "Quijote z=1 periodic-box calibration tests the reconstruction method, not "
                "eBOSS cut-sky selection transfer, ELG HOD realism, custom F+/F- neutrino wakes, "
                "or absolute Einstein-Vlasov amplitude."
            ),
        },
        "observed_galaxy_rows_used": False,
        "observed_odd_used": False,
        "absolute_EV_amplitude_calibrated": False,
        "covariance_inverse_used": False,
        "pvalue_or_detection_sigma": False,
    }
    atomic(OUT, state)

    print("E62_WSL_QUIJOTE_TRUTH_CALIBRATION_COMPLETE", flush=True)
    for axis in ("x", "y", "z"):
        m = metrics[axis]
        print(
            "AXIS", axis,
            "R_TRUTH", m["truth_pearson"],
            "SIGN_TRUTH", m["truth_sign_agreement"],
            "R_CUTOFF", m["cutoff_192_vs_256_pearson"],
            "SIGN_CUTOFF", m["cutoff_192_vs_256_sign_agreement"],
            "PASS", m["gate_pass"],
            flush=True,
        )
    print("DECISION", state["decision"]["interpretation"], flush=True)
    print("OBSERVED_ODD_USED", False, flush=True)
    print("REPORT", OUT, flush=True)


if __name__ == "__main__":
    main()
