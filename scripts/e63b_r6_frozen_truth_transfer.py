#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]

LOCK = ROOT/"eboss_workspace/cosmodc2/e63b_r6_final_pretruth_lock.json"
OUT = ROOT/"source_data/e63b_r6_outerrim_buffered_lightcone_truth_transfer_result.json"

EXPECTED_CACHE_SHA = \
    "74d0402bfb2dc87b7ba5f174fb43adb141db4dc21c7f7b1b8acc9a20eefa1f82"

EXPECTED_MEMBERSHIP_SHA = \
    "112836101c35286f803f99954969d2512e55a746c8fbeae84839b486134b1dd6"

RA0 = 62.5
DEC0 = -35.5

OUTER_RADIUS_DEG = 9.5
INNER_RADIUS_DEG = 2.5

SOURCE_Z0 = 0.75
SOURCE_Z1 = 1.16

PROBE_Z0 = 0.9
PROBE_Z1 = 1.0

MASS_FLOOR = 5.0e12

NSEL = 54893
NPROBE = 512
PROBE_SEED = 202610630001

OMEGA_M = 0.2648
OMEGA_L = 0.7352
H = 0.71
C_KM_S = 299792.458

N_TARGET = 1.65107e-4

R_SMOOTH = 16.0
RMAX_PRIMARY = 256.0
RMAX_CHECK = 192.0

TRUTH_R_MIN = 0.70
TRUTH_SIGN_MIN = 0.70
CUTOFF_R_MIN = 0.90
CUTOFF_SIGN_MIN = 0.90


def need(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def sha_file(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024*1024), b""):
            h.update(b)
    return h.hexdigest()


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def atomic(path, obj):
    raw = (
        json.dumps(
            obj,
            indent=2,
            sort_keys=True,
            allow_nan=False
        ) + "\n"
    ).encode()

    path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        dir=path.parent,
        prefix=".e63b_r6_",
        delete=False
    ) as f:
        tmp = Path(f.name)
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())

    try:
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def parse_bool(v):
    s = str(v).strip().lower()

    if s in ("true", "1", "t", "yes"):
        return True

    if s in ("false", "0", "f", "no"):
        return False

    raise RuntimeError(f"Unexpected bool {v!r}")


def chi_h(z):
    val = quad(
        lambda zz:
            1.0/math.sqrt(
                OMEGA_M*(1.0+zz)**3 + OMEGA_L
            ),
        0.0,
        float(z),
        epsabs=1e-10,
        epsrel=1e-10,
        limit=100
    )[0]

    return (C_KM_S/100.0)*val


def chi_mpc(z):
    return chi_h(z)/H


def solid_angle_circle(radius_deg):
    t = math.radians(radius_deg)
    return 2.0*math.pi*(1.0-math.cos(t))


def angular_sep_deg(ra, dec):
    r = np.deg2rad(np.asarray(ra, float))
    d = np.deg2rad(np.asarray(dec, float))

    r0 = math.radians(RA0)
    d0 = math.radians(DEC0)

    c = (
        np.sin(d)*math.sin(d0)
        + np.cos(d)*math.cos(d0)*np.cos(r-r0)
    )

    return np.rad2deg(
        np.arccos(np.clip(c, -1.0, 1.0))
    )


def position_angle_deg(ra, dec):
    r = np.deg2rad(np.asarray(ra, float))
    d = np.deg2rad(np.asarray(dec, float))

    r0 = math.radians(RA0)
    d0 = math.radians(DEC0)

    dr = r-r0

    y = np.sin(dr)*np.cos(d)

    x = (
        math.cos(d0)*np.sin(d)
        - math.sin(d0)*np.cos(d)*np.cos(dr)
    )

    return (
        np.rad2deg(np.arctan2(y, x)) + 360.0
    ) % 360.0


def sky_unit(ra, dec):
    r = np.deg2rad(np.asarray(ra, float))
    d = np.deg2rad(np.asarray(dec, float))

    return np.column_stack([
        np.cos(d)*np.cos(r),
        np.cos(d)*np.sin(r),
        np.sin(d)
    ])


def pearson(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)

    aa = a-np.mean(a)
    bb = b-np.mean(b)

    den = float(
        np.linalg.norm(aa)*np.linalg.norm(bb)
    )

    return float(aa@bb/den) if den > 0 else 0.0


def sign_agreement(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)

    nz = (a != 0) & (b != 0)

    need(
        np.count_nonzero(nz) >= 0.95*len(a),
        "Too many zero values"
    )

    return float(
        np.mean(np.sign(a[nz]) == np.sign(b[nz]))
    )


def vector_cosines(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)

    den = (
        np.linalg.norm(a, axis=1)
        * np.linalg.norm(b, axis=1)
    )

    m = den > 0

    need(
        np.count_nonzero(m) >= 0.95*len(a),
        "Too many zero vectors"
    )

    return (
        np.einsum("ij,ij->i", a[m], b[m])
        / den[m]
    )


def top_hat_kernel(delta):
    r2 = np.einsum("ij,ij->i", delta, delta)
    r = np.sqrt(r2)

    fac = np.empty_like(r)

    m = r < R_SMOOTH

    fac[m] = 1.0/(R_SMOOTH**3)
    fac[~m] = 1.0/np.maximum(r[~m]**3, 1e-300)

    return delta*fac[:, None], r


def reconstruct(source_pos, probe_pos, chunk):
    print(
        "E63B_R6_RECONSTRUCTION_TREE_START",
        "SOURCES", len(source_pos),
        flush=True
    )

    tree = cKDTree(source_pos, leafsize=32)

    primary = np.zeros(
        (len(probe_pos), 3),
        dtype=float
    )

    check = np.zeros_like(primary)

    done = 0

    for a0 in range(0, len(probe_pos), chunk):
        a1 = min(a0+chunk, len(probe_pos))

        nbs = tree.query_ball_point(
            probe_pos[a0:a1],
            RMAX_PRIMARY,
            workers=1
        )

        for local, js in enumerate(nbs):
            if not len(js):
                continue

            j = np.asarray(js, dtype=np.int64)

            delta = (
                source_pos[j]
                - probe_pos[a0+local]
            )

            K, r = top_hat_kernel(delta)

            primary[a0+local] = K.sum(axis=0)

            m = r <= RMAX_CHECK

            if np.any(m):
                check[a0+local] = K[m].sum(axis=0)

        done = a1

        if (
            done == len(probe_pos)
            or done % 32 == 0
        ):
            print(
                "E63B_R6_RECONSTRUCTION_PROGRESS",
                done,
                "/",
                len(probe_pos),
                flush=True
            )

    primary /= float(len(source_pos))
    check /= float(len(source_pos))

    return primary, check


def metrics(rec, chk, truth, pos):
    los = pos/np.linalg.norm(pos, axis=1)[:, None]

    rp = np.einsum("ij,ij->i", rec, los)
    rc = np.einsum("ij,ij->i", chk, los)
    tv = np.einsum("ij,ij->i", truth, los)

    out = {
        "LOS": {
            "truth_pearson":
                pearson(rp, tv),

            "truth_sign_agreement":
                sign_agreement(rp, tv),

            "cutoff_192_vs_256_pearson":
                pearson(rp, rc),

            "cutoff_192_vs_256_sign_agreement":
                sign_agreement(rp, rc),
        },

        "cartesian": {},
    }

    for j, axis in enumerate(("x", "y", "z")):
        out["cartesian"][axis] = {
            "truth_pearson":
                pearson(rec[:, j], truth[:, j]),

            "truth_sign_agreement":
                sign_agreement(
                    rec[:, j],
                    truth[:, j]
                ),

            "cutoff_192_vs_256_pearson":
                pearson(rec[:, j], chk[:, j]),

            "cutoff_192_vs_256_sign_agreement":
                sign_agreement(
                    rec[:, j],
                    chk[:, j]
                ),
        }

    vc = vector_cosines(rec, truth)

    out["vector_cosine"] = {
        "mean": float(np.mean(vc)),
        "median": float(np.median(vc)),
        "min": float(np.min(vc)),
        "max": float(np.max(vc)),
    }

    L = out["LOS"]

    L["gate_pass"] = bool(
        L["truth_pearson"] >= TRUTH_R_MIN
        and
        L["truth_sign_agreement"] >= TRUTH_SIGN_MIN
        and
        L["cutoff_192_vs_256_pearson"] >= CUTOFF_R_MIN
        and
        L["cutoff_192_vs_256_sign_agreement"] >= CUTOFF_SIGN_MIN
    )

    return out


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--run",
        action="store_true"
    )

    ap.add_argument(
        "--chunk",
        type=int,
        default=1
    )

    args = ap.parse_args()

    need(args.run, "Use --run")
    need(1 <= args.chunk <= 32, "Unsafe chunk")

    # --------------------------------------------------------
    # FINAL PRE-TRUTH LOCK
    # --------------------------------------------------------

    need(
        LOCK.is_file(),
        "Missing final R6 pre-truth lock"
    )

    lock = json.loads(LOCK.read_text())

    need(
        lock["stage"] ==
        "E63B_R6_FINAL_PRETRUTH_LOCK",
        "Wrong lock"
    )

    need(
        lock["truth_velocity_used_before_lock"] is False
        and
        lock["science_metrics_computed_before_lock"] is False
        and
        lock["observed_eBOSS_rows_used"] is False
        and
        lock["observed_odd_used"] is False,
        "Pre-truth guardrail changed"
    )

    fs = lock["frozen_science"]

    need(fs["NSEL"] == NSEL, "NSEL drift")
    need(fs["NPROBE"] == NPROBE, "NPROBE drift")
    need(fs["probe_seed"] == PROBE_SEED, "seed drift")

    need(
        fs["truth_pearson_min"] == TRUTH_R_MIN
        and
        fs["truth_sign_min"] == TRUTH_SIGN_MIN
        and
        fs["cutoff_pearson_min"] == CUTOFF_R_MIN
        and
        fs["cutoff_sign_min"] == CUTOFF_SIGN_MIN,
        "Science gate drift"
    )

    CACHE = ROOT/lock["retrieval"]["cache"]
    MEM = ROOT/lock["probe_boundary"]["membership"]

    need(CACHE.is_file(), "Missing frozen R6 cache")
    need(MEM.is_file(), "Missing boundary membership")

    cache_sha = sha_file(CACHE)
    mem_sha = sha_file(MEM)
    lock_sha = sha_file(LOCK)

    need(
        cache_sha == EXPECTED_CACHE_SHA,
        "Frozen R6 cache SHA mismatch"
    )

    need(
        mem_sha == EXPECTED_MEMBERSHIP_SHA,
        "Boundary membership SHA mismatch"
    )

    need(
        lock["retrieval"]["exact_cone_unique_centrals"]
        == 200493,
        "candidate-count drift"
    )

    # --------------------------------------------------------
    # MEMBERSHIP LOCK
    # --------------------------------------------------------

    membership = {}

    with MEM.open(newline="") as f:
        for row in csv.DictReader(f):
            gid = int(row["galaxy_id"])
            inside = bool(
                int(row["server_inside_probe_z"])
            )

            need(
                gid not in membership,
                "duplicate membership galaxy_id"
            )

            membership[gid] = inside

    need(
        len(membership) == 4,
        "boundary membership size drift"
    )

    need(
        sum(membership.values()) == 2,
        "boundary membership inside-count drift"
    )

    # --------------------------------------------------------
    # LOAD FROZEN R6 CACHE
    # --------------------------------------------------------

    rows = []

    with CACHE.open(newline="") as f:
        r = csv.DictReader(f)
        rows = list(r)

    need(
        len(rows) == 200493,
        f"R6 cache row drift: {len(rows)}"
    )

    gid = np.asarray(
        [int(x["galaxy_id"]) for x in rows],
        dtype=np.int64
    )

    halo = np.asarray(
        [int(x["halo_id"]) for x in rows],
        dtype=np.int64
    )

    ra = np.asarray(
        [float(x["ra_true"]) for x in rows]
    )

    dec = np.asarray(
        [float(x["dec_true"]) for x in rows]
    )

    z = np.asarray(
        [float(x["redshift_true"]) for x in rows]
    )

    mass = np.asarray(
        [float(x["halo_mass"]) for x in rows]
    )

    central = np.asarray(
        [parse_bool(x["is_central"]) for x in rows],
        dtype=bool
    )

    need(np.all(central), "Non-central row in R6 cache")

    need(
        len(np.unique(gid)) == len(gid),
        "duplicate galaxy_id"
    )

    need(
        len(np.unique(halo)) == len(halo),
        "duplicate central halo_id"
    )

    need(
        np.all(np.isfinite(ra))
        and np.all(np.isfinite(dec))
        and np.all(np.isfinite(z))
        and np.all(np.isfinite(mass)),
        "nonfinite basic cache values"
    )

    # --------------------------------------------------------
    # GEOMETRY / COVERAGE RECHECK
    # --------------------------------------------------------

    aang_all = angular_sep_deg(ra, dec)

    need(
        float(np.max(aang_all))
        <= OUTER_RADIUS_DEG + 1e-10,
        "R6 cache leaked outside exact cone"
    )

    pa = position_angle_deg(ra, dec)

    ann = (
        (aang_all >= 8.5)
        & (aang_all < 9.5)
    )

    coverage = []

    for j in range(24):
        lo = 15.0*j
        hi = 15.0*(j+1)

        coverage.append(
            int(
                np.count_nonzero(
                    ann
                    & (pa >= lo)
                    & (pa < hi)
                )
            )
        )

    need(
        min(coverage) >= 50,
        f"Coverage fail {coverage}"
    )

    need(
        min(coverage)
        ==
        lock["retrieval"]["coverage_min_sector_count"],
        "Coverage provenance drift"
    )

    print(
        "E63B_R6_FROZEN_COVERAGE_PASS",
        "MIN", min(coverage),
        flush=True
    )

    # --------------------------------------------------------
    # POSITION UNIT GUARD
    # --------------------------------------------------------

    take = np.linspace(
        0,
        len(rows)-1,
        1024,
        dtype=int
    )

    pos_guard = np.column_stack([
        [float(rows[i]["position_x"]) for i in take],
        [float(rows[i]["position_y"]) for i in take],
        [float(rows[i]["position_z"]) for i in take],
    ])

    rr = np.linalg.norm(pos_guard, axis=1)

    ch = np.asarray(
        [chi_mpc(z[i]) for i in take],
        dtype=float
    )

    ratio = float(np.median(rr/ch))

    if 0.97 <= ratio <= 1.03:
        scale = H
        unit_mode = "IRSA_GCR_Mpc_to_Mpc_over_h"

    elif H-0.025 <= ratio <= H+0.025:
        scale = 1.0
        unit_mode = "already_Mpc_over_h"

    else:
        raise RuntimeError(
            f"Unrecognized position units ratio={ratio}"
        )

    u = (
        pos_guard*scale
        /
        np.linalg.norm(
            pos_guard*scale,
            axis=1
        )[:, None]
    )

    us = sky_unit(
        ra[take],
        dec[take]
    )

    dots = np.einsum("ij,ij->i", u, us)

    need(
        float(np.median(dots)) >= 0.999
        and
        float(np.min(dots)) >= 0.995,
        "position/sky mismatch"
    )

    print(
        "E63B_R6_POSITION_UNIT_GUARD_PASS",
        unit_mode,
        "RATIO", ratio,
        flush=True
    )

    # --------------------------------------------------------
    # FROZEN TARGET DENSITY / SELECTION
    # --------------------------------------------------------

    omega = solid_angle_circle(OUTER_RADIUS_DEG)

    r0 = chi_h(SOURCE_Z0)
    r1 = chi_h(SOURCE_Z1)

    volume = omega*(r1**3-r0**3)/3.0

    nsel_calc = int(round(N_TARGET*volume))

    need(
        nsel_calc == NSEL,
        f"Frozen NSEL mismatch {nsel_calc}"
    )

    order = np.lexsort((gid, -mass))
    selected = order[:NSEL]

    minmass = float(np.min(mass[selected]))

    need(
        minmass >= 1.10*MASS_FLOOR,
        f"Retrieval floor too close: {minmass}"
    )

    selected_gid = gid[selected]

    print(
        "E63B_R6_FROZEN_TRACER_SELECTION_PASS",
        "NSEL", len(selected),
        "MIN_MASS", minmass,
        flush=True
    )

    # --------------------------------------------------------
    # PROBE MEMBERSHIP
    #
    # Interior serialized z values are direct.
    # Serialized 0.9000 / 1.0000 are resolved ONLY by the
    # frozen server-side membership lock.
    # --------------------------------------------------------

    aang = aang_all[selected]
    zang = z[selected]

    inner = aang <= INNER_RADIUS_DEG

    strict_z = (
        (zang > PROBE_Z0)
        & (zang < PROBE_Z1)
    )

    pmask = inner & strict_z

    boundary = (
        inner
        &
        (
            (zang == PROBE_Z0)
            |
            (zang == PROBE_Z1)
        )
    )

    boundary_local = np.flatnonzero(boundary)

    boundary_included = 0
    boundary_excluded = 0

    for ii in boundary_local:
        g = int(selected_gid[ii])

        need(
            g in membership,
            f"Boundary galaxy {g} missing membership lock"
        )

        if membership[g]:
            pmask[ii] = True
            boundary_included += 1
        else:
            pmask[ii] = False
            boundary_excluded += 1

    pool = np.flatnonzero(pmask)

    need(
        len(pool) >= NPROBE,
        f"Only {len(pool)} eligible probes"
    )

    print(
        "E63B_R6_PROBE_POOL_LOCK_APPLIED",
        "ELIGIBLE", len(pool),
        "SELECTED_BOUNDARY", len(boundary_local),
        "BOUNDARY_INCLUDED", boundary_included,
        "BOUNDARY_EXCLUDED", boundary_excluded,
        flush=True
    )

    # --------------------------------------------------------
    # DETERMINISTIC 512 PROBES
    # --------------------------------------------------------

    rng = np.random.default_rng(PROBE_SEED)

    pidx = np.sort(
        rng.choice(
            pool,
            size=NPROBE,
            replace=False
        )
    )

    probe_gid = selected_gid[pidx]

    # --------------------------------------------------------
    # NOW AND ONLY NOW LOAD THE VELOCITY TRUTH
    # --------------------------------------------------------

    pos_all = np.column_stack([
        np.asarray(
            [float(x["position_x"]) for x in rows]
        ),
        np.asarray(
            [float(x["position_y"]) for x in rows]
        ),
        np.asarray(
            [float(x["position_z"]) for x in rows]
        ),
    ]) * scale

    vel_all = np.column_stack([
        np.asarray(
            [float(x["velocity_x"]) for x in rows]
        ),
        np.asarray(
            [float(x["velocity_y"]) for x in rows]
        ),
        np.asarray(
            [float(x["velocity_z"]) for x in rows]
        ),
    ])

    need(
        np.isfinite(pos_all).all(),
        "nonfinite position"
    )

    need(
        np.isfinite(vel_all).all(),
        "nonfinite velocity"
    )

    source_pos = pos_all[selected]
    source_vel = vel_all[selected]

    ppos = source_pos[pidx]
    ptruth = source_vel[pidx]

    # --------------------------------------------------------
    # BUFFER GUARDS
    # --------------------------------------------------------

    pr = np.linalg.norm(ppos, axis=1)

    angular_buffer = (
        pr
        * math.sin(
            math.radians(
                OUTER_RADIUS_DEG-INNER_RADIUS_DEG
            )
        )
    )

    lower_buffer = pr-r0
    upper_buffer = r1-pr

    need(
        float(np.min(angular_buffer))
        > RMAX_PRIMARY,
        "Angular buffer below 256"
    )

    need(
        float(np.min(lower_buffer))
        > RMAX_PRIMARY,
        "Lower radial buffer below 256"
    )

    need(
        float(np.min(upper_buffer))
        > RMAX_PRIMARY,
        "Upper radial buffer below 256"
    )

    print(
        "E63B_R6_ALL_PRETRUTH_GUARDS_PASS",
        flush=True
    )

    print(
        "E63B_R6_RECONSTRUCTION_START",
        "SOURCES", len(source_pos),
        "PROBES", len(ppos),
        "CHUNK", args.chunk,
        "WORKERS", 1,
        flush=True
    )

    # --------------------------------------------------------
    # FIRST FROZEN TRUTH EXECUTION
    # --------------------------------------------------------

    rec, chk = reconstruct(
        source_pos,
        ppos,
        args.chunk
    )

    met = metrics(
        rec,
        chk,
        ptruth,
        ppos
    )

    passed = bool(
        met["LOS"]["gate_pass"]
    )

    result = {
        "stage":
            "E63B_R6_OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER",

        "date":
            "2026-10-14",

        "status":
            (
                "PASS_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED"
                if passed
                else
                "FAIL_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED"
            ),

        "frozen_pretruth_lock": {
            "path": str(LOCK.relative_to(ROOT)),
            "sha256": lock_sha,
        },

        "source": {
            "cache":
                str(CACHE.relative_to(ROOT)),

            "cache_sha256":
                cache_sha,

            "candidate_rows":
                len(rows),

            "selected_count":
                len(selected),

            "selected_galaxy_id_sha256":
                sha_bytes(
                    np.ascontiguousarray(
                        selected_gid
                    ).tobytes()
                ),
        },

        "geometry": {
            "center_ra_deg": RA0,
            "center_dec_deg": DEC0,
            "outer_radius_deg":
                OUTER_RADIUS_DEG,

            "outer_geometric_volume_Mpc_over_h_cubed":
                volume,

            "outer_annulus_24_sector_counts":
                coverage,

            "coverage_min":
                min(coverage),

            "position_unit_mode":
                unit_mode,

            "median_position_radius_over_chi_Mpc":
                ratio,

            "position_sky_direction_dot_median":
                float(np.median(dots)),

            "position_sky_direction_dot_min":
                float(np.min(dots)),
        },

        "tracer": {
            "target_number_density_h3_Mpc3":
                N_TARGET,

            "selected_count":
                len(selected),

            "achieved_number_density_h3_Mpc3":
                len(selected)/volume,

            "minimum_selected_halo_mass_Msun":
                minmass,

            "maximum_selected_halo_mass_Msun":
                float(np.max(mass[selected])),

            "central_only":
                True,

            "unit_number_weights":
                True,
        },

        "probes": {
            "count":
                NPROBE,

            "seed":
                PROBE_SEED,

            "inner_radius_deg":
                INNER_RADIUS_DEG,

            "redshift_true":
                [PROBE_Z0, PROBE_Z1],

            "eligible_selected_inner_count":
                len(pool),

            "boundary_membership_sha256":
                mem_sha,

            "selected_boundary_count":
                len(boundary_local),

            "selected_boundary_inside":
                boundary_included,

            "selected_boundary_excluded":
                boundary_excluded,

            "probe_galaxy_id_sha256":
                sha_bytes(
                    np.ascontiguousarray(
                        probe_gid
                    ).tobytes()
                ),

            "minimum_angular_boundary_distance_Mpc_over_h":
                float(np.min(angular_buffer)),

            "minimum_lower_radial_boundary_distance_Mpc_over_h":
                float(np.min(lower_buffer)),

            "minimum_upper_radial_boundary_distance_Mpc_over_h":
                float(np.min(upper_buffer)),
        },

        "reconstruction": {
            "smoothing_R_Mpc_over_h":
                R_SMOOTH,

            "primary_cutoff_Mpc_over_h":
                RMAX_PRIMARY,

            "check_cutoff_Mpc_over_h":
                RMAX_CHECK,

            "random_subtraction_used":
                False,

            "workers":
                1,

            "chunk":
                args.chunk,

            "same_kernel_as_E59_E62":
                True,
        },

        "gates": {
            "LOS_truth_pearson_min":
                TRUTH_R_MIN,

            "LOS_truth_sign_agreement_min":
                TRUTH_SIGN_MIN,

            "LOS_cutoff_pearson_min":
                CUTOFF_R_MIN,

            "LOS_cutoff_sign_agreement_min":
                CUTOFF_SIGN_MIN,
        },

        "metrics":
            met,

        "decision": {
            "buffered_lightcone_truth_gate_pass":
                passed,

            "interpretation":
                (
                    "OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_PASS; "
                    "NEXT_E64_SURVEY_MASK_HOD_TRUTH_TRANSFER"
                    if passed
                    else
                    "OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_FAIL; "
                    "DO_NOT_ADD_SURVEY_MASK_OR_OBSERVED_TAG"
                ),

            "important_limit":
                "E63B tests snapshot-to-lightcone transfer at the "
                "frozen E62 tracer density using central halos. "
                "It does not test eBOSS angular/radial selection, "
                "ELG HOD/photometric targeting, absolute "
                "Einstein-Vlasov amplitude, or observed significance."
        },

        "observed_eBOSS_rows_used":
            False,

        "observed_odd_used":
            False,

        "absolute_EV_amplitude_calibrated":
            False,

        "pvalue_or_detection_sigma":
            False,
    }

    atomic(OUT, result)

    L = met["LOS"]

    print()
    print("=== E63B R6 FROZEN TRUTH RESULT ===")

    print(
        "LOS",
        "R_TRUTH", L["truth_pearson"],
        "SIGN_TRUTH", L["truth_sign_agreement"],
        "R_CUTOFF", L["cutoff_192_vs_256_pearson"],
        "SIGN_CUTOFF",
        L["cutoff_192_vs_256_sign_agreement"],
        "PASS", L["gate_pass"],
        flush=True
    )

    print(
        "VECTOR_COS_MED",
        met["vector_cosine"]["median"],
        flush=True
    )

    print(
        "DECISION",
        result["decision"]["interpretation"],
        flush=True
    )

    print("OBSERVED_EBOSS_ROWS_USED", False)
    print("OBSERVED_ODD_USED", False)

    print(
        "REPORT",
        OUT,
        flush=True
    )

    print(
        "REPORT_SHA256",
        sha_file(OUT),
        flush=True
    )


if __name__ == "__main__":
    main()
