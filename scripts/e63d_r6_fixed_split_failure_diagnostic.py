#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os
import tempfile

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

RUNNER = ROOT/"scripts/e63b_r6_frozen_truth_transfer.py"

LOCK = ROOT/"eboss_workspace/cosmodc2/e63b_r6_final_pretruth_lock.json"

E63B = ROOT/"source_data/e63b_r6_outerrim_buffered_lightcone_truth_transfer_result.json"

E63C = ROOT/"source_data/e63c_r6_local_basis_failure_diagnostic.json"

OUT = ROOT/"source_data/e63d_r6_fixed_split_failure_diagnostic.json"


EXPECTED_CACHE_SHA = \
    "74d0402bfb2dc87b7ba5f174fb43adb141db4dc21c7f7b1b8acc9a20eefa1f82"

EXPECTED_MEM_SHA = \
    "112836101c35286f803f99954969d2512e55a746c8fbeae84839b486134b1dd6"


NSEL = 54893
NPROBE = 512
PROBE_SEED = 202610630001

RA0 = 62.5
DEC0 = -35.5

INNER_RADIUS_DEG = 2.5

PROBE_Z0 = 0.9
PROBE_Z1 = 1.0

# Frozen diagnostic splits.
Z_SPLIT = 0.95
DEC_SPLIT = -35.5

H = 0.71


def need(c, msg):
    if not c:
        raise RuntimeError(msg)


def sha_file(p):
    h = hashlib.sha256()

    with p.open("rb") as f:
        for b in iter(
            lambda: f.read(1024*1024),
            b""
        ):
            h.update(b)

    return h.hexdigest()


def atomic(path, obj):
    raw = (
        json.dumps(
            obj,
            indent=2,
            sort_keys=True,
            allow_nan=False
        ) + "\n"
    ).encode()

    path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with tempfile.NamedTemporaryFile(
        dir=path.parent,
        prefix=".e63d_",
        delete=False
    ) as f:
        q = Path(f.name)

        f.write(raw)
        f.flush()
        os.fsync(f.fileno())

    try:
        os.replace(q, path)

    finally:
        q.unlink(missing_ok=True)


def angular_sep_deg(ra, dec):
    r = np.deg2rad(
        np.asarray(ra, float)
    )

    d = np.deg2rad(
        np.asarray(dec, float)
    )

    r0 = math.radians(RA0)
    d0 = math.radians(DEC0)

    c = (
        np.sin(d)*math.sin(d0)
        +
        np.cos(d)*math.cos(d0)
        *
        np.cos(r-r0)
    )

    return np.rad2deg(
        np.arccos(
            np.clip(c, -1.0, 1.0)
        )
    )


def component_metrics(
    rec,
    chk,
    truth,
    basis,
    mask
):
    mask = np.asarray(mask, bool)

    need(
        np.count_nonzero(mask) >= 64,
        "Diagnostic subgroup too small"
    )

    rr = np.einsum(
        "ij,ij->i",
        rec[mask],
        basis[mask]
    )

    cc = np.einsum(
        "ij,ij->i",
        chk[mask],
        basis[mask]
    )

    tt = np.einsum(
        "ij,ij->i",
        truth[mask],
        basis[mask]
    )

    return {
        "n":
            int(np.count_nonzero(mask)),

        "truth_pearson":
            B.pearson(rr, tt),

        "truth_sign_agreement":
            B.sign_agreement(rr, tt),

        "cutoff_192_vs_256_pearson":
            B.pearson(rr, cc),

        "cutoff_192_vs_256_sign_agreement":
            B.sign_agreement(rr, cc),

        "reconstructed_mean":
            float(np.mean(rr)),

        "reconstructed_std":
            float(np.std(rr)),

        "truth_mean":
            float(np.mean(tt)),

        "truth_std":
            float(np.std(tt))
    }


def group_metrics(
    rec,
    chk,
    truth,
    bases,
    mask,
    zprobe,
    decprobe
):
    mask = np.asarray(mask, bool)

    n = int(np.count_nonzero(mask))

    need(
        n >= 64,
        f"Diagnostic subgroup too small: {n}"
    )

    return {
        "n": n,

        "z_true_mean":
            float(np.mean(zprobe[mask])),

        "z_true_min":
            float(np.min(zprobe[mask])),

        "z_true_max":
            float(np.max(zprobe[mask])),

        "dec_true_mean_deg":
            float(np.mean(decprobe[mask])),

        "dec_true_min_deg":
            float(np.min(decprobe[mask])),

        "dec_true_max_deg":
            float(np.max(decprobe[mask])),

        "local_basis": {
            name:
                component_metrics(
                    rec,
                    chk,
                    truth,
                    basis,
                    mask
                )
            for name, basis in bases.items()
        }
    }


def main():

    # --------------------------------------------------------
    # Import exact frozen E63B implementation.
    # --------------------------------------------------------

    spec = importlib.util.spec_from_file_location(
        "e63b_r6",
        RUNNER
    )

    global B
    B = importlib.util.module_from_spec(spec)

    spec.loader.exec_module(B)

    need(
        LOCK.is_file(),
        "Missing final pretruth lock"
    )

    need(
        E63B.is_file(),
        "Missing E63B result"
    )

    need(
        E63C.is_file(),
        "Missing completed E63C result"
    )

    lock = json.loads(
        LOCK.read_text()
    )

    e63b = json.loads(
        E63B.read_text()
    )

    e63c = json.loads(
        E63C.read_text()
    )

    # --------------------------------------------------------
    # Parent states must remain frozen.
    # --------------------------------------------------------

    need(
        e63b["status"]
        ==
        "FAIL_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED",
        "E63B parent status drift"
    )

    need(
        e63b["decision"][
            "buffered_lightcone_truth_gate_pass"
        ] is False,
        "E63B parent decision drift"
    )

    need(
        e63c["parent_E63B"][
            "original_fail_preserved"
        ] is True,
        "E63C provenance drift"
    )

    need(
        e63c["science_gate_reevaluated"]
        is False,
        "E63C unexpectedly reevaluated gate"
    )

    need(
        e63c["science_threshold_changed"]
        is False,
        "E63C threshold drift"
    )

    need(
        e63c["observed_eBOSS_rows_used"]
        is False
        and
        e63c["observed_odd_used"]
        is False,
        "Observation guardrail drift"
    )

    CACHE = ROOT/lock["retrieval"]["cache"]

    MEM = ROOT/lock[
        "probe_boundary"
    ]["membership"]

    need(
        sha_file(CACHE)
        ==
        EXPECTED_CACHE_SHA,
        "Frozen cache SHA drift"
    )

    need(
        sha_file(MEM)
        ==
        EXPECTED_MEM_SHA,
        "Frozen membership SHA drift"
    )

    # --------------------------------------------------------
    # Boundary membership.
    # --------------------------------------------------------

    membership = {}

    with MEM.open(newline="") as f:
        for r in csv.DictReader(f):

            gid = int(r["galaxy_id"])

            membership[gid] = bool(
                int(
                    r[
                        "server_inside_probe_z"
                    ]
                )
            )

    need(
        len(membership) == 4,
        "Membership count drift"
    )

    need(
        sum(membership.values()) == 2,
        "Membership inside-count drift"
    )

    # --------------------------------------------------------
    # Load exact frozen R6 cache.
    # --------------------------------------------------------

    with CACHE.open(newline="") as f:
        rows = list(
            csv.DictReader(f)
        )

    need(
        len(rows) == 200493,
        "R6 cache row-count drift"
    )

    gid = np.asarray([
        int(x["galaxy_id"])
        for x in rows
    ], dtype=np.int64)

    halo = np.asarray([
        int(x["halo_id"])
        for x in rows
    ], dtype=np.int64)

    ra = np.asarray([
        float(x["ra_true"])
        for x in rows
    ])

    dec = np.asarray([
        float(x["dec_true"])
        for x in rows
    ])

    z = np.asarray([
        float(x["redshift_true"])
        for x in rows
    ])

    mass = np.asarray([
        float(x["halo_mass"])
        for x in rows
    ])

    central = np.asarray([
        B.parse_bool(
            x["is_central"]
        )
        for x in rows
    ])

    need(
        np.all(central),
        "Non-central row in R6 cache"
    )

    need(
        len(np.unique(gid))
        ==
        len(gid),
        "Duplicate galaxy_id"
    )

    need(
        len(np.unique(halo))
        ==
        len(halo),
        "Duplicate halo_id"
    )

    # --------------------------------------------------------
    # Exact frozen NSEL.
    # --------------------------------------------------------

    order = np.lexsort(
        (gid, -mass)
    )

    selected = order[:NSEL]

    selected_gid = gid[selected]

    selected_sha = hashlib.sha256(
        np.ascontiguousarray(
            selected_gid
        ).tobytes()
    ).hexdigest()

    need(
        selected_sha
        ==
        e63b["source"][
            "selected_galaxy_id_sha256"
        ],
        "Selected source sample drift"
    )

    # --------------------------------------------------------
    # Exact frozen probe pool.
    # --------------------------------------------------------

    aang = angular_sep_deg(
        ra[selected],
        dec[selected]
    )

    zang = z[selected]

    inner = (
        aang <= INNER_RADIUS_DEG
    )

    pmask = (
        inner
        &
        (zang > PROBE_Z0)
        &
        (zang < PROBE_Z1)
    )

    boundary = (
        inner
        &
        (
            (zang == PROBE_Z0)
            |
            (zang == PROBE_Z1)
        )
    )

    for ii in np.flatnonzero(
        boundary
    ):
        g = int(
            selected_gid[ii]
        )

        need(
            g in membership,
            f"Missing boundary membership {g}"
        )

        pmask[ii] = membership[g]

    pool = np.flatnonzero(
        pmask
    )

    need(
        len(pool)
        ==
        e63b["probes"][
            "eligible_selected_inner_count"
        ],
        "Eligible probe pool drift"
    )

    rng = np.random.default_rng(
        PROBE_SEED
    )

    pidx = np.sort(
        rng.choice(
            pool,
            size=NPROBE,
            replace=False
        )
    )

    probe_gid = selected_gid[
        pidx
    ]

    probe_sha = hashlib.sha256(
        np.ascontiguousarray(
            probe_gid
        ).tobytes()
    ).hexdigest()

    need(
        probe_sha
        ==
        e63b["probes"][
            "probe_galaxy_id_sha256"
        ],
        "Frozen probe sample drift"
    )

    print(
        "E63D_FROZEN_SAMPLE_REPRODUCED",
        "NSEL", len(selected),
        "POOL", len(pool),
        "PROBES", len(pidx),
        flush=True
    )

    # --------------------------------------------------------
    # Position-unit convention exactly as parents.
    # --------------------------------------------------------

    take = np.linspace(
        0,
        len(rows)-1,
        1024,
        dtype=int
    )

    pos_guard = np.column_stack([
        [
            float(
                rows[i]["position_x"]
            )
            for i in take
        ],
        [
            float(
                rows[i]["position_y"]
            )
            for i in take
        ],
        [
            float(
                rows[i]["position_z"]
            )
            for i in take
        ]
    ])

    ch = np.asarray([
        B.chi_mpc(z[i])
        for i in take
    ])

    ratio = float(
        np.median(
            np.linalg.norm(
                pos_guard,
                axis=1
            )
            /
            ch
        )
    )

    if 0.97 <= ratio <= 1.03:

        scale = H

        unit_mode = \
            "IRSA_GCR_Mpc_to_Mpc_over_h"

    elif (
        H-0.025
        <= ratio
        <= H+0.025
    ):

        scale = 1.0

        unit_mode = \
            "already_Mpc_over_h"

    else:
        raise RuntimeError(
            f"Unknown position units {ratio}"
        )

    need(
        unit_mode
        ==
        e63b["geometry"][
            "position_unit_mode"
        ],
        "Position unit mode drift"
    )

    # --------------------------------------------------------
    # Exact same frozen positions and truth.
    # --------------------------------------------------------

    pos_all = np.column_stack([
        np.asarray([
            float(x["position_x"])
            for x in rows
        ]),
        np.asarray([
            float(x["position_y"])
            for x in rows
        ]),
        np.asarray([
            float(x["position_z"])
            for x in rows
        ])
    ]) * scale

    vel_all = np.column_stack([
        np.asarray([
            float(x["velocity_x"])
            for x in rows
        ]),
        np.asarray([
            float(x["velocity_y"])
            for x in rows
        ]),
        np.asarray([
            float(x["velocity_z"])
            for x in rows
        ])
    ])

    source_pos = pos_all[
        selected
    ]

    ppos = source_pos[
        pidx
    ]

    ptruth = vel_all[
        selected
    ][pidx]

    probe_z = z[
        selected
    ][pidx]

    probe_dec = dec[
        selected
    ][pidx]

    # --------------------------------------------------------
    # Same frozen E63B reconstruction.
    # --------------------------------------------------------

    print(
        "E63D_RECONSTRUCT_EXACT_FROZEN_E63B",
        "SOURCES", len(source_pos),
        "PROBES", len(ppos),
        flush=True
    )

    rec, chk = B.reconstruct(
        source_pos,
        ppos,
        chunk=1
    )

    # --------------------------------------------------------
    # Reproduce parent LOS exactly before any split.
    # --------------------------------------------------------

    parent = B.metrics(
        rec,
        chk,
        ptruth,
        ppos
    )

    old = e63b["metrics"]["LOS"]

    for k in [
        "truth_pearson",
        "truth_sign_agreement",
        "cutoff_192_vs_256_pearson",
        "cutoff_192_vs_256_sign_agreement"
    ]:
        need(
            abs(
                float(parent["LOS"][k])
                -
                float(old[k])
            )
            < 1e-12,
            f"Parent LOS mismatch {k}"
        )

    print(
        "E63D_PARENT_LOS_EXACT_REPRODUCTION_PASS",
        "R_TRUTH",
        parent["LOS"]["truth_pearson"],
        "SIGN_TRUTH",
        parent["LOS"][
            "truth_sign_agreement"
        ],
        flush=True
    )

    # --------------------------------------------------------
    # Same corrected local basis as completed E63C.
    # --------------------------------------------------------

    e_r = (
        ppos
        /
        np.linalg.norm(
            ppos,
            axis=1
        )[:, None]
    )

    phi = np.arctan2(
        e_r[:, 1],
        e_r[:, 0]
    )

    e_RA = np.column_stack([
        -np.sin(phi),
        np.cos(phi),
        np.zeros_like(phi)
    ])

    e_Dec = np.cross(
        e_r,
        e_RA
    )

    e_Dec /= np.linalg.norm(
        e_Dec,
        axis=1
    )[:, None]

    bases = {
        "radial_e_r": e_r,
        "tangential_e_RA": e_RA,
        "tangential_e_Dec": e_Dec
    }

    # --------------------------------------------------------
    # Full-sample E63C reproduction.
    # --------------------------------------------------------

    full_mask = np.ones(
        NPROBE,
        dtype=bool
    )

    full_local = {
        name:
            component_metrics(
                rec,
                chk,
                ptruth,
                basis,
                full_mask
            )
        for name, basis
        in bases.items()
    }

    for name in bases:
        oldc = e63c[
            "local_basis_metrics"
        ][name]

        for k in [
            "truth_pearson",
            "truth_sign_agreement",
            "cutoff_192_vs_256_pearson",
            "cutoff_192_vs_256_sign_agreement"
        ]:
            need(
                abs(
                    float(full_local[name][k])
                    -
                    float(oldc[k])
                )
                < 1e-12,
                f"E63C reproduction mismatch "
                f"{name} {k}"
            )

    print(
        "E63D_PARENT_E63C_LOCAL_BASIS_EXACT_REPRODUCTION_PASS",
        flush=True
    )

    # --------------------------------------------------------
    # PREDECLARED FIXED SPLITS.
    #
    # No optimization and no search for a favorable threshold.
    #
    # Redshift split is explicitly the serialized cache value.
    # Exact 0.9500 entries are assigned to HIGH by definition.
    # --------------------------------------------------------

    z_low = (
        probe_z < Z_SPLIT
    )

    z_high = (
        probe_z >= Z_SPLIT
    )

    dec_south = (
        probe_dec < DEC_SPLIT
    )

    dec_north = (
        probe_dec >= DEC_SPLIT
    )

    split_masks = {
        "z_low_lt_0p95":
            z_low,

        "z_high_ge_0p95":
            z_high,

        "dec_south_lt_minus35p5":
            dec_south,

        "dec_north_ge_minus35p5":
            dec_north
    }

    for name, mask in split_masks.items():
        need(
            np.count_nonzero(mask)
            >= 64,
            f"Too few probes in {name}: "
            f"{np.count_nonzero(mask)}"
        )

    exact_z_split_count = int(
        np.count_nonzero(
            probe_z == Z_SPLIT
        )
    )

    groups = {
        name:
            group_metrics(
                rec,
                chk,
                ptruth,
                bases,
                mask,
                probe_z,
                probe_dec
            )
        for name, mask
        in split_masks.items()
    }

    # --------------------------------------------------------
    # Simple descriptive contrasts.
    # No hypothesis-test p-values and no new gate.
    # --------------------------------------------------------

    def contrast(a, b):
        out = {}

        for comp in bases:
            A = groups[a][
                "local_basis"
            ][comp]

            Bv = groups[b][
                "local_basis"
            ][comp]

            out[comp] = {
                "delta_truth_pearson_B_minus_A":
                    float(
                        Bv["truth_pearson"]
                        -
                        A["truth_pearson"]
                    ),

                "delta_truth_sign_B_minus_A":
                    float(
                        Bv[
                            "truth_sign_agreement"
                        ]
                        -
                        A[
                            "truth_sign_agreement"
                        ]
                    ),

                "delta_cutoff_pearson_B_minus_A":
                    float(
                        Bv[
                            "cutoff_192_vs_256_pearson"
                        ]
                        -
                        A[
                            "cutoff_192_vs_256_pearson"
                        ]
                    )
            }

        return out

    contrasts = {
        "redshift_high_minus_low":
            contrast(
                "z_low_lt_0p95",
                "z_high_ge_0p95"
            ),

        "north_minus_south":
            contrast(
                "dec_south_lt_minus35p5",
                "dec_north_ge_minus35p5"
            )
    }

    result = {
        "stage":
            "E63D_R6_FIXED_SPLIT_POST_FAILURE_DIAGNOSTIC",

        "purpose":
            "Describe whether the frozen E63B lightcone-transfer "
            "failure varies across preregistered redshift and "
            "declination halves without changing estimator, "
            "sample, thresholds or science decision.",

        "parents": {
            "E63B_status":
                e63b["status"],

            "E63B_original_fail_preserved":
                True,

            "E63C_result":
                str(
                    E63C.relative_to(ROOT)
                ),

            "E63C_sha256":
                sha_file(E63C)
        },

        "frozen_inputs": {
            "cache_sha256":
                EXPECTED_CACHE_SHA,

            "membership_sha256":
                EXPECTED_MEM_SHA,

            "selected_galaxy_id_sha256":
                selected_sha,

            "probe_galaxy_id_sha256":
                probe_sha,

            "NSEL":
                NSEL,

            "NPROBE":
                NPROBE,

            "probe_seed":
                PROBE_SEED
        },

        "fixed_splits": {
            "redshift": {
                "variable":
                    "serialized redshift_true in frozen R6 cache",

                "threshold":
                    Z_SPLIT,

                "low":
                    "redshift_true < 0.95",

                "high":
                    "redshift_true >= 0.95",

                "exact_threshold_probe_count":
                    exact_z_split_count
            },

            "declination": {
                "variable":
                    "dec_true_deg",

                "threshold_deg":
                    DEC_SPLIT,

                "south":
                    "dec_true < -35.5 deg",

                "north":
                    "dec_true >= -35.5 deg"
            }
        },

        "parent_full_sample_local_basis":
            full_local,

        "groups":
            groups,

        "descriptive_contrasts":
            contrasts,

        "diagnostic_only":
            True,

        "pvalues_computed":
            False,

        "new_science_gate_defined":
            False,

        "science_gate_reevaluated":
            False,

        "science_threshold_changed":
            False,

        "estimator_changed":
            False,

        "probe_sample_changed":
            False,

        "observed_eBOSS_rows_used":
            False,

        "observed_odd_used":
            False,

        "E63B_fail_preserved":
            True
    }

    atomic(
        OUT,
        result
    )

    # --------------------------------------------------------
    # Human-readable output.
    # --------------------------------------------------------

    print()
    print(
        "=== E63D FIXED SPLIT COUNTS ==="
    )

    for name, mask in split_masks.items():
        print(
            name,
            "N",
            int(np.count_nonzero(mask))
        )

    print(
        "EXACT_Z_0.9500_PROBES",
        exact_z_split_count
    )

    for split_name, names in [
        (
            "REDSHIFT",
            [
                "z_low_lt_0p95",
                "z_high_ge_0p95"
            ]
        ),
        (
            "DECLINATION",
            [
                "dec_south_lt_minus35p5",
                "dec_north_ge_minus35p5"
            ]
        )
    ]:

        print()
        print(
            f"=== E63D {split_name} SPLIT ==="
        )

        for name in names:

            g = groups[name]

            print()
            print(
                name,
                "N", g["n"],
                "ZMEAN", g["z_true_mean"],
                "DECMEAN", g[
                    "dec_true_mean_deg"
                ]
            )

            for comp, v in g[
                "local_basis"
            ].items():

                print(
                    comp,
                    "R_TRUTH",
                    v["truth_pearson"],
                    "SIGN_TRUTH",
                    v[
                        "truth_sign_agreement"
                    ],
                    "R_CUTOFF",
                    v[
                        "cutoff_192_vs_256_pearson"
                    ],
                    "SIGN_CUTOFF",
                    v[
                        "cutoff_192_vs_256_sign_agreement"
                    ]
                )

    print()
    print(
        "=== E63D DESCRIPTIVE CONTRASTS ==="
    )

    for cname, c in contrasts.items():

        print()
        print(cname)

        for comp, v in c.items():

            print(
                comp,
                "DELTA_R_TRUTH",
                v[
                    "delta_truth_pearson_B_minus_A"
                ],
                "DELTA_SIGN_TRUTH",
                v[
                    "delta_truth_sign_B_minus_A"
                ],
                "DELTA_R_CUTOFF",
                v[
                    "delta_cutoff_pearson_B_minus_A"
                ]
            )

    print()
    print(
        "REPORT",
        OUT
    )

    print(
        "REPORT_SHA256",
        sha_file(OUT)
    )

    print(
        "E63D_R6_FIXED_SPLIT_DIAGNOSTIC_COMPLETE"
    )

    print(
        "ORIGINAL_E63B_FAIL_PRESERVED"
    )

    print(
        "NO_NEW_SCIENCE_GATE"
    )

    print(
        "NO_PVALUES"
    )

    print(
        "OBSERVED_EBOSS_REMAINS_SEALED"
    )


if __name__ == "__main__":
    main()
