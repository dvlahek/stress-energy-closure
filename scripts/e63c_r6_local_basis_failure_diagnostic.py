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


ROOT=Path(__file__).resolve().parents[1]

RUNNER=ROOT/"scripts/e63b_r6_frozen_truth_transfer.py"

LOCK=ROOT/"eboss_workspace/cosmodc2/e63b_r6_final_pretruth_lock.json"

E63B=ROOT/"source_data/e63b_r6_outerrim_buffered_lightcone_truth_transfer_result.json"

OUT=ROOT/"source_data/e63c_r6_local_basis_failure_diagnostic.json"

EXPECTED_CACHE_SHA="74d0402bfb2dc87b7ba5f174fb43adb141db4dc21c7f7b1b8acc9a20eefa1f82"
EXPECTED_MEM_SHA="112836101c35286f803f99954969d2512e55a746c8fbeae84839b486134b1dd6"

NSEL=54893
NPROBE=512
PROBE_SEED=202610630001

RA0=62.5
DEC0=-35.5

INNER_RADIUS_DEG=2.5

PROBE_Z0=0.9
PROBE_Z1=1.0

H=0.71


def need(c,msg):
    if not c:
        raise RuntimeError(msg)


def sha_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()


def atomic(path,obj):
    raw=(json.dumps(
        obj,
        indent=2,
        sort_keys=True,
        allow_nan=False
    )+"\n").encode()

    path.parent.mkdir(parents=True,exist_ok=True)

    with tempfile.NamedTemporaryFile(
        dir=path.parent,
        prefix=".e63c_",
        delete=False
    ) as f:
        q=Path(f.name)
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())

    try:
        os.replace(q,path)
    finally:
        q.unlink(missing_ok=True)


def parse_bool(v):
    s=str(v).strip().lower()

    if s in ("true","1","t","yes"):
        return True

    if s in ("false","0","f","no"):
        return False

    raise RuntimeError(f"Unexpected bool {v!r}")


def angular_sep_deg(ra,dec):
    r=np.deg2rad(np.asarray(ra,float))
    d=np.deg2rad(np.asarray(dec,float))

    r0=math.radians(RA0)
    d0=math.radians(DEC0)

    c=(
        np.sin(d)*math.sin(d0)
        + np.cos(d)*math.cos(d0)*np.cos(r-r0)
    )

    return np.rad2deg(
        np.arccos(np.clip(c,-1.0,1.0))
    )


def pearson(a,b):
    a=np.asarray(a,float)
    b=np.asarray(b,float)

    aa=a-np.mean(a)
    bb=b-np.mean(b)

    den=float(
        np.linalg.norm(aa)
        * np.linalg.norm(bb)
    )

    return float(aa@bb/den) if den>0 else 0.0


def sign_agreement(a,b):
    a=np.asarray(a,float)
    b=np.asarray(b,float)

    nz=(a!=0)&(b!=0)

    need(
        np.count_nonzero(nz)>=0.95*len(a),
        "Too many zero components"
    )

    return float(
        np.mean(
            np.sign(a[nz])==np.sign(b[nz])
        )
    )


def component_metrics(rec,chk,truth,basis):
    rr=np.einsum("ij,ij->i",rec,basis)
    cc=np.einsum("ij,ij->i",chk,basis)
    tt=np.einsum("ij,ij->i",truth,basis)

    return {
        "truth_pearson":
            pearson(rr,tt),

        "truth_sign_agreement":
            sign_agreement(rr,tt),

        "cutoff_192_vs_256_pearson":
            pearson(rr,cc),

        "cutoff_192_vs_256_sign_agreement":
            sign_agreement(rr,cc),

        "reconstructed_std":
            float(np.std(rr)),

        "truth_std":
            float(np.std(tt)),

        "reconstructed_mean":
            float(np.mean(rr)),

        "truth_mean":
            float(np.mean(tt))
    }


def main():

    # --------------------------------------------------------
    # Import exact frozen E63B implementation.
    # --------------------------------------------------------

    spec=importlib.util.spec_from_file_location(
        "e63b_r6",
        RUNNER
    )

    B=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(B)

    need(LOCK.is_file(),"Missing final pretruth lock")
    need(E63B.is_file(),"Missing frozen E63B result")

    lock=json.loads(LOCK.read_text())
    e63b=json.loads(E63B.read_text())

    need(
        e63b["status"]
        ==
        "FAIL_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED",
        "E63C requires frozen E63B FAIL parent"
    )

    need(
        e63b["decision"]["buffered_lightcone_truth_gate_pass"]
        is False,
        "Frozen E63B decision drift"
    )

    need(
        e63b["observed_eBOSS_rows_used"] is False
        and
        e63b["observed_odd_used"] is False,
        "Observed guardrail drift"
    )

    CACHE=ROOT/lock["retrieval"]["cache"]
    MEM=ROOT/lock["probe_boundary"]["membership"]

    need(sha_file(CACHE)==EXPECTED_CACHE_SHA,
         "R6 cache SHA drift")

    need(sha_file(MEM)==EXPECTED_MEM_SHA,
         "boundary membership SHA drift")

    membership={}

    with MEM.open(newline="") as f:
        for r in csv.DictReader(f):
            membership[int(r["galaxy_id"])]=bool(
                int(r["server_inside_probe_z"])
            )

    need(len(membership)==4,
         "membership count drift")

    # --------------------------------------------------------
    # Exact frozen cache.
    # --------------------------------------------------------

    with CACHE.open(newline="") as f:
        rows=list(csv.DictReader(f))

    need(len(rows)==200493,
         "R6 candidate count drift")

    gid=np.asarray([
        int(x["galaxy_id"])
        for x in rows
    ],dtype=np.int64)

    halo=np.asarray([
        int(x["halo_id"])
        for x in rows
    ],dtype=np.int64)

    ra=np.asarray([
        float(x["ra_true"])
        for x in rows
    ],dtype=float)

    dec=np.asarray([
        float(x["dec_true"])
        for x in rows
    ],dtype=float)

    z=np.asarray([
        float(x["redshift_true"])
        for x in rows
    ],dtype=float)

    mass=np.asarray([
        float(x["halo_mass"])
        for x in rows
    ],dtype=float)

    central=np.asarray([
        parse_bool(x["is_central"])
        for x in rows
    ])

    need(np.all(central),
         "noncentral row")

    need(len(np.unique(gid))==len(gid),
         "duplicate galaxy_id")

    need(len(np.unique(halo))==len(halo),
         "duplicate halo_id")

    # --------------------------------------------------------
    # Reproduce exact NSEL.
    # --------------------------------------------------------

    order=np.lexsort((gid,-mass))
    selected=order[:NSEL]

    selected_gid=gid[selected]

    selected_sha=hashlib.sha256(
        np.ascontiguousarray(
            selected_gid
        ).tobytes()
    ).hexdigest()

    need(
        selected_sha
        ==
        e63b["source"]["selected_galaxy_id_sha256"],
        "Frozen selected sample SHA mismatch"
    )

    aang=angular_sep_deg(
        ra[selected],
        dec[selected]
    )

    zang=z[selected]

    inner=aang<=INNER_RADIUS_DEG

    # Strict interior first.
    pmask=(
        inner
        &
        (zang>PROBE_Z0)
        &
        (zang<PROBE_Z1)
    )

    # Serialized boundaries resolved only through frozen membership.
    boundary=(
        inner
        &
        (
            (zang==PROBE_Z0)
            |
            (zang==PROBE_Z1)
        )
    )

    for ii in np.flatnonzero(boundary):
        g=int(selected_gid[ii])

        need(
            g in membership,
            f"Missing boundary membership for {g}"
        )

        pmask[ii]=membership[g]

    pool=np.flatnonzero(pmask)

    need(
        len(pool)
        ==
        e63b["probes"]["eligible_selected_inner_count"],
        "Eligible probe pool drift"
    )

    rng=np.random.default_rng(PROBE_SEED)

    pidx=np.sort(
        rng.choice(
            pool,
            size=NPROBE,
            replace=False
        )
    )

    probe_gid=selected_gid[pidx]

    probe_sha=hashlib.sha256(
        np.ascontiguousarray(
            probe_gid
        ).tobytes()
    ).hexdigest()

    need(
        probe_sha
        ==
        e63b["probes"]["probe_galaxy_id_sha256"],
        "Frozen 512-probe SHA mismatch"
    )

    print(
        "E63C_FROZEN_SAMPLE_REPRODUCED",
        "NSEL",len(selected),
        "POOL",len(pool),
        "PROBES",len(pidx),
        flush=True
    )

    # --------------------------------------------------------
    # Position unit mode must reproduce frozen E63B.
    # --------------------------------------------------------

    take=np.linspace(
        0,
        len(rows)-1,
        1024,
        dtype=int
    )

    pos_guard=np.column_stack([
        [float(rows[i]["position_x"]) for i in take],
        [float(rows[i]["position_y"]) for i in take],
        [float(rows[i]["position_z"]) for i in take]
    ])

    ch=np.asarray([
        B.chi_mpc(z[i])
        for i in take
    ])

    ratio=float(
        np.median(
            np.linalg.norm(pos_guard,axis=1)
            /
            ch
        )
    )

    if 0.97<=ratio<=1.03:
        scale=H
        unit_mode="IRSA_GCR_Mpc_to_Mpc_over_h"

    elif H-0.025<=ratio<=H+0.025:
        scale=1.0
        unit_mode="already_Mpc_over_h"

    else:
        raise RuntimeError(
            f"Position unit ratio {ratio}"
        )

    need(
        unit_mode
        ==
        e63b["geometry"]["position_unit_mode"],
        "Position unit mode drift"
    )

    # --------------------------------------------------------
    # Same frozen positions and truth vectors.
    # --------------------------------------------------------

    pos_all=np.column_stack([
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
    ])*scale

    vel_all=np.column_stack([
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

    source_pos=pos_all[selected]

    ppos=source_pos[pidx]
    ptruth=vel_all[selected][pidx]

    # --------------------------------------------------------
    # Exact same E63B reconstruction.
    # No new estimator, no retuning.
    # --------------------------------------------------------

    print(
        "E63C_RECONSTRUCT_EXACT_FROZEN_E63B",
        "SOURCES",len(source_pos),
        "PROBES",len(ppos),
        flush=True
    )

    rec,chk=B.reconstruct(
        source_pos,
        ppos,
        chunk=1
    )

    # --------------------------------------------------------
    # First consistency check: reproduce original LOS metrics.
    # --------------------------------------------------------

    base=B.metrics(
        rec,
        chk,
        ptruth,
        ppos
    )

    old=e63b["metrics"]["LOS"]
    now=base["LOS"]

    for k in [
        "truth_pearson",
        "truth_sign_agreement",
        "cutoff_192_vs_256_pearson",
        "cutoff_192_vs_256_sign_agreement"
    ]:
        need(
            abs(float(old[k])-float(now[k]))<1e-12,
            f"LOS reproduction mismatch {k}: "
            f"{old[k]} vs {now[k]}"
        )

    print(
        "E63C_PARENT_LOS_EXACT_REPRODUCTION_PASS",
        "R_TRUTH",now["truth_pearson"],
        "SIGN_TRUTH",now["truth_sign_agreement"],
        flush=True
    )

    # --------------------------------------------------------
    # Local spherical basis.
    #
    # e_r   : observer -> galaxy
    # e_RA  : increasing right ascension / east
    # e_Dec : increasing declination / north
    # --------------------------------------------------------

    # IMPORTANT:
    # Frozen E63B defines LOS directly from the physical position vector
    # ppos/|ppos|.  Use that exact direction here so the radial component
    # is mathematically identical to the parent E63B LOS observable.
    #
    # RA/Dec and the GCR position vector agree to extremely high accuracy,
    # but they are not bit-identical.  Building e_r independently from
    # serialized RA/Dec therefore introduced an artificial diagnostic
    # mismatch under the deliberately strict 1e-12 reproduction guard.

    e_r=ppos/np.linalg.norm(ppos,axis=1)[:,None]

    # Local east / increasing-RA direction obtained from the actual
    # position-vector azimuth.  This is perpendicular to e_r by
    # construction.
    phi=np.arctan2(e_r[:,1],e_r[:,0])

    e_RA=np.column_stack([
        -np.sin(phi),
        np.cos(phi),
        np.zeros_like(phi)
    ])

    # Local north / increasing-declination direction.  Cross-product
    # construction makes the three vectors an orthonormal right-handed
    # basis tied to the same physical position used by frozen E63B.
    e_Dec=np.cross(e_r,e_RA)

    e_Dec/=np.linalg.norm(
        e_Dec,
        axis=1
    )[:,None]

    # Orthonormality guard.
    need(
        np.max(
            np.abs(
                np.einsum("ij,ij->i",e_r,e_RA)
            )
        )<1e-12,
        "e_r/e_RA not orthogonal"
    )

    need(
        np.max(
            np.abs(
                np.einsum("ij,ij->i",e_r,e_Dec)
            )
        )<1e-12,
        "e_r/e_Dec not orthogonal"
    )

    need(
        np.max(
            np.abs(
                np.einsum("ij,ij->i",e_RA,e_Dec)
            )
        )<1e-12,
        "e_RA/e_Dec not orthogonal"
    )

    local={
        "radial_e_r":
            component_metrics(
                rec,chk,ptruth,e_r
            ),

        "tangential_e_RA":
            component_metrics(
                rec,chk,ptruth,e_RA
            ),

        "tangential_e_Dec":
            component_metrics(
                rec,chk,ptruth,e_Dec
            )
    }

    # Radial must be identical to E63B LOS.
    need(
        abs(
            local["radial_e_r"]["truth_pearson"]
            -
            old["truth_pearson"]
        )<1e-12,
        "Local radial != frozen LOS truth Pearson"
    )

    need(
        abs(
            local["radial_e_r"]["truth_sign_agreement"]
            -
            old["truth_sign_agreement"]
        )<1e-12,
        "Local radial != frozen LOS sign"
    )

    # --------------------------------------------------------
    # Additional diagnostic:
    # z axis decomposition relative to local basis.
    #
    # This tells us geometrically why global z can behave
    # differently without changing any science metric.
    # --------------------------------------------------------

    zhat=np.array([0.0,0.0,1.0])

    z_basis_geometry={
        "mean_abs_projection_on_e_r":
            float(np.mean(np.abs(e_r@zhat))),

        "mean_abs_projection_on_e_RA":
            float(np.mean(np.abs(e_RA@zhat))),

        "mean_abs_projection_on_e_Dec":
            float(np.mean(np.abs(e_Dec@zhat))),

        "median_projection_on_e_r":
            float(np.median(e_r@zhat)),

        "median_projection_on_e_Dec":
            float(np.median(e_Dec@zhat))
    }

    result={
        "stage":
            "E63C_R6_LOCAL_BASIS_POST_FAILURE_DIAGNOSTIC",

        "parent_E63B":{
            "result":
                str(E63B.relative_to(ROOT)),

            "status":
                e63b["status"],

            "gate_pass":
                False,

            "original_fail_preserved":
                True
        },

        "frozen_inputs":{
            "cache_sha256":
                EXPECTED_CACHE_SHA,

            "boundary_membership_sha256":
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

        "parent_LOS_reproduction":
            now,

        "local_basis_metrics":
            local,

        "global_cartesian_parent_metrics":
            e63b["metrics"]["cartesian"],

        "z_axis_local_basis_geometry":
            z_basis_geometry,

        "diagnostic_only":
            True,

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
            False
    }

    atomic(OUT,result)

    print()
    print("=== E63C LOCAL BASIS DIAGNOSTIC ===")

    for k,v in local.items():
        print(
            k,
            "R_TRUTH",v["truth_pearson"],
            "SIGN_TRUTH",v["truth_sign_agreement"],
            "R_CUTOFF",v["cutoff_192_vs_256_pearson"],
            "SIGN_CUTOFF",
            v["cutoff_192_vs_256_sign_agreement"]
        )

    print()
    print("=== GLOBAL z GEOMETRY ===")

    for k,v in z_basis_geometry.items():
        print(k,v)

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
        "E63C_R6_LOCAL_BASIS_DIAGNOSTIC_COMPLETE"
    )

    print(
        "ORIGINAL_E63B_FAIL_PRESERVED"
    )

    print(
        "NO_SCIENCE_GATE_REEVALUATION"
    )

    print(
        "OBSERVED_EBOSS_REMAINS_SEALED"
    )


if __name__=="__main__":
    main()
