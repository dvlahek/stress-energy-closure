#!/usr/bin/env python3

from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os
import shutil
import time

ROOT=Path(__file__).resolve().parents[1]

R5_DIR=ROOT/"eboss_workspace/cosmodc2/e63b_r5_tiles"
NEW_DIR=ROOT/"eboss_workspace/cosmodc2/e63b_r6_hybrid_missing"
OUT=ROOT/"eboss_workspace/cosmodc2/e63b_outerrim_central_candidates_R6.csv"
META=ROOT/"eboss_workspace/cosmodc2/e63b_r6_hybrid_retrieval_meta.json"

NEW_DIR.mkdir(parents=True,exist_ok=True)

TABLE="cosmodc2mockv1"

RA0=62.5
DEC0=-35.5
OUTER_RADIUS_DEG=9.5

SOURCE_Z0=0.75
SOURCE_Z1=1.16
MASS_FLOOR=5.0e12
TOP_CAP=600000

REQ=[
 "galaxy_id","halo_id","ra_true","dec_true","redshift_true",
 "position_x","position_y","position_z",
 "velocity_x","velocity_y","velocity_z",
 "halo_mass","is_central"
]

helper=ROOT/"scripts/e63_cosmodc2_outerrim_lightcone_source_preflight.py"
spec=importlib.util.spec_from_file_location("e63tap",helper)
P=importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)

def need(c,msg):
    if not c:
        raise RuntimeError(msg)

def sha_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def parse_bool(v):
    s=str(v).strip().lower()
    if s in ("true","1","t","yes"):
        return True
    if s in ("false","0","f","no"):
        return False
    raise RuntimeError(f"Unexpected boolean {v!r}")

def angular_sep_deg(ra,dec):
    r=math.radians(float(ra))
    d=math.radians(float(dec))
    r0=math.radians(RA0)
    d0=math.radians(DEC0)

    c=(
        math.sin(d)*math.sin(d0)
        +math.cos(d)*math.cos(d0)*math.cos(r-r0)
    )

    return math.degrees(
        math.acos(max(-1.0,min(1.0,c)))
    )

def position_angle_deg(ra,dec):
    r=math.radians(float(ra))
    d=math.radians(float(dec))
    r0=math.radians(RA0)
    d0=math.radians(DEC0)

    dr=r-r0
    y=math.sin(dr)*math.cos(d)
    x=(
        math.cos(d0)*math.sin(d)
        -math.sin(d0)*math.cos(d)*math.cos(dr)
    )

    return (
        math.degrees(math.atan2(y,x))+360.0
    )%360.0

# ------------------------------------------------------------
# Frozen missing geometry.
#
# Existing R5:
#   RA 42--68
#   Dec -50.6--31.4
#
# R6:
#   RA 50.5--74.5
#   Dec -45.3--25.7
#
# Missing union is represented by two NON-OVERLAPPING strips:
#
# East:
#   RA [68,74.5], Dec [-45.3,-25.7]
#   split into 4 Dec boxes.
#
# North:
#   RA [50.5,68], Dec [-31.4,-25.7]
#   split into 3 RA boxes.
# ------------------------------------------------------------

boxes=[]

# East strip: width=6.5, four boxes height=4.9.
for j in range(4):
    d0=-45.3 + j*4.9
    d1=d0+4.9

    boxes.append({
        "name":f"east_{j}",
        "ra0":68.0,
        "ra1":74.5,
        "dec0":d0,
        "dec1":d1,
    })

# North strip: height=5.7, three boxes.
north_width=(68.0-50.5)/3.0

for j in range(3):
    r0=50.5+j*north_width
    r1=50.5+(j+1)*north_width

    boxes.append({
        "name":f"north_{j}",
        "ra0":r0,
        "ra1":r1,
        "dec0":-31.4,
        "dec1":-25.7,
    })

queries=[]

print("E63B_R6_HYBRID_NEW_BOX_COUNT",len(boxes),flush=True)

# ------------------------------------------------------------
# Retrieve only the missing 7 boxes.
# ------------------------------------------------------------

for k,b in enumerate(boxes):

    rc=0.5*(b["ra0"]+b["ra1"])
    dc=0.5*(b["dec0"]+b["dec1"])
    width=b["ra1"]-b["ra0"]
    height=b["dec1"]-b["dec0"]

    cols=",".join(REQ)

    q=(
      f"SELECT TOP {TOP_CAP} {cols} "
      f"FROM {TABLE} WHERE "
      f"CONTAINS("
      f"POINT('ICRS',ra_true,dec_true),"
      f"BOX('ICRS',{rc:.10f},{dc:.10f},{width:.10f},{height:.10f})"
      f")=1 AND "
      f"(redshift_true>={SOURCE_Z0} "
      f"AND redshift_true<{SOURCE_Z1} "
      f"AND halo_mass>={MASS_FLOOR:.1f})"
    )

    qsha=hashlib.sha256(q.encode()).hexdigest()

    cp=NEW_DIR/f"{k:02d}_{b['name']}.csv"
    sp=NEW_DIR/f"{k:02d}_{b['name']}_tap_job.json"

    # Preserve a terminal infrastructure ERROR and allow an identical
    # frozen query to be resubmitted.
    if sp.exists() and not cp.exists():
        st=json.loads(sp.read_text())

        if (
            st.get("adql_sha256")
            and st["adql_sha256"] != qsha
        ):
            raise RuntimeError(
                f"ADQL SHA drift for {b['name']}"
            )

        if str(st.get("phase","")).upper()=="ERROR":
            jid=str(st.get("job_url","")).rstrip("/").split("/")[-1]
            archive=sp.with_name(
                sp.stem+
                f".ERROR_PRESERVED_job{jid}_"+
                time.strftime("%Y%m%d_%H%M%S")+
                ".json"
            )

            shutil.copy2(sp,archive)
            sp.unlink()

            print(
                "E63B_R6_HYBRID_ERROR_ARCHIVED",
                b["name"],archive,
                flush=True
            )

    if cp.exists():
        print(
            "E63B_R6_HYBRID_REUSE",
            k,b["name"],cp,
            flush=True
        )

    else:
        print(
            "E63B_R6_HYBRID_START",
            k,b["name"],
            "RA",rc,
            "DEC",dc,
            "W",width,
            "H",height,
            flush=True
        )

        raw=P.query_async(
            q,
            poll_seconds=5,
            max_wait_seconds=86400,
            state_path=sp,
        )

        tmp=cp.with_suffix(".tmp")
        tmp.write_text(raw)
        os.replace(tmp,cp)

        print(
            "E63B_R6_HYBRID_CACHED",
            k,b["name"],cp,
            flush=True
        )

    # Stream validation.
    n=0
    zupper=0

    with cp.open(newline="") as f:
        r=csv.DictReader(f)

        need(
            set(REQ)<=set(r.fieldnames or []),
            f"Bad columns in {cp}"
        )

        for row in r:
            n+=1

            z=float(row["redshift_true"])
            mass=float(row["halo_mass"])

            # Server query is strictly z < 1.16, but CSV serialization
            # may round valid internal values to 1.1600.
            need(
                z>=SOURCE_Z0 and z<=SOURCE_Z1,
                f"{b['name']} z leakage {z}"
            )

            if z==SOURCE_Z1:
                zupper+=1

            need(
                mass>=MASS_FLOOR,
                f"{b['name']} mass leakage {mass}"
            )

    need(
        n<TOP_CAP,
        f"{b['name']} hit TOP {TOP_CAP}; fail closed and subdivide"
    )

    print(
        "E63B_R6_HYBRID_ROWS",
        k,b["name"],n,
        flush=True
    )

    if zupper:
        print(
            "E63B_R6_HYBRID_SERIALIZED_UPPER_BOUNDARY",
            k,b["name"],zupper,
            "CSV_Z",SOURCE_Z1,
            "SERVER_CONTRACT redshift_true<1.16",
            flush=True
        )

    queries.append({
        **b,
        "center_ra":rc,
        "center_dec":dc,
        "width":width,
        "height":height,
        "rows":n,
        "query":q,
        "adql_sha256":qsha,
        "cache":str(cp.relative_to(ROOT)),
        "state":str(sp.relative_to(ROOT)),
    })

# ------------------------------------------------------------
# Merge:
#   all 16 existing R5 raw tiles
#   + only the 7 missing-strip boxes.
#
# Filter central + exact R6 cone WHILE STREAMING so memory stays small.
# ------------------------------------------------------------

r5_files=sorted(R5_DIR.glob("tile_*.csv"))

need(
    len(r5_files)==16,
    f"Expected 16 R5 raw tiles, found {len(r5_files)}"
)

new_files=[
    NEW_DIR/f"{k:02d}_{b['name']}.csv"
    for k,b in enumerate(boxes)
]

need(
    all(p.exists() for p in new_files),
    "Not all 7 hybrid missing-strip caches exist"
)

inputs=r5_files+new_files

by_gid={}

raw_rows=0
central_rows=0
exact_cone_central_rows_before_dedup=0
duplicate_count=0

for p in inputs:

    print(
        "E63B_R6_HYBRID_MERGE_INPUT",
        p.relative_to(ROOT),
        flush=True
    )

    with p.open(newline="") as f:
        r=csv.DictReader(f)

        need(
            set(REQ)<=set(r.fieldnames or []),
            f"Bad merge columns {p}"
        )

        for row in r:
            raw_rows+=1

            if not parse_bool(row["is_central"]):
                continue

            central_rows+=1

            ang=angular_sep_deg(
                row["ra_true"],
                row["dec_true"]
            )

            if ang>OUTER_RADIUS_DEG:
                continue

            exact_cone_central_rows_before_dedup+=1

            gid=int(float(row["galaxy_id"]))

            if gid in by_gid:
                duplicate_count+=1
            else:
                by_gid[gid]=row

need(
    len(by_gid)>=54893,
    f"Only {len(by_gid)} exact-cone unique centrals; "
    "below frozen NSEL=54893"
)

# Stable deterministic cache ordering.
gids=sorted(by_gid)

tmp=OUT.with_suffix(".tmp")
OUT.parent.mkdir(parents=True,exist_ok=True)

with tmp.open("w",newline="") as f:
    w=csv.DictWriter(
        f,
        fieldnames=REQ,
        extrasaction="ignore"
    )

    w.writeheader()

    for gid in gids:
        w.writerow(by_gid[gid])

os.replace(tmp,OUT)

# ------------------------------------------------------------
# Geometry-only real-catalog coverage guard.
# ------------------------------------------------------------

counts=[0]*24

for gid in gids:
    row=by_gid[gid]

    ang=angular_sep_deg(
        row["ra_true"],
        row["dec_true"]
    )

    if not (8.5<=ang<9.5):
        continue

    pa=position_angle_deg(
        row["ra_true"],
        row["dec_true"]
    )

    j=min(23,int(pa//15.0))
    counts[j]+=1

print()
print("=== E63B R6 HYBRID COVERAGE ===")

for j,n in enumerate(counts):
    print(
        "SECTOR",
        f"{j:02d}",
        f"PA=[{15*j},{15*(j+1)})",
        "N",n
    )

coverage_pass=min(counts)>=50

meta={
    "stage":"E63B_R6_HYBRID_RETRIEVAL",
    "truth_velocity_used":False,
    "science_metrics_computed":False,
    "observed_eBOSS_rows_used":False,
    "observed_odd_used":False,

    "r6_center":{
        "ra_deg":RA0,
        "dec_deg":DEC0,
        "outer_radius_deg":OUTER_RADIUS_DEG,
    },

    "reuse":{
        "r5_raw_tile_count":len(r5_files),
        "r5_rectangle_reuse_fraction":0.5171130952380952,
    },

    "new_retrieval":{
        "box_count":7,
        "boxes":queries,
    },

    "merge":{
        "raw_rows_streamed":raw_rows,
        "central_rows_streamed":central_rows,
        "exact_cone_central_rows_before_dedup":
            exact_cone_central_rows_before_dedup,
        "exact_cone_duplicate_rows":duplicate_count,
        "exact_cone_unique_centrals":len(gids),
        "cache":str(OUT.relative_to(ROOT)),
        "cache_sha256":sha_file(OUT),
    },

    "coverage":{
        "annulus_deg":[8.5,9.5],
        "sector_count":24,
        "minimum_required_per_sector":50,
        "counts":counts,
        "minimum_count":min(counts),
        "pass":coverage_pass,
    }
}

META.write_text(
    json.dumps(meta,indent=2)+"\n"
)

print()
print(
    "E63B_R6_HYBRID_MERGED_CACHE",
    OUT,
    "RAW_ROWS",raw_rows,
    "EXACT_CONE_UNIQUE_CENTRALS",len(gids),
    "DUPLICATES",duplicate_count,
    flush=True
)

print(
    "E63B_R6_HYBRID_COVERAGE_MIN",
    min(counts),
    "PASS",coverage_pass,
    flush=True
)

print(
    "E63B_R6_HYBRID_CACHE_SHA256",
    sha_file(OUT),
    flush=True
)

print(
    "E63B_R6_HYBRID_RETRIEVAL_COMPLETE_STOP_BEFORE_SCIENCE",
    flush=True
)

print("NO_TRUTH_VELOCITY_USED",flush=True)
print("NO_SCIENCE_METRICS_COMPUTED",flush=True)

need(
    coverage_pass,
    f"R6 real-catalog coverage guard FAIL counts={counts}"
)
