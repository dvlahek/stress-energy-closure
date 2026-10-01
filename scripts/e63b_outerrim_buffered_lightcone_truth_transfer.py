#!/usr/bin/env python3
"""
E63B: buffered OuterRim/CosmoDC2 lightcone truth transfer, MOCK/SIMULATION ONLY.

This stage follows:
  E62 periodic N-body truth PASS;
  E63 source/schema preflight PASS.

It keeps the E59/E62 reconstruction fixed and changes the geometry from a
periodic snapshot to an observer lightcone. Probes are explicitly buffered
from the queried cone/shell boundaries by >256 Mpc/h, so this stage tests
lightcone/evolution transfer before any explicit survey-mask subtraction.

The tracer number density is frozen to the exact E62 density:
165107 / (1000 Mpc/h)^3 = 1.65107e-4 h^3 Mpc^-3.

Observed eBOSS rows and the observed odd vector are never accessed.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import e63_cosmodc2_outerrim_lightcone_source_preflight as P

PREFLIGHT=ROOT/"source_data/e63_cosmodc2_outerrim_lightcone_source_preflight_result_compact_2026-10-01.json"
E62=ROOT/"source_data/e62_quijote_z1_truth_velocity_calibration_compact_summary_2026-10-01.json"
OUT=ROOT/"source_data/e63b_outerrim_buffered_lightcone_truth_transfer_result.json"
CACHE=ROOT/"eboss_workspace/cosmodc2/e63b_outerrim_central_candidates.csv"
LEGACY_CAPPED_CACHE=ROOT/"eboss_workspace/cosmodc2/e63b_outerrim_candidates.csv"

TABLE="cosmodc2mockv1"
RA0=55.0
DEC0=-41.0
OUTER_RADIUS_DEG=9.5
INNER_RADIUS_DEG=2.5
SOURCE_Z0=0.75
SOURCE_Z1=1.16
PROBE_Z0=0.9
PROBE_Z1=1.0
MASS_FLOOR=5.0e12
TOP_CAP=600000
NPROBE=512
PROBE_SEED=202610630001

OMEGA_M=0.2648
OMEGA_L=0.7352
H=0.71
C_KM_S=299792.458

N_TARGET=165107.0/1.0e9

R_SMOOTH=16.0
RMAX_PRIMARY=256.0
RMAX_CHECK=192.0

TRUTH_R_MIN=0.70
TRUTH_SIGN_MIN=0.70
CUTOFF_R_MIN=0.90
CUTOFF_SIGN_MIN=0.90

REQ=[
 "galaxy_id","halo_id","ra_true","dec_true","redshift_true",
 "position_x","position_y","position_z",
 "velocity_x","velocity_y","velocity_z",
 "halo_mass","is_central"
]

def need(c,msg):
    if not c:
        raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e63b_",delete=False) as f:
        q=Path(f.name)
        f.write(raw); f.flush(); os.fsync(f.fileno())
    try:
        os.replace(q,path)
    finally:
        q.unlink(missing_ok=True)

def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()

def chi_h(z):
    val=quad(lambda zz: 1.0/math.sqrt(OMEGA_M*(1+zz)**3+OMEGA_L),0,float(z),
             epsabs=1e-10,epsrel=1e-10,limit=100)[0]
    return (C_KM_S/100.0)*val

def chi_mpc(z):
    return chi_h(z)/H

def solid_angle_circle(radius_deg):
    t=math.radians(radius_deg)
    return 2.0*math.pi*(1.0-math.cos(t))

def angular_sep_deg(ra,dec):
    r=np.deg2rad(np.asarray(ra,float))
    d=np.deg2rad(np.asarray(dec,float))
    r0=math.radians(RA0); d0=math.radians(DEC0)
    c=np.sin(d)*math.sin(d0)+np.cos(d)*math.cos(d0)*np.cos(r-r0)
    return np.rad2deg(np.arccos(np.clip(c,-1.0,1.0)))

def position_angle_deg(ra,dec):
    r=np.deg2rad(np.asarray(ra,float)); d=np.deg2rad(np.asarray(dec,float))
    r0=math.radians(RA0); d0=math.radians(DEC0)
    dr=r-r0
    y=np.sin(dr)*np.cos(d)
    x=np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(dr)
    return (np.rad2deg(np.arctan2(y,x))+360.0)%360.0

def sky_unit(ra,dec):
    r=np.deg2rad(np.asarray(ra,float)); d=np.deg2rad(np.asarray(dec,float))
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])

def top_hat_kernel(delta):
    r2=np.einsum("ij,ij->i",delta,delta)
    r=np.sqrt(r2)
    fac=np.empty_like(r)
    m=r<R_SMOOTH
    fac[m]=1.0/(R_SMOOTH**3)
    fac[~m]=1.0/np.maximum(r[~m]**3,1e-300)
    return delta*fac[:,None],r

def reconstruct(source_pos,probe_pos,chunk=32):
    tree=cKDTree(source_pos,leafsize=32)
    primary=np.zeros((len(probe_pos),3),float)
    check=np.zeros((len(probe_pos),3),float)
    for a0 in range(0,len(probe_pos),chunk):
        a1=min(a0+chunk,len(probe_pos))
        nbs=tree.query_ball_point(probe_pos[a0:a1],RMAX_PRIMARY,workers=2)
        for local,js in enumerate(nbs):
            if not len(js):
                continue
            j=np.asarray(js,dtype=np.int64)
            delta=source_pos[j]-probe_pos[a0+local]
            K,r=top_hat_kernel(delta)
            primary[a0+local]=K.sum(axis=0)
            m=r<=RMAX_CHECK
            if np.any(m):
                check[a0+local]=K[m].sum(axis=0)
    primary/=float(len(source_pos))
    check/=float(len(source_pos))
    return primary,check

def pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    aa=a-np.mean(a); bb=b-np.mean(b)
    den=float(np.linalg.norm(aa)*np.linalg.norm(bb))
    return float(aa@bb/den) if den>0 else 0.0

def sign_agreement(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    nz=(a!=0)&(b!=0)
    need(np.count_nonzero(nz)>=0.95*len(a),"Too many zero values")
    return float(np.mean(np.sign(a[nz])==np.sign(b[nz])))

def vector_cosines(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    den=np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)
    m=den>0
    need(np.count_nonzero(m)>=0.95*len(a),"Too many zero vectors")
    return np.einsum("ij,ij->i",a[m],b[m])/den[m]

def query_string():
    cols=",".join(REQ)
    return (
      f"SELECT TOP {TOP_CAP} {cols} FROM {TABLE} WHERE "
      f"redshift_true>={SOURCE_Z0} AND redshift_true<{SOURCE_Z1} AND "
      f"halo_mass>={MASS_FLOOR:.1f} AND "
      f"is_central=1 AND "
      f"1=CONTAINS(POINT('ICRS',ra_true,dec_true),"
      f"CIRCLE('ICRS',{RA0},{DEC0},{OUTER_RADIUS_DEG}))"
    )

def load_or_query():
    q=query_string()
    CACHE.parent.mkdir(parents=True,exist_ok=True)
    if CACHE.is_file():
        raw=CACHE.read_text()
        print("E63B_REUSE_CACHED_CANDIDATES",CACHE,flush=True)
    else:
        print("E63B_TAP_QUERY_START",flush=True)
        raw=P.query_async(q,poll_seconds=5,max_wait_seconds=14400)
        tmp=CACHE.with_suffix(".tmp")
        tmp.write_text(raw)
        os.replace(tmp,CACHE)
        print("E63B_TAP_QUERY_CACHED",CACHE,flush=True)
    rows=P.parse_csv(raw)
    need(len(rows)<TOP_CAP,
         f"Candidate query hit TOP {TOP_CAP}; retrieval floor/cap must be revised before truth evaluation")
    return q,raw,rows

def parse_bool(v):
    s=str(v).strip().lower()
    if s in ("true","1","t","yes"): return True
    if s in ("false","0","f","no"): return False
    raise RuntimeError(f"Unexpected boolean value {v!r}")

def parse_rows(rows):
    n=len(rows)
    need(n>0,"No E63B candidate rows")
    out={}
    out["galaxy_id"]=np.asarray([int(float(r["galaxy_id"])) for r in rows],dtype=np.int64)
    out["halo_id"]=np.asarray([int(float(r["halo_id"])) for r in rows],dtype=np.int64)
    for k in ["ra_true","dec_true","redshift_true","position_x","position_y","position_z",
              "velocity_x","velocity_y","velocity_z","halo_mass"]:
        out[k]=np.asarray([float(r[k]) for r in rows],dtype=float)
        need(np.isfinite(out[k]).all(),f"Nonfinite {k}")
    out["is_central"]=np.asarray([parse_bool(r["is_central"]) for r in rows],dtype=bool)
    return out

def position_unit_guard(d):
    central=np.flatnonzero(d["is_central"])
    need(len(central)>=1024,"Too few central candidates for unit guard")
    take=central[np.linspace(0,len(central)-1,1024,dtype=int)]
    pos=np.column_stack([d["position_x"][take],d["position_y"][take],d["position_z"][take]])
    r=np.linalg.norm(pos,axis=1)
    ch=np.asarray([chi_mpc(z) for z in d["redshift_true"][take]],float)
    ratio=float(np.median(r/ch))
    if 0.97<=ratio<=1.03:
        scale=H
        mode="IRSA_GCR_Mpc_to_Mpc_over_h"
    elif H-0.025<=ratio<=H+0.025:
        scale=1.0
        mode="already_Mpc_over_h"
    else:
        raise RuntimeError(f"Unrecognized CosmoDC2 position unit geometry ratio={ratio}")
    posh=pos*scale
    u=posh/np.linalg.norm(posh,axis=1)[:,None]
    us=sky_unit(d["ra_true"][take],d["dec_true"][take])
    dots=np.einsum("ij,ij->i",u,us)
    need(float(np.median(dots))>=0.999 and float(np.min(dots))>=0.995,
         f"position/sky-direction mismatch median={np.median(dots)} min={np.min(dots)}")
    return scale,mode,ratio,float(np.median(dots)),float(np.min(dots))

def coverage_guard(d):
    m=d["is_central"]
    ang=angular_sep_deg(d["ra_true"][m],d["dec_true"][m])
    pa=position_angle_deg(d["ra_true"][m],d["dec_true"][m])
    ann=(ang>=8.5)&(ang<9.5)
    counts=[]
    for j in range(24):
        lo=15.0*j; hi=15.0*(j+1)
        counts.append(int(np.count_nonzero(ann&(pa>=lo)&(pa<hi))))
    need(min(counts)>=50,
         f"Outer-cone angular coverage guard failed sector counts={counts}")
    return counts

def metrics(rec,chk,truth,pos):
    los=pos/np.linalg.norm(pos,axis=1)[:,None]
    rp=np.einsum("ij,ij->i",rec,los)
    rc=np.einsum("ij,ij->i",chk,los)
    tv=np.einsum("ij,ij->i",truth,los)
    out={
      "LOS":{
        "truth_pearson":pearson(rp,tv),
        "truth_sign_agreement":sign_agreement(rp,tv),
        "cutoff_192_vs_256_pearson":pearson(rp,rc),
        "cutoff_192_vs_256_sign_agreement":sign_agreement(rp,rc),
      },
      "cartesian":{},
    }
    for j,a in enumerate(("x","y","z")):
        out["cartesian"][a]={
          "truth_pearson":pearson(rec[:,j],truth[:,j]),
          "truth_sign_agreement":sign_agreement(rec[:,j],truth[:,j]),
          "cutoff_192_vs_256_pearson":pearson(rec[:,j],chk[:,j]),
          "cutoff_192_vs_256_sign_agreement":sign_agreement(rec[:,j],chk[:,j]),
        }
    vc=vector_cosines(rec,truth)
    out["vector_cosine"]={
      "mean":float(np.mean(vc)),"median":float(np.median(vc)),
      "min":float(np.min(vc)),"max":float(np.max(vc))
    }
    L=out["LOS"]
    out["LOS"]["gate_pass"]=bool(
      L["truth_pearson"]>=TRUTH_R_MIN
      and L["truth_sign_agreement"]>=TRUTH_SIGN_MIN
      and L["cutoff_192_vs_256_pearson"]>=CUTOFF_R_MIN
      and L["cutoff_192_vs_256_sign_agreement"]>=CUTOFF_SIGN_MIN
    )
    return out

def self_test():
    need(abs(N_TARGET-0.000165107)<1e-15,"density constant drift")
    need("is_central=1" in query_string(),"central-only TAP retrieval predicate missing")
    need(chi_h(PROBE_Z0)-chi_h(SOURCE_Z0)>256.0,"lower radial buffer prereg invalid")
    need(chi_h(SOURCE_Z1)-chi_h(PROBE_Z1)>256.0,"upper radial buffer prereg invalid")
    need(chi_h(PROBE_Z0)*math.sin(math.radians(OUTER_RADIUS_DEG-INNER_RADIUS_DEG))>256.0,
         "angular buffer prereg invalid")
    d=np.array([[0.,0.,0.],[8.,0.,0.],[16.,0.,0.],[32.,0.,0.]])
    K,r=top_hat_kernel(d)
    need(np.allclose(K[0],0),"kernel zero fail")
    need(np.isclose(K[1,0],8/R_SMOOTH**3),"kernel inner fail")
    need(np.isclose(K[3,0],32/32**3),"kernel outer fail")
    print("E63B_SYNTHETIC_BUFFERED_LIGHTCONE_SELF_TEST_PASS",flush=True)

def main():
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=32)
    args=ap.parse_args()
    self_test()
    if args.self_test:return
    need(args.run,"Use --run for E63B")
    need(8<=args.chunk<=128,"Unsafe chunk")
    need(PREFLIGHT.is_file() and E62.is_file(),"Missing E63 preflight or E62 compact parent")
    pf=json.loads(PREFLIGHT.read_text()); e62=json.loads(E62.read_text())
    need(pf["status"]=="PASS_PUBLIC_OUTERRIM_LIGHTCONE_SOURCE_AVAILABLE"
         and pf["decision"]["source_preflight_pass"] is True,
         "E63 source preflight parent changed")
    need(e62["status"]=="PASS_TRUTH_CALIBRATION_QUANTIFIED"
         and e62["decision"]["truth_velocity_gate_pass"] is True,
         "E62 truth parent changed")
    need(pf["observed_eBOSS_rows_used"] is False and pf["observed_odd_used"] is False,
         "Observation guardrail changed")

    q,raw,rows=load_or_query()
    print("E63B_CANDIDATE_ROWS",len(rows),flush=True)
    d=parse_rows(rows)
    central=np.flatnonzero(d["is_central"])
    need(len(central)>0,"No central candidates")
    need(len(central)==len(rows),
         f"Server-side central-only retrieval leaked non-central rows: {len(rows)-len(central)}")
    need(len(np.unique(d["halo_id"][central]))==len(central),
         "Central sample has duplicate halo_id values")

    coverage=coverage_guard(d)
    scale,unit_mode,ratio,sky_med,sky_min=position_unit_guard(d)

    omega=solid_angle_circle(OUTER_RADIUS_DEG)
    r0=chi_h(SOURCE_Z0); r1=chi_h(SOURCE_Z1)
    volume=omega*(r1**3-r0**3)/3.0
    nsel=int(round(N_TARGET*volume))
    need(len(central)>=nsel,
         f"Only {len(central)} central candidates for frozen target {nsel}")

    mass=d["halo_mass"][central]
    gid=d["galaxy_id"][central]
    order=np.lexsort((gid,-mass))
    selected=central[order[:nsel]]
    minmass=float(np.min(d["halo_mass"][selected]))
    need(minmass>=1.10*MASS_FLOOR,
         f"Retrieval floor too close to science threshold: min selected mass={minmass}")

    pos=np.column_stack([d["position_x"][selected],d["position_y"][selected],d["position_z"][selected]])*scale
    vel=np.column_stack([d["velocity_x"][selected],d["velocity_y"][selected],d["velocity_z"][selected]])
    zang=d["redshift_true"][selected]
    aang=angular_sep_deg(d["ra_true"][selected],d["dec_true"][selected])

    pmask=(aang<=INNER_RADIUS_DEG)&(zang>=PROBE_Z0)&(zang<PROBE_Z1)
    pool=np.flatnonzero(pmask)
    need(len(pool)>=NPROBE,
         f"Only {len(pool)} selected inner probes < frozen {NPROBE}")
    rng=np.random.default_rng(PROBE_SEED)
    pidx=np.sort(rng.choice(pool,size=NPROBE,replace=False))
    ppos=pos[pidx]; ptruth=vel[pidx]

    pr=np.linalg.norm(ppos,axis=1)
    angular_buffer=pr*math.sin(math.radians(OUTER_RADIUS_DEG-INNER_RADIUS_DEG))
    lower_buffer=pr-r0
    upper_buffer=r1-pr
    need(float(np.min(angular_buffer))>RMAX_PRIMARY,
         f"Angular buffer below 256: {np.min(angular_buffer)}")
    need(float(np.min(lower_buffer))>RMAX_PRIMARY,
         f"Lower radial buffer below 256: {np.min(lower_buffer)}")
    need(float(np.min(upper_buffer))>RMAX_PRIMARY,
         f"Upper radial buffer below 256: {np.min(upper_buffer)}")

    print("E63B_RECONSTRUCTION_START SOURCES",len(pos),"PROBES",len(ppos),flush=True)
    rec,chk=reconstruct(pos,ppos,args.chunk)
    met=metrics(rec,chk,ptruth,ppos)
    passed=bool(met["LOS"]["gate_pass"])

    result={
      "stage":"E63B_OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER",
      "date":"2026-10-01",
      "status":"PASS_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED" if passed else "FAIL_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED",
      "source":{
        "table":TABLE,"tap_endpoint":"https://irsa.ipac.caltech.edu/TAP",
        "query":q,"cache":str(CACHE),"cache_sha256":sha256_bytes(raw.encode()),
        "legacy_all_galaxy_capped_cache":str(LEGACY_CAPPED_CACHE),
        "server_side_central_filter":True,
        "retrieval_complete_below_TOP_cap":bool(len(rows)<TOP_CAP),
        "candidate_rows":len(rows),"central_candidate_rows":len(central),
        "outer_radius_deg":OUTER_RADIUS_DEG,
        "source_redshift_true":[SOURCE_Z0,SOURCE_Z1],
        "candidate_mass_floor_Msun":MASS_FLOOR
      },
      "geometry":{
        "outerrim_cosmology":{"Omega_m":OMEGA_M,"Omega_Lambda":OMEGA_L,"h":H},
        "position_unit_mode":unit_mode,
        "median_position_radius_over_chi_Mpc":ratio,
        "position_sky_direction_dot_median":sky_med,
        "position_sky_direction_dot_min":sky_min,
        "outer_cone_solid_angle_sr":omega,
        "outer_geometric_volume_Mpc_over_h_cubed":volume,
        "outer_annulus_24_sector_central_candidate_counts":coverage
      },
      "tracer":{
        "target_number_density_h3_Mpc3":N_TARGET,
        "selected_count":len(selected),
        "achieved_number_density_h3_Mpc3":len(selected)/volume,
        "minimum_selected_halo_mass_Msun":minmass,
        "maximum_selected_halo_mass_Msun":float(np.max(d["halo_mass"][selected])),
        "central_only":True,"unit_number_weights":True
      },
      "probes":{
        "count":NPROBE,"seed":PROBE_SEED,
        "inner_radius_deg":INNER_RADIUS_DEG,
        "redshift_true":[PROBE_Z0,PROBE_Z1],
        "eligible_selected_inner_count":len(pool),
        "minimum_angular_boundary_distance_Mpc_over_h":float(np.min(angular_buffer)),
        "minimum_lower_radial_boundary_distance_Mpc_over_h":float(np.min(lower_buffer)),
        "minimum_upper_radial_boundary_distance_Mpc_over_h":float(np.min(upper_buffer)),
        "selected_source_indices_sha256":sha256_bytes(np.ascontiguousarray(selected).tobytes()),
        "probe_selected_indices_sha256":sha256_bytes(np.ascontiguousarray(pidx).tobytes())
      },
      "reconstruction":{
        "smoothing_R_Mpc_over_h":R_SMOOTH,
        "primary_cutoff_Mpc_over_h":RMAX_PRIMARY,
        "check_cutoff_Mpc_over_h":RMAX_CHECK,
        "random_subtraction_used":False,
        "same_kernel_as_E59_E62":True
      },
      "gates":{
        "LOS_truth_pearson_min":TRUTH_R_MIN,
        "LOS_truth_sign_agreement_min":TRUTH_SIGN_MIN,
        "LOS_cutoff_pearson_min":CUTOFF_R_MIN,
        "LOS_cutoff_sign_agreement_min":CUTOFF_SIGN_MIN
      },
      "metrics":met,
      "decision":{
        "buffered_lightcone_truth_gate_pass":passed,
        "interpretation":(
          "OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_PASS; NEXT_E64_SURVEY_MASK_HOD_TRUTH_TRANSFER"
          if passed else
          "OUTERRIM_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_FAIL; DO_NOT_ADD_SURVEY_MASK_OR_OBSERVED_TAG"
        ),
        "important_limit":"E63B isolates snapshot-to-lightcone transfer at the E62 tracer density using central halos. It does not yet test eBOSS angular/radial selection, ELG HOD/photometric targeting, or absolute Einstein-Vlasov amplitude."
      },
      "observed_eBOSS_rows_used":False,
      "observed_odd_used":False,
      "absolute_EV_amplitude_calibrated":False,
      "pvalue_or_detection_sigma":False
    }
    atomic(OUT,result)
    print("E63B_BUFFERED_LIGHTCONE_COMPLETE",flush=True)
    L=met["LOS"]
    print("LOS R_TRUTH",L["truth_pearson"],
          "SIGN_TRUTH",L["truth_sign_agreement"],
          "R_CUTOFF",L["cutoff_192_vs_256_pearson"],
          "SIGN_CUTOFF",L["cutoff_192_vs_256_sign_agreement"],
          "PASS",L["gate_pass"],flush=True)
    print("VECTOR_COS_MED",met["vector_cosine"]["median"],flush=True)
    print("DECISION",result["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
