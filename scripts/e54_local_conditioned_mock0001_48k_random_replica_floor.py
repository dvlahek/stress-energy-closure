#!/usr/bin/env python3
"""
E54: finite-random floor from fresh 48k random replicas at fixed mock0001 galaxies.

The E53 Stage-2 result showed that the conditioned F+/F- fingerprint is stable,
but 4800R -> 48000R still changes the marked-LS background materially.

E54 therefore holds EVERYTHING scientific fixed:
  - mock0001 full eligible galaxies,
  - E51 four technical sign fields,
  - E19 F+/F- injection basis,
  - pair geometry, s/mu bins, orientation and projection,

and changes only the Monte-Carlo random catalogue. Four fresh 48k random
subsamples are drawn independently (within each replicate, without replacement)
from the same frozen source-eligible random pool for each cap/tracer.

These are technical random replicas, not new cosmological mocks and not a
physical velocity reconstruction. Replicas can overlap because the source pool
is finite; overlap fractions are reported explicitly.

NO observed galaxy rows. NO observed odd vector. NO covariance inverse/p-values.
"""
from __future__ import annotations

import argparse, json, math, os, sys, tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import run_eboss_dr16_full_eligible_mock0001_sparse as V1
import e53_local_conditioned_mock0001_sample_scaling as E53

PARENT=ROOT/"source_data/e53_conditioned_mock0001_sample_scaling.json"
OUT=ROOT/"source_data/e54_conditioned_mock0001_48k_random_replica_floor.json"

CAPS=("NGC","SGC")
TRACERS=V1.TRACERS
ROLES=V1.ROLES
FIELDS=E53.FIELDS
FIELD_NAMES=E53.FIELD_NAMES
TERMS=("D1D2","D1R2","R1D2","R1R2")

REPS=("R1","R2","R3","R4")
N_RANDOM=48000
SEED_ROOT=202609305400

CAL_SD_MEDIAN_MAX=0.01
CAL_SD_FIELD_MAX=0.02
ANGLE_SPREAD_MAX_RAD=0.002
DIFF_RATIO_SPREAD_MAX=0.01

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e54_",delete=False) as f:
        t=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(t,path)
    finally: t.unlink(missing_ok=True)

def rep_seed(rep,cap,tracer):
    return (SEED_ROOT + 10000*(REPS.index(rep)+1)
            + 1000*CAPS.index(cap) + 100*TRACERS.index(tracer))

def full_random_pool(path,cap,tracer,oldmeta):
    with fits.open(path,memmap=False) as hdus:
        hdus.verify("exception")
        tabs=[h for h in hdus if isinstance(h,fits.BinTableHDU)]
        need(len(tabs)==1 and int(tabs[0].header["NAXIS2"])==oldmeta["input_FITS_rows"],
             "Frozen random source BINTABLE/header differs")
        tab=tabs[0]
        needed={"RA","DEC","Z",*V1.WEIGHT_COLUMNS}
        if tracer=="eBOSS_ELG": needed.add("chunk")
        need(needed.issubset(set(tab.columns.names)),
             "Random source missing frozen columns")
        d=tab.data
        ra,dec,z=(np.asarray(d[n],dtype="f8") for n in ("RA","DEC","Z"))
        need(all(np.isfinite(v).all() for v in (ra,dec,z)),
             "Nonfinite random coordinates")
        candidate=(z>=.9)&(z<1.)
        weights,retained=V1.validated_weight_product({
            n:np.asarray(d[n][candidate],dtype="f8") for n in V1.WEIGHT_COLUMNS
        })
        eligible=np.flatnonzero(candidate)[retained]
        valid_weight=weights[retained]
        need(int(candidate.sum())==oldmeta["candidate_highz_rows_before_weight_validation"]
             and len(eligible)==oldmeta["eligible_highz_rows_after_fixed_weight_gate"]
             and np.isfinite(valid_weight).all() and np.all(valid_weight>0),
             "Frozen source-eligible random pool drift")
        return (ra[eligible].copy(),dec[eligible].copy(),
                z[eligible].copy(),valid_weight.copy())

def sample_pool(pool,seed,n=N_RANDOM):
    need(len(pool[0])>=n,"Eligible random pool smaller than requested replica")
    pos=np.sort(np.random.default_rng(seed).choice(len(pool[0]),size=n,replace=False))
    cat=tuple(np.asarray(v[pos],dtype="f8").copy() for v in pool)
    return cat,pos

def term_to_json(rec):
    return E53.encode_term(rec)

def summary_metrics(summary):
    out={
      "postwindow_angle_rad":float(summary["postwindow_angle_rad"]),
      "difference_norm_over_plus_norm":float(summary["difference_norm_over_plus_norm"]),
      "shape_residual_norm_over_plus_norm":float(summary["shape_residual_norm_over_plus_norm"]),
      "rho_minus_on_plus":float(summary["rho_minus_on_plus"]),
      "by_tag_field":{}
    }
    for f in FIELD_NAMES:
        z=summary["by_tag_field"][f]
        out["by_tag_field"][f]={
          "baseline_L2":float(z["baseline_L2"]),
          "calibrated_diff":float(z["calibrated_Fminus_minus_Fplus_coordinate"]),
          "plus_lambda":float(z["plus_apparent_lambda"]),
          "minus_lambda":float(z["minus_apparent_lambda"]),
          "shape_only":float(z["shape_only_coordinate"]),
        }
    return out

def descriptive(vals):
    a=np.asarray(vals,float)
    return {
      "n":int(len(a)),
      "mean":float(np.mean(a)),
      "median":float(np.median(a)),
      "sample_sd":float(np.std(a,ddof=1)),
      "MAD":float(np.median(np.abs(a-np.median(a)))),
      "min":float(np.min(a)),
      "max":float(np.max(a)),
      "naive_se_if_independent":float(np.std(a,ddof=1)/math.sqrt(len(a))),
    }

def aggregate(case):
    reps=[case["replicas"][r]["summary"] for r in REPS]
    out={
      "postwindow_angle_rad":descriptive([z["postwindow_angle_rad"] for z in reps]),
      "difference_norm_over_plus_norm":descriptive(
          [z["difference_norm_over_plus_norm"] for z in reps]),
      "by_tag_field":{}
    }
    for f in FIELD_NAMES:
        out["by_tag_field"][f]={}
        for key in ("baseline_L2","calibrated_diff","plus_lambda","minus_lambda","shape_only"):
            out["by_tag_field"][f][key]=descriptive(
                [z["by_tag_field"][f][key] for z in reps])
    angle_spread=out["postwindow_angle_rad"]["max"]-out["postwindow_angle_rad"]["min"]
    diff_spread=(out["difference_norm_over_plus_norm"]["max"]
                 -out["difference_norm_over_plus_norm"]["min"])
    cal_sds=[out["by_tag_field"][f]["calibrated_diff"]["sample_sd"] for f in FIELD_NAMES]
    out["engineering_gate"]={
      "angle_spread_rad":float(angle_spread),
      "difference_ratio_spread":float(diff_spread),
      "calibrated_diff_sd_by_field":dict(zip(FIELD_NAMES,map(float,cal_sds))),
      "calibrated_diff_sd_median":float(np.median(cal_sds)),
      "calibrated_diff_sd_max":float(np.max(cal_sds)),
      "geometry_stable":bool(angle_spread<=ANGLE_SPREAD_MAX_RAD and
                             diff_spread<=DIFF_RATIO_SPREAD_MAX),
      "48k_random_floor_small_on_unit_lambda_scale":bool(
          np.median(cal_sds)<=CAL_SD_MEDIAN_MAX and
          np.max(cal_sds)<=CAL_SD_FIELD_MAX),
      "thresholds":{
        "cal_sd_median_max":CAL_SD_MEDIAN_MAX,
        "cal_sd_field_max":CAL_SD_FIELD_MAX,
        "angle_spread_max_rad":ANGLE_SPREAD_MAX_RAD,
        "difference_ratio_spread_max":DIFF_RATIO_SPREAD_MAX,
      }
    }
    return out

def overlap_summary(pos_by_rep):
    out={}
    for tr in TRACERS:
        d={}
        for i,a in enumerate(REPS):
            for b in REPS[i+1:]:
                A=pos_by_rep[a][tr]; B=pos_by_rep[b][tr]
                d[a+"__"+b]=float(len(np.intersect1d(A,B,assume_unique=True))/N_RANDOM)
        out[tr]={
          "pairwise_fractional_overlap":d,
          "median_pairwise_fractional_overlap":float(np.median(list(d.values()))),
        }
    return out

def synthetic_self_test():
    need(len(REPS)==4 and len(set(REPS))==4,"Replica contract corrupted")
    seeds=[rep_seed(r,c,t) for r in REPS for c in CAPS for t in TRACERS]
    need(len(seeds)==len(set(seeds)),"Replica seeds not unique")
    rng=np.random.default_rng(54001)
    pool=tuple(rng.normal(size=90000) for _ in range(4))
    _,p1=sample_pool(pool,rep_seed("R1","NGC",TRACERS[0]))
    _,p2=sample_pool(pool,rep_seed("R2","NGC",TRACERS[0]))
    need(len(p1)==N_RANDOM and len(np.unique(p1))==N_RANDOM,
         "Replica sampling without-replacement failed")
    ov=len(np.intersect1d(p1,p2,assume_unique=True))/N_RANDOM
    need(0<ov<1,"Replica overlap diagnostic failed")
    print("E54_SYNTHETIC_RANDOM_REPLICA_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()

    synthetic_self_test()
    if args.self_test: return
    need(args.run,"Use --run for the mock-only E54 calculation")
    need(PARENT.is_file(),"Missing local E53 Stage2 parent JSON")
    parent=json.loads(PARENT.read_text())
    need(parent["status"]=="PASS_FULLD_4800_AND_48000_CONDITIONED_SAMPLE_SCALING",
         "E54 requires completed E53 Stage2 parent")
    need(parent["observed_odd_used"] is False and
         parent["observed_galaxy_rows_used"] is False,
         "E53 parent guardrail changed")

    p,a02,e4=V1.source_only_gate()
    original_protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,p0=A02.preflight(
        original_protocol,require_local=True)
    paths,total=A02.resolve_72_sources(
        original_protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state.get("stage")=="E54_CONDITIONED_MOCK0001_48K_RANDOM_REPLICA_FLOOR"
             and state.get("observed_odd_used") is False,
             "Existing E54 checkpoint incompatible")
    else:
        state={
          "stage":"E54_CONDITIONED_MOCK0001_48K_RANDOM_REPLICA_FLOOR",
          "date":"2026-09-30",
          "status":"INCOMPLETE",
          "mock_id":1,
          "random_size_per_tracer":N_RANDOM,
          "fresh_replicas":list(REPS),
          "seed_root":SEED_ROOT,
          "same_E51_fields":True,
          "same_E19_injection":True,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,
          "observed_odd_used":False,
          "physical_velocity_reconstruction_used":False,
          "pvalue_or_detection_sigma":False,
          "covariance_inverse_used":False,
        }
        atomic(OUT,state)

    for cap in CAPS:
        ck=f"0001/{cap}"
        case=state["cases"].setdefault(ck,{
          "cap":cap,"status":"INCOMPLETE","replicas":{},
          "random_pool_rows":{},"random_replica_positions_SHA256":{}
        })
        metadata=a02["cases"][ck]["input_sample_diagnostics"]
        prior=e4["cases"][ck]
        cats={}; info={}; pools={}

        for tracer in TRACERS:
            for role in ROLES:
                label=tracer+"_"+role
                src=paths[A02.key(1,cap,tracer)+"/"+role]
                cats[(tracer,role)],info[label]=V1.load_original_full(
                    src,cap,tracer,role,metadata[label],prior,p0)
                need(info[label]["full_eligible_rows"]==
                     parent["cases"][ck]["source_full_eligible_diagnostics"][label]["full_eligible_rows"],
                     "E53 full-eligible row count drift "+ck+"/"+label)
                if role=="ran":
                    pools[tracer]=full_random_pool(src,cap,tracer,metadata[label])
                    need(len(pools[tracer][0])==info[label]["full_eligible_rows"],
                         "Replica pool length drift")
                    case["random_pool_rows"][tracer]=len(pools[tracer][0])

        dist=V1.cached_exact_distance(cats)

        parent_dd=parent["cases"][ck]["levels"]["nested_48000"]["terms"]["D1D2"]
        dl=cats["eBOSS_LRG","dat"]["full"]
        de=cats["eBOSS_ELG","dat"]["full"]

        pos_by_rep={}
        for rep in REPS:
            rr=case["replicas"].setdefault(rep,{"status":"INCOMPLETE","terms":{}})
            repcats={}
            pos_by_rep[rep]={}
            for tracer in TRACERS:
                seed=rep_seed(rep,cap,tracer)
                cat,pos=sample_pool(pools[tracer],seed)
                repcats[tracer]=cat
                pos_by_rep[rep][tracer]=pos
                case["random_replica_positions_SHA256"].setdefault(rep,{})[tracer]=V1.E3.arr_sha(pos)
                rr.setdefault("seeds",{})[tracer]=seed
                rr.setdefault("catalogue_SHA256",{})[tracer]=V1.catalogue_sha(cat)

            samples={
              "D1R2":(dl,repcats["eBOSS_ELG"]),
              "R1D2":(repcats["eBOSS_LRG"],de),
              "R1R2":(repcats["eBOSS_LRG"],repcats["eBOSS_ELG"]),
            }
            for term in ("D1R2","R1D2","R1R2"):
                if term in rr["terms"]:
                    print("E54_REUSE_TERM",ck,rep,term,flush=True)
                    continue
                rec=E53.sparse_marked(*samples[term],dist,chunk=args.chunk,
                                      max_candidates=args.max_candidate_pairs_per_term)
                rr["terms"][term]=term_to_json(rec)
                atomic(OUT,state)
                print("E54_TERM_PASS",ck,rep,term,
                      "CAND",rec["meta"]["candidate_neighbour_pairs"],
                      "ACC",rec["meta"]["accepted_pairs"],flush=True)

            termmap={"D1D2":parent_dd,**rr["terms"]}
            summary=E53.summarize_level(termmap)
            rr["summary"]=summary_metrics(summary)
            rr["status"]="complete"
            atomic(OUT,state)
            print("E54_REPLICA_PASS",ck,rep,
                  "ANGLE",rr["summary"]["postwindow_angle_rad"],flush=True)

        case["fresh_replica_overlap"]=overlap_summary(pos_by_rep)
        case["aggregate"]=aggregate(case)
        case["E53_nested_48000_reference"]=summary_metrics(
            parent["cases"][ck]["levels"]["nested_48000"]["summary"])
        case["status"]="complete"
        atomic(OUT,state)

    gates=[state["cases"][f"0001/{c}"]["aggregate"]["engineering_gate"] for c in CAPS]
    all_geom=all(g["geometry_stable"] for g in gates)
    all_floor=all(g["48k_random_floor_small_on_unit_lambda_scale"] for g in gates)
    state["decision"]={
      "geometry_stable_all_caps":all_geom,
      "48k_random_floor_small_all_caps":all_floor,
      "engineering_interpretation":(
        "48K_REPLICA_FLOOR_SMALL_ENOUGH_FOR_NINEMOCK_ENGINEERING_EXPANSION"
        if all_geom and all_floor else
        "48K_REPLICA_FLOOR_STILL_MATERIAL; USE_REPLICA_AVERAGING_OR_HIGHER_DENSITY_BEFORE_NINEMOCK_EXPANSION"
      ),
      "important_limit":"This gate is relative to unit technical lambda, not a physical eBOSS detection threshold."
    }
    state["status"]="PASS_RANDOM_REPLICA_FLOOR_QUANTIFIED"
    atomic(OUT,state)

    print("E54_WSL_RANDOM_REPLICA_FLOOR_PASS",flush=True)
    for cap in CAPS:
        a=state["cases"][f"0001/{cap}"]["aggregate"]
        g=a["engineering_gate"]
        print("CAP",cap,
              "ANGLE_SPREAD",g["angle_spread_rad"],
              "DIFF_RATIO_SPREAD",g["difference_ratio_spread"],
              "CAL_SD_MED",g["calibrated_diff_sd_median"],
              "CAL_SD_MAX",g["calibrated_diff_sd_max"],
              "GEOM_OK",g["geometry_stable"],
              "RANDOM_FLOOR_OK",g["48k_random_floor_small_on_unit_lambda_scale"],
              flush=True)
        for f in FIELD_NAMES:
            print("TAG",f,
                  "CAL_SD",a["by_tag_field"][f]["calibrated_diff"]["sample_sd"],
                  "CAL_MEAN",a["by_tag_field"][f]["calibrated_diff"]["mean"],
                  "L2_MEAN",a["by_tag_field"][f]["baseline_L2"]["mean"],flush=True)
    print("DECISION",state["decision"]["engineering_interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
