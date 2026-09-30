#!/usr/bin/env python3
"""
E55: full-eligible nine-mock conditioned estimator with two fresh 48k random replicas.

Purpose
-------
E54 showed that one 48k random realization has a sufficiently small engineering
floor on the unit technical-lambda scale. E55 now measures the *mock-to-mock*
background of the same frozen conditioned estimator using all nine pre-fixed
eBOSS EZmock IDs.

For each mock/cap:
  - use ALL frozen source-eligible LRG and ELG galaxies;
  - use two predeclared fresh 48k random replicas (R1,R2) per tracer;
  - keep the E51 four technical sign fields and E19 F+/F- basis unchanged;
  - average R1/R2 only after each full marked-LS estimate is complete;
  - compress to the one-dimensional calibrated F--F+ coordinate.

Mock0001 R1/R2 seed definitions are exactly the E54 seeds. If the local E54
result exists, its completed R1/R2 summaries are reused verbatim and checked.

NO observed galaxy rows. NO observed odd vector. NO physical velocity
reconstruction. NO inverse covariance. NO p-values or detection significance.
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

E54=ROOT/"source_data/e54_conditioned_mock0001_48k_random_replica_floor.json"
OUT=ROOT/"source_data/e55_conditioned_full_eligible_ninemock_two48k_result.json"

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES
FIELDS=E53.FIELDS
FIELD_NAMES=E53.FIELD_NAMES
REPS=("R1","R2")
N_RANDOM=48000
SEED_ROOT=202609305400

RANDOM_TO_BETWEEN_MEDIAN_MAX=0.50
RANDOM_TO_BETWEEN_FIELD_MAX=1.00

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e55_",delete=False) as f:
        t=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(t,path)
    finally: t.unlink(missing_ok=True)

def rep_seed(mid,rep,cap,tracer):
    return (SEED_ROOT
            + 1_000_000*IDS.index(mid)
            + 10_000*(REPS.index(rep)+1)
            + 1_000*CAPS.index(cap)
            + 100*TRACERS.index(tracer))

def catalogue_sha(cat):
    return V1.catalogue_sha(cat)

def load_full_verified(path,*,mid,cap,tracer,role,expected_rows,
                       archived_info,origproto):
    fixed,info=A02.sample_catalogue(
        path,expected_rows=expected_rows,cap=cap,tracer=tracer,
        role=role,expected_highz=None,p=origproto)
    need(info["selected_array_SHA256"]==archived_info["selected_array_SHA256"],
         f"A02 selected SHA drift {mid}/{cap}/{tracer}/{role}")
    need(info["eligible_highz_rows_after_fixed_weight_gate"]==
         archived_info["eligible_highz_rows_after_fixed_weight_gate"],
         "A02 eligible count drift")

    with fits.open(path,memmap=False) as hdus:
        hdus.verify("exception")
        tabs=[h for h in hdus if isinstance(h,fits.BinTableHDU)]
        need(len(tabs)==1 and int(tabs[0].header["NAXIS2"])==expected_rows,
             "Frozen FITS row count drift")
        tab=tabs[0]
        needed={"RA","DEC","Z",*V1.WEIGHT_COLUMNS}
        if tracer=="eBOSS_ELG": needed.add("chunk")
        need(needed.issubset(set(tab.columns.names)),"Missing frozen columns")
        d=tab.data
        ra,dec,z=(np.asarray(d[n],dtype="f8") for n in ("RA","DEC","Z"))
        need(all(np.isfinite(v).all() for v in (ra,dec,z)),
             "Nonfinite source coordinates")
        cand=(z>=.9)&(z<1.)
        weights,retained=V1.validated_weight_product({
          n:np.asarray(d[n][cand],dtype="f8") for n in V1.WEIGHT_COLUMNS
        })
        eligible=np.flatnonzero(cand)[retained]
        valid_weight=weights[retained]
        need(len(eligible)==archived_info["eligible_highz_rows_after_fixed_weight_gate"]
             and np.isfinite(valid_weight).all() and np.all(valid_weight>0),
             "Full eligible source population drift")
        full=(ra[eligible].copy(),dec[eligible].copy(),
              z[eligible].copy(),valid_weight.copy())
    return fixed,full,info

def sample_random(pool,seed):
    need(len(pool[0])>=N_RANDOM,"Eligible random pool < 48k")
    pos=np.sort(np.random.default_rng(seed).choice(
        len(pool[0]),size=N_RANDOM,replace=False))
    cat=tuple(np.asarray(v[pos],dtype="f8").copy() for v in pool)
    return cat,pos

def exact_distance_cache(cats):
    zz=np.unique(np.concatenate([np.asarray(c[2],dtype="f8") for c in cats]))
    rr=V1.comoving_mpc_over_h(zz,V1.PRIMARY_GEOMETRY)
    need(np.isfinite(rr).all() and np.all(np.diff(rr)>0),
         "Invalid exact distance cache")
    def dist(z):
        vals=np.asarray(z,dtype="f8")
        ii=np.searchsorted(zz,vals)
        need(np.all(ii<len(zz)) and np.array_equal(zz[ii],vals),
             "Distance cache miss")
        return rr[ii]
    return dist

def compact_summary(summary):
    return {
      "postwindow_angle_rad":float(summary["postwindow_angle_rad"]),
      "difference_norm_over_plus_norm":float(summary["difference_norm_over_plus_norm"]),
      "shape_residual_norm_over_plus_norm":float(summary["shape_residual_norm_over_plus_norm"]),
      "rho_minus_on_plus":float(summary["rho_minus_on_plus"]),
      "by_tag_field":{
        f:{
          "baseline_L2":float(summary["by_tag_field"][f]["baseline_L2"]),
          "calibrated_diff":float(summary["by_tag_field"][f]["calibrated_Fminus_minus_Fplus_coordinate"]),
          "plus_lambda":float(summary["by_tag_field"][f]["plus_apparent_lambda"]),
          "minus_lambda":float(summary["by_tag_field"][f]["minus_apparent_lambda"]),
          "shape_only":float(summary["by_tag_field"][f]["shape_only_coordinate"]),
        } for f in FIELD_NAMES
      }
    }

def mean_pair(a,b):
    return .5*(float(a)+float(b))

def describe(v):
    a=np.asarray(v,float)
    med=float(np.median(a))
    av=np.sort(np.abs(a)); n=len(a)
    k83=6
    k94=8
    return {
      "n":n,"mean":float(np.mean(a)),"median":med,
      "sample_sd":float(np.std(a,ddof=1)),
      "MAD":float(np.median(np.abs(a-med))),
      "min":float(np.min(a)),"max":float(np.max(a)),
      "descriptive_abs_threshold_83pct":float(av[k83-1]),
      "descriptive_abs_threshold_94pct":float(av[k94-1]),
      "values":[float(x) for x in a],
    }

def aggregate(state):
    out={}
    for cap in CAPS:
        cs=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        capout={"by_tag_field":{}}
        angle_means=[mean_pair(c["replicas"]["R1"]["summary"]["postwindow_angle_rad"],
                               c["replicas"]["R2"]["summary"]["postwindow_angle_rad"])
                     for c in cs]
        diff_means=[mean_pair(c["replicas"]["R1"]["summary"]["difference_norm_over_plus_norm"],
                              c["replicas"]["R2"]["summary"]["difference_norm_over_plus_norm"])
                    for c in cs]
        capout["postwindow_angle_mock_means"]=describe(angle_means)
        capout["difference_norm_over_plus_mock_means"]=describe(diff_means)
        ratios=[]
        for f in FIELD_NAMES:
            mockmeans=[]; diffs=[]
            for c in cs:
                a=c["replicas"]["R1"]["summary"]["by_tag_field"][f]["calibrated_diff"]
                b=c["replicas"]["R2"]["summary"]["by_tag_field"][f]["calibrated_diff"]
                mockmeans.append(mean_pair(a,b))
                diffs.append(float(a)-float(b))
            between=float(np.std(mockmeans,ddof=1))
            within=float(math.sqrt(np.mean(np.asarray(diffs,float)**2)/2.0))
            ratio=(within/between if between>0 else None)
            if ratio is not None: ratios.append(ratio)
            capout["by_tag_field"][f]={
              "mock_mean_calibrated_diff":describe(mockmeans),
              "R1_minus_R2":describe(diffs),
              "diagnostic_single_rep_random_scale_from_pair_differences":within,
              "between_mock_sd_of_two_rep_mean":between,
              "random_scale_over_between_mock_sd":ratio,
            }
        capout["engineering_gate"]={
          "random_scale_over_between_by_field":{
            f:capout["by_tag_field"][f]["random_scale_over_between_mock_sd"]
            for f in FIELD_NAMES},
          "median_ratio":float(np.median(ratios)),
          "max_ratio":float(np.max(ratios)),
          "median_max_allowed":RANDOM_TO_BETWEEN_MEDIAN_MAX,
          "field_max_allowed":RANDOM_TO_BETWEEN_FIELD_MAX,
        }
        capout["engineering_gate"]["random_integration_subdominant"]=bool(
            capout["engineering_gate"]["median_ratio"]<=RANDOM_TO_BETWEEN_MEDIAN_MAX
            and capout["engineering_gate"]["max_ratio"]<=RANDOM_TO_BETWEEN_FIELD_MAX)
        out[cap]=capout
    return out

def synthetic_self_test():
    seeds=[rep_seed(mid,r,c,t) for mid in IDS for r in REPS for c in CAPS for t in TRACERS]
    need(len(seeds)==len(set(seeds)),"E55 seed collision")
    need(rep_seed(1,"R1","NGC","eBOSS_LRG")==202609315400,
         "E55 mock0001 R1 no longer replays E54")
    need(rep_seed(1,"R2","NGC","eBOSS_LRG")==202609325400,
         "E55 mock0001 R2 no longer replays E54")
    print("E55_SYNTHETIC_SEED_AND_AGGREGATION_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()
    synthetic_self_test()
    if args.self_test: return
    need(args.run,"Use --run for E55 mock-only calculation")

    p,a02,e4=V1.source_only_gate()
    protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        protocol,require_local=True)
    paths,total=A02.resolve_72_sources(
        protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=a02
    need(archive["completed_cases"]==18 and archive["observed_odd_data_vector_read"] is False,
         "A02 parent changed")

    e54=None
    if E54.is_file():
        e54=json.loads(E54.read_text())
        need(e54["status"]=="PASS_RANDOM_REPLICA_FLOOR_QUANTIFIED"
             and e54["decision"]["48k_random_floor_small_all_caps"] is True,
             "Local E54 parent failed")

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state.get("stage")=="E55_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_TWO48K",
             "Existing E55 checkpoint incompatible")
    else:
        state={
          "stage":"E55_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_TWO48K",
          "date":"2026-09-30","status":"INCOMPLETE",
          "mock_ids":list(IDS),"caps":list(CAPS),"replicas":list(REPS),
          "random_size_per_tracer":N_RANDOM,"seed_root":SEED_ROOT,
          "same_E51_fields":True,"same_E19_injection":True,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},"observed_galaxy_rows_used":False,
          "observed_odd_used":False,"physical_velocity_reconstruction_used":False,
          "covariance_inverse_used":False,"pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        case=state["cases"].setdefault(ck,{
          "status":"INCOMPLETE","id":mid,"cap":cap,"replicas":{},
          "full_eligible_rows":{},"random_pool_rows":{}
        })
        if case.get("status")=="complete":
            print("E55_REUSE_CASE",ck,flush=True); continue

        if mid==1 and e54 is not None:
            ec=e54["cases"][ck]
            for rep in REPS:
                expected_seed={tr:rep_seed(mid,rep,cap,tr) for tr in TRACERS}
                need(ec["replicas"][rep]["seeds"]==expected_seed,
                     "E54 seed replay changed "+ck+"/"+rep)
                case["replicas"][rep]={
                  "status":"complete_reused_E54",
                  "seeds":expected_seed,
                  "summary":compact_summary(ec["replicas"][rep]["summary"])
                }
            case["random_pool_rows"]=ec["random_pool_rows"]
            case["fresh_replica_overlap"]=ec["fresh_replica_overlap"]
            case["status"]="complete"
            atomic(OUT,state)
            print("E55_REUSED_E54_MOCK0001",ck,flush=True)
            continue

        archived=archive["cases"][ck]
        prior=next((z for z in old0001["cases"] if z["cap"]==cap),None) if mid==1 else None
        full_data={}; random_pool={}; fullinfo={}
        for tracer in TRACERS:
          k=A02.key(mid,cap,tracer)
          for role in ROLES:
            src=paths[k+"/"+role]
            rows=A02.source_header_rows(src,prior,tracer,role)
            fixed,full,info=load_full_verified(
                src,mid=mid,cap=cap,tracer=tracer,role=role,
                expected_rows=rows,
                archived_info=archived["input_sample_diagnostics"][tracer+"_"+role],
                origproto=origproto)
            if role=="dat":
                full_data[tracer]=full
                case["full_eligible_rows"][tracer]=len(full[0])
            else:
                random_pool[tracer]=full
                case["random_pool_rows"][tracer]=len(full[0])
            fullinfo[tracer+"_"+role]={
              "eligible_rows":len(full[0]),
              "A02_selected_array_SHA256":info["selected_array_SHA256"],
            }
        case["source_replay"]=fullinfo

        repcats={}; reppos={}
        for rep in REPS:
            repcats[rep]={}; reppos[rep]={}
            for tracer in TRACERS:
                cat,pos=sample_random(random_pool[tracer],rep_seed(mid,rep,cap,tracer))
                repcats[rep][tracer]=cat; reppos[rep][tracer]=pos

        dist=exact_distance_cache(
            list(full_data.values())+
            [repcats[r][t] for r in REPS for t in TRACERS])

        if "D1D2" not in case:
            dd=E53.sparse_marked(full_data["eBOSS_LRG"],full_data["eBOSS_ELG"],
                                 dist,chunk=args.chunk,
                                 max_candidates=args.max_candidate_pairs_per_term)
            case["D1D2"]=E53.encode_term(dd)
            atomic(OUT,state)
            print("E55_DD_PASS",ck,"ACC",dd["meta"]["accepted_pairs"],flush=True)

        for rep in REPS:
            rr=case["replicas"].setdefault(rep,{
              "status":"INCOMPLETE",
              "seeds":{t:rep_seed(mid,rep,cap,t) for t in TRACERS},
              "terms":{}
            })
            samples={
              "D1R2":(full_data["eBOSS_LRG"],repcats[rep]["eBOSS_ELG"]),
              "R1D2":(repcats[rep]["eBOSS_LRG"],full_data["eBOSS_ELG"]),
              "R1R2":(repcats[rep]["eBOSS_LRG"],repcats[rep]["eBOSS_ELG"]),
            }
            for term,(a,b) in samples.items():
                if term in rr["terms"]:
                    print("E55_REUSE_TERM",ck,rep,term,flush=True); continue
                z=E53.sparse_marked(a,b,dist,chunk=args.chunk,
                                    max_candidates=args.max_candidate_pairs_per_term)
                rr["terms"][term]=E53.encode_term(z)
                atomic(OUT,state)
                print("E55_TERM_PASS",ck,rep,term,
                      "ACC",z["meta"]["accepted_pairs"],flush=True)
            sm=E53.summarize_level({"D1D2":case["D1D2"],**rr["terms"]})
            rr["summary"]=compact_summary(sm)
            rr["status"]="complete"
            atomic(OUT,state)
            print("E55_REPLICA_PASS",ck,rep,
                  "ANGLE",rr["summary"]["postwindow_angle_rad"],flush=True)

        case["fresh_replica_overlap"]={}
        for tr in TRACERS:
            A=reppos["R1"][tr]; B=reppos["R2"][tr]
            case["fresh_replica_overlap"][tr]=float(
                len(np.intersect1d(A,B,assume_unique=True))/N_RANDOM)
        case["status"]="complete"
        atomic(OUT,state)
        print("E55_CASE_PASS",ck,flush=True)

    need(len(state["cases"])==18 and
         all(c["status"]=="complete" for c in state["cases"].values()),
         "Not all E55 cases complete")
    state["aggregate"]=aggregate(state)
    gates=[state["aggregate"][cap]["engineering_gate"]["random_integration_subdominant"]
           for cap in CAPS]
    state["decision"]={
      "random_integration_subdominant_all_caps":all(gates),
      "interpretation":(
        "TWO48K_AVERAGED_NINEMOCK_BACKGROUND_READY_FOR_NEXT_PHYSICAL_TAG_GATE"
        if all(gates) else
        "RANDOM_INTEGRATION_NOT_SUBDOMINANT; INCREASE_REPLICA_AVERAGING_BEFORE_PHYSICAL_TAG_GATE"
      ),
      "important_limit":"Nine mocks give descriptive 1D background only; no robust tail probability, detection or exclusion."
    }
    state["status"]="PASS_NINEMOCK_FULL_ELIGIBLE_TWO48K_BACKGROUND_QUANTIFIED"
    atomic(OUT,state)

    print("E55_WSL_NINEMOCK_PASS",flush=True)
    for cap in CAPS:
        a=state["aggregate"][cap]
        print("CAP",cap,
              "ANGLE_SD",a["postwindow_angle_mock_means"]["sample_sd"],
              "RANDOM_RATIO_MED",a["engineering_gate"]["median_ratio"],
              "RANDOM_RATIO_MAX",a["engineering_gate"]["max_ratio"],
              "RANDOM_SUBDOM",a["engineering_gate"]["random_integration_subdominant"],
              flush=True)
        for f in FIELD_NAMES:
            z=a["by_tag_field"][f]
            print("TAG",f,
                  "MOCK_SD",z["between_mock_sd_of_two_rep_mean"],
                  "RANDOM_SCALE",z["diagnostic_single_rep_random_scale_from_pair_differences"],
                  "RATIO",z["random_scale_over_between_mock_sd"],
                  "A83",z["mock_mean_calibrated_diff"]["descriptive_abs_threshold_83pct"],
                  "A94",z["mock_mean_calibrated_diff"]["descriptive_abs_threshold_94pct"],
                  flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
