#!/usr/bin/env python3
"""
E58: seven-replica 48k random-integration closure.

E57 passed SGC but NGC missed the frozen median final-mean precision gate by
0.001413 (0.251413 > 0.25). E58 changes only numerical integration precision:
reuse R1-R6 and DD from E57, add one deterministic R7 per mock/cap/tracer,
keep the E56/E57 gate unchanged.

R=7 is the minimum integer implied by the E57 NGC diagnostic under standard
1/sqrt(R) Monte-Carlo scaling:
    ceil(6 * (0.2514132726578298 / 0.25)^2) = 7.

No observed galaxy rows, no observed odd vector, no physical velocity
reconstruction, no covariance inverse, no p-values or detection claim.
"""
from __future__ import annotations
import argparse, json, math, os, sys, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import run_eboss_dr16_full_eligible_mock0001_sparse as V1
import e53_local_conditioned_mock0001_sample_scaling as E53
import e55_local_conditioned_full_eligible_ninemock_two48k as E55
import e57_local_conditioned_full_eligible_ninemock_six48k as E57

PARENT=ROOT/"source_data/e57_conditioned_full_eligible_ninemock_six48k_result.json"
OUT=ROOT/"source_data/e58_conditioned_full_eligible_ninemock_seven48k_result.json"

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES
FIELD_NAMES=E53.FIELD_NAMES
REPS=("R1","R2","R3","R4","R5","R6","R7")
OLD_REPS=REPS[:-1]
NEW_REP="R7"
N_RANDOM=48000
SEED_ROOT=202609305400
MEDIAN_MAX=0.25
FIELD_MAX=0.50

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e58_",delete=False) as f:
        t=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(t,path)
    finally: t.unlink(missing_ok=True)

def rep_seed(mid,rep,cap,tracer):
    return (SEED_ROOT
            + 1_000_000*IDS.index(mid)
            + 10_000*(REPS.index(rep)+1)
            + 1_000*CAPS.index(cap)
            + 100*TRACERS.index(tracer))

def describe(v):
    a=np.asarray(v,float); med=float(np.median(a)); av=np.sort(np.abs(a))
    return {
      "n":len(a),"mean":float(np.mean(a)),"median":med,
      "sample_sd":float(np.std(a,ddof=1)),
      "MAD":float(np.median(np.abs(a-med))),
      "min":float(np.min(a)),"max":float(np.max(a)),
      "descriptive_abs_threshold_83pct":float(av[5]),
      "descriptive_abs_threshold_94pct":float(av[7]),
      "values":[float(x) for x in a],
    }

def pooled(matrix):
    a=np.asarray(matrix,float)
    need(a.shape==(len(IDS),len(REPS)),"Unexpected E58 matrix shape")
    m=np.mean(a,axis=1,keepdims=True)
    single=math.sqrt(float(np.sum((a-m)**2))/(len(IDS)*(len(REPS)-1)))
    return single,single/math.sqrt(len(REPS))

def aggregate(state):
    out={}
    for cap in CAPS:
        cases=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        z={"by_tag_field":{}}
        angles=[float(np.mean([c["replicas"][r]["summary"]["postwindow_angle_rad"]
                              for r in REPS])) for c in cases]
        diffs=[float(np.mean([c["replicas"][r]["summary"]["difference_norm_over_plus_norm"]
                             for r in REPS])) for c in cases]
        z["postwindow_angle_sevenrep_mock_means"]=describe(angles)
        z["difference_norm_over_plus_sevenrep_mock_means"]=describe(diffs)
        ratios=[]
        for f in FIELD_NAMES:
            matrix=[]; means=[]
            for c in cases:
                vals=[float(c["replicas"][r]["summary"]["by_tag_field"][f]["calibrated_diff"])
                      for r in REPS]
                matrix.append(vals); means.append(float(np.mean(vals)))
            single,eff=pooled(matrix)
            between=float(np.std(means,ddof=1))
            ratio=eff/between
            ratios.append(ratio)
            z["by_tag_field"][f]={
              "sevenrep_mock_mean_calibrated_diff":describe(means),
              "pooled_single_rep_random_sd":float(single),
              "effective_random_se_sevenrep_mean":float(eff),
              "between_mock_sd_sevenrep_mean":between,
              "effective_random_over_between_mock_sd":float(ratio),
              "descriptive_random_deconvolved_between_mock_sd":
                  float(math.sqrt(max(0.0,between*between-eff*eff))),
              "replica_matrix":[[float(x) for x in row] for row in matrix],
            }
        med=float(np.median(ratios)); mx=float(np.max(ratios))
        z["engineering_gate"]={
          "effective_random_over_between_by_field":{
            f:z["by_tag_field"][f]["effective_random_over_between_mock_sd"]
            for f in FIELD_NAMES},
          "median_ratio":med,"max_ratio":mx,
          "median_max_allowed":MEDIAN_MAX,"field_max_allowed":FIELD_MAX,
          "sevenrep_random_mean_subdominant":bool(med<=MEDIAN_MAX and mx<=FIELD_MAX)
        }
        out[cap]=z
    return out

def self_test():
    seeds=[rep_seed(mid,r,c,t) for mid in IDS for r in REPS for c in CAPS for t in TRACERS]
    need(len(seeds)==len(set(seeds)),"E58 seed collision")
    need(rep_seed(1,"R7","NGC","eBOSS_LRG")==202609375400,"R7 NGC LRG seed mismatch")
    need(rep_seed(1,"R7","SGC","eBOSS_ELG")==202609376500,"R7 SGC ELG seed mismatch")
    required=math.ceil(6*(0.2514132726578298/0.25)**2)
    need(required==7,"E58 minimum-R derivation changed")
    print("E58_SYNTHETIC_SEVEN_REPLICA_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()
    self_test()
    if args.self_test: return
    need(args.run,"Use --run for E58")
    need(PARENT.is_file(),"Missing E57 result")
    parent=json.loads(PARENT.read_text())
    need(parent["status"]=="PASS_NINEMOCK_FULL_ELIGIBLE_SIX48K_BACKGROUND_QUANTIFIED",
         "E57 parent incomplete")
    need(parent["decision"]["sixrep_random_mean_subdominant_all_caps"] is False,
         "E57 decision unexpectedly changed")
    need(parent["aggregate"]["NGC"]["engineering_gate"]["sixrep_random_mean_subdominant"] is False,
         "E57 NGC parent unexpectedly passed")
    need(parent["aggregate"]["SGC"]["engineering_gate"]["sixrep_random_mean_subdominant"] is True,
         "E57 SGC parent unexpectedly failed")
    need(parent["observed_odd_used"] is False and parent["observed_galaxy_rows_used"] is False,
         "Parent observation guardrail changed")

    p,a02,e4=V1.source_only_gate()
    protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        protocol,require_local=True)
    paths,total=A02.resolve_72_sources(protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=a02
    need(archive["completed_cases"]==18 and archive["observed_odd_data_vector_read"] is False,
         "A02 parent changed")

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state["stage"]=="E58_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_SEVEN48K",
             "Existing E58 checkpoint incompatible")
    else:
        state={
          "stage":"E58_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_SEVEN48K",
          "date":"2026-09-30","status":"INCOMPLETE",
          "mock_ids":list(IDS),"caps":list(CAPS),"replicas":list(REPS),
          "new_replicas":[NEW_REP],"random_size_per_tracer":N_RANDOM,
          "seed_root":SEED_ROOT,"E57_R1_R6_reused":True,
          "gate_carried_forward_unchanged":True,
          "same_E51_fields":True,"same_E19_injection":True,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,"observed_odd_used":False,
          "physical_velocity_reconstruction_used":False,
          "covariance_inverse_used":False,"pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        case=state["cases"].setdefault(ck,{"status":"INCOMPLETE","id":mid,"cap":cap,"replicas":{}})
        if case.get("status")=="complete":
            print("E58_REUSE_CASE",ck,flush=True); continue

        pc=parent["cases"][ck]
        need(pc["status"]=="complete","E57 case incomplete "+ck)
        need("D1D2" in pc,"E57 D1D2 missing "+ck)
        case["D1D2"]=pc["D1D2"]
        for rep in OLD_REPS:
            expected={tr:rep_seed(mid,rep,cap,tr) for tr in TRACERS}
            need(pc["replicas"][rep]["seeds"]==expected,"E57 seed mismatch "+ck+"/"+rep)
            case["replicas"][rep]={
              "status":"complete_reused_E57","seeds":expected,
              "summary":E57.compact_from_any(pc["replicas"][rep]["summary"])
            }

        archived=archive["cases"][ck]
        full_data={}; random_pool={}
        prior=None
        for tracer in TRACERS:
          k=A02.key(mid,cap,tracer)
          for role in ROLES:
            src=paths[k+"/"+role]
            rows=A02.source_header_rows(src,prior,tracer,role)
            fixed,full,info=E55.load_full_verified(
                src,mid=mid,cap=cap,tracer=tracer,role=role,
                expected_rows=rows,
                archived_info=archived["input_sample_diagnostics"][tracer+"_"+role],
                origproto=origproto)
            if role=="dat": full_data[tracer]=full
            else: random_pool[tracer]=full

        repcats={}
        for tracer in TRACERS:
            cat,pos=E55.sample_random(random_pool[tracer],rep_seed(mid,NEW_REP,cap,tracer))
            repcats[tracer]=cat
        dist=E55.exact_distance_cache(list(full_data.values())+list(repcats.values()))

        rr=case["replicas"].setdefault(NEW_REP,{
          "status":"INCOMPLETE",
          "seeds":{t:rep_seed(mid,NEW_REP,cap,t) for t in TRACERS},
          "terms":{}
        })
        samples={
          "D1R2":(full_data["eBOSS_LRG"],repcats["eBOSS_ELG"]),
          "R1D2":(repcats["eBOSS_LRG"],full_data["eBOSS_ELG"]),
          "R1R2":(repcats["eBOSS_LRG"],repcats["eBOSS_ELG"]),
        }
        for term,(a,b) in samples.items():
            if term in rr["terms"]:
                print("E58_REUSE_TERM",ck,term,flush=True); continue
            z=E53.sparse_marked(a,b,dist,chunk=args.chunk,
                                max_candidates=args.max_candidate_pairs_per_term)
            rr["terms"][term]=E53.encode_term(z)
            atomic(OUT,state)
            print("E58_TERM_PASS",ck,term,"ACC",z["meta"]["accepted_pairs"],flush=True)

        sm=E53.summarize_level({"D1D2":case["D1D2"],**rr["terms"]})
        rr["summary"]=E57.compact_from_any(sm)
        rr["status"]="complete"
        case["status"]="complete"
        atomic(OUT,state)
        print("E58_CASE_PASS",ck,"ANGLE",rr["summary"]["postwindow_angle_rad"],flush=True)

    need(len(state["cases"])==18 and all(c["status"]=="complete" for c in state["cases"].values()),
         "Not all E58 cases complete")
    state["aggregate"]=aggregate(state)
    gates=[state["aggregate"][cap]["engineering_gate"]["sevenrep_random_mean_subdominant"]
           for cap in CAPS]
    state["decision"]={
      "sevenrep_random_mean_subdominant_all_caps":bool(all(gates)),
      "interpretation":(
        "SEVENREP_RANDOM_INTEGRATION_GATE_CLOSED; ADVANCE_TO_PHYSICAL_VELOCITY_TAG_AND_ABSOLUTE_AMPLITUDE_CALIBRATION"
        if all(gates) else
        "SEVENREP_RANDOM_INTEGRATION_STILL_MATERIAL; STOP_REPLICA_CHASING_AND_USE_DETERMINISTIC_OR_HIGHER_DENSITY_RANDOM_INTEGRATION"
      ),
      "important_limit":"Passing closes only conditional random-subset Monte-Carlo precision. It does not calibrate a physical velocity tag, absolute eBOSS EV amplitude, finite-parent-pool common bias, covariance inference, or significance."
    }
    state["status"]="PASS_NINEMOCK_FULL_ELIGIBLE_SEVEN48K_BACKGROUND_QUANTIFIED"
    atomic(OUT,state)

    print("E58_WSL_SEVENREP_NINEMOCK_PASS",flush=True)
    for cap in CAPS:
        g=state["aggregate"][cap]["engineering_gate"]
        print("CAP",cap,"EFF_RATIO_MED",g["median_ratio"],"EFF_RATIO_MAX",g["max_ratio"],
              "SEVENREP_RANDOM_SUBDOM",g["sevenrep_random_mean_subdominant"],flush=True)
        for f in FIELD_NAMES:
            z=state["aggregate"][cap]["by_tag_field"][f]
            print("TAG",f,"MOCK_SD",z["between_mock_sd_sevenrep_mean"],
                  "SEVENREP_RAND_SE",z["effective_random_se_sevenrep_mean"],
                  "RATIO",z["effective_random_over_between_mock_sd"],flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
