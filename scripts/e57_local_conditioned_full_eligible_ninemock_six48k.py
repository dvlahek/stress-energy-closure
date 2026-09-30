#!/usr/bin/env python3
"""
E57: six-replica 48k random integration across all nine full-eligible mocks.

E56 retained the frozen estimator and four 48k random replicas. Both caps
missed the frozen median precision target 0.25 narrowly, while both per-field
maxima were already below 0.50.

E57 changes ONLY numerical random-integration precision:

  * R1-R4 summaries are reused verbatim from E56;
  * R5/R6 are added with deterministic preregistered seeds;
  * the same all-source-eligible mock galaxies are used;
  * E51 tag fields A-D and E19 F+/F- basis remain unchanged;
  * observed galaxy rows and observed odd vector remain SEALED.

R=6 is selected as the minimum integer implied by standard 1/sqrt(R)
Monte-Carlo scaling from the E56 random-precision diagnostic while KEEPING
the E56 final-mean gate fixed. This is an engineering resource decision, not
a science cut or an observed-data tuning step.

For each cap/tag field and mock i:
  mbar_i = mean over R1..R6 calibrated F--F+ coordinates.

The pooled within-mock single-replica variance is
  s_R^2 = sum_i sum_r (m_ir-mbar_i)^2 / [Nmock*(R-1)].

The effective random-integration error of the six-replica mean is
  sigma_Rbar = s_R/sqrt(6).

Frozen gate, unchanged from E56:
  median_field(sigma_Rbar/s_mock) <= 0.25
  max_field(sigma_Rbar/s_mock)    <= 0.50
in EACH cap.

Replica averaging estimates conditional random-subset Monte-Carlo noise only.
It does not diagnose a common finite-parent-pool discretization bias.

No physical velocity reconstruction, no covariance inverse, no p-values,
no detection/exclusion claim.
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

E54=ROOT/"source_data/e54_conditioned_mock0001_48k_random_replica_floor.json"
E56_RESULT=ROOT/"source_data/e56_conditioned_full_eligible_ninemock_four48k_result.json"
OUT=ROOT/"source_data/e57_conditioned_full_eligible_ninemock_six48k_result.json"

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES
FIELD_NAMES=E53.FIELD_NAMES
REPS=("R1","R2","R3","R4","R5","R6")
OLD_REPS=("R1","R2","R3","R4")
NEW_REPS=("R5","R6")
N_RANDOM=48000
SEED_ROOT=202609305400

MEDIAN_EFFECTIVE_RATIO_MAX=0.25
FIELD_EFFECTIVE_RATIO_MAX=0.50

def need(c,msg):
    if not c:
        raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e57_",delete=False) as f:
        t=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try:
        os.replace(t,path)
    finally:
        t.unlink(missing_ok=True)

def rep_seed(mid,rep,cap,tracer):
    return (SEED_ROOT
            + 1_000_000*IDS.index(mid)
            + 10_000*(REPS.index(rep)+1)
            + 1_000*CAPS.index(cap)
            + 100*TRACERS.index(tracer))

def compact_from_any(summary):
    z=summary["by_tag_field"][FIELD_NAMES[0]]
    if "calibrated_diff" in z:
        return {
          "postwindow_angle_rad":float(summary["postwindow_angle_rad"]),
          "difference_norm_over_plus_norm":float(summary["difference_norm_over_plus_norm"]),
          "shape_residual_norm_over_plus_norm":float(summary["shape_residual_norm_over_plus_norm"]),
          "rho_minus_on_plus":float(summary["rho_minus_on_plus"]),
          "by_tag_field":{
            f:{
              "baseline_L2":float(summary["by_tag_field"][f]["baseline_L2"]),
              "calibrated_diff":float(summary["by_tag_field"][f]["calibrated_diff"]),
              "plus_lambda":float(summary["by_tag_field"][f]["plus_lambda"]),
              "minus_lambda":float(summary["by_tag_field"][f]["minus_lambda"]),
              "shape_only":float(summary["by_tag_field"][f]["shape_only"]),
            } for f in FIELD_NAMES
          }
        }
    return E55.compact_summary(summary)

def describe(v):
    a=np.asarray(v,float); med=float(np.median(a)); av=np.sort(np.abs(a))
    return {
      "n":int(len(a)),
      "mean":float(np.mean(a)),
      "median":med,
      "sample_sd":float(np.std(a,ddof=1)),
      "MAD":float(np.median(np.abs(a-med))),
      "min":float(np.min(a)),
      "max":float(np.max(a)),
      "descriptive_abs_threshold_83pct":float(av[5]),
      "descriptive_abs_threshold_94pct":float(av[7]),
      "values":[float(x) for x in a],
    }

def pooled_random_scale(replica_values_by_mock):
    a=np.asarray(replica_values_by_mock,float)
    need(a.shape==(len(IDS),len(REPS)),"Unexpected E57 replica matrix shape")
    means=np.mean(a,axis=1,keepdims=True)
    sse=float(np.sum((a-means)**2))
    dof=len(IDS)*(len(REPS)-1)
    single=math.sqrt(sse/dof)
    effective=single/math.sqrt(len(REPS))
    return single,effective

def aggregate(state):
    out={}
    for cap in CAPS:
        cases=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        capout={"by_tag_field":{}}
        angle_mockmeans=[]; diff_mockmeans=[]
        for c in cases:
            angle_mockmeans.append(float(np.mean([
                c["replicas"][r]["summary"]["postwindow_angle_rad"] for r in REPS])))
            diff_mockmeans.append(float(np.mean([
                c["replicas"][r]["summary"]["difference_norm_over_plus_norm"] for r in REPS])))
        capout["postwindow_angle_sixrep_mock_means"]=describe(angle_mockmeans)
        capout["difference_norm_over_plus_sixrep_mock_means"]=describe(diff_mockmeans)

        ratios=[]
        for f in FIELD_NAMES:
            matrix=[]; mockmeans=[]
            for c in cases:
                vals=[float(c["replicas"][r]["summary"]["by_tag_field"][f]["calibrated_diff"])
                      for r in REPS]
                matrix.append(vals); mockmeans.append(float(np.mean(vals)))
            single,effective=pooled_random_scale(matrix)
            between=float(np.std(np.asarray(mockmeans,float),ddof=1))
            ratio=(effective/between if between>0 else None)
            need(ratio is not None and np.isfinite(ratio),
                 "Invalid effective-random/between-mock ratio")
            ratios.append(ratio)
            intrinsic_var=max(0.0,between*between-effective*effective)
            capout["by_tag_field"][f]={
              "sixrep_mock_mean_calibrated_diff":describe(mockmeans),
              "pooled_single_rep_random_sd":float(single),
              "effective_random_se_sixrep_mean":float(effective),
              "between_mock_sd_sixrep_mean":between,
              "effective_random_over_between_mock_sd":float(ratio),
              "descriptive_random_deconvolved_between_mock_sd":float(math.sqrt(intrinsic_var)),
              "replica_matrix":[[float(v) for v in row] for row in matrix],
            }

        capout["engineering_gate"]={
          "effective_random_over_between_by_field":{
             f:capout["by_tag_field"][f]["effective_random_over_between_mock_sd"]
             for f in FIELD_NAMES},
          "median_ratio":float(np.median(ratios)),
          "max_ratio":float(np.max(ratios)),
          "median_max_allowed":MEDIAN_EFFECTIVE_RATIO_MAX,
          "field_max_allowed":FIELD_EFFECTIVE_RATIO_MAX,
        }
        capout["engineering_gate"]["sixrep_random_mean_subdominant"]=bool(
            capout["engineering_gate"]["median_ratio"]<=MEDIAN_EFFECTIVE_RATIO_MAX
            and capout["engineering_gate"]["max_ratio"]<=FIELD_EFFECTIVE_RATIO_MAX)
        out[cap]=capout
    return out

def synthetic_self_test():
    seeds=[rep_seed(mid,r,c,t)
           for mid in IDS for r in REPS for c in CAPS for t in TRACERS]
    need(len(seeds)==len(set(seeds)),"E57 seed collision")
    expected={
      ("R5","NGC","eBOSS_LRG"):202609355400,
      ("R6","NGC","eBOSS_LRG"):202609365400,
      ("R5","SGC","eBOSS_ELG"):202609356500,
      ("R6","SGC","eBOSS_ELG"):202609366500,
    }
    for k,v in expected.items():
        need(rep_seed(1,*k)==v,"E57 mock0001 seed replay failed")
    z=np.arange(54,dtype=float).reshape(9,6)
    s,e=pooled_random_scale(z)
    need(np.isfinite(s) and abs(e-s/math.sqrt(6))<1e-15,
         "E57 sqrt(6) scaling failed")
    print("E57_SYNTHETIC_SIX_REPLICA_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()

    synthetic_self_test()
    if args.self_test:
        return
    need(args.run,"Use --run for E57 mock-only calculation")
    need(E54.is_file() and E56_RESULT.is_file(),
         "E57 requires completed local E54 and E56 result JSONs")

    e54=json.loads(E54.read_text())
    e56=json.loads(E56_RESULT.read_text())
    need(e54["status"]=="PASS_RANDOM_REPLICA_FLOOR_QUANTIFIED","E54 parent incomplete")
    need(e56["status"]=="PASS_NINEMOCK_FULL_ELIGIBLE_FOUR48K_BACKGROUND_QUANTIFIED",
         "E56 parent incomplete")
    need(e56["decision"]["fourrep_random_mean_subdominant_all_caps"] is False,
         "E56 decision unexpectedly changed")
    for parent in (e54,e56):
        need(parent["observed_odd_used"] is False
             and parent["observed_galaxy_rows_used"] is False,
             "Parent observation guardrail changed")

    p,a02,e4=V1.source_only_gate()
    protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        protocol,require_local=True)
    paths,total=A02.resolve_72_sources(
        protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=a02
    need(archive["completed_cases"]==18
         and archive["observed_odd_data_vector_read"] is False,
         "A02 parent changed")

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state.get("stage")=="E57_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_SIX48K",
             "Existing E57 checkpoint incompatible")
    else:
        state={
          "stage":"E57_CONDITIONED_FULL_ELIGIBLE_NINEMOCK_SIX48K",
          "date":"2026-09-30","status":"INCOMPLETE",
          "mock_ids":list(IDS),"caps":list(CAPS),"replicas":list(REPS),
          "new_replicas":list(NEW_REPS),
          "random_size_per_tracer":N_RANDOM,"seed_root":SEED_ROOT,
          "E56_R1_R4_reused":True,
          "same_E51_fields":True,"same_E19_injection":True,
          "gate_carried_forward_unchanged_from_E56":True,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,
          "observed_odd_used":False,
          "physical_velocity_reconstruction_used":False,
          "covariance_inverse_used":False,
          "pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        case=state["cases"].setdefault(ck,{
          "status":"INCOMPLETE","id":mid,"cap":cap,"replicas":{}
        })
        if case.get("status")=="complete":
            print("E57_REUSE_CASE",ck,flush=True); continue

        pc=e56["cases"][ck]
        need(pc["status"]=="complete","E56 case incomplete "+ck)
        for rep in OLD_REPS:
            expected={tr:rep_seed(mid,rep,cap,tr) for tr in TRACERS}
            need(pc["replicas"][rep]["seeds"]==expected,
                 "E56 seed mismatch "+ck+"/"+rep)
            case["replicas"][rep]={
              "status":"complete_reused_E56",
              "seeds":expected,
              "summary":compact_from_any(pc["replicas"][rep]["summary"])
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

        repcats={}; reppos={}
        for rep in NEW_REPS:
            repcats[rep]={}; reppos[rep]={}
            for tracer in TRACERS:
                cat,pos=E55.sample_random(
                    random_pool[tracer],rep_seed(mid,rep,cap,tracer))
                repcats[rep][tracer]=cat; reppos[rep][tracer]=pos

        dist=E55.exact_distance_cache(
            list(full_data.values())+
            [repcats[r][t] for r in NEW_REPS for t in TRACERS])

        if "D1D2" in pc:
            case["D1D2"]=pc["D1D2"]
        elif "D1D2" not in case:
            dd=E53.sparse_marked(
                full_data["eBOSS_LRG"],full_data["eBOSS_ELG"],dist,
                chunk=args.chunk,
                max_candidates=args.max_candidate_pairs_per_term)
            case["D1D2"]=E53.encode_term(dd)
            e54r1=e54["cases"][ck]["replicas"]["R1"]
            replay=compact_from_any(E53.summarize_level(
                {"D1D2":case["D1D2"],**e54r1["terms"]}))
            archived_summary=compact_from_any(e54r1["summary"])
            need(abs(replay["postwindow_angle_rad"]-
                     archived_summary["postwindow_angle_rad"])<5e-12,
                 "mock0001 DD replay angle mismatch "+ck)
            for f in FIELD_NAMES:
                need(abs(replay["by_tag_field"][f]["calibrated_diff"]-
                         archived_summary["by_tag_field"][f]["calibrated_diff"])<5e-12,
                     "mock0001 DD replay calibrated mismatch "+ck+"/"+f)
            atomic(OUT,state)
            print("E57_MOCK0001_DD_REPLAY_PASS",ck,flush=True)

        for rep in NEW_REPS:
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
                    print("E57_REUSE_TERM",ck,rep,term,flush=True); continue
                z=E53.sparse_marked(
                    a,b,dist,chunk=args.chunk,
                    max_candidates=args.max_candidate_pairs_per_term)
                rr["terms"][term]=E53.encode_term(z)
                atomic(OUT,state)
                print("E57_TERM_PASS",ck,rep,term,
                      "ACC",z["meta"]["accepted_pairs"],flush=True)
            sm=E53.summarize_level({"D1D2":case["D1D2"],**rr["terms"]})
            rr["summary"]=compact_from_any(sm)
            rr["status"]="complete"
            atomic(OUT,state)
            print("E57_REPLICA_PASS",ck,rep,
                  "ANGLE",rr["summary"]["postwindow_angle_rad"],flush=True)

        case["R5_R6_overlap"]={}
        for tr in TRACERS:
            case["R5_R6_overlap"][tr]=float(
                len(np.intersect1d(reppos["R5"][tr],reppos["R6"][tr],
                                   assume_unique=True))/N_RANDOM)
        case["status"]="complete"
        atomic(OUT,state)
        print("E57_CASE_PASS",ck,flush=True)

    need(len(state["cases"])==18
         and all(c["status"]=="complete" for c in state["cases"].values()),
         "Not all E57 cases complete")
    state["aggregate"]=aggregate(state)
    gates=[state["aggregate"][cap]["engineering_gate"]["sixrep_random_mean_subdominant"]
           for cap in CAPS]
    state["decision"]={
      "sixrep_random_mean_subdominant_all_caps":bool(all(gates)),
      "interpretation":(
        "SIXREP_RANDOM_INTEGRATION_GATE_CLOSED; ADVANCE_TO_PHYSICAL_VELOCITY_TAG_AND_ABSOLUTE_AMPLITUDE_CALIBRATION"
        if all(gates) else
        "SIXREP_RANDOM_INTEGRATION_STILL_MATERIAL; DO_NOT_ADVANCE_PHYSICAL_TAG_GATE"
      ),
      "important_limit":"Passing would close only conditional random-subset Monte-Carlo precision. It would not calibrate a physical velocity tag, an absolute eBOSS signal, finite-parent-pool common bias, or statistical significance."
    }
    state["status"]="PASS_NINEMOCK_FULL_ELIGIBLE_SIX48K_BACKGROUND_QUANTIFIED"
    atomic(OUT,state)

    print("E57_WSL_SIXREP_NINEMOCK_PASS",flush=True)
    for cap in CAPS:
        a=state["aggregate"][cap]; g=a["engineering_gate"]
        print("CAP",cap,
              "EFF_RATIO_MED",g["median_ratio"],
              "EFF_RATIO_MAX",g["max_ratio"],
              "SIXREP_RANDOM_SUBDOM",g["sixrep_random_mean_subdominant"],
              flush=True)
        for f in FIELD_NAMES:
            z=a["by_tag_field"][f]
            print("TAG",f,
                  "MOCK_SD",z["between_mock_sd_sixrep_mean"],
                  "SINGLE_RAND_SD",z["pooled_single_rep_random_sd"],
                  "SIXREP_RAND_SE",z["effective_random_se_sixrep_mean"],
                  "RATIO",z["effective_random_over_between_mock_sd"],
                  "A83",z["sixrep_mock_mean_calibrated_diff"]["descriptive_abs_threshold_83pct"],
                  "A94",z["sixrep_mock_mean_calibrated_diff"]["descriptive_abs_threshold_94pct"],
                  flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
