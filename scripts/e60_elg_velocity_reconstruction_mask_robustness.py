#!/usr/bin/env python3
"""
E60: ELG velocity-reconstruction random-mask robustness, MOCK ONLY.

Parent E59 selected the eBOSS ELG density reconstruction by a prospectively
frozen split-half screen. E59 used one fixed E58 R1 48k random selection field
and explicitly left reconstruction-mask robustness open.

E60 closes that numerical/survey-selection issue before any truth-velocity
calibration:
  * same 9 fixed EZmocks and NGC/SGC caps;
  * same full eligible mock galaxies and same E59 512 DD midpoint probes;
  * ELG reconstruction only (selected by E59, no post-result tracer search);
  * exact E58 R1--R7 48k ELG random-subset seed contract;
  * same R=16 Mpc/h top-hat-smoothed real-space linear velocity kernel and
    256 Mpc/h primary cutoff;
  * compare all 21 random-replica pairs and leave-one-out consensus fields;
  * quantify replica RMS relative to the seven-replica mean field;
  * on mock0001 in each cap, additionally compare the seven-replica mean to
    the FULL eligible ELG random pool.

No observed galaxies, no observed odd vector, no true velocities, no inverse
covariance, no p-values, no detection significance, no F-state retuning.

If E60 passes, the next stage is truth-labelled N-body/lightcone calibration.
If E60 fails, do not proceed to truth/observations; replace the reconstruction
selection integral with full-pool or a deterministic denser scheme.
"""
from __future__ import annotations

import argparse, json, math, os, sys, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import run_eboss_dr16_full_eligible_mock0001_sparse as V1
import e55_local_conditioned_full_eligible_ninemock_two48k as E55
import e59_mock_velocity_tag_repeatability as E59

E58=ROOT/"source_data/e58_conditioned_full_eligible_ninemock_seven48k_result.json"
E59_RESULT=ROOT/"source_data/e59_mock_velocity_tag_repeatability.json"
OUT=ROOT/"source_data/e60_elg_velocity_reconstruction_mask_robustness.json"

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES
SELECTED_TRACER="eBOSS_ELG"
REPS=("R1","R2","R3","R4","R5","R6","R7")
N_RANDOM=48000
E58_SEED_ROOT=202609305400

CAP_MEDIAN_PAIRWISE_PEARSON_MIN=0.95
CAP_MEDIAN_PAIRWISE_SIGN_MIN=0.90
CAP_MEDIAN_RANDOM_RMS_RATIO_MAX=0.25
CAP_MAX_RANDOM_RMS_RATIO_MAX=0.50

FULLPOOL_PEARSON_MIN=0.98
FULLPOOL_SIGN_MIN=0.95
FULLPOOL_RMS_DIFF_RATIO_MAX=0.25

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e60_",delete=False) as f:
        q=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(q,path)
    finally: q.unlink(missing_ok=True)

def rep_seed(mid,rep,cap,tracer=SELECTED_TRACER):
    return (E58_SEED_ROOT
            +1_000_000*IDS.index(mid)
            +10_000*(REPS.index(rep)+1)
            +1_000*CAPS.index(cap)
            +100*TRACERS.index(tracer))

def sample_random(pool,seed):
    need(len(pool[0])>=N_RANDOM,"Eligible random pool <48k")
    pos=np.sort(np.random.default_rng(seed).choice(
        len(pool[0]),size=N_RANDOM,replace=False))
    cat=tuple(np.asarray(v[pos],dtype="f8").copy() for v in pool)
    return cat,pos

def pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a-np.mean(a); b=b-np.mean(b)
    den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(a@b/den) if den>0 else 0.0

def sign_agreement(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    nz=(a!=0)&(b!=0)
    need(np.count_nonzero(nz)>=0.95*len(a),"Too many exact-zero LOS values")
    return float(np.mean(np.sign(a[nz])==np.sign(b[nz])))

def describe(v):
    a=np.asarray(v,float)
    need(np.isfinite(a).all() and len(a)>0,"Bad E60 metric vector")
    return {
      "n":len(a),"mean":float(np.mean(a)),"median":float(np.median(a)),
      "min":float(np.min(a)),"max":float(np.max(a)),
      "sample_sd":float(np.std(a,ddof=1)) if len(a)>1 else 0.0,
    }

def full_data_kernel(cat,probes,chunk):
    x,u,w=E59.cat_xyz(cat)
    from scipy.spatial import cKDTree
    f=E59.weighted_full_field(cKDTree(x,leafsize=32),x,w,probes,
                              E59.RMAX_PRIMARY,chunk)
    sw=float(np.sum(w))
    need(sw>0,"Bad data weight sum")
    return {k:f[k]/sw for k in ("primary","check")},sw

def random_kernel(cat,probes,chunk):
    x,u,w=E59.cat_xyz(cat)
    from scipy.spatial import cKDTree
    f=E59.weighted_full_field(cKDTree(x,leafsize=32),x,w,probes,
                              E59.RMAX_PRIMARY,chunk)
    sw=float(np.sum(w))
    need(sw>0,"Bad random weight sum")
    return {k:f[k]/sw for k in ("primary","check")},sw

def reconstruction_los(D,R,probes):
    return E59.los(D["primary"]-R["primary"],probes)

def case_metrics(rep_los):
    A=np.stack([np.asarray(rep_los[r],float) for r in REPS],axis=0)
    pair_r=[]; pair_s=[]
    for i in range(len(REPS)):
      for j in range(i+1,len(REPS)):
        pair_r.append(pearson(A[i],A[j]))
        pair_s.append(sign_agreement(A[i],A[j]))
    mean=np.mean(A,axis=0)
    loo_r=[]; loo_s=[]
    for i in range(len(REPS)):
        loo=(np.sum(A,axis=0)-A[i])/(len(REPS)-1)
        loo_r.append(pearson(A[i],loo))
        loo_s.append(sign_agreement(A[i],loo))
    random_rms=float(np.sqrt(np.mean(np.var(A,axis=0,ddof=1))))
    mean_rms=float(np.sqrt(np.mean(mean*mean)))
    ratio=random_rms/mean_rms if mean_rms>0 else math.inf
    return {
      "pairwise_pearson":describe(pair_r),
      "pairwise_sign_agreement":describe(pair_s),
      "replica_vs_leave_one_out_pearson":describe(loo_r),
      "replica_vs_leave_one_out_sign_agreement":describe(loo_s),
      "random_replica_rms":random_rms,
      "sevenrep_mean_field_rms":mean_rms,
      "random_rms_over_sevenrep_mean_rms":ratio,
      "sevenrep_mean_los":mean.tolist(),
    }

def aggregate(state):
    out={}
    for cap in CAPS:
        cs=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        pr=[c["metrics"]["pairwise_pearson"]["median"] for c in cs]
        ps=[c["metrics"]["pairwise_sign_agreement"]["median"] for c in cs]
        rr=[c["metrics"]["random_rms_over_sevenrep_mean_rms"] for c in cs]
        out[cap]={
          "case_median_pairwise_pearson":describe(pr),
          "case_median_pairwise_sign_agreement":describe(ps),
          "random_rms_over_mean_field_rms":describe(rr),
        }
        out[cap]["gate"]={
          "median_pairwise_pearson_min":CAP_MEDIAN_PAIRWISE_PEARSON_MIN,
          "median_pairwise_sign_min":CAP_MEDIAN_PAIRWISE_SIGN_MIN,
          "median_random_rms_ratio_max":CAP_MEDIAN_RANDOM_RMS_RATIO_MAX,
          "max_random_rms_ratio_max":CAP_MAX_RANDOM_RMS_RATIO_MAX,
          "pass":bool(
            out[cap]["case_median_pairwise_pearson"]["median"]>=CAP_MEDIAN_PAIRWISE_PEARSON_MIN
            and out[cap]["case_median_pairwise_sign_agreement"]["median"]>=CAP_MEDIAN_PAIRWISE_SIGN_MIN
            and out[cap]["random_rms_over_mean_field_rms"]["median"]<=CAP_MEDIAN_RANDOM_RMS_RATIO_MAX
            and out[cap]["random_rms_over_mean_field_rms"]["max"]<=CAP_MAX_RANDOM_RMS_RATIO_MAX
          )
        }
    return out

def self_test():
    seeds=[rep_seed(mid,r,c) for mid in IDS for r in REPS for c in CAPS]
    need(len(seeds)==len(set(seeds)),"E60 E58-seed collision")
    need(rep_seed(1,"R1","NGC")==202609315500,"E60 R1 NGC ELG seed drift")
    need(rep_seed(1,"R7","SGC")==202609376500,"E60 R7 SGC ELG seed drift")
    x=np.array([1.,-1.,2.,-3.]); y=np.array([2.,-2.,1.,-4.])
    need(sign_agreement(x,y)==1.0 and pearson(x,y)>0.9,"E60 metric self-test fail")
    print("E60_SYNTHETIC_MASK_ROBUSTNESS_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    args=ap.parse_args()
    self_test()
    if args.self_test:return
    need(args.run,"Use --run for E60 mock-only calculation")
    need(16<=args.chunk<=256,"Unsafe E60 chunk")

    need(E58.is_file() and E59_RESULT.is_file(),"Missing local E58/E59 parent")
    e58=json.loads(E58.read_text()); e59=json.loads(E59_RESULT.read_text())
    need(e58["status"]=="PASS_NINEMOCK_FULL_ELIGIBLE_SEVEN48K_BACKGROUND_QUANTIFIED"
         and e58["decision"]["sevenrep_random_mean_subdominant_all_caps"] is True,
         "E58 parent gate changed")
    need(e59["status"]=="PASS_MOCK_VELOCITY_TAG_REPEATABILITY_QUANTIFIED"
         and e59["decision"]["selected_tracer_reconstruction_for_truth_calibration"]==SELECTED_TRACER
         and e59["candidate_reconstructions"][SELECTED_TRACER]["repeatability_screen_pass"] is True,
         "E59 selected-tracer parent changed")
    need(e59["observed_odd_used"] is False and e59["observed_galaxy_rows_used"] is False
         and e59["physical_velocity_truth_correlation_calibrated"] is False,
         "E59 guardrail changed")

    p,a02,e4=V1.source_only_gate()
    protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        protocol,require_local=True)
    paths,total=A02.resolve_72_sources(protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=a02
    need(archive["completed_cases"]==18 and archive["observed_odd_data_vector_read"] is False,
         "A02 source parent changed")

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state.get("stage")=="E60_ELG_VELOCITY_RECONSTRUCTION_MASK_ROBUSTNESS",
             "Existing E60 checkpoint incompatible")
    else:
        state={
          "stage":"E60_ELG_VELOCITY_RECONSTRUCTION_MASK_ROBUSTNESS",
          "date":"2026-09-30","status":"INCOMPLETE",
          "selected_tracer_fixed_from_E59":SELECTED_TRACER,
          "mock_ids":list(IDS),"caps":list(CAPS),"replicas":list(REPS),
          "random_size_per_replica":N_RANDOM,
          "same_E59_kernel":True,"same_E59_probe_policy":True,
          "full_pool_anchor_mock_id":1,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,"observed_odd_used":False,
          "true_velocity_used":False,"absolute_EV_amplitude_calibrated":False,
          "covariance_inverse_used":False,"pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        case=state["cases"].setdefault(ck,{
          "id":mid,"cap":cap,"status":"INCOMPLETE","replicas":{}
        })
        if case.get("status")=="complete":
            print("E60_REUSE_CASE",ck,flush=True); continue

        archived=archive["cases"][ck]
        prior=next((z for z in old0001["cases"] if z["cap"]==cap),None) if mid==1 else None
        full_data={}; random_pool={}; replay={}
        for tr in TRACERS:
          k=A02.key(mid,cap,tr)
          for role in ROLES:
            src=paths[k+"/"+role]
            rows=A02.source_header_rows(src,prior,tr,role)
            fixed,full,info=E55.load_full_verified(
                src,mid=mid,cap=cap,tracer=tr,role=role,expected_rows=rows,
                archived_info=archived["input_sample_diagnostics"][tr+"_"+role],
                origproto=origproto)
            (full_data if role=="dat" else random_pool)[tr]=full
            replay[tr+"_"+role]={
              "eligible_rows":len(full[0]),
              "A02_selected_array_SHA256":info["selected_array_SHA256"]
            }

        probes,pmeta=E59.pair_probes(full_data["eBOSS_LRG"],full_data["eBOSS_ELG"],
                                     mid,cap,E59.NPROBE)
        parent=e59["cases"][ck]
        need(pmeta["probe_cartesian_SHA256"]==parent["probe_midpoints"]["probe_cartesian_SHA256"],
             "E59 probe midpoint replay changed "+ck)

        D,dwsum=full_data_kernel(full_data[SELECTED_TRACER],probes,args.chunk)
        rep_los={}
        for rep in REPS:
            if rep in case["replicas"] and "los" in case["replicas"][rep]:
                rep_los[rep]=case["replicas"][rep]["los"]
                print("E60_REUSE_REPLICA",ck,rep,flush=True)
                continue
            cat,pos=sample_random(random_pool[SELECTED_TRACER],rep_seed(mid,rep,cap))
            R,rwsum=random_kernel(cat,probes,args.chunk)
            z=reconstruction_los(D,R,probes)
            rec={
              "seed":rep_seed(mid,rep,cap),
              "random_rows":len(cat[0]),
              "selected_positions_SHA256":E59.arr_sha(pos),
              "selected_catalogue_SHA256":E55.catalogue_sha(cat),
              "random_weight_sum":rwsum,
              "los":z.tolist()
            }
            if rep=="R1":
                pR1=parent["fixed_random_mask"][SELECTED_TRACER]
                need(rec["seed"]==pR1["E58_R1_seed"]
                     and rec["selected_positions_SHA256"]==pR1["selected_positions_SHA256"]
                     and rec["selected_catalogue_SHA256"]==pR1["selected_catalogue_SHA256"],
                     "E59 R1 random subset replay changed "+ck)
                old=np.asarray(parent["tracers"][SELECTED_TRACER]["full_primary_los"],float)
                gap=float(np.max(np.abs(z-old)))
                need(gap<=5e-13,"E59 R1 reconstructed LOS replay drift "+ck)
                rec["E59_R1_max_abs_LOS_replay_gap"]=gap
            case["replicas"][rep]=rec
            rep_los[rep]=rec["los"]
            case["source_replay"]=replay
            case["probe_midpoints"]=pmeta
            case["data_weight_sum"]=dwsum
            atomic(OUT,state)
            print("E60_REPLICA_PASS",ck,rep,flush=True)

        case["metrics"]=case_metrics(rep_los)

        if mid==1:
            Rfull,rwsum=random_kernel(random_pool[SELECTED_TRACER],probes,args.chunk)
            zfull=reconstruction_los(D,Rfull,probes)
            zmean=np.asarray(case["metrics"]["sevenrep_mean_los"],float)
            diff=zmean-zfull
            denom=float(np.sqrt(np.mean(zfull*zfull)))
            case["full_pool_anchor"]={
              "full_eligible_random_rows":len(random_pool[SELECTED_TRACER][0]),
              "full_eligible_random_weight_sum":rwsum,
              "sevenmean_vs_fullpool_pearson":pearson(zmean,zfull),
              "sevenmean_vs_fullpool_sign_agreement":sign_agreement(zmean,zfull),
              "rms_difference_over_fullpool_field_rms":(
                  float(np.sqrt(np.mean(diff*diff))/denom) if denom>0 else math.inf
              )
            }
        case["status"]="complete"
        atomic(OUT,state)
        print("E60_CASE_PASS",ck,
              "PAIR_R_MED",case["metrics"]["pairwise_pearson"]["median"],
              "SIGN_MED",case["metrics"]["pairwise_sign_agreement"]["median"],
              "RMS_RATIO",case["metrics"]["random_rms_over_sevenrep_mean_rms"],flush=True)

    need(len(state["cases"])==18 and all(c["status"]=="complete" for c in state["cases"].values()),
         "Not all E60 cases complete")
    state["aggregate"]=aggregate(state)

    anchors={}
    for cap in CAPS:
        a=state["cases"][f"0001/{cap}"]["full_pool_anchor"]
        anchors[cap]={
          **a,
          "pass":bool(
            a["sevenmean_vs_fullpool_pearson"]>=FULLPOOL_PEARSON_MIN
            and a["sevenmean_vs_fullpool_sign_agreement"]>=FULLPOOL_SIGN_MIN
            and a["rms_difference_over_fullpool_field_rms"]<=FULLPOOL_RMS_DIFF_RATIO_MAX
          )
        }
    state["full_pool_anchor_summary"]=anchors

    capok=all(state["aggregate"][cap]["gate"]["pass"] for cap in CAPS)
    anchorok=all(anchors[cap]["pass"] for cap in CAPS)
    state["decision"]={
      "mask_robustness_pass":bool(capok and anchorok),
      "interpretation":(
        "ELG_RECONSTRUCTION_MASK_ROBUSTNESS_PASS; NEXT_TRUTH_VELOCITY_CALIBRATION"
        if capok and anchorok else
        "ELG_RECONSTRUCTION_MASK_ROBUSTNESS_FAIL; USE_FULL_POOL_OR_DENSER_DETERMINISTIC_SELECTION_BEFORE_TRUTH"
      ),
      "important_limit":(
        "E60 tests numerical/survey-selection random-mask robustness only. "
        "It does not measure correlation with true halo, matter, baryon, or neutrino-CDM relative velocity."
      )
    }
    state["status"]="PASS_ELG_MASK_ROBUSTNESS_QUANTIFIED"
    atomic(OUT,state)

    print("E60_WSL_ELG_MASK_ROBUSTNESS_PASS",flush=True)
    for cap in CAPS:
        a=state["aggregate"][cap]
        print("CAP",cap,
              "PAIR_R_MED",a["case_median_pairwise_pearson"]["median"],
              "PAIR_SIGN_MED",a["case_median_pairwise_sign_agreement"]["median"],
              "RANDOM_RMS_RATIO_MED",a["random_rms_over_mean_field_rms"]["median"],
              "RANDOM_RMS_RATIO_MAX",a["random_rms_over_mean_field_rms"]["max"],
              "CAP_GATE",a["gate"]["pass"],flush=True)
        b=anchors[cap]
        print("FULLPOOL",cap,
              "R",b["sevenmean_vs_fullpool_pearson"],
              "SIGN",b["sevenmean_vs_fullpool_sign_agreement"],
              "RMS_RATIO",b["rms_difference_over_fullpool_field_rms"],
              "PASS",b["pass"],flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
