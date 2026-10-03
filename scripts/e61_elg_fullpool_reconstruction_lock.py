#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,math,os,sys,tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import run_eboss_dr16_full_eligible_mock0001_sparse as V1
import e55_local_conditioned_full_eligible_ninemock_two48k as E55
import e59_mock_velocity_tag_repeatability as E59
import e60_elg_velocity_reconstruction_mask_robustness as E60

E59_RESULT=ROOT/"source_data/e59_mock_velocity_tag_repeatability.json"
E60_RESULT=ROOT/"source_data/e60_elg_velocity_reconstruction_mask_robustness.json"
OUT=ROOT/"source_data/e61_elg_fullpool_reconstruction_lock.json"
IDS=A02.IDS; CAPS=A02.CAPS; TRACERS=A02.TRACERS; ROLES=A02.ROLES
TRACER="eBOSS_ELG"

CAP_R=0.98; CAP_S=0.95; CAP_D=0.25
CASE_R=0.95; CASE_S=0.90; CASE_D=0.50

def need(c,m):
    if not c: raise RuntimeError(m)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e61_",delete=False) as f:
        q=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(q,path)
    finally: q.unlink(missing_ok=True)

def pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a-a.mean(); b=b-b.mean()
    d=float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(a@b/d) if d>0 else 0.0

def sign_agree(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    m=(a!=0)&(b!=0)
    need(m.sum()>=.95*len(a),"too many zero LOS values")
    return float(np.mean(np.sign(a[m])==np.sign(b[m])))

def compare(s,f):
    s=np.asarray(s,float); f=np.asarray(f,float)
    rf=float(np.sqrt(np.mean(f*f))); rd=float(np.sqrt(np.mean((s-f)**2)))
    return {"pearson":pearson(s,f),"sign_agreement":sign_agree(s,f),
            "rms_difference_ratio":rd/rf if rf>0 else math.inf,
            "fullpool_field_rms":rf,
            "sevenrep_mean_field_rms":float(np.sqrt(np.mean(s*s)))}

def desc(v):
    a=np.asarray(v,float)
    return {"n":len(a),"mean":float(a.mean()),"median":float(np.median(a)),
            "min":float(a.min()),"max":float(a.max()),
            "sample_sd":float(a.std(ddof=1))}

def aggregate(state):
    out={}
    for cap in CAPS:
        cs=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        r=desc([c["comparison"]["pearson"] for c in cs])
        s=desc([c["comparison"]["sign_agreement"] for c in cs])
        d=desc([c["comparison"]["rms_difference_ratio"] for c in cs])
        ok=(r["median"]>=CAP_R and s["median"]>=CAP_S and d["median"]<=CAP_D
            and r["min"]>=CASE_R and s["min"]>=CASE_S and d["max"]<=CASE_D)
        out[cap]={"pearson":r,"sign_agreement":s,"rms_difference_ratio":d,
                  "sevenrep_surrogate_gate":{
                    "cap_median_pearson_min":CAP_R,"cap_median_sign_min":CAP_S,
                    "cap_median_rms_ratio_max":CAP_D,"case_pearson_min":CASE_R,
                    "case_sign_min":CASE_S,"case_rms_ratio_max":CASE_D,
                    "pass":bool(ok)}}
    return out

def self_test():
    x=np.array([1.,-2.,3.,4.]); y=np.array([1.1,-1.9,2.9,4.1])
    z=compare(x,y)
    need(z["pearson"]>.99 and z["sign_agreement"]==1.0,"comparison self-test")
    print("E61_SYNTHETIC_FULLPOOL_LOCK_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    a=ap.parse_args()
    self_test()
    if a.self_test:return
    need(a.run,"use --run")
    need(16<=a.chunk<=256,"unsafe chunk")
    e59=json.loads(E59_RESULT.read_text()); e60=json.loads(E60_RESULT.read_text())
    need(e59["decision"]["selected_tracer_reconstruction_for_truth_calibration"]==TRACER,
         "E59 selected tracer changed")
    need(e60["status"]=="PASS_ELG_MASK_ROBUSTNESS_QUANTIFIED"
         and e60["decision"]["mask_robustness_pass"] is False,"E60 parent changed")
    need(e60["observed_odd_used"] is False and e60["observed_galaxy_rows_used"] is False,
         "observation guardrail changed")

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
        need(state["stage"]=="E61_ELG_FULLPOOL_RECONSTRUCTION_LOCK","checkpoint mismatch")
    else:
        state={"stage":"E61_ELG_FULLPOOL_RECONSTRUCTION_LOCK","date":"2026-09-30",
               "status":"INCOMPLETE","selected_tracer":TRACER,
               "full_eligible_random_pool_is_final_selection_operator":True,
               "sevenrep_mean_is_surrogate_diagnostic_only":True,
               "same_E59_kernel":True,"same_E59_probe_policy":True,
               "cases":{},"observed_galaxy_rows_used":False,"observed_odd_used":False,
               "true_velocity_used":False,"absolute_EV_amplitude_calibrated":False,
               "covariance_inverse_used":False,"pvalue_or_detection_sigma":False}
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        c=state["cases"].setdefault(ck,{"id":mid,"cap":cap,"status":"INCOMPLETE"})
        if c.get("status")=="complete":
            print("E61_REUSE_CASE",ck,flush=True); continue
        archived=archive["cases"][ck]
        prior=next((z for z in old0001["cases"] if z["cap"]==cap),None) if mid==1 else None
        data={}; ran={}
        for tr in TRACERS:
          k=A02.key(mid,cap,tr)
          for role in ROLES:
            src=paths[k+"/"+role]
            rows=A02.source_header_rows(src,prior,tr,role)
            fixed,full,info=E55.load_full_verified(
                src,mid=mid,cap=cap,tracer=tr,role=role,expected_rows=rows,
                archived_info=archived["input_sample_diagnostics"][tr+"_"+role],
                origproto=origproto)
            (data if role=="dat" else ran)[tr]=full

        probes,pmeta=E59.pair_probes(data["eBOSS_LRG"],data["eBOSS_ELG"],mid,cap,E59.NPROBE)
        need(pmeta["probe_cartesian_SHA256"]==
             e59["cases"][ck]["probe_midpoints"]["probe_cartesian_SHA256"],
             "E59 probe replay changed "+ck)
        D,dws=E60.full_data_kernel(data[TRACER],probes,a.chunk)
        R,rws=E60.random_kernel(ran[TRACER],probes,a.chunk)
        full=E60.reconstruction_los(D,R,probes)
        seven=np.asarray(e60["cases"][ck]["metrics"]["sevenrep_mean_los"],float)
        cmp=compare(seven,full)

        if mid==1:
            old=e60["full_pool_anchor_summary"][cap]
            need(abs(cmp["pearson"]-old["sevenmean_vs_fullpool_pearson"])<5e-13,
                 "E60 Pearson replay drift")
            need(abs(cmp["sign_agreement"]-old["sevenmean_vs_fullpool_sign_agreement"])<5e-13,
                 "E60 sign replay drift")
            need(abs(cmp["rms_difference_ratio"]-old["rms_difference_over_fullpool_field_rms"])<5e-13,
                 "E60 RMS replay drift")

        c.update({"probe_midpoints":pmeta,
                  "full_eligible_ELG_random_rows":len(ran[TRACER][0]),
                  "full_eligible_ELG_random_weight_sum":rws,"ELG_data_weight_sum":dws,
                  "comparison":cmp,"fullpool_reconstruction_los":full.tolist(),
                  "status":"complete"})
        atomic(OUT,state)
        print("E61_CASE_PASS",ck,"R",cmp["pearson"],"SIGN",cmp["sign_agreement"],
              "RMS",cmp["rms_difference_ratio"],flush=True)

    state["aggregate"]=aggregate(state)
    surrogate=all(state["aggregate"][c]["sevenrep_surrogate_gate"]["pass"] for c in CAPS)
    state["decision"]={
      "full_pool_reconstruction_operator_locked":True,
      "sevenrep_surrogate_valid_all_caps":bool(surrogate),
      "interpretation":("FULLPOOL_ELG_RECONSTRUCTION_LOCKED; SEVENREP_SURROGATE_VALID; NEXT_TRUTH_VELOCITY_CALIBRATION"
        if surrogate else
        "FULLPOOL_ELG_RECONSTRUCTION_LOCKED; SEVENREP_SURROGATE_REJECTED; NEXT_TRUTH_VELOCITY_CALIBRATION_MUST_USE_FULLPOOL_OPERATOR"),
      "important_limit":"No truth-velocity calibration yet."}
    state["status"]="PASS_ELG_FULLPOOL_RECONSTRUCTION_LOCK_QUANTIFIED"
    atomic(OUT,state)

    print("E61_WSL_FULLPOOL_LOCK_PASS",flush=True)
    for cap in CAPS:
        z=state["aggregate"][cap]
        print("CAP",cap,"R_MED",z["pearson"]["median"],"R_MIN",z["pearson"]["min"],
              "SIGN_MED",z["sign_agreement"]["median"],"SIGN_MIN",z["sign_agreement"]["min"],
              "RMS_MED",z["rms_difference_ratio"]["median"],"RMS_MAX",z["rms_difference_ratio"]["max"],
              "SURROGATE_PASS",z["sevenrep_surrogate_gate"]["pass"],flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
