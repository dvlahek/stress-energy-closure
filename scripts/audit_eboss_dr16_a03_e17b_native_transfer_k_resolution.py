#!/usr/bin/env python3
"""E17B preregistered CLASS native k-resolution stress of frozen E17A both-leg linear transfers.

Fresh CLASS at two numerical resolutions only. No new distribution, triangle, halo kernel,
galaxy estimator, observed odd, catalogue, mock, science cut, seed or physical B.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"code"), str(ROOT/"scripts")]
import audit_eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk as e17
import audit_eboss_dr16_a03_e14_multik_direct_vTk_source_shape as e14

PRE = ROOT/"source_data/eboss_dr16_a03_e17b_native_transfer_k_resolution_prereg_2026-09-27.json"
PRE_BLOB = "d0389e336de9e9964cade589fc17594a595f2e17"
ARCH = ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27"
CSV = ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
OUT = ROOT/"eboss_workspace/a03_physics_source/e17b_native_transfer"
STATES=("FD","plus","minus")
ZS=(.945,.95,.955)
FIELDS=("Pcb","theta","dcb")
OMEGA_B=.02237
OMEGA_CDM=.1200

def require(b,s):
    if not b: raise RuntimeError("E17B_FAIL_CLOSED: "+s)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def packed(obj):return (json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
def save_once(path,obj):
    raw=packed(obj);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"immutable output collision at "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17b_",delete=False) as f:
            tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17B_SAVED",str(path),"SHA256",sha(raw),flush=True)
def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"preregistration blob changed")
    p=json.loads(PRE.read_bytes())
    require(p["status"]=="E17B_PROSPECTIVELY_REGISTERED_BEFORE_ANY_REFINED_CLASS_EXECUTION"
        and p["absolute_STOP"]["no_observed_odd_read"] and p["absolute_STOP"]["no_main_mutation"],
        "scope or observed STOP changed")
    a=p["parent"]
    require(blob((ARCH/"archive_manifest.json").read_bytes())==a["E17A_archive_manifest_git_blob"],
        "E17A archive manifest git blob changed")
    for fn,expected in [
        ("e17_original_joint_three_state_CLASS_both_short_legs_source_only.json",a["E17A_original_joint_sha256"]),
        ("e17_independent_fresh_CLASS_scalar_replay.json",a["E17A_independent_sha256"])]:
        require(sha((ARCH/fn).read_bytes())==expected,"E17A joint/independent missing: "+fn)
    for st in STATES:
        require(sha((ARCH/("e17_original_both_legs_"+st+".json")).read_bytes())==
                a["E17A_original_per_state_sha256"][st],"E17A original state changed: "+st)
    require(sha(CSV.read_bytes())==a["E8_original_csv_sha256"]
            and sha(e17.E16RAW.read_bytes())==a["E16_original_sha256"]
            and blob(e17.PRE.read_bytes())==a["E17A_protocol_git_blob"],
            "original E8/E16/E17A parent drift")
    require(p["frozen"]["states"]==list(STATES) and p["frozen"]["tracer_order"]==["LRG","ELG"]
            and p["frozen"]["geometry_count"]==576 and p["frozen"]["z"]==list(ZS)
            and p["precision"]["tier_high"]=={"k_per_decade_for_pk":40,"k_per_decade_for_bao":140}
            and p["precision"]["tier_ultra"]=={"k_per_decade_for_pk":80,"k_per_decade_for_bao":280},
            "predeclared precision/geometry altered")
    e17.gate()
    return p
def flatten_original(state):
    x=json.loads((ARCH/("e17_original_both_legs_"+state+".json")).read_bytes())
    require(x["state"]==state and x["CLASS_transfer_grid"]["N"]==81
            and len(x["rows_original_576_geometry_three_z_two_legs"])==576
            and x["observed_odd_read"] is False,
            "baseline 81-node state/observed provenance drift")
    rows=x["rows_original_576_geometry_three_z_two_legs"]
    fresh=e17.geometry()
    require(len(fresh)==len(rows),"original E16 geometry count changed")
    legs=[]
    for n,(old,new) in enumerate(zip(rows,fresh)):
        require((old["k"],old["K"],old["mu_s"],old["mu_L"],old["phi"])==
                (new["k"],new["K"],new["mu_s"],new["mu_L"],new["phi"]),
                "original E17A 48 orientation order changed")
        for lab in ("k1_h_Mpc","k2_h_Mpc"):
            require(abs(old[lab]-new[lab])<1e-14,"original E17A geometry leg changed")
            legs.append(old[lab])
    values={}
    for z in ZS:
        d={"Pcb":[],"theta":[],"dcb":[]}
        for row in rows:
            for l in ("leg1","leg2"):
                o=row["z_nodes"][str(z)][l]
                d["Pcb"].append(o["P_cb_CLASS_Mpc3"])
                d["theta"].append(o["theta_rel_DIRECT_vTk_per_R_Mpc_inv"])
                d["dcb"].append(o["delta_cb_Newtonian_per_R"])
        values[str(z)]=d
    return np.asarray(legs,float),values,x

def evaluated_tier(cosmo,h,legs,setting):
    actual={}
    datasets={}
    for z in ZS:
        t=cosmo.get_transfer(z=z,output_format="class")
        require({"k (h/Mpc)","t_ncdm[0]","t_cdm","d_b","d_cdm"}.issubset(t),
                "refined CLASS direct vTk/CB transfer columns missing")
        kin=np.asarray(t["k (h/Mpc)"],float)
        require(np.isfinite(kin).all() and np.all(np.diff(kin)>0)
                and kin[0]<min(legs) and kin[-1]>max(legs),
                "refined CLASS does not bracket every actual leg")
        if actual:require(np.array_equal(kin,actual["k_h_Mpc"]),
                          "refined native transfer grid changed across z")
        else:actual={"k_h_Mpc":kin.tolist(),"N":len(kin)}
        i=np.searchsorted(kin,legs,side="right")
        require((i>0).all() and (i<len(kin)).all(),"refined native transfer extrapolation")
        tnu=np.asarray(t["t_ncdm[0]"],float)
        tcdm=np.asarray(t["t_cdm"],float)
        db=np.asarray(t["d_b"],float)
        dc=np.asarray(t["d_cdm"],float)
        delta=(OMEGA_B*db+OMEGA_CDM*dc)/(OMEGA_B+OMEGA_CDM)
        values={
            "Pcb":[float(cosmo.pk_cb_lin(float(k)*h,z)) for k in legs],
            "theta":np.interp(legs,kin,tnu-tcdm).tolist(),
            "dcb":np.interp(legs,kin,delta).tolist()}
        require(all(np.isfinite(values[f]).all() for f in FIELDS)
                and min(values["Pcb"])>0,"refined CLASS field nonfinite/nonpositive")
        datasets[str(z)]={"native":{
                "t_ncdm[0]":tnu.tolist(),"t_cdm":tcdm.tolist(),
                "d_b":db.tolist(),"d_cdm":dc.tolist()},
                "actual_leg_samples":values}
    return {"precision_parameters":setting,"native_grid":actual,"z_nodes":datasets}

def metrics(tier_a,tier_b,field,scale_floor):
    diffs=[];worst=None
    for zi,z in enumerate(ZS):
        va=tier_a[str(z)][field];vb=tier_b[str(z)][field]
        require(len(va)==1152 and len(vb)==1152,"one original E16 two-leg sample missing")
        for j,(a,b) in enumerate(zip(va,vb)):
            v=abs(a-b)/max(abs(a),abs(b),scale_floor)
            diffs.append(v)
            if worst is None or v>worst["scaled_gap"]:
                worst={"scaled_gap":float(v),"absolute_gap":float(abs(a-b)),
                        "z":z,"geometry_index":j//2,"leg_index":1+j%2,
                        "original_pair_index":(j//2)//48,
                        "old_value":float(a),"new_value":float(b)}
    return {"N":len(diffs),"max_scaled_gap":max(diffs),"p95_scaled_gap":
            float(np.quantile(diffs,.95)),"rms_scaled_gap":
            float(math.sqrt(np.mean(np.square(diffs)))),"worst":worst}

def run_state(st,out):
    p=gate();require(st in STATES,"unknown original F state")
    import class_response_optimize as cro
    import wake_two_tracer_fisher as base
    from classy import Class
    require(cro.CLASS_COMMIT==p["parent"]["original_CLASS_commit"]
        and cro.OMEGA_B==OMEGA_B and cro.OMEGA_CDM==OMEGA_CDM
        and base.NS==.9649,"original CLASS cosmology changed")
    legs,baseline,original=flatten_original(st)
    arr=np.genfromtxt(CSV,delimiter=",",names=True)
    q=np.asarray(arr["q_dimensionless"],float)
    f=np.asarray(arr[e14.COLS[st]],float)
    require(len(q)==4000 and np.min(f)>0,"original PSD invalid")
    out.mkdir(parents=True,exist_ok=True)
    psd=out/("e17b_original_E13_"+st+".dat")
    source="\n".join("%.14e %.14e"%(x,y) for x,y in zip(q,f))+"\n"
    if psd.exists():require(psd.read_text()==source,"original PSD bytes changed")
    else:psd.write_text(source)
    require(sha(psd.read_bytes())==original["original_E13_PSD_SHA256"],
            "original E17A PSD is not byte-exact")
    base.ZBINS=np.asarray([.9,1.],float)
    params=base.class_params(psd,.06)
    params["output"]="mPk,dTk,vTk"
    require(params["gauge"]=="newtonian" and params["P_k_max_h/Mpc"]>=.105
            and params["z_max_pk"]>=.955,"baseline CLASS range or gauge changed")
    computed={}
    for tier in ("tier_high","tier_ultra"):
        precise=dict(params);precise.update(p["precision"][tier])
        cosm=Class();cosm.set(precise);cosm.compute()
        try:
            h=float(cosm.h())
            require(abs(h-original["h"])<1e-12,"refined CLASS h changed")
            print("E17B_FRESH_CLASS_TIER",st,tier,flush=True)
            computed[tier]=evaluated_tier(cosm,h,legs,p["precision"][tier])
            computed[tier]["E14_three_short_k_Pcb_Mpc3"]={str(k):float(cosm.pk_cb_lin(k*h,.95)) for k in e17.KS}
        finally:cosm.struct_cleanup();cosm.empty()
    nh=computed["tier_high"]["native_grid"]["N"]
    nu=computed["tier_ultra"]["native_grid"]["N"]
    sizes={"baseline":original["CLASS_transfer_grid"]["N"],
           "tier_high":nh,"tier_ultra":nu}
    require(nh>sizes["baseline"] and nu>nh,
            "refined CLASS settings failed to increase native transfer-k node count")
    summary={}
    for field in FIELDS:
        high={z:computed["tier_high"]["z_nodes"][z]["actual_leg_samples"] for z in map(str,ZS)}
        ultra={z:computed["tier_ultra"]["z_nodes"][z]["actual_leg_samples"] for z in map(str,ZS)}
        scale_floor=(1e-6*max(abs(v) for z in ultra.values() for v in z[field])
                     if field!="Pcb" else 1e-30)
        summary[field]={"baseline_vs_ultra":metrics(baseline,ultra,field,scale_floor),
                        "high_vs_ultra":metrics(high,ultra,field,scale_floor),
                        "scale_floor_for_near_zero":scale_floor}
    candidate=p["prereg_engineering_diagnostics_NOT_physical_error_budget"]
    archived_anchors=json.loads((e17.E14DIR/("e14_short_Pcb_"+st+".json")).read_bytes())["CLASS_Pcb_short_per_state_Mpc3"]
    qa={"grid_refinement":sizes,"max_original_anchor_Pcb_relative_gap":
        max(abs(computed[tier]["E14_three_short_k_Pcb_Mpc3"][str(k)]-
                archived_anchors[str(k)])/archived_anchors[str(k)]
            for tier in ("tier_high","tier_ultra") for k in e17.KS),
        "baseline_ultra_transfer_candidate":
        max(summary[x]["baseline_vs_ultra"]["max_scaled_gap"] for x in ("theta","dcb"))<=
            candidate["max_baseline_ultra_scaled_transfer_gap_CANDIDATE"],
        "high_ultra_transfer_candidate":
        max(summary[x]["high_vs_ultra"]["max_scaled_gap"] for x in ("theta","dcb"))<=
            candidate["max_high_ultra_scaled_transfer_gap_CANDIDATE"]}
    qa["CLASS_Pcb_anchor_candidate"]=qa["max_original_anchor_Pcb_relative_gap"]<=candidate[
        "max_original_anchor_Pcb_relative_gap"]
    qa["engineering_candidate_not_physical_certification"]=all(
        qa[x] for x in ("baseline_ultra_transfer_candidate",
                      "high_ultra_transfer_candidate","CLASS_Pcb_anchor_candidate"))
    result={"date":"2026-09-27","status":
        "E17B_REFINED_NATIVE_CLASS_LINEAR_TRANSFER_DIAGNOSTIC_COMPLETE_FULL_PHYSICAL_B_BLOCKED",
        "state":st,"prospective_protocol_git_blob":PRE_BLOB,
        "original_E17A_report_sha256":p["parent"]["E17A_original_per_state_sha256"][st],
        "original_E13_PSD_sha256":sha(psd.read_bytes()),
        "CLASS_commit":cro.CLASS_COMMIT,
        "geometry_count":576,"z_nodes":list(ZS),"both_legs_per_z":1152,
        "omega_b":OMEGA_B,"omega_cdm":OMEGA_CDM,
        "native_CLASS_evaluations_by_fixed_tier":computed,"convergence_diagnostics":summary,
        "QA":qa,"error_scale_definition":p["precision"]["scaled_error"],
        "provenance":"CLASS direct Pcb each actual leg; CLASS direct vTk and d_cb native transfer interpolation within verified refined grids. Neither the field nor the tested resolution proves a retarded halo/tracer response.",
        "original_R16_rank_not_replaced_by_short_mode":True,
        "observed_odd_read":False,"new_mock_catalogue_download_science_cut_seed":False,
        "physical_finite_K_B_xi_SNR_significance":None,"PR_draft_main_untouched":True}
    save_once(out/("e17b_refined_native_CLASS_"+st+".json"),result)
    print("E17B_STATE",st,"GRID",sizes,"THETA_BASELINE_ULTRA",
        summary["theta"]["baseline_vs_ultra"]["max_scaled_gap"],
        "THETA_HIGH_ULTRA",summary["theta"]["high_vs_ultra"]["max_scaled_gap"],
        "CANDIDATE",qa["engineering_candidate_not_physical_certification"],flush=True)

def aggregate(out):
    p=gate()
    reports={}
    for st in STATES:
        f=out/("e17b_refined_native_CLASS_"+st+".json")
        raw=f.read_bytes();j=json.loads(raw)
        require(j["state"]==st and j["prospective_protocol_git_blob"]==PRE_BLOB
                and j["geometry_count"]==576 and not j["observed_odd_read"]
                and j["physical_finite_K_B_xi_SNR_significance"] is None,
                "original state E17B report invalid")
        reports[st]={"sha256":sha(raw),"size_bytes":len(raw),
            "QA":j["QA"],"diagnostics":{fld:{
                k:j["convergence_diagnostics"][fld][k]["max_scaled_gap"]
                for k in ("baseline_vs_ultra","high_vs_ultra")} for fld in FIELDS}}
    result={"date":"2026-09-27","status":
        "E17B_THREE_STATE_NATIVE_TRANSFER_K_REFINEMENT_DIAGNOSTIC_COMPLETE_PHYSICAL_B_BLOCKED",
        "prospective_protocol_git_blob":PRE_BLOB,"original_E17A_joint_sha256":
        p["parent"]["E17A_original_joint_sha256"],"state_reports":reports,
        "all_three_state_qa_candidate":all(reports[s]["QA"][
            "engineering_candidate_not_physical_certification"] for s in STATES),
        "interpretation":"Resolution stability is a conditional engineering check, not certified transfer interpolation truncation bound or retarded physical three-point halo coupling.",
        "E14_E15_original_72_24_unmodified":True,"observed_odd_read":False,
        "physical_finite_K_bispectrum":"BLOCKED","main_untouched_PR_draft":True}
    save_once(out/"e17b_original_joint_native_CLASS_transfer_resolution.json",result)
    print("E17B_JOINT",result["all_three_state_qa_candidate"],flush=True)

def check_result(out):
    a=json.loads((out/"e17b_original_joint_native_CLASS_transfer_resolution.json").read_bytes())
    require(a["all_three_state_qa_candidate"],
            "predeclared engineering resolution candidate not met; archived reports are authoritative")
    print("E17B_PREREG_ENGINEERING_RESOLUTION_CANDIDATE_PASS_NOT_PHYSICAL_B",flush=True)

if __name__=="__main__":
    ar=argparse.ArgumentParser()
    g=ar.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--state",choices=STATES)
    g.add_argument("--aggregate",action="store_true")
    g.add_argument("--check",action="store_true")
    ar.add_argument("--output-dir",type=Path,default=OUT)
    a=ar.parse_args()
    if a.self_test:
        gate();print("E17B_PREREG_PARENT_SHA_GEOMETRY_OBSERVED_STOP_PASS",flush=True)
    elif a.state:run_state(a.state,a.output_dir)
    elif a.aggregate:aggregate(a.output_dir)
    else:check_result(a.output_dir)
