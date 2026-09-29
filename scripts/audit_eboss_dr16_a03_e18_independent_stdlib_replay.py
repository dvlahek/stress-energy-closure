#!/usr/bin/env python3
"""Independent stdlib audit of RETROSPECTIVE existing E4 9-mock 24D summary.

Reads four small tracked Git JSONs only. No FITS, observed data, WSL,
new randoms, covariance inversion, p values or physical inference.
"""
import copy
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/"source_data"
INPUTS={
"protocol":("eboss_dr16_a03_e18_posthoc_nine_mock_paired_random_24d_stability_protocol_2026-09-29.json","9b13de770d4c9b051f8b6af7c5a6e2a342df77e3"),
"original":("eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json","32246a495506deac5b9dfe71bea51c3db33202b1"),
"erratum":("eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json","0baf9c5a009bf8f5e78f5725e44395bad6f91715"),
"result":("eboss_dr16_a03_e18_posthoc_nine_mock_24d_paired_random_vs_mock_scatter_descriptive_result_2026-09-29.json","86ea80bca942143d20e93123aa7d8d0892431643"),
}
IDS=[1,125,250,375,500,625,750,875,1000]
CAPS=["NGC","SGC"]
LEVELS=["full_1200_original","nested_2400","nested_4800"]
SIZE={"full_1200_original":1200,"nested_2400":2400,"nested_4800":4800}
CMP=[
("full_1200_original","nested_2400","nested_2400_minus_full_1200_original"),
("nested_2400","nested_4800","nested_4800_minus_nested_2400"),
("full_1200_original","nested_4800","nested_4800_minus_full_1200_original"),
]
GROUPS={
"NGC_ell1":range(0,6),"NGC_ell3":range(6,12),
"SGC_ell1":range(12,18),"SGC_ell3":range(18,24),
"NGC_12D":range(12),"SGC_12D":range(12,24),
"joint_24D":range(24)
}

def need(condition,message):
    if not condition:
        raise ValueError("E18_INDEPENDENT_SOURCE_ONLY_FAIL_CLOSED: "+message)

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def inputs():
    obj={}
    for label,(name,sha) in INPUTS.items():
        p=ROOT/name
        need(p.is_file() and not p.is_symlink(),"source file missing")
        raw=p.read_bytes()
        need(blob(raw)==sha,"input Git blob drift: "+label)
        obj[label]=json.loads(raw)
    p,e,err,res=(obj[x] for x in ("protocol","original","erratum","result"))
    need(p["classification"].startswith("EXPLICITLY_POSTHOC") and
         p["descriptive_outputs_before_computation"]["no_posthoc_acceptance_threshold"] and
         all(p["absolute_STOP"].values()),"not original retrospective protocol")
    need(e["completed_cases"]==18 and e["failed_cases"]==0 and
         not e["observed_galaxy_rows_read"] and not e["observed_random_rows_read"] and
         not e["observed_odd_data_vector_read"] and not e["new_science_selection_applied"],
         "original E4 observed/status contract")
    need(err["original_report_git_blob_sha1"]==INPUTS["original"][1] and
         err["original_report_sha256"]==p["parent_original_E4"]["original_user_JSON_sha256"] and
         err["no_observed_rows_or_odd_read"] and err["pilot_joint_dim"]==24,
         "original independent erratum provenance")
    need(e["fixed_ids"]==IDS and err["fixed_original_IDs"]==IDS and
         e["caps"]==CAPS and err["caps"]==CAPS and
         e["random_levels"]==LEVELS and len(e["cases"])==18,
         "original matched fixed cohort mismatch")
    need(res["observed_odd_SEALED"] and res["nine_ID_centered_covariance_rank_upper_bound"]==8 and
         res["E7_48k_fullgalaxy_different_DD_excluded"] and
         res["not_inferential_covariance_or_wake_physics"] and
         res["original_E4_galaxy_sample_identical_across_nested_random_levels"],
         "reported result overclaim")
    return e,err,res

def original_vectors(e):
    vs={}
    for i in IDS:
        k=f"{i:04d}"
        vs[k]={}
        for level in LEVELS:
            vector=[]
            for cap in CAPS:
                case=e["cases"][k+"/"+cap]
                c=case["levels"][level]
                b=c["baseline_original"]
                need(case["id"]==i and case["cap"]==cap and case["status"]=="complete",
                     "mock ID/cap has been substituted")
                need(c["LRG_R_rows"]==SIZE[level] and c["ELG_R_rows"]==SIZE[level] and
                     b["status"]=="complete_full_RR_support" and
                     b["RR_supported_cells"]==144 and
                     b["RR_missing_indices"]==[] and
                     b["no_unsupported_RR_zero_fill_or_projection"],
                     "RR support or nested random count changed")
                fixed_DD=case["levels"]["full_1200_original"]["baseline_original"]["forward_pair_terms"]["D1D2"]["weighted_histogram_SHA256"]
                need(b["forward_pair_terms"]["D1D2"]["weighted_histogram_SHA256"]==fixed_DD,
                     "not same DD mock galaxy sample")
                for ell in ("1","3"):
                    v=b["original_physical_domain_multipoles"][ell]["values_by_fixed_s_bin"]
                    need(len(v)==6 and all(math.isfinite(x) for x in v),"invalid six-bin odd")
                    vector.extend(v)
            need(len(vector)==24,"frozen 24D order wrong")
            vs[k][level]=vector
    return vs

def close(a,b):
    need(a is not None and b is not None and math.isfinite(a) and math.isfinite(b) and
         math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),"descriptive arithmetic mismatch")

def verify(e,err,res):
    v=original_vectors(e)
    for low,high,key in CMP:
        recorded=res["comparisons"][key]
        need(recorded["low_random_level"]==low and recorded["high_random_level"]==high,
             "wrong paired comparison")
        high_rows=[v[f"{ident:04d}"][high] for ident in IDS]
        diff_rows=[]
        for ident in IDS:
            k=f"{ident:04d}"
            delta=[z-y for z,y in zip(v[k][high],v[k][low])]
            archived=err["joint_24D_drift_by_fixed_mock_ID"][k][key]["signed_24D_pilot_shift"]
            need(len(archived)==24 and all(math.isclose(z,y,rel_tol=1e-12,abs_tol=1e-12) for z,y in zip(delta,archived)),
                 "original independent signed E4 24D erratum mismatch")
            diff_rows.append(delta)
            r=recorded["all_9_fixed_IDs"][k]
            close(math.sqrt(sum(x*x for x in delta)),r["paired_shift_joint_24D_L2"])
            close(max(abs(x) for x in delta),r["paired_shift_max_abs_component"])
        need(len(recorded["all_9_fixed_IDs"])==9 and len(recorded["all_24_components"])==24,
             "fixed mock or full component loss")
        for group,ii in GROUPS.items():
            count=len(ii)
            delta_sq=sum(diff_rows[i][j]**2 for i in range(9) for j in ii)
            means={j:sum(high_rows[i][j] for i in range(9))/9 for j in ii}
            centered_sq=sum((high_rows[i][j]-means[j])**2 for i in range(9) for j in ii)
            paired=math.sqrt(delta_sq/(9*count))
            between=math.sqrt(centered_sq/(8*count))
            r=recorded["by_group"][group]
            close(paired,r["paired_RMS"])
            close(between,r["between_ID_high_level_centered_RMS"])
            close(paired/between if between else None,r["paired_over_between_descriptive_ratio"])
        for j,r in enumerate(recorded["all_24_components"]):
            delta=[diff_rows[i][j] for i in range(9)]
            h=[high_rows[i][j] for i in range(9)]
            hmean=sum(h)/9
            paired=math.sqrt(sum(x*x for x in delta)/9)
            between=math.sqrt(sum((x-hmean)**2 for x in h)/8)
            need(r["index"]==j and r["cap"]==("NGC" if j<12 else "SGC") and
                 r["ell"]==(1 if j%12<6 else 3) and
                 r["s_Mpc_per_h"]==(30,50,70,90,110,130)[j%6],"component indexing drift")
            close(sum(delta)/9,r["signed_paired_mean"])
            close(paired,r["paired_RMS"])
            close(between,r["between_ID_high_level_sd"])
            close(paired/between if between else None,r["paired_over_between_descriptive_ratio"])
    return True

def reject(label,func):
    try:func()
    except (ValueError,KeyError,TypeError,IndexError):
        print("E18_INDEPENDENT_NEGATIVE_REJECT",label,flush=True)
    else:raise AssertionError("accepted synthetic negative: "+label)

def main():
    e,err,res=inputs()
    verify(e,err,res)
    changed=copy.deepcopy(e)
    changed["cases"]["0001/NGC"]["levels"]["nested_4800"]["baseline_original"]["RR_supported_cells"]=143
    reject("SILENT_RR_HOLE",lambda:verify(changed,err,res))
    changed=copy.deepcopy(e)
    changed["cases"]["1000/SGC"]["levels"]["nested_2400"]["baseline_original"]["forward_pair_terms"]["D1D2"]["weighted_histogram_SHA256"]="changedDD"
    reject("GALAXY_DD_CHANGE",lambda:verify(changed,err,res))
    changed=copy.deepcopy(err)
    changed["joint_24D_drift_by_fixed_mock_ID"]["0001"]["nested_4800_minus_nested_2400"]["signed_24D_pilot_shift"][0]=100
    reject("ERRATUM_SIGNED_DELTA_TAMPER",lambda:verify(e,changed,res))
    changed=copy.deepcopy(res)
    changed["comparisons"]["nested_4800_minus_nested_2400"]["by_group"]["joint_24D"]["paired_RMS"]=0
    reject("REPORTED_NUMERIC_DRIFT",lambda:verify(e,err,changed))
    changed=copy.deepcopy(e)
    del changed["cases"]["0125/SGC"]
    reject("MISSING_MATCHED_MOCK_CAP",lambda:verify(changed,err,res))
    print("E18_FIVE_INDEPENDENT_NEGATIVE_CONTROLS_PASS",flush=True)
    for name in ("nested_2400_minus_full_1200_original",
                 "nested_4800_minus_nested_2400",
                 "nested_4800_minus_full_1200_original"):
        g=res["comparisons"][name]["by_group"]["joint_24D"]
        print("E18_INDEPENDENT_REPLAY",name,
              "paired_RMS",format(g["paired_RMS"],".12g"),
              "between_ID_RMS",format(g["between_ID_high_level_centered_RMS"],".12g"),
              "DESCRIPTIVE_RATIO_NOT_SIGMA",format(g["paired_over_between_descriptive_ratio"],".12g"),
              flush=True)
    print("E18_ALL_9_IDS_3_LEVELS_24D_ORIGINAL_ERRATUM_INDEPENDENT_STD_LIB_PASS",flush=True)
    print("E18_NO_OBSERVED_ODD_NO_INFERENTIAL_COVARIANCE_NO_LOCAL_WSL",flush=True)

if __name__=="__main__":
    main()
