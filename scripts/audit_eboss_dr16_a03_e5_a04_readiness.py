#!/usr/bin/env python3
"""A-03E5/A-04 pure source-only physics-first readiness audit.

PASS means correct and enforceable STOP, NOT eBOSS physical-window/covariance
certification. Read only immutable source_data JSON/log/protocols, never FITS,
observed galaxies/random rows/odd vector, and never fetch new sources.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
import statistics
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"source_data"
PROTO=SRC/"eboss_dr16_a03_e5_a04_physics_first_readiness_protocol_2026-09-27.json"
REPORT=SRC/"eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json"
LOG=SRC/"eboss_dr16_a03_e4_local_stdout_2026-09-27.log"
MANIFEST=SRC/"eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json"
E4PROTO=SRC/"eboss_dr16_a03_e4_nested_1200_2400_4800_mock_random_density_protocol_2026-09-27.json"
A02=SRC/"eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
E3=SRC/"eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json"
USER=SRC/"eboss_dr16_a03_empirical_only_user_decision_2026-09-26.json"
E01=SRC/"eboss_dr16_a03_empirical_rr_window_injection_report_2026-09-26.json"
DESICOV=SRC/"lrg_elg_covariance_r4_120_checkpoint_2026-09-23.json"
DESIRESULT=SRC/"lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json"
PROTO_SHA1="0c81483bd9c70daf7cc6341f35a857761ab964d6"
E4_SHA="ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f"
E4_LOG_SHA="32d93ec0aae6214881ba7852dfe148d2e82a284de5389aeba1dd2c2ba06798a4"
E4_MAN_SHA1="8e73c473319963c2c1d8a1c9887226da889d2126"
E4_PROTO_SHA1="b9e040e026cb2e6f6c1ad896f556fa29e295a6a1"
A02_SHA="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
E3_SHA="1e52f95be973aa6616337c2427fde60725a0df31ce7d3a818e99d1dd1fb0fb35"
USER_SHA1="203c6550b320c731576359c210fa045210edb9c3"
E01_SHA1="35cb92b3ce2b870ae6d984f0e9692266b3e025ac"
DESICOV_SHA1="e1c11ab41fe459fa261c59e731f9bffa7b0ccd42"
DESIRESULT_SHA1="4babbd9e23ba541a29caf5fb8fa5c3f22bca8388"
IDS=(1,125,250,375,500,625,750,875,1000)
CAPS=("NGC","SGC")
LEVELS=("full_1200_original","nested_2400","nested_4800")
PASS_E4="A03E4_NINE_MOCK_NESTED_1200_2400_4800_RANDOM_DENSITY_DESCRIPTIVE_ONLY"

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def read(p):
    return json.loads(p.read_bytes())

def check(ok,why):
    if not ok:raise ValueError(why)

def immutable():
    checks=((PROTO,"blob",PROTO_SHA1,None),
        (REPORT,"sha",E4_SHA,444242),(LOG,"sha",E4_LOG_SHA,11200),
        (MANIFEST,"blob",E4_MAN_SHA1,None),(E4PROTO,"blob",E4_PROTO_SHA1,None),
        (A02,"sha",A02_SHA,None),(E3,"sha",E3_SHA,None),
        (USER,"blob",USER_SHA1,None),(E01,"blob",E01_SHA1,None),
        (DESICOV,"blob",DESICOV_SHA1,None),(DESIRESULT,"blob",DESIRESULT_SHA1,None))
    for path,kind,expected,n in checks:
        raw=path.read_bytes()
        check((n is None or len(raw)==n) and
            (blob(raw) if kind=="blob" else sha(raw))==expected,
            "Frozen original source SHA/size changed "+str(path))
    return (read(x) for x in (PROTO,REPORT,MANIFEST,E4PROTO,USER,E01,DESICOV,DESIRESULT))

def audit(p,e4,m,e4p,user,e01,desicov,desiresult):
    parent=p["parents"];a03=p["A03_E5_prospective_requirements"];a04=p["A04_eboss_prospective_requirements"]
    check(p["date"]=="2026-09-27" and
        p["scope"]=="PLANNING_AND_INDEPENDENT_SOURCE_ONLY_READINESS_NOT_NEW_FITS_OR_INFERENCE" and
        p["registration_honesty"].startswith("E4 already measured;") and
        p["current_decision"]=="STOP_A03_PHYSICAL_AND_A04_EBOSS_INFERENCE_OBSERVED_ODD_SEALED" and
        p["branch"]=="audit/eboss-elg-bit8-ra-orientation-20260925" and
        p["pr_number"]==1,
        "Original after-E4 planning vs prospective science distinction lost")
    for k,actual in (
        ("E4_report_sha256",E4_SHA),("E4_report_exact_bytes",444242),
        ("E4_report_git_blob_sha1","32246a495506deac5b9dfe71bea51c3db33202b1"),
        ("E4_stdout_sha256",E4_LOG_SHA),("E4_stdout_exact_bytes",11200),
        ("E4_stdout_git_blob_sha1","31eb941a5a0fb65a8efca2a25530099b0bcc057d"),
        ("E4_manifest_git_blob_sha1",E4_MAN_SHA1),
        ("E4_prereg_git_blob_sha1",E4_PROTO_SHA1),
        ("A02_report_sha256",A02_SHA),("E3_report_sha256",E3_SHA),
        ("A03_user_decision_git_blob_sha1",USER_SHA1),
        ("A03_E0E1_report_git_blob_sha1",E01_SHA1),
        ("DESI_120mock_covariance_git_blob_sha1",DESICOV_SHA1),
        ("DESI_windowed_result_git_blob_sha1",DESIRESULT_SHA1)):
        check(parent[k]==actual,"Retrospective frozen parent identity altered: "+k)
    for flag in ("observed_galaxy_rows_read","observed_random_rows_read",
        "observed_odd_data_vector_read","new_science_selection_applied",
        "empirical_pair_window_physically_validated","EBOSS_18D_covariance_calibrated",
        "author_contact_approved","new_download_authorized",
        "main_mutation_permitted","unblinding_authorized"):
        check(p[flag] is False,"A03E5/A04 STOP guard changed "+flag)
    check(a03["status"].startswith("BLOCKED_") and
        a03["numeric_error_budget_now"]==
        "NOT_DEFINED_UNTIL_FINAL_OBSERVABLE_THEORY_AND_INDEPENDENT_COVARIANCE_FROZEN" and
        a03["new_random_density_or_seed_run_now"]=="NOT_AUTHORIZED" and
        a03["physical_pair_window_certified_now"] is False and
        len(a03["must_freeze_BEFORE_future_physical_measurement"])==5,
        "Physics-derived numerical budget and preregistration guard unexpectedly bypassed")
    check(a04["status"].startswith("BLOCKED_") and
        a04["eboss_final_18D"].startswith("NOT_YET_DEFINED;") and
        a04["DESI_18D"]=="SEPARATE_DESI_DR1_3Zx6S_DIPOLE_NOT_TRANSFERABLE_TO_EBOSS" and
        a04["pilot_geometry"]["redshift"]==[.9,1.] and
        a04["pilot_geometry"]["s_edges_mpc_h"]==[20,40,60,80,100,120,140] and
        a04["pilot_geometry"]["odd_ells"]==[1,3] and
        a04["pilot_geometry"]["caps"]==["NGC","SGC"] and
        a04["pilot_geometry"]["components_per_cap"]==12 and
        a04["pilot_geometry"]["joint_components"]==24 and
        a04["current_mock_cohort"]["distinct_mock_realizations"]==9 and
        a04["current_mock_cohort"]["covariance_centered_rank_upper_bound"]==8 and
        a04["current_mock_cohort"]["can_invert_p12"] is False and
        a04["current_mock_cohort"]["can_invert_p18"] is False and
        a04["current_mock_cohort"]["can_invert_p24"] is False and
        a04["large_download_authorized"] is False,
        "eBOSS observable dimension, covariance rank or DESI-isolation gate changed")
    check(desicov["scope"]=="DESI DR1 LRG–ELG covariance from 120 mock realizations" and
        desicov["n_mocks"]==120 and desicov["primary_vector"]=="18D dipole" and
        desicov["frozen_grid"]["z_edges"]==[.8,.9,1.,1.1] and
        len(desicov["frozen_grid"]["separation_centers_Mpc_over_h"])==6 and
        desiresult["scope"]=="Final RR-window-convolved DESI DR1 LRGxELG 120-mock finite-ensemble calibration",
        "External DESI-only 18D/120-mock result improperly relabelled")
    check(user["decision"].startswith("Do not contact Raichoor") and
        user["author_contact_approved"] is False and
        user["email_sent"] is False and
        user["additional_mock_sources_download_authorized"] is False and
        user["observed_odd_data_vector_read"] is False and
        e01["completed_cases"]==18 and
        e01["observed_odd_data_vector_read"] is False and
        e01["physical_pair_window_certified"] is False,
        "Prior user empirical-only no-author-contact or E0/E1 physical limitation changed")
    check(e4["status"]==PASS_E4 and e4["completed_cases"]==18 and
        e4["failed_cases"]==0 and e4["errors"]==[] and
        e4["fixed_ids"]==list(IDS) and e4["caps"]==list(CAPS) and
        e4["random_levels"]==list(LEVELS) and
        e4["observed_odd_data_vector_read"] is False and
        e4["physical_empirical_pair_window_certified"] is False and
        e4["inferential_18D_covariance_computed"] is False and
        e4p["fixed_cases"]==18 and
        e4p["physical_empirical_pair_window_certified"] is False,
        "Original frozen E4 mock-only scope or original prereg changed")
    bycap={c:[] for c in CAPS}
    count=0
    for mid in IDS:
        for cap in CAPS:
            key=f"{mid:04d}/{cap}";case=e4["cases"][key]
            check(case["status"]=="complete" and
                case["no_observed_odd_data_read"] is True and
                case["not_a_physical_window_or_inferential_covariance"] is True,
                "Original fixed E4 case/scope not complete "+key)
            for level in LEVELS:
                record=case["levels"][level]["baseline_original"]
                check(record["RR_supported_cells"]==144 and
                    record["no_unsupported_RR_zero_fill_or_projection"] is True,
                    "Original E4 RR-support evidence changed "+key+"/"+level)
                count+=1
            v=case["pairwise_common_support_density_comparisons"]["nested_4800_vs_nested_2400"]
            check(v["common_supported_cells"]==144 and
                math.isfinite(v["relative_L1_xi_shift_on_common_support"]) and
                v["relative_L1_xi_shift_on_common_support"]>=0,
                "Original E4 preregistered common-support L1 invalid "+key)
            bycap[cap].append(v["relative_L1_xi_shift_on_common_support"])
    check(count==54 and len(e4["cases"])==18 and
        len(e4["per_case_checkpoint_SHA256"])==18 and
        m["full_RR_baseline_scenario_records"]==54 and
        m["source_only_full_gzip_rehash_independently_repeated_outside_WSL"] is False and
        m["not_a_detection_or_exclusion"] is True,
        "Frozen E4 upload/manifest independent audit scope changed")
    for cap in CAPS:
        n=statistics.median(bycap[cap])
        s=p["E4_retrospective_status"]["relative_L1_xi_median_2400_to_4800"][cap]
        t=m["per_cap_descriptive_summaries"][cap]["relative_L1_xi_shift_on_common_support_by_density_comparison"]["nested_4800_vs_nested_2400"]["median"]
        check(math.isclose(n,s,rel_tol=0,abs_tol=1e-14) and
            math.isclose(n,t,rel_tol=0,abs_tol=1e-14),
            "E4 original nonconvergence output was sanitized or changed "+cap)
    check(p["E4_retrospective_status"]["finite_random_convergence"]=="NOT_DEMONSTRATED" and
        p["E4_retrospective_status"]["RR_full_support_certifies_estimator_convergence"] is False if
        "RR_full_support_certifies_estimator_convergence" in p["E4_retrospective_status"] else
        p["E4_retrospective_status"]["finite_random_convergence"]=="NOT_DEMONSTRATED",
        "Finite random nonconvergence was misreported")
    return {c:statistics.median(bycap[c]) for c in CAPS}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    args_loaded=tuple(immutable())
    med=audit(*args_loaded)
    if args.self_test:
        p,e4,m,e4p,user,e01,desicov,desiresult=args_loaded
        for field in ("observed_odd_data_vector_read","unblinding_authorized"):
            q=copy.deepcopy(p);q[field]=True
            try:audit(q,e4,m,e4p,user,e01,desicov,desiresult)
            except ValueError:pass
            else:raise AssertionError("Observed unblinding/odd negative control did not STOP "+field)
        wrong=copy.deepcopy(p);wrong["A04_eboss_prospective_requirements"]["pilot_geometry"]["joint_components"]=18
        try:audit(wrong,e4,m,e4p,user,e01,desicov,desiresult)
        except ValueError:pass
        else:raise AssertionError("Illegally reduced eBOSS pilot from 24D to 18D")
        wrong=copy.deepcopy(p);wrong["E4_retrospective_status"]["finite_random_convergence"]="DEMONSTRATED"
        try:audit(wrong,e4,m,e4p,user,e01,desicov,desiresult)
        except ValueError:pass
        else:raise AssertionError("Posthoc E4 convergence fiction accepted")
        wrong=copy.deepcopy(desicov);wrong["scope"]="eBOSS DR16"
        try:audit(p,e4,m,e4p,user,e01,wrong,desiresult)
        except ValueError:pass
        else:raise AssertionError("DESI covariance relabelled eBOSS")
    print("A03E5_A04_SOURCE_ONLY_READINESS_PARENT_AND_NEGATIVE_CONTROL_AUDIT_OK",flush=True)
    print("E4_ORIGINAL_FINITE_RANDOM_L1_2400_TO_4800_MEDIAN",med,flush=True)
    print("PHYSICAL_A03_STOP E_BOSS_18D_A04_STOP OBSERVED_ODD_SEALED",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
