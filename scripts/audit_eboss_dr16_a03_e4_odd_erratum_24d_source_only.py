#!/usr/bin/env python3
"""Independent source-only validation of the separate E4 ell=1,3 erratum.

Original E4 JSON and its derivative original manifest are IMMUTABLE. Recompute
18 x 2 ell x 3 density comparisons from archived six-bin multipoles, both
corrected per-cap nine-mock summaries, and the descriptive 24D joint pilot
drift. This is RETROSPECTIVE; no FITS, observed rows, new randoms or inference.
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
ORIGINAL=SRC/"eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json"
MANIFEST=SRC/"eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json"
ERRATUM=SRC/"eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json"
ORIGINAL_SHA256="ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f"
ORIGINAL_BLOB="32246a495506deac5b9dfe71bea51c3db33202b1"
MANIFEST_BLOB="8e73c473319963c2c1d8a1c9887226da889d2126"
ERRATUM_BLOB="0baf9c5a009bf8f5e78f5725e44395bad6f91715"
IDS=(1,125,250,375,500,625,750,875,1000)
CAPS=("NGC","SGC")
ELLS=("1","3")
LEVELS=("full_1200_original","nested_2400","nested_4800")
COMPARISONS=(("nested_2400","full_1200_original"),
             ("nested_4800","nested_2400"),
             ("nested_4800","full_1200_original"))
S_CENTERS=(30,50,70,90,110,130)
PASS="A03E4_NINE_MOCK_NESTED_1200_2400_4800_RANDOM_DENSITY_DESCRIPTIVE_ONLY"

def require(ok,message):
    if not ok:raise ValueError(message)

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def load(p):
    return json.loads(p.read_bytes())

def close(a,b,label):
    require(isinstance(a,(int,float)) and isinstance(b,(int,float)) and
        math.isfinite(a) and math.isfinite(b) and
        math.isclose(a,b,rel_tol=1e-13,abs_tol=1e-13),
        "Corrected E4 odd diagnostics source-only mismatch: "+label)

def summary(values):
    require(len(values)==9 and all(math.isfinite(x) and x>=0 for x in values),
        "E4 all-nine summary needs nine nonnegative finite values")
    return {"min":min(values),"median":statistics.median(values),"max":max(values)}

def compare_summary(actual,expected,where):
    require(set(actual)=={"min","median","max"},where+" malformed summary keys")
    for k,v in expected.items():close(actual[k],v,where+"/"+k)

def audit(original,legacy,err):
    require(original["status"]==PASS and original["completed_cases"]==18 and
        original["failed_cases"]==0 and original["errors"]==[] and
        original["fixed_ids"]==list(IDS) and original["caps"]==list(CAPS) and
        original["random_levels"]==list(LEVELS) and
        original["observed_odd_data_vector_read"] is False and
        original["physical_empirical_pair_window_certified"] is False and
        original["inferential_18D_covariance_computed"] is False,
        "Original immutable E4 18-case and sealed scope failed")
    keys=tuple(f"{mid:04d}/{cap}" for mid in IDS for cap in CAPS)
    require(tuple(original["cases"])==keys and
        tuple(original["per_case_checkpoint_SHA256"])==keys and
        set(err["per_case"])==set(keys) and
        err["fixed_original_IDs"]==list(IDS) and
        err["caps"]==list(CAPS) and
        err["ells"]==list(ELLS) and
        err["s_centers_mpc_h"]==list(S_CENTERS) and
        err["pilot_joint_dim"]==24 and
        err["pilot_joint_vector_order"]==
            "NGC ell1 s=30..130; NGC ell3 s=30..130; SGC ell1 s=30..130; SGC ell3 s=30..130",
        "E4 erratum frozen original case/24D pilot ordering differs")
    require(err["original_report_sha256"]==ORIGINAL_SHA256 and
        err["original_report_git_blob_sha1"]==ORIGINAL_BLOB and
        err["original_unchanged_manifest_git_blob_sha1"]==MANIFEST_BLOB and
        err["after_the_original_E4_trial"] is True and
        err["not_prospectively_registered_E4_acceptance"] is True and
        err["scope"]=="RETROSPECTIVE_FINITE_RANDOM_MOCK_ONLY_NO_INFERENCE" and
        err["RR_supported_144_of_144_at_original_1200"] is True and
        err["original_E4_nonconvergence_not_reclassified"] is True and
        err["no_observed_rows_or_odd_read"] is True and
        err["no_new_source_or_download"] is True and
        err["no_new_random_seed_cut_or_threshold"] is True and
        err["DESI_120_mock_18D_covariance_not_used"] is True and
        err["physical_eBOSS_pair_window_certified"] is False and
        err["eBOSS_inferential_covariance_computed"] is False,
        "E4 erratum original fingerprints, timing or anti-inference flags changed")
    label_pairs={hi+"_minus_"+lo:(hi,lo) for hi,lo in COMPARISONS}
    require(err["comparisons"]==list(label_pairs),
        "E4 original three density comparisons changed")
    require(tuple(legacy["per_cap_descriptive_summaries"])==CAPS,
        "Original legacy E4 manifest cap list changed")
    for cap in CAPS:
        bad=legacy["per_cap_descriptive_summaries"][cap][
            "odd_ell_median_case_max_abs_difference_by_density_comparison"]
        require(bad.get("1")=={} and bad.get("3")=={} and
            set(bad)=={"1","3",*label_pairs},
            "Historical E4 derivative manifest bug differs; do not silently patch original")
    validated_vectors=0
    for mid in IDS:
        midstr=f"{mid:04d}"
        for cap in CAPS:
            key=midstr+"/"+cap
            case=original["cases"][key]
            corr=err["per_case"][key]
            require(case["status"]=="complete" and
                corr["original_E4_checkpoint_SHA256"]==
                original["per_case_checkpoint_SHA256"][key] and
                set(corr["by_ell_and_comparison"])==set(ELLS),
                "E4 original case/checkpoint source misalignment "+key)
            for ell in ELLS:
                require(set(corr["by_ell_and_comparison"][ell])==set(label_pairs),
                    "E4 corrected ell comparison keys absent "+key+"/"+ell)
                for label,(hi,lo) in label_pairs.items():
                    x=case["levels"][lo]["baseline_original"][
                        "original_physical_domain_multipoles"][ell]["values_by_fixed_s_bin"]
                    y=case["levels"][hi]["baseline_original"][
                        "original_physical_domain_multipoles"][ell]["values_by_fixed_s_bin"]
                    require(len(x)==len(y)==6 and all(
                        math.isfinite(z) for z in x+y),
                        "Original E4 archived ell values are non-finite or wrong grid "+key)
                    expected=[a-b for a,b in zip(y,x)]
                    q=corr["by_ell_and_comparison"][ell][label]
                    require(len(q["signed_six_s_bin_shift"])==6,
                        "Erratum E4 signed six-bin vector shape mismatch "+key)
                    for i,(a,b) in enumerate(zip(q["signed_six_s_bin_shift"],expected)):
                        close(a,b,key+"/"+ell+"/"+label+"/"+str(i))
                    absolute=sorted(abs(v) for v in expected)
                    close(q["max_abs_over_s"],max(absolute),
                          key+"/"+ell+"/"+label+"/maximum")
                    close(q["median_abs_over_s"],statistics.median(absolute),
                          key+"/"+ell+"/"+label+"/median")
                    validated_vectors+=1
            for level in LEVELS:
                require(case["levels"][level]["baseline_original"]["RR_supported_cells"]==144,
                        "Original E4 RR support changed "+key+"/"+level)
        for label in label_pairs:
            joint=[]
            for cap in CAPS:
                for ell in ELLS:
                    joint+=err["per_case"][midstr+"/"+cap][
                        "by_ell_and_comparison"][ell][label]["signed_six_s_bin_shift"]
            require(len(joint)==24,"Wrong 24D diagnostic dimension")
            q=err["joint_24D_drift_by_fixed_mock_ID"][midstr][label]
            require(len(q["signed_24D_pilot_shift"])==24,
                    "Erratum 24D vector missing "+midstr+"/"+label)
            for i,(a,b) in enumerate(zip(q["signed_24D_pilot_shift"],joint)):
                close(a,b,midstr+"/"+label+"/24D/"+str(i))
            close(q["maximum_abs_component"],max(abs(x) for x in joint),
                  midstr+"/"+label+"/24D max abs")
            close(q["euclidean_norm"],math.hypot(*joint),
                  midstr+"/"+label+"/24D L2")
    require(validated_vectors==18*2*3,
            "All 108 per-case ell/density six-bin vectors must be verified")
    for cap in CAPS:
        c=err["per_cap_corrected_odd_diagnostics"][cap]
        require(c["fixed_mock_count"]==9 and
            set(c["by_ell_and_comparison"])==set(ELLS),
            "Corrected cap does not contain both independently scoped ell values "+cap)
        for ell in ELLS:
            for label in label_pairs:
                rows=[err["per_case"][f"{mid:04d}/{cap}"][
                    "by_ell_and_comparison"][ell][label] for mid in IDS]
                result=c["by_ell_and_comparison"][ell][label]
                compare_summary(result["per_mock_max_abs_over_s"],
                    summary([v["max_abs_over_s"] for v in rows]),
                    cap+"/"+ell+"/"+label)
                expect=[statistics.median(
                    [abs(v["signed_six_s_bin_shift"][i]) for v in rows])
                    for i in range(6)]
                require(len(result["per_s_bin_median_abs_across_nine"])==6,
                    "Six separation bins missing from corrected median "+cap+"/"+ell)
                for i,(a,b) in enumerate(zip(
                    result["per_s_bin_median_abs_across_nine"],expect)):
                    close(a,b,cap+"/"+ell+"/"+label+"/s"+str(i))
                if ell=="3":
                    previous=legacy["per_cap_descriptive_summaries"][cap][
                        "odd_ell_median_case_max_abs_difference_by_density_comparison"][label]
                    compare_summary(previous,
                        result["per_mock_max_abs_over_s"],
                        "Legacy unscoped ell3 summary preserved "+cap+"/"+label)
    for label in label_pairs:
        rows=[err["joint_24D_drift_by_fixed_mock_ID"][f"{mid:04d}"][label]
              for mid in IDS]
        q=err["joint_24D_all_nine_descriptive"][label]
        compare_summary(q["maximum_abs_component"],
            summary([v["maximum_abs_component"] for v in rows]),"joint 24D Linf/"+label)
        compare_summary(q["euclidean_norm"],
            summary([v["euclidean_norm"] for v in rows]),"joint 24D L2/"+label)
    return validated_vectors

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    raw_original=ORIGINAL.read_bytes()
    raw_manifest=MANIFEST.read_bytes()
    raw_erratum=ERRATUM.read_bytes()
    require(hashlib.sha256(raw_original).hexdigest()==ORIGINAL_SHA256 and
        blob(raw_original)==ORIGINAL_BLOB and
        blob(raw_manifest)==MANIFEST_BLOB and
        blob(raw_erratum)==ERRATUM_BLOB,
        "Source-only E4 original report, manifest or independent erratum bytes mutated")
    original,manifest,erratum=(json.loads(raw) for raw in
        (raw_original,raw_manifest,raw_erratum))
    vectors=audit(original,manifest,erratum)
    if args.self_test:
        bad=copy.deepcopy(erratum)
        bad["per_cap_corrected_odd_diagnostics"]["NGC"]["by_ell_and_comparison"]["1"][
            "nested_4800_minus_nested_2400"]["per_mock_max_abs_over_s"]["median"]+=.001
        try:audit(original,manifest,bad)
        except ValueError:pass
        else:raise AssertionError("Tampered missing original ell1 replacement accepted")
        bad=copy.deepcopy(erratum)
        bad["per_case"]["0001/SGC"]["by_ell_and_comparison"]["3"][
            "nested_4800_minus_nested_2400"]["signed_six_s_bin_shift"][0]=0.
        try:audit(original,manifest,bad)
        except ValueError:pass
        else:raise AssertionError("Tampered per-case E4 ell3 signed source accepted")
        bad=copy.deepcopy(erratum)
        bad["pilot_joint_dim"]=18
        try:audit(original,manifest,bad)
        except ValueError:pass
        else:raise AssertionError("Illegitimate eBOSS 24D->18D diagnostic change accepted")
        bad=copy.deepcopy(erratum)
        bad["no_observed_rows_or_odd_read"]=False
        try:audit(original,manifest,bad)
        except ValueError:pass
        else:raise AssertionError("Tampered observed-odd source-only guard accepted")
    print("E4_ORIGINAL_MANIFEST_ELL1_TRUNCATION_ERRATUM_SOURCE_ONLY_PASS",
          vectors,"SIX_BIN_ELL_DENSITY_VECTORS",
          "9_JOINT_24D_PILOT_VECTORS_PER_COMPARISON",flush=True)
    print("NO_NEW_MOCK_OR_OBSERVED_FITS_PHYSICAL_A03_A04_OPEN_ODD_SEALED",
          flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
