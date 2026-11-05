#!/usr/bin/env python3
"""E4 original-upload SOURCE-ONLY independent archive audit; NO FITS or observed data.

All uploaded original report/log bytes are immutable. Recompute preregistered
PCG64 complements using only frozen A02 eligible-source counts; compare all
original A02/E2/E3 sample/pair/xi SHA and full 72-source stdout SHA evidence.
The WSL, not this source-only CI, rehashed original 72 local gzip files.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import statistics
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/"source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json"
LOG=ROOT/"source_data/eboss_dr16_a03_e4_local_stdout_2026-09-27.log"
MANIFEST=ROOT/"source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json"
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e4_nested_1200_2400_4800_mock_random_density_protocol_2026-09-27.json"
RUNNER=ROOT/"scripts/audit_eboss_dr16_a03_e4_nested_random_density.py"
A02=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
E2=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json"
E3=ROOT/"source_data/eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json"
E3MAN=ROOT/"source_data/eboss_dr16_a03_e3_mock_galaxy_nested_random_split_synthetic_dd_uploaded_manifest_2026-09-26.json"
G_SOURCE=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_report_2026-09-26.json"
R_SOURCE=ROOT/"source_data/eboss_dr16_nine_ezmock_matched_random_sha_report_2026-09-26.json"
REPORT_SHA="ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f"
LOG_SHA="32d93ec0aae6214881ba7852dfe148d2e82a284de5389aeba1dd2c2ba06798a4"
E3_SHA="1e52f95be973aa6616337c2427fde60725a0df31ce7d3a818e99d1dd1fb0fb35"
A02_SHA="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
E2_SHA="6ef86f904bb0cdc407adabe2e125146ce4d0671ea1ac932e13d11fc481cd596a"
PROTO_BLOB="b9e040e026cb2e6f6c1ad896f556fa29e295a6a1"
RUNNER_BLOB="fda4a5eda420a9f69bd4d393bae02419155958f8"
E3MAN_BLOB="3674a0345e22391143fd7b3c10a83a871e0b72cf"
IDS=(1,125,250,375,500,625,750,875,1000)
CAPS=("NGC","SGC")
TRACERS=("eBOSS_LRG","eBOSS_ELG")
ROLES=("dat","ran")
LEVELS=("full_1200_original","nested_2400","nested_4800")
TERMS=("D1D2","D1R2","R1D2","R1R2")
PASS="A03E4_NINE_MOCK_NESTED_1200_2400_4800_RANDOM_DENSITY_DESCRIPTIVE_ONLY"
FIXED_KEYS=tuple(f"{mid:04d}/{cap}" for mid in IDS for cap in CAPS)
SOURCE_LINE=re.compile(r"^A02_REHASH_72_INPUT_SHA_OK (\d{4})/(NGC|SGC)/(eBOSS_LRG|eBOSS_ELG) (dat|ran) (\d+) ([0-9a-f]{64})$")
SHA_PATTERN=re.compile(r"^[a-f0-9]{64}$")

def require(condition,message):
    if not condition:
        raise ValueError(message)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode("ascii")+b"\0"+raw).hexdigest()

def arr_sha(a):
    return sha(np.ascontiguousarray(a).tobytes())

def load(p):
    return json.loads(p.read_bytes())

def exact_uploaded_parent_gate():
    for p,n,expected,digest in ((REPORT,444242,REPORT_SHA,"32246a495506deac5b9dfe71bea51c3db33202b1"),
                                (LOG,11200,LOG_SHA,"31eb941a5a0fb65a8efca2a25530099b0bcc057d")):
        b=p.read_bytes()
        require(len(b)==n and sha(b)==expected and blob(b)==digest,
                "Exact original uploaded E4 report or WSL stdout archive mutated: "+str(p))
    for p,expected in ((A02,A02_SHA),(E2,E2_SHA),(E3,E3_SHA)):
        require(sha(p.read_bytes())==expected,"Exact parent A02/E2/E3 archive bytes changed: "+str(p))
    require(blob(PROTOCOL.read_bytes())==PROTO_BLOB and
            sha(PROTOCOL.read_bytes())==load(REPORT)["protocol_sha256"],
            "E4 preregistered protocol Git blob/full SHA changed")
    require(blob(RUNNER.read_bytes())==RUNNER_BLOB,
            "E4 original executable runner Git blob changed")
    require(blob(E3MAN.read_bytes())==E3MAN_BLOB,
            "Original E3 uploaded manifest Git blob changed")
    manifest=load(MANIFEST)
    require(manifest["original_uploaded_sha256"]==REPORT_SHA and
            manifest["original_uploaded_git_blob_sha1"]=="32246a495506deac5b9dfe71bea51c3db33202b1" and
            manifest["original_uploaded_exact_bytes"]==444242 and
            manifest["original_stdout_sha256"]==LOG_SHA and
            manifest["original_stdout_git_blob_sha1"]=="31eb941a5a0fb65a8efca2a25530099b0bcc057d" and
            manifest["original_stdout_exact_bytes"]==11200 and
            manifest["original_E4_protocol_git_blob_sha1"]==PROTO_BLOB and
            manifest["E4_runner_git_blob_sha1"]==RUNNER_BLOB and
            manifest["parent_E3_report_sha256"]==E3_SHA and
            manifest["parent_A02_report_sha256"]==A02_SHA and
            manifest["parent_A03E2_report_sha256"]==E2_SHA and
            manifest["source_only_full_gzip_rehash_independently_repeated_outside_WSL"] is False and
            manifest["observed_odd_data_vector_read"] is False and
            manifest["new_mock_downloads"] is False,
            "Original E4 independent manifest SHA/scope drift")
    return manifest

def check_log(lines,g,r):
    require(len(lines)==95,"Original WSL log line count changed")
    count=0;size=0
    for mid in IDS:
        for cap in CAPS:
            for tracer in TRACERS:
                for role in ROLES:
                    row=lines[count]
                    m=SOURCE_LINE.fullmatch(row)
                    key=f"{mid:04d}/{cap}/{tracer}"
                    require(m is not None and m.group(1)==f"{mid:04d}" and
                        m.group(2)==cap and m.group(3)==tracer and m.group(4)==role,
                        "Wrong original pre-FITS 72 source SHA row/order "+str(count))
                    parent=(g["existing_frozen_0001_sources"][key] if mid==1 else
                            g["new_galaxy_sources"][key]) if role=="dat" else r["verified_random_sources"][key]
                    expected=(parent["full_compressed_sha256_reverified"] if mid==1 else
                              parent["full_compressed_sha256_first_seen"]) if role=="dat" else parent["fully_reverified_compressed_sha256"]
                    require(int(m.group(5))==parent["compressed_bytes"] and m.group(6)==expected,
                            "Original full gzip SHA or bytes differ from pinned archived source: "+key+"/"+role)
                    size+=int(m.group(5));count+=1
    require(count==72 and size==1860198719 and
        lines[72]=="A02_ALL_72_COMPLETE_GZIP_SHA_OK 72 1860198719",
        "Pre-FITS 72 original compressed source identity not complete")
    for i,key in enumerate(FIXED_KEYS,73):
        require(lines[i]==("A03E4_FIXED_MOCK_CASE "+key+
            " RR_SUPPORT {'full_1200_original': 144, 'nested_2400': 144, 'nested_4800': 144}"),
            "Original E4 fixed case order/full RR console status changed "+key)
    require(lines[91]=="A03E4_NESTED_1200_2400_4800_MOCK_ONLY "+PASS and
        lines[93]=="COMPLETED_CASES 18" and lines[94]=="OBSERVED_ODD_DATA_READ False",
        "Original local E4 terminal status/observed-odd guard changed")
    return count,size

def audit_model(j,lines,manifest,a02,e2,e3,e3man,g,r):
    require(j["status"]==PASS and j["completed_cases"]==18 and
        j["failed_cases"]==0 and j["errors"]==[] and
        j["fixed_ids"]==list(IDS) and j["caps"]==list(CAPS) and
        j["random_levels"]==list(LEVELS) and
        list(j["cases"])==list(FIXED_KEYS) and
        list(j["per_case_checkpoint_SHA256"])==list(FIXED_KEYS),
        "E4 fixed 18-case PASS/order/checkpoint scope changed")
    require(j["parent_E3_archive_sha256"]==E3_SHA and
        j["parent_A02_archive_sha256"]==A02_SHA and
        j["total_original_72_compressed_bytes"]==1860198719 and
        j["all_72_full_gzip_SHA_rehashed_before_any_E4_FITS_row"] is True and
        j["observed_galaxy_rows_read"] is False and
        j["observed_random_rows_read"] is False and
        j["observed_odd_data_vector_read"] is False and
        j["new_science_selection_applied"] is False and
        j["physical_empirical_pair_window_certified"] is False and
        j["inferential_18D_covariance_computed"] is False and
        j["not_a_detection_or_exclusion"] is True,
        "E4 physical scope/source preflight/observed seal changed")
    require(e3["status"]=="A03E3_NINE_MOCK_NESTED_RANDOM_AND_SYNTHETIC_DD_ODD_DESCRIPTIVE_ONLY"
        and e3["errors"]==[] and e3["completed_cases"]==18 and
        list(e3["cases"])==list(FIXED_KEYS) and
        e3man["original_uploaded_sha256"]==E3_SHA and
        e2["completed_cases"]==18 and e2["errors"]==[] and
        a02["completed_cases"]==18 and a02["errors"]==[],
        "Original A02/E2/E3 18-case parent scope no longer closed")
    count,size=check_log(lines,g,r)
    metrics={cap:{x:[] for x in ("nested_2400_vs_full_1200_original",
        "nested_4800_vs_nested_2400","nested_4800_vs_full_1200_original")} for cap in CAPS}
    pairchecks=0
    for key in FIXED_KEYS:
        case=j["cases"][key];a=a02["cases"][key];e=e2["cases"][key];old=e3["cases"][key];m=manifest["per_fixed_case_source_and_output_fingerprints"][key]
        mid=int(key[:4]);cap=key[5:]
        raw=json.dumps(case,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
        require(case["status"]=="complete" and case["id"]==mid and case["cap"]==cap and
            case["original_A02_exact_input_metadata_replayed"] is True and
            case["parent_E3_original_1200_random_catalogue_SHA_replayed"] is True and
            sha(raw)==j["per_case_checkpoint_SHA256"][key]==m["per_case_checkpoint_SHA256"],
            "E4 original per-case checkpoint/parent flag changed "+key)
        for ti,tracer in enumerate(TRACERS):
            for role in ROLES:
                field=tracer+"_"+role;actual=case["original_A02_all_four_selected_catalogue_SHA256"][field]
                require(actual==a["input_sample_diagnostics"][field]["selected_array_SHA256"]
                    ==e["original_selected_array_SHA256"][field]
                    ==m["original_A02_all_four_selected_catalogue_SHA256"][field],
                    "Original A02/E2 exact fixed sample SHA changed "+key+"/"+field)
            ext=case["fixed_random_extensions_by_tracer"][tracer]
            oldran=a["input_sample_diagnostics"][tracer+"_ran"]
            eligible=ext["eligible_after_original_A02_gate"]
            require(eligible==oldran["eligible_highz_rows_after_fixed_weight_gate"] and eligible>=4800,
                "Original A02 random eligible count changed "+key+"/"+tracer)
            seed=93127+10000*CAPS.index(cap)+100*ti+1
            original=np.sort(np.random.default_rng(seed).choice(eligible,size=1200,replace=False))
            xseed=20502027+100000*mid+1000*CAPS.index(cap)+100*ti
            avail=np.setdiff1d(np.arange(eligible,dtype=np.int64),original,assume_unique=True)
            extra=avail[np.random.default_rng(xseed).choice(len(avail),size=3600,replace=False)]
            selected=(original,np.sort(np.concatenate((original,extra[:1200]))),
                      np.sort(np.concatenate((original,extra))))
            require(ext["fixed_supplement_PCG64_seed"]==xseed and
                ext["original_1200_eligible_positions_SHA256"]==arr_sha(original) and
                ext["supplement_3600_UNSORTED_rng_order_SHA256"]==arr_sha(extra) and
                ext["supplement_disjoint_from_original_1200"] is True and
                ext["nested_1200_within_2400_within_4800"] is True and
                ext["fixed_supplement_PCG64_seed"]==m["per_tracer_extension_fingerprints"][tracer]["fixed_supplement_PCG64_seed"],
                "Independent PCG64 frozen extra seed/original index SHA mismatch "+key+"/"+tracer)
            for level,num,idx in zip(LEVELS,(1200,2400,4800),selected):
                t=ext["levels"][level]
                require(t["n"]==num and t["selected_eligible_positions_SHA256"]==arr_sha(idx) and
                    t["selected_eligible_positions_SHA256"]==
                    m["per_tracer_extension_fingerprints"][tracer]["selected_eligible_positions_SHA256_by_level"][level] and
                    t["four_vector_catalogue_SHA256"]==
                    m["per_tracer_extension_fingerprints"][tracer]["four_vector_catalogue_SHA256_by_level"][level] and
                    SHA_PATTERN.fullmatch(t["selected_original_FITS_row_indices_SHA256"]) is not None,
                    "Fixed nested 1200/2400/4800 source index SHA or catalogue digest mismatch "+key+"/"+tracer+"/"+level)
                if tracer=="eBOSS_ELG":
                    require(sum(t["exact_selected_ELG_chunk_counts"].values())==num,
                        "Full ELG exact chunk count lost "+key+"/"+level)
            require(ext["levels"]["full_1200_original"]["four_vector_catalogue_SHA256"]==
                old["random_permutation_and_levels"][tracer]["per_level"]["full_1200"]["subcatalogue_four_vector_SHA256"],
                "Original E3 full 1200 random vector SHA not reproduced "+key+"/"+tracer)
        original_xi=e["scenarios"]["baseline_original"]["forward_xi_grid_sha256"]
        dd_ref=None
        for level in LEVELS:
            entry=case["levels"][level];rec=entry["baseline_original"]
            require(rec["RR_supported_cells"]==144 and rec["RR_missing_indices"]==[] and
                rec["status"]=="complete_full_RR_support" and
                rec["original_physical_domain_multipoles"] is not None and
                rec["no_unsupported_RR_zero_fill_or_projection"] is True and
                rec["max_abs_xi_reverse_on_supported_cells"]<1e-8 and
                entry["LRG_R_rows"]==entry["ELG_R_rows"]==
                {"full_1200_original":1200,"nested_2400":2400,"nested_4800":4800}[level] and
                rec["full_6x24_xi_SHA256"]==m["baseline_full_6x24_xi_SHA256_by_density"][level] and
                m["RR_supported_cells_by_density"][level]==144,
                "Original full RR/full multipole/forward reverse closure changed "+key+"/"+level)
            require(rec["RR_support_mask_sha256"]==
                arr_sha(np.ones((6,24),dtype=np.uint8)),
                "Original positive RR support mask SHA changed "+key+"/"+level)
            for side in ("forward_pair_terms","reverse_pair_terms"):
                for term in TERMS:
                    t=rec[side][term]
                    require(t["normalization_relative_residual"]<1e-12 and
                            t["accepted_pairs"]>=0 and
                            SHA_PATTERN.fullmatch(t["weighted_histogram_SHA256"]) is not None,
                            "E4 pair normalizations/weighted SHA invalid "+key+"/"+level+"/"+side+"/"+term)
                    if level=="full_1200_original":
                        parent=e["scenarios"]["baseline_original"][side][term]
                        require(t["weighted_histogram_SHA256"]==parent["weighted_histogram_SHA256"] and
                            t["accepted_pairs"]==parent["accepted_pairs"] and
                            math.isclose(t["independently_normalized_pair_weight"],
                                         parent["independently_normalized_pair_weight"],
                                         rel_tol=1e-12,abs_tol=1e-12),
                            "Original E2 exact full 1200 weighted pair SHA/count/norm changed "+key+"/"+side+"/"+term)
                    pairchecks+=1
            require(rec["forward_pair_terms"]["D1D2"]["accepted_pairs"]==
                rec["reverse_pair_terms"]["D1D2"]["accepted_pairs"] and
                rec["forward_pair_terms"]["D1R2"]["accepted_pairs"]==
                rec["reverse_pair_terms"]["R1D2"]["accepted_pairs"] and
                rec["forward_pair_terms"]["R1D2"]["accepted_pairs"]==
                rec["reverse_pair_terms"]["D1R2"]["accepted_pairs"] and
                rec["forward_pair_terms"]["R1R2"]["accepted_pairs"]==
                rec["reverse_pair_terms"]["R1R2"]["accepted_pairs"],
                "True forward/reverse accepted cross-pair count mismatch "+key+"/"+level)
            dd=(rec["forward_pair_terms"]["D1D2"],rec["reverse_pair_terms"]["D1D2"])
            if dd_ref is not None:
                require(dd==dd_ref,"Original DD changed when only random density changed "+key+"/"+level)
            dd_ref=dd
            for ell in ("0","1","2","3"):
                v=rec["original_physical_domain_multipoles"][ell]["values_by_fixed_s_bin"]
                require(len(v)==6 and all(math.isfinite(x) for x in v),
                        "Full 144 RR multipoles not finite 6 bins "+key+"/"+level+"/"+ell)
            if level=="full_1200_original":
                require(rec["full_6x24_xi_SHA256"]==original_xi==
                    old["levels"]["full_1200"]["scenarios"]["baseline_original"]["full_6x24_xi_SHA256"]==
                    e3man["per_fixed_case_source_and_output_fingerprints"][key]
                    ["original_E2_full_1200_xi_SHA256_by_scenario"]["baseline_original"],
                    "E4 full 1200 baseline E2/E3 exact xi SHA changed "+key)
        for name in metrics[cap]:
            v=case["pairwise_common_support_density_comparisons"][name]
            require(v["common_supported_cells"]==144 and v["no_common_support_warning"] is False and
                    all(math.isfinite(v[k]) and v[k]>=0 for k in
                        ("max_abs_xi_difference_on_common_cells",
                         "median_abs_xi_difference_on_common_cells",
                         "relative_L1_xi_shift_on_common_support")),
                    "Predeclared xi common support density diagnostic is incomplete "+key+"/"+name)
            metrics[cap][name].append(v["relative_L1_xi_shift_on_common_support"])
    require(pairchecks==432,"18x3x8 original pair checks incomplete")
    for cap in CAPS:
        p=manifest["per_cap_descriptive_summaries"][cap]
        require(p["case_count"]==9 and p["RR_support_all_three_levels_in_all_nine_cases"] is True,
                "Original cap/coverage summary changed "+cap)
        for name,vals in metrics[cap].items():
            s=p["relative_L1_xi_shift_on_common_support_by_density_comparison"][name]
            require(s["min"]==min(vals) and s["median"]==statistics.median(vals) and
                    s["max"]==max(vals),
                    "Manifest predeclared per-cap all-nine L1 summaries disagree "+cap+"/"+name)
    return count,pairchecks

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    m=exact_uploaded_parent_gate()
    data=load(REPORT);a02=load(A02);e2=load(E2);e3=load(E3);e3man=load(E3MAN)
    g=load(G_SOURCE);r=load(R_SOURCE)
    lines=LOG.read_bytes().decode("ascii").splitlines()
    source_count,pairs=audit_model(data,lines,m,a02,e2,e3,e3man,g,r)
    if args.self_test:
        for field in ("observed_odd_data_vector_read","physical_empirical_pair_window_certified"):
            bad=copy.deepcopy(data);bad[field]=True
            try:audit_model(bad,lines,m,a02,e2,e3,e3man,g,r)
            except ValueError:pass
            else:raise AssertionError("Negative observed/physical-scope control unexpectedly passed "+field)
        tampered=copy.deepcopy(data)
        tampered["cases"]["0001/NGC"]["original_A02_all_four_selected_catalogue_SHA256"]["eBOSS_ELG_ran"]="0"*64
        try:audit_model(tampered,lines,m,a02,e2,e3,e3man,g,r)
        except ValueError:pass
        else:raise AssertionError("Tampered original A02 sample SHA unexpectedly accepted")
        wrong=list(lines);wrong[0]=wrong[0].replace("67862ecb","00000000",1)
        try:audit_model(data,wrong,m,a02,e2,e3,e3man,g,r)
        except ValueError:pass
        else:raise AssertionError("Tampered original full gzip source SHA line unexpectedly accepted")
    require(m["full_RR_baseline_scenario_records"]==54 and
            m["independent_numpy_PCG64_all_original_and_extra_seeds_and_1200_2400_4800_index_SHA_checked"]==36 and
            m["source_72_log_entries_independently_matched_archived_galaxy_random_SHA_and_bytes"]==source_count and
            m["all_forward_reverse_pair_terms_checked"]==pairs,
            "Independent uploaded E4 manifest expected-source audit totals changed")
    print("A03E4_ORIGINAL_UPLOADED_ARCHIVE_SOURCE_ONLY_18_CASES_OK",
          "72_ORIGINAL_SOURCE_SHA_LINES",source_count,
          "36_EXACT_PCG64_EXTENSION_REPLAYS","54_FULL_RR_CASE_LEVELS",
          "432_FORWARD_REVERSE_PAIR_TERMS",pairs,flush=True)
    print("E4_INTERPRETATION_FINITE_RANDOM_NOT_CONVERGED_PHYSICAL_A03_AND_A04_OPEN_ODD_SEALED",flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
