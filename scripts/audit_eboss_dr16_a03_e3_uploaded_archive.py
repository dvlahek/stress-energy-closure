#!/usr/bin/env python3
"""Independent source-only A-03E3 original-upload archive audit; NO FITS/observed reads."""
from __future__ import annotations
import argparse
import ast
import copy
import hashlib
import json
import re
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent.parent
SOURCE=ROOT/"source_data"
REPORT=SOURCE/"eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json"
LOG=SOURCE/"eboss_dr16_a03_e3_local_stdout_2026-09-26.log"
MANIFEST=SOURCE/"eboss_dr16_a03_e3_mock_galaxy_nested_random_split_synthetic_dd_uploaded_manifest_2026-09-26.json"
PROTOCOL=SOURCE/"eboss_dr16_a03_e3_nested_mock_random_split_synthetic_dd_injection_protocol_2026-09-26.json"
E2=SOURCE/"eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json"
A02=SOURCE/"eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
REPORT_SHA="1e52f95be973aa6616337c2427fde60725a0df31ce7d3a818e99d1dd1fb0fb35"
REPORT_BLOB="fa32a0f1546307e7bcdfe3a7911e613b83bb4f53"
LOG_SHA="869d46781871f6d88b5f487618e6d1e0c823b9f35d72c0edf98b673f9684e1d6"
LOG_BLOB="6bb9266abd1aa00087d4150e6e89ac816b205ebe"
MANIFEST_BLOB="3674a0345e22391143fd7b3c10a83a871e0b72cf"
PROTOCOL_SHA="bb4cd60ff84a167a09d91324f8a2fce8948f15b86766ae82f4bed121df384a0d"
E2_SHA="6ef86f904bb0cdc407adabe2e125146ce4d0671ea1ac932e13d11fc481cd596a"
A02_SHA="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
IDS=(1,125,250,375,500,625,750,875,1000)
CAPS=("NGC","SGC")
TRACERS=("eBOSS_LRG","eBOSS_ELG")
LEVELS=("half_A_600","half_B_600","three_quarter_A_900","full_1200")
SCENARIOS=("baseline_original","plus_5pct_ELG_R_z_ramp","minus_5pct_ELG_R_z_ramp")
TERMS=("D1D2","D1R2","R1D2","R1R2")
MIRROR={"D1D2":"D1D2","D1R2":"R1D2","R1D2":"D1R2","R1R2":"R1R2"}

def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def arrsha(a):return sha(np.ascontiguousarray(a).tobytes())
def require(ok,msg):
    if not ok:raise ValueError(msg)
def obj(p):return json.loads(p.read_bytes())

def audit(j,m,e2,a02,lines):
    keys=[f"{mid:04d}/{cap}" for mid in IDS for cap in CAPS]
    require(j["status"]=="A03E3_NINE_MOCK_NESTED_RANDOM_AND_SYNTHETIC_DD_ODD_DESCRIPTIVE_ONLY"
        and j["completed_cases"]==18 and j["failed_cases"]==0 and j["errors"]==[]
        and list(j["cases"])==keys and list(j["per_case_checkpoint_SHA256"])==keys
        and j["protocol_sha256"]==PROTOCOL_SHA and j["parent_E2_archive_sha256"]==E2_SHA
        and j["parent_A02_archive_sha256"]==A02_SHA
        and j["fixed_ids"]==list(IDS) and j["caps"]==list(CAPS)
        and j["random_levels"]==list(LEVELS) and j["weight_scenarios"]==list(SCENARIOS)
        and j["DD_pair_injection_amplitudes"]==[.02,-.02]
        and j["all_72_full_gzip_SHA_rehashed_before_any_E3_FITS_row"] is True
        and j["total_original_72_compressed_bytes"]==1860198719
        and all(j[x] is False for x in (
            "observed_galaxy_rows_read","observed_random_rows_read","observed_odd_data_vector_read",
            "new_science_selection_applied","physical_empirical_pair_window_certified",
            "physical_wake_injection_recovery_certified","inferential_18D_covariance_computed"))
        and j["not_a_detection_or_exclusion"] is True,"E3 fixed IDs/scope/parent/report changed")
    require(e2["completed_cases"]==18 and e2["errors"]==[] and list(e2["cases"])==keys
        and a02["completed_cases"]==18 and a02["errors"]==[]
        and m["number_of_cases"]==18 and m["failed_cases"]==0
        and m["original_uploaded_exact_bytes"]==1848217 and m["original_uploaded_sha256"]==REPORT_SHA
        and m["original_uploaded_git_blob_sha1"]==REPORT_BLOB
        and m["archived_report_git_blob_sha1"]==REPORT_BLOB
        and m["archive_byte_identical_to_uploaded_source"] is True
        and m["original_stdout_sha256"]==LOG_SHA
        and m["original_stdout_git_blob_sha1"]==LOG_BLOB
        and m["stdout_archive_byte_identical_to_uploaded_source"] is True
        and m["registered_E3_protocol_sha256"]==PROTOCOL_SHA
        and m["exact_parent_E2_archive_sha256"]==E2_SHA
        and m["exact_parent_A02_archive_sha256"]==A02_SHA
        and list(m["per_fixed_case_source_and_output_fingerprints"])==keys,
        "Independent E3 original upload manifest or archived parent mismatch")
    require(len(lines)==95 and lines[72]=="A02_ALL_72_COMPLETE_GZIP_SHA_OK 72 1860198719",
            "All-72-SHA-before-FITS log sequence changed")
    totals=0
    for ix,line in enumerate(lines[:72]):
        expected_key=keys[ix//4]+"/"+TRACERS[(ix%4)//2]
        expected_role=("dat","ran")[ix%2]
        hit=re.fullmatch(r"A02_REHASH_72_INPUT_SHA_OK (\d{4}/(?:NGC|SGC)/(?:eBOSS_LRG|eBOSS_ELG)) (dat|ran) (\d+) ([0-9a-f]{64})",line)
        require(hit is not None and hit.group(1)==expected_key and hit.group(2)==expected_role,
                "Original 72-SHA log order or source identity changed")
        totals+=int(hit.group(3))
    require(totals==1860198719 and lines[91].startswith("A03E3_NESTED_MOCK_RANDOM_SPLIT_AND_DD_ODD_INJECTION ")
        and lines[93]=="COMPLETED_CASES 18" and lines[94]=="OBSERVED_ODD_DATA_READ False",
        "Log final completion or source byte totals differ")
    scenario_count=sparse=inject_count=0
    maxres={k:0. for k in ("forward_analytic","reverse_analytic","mirror_increment","mirror_injected_xi")}
    for key in keys:
        c=j["cases"][key];old=e2["cases"][key];prior=a02["cases"][key]
        f=m["per_fixed_case_source_and_output_fingerprints"][key]
        raw=json.dumps(c,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
        require(sha(raw)==j["per_case_checkpoint_SHA256"][key]==f["per_case_checkpoint_SHA256"],
                "Per-case JSON checkpoint bytes changed: "+key)
        require(c["status"]=="complete" and c["original_e2_selected_array_sha_replayed"] is True
            and c["full_original_e2_all_three_pair_and_xi_SHA_replayed"] is True
            and c["no_observed_rows_read"] is True and c["no_physical_wake_injection_or_covariance"] is True,
            "E3 parent/scope case changed: "+key)
        mid=int(key[:4]);cap=key[-3:]
        require(c["id"]==mid and c["cap"]==cap,"Case ID/cap changed")
        for ti,tr in enumerate(TRACERS):
            rec=c["random_permutation_and_levels"][tr]
            seed=20402026+100000*mid+1000*CAPS.index(cap)+100*ti
            perm=np.random.default_rng(seed).permutation(1200)
            subsets=(np.sort(perm[:600]),np.sort(perm[600:]),np.sort(perm[:900]),np.arange(1200,dtype=np.int64))
            require(rec["fixed_numpy_PCG64_seed"]==seed==f["fixed_random_seeds"][tr]
                and arrsha(perm)==rec["permutation_original_1200_index_SHA256"]
                and rec["half_A_and_B_disjoint_and_cover_original_1200"] is True
                and rec["half_A_nested_in_900"] is True
                and not np.intersect1d(subsets[0],subsets[1]).size
                and np.array_equal(np.union1d(subsets[0],subsets[1]),subsets[3])
                and np.all(np.isin(subsets[0],subsets[2])),"Random seed/nested split changed: "+key+tr)
            for lev,idx,n in zip(LEVELS,subsets,(600,600,900,1200)):
                p=rec["per_level"][lev]
                require(p["n"]==n and p["selected_original_1200_index_SHA256"]==arrsha(idx),
                        "Original subset index SHA changed: "+key+tr+lev)
            require(rec["per_level"]["full_1200"]["subcatalogue_four_vector_SHA256"]==
                old["original_selected_array_SHA256"][tr+"_ran"]==
                prior["input_sample_diagnostics"][tr+"_ran"]["selected_array_SHA256"],
                "Original full-1200 random catalogue SHA changed: "+key+tr)
        for tr in TRACERS:
            for role in ("dat","ran"):
                name=tr+"_"+role
                require(old["original_selected_array_SHA256"][name]==
                        prior["input_sample_diagnostics"][name]["selected_array_SHA256"],
                        "Frozen A02 sample SHA changed: "+key+name)
        supports={}
        for lev in LEVELS:
            level=c["levels"][lev]
            supports[lev]=level["scenarios"]["baseline_original"]["RR_supported_cells"]
            require(level["LRG_R_rows"]==level["ELG_R_rows"]==
                    {"half_A_600":600,"half_B_600":600,"three_quarter_A_900":900,"full_1200":1200}[lev],
                    "Unexpected random rows: "+key+lev)
            require(list(level["scenarios"])==list(SCENARIOS if lev!="three_quarter_A_900" else SCENARIOS[:1]),
                    "Unregistered E2 stress scenario: "+key+lev)
            base=level["scenarios"]["baseline_original"]
            for sc,v in level["scenarios"].items():
                scenario_count+=1; cells=v["RR_supported_cells"]
                mask=np.ones((6,24),dtype=np.uint8)
                for idx in v["RR_missing_indices"]:
                    require(len(idx)==2 and 0<=idx[0]<6 and 0<=idx[1]<24,"Bad sparse RR index")
                    mask[idx[0],idx[1]]=0
                require(cells==int(mask.sum()) and len(v["RR_missing_indices"])==144-cells
                    and arrsha(mask)==v["RR_support_mask_sha256"]
                    and v["no_unsupported_RR_zero_fill_or_projection"] is True
                    and v["max_abs_xi_reverse_on_supported_cells"]<1e-8,
                    "RR support/signed mirror failed: "+key+lev+sc)
                full=cells==144
                if not full:sparse+=1
                require((v["original_physical_domain_multipoles"] is not None)==full
                    and (v["full_6x24_xi_SHA256"] is not None)==full
                    and v["status"]==("complete_full_RR_support" if full else "complete_sparse_RR_support_no_full_multipoles"),
                    "Incomplete RR incorrectly projected: "+key+lev+sc)
                if full:
                    require(all(len(v["original_physical_domain_multipoles"][str(ell)]["values_by_fixed_s_bin"])==6
                        for ell in range(4)),"Missing complete projected multipoles")
                for name,rname in MIRROR.items():
                    a=v["forward_pair_terms"][name];b=v["reverse_pair_terms"][rname]
                    require(a["accepted_pairs"]==b["accepted_pairs"]
                        and abs(a["independently_normalized_pair_weight"]-b["independently_normalized_pair_weight"])<1e-8
                        and a["normalization_relative_residual"]<1e-10
                        and b["normalization_relative_residual"]<1e-10,
                        "Four independently normalized LS terms differ: "+key+lev+sc)
                if lev=="full_1200":
                    parent=old["scenarios"][sc]
                    require(cells==144 and v["full_6x24_xi_SHA256"]==parent["forward_xi_grid_sha256"]==
                        f["original_E2_full_1200_xi_SHA256_by_scenario"][sc],
                        "Full-1200 exact original E2 xi fingerprint changed: "+key+sc)
                    for side in ("forward_pair_terms","reverse_pair_terms"):
                        for term in TERMS:
                            a=v[side][term];b=parent[side][term]
                            require(a["weighted_histogram_SHA256"]==b["weighted_histogram_SHA256"]
                                and a["accepted_pairs"]==b["accepted_pairs"],
                                "Full-1200 original E2 weighted pair SHA differs: "+key+sc+side+term)
                    if sc!="baseline_original":
                        require(v["ELG_R_weight_stress"]["modified_selected_weight_SHA256"]==
                            parent["ELG_R_weight_stress"]["modified_selected_ELG_R_weight_SHA256"],
                            "Full-1200 original E2 changed ELG random SHA differs")
                if sc!="baseline_original":
                    require(v["RR_support_mask_sha256"]==base["RR_support_mask_sha256"]
                        and v["ELG_R_weight_stress"]["only_ELG_RANDOM_weight_modified"] is True
                        and v["ELG_R_weight_stress"]["factor_min"]>=.95-1e-12
                        and v["ELG_R_weight_stress"]["factor_max"]<=1.05+1e-12,
                        "Unexpected ELG random-only stress or RR support changed")
                    for side,unaffected in (("forward_pair_terms",("D1D2","R1D2")),
                                            ("reverse_pair_terms",("D1D2","D1R2"))):
                        for term in unaffected:
                            require(v[side][term]==base[side][term],
                                    "Unexpected non-ELG random pair change")
            for amp in ("0.02","-0.02"):
                z=level["DD_pair_level_synthetic_odd_injections"][amp]
                inject_count+=1
                require(z["amplitude"]==float(amp)
                    and z["post_binning_DD_only_modification"] is True
                    and z["original_DD_pair_count_and_normalization_unchanged"] is True
                    and z["not_a_physical_galaxy_wake_injection"] is True
                    and z["positive_RR_supported_cell_count"]==base["RR_supported_cells"],
                    "Unregistered/survey DD injection detected: "+key+lev+amp)
                for tag,field,threshold in (
                    ("forward_analytic","max_abs_analytic_xi_increment_residual_on_supported_cells",1e-11),
                    ("reverse_analytic","max_abs_reverse_analytic_xi_increment_residual_on_supported_cells",1e-11),
                    ("mirror_increment","max_abs_forward_reverse_injected_increment_mirror_residual_on_supported_cells",1e-8),
                    ("mirror_injected_xi","max_abs_forward_reverse_injected_xi_residual_on_supported_cells",1e-8)):
                    v=z[field]
                    require(0<=v<threshold,"Signed DD estimator injection failed: "+key+lev+amp+field)
                    maxres[tag]=max(maxres[tag],v)
                full=base["RR_supported_cells"]==144
                projected=z["full_ell0to3_injected_increment_only_if_144_supported"]
                require((projected is not None)==full and
                    z["unsupported_RR_cells_have_no_injected_xi_or_full_multipole"] is (not full),
                    "Sparse DD injection projected or full RR injection skipped")
                if full:
                    for ell in map(str,range(4)):
                        u=projected[ell]["delta_injected_minus_baseline_by_s_bin"]
                        v=projected[ell]["expected_finite_bin_pair_analytic_increment_by_s_bin"]
                        require(len(u)==len(v)==6 and max(abs(a-b) for a,b in zip(u,v))<1e-11,
                                "Synthetic complete 6-bin multipole recovery failed")
        require(supports==f["RR_supported_baseline_by_level"],
                "Per-case fixed RR support differs from separately archived manifest")
        entry=lines[73+keys.index(key)]
        mline=re.fullmatch(r"A03E3_FIXED_MOCK_CASE (\d{4}/(?:NGC|SGC)) RR_SUPPORT (.*)",entry)
        require(mline is not None and mline.group(1)==key and ast.literal_eval(mline.group(2))==supports,
                "Exact stdout per-case RR support or order changed")
    require(scenario_count==m["number_of_scenario_records"]==180
        and sparse==m["number_of_sparse_scenario_records_without_full_multipoles"]==56
        and inject_count==m["number_of_synthetic_DD_injections"]==144,
        "E3 18-case/180-scenario/144-injection summary drift")
    for tag,v in maxres.items():
        require(v==m["max_synthetic_DD_residuals"][tag],
                "Manifest vs actual DD injection residual maximum changed: "+tag)
    return scenario_count,sparse,inject_count,maxres

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    r,l,m,p,e,a=(q.read_bytes() for q in (REPORT,LOG,MANIFEST,PROTOCOL,E2,A02))
    require(len(r)==1848217 and sha(r)==REPORT_SHA and blob(r)==REPORT_BLOB
        and len(l)==11544 and sha(l)==LOG_SHA and blob(l)==LOG_BLOB
        and blob(m)==MANIFEST_BLOB and sha(p)==PROTOCOL_SHA
        and sha(e)==E2_SHA and sha(a)==A02_SHA,
        "Exact original archive/log/manifest/protocol/parent bytes or SHA changed")
    j=json.loads(r); manifest=json.loads(m);parent=json.loads(e);a02=json.loads(a)
    result=audit(j,manifest,parent,a02,l.decode().splitlines())
    if args.self_test:
        modified=copy.deepcopy(j)
        modified["observed_odd_data_vector_read"]=True
        try:audit(modified,manifest,parent,a02,l.decode().splitlines())
        except ValueError:pass
        else:raise AssertionError("Observed-odd flag tamper was accepted")
        modified=copy.deepcopy(j)
        modified["cases"]["0001/NGC"]["levels"]["half_A_600"]["scenarios"]["baseline_original"]["RR_missing_indices"]=[]
        try:audit(modified,manifest,parent,a02,l.decode().splitlines())
        except ValueError:pass
        else:raise AssertionError("Sparse RR mask tamper was accepted")
    print("A03E3_ORIGINAL_UPLOAD_1848217_BYTE_18_CASE_ARCHIVE_SOURCE_ONLY_OK",flush=True)
    print("SCENARIOS",result[0],"SPARSE_NO_PROJECTION",result[1],"DD_INJECTIONS",result[2],flush=True)
    print("MAX_SYNTHETIC_DD_INCREMENT_RESIDUALS",result[3],flush=True)
    print("OBSERVED_ODD_READ",False,"NO_FITS",True,flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
