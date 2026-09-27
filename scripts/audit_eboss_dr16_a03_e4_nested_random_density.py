#!/usr/bin/env python3
"""A-03E4: frozen MOCK-ONLY 1200 -> 2400 -> 4800 nested matched random-density test.

Rehash the 72 original compressed full mock-galaxy/random files before ANY
new E4 FITS header/row. Use the EXACT original A02 600D/1200R sample and
select only extra RANDOM rows from the eligible original same-source complement.
All nine original IDs and both caps are retained. NO observed data or new URLs.
This reports finite-random sensitivity, NOT a survey window or detection.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection as E3
from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import sample_catalogue
from eboss_dr16_fiducial import WEIGHT_COLUMNS, validated_weight_product

ROOT=A02.ROOT
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e4_nested_1200_2400_4800_mock_random_density_protocol_2026-09-27.json"
E3_REPORT=ROOT/"source_data/eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json"
E3_MANIFEST=ROOT/"source_data/eboss_dr16_a03_e3_mock_galaxy_nested_random_split_synthetic_dd_uploaded_manifest_2026-09-26.json"
E3_STDOUT=ROOT/"source_data/eboss_dr16_a03_e3_local_stdout_2026-09-26.log"
E2_REPORT=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json"
A02_REPORT=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
A02_RUNNER=ROOT/"scripts/audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport.py"
A02_SAMPLER=ROOT/"scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py"
E3_RUNNER=ROOT/"scripts/audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection.py"
PASS="A03E4_NINE_MOCK_NESTED_1200_2400_4800_RANDOM_DENSITY_DESCRIPTIVE_ONLY"
STOP="A03E4_NINE_MOCK_NESTED_RANDOM_DENSITY_INCOMPLETE_STOP"
LEVELS=("full_1200_original","nested_2400","nested_4800")
EXTRA_SEED_ROOT=20502027
TOTAL_SOURCE_BYTES=1860198719

def sha(raw):return hashlib.sha256(raw).hexdigest()
def file_sha(path):return sha(path.read_bytes())
def arr_sha(a):return E3.arr_sha(a)
def blob(path):return E3.E2.blob(path)
def assert_true(ok,msg):
    if not ok:raise ValueError(msg)

def source_gate(p):
    assert_true(
        p["registered_before"]=="Any first E4 FITS header, row or E4 extra random index draw; E3 original 18/18 completed and exact uploaded archive independently audited"
        and p["user_decision"]=="EMPIRICAL_ONLY_DO_NOT_CONTACT_AUTHORS"
        and p["same_fixed_mock_ids"]==list(A02.IDS)
        and p["caps"]==list(A02.CAPS) and p["tracers"]==list(A02.TRACERS)
        and p["roles"]==list(A02.ROLES) and p["fixed_cases"]==18
        and p["random_extension"]["levels_in_order"]==list(LEVELS)
        and p["random_extension"]["supplement_seed_formula"]==
            "20502027+100000*mockID+1000*capIndex+100*tracerIndex; numpy.random.default_rng(seed).choice(len(eligible_complement),size=3600,replace=False)"
        and p["parent_E3_report_sha256"]=="1e52f95be973aa6616337c2427fde60725a0df31ce7d3a818e99d1dd1fb0fb35"
        and p["parent_E3_report_git_blob_sha1"]=="fa32a0f1546307e7bcdfe3a7911e613b83bb4f53"
        and p["parent_E3_uploaded_manifest_git_blob_sha1"]=="3674a0345e22391143fd7b3c10a83a871e0b72cf"
        and p["parent_E3_protocol_sha256"]=="bb4cd60ff84a167a09d91324f8a2fce8948f15b86766ae82f4bed121df384a0d"
        and p["parent_A02_report_sha256"]=="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
        and p["parent_A03E2_report_sha256"]=="6ef86f904bb0cdc407adabe2e125146ce4d0671ea1ac932e13d11fc481cd596a"
        and p["scenario"].startswith("baseline_original only for 1200,2400,4800")
        and p["observed_galaxy_rows_read"] is False and p["observed_random_rows_read"] is False
        and p["observed_odd_data_vector_read"] is False
        and p["new_science_selection_applied"] is False
        and p["author_contact_permitted"] is False
        and p["new_downloads_permitted"] is False
        and p["main_mutation_permitted"] is False
        and p["physical_empirical_pair_window_certified"] is False
        and p["inferential_18D_covariance_computed"] is False
        and p["not_a_detection_or_exclusion"] is True,
        "E4 preregistration/parent/no-observed contract changed")
    assert_true(
        E3_REPORT.stat().st_size==p["parent_E3_report_exact_bytes"]==1848217
        and file_sha(E3_REPORT)==p["parent_E3_report_sha256"]
        and blob(E3_REPORT)==p["parent_E3_report_git_blob_sha1"]
        and file_sha(E3_STDOUT)==p["parent_E3_stdout_sha256"]
        and blob(E3_MANIFEST)==p["parent_E3_uploaded_manifest_git_blob_sha1"]
        and file_sha(E3.PROTOCOL)==p["parent_E3_protocol_sha256"]
        and blob(E3.PROTOCOL)==p["parent_E3_protocol_git_blob_sha1"]
        and file_sha(E2_REPORT)==p["parent_A03E2_report_sha256"]
        and file_sha(A02_REPORT)==p["parent_A02_report_sha256"]
        and blob(A02_RUNNER)==p["parent_A02_implementation_git_blob_sha1"]
        and blob(A02_SAMPLER)==p["parent_A02_sampler_git_blob_sha1"]
        and blob(E3_RUNNER)==p["parent_E3_runner_git_blob_sha1"],
        "E4 frozen original E3/E2/A02 uploaded bytes/protocol/source drift")
    parent=A02.load(E3_REPORT)
    manifest=A02.load(E3_MANIFEST)
    expected=[f"{mid:04d}/{cap}" for mid in A02.IDS for cap in A02.CAPS]
    assert_true(parent["status"]==E3.PASS and parent["completed_cases"]==18
        and parent["failed_cases"]==0 and parent["errors"]==[]
        and list(parent["cases"])==expected
        and manifest["original_uploaded_sha256"]==p["parent_E3_report_sha256"]
        and parent["observed_odd_data_vector_read"] is False
        and parent["all_72_full_gzip_SHA_rehashed_before_any_E3_FITS_row"] is True
        and all(parent["cases"][key]["status"]=="complete"
                and parent["cases"][key]["full_original_e2_all_three_pair_and_xi_SHA_replayed"] is True
                and parent["cases"][key]["no_observed_rows_read"] is True
                for key in expected),
        "Original 18/18 E3 report/source-only acceptance changed")
    a02p,a02,e2,_=E3.exact_gate(A02.load(E3.PROTOCOL))
    return parent,a02p,a02,e2

def supplement_positions(eligible_count,original_positions,mid,cap,tracer):
    if mid not in A02.IDS or cap not in A02.CAPS or tracer not in A02.TRACERS:
        raise ValueError("Unregistered E4 ID, cap or tracer")
    original_positions=np.asarray(original_positions,dtype=np.int64)
    if (original_positions.shape!=(1200,)
        or np.any(original_positions<0) or np.any(original_positions>=eligible_count)
        or not np.array_equal(original_positions,np.unique(original_positions))):
        raise ValueError("Original A02 sorted 1200 eligible positions changed")
    available=np.setdiff1d(np.arange(eligible_count,dtype=np.int64),
                           original_positions,assume_unique=True)
    if len(available)<3600:
        raise ValueError("Insufficient unchanged eligible random complement for fixed 4800")
    seed=EXTRA_SEED_ROOT+100000*mid+1000*A02.CAPS.index(cap)+100*A02.TRACERS.index(tracer)
    extra=available[np.random.default_rng(seed).choice(len(available),size=3600,replace=False)]
    positions={
        "full_1200_original":original_positions,
        "nested_2400":np.sort(np.concatenate((original_positions,extra[:1200]))),
        "nested_4800":np.sort(np.concatenate((original_positions,extra)))}
    if (any(len(np.unique(positions[lev]))!=n for lev,n in zip(LEVELS,(1200,2400,4800)))
        or not np.all(np.isin(positions["full_1200_original"],positions["nested_2400"]))
        or not np.all(np.isin(positions["nested_2400"],positions["nested_4800"]))):
        raise ValueError("Fixed 1200/2400/4800 random supplement is not strictly nested")
    return seed,extra,positions

def extend_random(path,original,meta,mid,cap,tracer,expected_rows,origproto):
    high=origproto["exact_fixed_highz_bin"]
    seed_root=origproto["sampling"]["seed_root"]
    oldseed=seed_root+10000*A02.CAPS.index(cap)+100*A02.TRACERS.index(tracer)+1
    assert_true(oldseed==meta["seed"] and meta["selected_rows"]==1200,
                "A02 original 1200 seed/count differs before E4 extension")
    with fits.open(path,memmap=False) as hdus:
        hdus.verify("exception")
        tabs=[h for h in hdus if isinstance(h,fits.BinTableHDU)]
        if len(tabs)!=1 or int(tabs[0].header["NAXIS2"])!=expected_rows:
            raise ValueError("Source-pinned E4 random BINTABLE/header differs")
        table=tabs[0]
        required={"RA","DEC","Z",*WEIGHT_COLUMNS}
        if tracer=="eBOSS_ELG":required.add("chunk")
        if not required.issubset(table.columns.names):
            raise ValueError("Source-pinned random missing required original fields")
        data=table.data
        z=np.asarray(data["Z"],dtype="f8")
        ra=np.asarray(data["RA"],dtype="f8")
        dec=np.asarray(data["DEC"],dtype="f8")
        if (len(z)!=expected_rows or not np.isfinite(z).all()
            or not np.isfinite(ra).all() or not np.isfinite(dec).all()
            or np.any(ra<0) or np.any(ra>=360)
            or np.any(dec<-90) or np.any(dec>90)):
            raise ValueError("Changed original random source coordinate/value gate")
        candidate=(z>=high[0])&(z<high[1])
        weights,retained=validated_weight_product({
            name:np.asarray(data[name][candidate],dtype="f8") for name in WEIGHT_COLUMNS})
        eligible=np.flatnonzero(candidate)[retained]
        eligible_weight=weights[retained]
        if (len(eligible)!=meta["eligible_highz_rows_after_fixed_weight_gate"]
            or int(candidate.sum())!=meta["candidate_highz_rows_before_weight_validation"]
            or int(candidate.sum())-len(eligible)!=meta["numerical_zero_weight_excluded_rows"]):
            raise ValueError("E4 re-opened full random source eligible count differs from A02")
        original_positions=np.sort(np.random.default_rng(oldseed).choice(
            len(eligible),size=1200,replace=False))
        expected_original=tuple(np.asarray(v,dtype="f8").copy() for v in (
            ra[eligible[original_positions]],dec[eligible[original_positions]],
            z[eligible[original_positions]],eligible_weight[original_positions]))
        if (E3.catalogue_sha(expected_original)!=meta["selected_array_SHA256"]
            or any(not np.array_equal(x,y) for x,y in zip(expected_original,original))):
            raise ValueError("Original A02 1200 random was reselected/altered")
        seed,extra,positions=supplement_positions(
            len(eligible),original_positions,mid,cap,tracer)
        catalogs,summary={},{}
        for level in LEVELS:
            inds=positions[level]
            cat=tuple(np.asarray(v,dtype="f8").copy() for v in (
                ra[eligible[inds]],dec[eligible[inds]],z[eligible[inds]],eligible_weight[inds]))
            if (not all(len(a)==len(inds) and np.isfinite(a).all() for a in cat)
                or np.any(cat[3]<=0) or np.any(cat[2]<.9) or np.any(cat[2]>=1.)):
                raise ValueError("Extended E4 random vector violates original A02 gates")
            catalogs[level]=cat
            summ={"n":len(inds),"selected_eligible_positions_SHA256":arr_sha(inds),
                "selected_original_FITS_row_indices_SHA256":arr_sha(eligible[inds]),
                "four_vector_catalogue_SHA256":E3.catalogue_sha(cat)}
            if tracer=="eBOSS_ELG":
                labels=np.asarray(data["chunk"])[eligible[inds]]
                def label(x):
                    if isinstance(x,(bytes,np.bytes_)):x=x.decode("ascii",errors="strict")
                    v=str(x).strip()
                    if not v or len(v)>128:raise ValueError("Invalid original exact ELG chunk label")
                    return v
                unique,count=np.unique(labels,return_counts=True)
                summ["exact_selected_ELG_chunk_counts"]={
                    label(v):int(n) for v,n in zip(unique,count)}
                if sum(summ["exact_selected_ELG_chunk_counts"].values())!=len(inds):
                    raise ValueError("ELG exact chunk descriptive count loss")
            summary[level]=summ
        if (summary["full_1200_original"]["four_vector_catalogue_SHA256"]!=
            meta["selected_array_SHA256"] or not np.array_equal(
                positions["full_1200_original"],original_positions)):
            raise ValueError("Original 1200R SHA or row order no longer matches")
        return catalogs,{"fixed_supplement_PCG64_seed":seed,
            "eligible_after_original_A02_gate":len(eligible),
            "original_1200_eligible_positions_SHA256":arr_sha(original_positions),
            "supplement_3600_UNSORTED_rng_order_SHA256":arr_sha(extra),
            "supplement_disjoint_from_original_1200":True,
            "nested_1200_within_2400_within_4800":True,"levels":summary}

def comparisons(left,right):
    rec=E3.compare_supported(left,right)
    mask=left["support"]&right["support"]
    if not np.any(mask):
        rec["relative_L1_xi_shift_on_common_support"]=None
    else:
        numerator=float(np.abs(left["xi"][mask]-right["xi"][mask]).sum())
        denominator=float(np.abs(right["xi"][mask]).sum())
        rec["relative_L1_xi_shift_on_common_support"]=(
            float(numerator/denominator) if denominator>0 else None)
    rec["relative_L1_denominator_convention"]="sum abs(xi of right-hand comparison baseline on exactly common positive-RR cells); null if zero"
    return rec

def one_case(mid,cap,cats,meta,paths,prior,e2case,e3case,origproto):
    key=f"{mid:04d}/{cap}"
    for tracer in A02.TRACERS:
        for role in A02.ROLES:
            field=tracer+"_"+role
            assert_true(meta[field]==prior["input_sample_diagnostics"][field],
                        "Original A02 deterministic 600D/1200R source metadata changed: "+key+field)
    bytracer={}
    extras={}
    for tracer in A02.TRACERS:
        path=paths[A02.key(mid,cap,tracer)+"/ran"]
        saved=meta[tracer+"_ran"]
        samples,diagnostic=extend_random(
            path,cats[tracer,"ran"],saved,mid,cap,tracer,
            saved["input_FITS_rows"],origproto)
        assert_true(samples["full_1200_original"] is not None and
            diagnostic["levels"]["full_1200_original"]["four_vector_catalogue_SHA256"]==
            e3case["random_permutation_and_levels"][tracer]["per_level"]["full_1200"]["subcatalogue_four_vector_SHA256"],
            "Parent E3 original 1200 four-vector SHA changed: "+key+tracer)
        bytracer[tracer]=samples
        extras[tracer]=diagnostic
    output={"status":"complete","id":mid,"cap":cap,
        "original_A02_all_four_selected_catalogue_SHA256":{
            k:v["selected_array_SHA256"] for k,v in meta.items()},
        "original_A02_exact_input_metadata_replayed":True,
        "parent_E3_original_1200_random_catalogue_SHA_replayed":True,
        "fixed_random_extensions_by_tracer":extras,
        "levels":{},"pairwise_common_support_density_comparisons":{},
        "no_observed_galaxy_or_random_rows_read":True,
        "no_observed_odd_data_read":True,
        "not_a_physical_window_or_inferential_covariance":True}
    internal={}
    baseline={}
    for level in LEVELS:
        rl=bytracer["eBOSS_LRG"][level]
        re=bytracer["eBOSS_ELG"][level]
        rec,arrays=E3.calculate(cats,rl,re)
        if level=="full_1200_original":
            E3.compare_e2_level(e2case,"baseline_original",rec,
                                {"modified_selected_weight_SHA256":arr_sha(re[3])})
            parent=e3case["levels"]["full_1200"]["scenarios"]["baseline_original"]
            assert_true(rec["full_6x24_xi_SHA256"]==parent["full_6x24_xi_SHA256"]
                and rec["RR_support_mask_sha256"]==parent["RR_support_mask_sha256"],
                "E4 original 1200 xi or positive RR mask differs from E3: "+key)
        else:
            for side in ("forward_pair_terms","reverse_pair_terms"):
                assert_true(rec[side]["D1D2"]==
                    baseline["full_1200_original"][side]["D1D2"],
                    "Changing RANDOM density changed original mock galaxy-galaxy DD pairs: "+key)
        output["levels"][level]={"LRG_R_rows":len(rl[0]),"ELG_R_rows":len(re[0]),
            "baseline_original":rec,"not_a_survey_selection_or_wake_inference":True}
        internal[level]=arrays
        baseline[level]=rec
    for l,r in (("nested_2400","full_1200_original"),
                ("nested_4800","nested_2400"),
                ("nested_4800","full_1200_original")):
        output["pairwise_common_support_density_comparisons"][l+"_vs_"+r]=comparisons(
            internal[l],internal[r])
    return output

def self_test(p):
    parent,a02p,a02,e2=source_gate(p)
    for mid,cap,tracer in ((1,"NGC","eBOSS_LRG"),(125,"SGC","eBOSS_ELG")):
        original=np.sort(np.random.default_rng(93127).choice(8000,size=1200,replace=False))
        seed,extra,levels=supplement_positions(8000,original,mid,cap,tracer)
        seed2,extra2,levels2=supplement_positions(8000,original,mid,cap,tracer)
        assert_true(seed==seed2 and np.array_equal(extra,extra2) and
            all(np.array_equal(levels[x],levels2[x]) for x in LEVELS)
            and len(np.intersect1d(original,extra))==0
            and all(len(levels[k])==n for k,n in zip(LEVELS,(1200,2400,4800))),
            "E4 pure synthetic nested random source-only regression failed")
    for fake,original in ((4000,np.arange(1200)),(8000,np.ones(1200,dtype=np.int64))):
        try:supplement_positions(fake,original,1,"NGC","eBOSS_LRG")
        except ValueError:pass
        else:raise AssertionError("E4 invalid eligible complement or duplicate original accepted")
    q=copy.deepcopy(p);q["observed_odd_data_vector_read"]=True
    try:source_gate(q)
    except ValueError:pass
    else:raise AssertionError("Tampered observed-odd-access prereg flag was accepted")
    assert_true(len(parent["cases"])==len(a02["cases"])==len(e2["cases"])==18,
                "Original source-only 18-case parent gate incomplete")
    print("A03E4_FROZEN_SOURCE_ONLY_1200_2400_4800_SYNTHETIC_SELF_TEST_OK",flush=True)

def run(p):
    parent,a02p,a02,e2=source_gate(p)
    dest=ROOT/p["local_output"]
    phash=file_sha(PROTOCOL)
    if dest.exists():
        state=A02.load(dest)
        if state.get("status")==STOP and state.get("errors"):
            raise ValueError("Existing E4 STOP with errors: preserve exact original bytes and stdout; no rerun")
        assert_true(state.get("protocol_sha256")==phash
            and state.get("parent_E3_archive_sha256")==p["parent_E3_report_sha256"]
            and state.get("parent_A02_archive_sha256")==p["parent_A02_report_sha256"]
            and state.get("observed_odd_data_vector_read") is False
            and state.get("status") in (STOP,PASS)
            and set(state.get("cases",{})).issubset(parent["cases"]),
            "Prior E4 checkpoint differs from frozen protocol/parent")
        for key,c in state["cases"].items():
            raw=json.dumps(c,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
            assert_true(c["status"]=="complete" and
                sha(raw)==state.get("per_case_checkpoint_SHA256",{}).get(key),
                "Prior successful E4 per-case checkpoint SHA changed: "+key)
    else:
        state=None
    pre=A02.preflight(a02p,require_local=True)
    g,r,gman,rman,ref,original,gproto,rproto,rawproto,origproto=pre
    paths,total=A02.resolve_72_sources(a02p,g,r,gman,rman,ref,gproto,rproto,rawproto)
    assert_true(total==TOTAL_SOURCE_BYTES,"Original 72 gzip source byte count drift")
    if state is None:
        state={"status":STOP,"protocol_sha256":phash,
            "parent_E3_archive_sha256":p["parent_E3_report_sha256"],
            "parent_A02_archive_sha256":p["parent_A02_report_sha256"],
            "all_72_full_gzip_SHA_rehashed_before_any_E4_FITS_row":True,
            "total_original_72_compressed_bytes":total,
            "fixed_ids":list(A02.IDS),"caps":list(A02.CAPS),"random_levels":list(LEVELS),
            "cases":{},"per_case_checkpoint_SHA256":{},"errors":[],
            "observed_galaxy_rows_read":False,"observed_random_rows_read":False,
            "observed_odd_data_vector_read":False,"new_science_selection_applied":False,
            "physical_empirical_pair_window_certified":False,
            "inferential_18D_covariance_computed":False,
            "not_a_detection_or_exclusion":True}
        A02.atomic(dest,state)
    if state["status"]==PASS:
        assert_true(len(state["cases"])==18,"Prior complete E4 checkpoint omitted fixed cases")
        print("A03E4_ALREADY_COMPLETE_AFTER_FULL_72_REHASH",flush=True)
        return state
    for mid in A02.IDS:
        for cap in A02.CAPS:
            key=f"{mid:04d}/{cap}"
            if key in state["cases"]:
                print("A03E4_REUSED_IMMUTABLE_SUCCESS",key,flush=True)
                continue
            if mid!=1 and any(state["cases"].get(f"0001/{c}",{}).get("status")!="complete"
                              for c in A02.CAPS):
                raise ValueError("Both original fixed ID0001 cap cases must complete before other IDs")
            cats,metadata={},{}
            try:
                previous=next((x for x in original["cases"] if x["cap"]==cap),None) if mid==1 else None
                for tracer in A02.TRACERS:
                    for role in A02.ROLES:
                        path=paths[A02.key(mid,cap,tracer)+"/"+role]
                        rows=A02.source_header_rows(path,previous,tracer,role)
                        values,meta=sample_catalogue(
                            path,expected_rows=rows,cap=cap,tracer=tracer,
                            role=role,expected_highz=None,p=origproto)
                        cats[tracer,role]=values
                        metadata[tracer+"_"+role]=meta
                current=one_case(mid,cap,cats,metadata,paths,a02["cases"][key],
                                 e2["cases"][key],parent["cases"][key],origproto)
            except Exception as exc:
                state["cases"][key]={"status":"incomplete_stop","id":mid,"cap":cap,
                    "errors":[str(exc)],"source_roles_replayed_before_failure":list(metadata),
                    "no_post_result_reselection":True}
                state["errors"].append(key+": "+str(exc))
                A02.atomic(dest,state)
                raise ValueError("E4 fixed case failed, preserving original checkpoint: "+key+
                    " "+str(exc)) from exc
            state["cases"][key]=current
            state["per_case_checkpoint_SHA256"][key]=sha(
                json.dumps(current,sort_keys=True,separators=(",",":"),allow_nan=False).encode())
            A02.atomic(dest,state)
            print("A03E4_FIXED_MOCK_CASE",key,"RR_SUPPORT",{
                level:current["levels"][level]["baseline_original"]["RR_supported_cells"]
                for level in LEVELS},flush=True)
    state["completed_cases"]=sum(x.get("status")=="complete" for x in state["cases"].values())
    state["failed_cases"]=sum(x.get("status")!="complete" for x in state["cases"].values())
    state["status"]=PASS if state["completed_cases"]==18 and state["failed_cases"]==0 and not state["errors"] else STOP
    A02.atomic(dest,state)
    print("A03E4_NESTED_1200_2400_4800_MOCK_ONLY",state["status"],flush=True)
    print("REPORT",dest,flush=True)
    print("COMPLETED_CASES",state["completed_cases"],flush=True)
    print("OBSERVED_ODD_DATA_READ",False,flush=True)
    return state

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    p=A02.load(PROTOCOL)
    if args.self_test:
        self_test(p);return 0
    try:
        out=run(p)
        return 0 if out["status"]==PASS else 2
    except Exception as exc:
        path=ROOT/p["local_output"]
        if path.exists():
            try:
                old=A02.load(path)
                if old.get("status")==STOP and old.get("errors"):
                    print("A03E4_PRESERVED_ORIGINAL_STOP_NO_RESERIALIZATION",path,flush=True)
                    print("A03E4_ORIGINAL_STOP_ERRORS",chr(10).join(old["errors"]),flush=True)
                    return 2
            except Exception:pass
        print("A03E4_FAILED_CLOSED",str(exc),flush=True)
        return 2

if __name__=="__main__":
    raise SystemExit(main())
