#!/usr/bin/env python3
"""SOURCE-PINNED SGC-only recovery after original full-D mock0001 SGC parity STOP.

Original v1 NGC is reused byte-for-byte from its immutable SHA-pinned JSON.
Original v1 SGC failed independent sparse reverse counting at R1D2: its
mismatch magnitude was NOT checkpointed and its cause is not presumed.
For the production-scale pair estimator we count each forward unordered
tracer pair exactly once, then construct its reversed orientation by the
mathematical signed-mu mirror, preserving the original separate LS
normalizations. This does NOT claim an independent full-size reverse count.
Small original 600D/1200R forward/reverse independent dense SHA replay and
sparse-vs-dense closure remain mandatory. SGC forward pair histograms are
atomically checkpointed individually so WSL termination does not lose work.
The default --run performs ONLY SGC original 4800R; 48000R is opt-in.
No observed galaxy, observed random, observed odd, source download or new cut.
"""
from __future__ import annotations

import argparse
import fcntl
import gc
import hashlib
import json
import math
import tempfile
from pathlib import Path

import numpy as np

import run_eboss_dr16_full_eligible_mock0001_sparse as V1

ROOT=V1.ROOT
PROTOCOL=ROOT/"source_data/eboss_dr16_full_eligible_mock0001_sgc_resume_v2_after_failure_protocol_2026-09-27.json"
PROTOCOL_BLOB="c2ff44c037979c4695523dcc4027072928bf57d4"
V21_CONTRACT=ROOT/"source_data/eboss_dr16_mock0001_sgc_48000_chunked_v21_engineering_protocol_2026-09-27.json"
V21_BLOB="20481d7daffcb49fad0f2b554c67925b29846981"
FIRST_SLICE_ROWS=4000
V1_BLOB="793a69301afc41c5254a78b2d56fe40c0cba4876"
CAP="SGC"
KEY="0001/SGC"
NGC="0001/NGC"
SUCCESS="EBOSS_FULL_ELIGIBLE_MOCK0001_SGC_RECOVERY_V2_ENGINEERING_ONLY"
STAGED="EBOSS_FULL_ELIGIBLE_MOCK0001_SGC_4800_ONLY_STAGED_STOP"
STOP="EBOSS_FULL_ELIGIBLE_MOCK0001_SGC_RECOVERY_INCOMPLETE_STOP"
BUDGET=60_000_000


def need(ok,msg):
    if not ok:raise ValueError(msg)


def jsha(obj):
    return V1.sha(json.dumps(obj,sort_keys=True,separators=(",",":"),
                             allow_nan=False).encode())


def verify_array(rec, *, pair=True):
    vals=rec["histogram_6x24"] if pair else rec["xi_6x24"]
    a=np.ascontiguousarray(np.asarray(vals,dtype="f8"))
    need(a.shape==(6,24) and np.isfinite(a).all() and
         (not pair or np.all(a>=0)) and
         V1.sha(a.tobytes())==rec["histogram_SHA256" if pair else "xi_SHA256"],
         "Frozen sparse histogram/xi SHA, shape, finite or positivity mismatch")
    if pair:
        need(isinstance(rec["accepted_pairs"],int) and rec["accepted_pairs"]>=0
             and rec["pair_normalization"]>0 and
             math.isfinite(rec["pair_normalization"]),
             "Frozen weighted pair accepted count or normalization invalid")
    return a


def guards():
    raw=PROTOCOL.read_bytes()
    need(V1.blob(raw)==PROTOCOL_BLOB, "V2 postmortem protocol Git blob changed")
    p=json.loads(raw)
    v21raw=V21_CONTRACT.read_bytes()
    c=json.loads(v21raw)
    need(V1.blob(v21raw)==V21_BLOB and
         c["source_branch"]==p["source_branch"] and
         c["original_v2_runner_git_blob_sha1"]=="808b253c4998335b67a4117bac28facc170892e1" and
         c["original_v2_4800_stage_uploaded_sha256"]=="d3f298a9550981c0ece7d6c9e7aee0d3fb5cebaf0d80cdd39805cda3bb08e238" and
         c["scope"]["cap"]==CAP and c["scope"]["level"]=="nested_48000" and
         c["scope"]["first_catalogue_chunk_rows"]==FIRST_SLICE_ROWS and
         c["scope"]["terms"]==["D1R2","R1D2","R1R2"] and
         c["scope"]["unchanged_4800_checkpoint"] is True and
         c["scope"]["observed_odd_data_read"] is False and
         c["scope"]["new_download"] is False and
         c["scope"]["new_science_selection"] is False and
         c["scope"]["main_change"] is False,
         "Pre-48000 post-SGC4800 source-pinned engineering slice protocol changed")
    need(V1.blob(Path(V1.__file__).read_bytes())==V1_BLOB and
         p["original_v1_runner_blob_sha1"]==V1_BLOB and
         p["source_branch"]=="audit/eboss-elg-bit8-ra-orientation-20260925" and
         p["repair_scope"]["cap"]==CAP and
         p["repair_scope"]["old_NGC"]=="READ_ONLY_REUSE_V1_NGC_NO_RECOMPUTATION" and
         p["repair_scope"]["original_v1_report"]=="READ_ONLY_BYTE_FOR_BYTE_PRESERVED" and
         p["algorithm"].startswith("Compute each SGC forward LS term only once") and
         p["large_pair_chunk"]==32 and p["checkpoint_every_pair"] is True and
         p["parent_v1_SGC_error"]==V1_STOP_REASON,
         "V1 source or V2 source-pinned SGC-only repair contract changed")
    for flag in ("observed_galaxy_rows_read","observed_random_rows_read",
                 "observed_odd_data_vector_read","new_fits_download",
                 "new_mock_ids","author_contact","change_main","new_physics_cut",
                 "physical_A03_window_certified","eBOSS_A04_inference_performed"):
        need(p[flag] is False,"Forbidden V2 scope: "+flag)
    # Source-only parental JSON checks and exact Git code provenance. No FITS.
    orig,a02,e4=V1.source_only_gate()
    need(orig["fixed_id"]==[1] and
         p["original_original_E4_report_sha256"]==V1.OLD_E4_SHA and
         p["original_72_total_source_bytes"]==V1.TOTAL_SOURCE_BYTES and
         p["parent_v1_state"]==V1.STOP,
         "V1 source gate or frozen original v1 failure scope changed")
    return p,a02,e4


V1_STOP_REASON="Full-D sparse forward reverse four-term parity failed: R1D2"


def immutable_original(p):
    original=ROOT/"eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_report.json"
    raw=original.read_bytes()
    need(len(raw)==p["parent_v1_full_report_exact_bytes"] and
         V1.sha(raw)==p["parent_v1_full_report_sha256"],
         "Original v1 NGC-success/SGC-STOP report bytes changed. Do not delete/overwrite.")
    j=json.loads(raw)
    need(j["status"]==p["parent_v1_state"] and j["errors"]==[
        KEY+": "+V1_STOP_REASON] and
        j["all_72_full_gzip_rehashed_before_any_FITS"] is True and
        j["source_72_total_compressed_bytes"]==V1.TOTAL_SOURCE_BYTES and
        j["observed_odd_data_vector_read"] is False and
        j["observed_galaxy_rows_read"] is False and
        j["observed_random_rows_read"] is False and
        j["cases"][KEY]["status"]=="incomplete_stop" and
        j["cases"][KEY]["error"]==V1_STOP_REASON and
        jsha(j["cases"][NGC])==j["case_canonical_SHA256"][NGC]==
        p["parent_v1_NGC_canonical_case_sha256"],
        "Original v1 report NGC SHA/scope or exact original SGC error differs")
    for stage in V1.LEVELS:
        d=j["cases"][NGC]["levels"][stage]
        need(d["RR_supported_cells"]==144 and
             d["no_observed_galaxies_or_odd_read"] is True,
             "Original v1 NGC RR support/observed guard changed")
        verify_array(d,pair=False)
        for side in ("forward","reverse"):
            for term in V1.TERMS:
                verify_array(d[side][term])
    return j


def fwd_samples(cats,level):
    dl=cats[("eBOSS_LRG","dat")]["full"]
    de=cats[("eBOSS_ELG","dat")]["full"]
    rl=cats[("eBOSS_LRG","ran")][level]
    re=cats[("eBOSS_ELG","ran")][level]
    return {"D1D2":(dl,de),"D1R2":(dl,re),
            "R1D2":(rl,de),"R1R2":(rl,re)}


def source_signature(info):
    return {
        label:(record["full_eligible_four_vector_SHA256"]
               if label.endswith("_dat")
               else record["nested_48000_four_vector_SHA256"])
        for label,record in sorted(info.items())
    }


def validate_progress(progress,*,synthetic_slice_rows=None):
    # Production accepts only frozen 4159/48000 first- and 15860/48000
    # second-catalogue populations. Source-only synthetic exercises 37x43.
    slice_rows=(FIRST_SLICE_ROWS if synthetic_slice_rows is None
                else synthetic_slice_rows)
    allowed_first=((4159,48000) if synthetic_slice_rows is None else (37,))
    allowed_second=((15860,48000) if synthetic_slice_rows is None else (43,))
    sources=progress.get("sources")
    empty_uninitialized=(sources is None and not progress.get("levels")
                         and not progress.get("assembled_levels"))
    pinned=(isinstance(sources,dict) and len(sources)==4
            and all(isinstance(v,str) and len(v)==64 for v in sources.values()))
    need(progress.get("cap")==CAP and (empty_uninitialized or pinned) and
         set(progress["levels"]).issubset(V1.LEVELS),
         "Existing v2 progress lacks valid original SGC source or level guards")
    for stage,termmap in progress["levels"].items():
        need(set(termmap).issubset(V1.TERMS),
             "Existing v2 progress contains forbidden forward term")
        for term,rec in termmap.items():
            verify_array(rec)
            need(rec.get("forward_counted_not_mirrored") is True and
                 rec["candidate_neighbour_pairs"] >= rec["accepted_pairs"],
                 "Stored checkpoint is not genuinely forward-counted")
    partial=progress.get("partial_48000",{})
    need(isinstance(partial,dict) and set(partial).issubset({"D1R2","R1D2","R1R2"}),
         "Forbidden in-progress 48000R first-catalogue term")
    for term,parts in partial.items():
        need(isinstance(parts,dict),"48000R source-pinned first-row partials malformed")
        for startkey,rec in parts.items():
            need(startkey.isdecimal() and str(int(startkey))==startkey and
                 int(startkey)%slice_rows==0 and
                 int(startkey)>=0 and int(startkey)<48000 and
                 rec["first_row_start"]==int(startkey) and
                 int(startkey)<rec["first_row_stop"]<=48000 and
                 rec["first_row_stop"]-int(startkey)<=slice_rows and
                 rec["first_total_size"] in allowed_first and
                 rec["first_row_stop"]<=rec["first_total_size"] and
                 rec["second_total_size"] in allowed_second and
                 rec.get("forward_counted_not_mirrored") is True and
                 rec["candidate_neighbour_pairs"]>=rec["accepted_pairs"],
                 "Untrusted 48000R per-slice row bounds/counts detected")
            verify_array(rec)
    assembled=progress.get("assembled_levels",{})
    need(set(assembled).issubset(progress["levels"]),
         "Assembled SGC level not backed by stored forward pair checkpoints")
    for stage,rec in assembled.items():
        verify_array(rec,pair=False)
        need(rec["RR_supported_cells"]==144 and
             rec["reverse_source"]==
             "EXACT_GEOMETRIC_SIGNED_MU_MIRROR_NOT_INDEPENDENT_RECOUNT" and
             all(term in progress["levels"][stage] for term in V1.TERMS),
             "Assembled SGC level missing full RR support, mirror provenance or four terms")
        reconstructed=close_mirrored_level(progress["levels"][stage])
        need(reconstructed["xi_SHA256"]==rec["xi_SHA256"],
             "Tampered assembled SGC xi differs from independently rebuilt stored pairs")
    return progress


def init_state(dest,p,old):
    if dest.exists():
        s=json.loads(dest.read_bytes())
        need(s["parent_v1_full_report_sha256"]==
             p["parent_v1_full_report_sha256"] and
             s["parent_v1_NGC_canonical_case_sha256"]==
             p["parent_v1_NGC_canonical_case_sha256"] and
             s["observed_odd_data_vector_read"] is False and
             jsha(s["cases"][NGC])==p["parent_v1_NGC_canonical_case_sha256"] and
             s["status"] in (SUCCESS,STAGED,STOP),
             "Existing v2 checkpoint scope or original NGC parent changed")
        validate_progress(s["sgc_progress"])
        return s
    s={
        "status":STOP,
        "source_protocol_git_blob":PROTOCOL_BLOB,
        "parent_v1_full_report_sha256":p["parent_v1_full_report_sha256"],
        "parent_v1_NGC_canonical_case_sha256":p["parent_v1_NGC_canonical_case_sha256"],
        "parent_original_v1_SGC_failure":V1_STOP_REASON,
        "original_v1_report_unchanged":True,
        "new_SGC_large_pair_reverse_method":"GEOMETRIC_MIRROR_FROM_COUNTED_FORWARD_NOT_INDEPENDENT",
        "original_600D_1200R_independent_both_orientations_required":True,
        "mock_id":1,"caps":["NGC","SGC"],
        "cases":{NGC:old["cases"][NGC]},
        "sgc_progress":{"cap":CAP,"sources":None,"levels":{}},
        "error_history":[],
        "errors":[],
        "observed_galaxy_rows_read":False,
        "observed_random_rows_read":False,
        "observed_odd_data_vector_read":False,
        "physical_window_certified":False,
        "inferential_covariance_computed":False,
        "wake_detection_or_exclusion_computed":False,
        "new_source_download":False
    }
    V1.A02.atomic(dest,s)
    return s



def chunked_48000_forward_pair(s,dest,term,first,second,distance,
                               *,first_slice_rows=FIRST_SLICE_ROWS):
    """Pair-preserving bounded SGC 48k counting with per-slice SHA checkpoints.

    Disjoint consecutive first-catalogue slices each see the COMPLETE second
    catalogue. Sum in frozen increasing first-index order, and normalize once
    by the original full input-catalogue weight sums. This changes only the
    floating addition tree, not the pair set, selection, or physics.
    """
    need(term in ("D1R2","R1D2","R1R2") and
         first_slice_rows>=1 and first_slice_rows<=FIRST_SLICE_ROWS and
         len(first[0])>0 and len(second[0])>0,
         "Illegal 48000R stage, first-slice size or empty source")
    partial=s["sgc_progress"].setdefault("partial_48000",{})
    saved=partial.setdefault(term,{})
    accumulated=np.zeros((6,24),dtype="f8")
    accepted=candidates=0
    n1,n2=len(first[0]),len(second[0])
    sum1=float(np.sum(first[3],dtype="f8"))
    sum2=float(np.sum(second[3],dtype="f8"))
    full_norm=sum1*sum2
    need(math.isfinite(full_norm) and full_norm>0,
         "SGC 48k original whole-catalogue weighted pair normalization invalid")
    for start in range(0,n1,first_slice_rows):
        stop=min(start+first_slice_rows,n1)
        key=str(start)
        if key in saved:
            rec=saved[key]
            verify_array(rec)
            need(rec["first_row_start"]==start and
                 rec["first_row_stop"]==stop and
                 rec["first_total_size"]==n1 and rec["second_total_size"]==n2 and
                 np.isclose(rec["pair_normalization"],
                            float(np.sum(first[3][start:stop],dtype="f8"))*sum2,
                            rtol=1e-13,atol=0),
                 "Saved 48000R first-row slice no longer matches original source weights")
            h=np.asarray(rec["histogram_6x24"],dtype="f8")
            print("V21_REUSE_SHA_VERIFIED_FIRST_SLICE",term,start,stop,flush=True)
        else:
            cat1=tuple(np.asarray(col[start:stop],dtype="f8") for col in first)
            h,meta=V1.kdtree_rr(
                cat1,second,V1.S_EDGES,V1.MU_EDGES,distance,.05,
                chunk=32,max_neighbour_pairs=BUDGET)
            rec={"accepted_pairs":int(meta["accepted_weighted_pair_count"]),
                 "candidate_neighbour_pairs":int(meta["candidate_neighbour_pairs"]),
                 "pair_normalization":float(meta["pair_normalization"]),
                 "histogram_6x24":h.tolist(),
                 "histogram_SHA256":V1.sha(np.ascontiguousarray(h).tobytes()),
                 "forward_counted_not_mirrored":True,
                 "first_row_start":start,"first_row_stop":stop,
                 "first_total_size":n1,"second_total_size":n2}
            verify_array(rec)
            need(rec["pair_normalization"]>0 and
                 np.isclose(rec["pair_normalization"],
                            float(np.sum(first[3][start:stop],dtype="f8"))*sum2,
                            rtol=1e-13,atol=0),
                 "Saved first-row slice weighted normalization invalid")
            saved[key]=rec
            V1.A02.atomic(dest,s)
            print("V21_ATOMIC_FIRST_SLICE_SAVED",term,start,stop,
                  "accepted",rec["accepted_pairs"],
                  "candidates",rec["candidate_neighbour_pairs"],flush=True)
        accumulated+=h
        accepted+=rec["accepted_pairs"]
        candidates+=rec["candidate_neighbour_pairs"]
        need(candidates<=BUDGET and np.isfinite(accumulated).all() and
             np.all(accumulated>=0),"Cumulative SGC 48k candidate budget or histogram invalid")
        del h
        gc.collect()
    return accumulated,{
        "accepted_weighted_pair_count":accepted,
        "candidate_neighbour_pairs":candidates,
        "pair_normalization":full_norm,
        "source_first_slices":len(saved),
        "first_slice_rows":first_slice_rows,
        "sum_first_weights":sum1,
        "sum_second_weights":sum2
    }


def checkpointed_forward(s,dest,cats,level,distance,old_e4):
    progress=s["sgc_progress"]
    stages=progress["levels"]
    if level not in stages:stages[level]={}
    store=stages[level]
    pairs=fwd_samples(cats,level)
    for term in V1.TERMS:
        if term in store:
            verify_array(store[term])
            print("V2_REUSE_VALIDATED_SGC_FORWARD_PAIR",level,term,flush=True)
            continue
        if level=="nested_48000" and term=="D1D2":
            before=stages["original_4800"]["D1D2"]
            store[term]=dict(before,derived_from_identical_galaxy_pair_stage="original_4800")
            V1.A02.atomic(dest,s)
            print("V2_REUSE_UNCHANGED_FULL_GALAXY_DD",level,flush=True)
            continue
        first,second=pairs[term]
        if level=="nested_48000":
            h,meta=chunked_48000_forward_pair(
                s,dest,term,first,second,distance)
        else:
            h,meta=V1.kdtree_rr(first,second,V1.S_EDGES,V1.MU_EDGES,
                                distance,.05,chunk=32,
                                max_neighbour_pairs=BUDGET)
        need(h.shape==(6,24) and np.isfinite(h).all() and np.all(h>=0) and
             meta["pair_normalization"]>0,
             "SGC sparse forward weighted pair histogram invalid")
        rec={
            "accepted_pairs":int(meta["accepted_weighted_pair_count"]),
            "candidate_neighbour_pairs":int(meta["candidate_neighbour_pairs"]),
            "pair_normalization":float(meta["pair_normalization"]),
            "histogram_6x24":h.tolist(),
            "histogram_SHA256":V1.sha(np.ascontiguousarray(h).tobytes()),
            "forward_counted_not_mirrored":True,
            "source_pairs":[len(first[0]),len(second[0])],
        }
        if level=="nested_48000":
            rec["first_catalogue_slice_rows"]=FIRST_SLICE_ROWS
            rec["slice_count"]=meta["source_first_slices"]
        verify_array(rec)
        if level=="original_4800" and term=="R1R2":
            old=old_e4["levels"]["nested_4800"]["baseline_original"]["forward_pair_terms"]["R1R2"]
            need(rec["accepted_pairs"]==old["accepted_pairs"] and
                 np.isclose(float(h.sum()),old["total_weighted_pairs_in_fixed_s_mu_bins"],
                            atol=2e-8,rtol=5e-10),
                 "Original SGC E4 exact same 4800R RR pair population differs")
        store[term]=rec
        V1.A02.atomic(dest,s)
        print("V2_SGC_FORWARD_PAIR_SAVED",level,term,
              "accepted",rec["accepted_pairs"],
              "candidates",rec["candidate_neighbour_pairs"],flush=True)
        del h
        gc.collect()
    return store


def close_mirrored_level(store):
    need(set(store)==set(V1.TERMS),
         "SGC level must contain all four forward paired terms")
    fwd={term:verify_array(store[term]) for term in V1.TERMS}
    norms={term:store[term]["pair_normalization"] for term in V1.TERMS}
    # Reversing the physical tracer order exchanges D1R2 <-> R1D2, flips
    # signed mu, and preserves the exact same pair weights/normalization.
    reversed_h={}
    reverse_records={}
    reverse_norms={}
    for fterm,rterm in V1.MAPPING.items():
        arr=np.ascontiguousarray(fwd[fterm][:,::-1])
        reversed_h[rterm]=arr
        reverse_norms[rterm]=norms[fterm]
        reverse_records[rterm]={
            "accepted_pairs":store[fterm]["accepted_pairs"],
            "candidate_neighbour_pairs":store[fterm]["candidate_neighbour_pairs"],
            "pair_normalization":norms[fterm],
            "histogram_6x24":arr.tolist(),
            "histogram_SHA256":V1.sha(arr.tobytes()),
            "geometric_mirror_of_forward":fterm,
            "independently_recounted":False,
        }
    xi,support=V1.cross_landy_szalay(fwd,norms)
    rev,rsupport=V1.cross_landy_szalay(reversed_h,reverse_norms)
    need(int(support.sum())==144 and
         np.array_equal(support,rsupport[:,::-1]) and
         np.isfinite(xi).all() and
         np.max(np.abs(xi-rev[:,::-1]))<5e-10,
         "Derived physical tracer mirror failed original source geometry/LS algebra")
    out={"forward":store,"reverse":reverse_records,
         "reverse_source":"EXACT_GEOMETRIC_SIGNED_MU_MIRROR_NOT_INDEPENDENT_RECOUNT",
         "xi_6x24":xi.tolist(),
         "xi_SHA256":V1.sha(np.ascontiguousarray(xi).tobytes()),
         "RR_supported_cells":144,
         "max_abs_reverse_xi":float(np.max(np.abs(xi-rev[:,::-1]))),
         "odd_and_even_multipoles":V1.project(xi,support,V1.MU_EDGES,(0,1,2,3)),
         "no_observed_galaxies_or_odd_read":True}
    verify_array(out,pair=False)
    return out


def synthetic_test():
    V1.synthetic_test()
    rng=np.random.default_rng(27112026)
    arrays={term:rng.uniform(.05,1.5,size=(6,24)) for term in V1.TERMS}
    norm={"D1D2":11.,"D1R2":12.,"R1D2":13.,"R1R2":14.}
    reversed_h={r:arrays[f][:,::-1].copy() for f,r in V1.MAPPING.items()}
    reversed_n={r:norm[f] for f,r in V1.MAPPING.items()}
    x,s=V1.cross_landy_szalay(arrays,norm)
    y,rs=V1.cross_landy_szalay(reversed_h,reversed_n)
    need(s.all() and rs.all() and
         np.allclose(x,y[:,::-1],atol=2e-14,rtol=0),
         "Synthetic distinct DR/RD normalizations not preserved on mirror")
    positive={}
    for term in V1.TERMS:
        h=np.ascontiguousarray(arrays[term],dtype="f8")
        positive[term]={"histogram_6x24":h.tolist(),
                        "histogram_SHA256":V1.sha(h.tobytes()),
                        "accepted_pairs":100,"candidate_neighbour_pairs":101,
                        "pair_normalization":norm[term],
                        "forward_counted_not_mirrored":True}
    assembled=close_mirrored_level(positive)
    need(assembled["RR_supported_cells"]==144 and
         assembled["reverse"]["R1D2"]["geometric_mirror_of_forward"]=="D1R2",
         "Synthetic full 4-term checkpointed physical reverse construction failed")
    progress={"cap":CAP,"sources":{key:"1"*64 for key in
             ("eBOSS_LRG_dat","eBOSS_LRG_ran","eBOSS_ELG_dat","eBOSS_ELG_ran")},
             "levels":{"original_4800":positive},
             "assembled_levels":{"original_4800":assembled}}
    validate_progress(progress)
    # This is an actual run of the 48k streaming algorithm on SYNTHETIC
    # catalogues only, with 11-row first-catalogue slices and atomic JSON.
    tiny1=(rng.uniform(180,182,37),rng.uniform(5,6,37),
           rng.uniform(.91,.99,37),rng.uniform(.6,1.7,37))
    tiny2=(rng.uniform(180,182,43),rng.uniform(5,6,43),
           rng.uniform(.91,.99,43),rng.uniform(.7,1.4,43))
    dist=lambda z: np.asarray(z,dtype="f8")*2800.
    direct,dmeta=V1.kdtree_rr(
        tiny1,tiny2,V1.S_EDGES,V1.MU_EDGES,dist,.05,
        chunk=11,max_neighbour_pairs=BUDGET)
    with tempfile.TemporaryDirectory() as folder:
        saved_state={"sgc_progress":{"cap":CAP,
                     "sources":progress["sources"],"levels":{},
                     "partial_48000":{}}}
        out=Path(folder)/"synthetic_atomic_48k_slices.json"
        h,meta=chunked_48000_forward_pair(
            saved_state,out,"R1R2",tiny1,tiny2,dist,first_slice_rows=11)
        need(meta["source_first_slices"]==4 and
             meta["accepted_weighted_pair_count"]==dmeta["accepted_weighted_pair_count"] and
             meta["candidate_neighbour_pairs"]==dmeta["candidate_neighbour_pairs"] and
             np.isclose(meta["pair_normalization"],dmeta["pair_normalization"],
                        rtol=1e-14,atol=0) and
             np.allclose(h,direct,rtol=1e-12,atol=1e-10),
             "Synthetic 4-slice SGC weighted sparse full-vs-sliced exact pair closure failed")
        stored=json.loads(out.read_bytes())
        validate_progress(stored["sgc_progress"],synthetic_slice_rows=11)
        # Replay from persistent storage, not in-memory mutated state.
        h2,m2=chunked_48000_forward_pair(
            stored,out,"R1R2",tiny1,tiny2,dist,first_slice_rows=11)
        need(np.array_equal(h,h2) and
             m2["accepted_weighted_pair_count"]==meta["accepted_weighted_pair_count"],
             "Stored SGC 48k per-slice resume changes completed histograms")
        stored["sgc_progress"]["partial_48000"]["R1R2"]["0"]["histogram_SHA256"]="0"*64
        try:validate_progress(stored["sgc_progress"],synthetic_slice_rows=11)
        except ValueError:pass
        else:raise AssertionError("Tampered 48k partial slice was accepted")
    print("EBOSS_SGC_V21_FOUR_SYNTHETIC_DISJOINT_PAIR_SLICES_AND_RESUME_OK",
          "NO_FITS NO_OBSERVED",flush=True)
    bad={"histogram_6x24":np.ones((6,24)).tolist(),
         "histogram_SHA256":"0"*64,"accepted_pairs":7,
         "pair_normalization":1.}
    try:verify_array(bad)
    except ValueError:pass
    else:raise AssertionError("Tampered SGC intermediate histogram was accepted")
    print("EBOSS_FULL_MOCK_SGC_V2_SOURCE_ONLY_SYNTHETIC_MIRROR_AND_CHECKPOINT_OK",
          "NO_FITS NO_OBSERVED ODD_SEALED",flush=True)


def run(p,a02,e4,*,include_48000):
    old=immutable_original(p)
    v2=ROOT/p["repair_scope"]["new_v2_report"]
    v2.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive process lock is OS-released even if WSL kills the Python PID.
    with (v2.parent/(v2.name+".lock")).open("a+") as lock:
        try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("Another SGC V2 recovery is running; refuse concurrent writes") from exc
        state=init_state(v2,p,old)
        if state["status"]==SUCCESS:
            print("V2_COMPLETE_REUSING_SAVED_SGC_AND_NGC_NO_RECOMPUTE",v2,flush=True)
            return 0
        # Previous source references are fully rehashed before any new FITS
        # even on a resumed checkpoint after sudden WSL process termination.
        original_protocol=V1.A02.load(V1.A02.PROTOCOL)
        g,r,gman,rman,ref,original,gproto,rproto,rawproto,p0=V1.A02.preflight(
            original_protocol,require_local=True)
        paths,total=V1.A02.resolve_72_sources(
            original_protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
        need(len(paths)==72 and total==V1.TOTAL_SOURCE_BYTES,
             "All 72 original full compressed sources must rehash before SGC FITS")
        state["all_72_original_gzip_rehashed_before_any_SGC_FITS"]=True
        V1.A02.atomic(v2,state)
        cats,info={},{}
        for tracer in V1.TRACERS:
            for role in V1.ROLES:
                label=tracer+"_"+role
                path=paths[V1.A02.key(1,CAP,tracer)+"/"+role]
                cats[(tracer,role)],info[label]=V1.load_original_full(
                    path,CAP,tracer,role,a02["cases"][KEY]["input_sample_diagnostics"][label],
                    e4["cases"][KEY],p0)
                print("V2_SGC_FULL_SOURCE_VERIFIED",label,
                      info[label]["full_eligible_rows"],flush=True)
        fingerprint=source_signature(info)
        progress=state["sgc_progress"]
        if progress["sources"] is None:
            progress["sources"]=fingerprint
        need(progress["sources"]==fingerprint,
             "Changed SGC full source/sample SHA since stored partial V2 pair checkpoint")
        V1.A02.atomic(v2,state)
        dist=V1.cached_exact_distance(cats)
        benchmark=V1.validate_sparse_against_original_dense(
            cats,a02["cases"][KEY],dist)
        need(len(benchmark)==8 and
             all(value["sparse_dense_max_abs_weighted_hist_error"]<2e-9
                 for value in benchmark.values()),
             "SGC original eight-pair independent dense/sparse SHA replay failed")
        progress["original_eight_dense_sparse_independent_closure"]=benchmark
        V1.A02.atomic(v2,state)
        levels=("original_4800","nested_48000") if include_48000 else ("original_4800",)
        try:
            for stage in levels:
                store=checkpointed_forward(
                    state,v2,cats,stage,dist,e4["cases"][KEY])
                assembled=close_mirrored_level(store)
                if stage=="original_4800":
                    assembled["original_E4_same_RR_4800_count_and_sum_replayed"]=True
                progress.setdefault("assembled_levels",{})[stage]=assembled
                V1.A02.atomic(v2,state)
                print("V2_SGC_LEVEL_VALIDATED",stage,"RR_SUPPORT",assembled["RR_supported_cells"],
                      "REVERSE_GEOMETRIC_NOT_INDEPENDENT",flush=True)
            if not include_48000:
                state["status"]=STAGED
                state["errors"]=[]
                V1.A02.atomic(v2,state)
                print("V2_SGC_4800_STAGED_STOP_SAVED",v2,flush=True)
                return 0
            c={"mock_id":1,"cap":CAP,
               "original_four_source_SHA_replayed":True,
               "all_original_600D_1200R_eight_weighted_hist_SHAs_replayed":True,
               "sparse_dense_eight_pair_closure":benchmark,
               "source_full_eligible_diagnostics":info,
               "levels":{stage:progress["assembled_levels"][stage] for stage in V1.LEVELS},
               "large_pair_reverse_validity":"DERIVED_EXACT_SIGNED_MU_MIRROR_NOT_INDEPENDENT",
               "original_v1_SGC_reversal_failure_unresolved":V1_STOP_REASON}
            a=np.asarray(c["levels"]["nested_48000"]["xi_6x24"])
            b=np.asarray(c["levels"]["original_4800"]["xi_6x24"])
            denom=float(np.abs(a).sum())
            c["full_galaxy_4800_to_48000_relative_xi_L1"]=(
                float(np.abs(a-b).sum()/denom) if denom else None)
            c["full_galaxy_4800_to_48000_odd_multipole_differences"]={
                str(ell):[
                    float(u-v) for u,v in zip(
                        c["levels"]["nested_48000"]["odd_and_even_multipoles"][str(ell)]["values_by_fixed_s_bin"],
                        c["levels"]["original_4800"]["odd_and_even_multipoles"][str(ell)]["values_by_fixed_s_bin"])]
                for ell in (1,3)}
            c["scope"]="ONE_FULL_GALAXY_MOCK_ONLY_REPAIRED_ENGINEERING_NOT_PHYSICAL_SENSITIVITY"
            state["cases"][KEY]=c
            state["case_canonical_SHA256"]={NGC:jsha(state["cases"][NGC]),KEY:jsha(c)}
            state["status"]=SUCCESS
            if state["errors"]:
                state["error_history"].extend(state["errors"])
            state["errors"]=[]
            state["completed_cases"]=2
            V1.A02.atomic(v2,state)
            print("V2_SGC_COMPLETE_FROZEN_NGC_REUSED",
                  "SGC_XI_L1_4800_TO_48000",
                  c["full_galaxy_4800_to_48000_relative_xi_L1"],
                  "V1_INDEPENDENT_SGC_FAILURE_REMAINS_IN_IMMUTABLE_PARENT",
                  flush=True)
            print("V2_REPORT",v2,flush=True)
            return 0
        except Exception as exc:
            err=KEY+": "+str(exc)
            state["errors"]=[err]
            state["error_history"].append(err)
            state["status"]=STOP
            V1.A02.atomic(v2,state)
            print("V2_SGC_FAIL_CLOSED_CHECKPOINT_PRESERVED",err,flush=True)
            raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--include-48000",action="store_true",
                    help="Explicitly run sparse 48000R only after 4800R checkpoint")
    args=ap.parse_args()
    p,a02,e4=guards()
    if args.self_test:
        need(not args.run and not args.include_48000,
             "Self-test cannot read any original v1 report or WSL FITS")
        synthetic_test()
        return 0
    need(args.run,"Real local SGC must be requested with --run")
    return run(p,a02,e4,include_48000=args.include_48000)


if __name__=="__main__":
    raise SystemExit(main())
