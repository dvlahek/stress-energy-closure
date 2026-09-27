#!/usr/bin/env python3
"""Full-eligible MOCK0001 eBOSS LRGxELG sparse-pair diagnostic, NO OBSERVED DATA.

After complete re-SHA of the original 72 local gzip sources, replay original
600D/1200R sample hashes and 4800R E4 catalogues. Reproduce original dense
600D/1200R weighted pairs with independently scalar-tested sparse cKDTree
pairs; then use ALL source-eligible mock galaxies with ORIGINAL 4800 and
prespecified extended 48000 randoms per tracer/cap. No download or inference.

The purpose is to replace an underpowered 600D implementation pilot with a
full-mock estimator *engineering* test, not to certify the physical window.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import audit_eboss_dr16_a03_e4_nested_random_density as E4
import audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection as E3
from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import (
    oriented_pair_terms, rr_histogram,
)
from audit_eboss_dr16_rr_kdtree_pilot import kdtree_rr
from audit_eboss_dr16_ezmock0001_galaxy_odd_projection import project
from eboss_dr16_fiducial import (
    PRIMARY_GEOMETRY, WEIGHT_COLUMNS, comoving_mpc_over_h,
    validated_weight_product,
)
from check_eboss_cross_ls_synthetic import cross_landy_szalay

ROOT = A02.ROOT
PROTOCOL = ROOT / "source_data/eboss_dr16_full_eligible_mock_galaxy_4800_48000_random_pilot_protocol_2026-09-27.json"
OLD_A02 = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
OLD_E4 = ROOT / "source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json"
OLD_MANIFEST = ROOT / "source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json"
PROTOCOL_BLOB = "e53986db985f8bfcbe8e692eb9c3aed2d76ed5ed"
OLD_A02_SHA = "15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
OLD_E4_SHA = "ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f"
OLD_MANIFEST_BLOB = "8e73c473319963c2c1d8a1c9887226da889d2126"
OLD_E4_RUNNER_BLOB = "fda4a5eda420a9f69bd4d393bae02419155958f8"
SPARSE_RUNNER_BLOB = "eBOSS existing sparse source script source-pin from audit_eboss_dr16_rr_kdtree_pilot.py"
PASS = "EBOSS_FULL_ELIGIBLE_MOCK0001_BOTH_CAPS_4800_48000_SPARSE_DESCRIPTIVE_ONLY"
STOP = "EBOSS_FULL_ELIGIBLE_MOCK0001_INCOMPLETE_STOP"
SIZES = (4800, 48000)
LEVELS = ("original_4800","nested_48000")
EXTRA = 43200
NEW_SEED_ROOT = 2026092707
CAPS = A02.CAPS
TRACERS = A02.TRACERS
ROLES = A02.ROLES
MAPPING = A02.MAPPING
S_EDGES = A02.S_EDGES
MU_EDGES = A02.MU_EDGES
TERMS = ("D1D2","D1R2","R1D2","R1R2")
TOTAL_SOURCE_BYTES = 1860198719


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def blob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def guard(flag, msg):
    if not flag:
        raise ValueError(msg)


def catalogue_sha(cat):
    return E3.catalogue_sha(cat)


def source_only_gate():
    raw = PROTOCOL.read_bytes()
    guard(blob(raw) == PROTOCOL_BLOB, "Registered engineering protocol Git blob changed")
    p = json.loads(raw)
    a02raw = OLD_A02.read_bytes()
    e4raw = OLD_E4.read_bytes()
    guard(sha(a02raw) == OLD_A02_SHA and sha(e4raw) == OLD_E4_SHA,
          "Original parent 18-case source SHA changed")
    guard(blob(OLD_MANIFEST.read_bytes()) == OLD_MANIFEST_BLOB and
          blob(E4.__file__ and Path(E4.__file__).read_bytes()) == OLD_E4_RUNNER_BLOB,
          "Original E4 runner/manifest bytes changed")
    guard(p["branch"] == "audit/eboss-elg-bit8-ra-orientation-20260925" and
          p["fixed_id"] == [1] and p["caps"] == list(CAPS) and
          p["tracer_order"] == list(TRACERS) and p["roles"] == list(ROLES) and
          p["redshift"] == [.9,1.] and
          p["geometry"]["s_edges_mpc_h"] == S_EDGES.tolist() and
          p["geometry"]["signed_mu_bins"] == 24 and
          p["geometry"]["fiducial"] == PRIMARY_GEOMETRY and
          p["randoms"]["nested_levels"] == list(SIZES) and
          p["randoms"]["new_extra_count_each_tracer_cap"] == EXTRA and
          p["galaxies"]["galaxy_resampling"] is False and
          p["source_hash_requirements"]["all_72_full_original_gzip_SHA_reverified_before_any_new_FITS"] is True and
          p["execution"]["no_new_input_download"] is True and
          p["execution"]["no_observed_catalogue_row_read"] is True and
          p["interpretation"]["cannot_establish"].startswith("Historical ELG"),
          "Full-D engineering scope, fiducial or 48000R seed contract changed")
    for k in ("observed_galaxy_rows_read","observed_random_rows_read",
              "observed_odd_data_vector_read","new_science_cut",
              "additional_mock_downloads","author_contact","main_change",
              "unblinding_authorized"):
        guard(p[k] is False, "Forbidden full-D scope/observed data access: "+k)
    a02 = json.loads(a02raw)
    e4 = json.loads(e4raw)
    guard(a02["completed_cases"] == e4["completed_cases"] == 18 and
          a02["observed_odd_data_vector_read"] is False and
          e4["observed_odd_data_vector_read"] is False and
          a02["cases"]["0001/NGC"]["status"].startswith("oriented_pair") and
          e4["cases"]["0001/SGC"]["status"] == "complete",
          "Original source-only paired mock cohort incomplete")
    return p, a02, e4


def original_seed(p0, cap, tracer, role):
    return (p0["sampling"]["seed_root"] + 10000 * CAPS.index(cap)
            + 100 * TRACERS.index(tracer) + int(role == "ran"))


def copy_four(ra, dec, z, weights, eligible, positions):
    indices = eligible[positions]
    return tuple(np.asarray(v, dtype="f8").copy()
                 for v in (ra[indices], dec[indices],
                           z[indices], weights[positions]))


def load_original_full(path, cap, tracer, role, oldmeta, e4record, p0):
    # Must only be called after A02.resolve_72_sources() returns.
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tabs = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        guard(len(tabs) == 1 and
              int(tabs[0].header["NAXIS2"]) == oldmeta["input_FITS_rows"],
              "Frozen source full mock BINTABLE/header differs")
        tab = tabs[0]
        cols = set(tab.columns.names)
        needed = {"RA","DEC","Z",*WEIGHT_COLUMNS}
        if tracer == "eBOSS_ELG":
            needed.add("chunk")
        guard(needed.issubset(cols), "Original full mock source missing exact columns")
        d = tab.data
        ra,dec,z = (np.asarray(d[name],dtype="f8") for name in ("RA","DEC","Z"))
        guard(all(np.isfinite(v).all() for v in (ra,dec,z)) and
              np.all((ra>=0)&(ra<360)) and
              np.all((dec>=-90)&(dec<=90)),
              "Full mock source coordinate/nonfinite gate differs")
        candidate = (z>=.9)&(z<1.)
        weights, retained = validated_weight_product({
            name: np.asarray(d[name][candidate],dtype="f8")
            for name in WEIGHT_COLUMNS})
        eligible = np.flatnonzero(candidate)[retained]
        valid_weight = weights[retained]
        guard(int(candidate.sum()) == oldmeta["candidate_highz_rows_before_weight_validation"] and
              len(eligible) == oldmeta["eligible_highz_rows_after_fixed_weight_gate"] and
              int(candidate.sum())-len(eligible) == oldmeta["numerical_zero_weight_excluded_rows"] and
              np.isfinite(valid_weight).all() and np.all(valid_weight>0),
              "Frozen original full mock source eligible population/weights drift")
        nold = 600 if role=="dat" else 1200
        seed = original_seed(p0,cap,tracer,role)
        guard(oldmeta["seed"] == seed and oldmeta["selected_rows"] == nold,
              "Original 600D/1200R source sample seed contract differs")
        origpos = np.sort(np.random.default_rng(seed).choice(
            len(eligible),size=nold,replace=False))
        original = copy_four(ra,dec,z,valid_weight,eligible,origpos)
        guard(catalogue_sha(original) == oldmeta["selected_array_SHA256"] and
              catalogue_sha(original) == e4record[
                  "original_A02_all_four_selected_catalogue_SHA256"][tracer+"_"+role],
              "Full mock re-read does not reproduce original 600D/1200R SHA")
        full = copy_four(ra,dec,z,valid_weight,eligible,
                         np.arange(len(eligible),dtype=np.int64)) if role=="dat" else None
        info = {"full_eligible_rows":int(len(eligible)),
                "original_fixed_subsample_SHA256":catalogue_sha(original),
                "original_selected_rows":nold,
                "original_source_eligible_count_replayed":True}
        if role=="dat":
            info["full_eligible_four_vector_SHA256"]=catalogue_sha(full)
            if tracer=="eBOSS_ELG":
                labels=np.asarray(d["chunk"])[eligible]
                info["full_eligible_original_chunk_counts"]={
                    (k.decode("ascii") if isinstance(k,(bytes,np.bytes_)) else str(k)).strip():int(n)
                    for k,n in zip(*np.unique(labels,return_counts=True))}
            return {"original":original,"full":full}, info
        newseed = (NEW_SEED_ROOT+100000+1000*CAPS.index(cap)+
                   100*TRACERS.index(tracer))
        seed4,extra4,pos4 = E4.supplement_positions(
            len(eligible),origpos,1,cap,tracer)
        old4800=pos4["nested_4800"]
        orig4=copy_four(ra,dec,z,valid_weight,eligible,old4800)
        previous=e4record["fixed_random_extensions_by_tracer"][tracer]
        guard(seed4 == previous["fixed_supplement_PCG64_seed"] and
              E3.arr_sha(old4800)==previous["levels"]["nested_4800"][
                  "selected_eligible_positions_SHA256"] and
              catalogue_sha(orig4)==previous["levels"]["nested_4800"][
                  "four_vector_catalogue_SHA256"],
              "Original E4 4800R source vectors or indices changed")
        complement=np.setdiff1d(np.arange(len(eligible),dtype=np.int64),
                                old4800,assume_unique=True)
        guard(len(complement)>=EXTRA,"Not enough original eligible same-source randoms")
        chosen=complement[np.random.default_rng(newseed).choice(
            len(complement),size=EXTRA,replace=False)]
        fullpos=np.sort(np.r_[old4800,chosen])
        guard(len(fullpos)==48000 and len(np.unique(fullpos))==48000 and
              np.all(np.isin(old4800,fullpos)),
              "Original nested 4800R not preserved in 48000R")
        big=copy_four(ra,dec,z,valid_weight,eligible,fullpos)
        info.update({"original_E4_4800_catalogue_SHA256":catalogue_sha(orig4),
                     "new_43200_complement_seed":newseed,
                     "new_43200_original_eligible_indices_SHA256":E3.arr_sha(chosen),
                     "nested_48000_eligible_indices_SHA256":E3.arr_sha(fullpos),
                     "nested_48000_four_vector_SHA256":catalogue_sha(big)})
        return {"original":original,"original_4800":orig4,"nested_48000":big},info


def validate_sparse_against_original_dense(cats, oldcase, distance):
    dl,de=cats["eBOSS_LRG","dat"]["original"],cats["eBOSS_ELG","dat"]["original"]
    rl,re=cats["eBOSS_LRG","ran"]["original"],cats["eBOSS_ELG","ran"]["original"]
    originals={"D1D2":(dl,de),"D1R2":(dl,re),
               "R1D2":(rl,de),"R1R2":(rl,re)}
    results={}
    for side,samples in (
        ("forward_pair_terms",originals),
        ("reverse_pair_terms",{
            "D1D2":(de,dl),"D1R2":(de,rl),
            "R1D2":(re,dl),"R1R2":(re,rl)
        })):
        for term,(c1,c2) in samples.items():
            dense,meta=rr_histogram(c1,c2,S_EDGES,MU_EDGES,.05,distance,block=128)
            expected=oldcase["cross_ls"][side][term]
            guard(meta["accepted_pairs"] == expected["accepted_pairs"] and
                  sha(np.ascontiguousarray(dense).tobytes())==
                  expected["weighted_histogram_SHA256"],
                  "Original 600D/1200R dense pair histogram SHA reproduction failed: "+side+"/"+term)
            sparse,smeta=kdtree_rr(
                c1,c2,S_EDGES,MU_EDGES,distance,.05,chunk=64,
                max_neighbour_pairs=25000000)
            guard(smeta["accepted_weighted_pair_count"]==meta["accepted_pairs"] and
                  np.allclose(sparse,dense,atol=2e-9,rtol=2e-11) and
                  np.isclose(smeta["pair_normalization"],
                             meta["pair_normalization"],atol=0,rtol=1e-12),
                  "Sparse/dense exact original-weighted pair closure failed: "+side+"/"+term)
            results[side+"/"+term]={"accepted_pairs":meta["accepted_pairs"],
                         "original_dense_exact_SHA256":expected["weighted_histogram_SHA256"],
                         "sparse_dense_max_abs_weighted_hist_error":
                           float(np.max(np.abs(sparse-dense)))}
    return results


def cached_exact_distance(cats):
    # Compute original fiducial distance for each *exact* supplied redshift
    # only once per cap; no interpolated/fitted cosmological calibration.
    full=[cats[(tr,"dat")]["full"][2] for tr in TRACERS]
    full += [cats[(tr,"ran")]["nested_48000"][2] for tr in TRACERS]
    zz=np.unique(np.concatenate(full))
    rr=comoving_mpc_over_h(zz,PRIMARY_GEOMETRY)
    guard(np.isfinite(rr).all() and np.all(np.diff(rr)>0),
          "Primary eBOSS exact input-redshift distance mapping invalid")
    def distance(z):
        vals=np.asarray(z,dtype="f8")
        ii=np.searchsorted(zz,vals)
        guard(np.all(ii < len(zz)) and np.array_equal(zz[ii],vals),
              "Distance cache hit with unverified/redshift-resampled row")
        return rr[ii]
    return distance


def full_pairs(cats,level,distance,*,max_candidates):
    dl,de=cats["eBOSS_LRG","dat"]["full"],cats["eBOSS_ELG","dat"]["full"]
    rl,re=cats["eBOSS_LRG","ran"][level],cats["eBOSS_ELG","ran"][level]
    out={}
    normalized={}
    for side,samples in (
        ("forward",{"D1D2":(dl,de),"D1R2":(dl,re),
                    "R1D2":(rl,de),"R1R2":(rl,re)}),
        ("reverse",{"D1D2":(de,dl),"D1R2":(de,rl),
                    "R1D2":(re,dl),"R1R2":(re,rl)})):
        out[side]={}
        arrays={}
        norms={}
        for name,(first,second) in samples.items():
            h,meta=kdtree_rr(first,second,S_EDGES,MU_EDGES,
                             distance,.05,chunk=64,
                             max_neighbour_pairs=max_candidates)
            denom=meta["pair_normalization"]
            guard(h.shape==(6,24) and
                  np.all(np.isfinite(h)) and np.all(h>=0) and
                  denom>0, "Nonfinite sparse full-D pair histogram")
            out[side][name]={"accepted_pairs":meta["accepted_weighted_pair_count"],
                "candidate_neighbour_pairs":meta["candidate_neighbour_pairs"],
                "pair_normalization":denom,
                "histogram_6x24":h.tolist(),
                "histogram_SHA256":sha(np.ascontiguousarray(h).tobytes())}
            arrays[name]=h
            norms[name]=denom
            print("FULL_MOCK_PAIR",level,side,name,
                  "accepted",meta["accepted_weighted_pair_count"],
                  "candidates",meta["candidate_neighbour_pairs"],flush=True)
        xi,support=cross_landy_szalay(arrays,norms)
        normalized[side]={"xi":xi,"support":support,
                          "normalized_pair_hist":{k:arrays[k]/norms[k] for k in TERMS}}
    fw,rv=normalized["forward"],normalized["reverse"]
    guard(np.array_equal(fw["support"],rv["support"][:,::-1]),
          "Full-D sparse pair reverse positive RR masks differ")
    for term,reverse in MAPPING.items():
        normed=normalized["forward"]["normalized_pair_hist"][term]
        rev=normalized["reverse"]["normalized_pair_hist"][reverse][:,::-1]
        guard(np.allclose(normed,rev,atol=1e-11,rtol=2e-10) and
              out["forward"][term]["accepted_pairs"]==
              out["reverse"][reverse]["accepted_pairs"],
              "Full-D sparse forward reverse four-term parity failed: "+term)
    guard(fw["xi"].shape==(6,24) and
          np.count_nonzero(fw["support"])==144 and
          np.isfinite(fw["xi"]).all() and
          np.allclose(fw["xi"],rv["xi"][:,::-1],atol=5e-9,rtol=2e-9),
          "Full-D full-RR-support or true independent signed-mu xi reversal failed")
    out["xi_6x24"]=fw["xi"].tolist()
    out["xi_SHA256"]=sha(np.ascontiguousarray(fw["xi"]).tobytes())
    out["RR_supported_cells"]=int(np.count_nonzero(fw["support"]))
    out["max_abs_reverse_xi"]=float(np.max(
        np.abs(fw["xi"]-rv["xi"][:,::-1])))
    out["odd_and_even_multipoles"]=project(fw["xi"],fw["support"],MU_EDGES,(0,1,2,3))
    out["no_observed_galaxies_or_odd_read"]=True
    return out


def synthetic_test():
    rng=np.random.default_rng(20260927)
    n,m=37,43
    one=(rng.uniform(180,182,n),rng.uniform(5,6,n),
         rng.uniform(.91,.99,n),rng.uniform(.6,1.9,n))
    two=(rng.uniform(180,182,m),rng.uniform(5,6,m),
         rng.uniform(.91,.99,m),rng.uniform(.7,1.4,m))
    testdist=lambda z: np.asarray(z,dtype="f8")*2800.
    edges=np.array([0.,20.,40.,60.,80.,100.,140.])
    direct,old=rr_histogram(one,two,edges,MU_EDGES,.05,testdist,block=16)
    sparse,new=kdtree_rr(one,two,edges,MU_EDGES,testdist,.05,chunk=11,
                         max_neighbour_pairs=100000)
    guard(np.allclose(direct,sparse,rtol=1e-11,atol=1e-9) and
          old["accepted_pairs"]==new["accepted_weighted_pair_count"],
          "Synthetic sparse full-D pair kernel differs from original dense")
    malformed=np.arange(1200,dtype=np.int64)
    try:
        E4.supplement_positions(4000,malformed,1,"NGC","eBOSS_LRG")
    except ValueError:
        pass
    else:
        raise AssertionError("Original E4 eligible complement fail-closed guard broken")
    print("FULL_MOCK_SPARSE_VS_ORIGINAL_DENSE_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()
    p,a02,e4=source_only_gate()
    if args.self_test:
        synthetic_test()
        return 0
    guard(args.run and args.max_candidate_pairs_per_term>=1000000,
          "Only --self-test or explicit --run mock-only is supported")
    dest=ROOT/p["execution"]["output_json"]
    guard(not dest.exists(), "Original full-D output exists; preserve previous success/failure bytes")
    dest.parent.mkdir(parents=True,exist_ok=True)
    # Full source provenance BEFORE FIRST new FITS header or row.
    original_protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,original,gproto,rproto,rawproto,p0= A02.preflight(
        original_protocol,require_local=True)
    paths,total=A02.resolve_72_sources(
        original_protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    guard(len(paths)==72 and total==TOTAL_SOURCE_BYTES,
          "Full 72 source rehash/source byte count incomplete")
    state={"status":STOP,"source_protocol_blob":PROTOCOL_BLOB,
           "original_A02_sha256":OLD_A02_SHA,"original_E4_sha256":OLD_E4_SHA,
           "all_72_full_gzip_rehashed_before_any_FITS":True,
           "source_72_total_compressed_bytes":total,
           "full_eligible_galaxy_no_new_cut":True,
           "mock_id":1,"caps":list(CAPS),"levels":list(LEVELS),
           "cases":{},"errors":[],"observed_galaxy_rows_read":False,
           "observed_random_rows_read":False,
           "observed_odd_data_vector_read":False,
           "physical_window_certified":False,
           "inferential_covariance_computed":False,
           "wake_detection_or_exclusion_computed":False,
           "new_source_download":False}
    A02.atomic(dest,state)
    for cap in CAPS:
        key=f"0001/{cap}"
        metadata=a02["cases"][key]["input_sample_diagnostics"]
        prior=e4["cases"][key]
        cats={}
        info={}
        try:
            for tracer in TRACERS:
                for role in ROLES:
                    k=tracer+"_"+role
                    src=paths[A02.key(1,cap,tracer)+"/"+role]
                    cats[(tracer,role)],info[k]=load_original_full(
                        src,cap,tracer,role,metadata[k],prior,p0)
                    print("FULL_MOCK_INPUT_SHA_OK",key,k,
                          "full_eligible",info[k]["full_eligible_rows"],flush=True)
            # The production cache is built on exactly the same full-D/48k
            # catalogue redshifts and reused by both 4800 and 48000 runs.
            dist=cached_exact_distance(cats)
            benchmark=validate_sparse_against_original_dense(
                cats,a02["cases"][key],dist)
            print("FULL_MOCK_SPARSE_DENSE_8PAIR_ORIGINAL_SHA_REPLAY_OK",
                  key,flush=True)
            case={"mock_id":1,"cap":cap,"original_four_source_SHA_replayed":True,
                  "all_original_600D_1200R_eight_weighted_hist_SHAs_replayed":True,
                  "sparse_dense_eight_pair_closure":benchmark,
                  "source_full_eligible_diagnostics":info,"levels":{}}
            for level in LEVELS:
                levelout=full_pairs(
                    cats,level,dist,max_candidates=args.max_candidate_pairs_per_term)
                if level=="original_4800":
                    oldrr=prior["levels"]["nested_4800"]["baseline_original"][
                        "forward_pair_terms"]["R1R2"]
                    got=levelout["forward"]["R1R2"]
                    guard(got["accepted_pairs"]==oldrr["accepted_pairs"] and
                          np.isclose(np.sum(got["histogram_6x24"]),
                                     oldrr["total_weighted_pairs_in_fixed_s_mu_bins"],
                                     rtol=5e-10,atol=2e-8),
                          "Same original E4 4800R source RR accepted counts/weight changed")
                    levelout["original_E4_same_RR_4800_count_and_sum_replayed"]=True
                case["levels"][level]=levelout
            left=np.asarray(case["levels"]["nested_48000"]["xi_6x24"])
            right=np.asarray(case["levels"]["original_4800"]["xi_6x24"])
            denom=float(np.sum(np.abs(left)))
            case["full_galaxy_4800_to_48000_relative_xi_L1"] = (
                float(np.sum(np.abs(left-right))/denom) if denom>0 else None)
            case["full_galaxy_4800_to_48000_odd_multipole_differences"]={
                str(ell):[float(x-y) for x,y in zip(
                    case["levels"]["nested_48000"]["odd_and_even_multipoles"][str(ell)]["values_by_fixed_s_bin"],
                    case["levels"]["original_4800"]["odd_and_even_multipoles"][str(ell)]["values_by_fixed_s_bin"])]
                for ell in (1,3)}
            case["scope"]="ONE_FULL_GALAXY_MOCK_ONLY_ENGINEERING_NOT_PHYSICAL_SENSITIVITY"
            state["cases"][key]=case
            state["case_canonical_SHA256"]=state.get("case_canonical_SHA256",{})
            state["case_canonical_SHA256"][key]=sha(json.dumps(
                case,sort_keys=True,separators=(",",":"),allow_nan=False).encode())
            A02.atomic(dest,state)
            print("FULL_MOCK_CAP_COMPLETE",key,"xi_L1_4800_TO_48000",
                  case["full_galaxy_4800_to_48000_relative_xi_L1"],flush=True)
        except Exception as exc:
            state["errors"].append(key+": "+str(exc))
            state["cases"][key]={"status":"incomplete_stop",
                                "error":str(exc),"original_inputs_checked":list(info)}
            A02.atomic(dest,state)
            print("FULL_MOCK_FAIL_CLOSED",key,str(exc),flush=True)
            raise
    state["completed_cases"]=2
    state["status"]=PASS if not state["errors"] else STOP
    A02.atomic(dest,state)
    print("EBOSS_FULL_ELIGIBLE_MOCK0001",state["status"],flush=True)
    print("REPORT",dest,flush=True)
    print("OBSERVED_ODD_DATA_READ",False,flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
