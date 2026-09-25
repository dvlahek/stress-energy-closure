#!/usr/bin/env python3
"""Data-blind, exact-chunk random-only audit of three official ELG extra polygons.

This checks sample membership of published DR16 ELG *clustering RANDOMS*
in three official additional MANGLE polygons. Inputs and diagnostic
strata are frozen. No veto is reapplied to released clustering samples.
ELG BRICKMASK images, its extra HEALPix bit 8, LRG masks, physical
LRG×ELG window and galaxy odd signal are deliberately out of scope.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from aggregate_eboss_dr16_9mock_window import IDS, CAPS
from audit_eboss_dr16_9mock_fine_weighted_nz import (
    acquire, chunk_label, json_write_atomic, preflight,
)
from audit_eboss_dr16_official_polygon_sha import (
    digest_file, stream_pinned,
)
from eboss_dr16_fiducial import WEIGHT_COLUMNS, validated_weight_product
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_extra_random_membership_protocol_2026-09-25.json"
EVIDENCE = ROOT / "source_data/eboss_dr16_elg_three_extra_polygon_sha_2026-09-25.json"
UPSTREAM = ROOT / "scripts/eBOSS_ELG_extra.py"
SAMPLE_PER_STRATUM = 2500
SEED_BASE = 91307
BINS = ((0.6,0.7),(0.7,0.8),(0.8,0.9),(0.9,1.0))
CHUNKS = {"NGC": ("eboss23", "eboss25"),
          "SGC": ("eboss21", "eboss22")}
MASK_WORDS = (
    ("ELG_centerpost.ply", 1 << 9),
    ("ELG_TDSSFES_62arcsec.pix.snap.balk.ply", 1 << 10),
    ("ebosselg_badphot.26Aug2019.ply", 1 << 11),
)
EXPECTED_UPSTREAM_BLOB = "d66c5e6ed6a8260b5e8dcc840795485345392ced"


def git_blob_sha(path):
    b = path.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode() + b).hexdigest()


def validate_protocol_and_evidence(protocol, evidence, pinned_manifest):
    if (protocol.get("sample_per_stratum") != SAMPLE_PER_STRATUM
            or protocol.get("seed_base") != SEED_BASE
            or protocol.get("candidate_z_bins")
            != [list(x) for x in BINS]
            or tuple(protocol.get("mock_ids",[])) != IDS
            or tuple(protocol.get("caps",[])) != CAPS
            or {k:tuple(v) for k,v in protocol.get("chunks_by_cap",{}).items()} != CHUNKS
            or protocol.get("observed_odd_data_vector_read") is not False
            or evidence.get("status") != "official_elg_extra_polygon_bytes_pinned_only"
            or evidence.get("workflow_run") != 36100783598
            or evidence.get("source_upstream_commit")
            not in (None, "b9eb684a579b56ec3dbdb46549224be7e3fa2830")
            or evidence.get("observed_odd_data_vector_read") is not False
            or git_blob_sha(UPSTREAM) != EXPECTED_UPSTREAM_BLOB):
        raise ValueError("Fixed three-polygon membership registration invalid")
    keys = {
        row["filename"]: row
        for row in evidence.get("products", [])
    }
    report = {
        row["filename"]: row
        for row in pinned_manifest.get("products", [])
    }
    if (pinned_manifest.get("status")
            != "official_elg_extra_polygon_bytes_pinned_only"
            or pinned_manifest.get("observed_odd_data_vector_read") is not False
            or pinned_manifest.get("physical_joint_mask_certified") is not False
            or set(keys) != {x[0] for x in MASK_WORDS}
            or set(report) != set(keys)):
        raise ValueError("Pinned official extra-veto polygon manifest invalid")
    for name, bit in MASK_WORDS:
        x, y = keys[name], report[name]
        if (x["bit"] != bit.bit_length()-1
                or x["sha256"] != y["sha256"]
                or x["bytes"] != y["bytes"]
                or y["role"] != "elg_extra_polygon"
                or y["source_url"] != evidence["source_root"] + name):
            raise ValueError("Official ELG polygon version, bit or SHA changed: "+name)
    return keys, report


def acquire_polygon_sources(evidence, report, polygon_dir, *, timeout,
                            no_download):
    paths = {}
    root = evidence["source_root"]
    for name,_ in MASK_WORDS:
        entry = report[name]
        path = polygon_dir / name
        if path.exists():
            sha, nbytes = digest_file(path)
            if sha != entry["sha256"] or nbytes != entry["bytes"]:
                raise ValueError("Cached official ELG polygon SHA256 mismatch: "+name)
        else:
            if no_download:
                raise FileNotFoundError("Missing pinned official extra mask: "+str(path))
            sha,nbytes = stream_pinned(entry["source_url"],path,
                                      root,name,timeout)
            if sha != entry["sha256"] or nbytes != entry["bytes"]:
                path.unlink(missing_ok=True)
                raise ValueError("Official extra polygon bytes differ from pinned source: "+name)
        paths[name] = path
        print("ELG_EXTRA_POLYGON_SHA_OK",name,entry["sha256"],nbytes,flush=True)
    return paths


def stratum_seed(kind, mid, cap, chunk, zi):
    mid = 0 if kind == "observed" else mid
    return (SEED_BASE + 100_000_000 * CAPS.index(cap)
            + 10000 * mid + 100 * CHUNKS[cap].index(chunk) + zi)


def probe_membership(ra, dec, weights, polygons):
    if (len(ra) != len(dec) or len(ra) != len(weights) or not len(ra)
            or not np.isfinite(ra).all() or not np.isfinite(dec).all()
            or not np.isfinite(weights).all() or np.any(weights<=0)):
        raise ValueError("Invalid selected random-only polygon sample")
    bits = np.zeros(len(ra), dtype="i8")
    by_polygon = {}
    for name, bit in MASK_WORDS:
        ids = np.asarray(polygons[name].polyid(ra,dec))
        if ids.shape != (len(ra),):
            raise ValueError("Unexpected official polygon membership shape")
        inside = ids != -1
        bits[inside] |= bit
        by_polygon[name] = {
            "source_bit": int(bit),
            "sample_inside_rows": int(np.count_nonzero(inside)),
            "sample_inside_row_fraction": float(np.mean(inside)),
            "sample_inside_weight": float(np.sum(weights[inside],dtype="f8")),
            "sample_inside_weight_fraction": float(
                np.sum(weights[inside],dtype="f8") / np.sum(weights,dtype="f8")),
        }
    union = bits != 0
    values,count = np.unique(bits,return_counts=True)
    hist = {str(int(b)):int(n) for b,n in zip(values,count)}
    weighted_hist = {
        str(int(b)):float(np.sum(weights[bits==b],dtype="f8"))
        for b in values
    }
    if sum(hist.values()) != len(ra) or sum(hist.values()) != len(weights):
        raise ValueError("ELG extra polygon mask partition inconsistent")
    return {
        "polygon_diagnostic_membership": by_polygon,
        "union_sample_inside_rows": int(np.count_nonzero(union)),
        "union_sample_inside_weight": float(np.sum(weights[union],dtype="f8")),
        "union_sample_inside_row_fraction": float(np.mean(union)),
        "union_sample_inside_weight_fraction": float(
            np.sum(weights[union],dtype="f8") / np.sum(weights,dtype="f8")),
        "membership_maskword_counts": hist,
        "membership_maskword_weight_sums": weighted_hist,
    }


def random_strata(path, *, cap, kind, mid, expected_rows, polygons):
    selected = {}
    with fits.open(path,memmap=False) as hdus:
        hdus.verify("exception")
        tables=[h for h in hdus if isinstance(h,fits.BinTableHDU)]
        if len(tables)!=1:
            raise ValueError("Expected exactly one released random FITS BINTABLE")
        h=tables[0]
        if expected_rows is not None and h.header["NAXIS2"] != expected_rows:
            raise ValueError("Released ELG random FITS header row count changed")
        if not {"RA","DEC","Z","chunk",*WEIGHT_COLUMNS}.issubset(h.columns.names):
            raise ValueError("Expected published DR16 ELG chunk and weight columns")
        data=h.data
        z=np.asarray(data["Z"],dtype="f8")
        ra=np.asarray(data["RA"],dtype="f8")
        dec=np.asarray(data["DEC"],dtype="f8")
        if (not np.isfinite(z).all() or not np.isfinite(ra).all()
                or not np.isfinite(dec).all() or np.any(ra<0)
                or np.any(ra>=360) or np.any(dec<-90) or np.any(dec>90)):
            raise ValueError("Invalid published ELG random FITS coordinates")
        candidate=(z>=BINS[0][0])&(z<BINS[-1][1])
        fields={n:np.asarray(data[n][candidate],dtype="f8")
                for n in WEIGHT_COLUMNS}
        weight,gate=validated_weight_product(fields)
        zv=z[candidate][gate]
        rav=ra[candidate][gate]
        decv=dec[candidate][gate]
        wv=weight[gate]
        labels=np.asarray([
            chunk_label(v) for v in data["chunk"][candidate][gate]
        ])
        if set(np.unique(labels)) != set(CHUNKS[cap]):
            raise ValueError("Unexpected or missing exact published ELG chunk label in "+cap)
        for chunk in CHUNKS[cap]:
            for zi,(zlo,zhi) in enumerate(BINS):
                rows=np.flatnonzero((labels==chunk)&(zv>=zlo)&(zv<zhi))
                if len(rows)<SAMPLE_PER_STRATUM:
                    raise ValueError(
                        "Fixed source sample below 2500 in "+cap+"/"+chunk+"/"+str(zi))
                seed=stratum_seed(kind,mid,cap,chunk,zi)
                subset=rows[np.random.default_rng(seed).permutation(len(rows))[
                    :SAMPLE_PER_STRATUM]]
                item=probe_membership(
                    rav[subset],decv[subset],wv[subset],polygons)
                selected[f"{chunk}_z{zi}"]={
                    "cap":cap,"chunk":chunk,"zbin":[zlo,zhi],
                    "retained_parent_rows":int(len(rows)),
                    "seed":int(seed),"sampled_rows":SAMPLE_PER_STRATUM,
                    "sample_weight_sum":float(np.sum(wv[subset],dtype="f8")),
                    **item,
                }
        retained_count=len(zv)
    if len(selected)!=8:
        raise ValueError("Expected eight exact chunk-by-z selection strata")
    print("ELG_EXTRA_RANDOM_MEMBERSHIP_INPUT_OK",kind,mid,cap,
          retained_count,len(selected),flush=True)
    return selected,retained_count


def load_polygons(paths):
    import pymangle
    polygons={name:pymangle.Mangle(str(path)) for name,path in paths.items()}
    if set(polygons)!={x[0] for x in MASK_WORDS}:
        raise ValueError("Incomplete official extra mask polygon readers")
    return polygons


def self_test():
    p=json.loads(PROTOCOL.read_text())
    evidence=json.loads(EVIDENCE.read_text())
    assert p["sample_per_stratum"]==SAMPLE_PER_STRATUM
    assert len(evidence["products"])==3
    assert git_blob_sha(UPSTREAM)==EXPECTED_UPSTREAM_BLOB
    class Toy:
        def __init__(self,bound):
            self.bound=bound
        def polyid(self,ra,dec):
            return np.where(ra>self.bound,12,-1)
    toy={
        name:Toy(bound)
        for (name,_),bound in zip(MASK_WORDS,(0.5,1.5,2.5))
    }
    r=probe_membership(
        np.array([0.,1.,2.,3.]),np.zeros(4),np.ones(4),toy)
    assert r["union_sample_inside_rows"]==3
    assert r["polygon_diagnostic_membership"][
        "ELG_centerpost.ply"]["sample_inside_rows"]==3
    assert r["membership_maskword_counts"]=={
        "0":1,"512":1,"1536":1,"3584":1}
    assert stratum_seed("observed",0,"SGC","eboss21",3)!=stratum_seed(
        "mock",1,"SGC","eboss21",3)
    print("EBOSS_ELG_EXTRA_RANDOM_MEMBERSHIP_SELF_TEST_OK",flush=True)


def run(args):
    proto=json.loads(PROTOCOL.read_text())
    evidence=json.loads(EVIDENCE.read_text())
    pinned=json.loads(args.polygon_manifest.read_text())
    expected,_=validate_protocol_and_evidence(proto,evidence,pinned)
    paths=acquire_polygon_sources(
        evidence,_,args.polygon_dir,
        timeout=args.timeout,no_download=args.no_download)
    polygons=load_polygons(paths)
    mocksha,obssha,obs_counts,_=preflight(args)
    mids=IDS if args.mock_ids=="all" else (
        () if args.mock_ids=="none" else
        tuple(int(x.strip()) for x in args.mock_ids.split(",")))
    caps=CAPS if args.caps=="all" else tuple(
        x.strip() for x in args.caps.split(","))
    if (len(caps)==0 or any(c not in CAPS for c in caps)
            or len(set(caps))!=len(caps)
            or any(mid not in IDS for mid in mids)
            or len(set(mids))!=len(mids)):
        raise ValueError("Only fixed mock IDs and both declared caps are valid")
    out=args.out_dir
    report={
        "study":"Data-blind published extra ELG polygon random membership",
        "status":"elg_extra_polygon_membership_incomplete",
        "protocol":str(PROTOCOL.relative_to(ROOT)),
        "polygon_sha_evidence":str(EVIDENCE.relative_to(ROOT)),
        "sample_per_chunk_redshift_stratum":SAMPLE_PER_STRATUM,
        "requested_caps":list(caps),"requested_mock_ids":list(mids),
        "cases":[],"errors":[],
        "official_elg_bit8_checked":False,
        "full_elg_brickmask_checked":False,
        "lrg_official_veto_checked":False,
        "joint_pair_window_certified":False,
        "observed_galaxy_data_read":False,
        "mock_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
    }
    for cap in caps:
        for kind,mid in [("observed",0)]+[("mock",i) for i in mids]:
            source_sha=obssha[cap,"ELG"] if kind=="observed" else mocksha[mid,cap,"ELG"]
            source_name=(RANDOMS["ELG",cap][0] if kind=="observed" else
                         Path(mock_path("eBOSS_ELG",cap,"ran",mid)).name)
            audit_path=out/f"{kind}_{mid:04d}_{cap}_elg_extra_membership.json"
            if audit_path.exists():
                record=json.loads(audit_path.read_text())
                if (record.get("status")!="elg_extra_polygon_random_case_complete"
                    or record.get("source_sha256")!=source_sha
                    or record.get("polygon_sha256")
                    !={k:v["sha256"] for k,v in expected.items()}
                    or record.get("observed_odd_data_vector_read") is not False):
                    raise ValueError("Stored mask membership checkpoint input differs: "+str(audit_path))
            else:
                local=(args.observed_cache_dir if kind=="observed"
                       else args.mock_cache_dir)/source_name
                relative=None if kind=="observed" else mock_path(
                    "eBOSS_ELG",cap,"ran",mid)
                existed=local.is_file()
                local,digest,size=acquire(
                    local,source_sha,observed=kind=="observed",
                    filename=source_name if kind=="observed" else None,
                    url=MOCK_BASE+relative if relative else None,
                    timeout=args.timeout,no_download=args.no_download)
                expected_rows=RANDOMS["ELG",cap][1] if kind=="observed" else None
                try:
                    strata,nrows=random_strata(
                        local,cap=cap,kind=kind,mid=mid,
                        expected_rows=expected_rows,polygons=polygons)
                    record={
                        "status":"elg_extra_polygon_random_case_complete",
                        "kind":kind,"mock_id":mid if kind=="mock" else None,
                        "cap":cap,"source_sha256":digest,
                        "source_file_bytes":size,
                        "source_filename":source_name,
                        "retained_candidate_rows":nrows,
                        "polygon_sha256":{
                            k:v["sha256"] for k,v in expected.items()},
                        "strata":strata,
                        "observed_galaxy_data_read":False,
                        "mock_galaxy_data_read":False,
                        "observed_odd_data_vector_read":False,
                        "physical_joint_mask_certified":False,
                    }
                    json_write_atomic(audit_path,record)
                finally:
                    if (kind=="mock" and args.purge_new_mock_fits
                            and not existed):
                        local.unlink(missing_ok=True)
            report["cases"].append({
                "kind":kind,"mock_id":mid if kind=="mock" else None,
                "cap":cap,"record_path":str(audit_path),
                "strata_count":len(record["strata"]),
                "union_inside_total_sampled_rows":sum(
                    x["union_sample_inside_rows"]
                    for x in record["strata"].values()),
                "sampled_total_rows":sum(
                    x["sampled_rows"] for x in record["strata"].values()),
            })
            report["status"]="elg_extra_polygon_membership_checkpoint"
            json_write_atomic(out/"elg_extra_membership_summary.json",report)
            gc.collect()
    full=(set(caps)==set(CAPS) and set(mids)==set(IDS))
    report["status"]=("elg_extra_polygon_random_9mock_cohort_complete"
                      if full else "elg_extra_polygon_random_pilot_complete")
    report["note"]=(
        "Polygon membership on fixed-seed sampled released ELG clustering "
        "RANDOMS only. No veto reapplication; official extra HEALPix bit8, "
        "all brickmask images and full LRG×ELG physical window remain open.")
    json_write_atomic(out/"elg_extra_membership_summary.json",report)
    print("EBOSS_ELG_EXTRA_RANDOM_MEMBERSHIP",report["status"],
          len(report["cases"]),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ensemble-json",type=Path)
    ap.add_argument("--polygon-manifest",type=Path)
    ap.add_argument("--polygon-dir",type=Path,
                    default=Path("eboss_workspace/official_mask_inventory/official_polygons"))
    ap.add_argument("--observed-cache-dir",type=Path,
                    default=Path("eboss_workspace/local_rr/fits"))
    ap.add_argument("--mock-cache-dir",type=Path,
                    default=Path("eboss_workspace/local_nz/mock_fits"))
    ap.add_argument("--out-dir",type=Path,
                    default=Path("eboss_workspace/official_mask_inventory/membership"))
    ap.add_argument("--caps",default="SGC",
                    help="One of NGC, SGC or all")
    ap.add_argument("--mock-ids",default="none",
                    help="none, comma-separated fixed nine IDs, or all")
    ap.add_argument("--timeout",type=float,default=120)
    ap.add_argument("--no-download",action="store_true")
    ap.add_argument("--purge-new-mock-fits",action="store_true")
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (not args.ensemble_json or not args.polygon_manifest
            or args.timeout<=0):
        ap.error("Require --ensemble-json, --polygon-manifest and positive timeout")
    args.out_dir.mkdir(parents=True,exist_ok=True)
    try:
        run(args)
    except Exception as exc:
        json_write_atomic(args.out_dir/"elg_extra_membership_failure.json",{
            "status":"elg_extra_polygon_random_membership_failed",
            "error":str(exc),"observed_galaxy_data_read":False,
            "observed_odd_data_vector_read":False,
            "physical_joint_mask_certified":False,
        })
        print("EBOSS_ELG_EXTRA_RANDOM_MEMBERSHIP_FAILED",
              type(exc).__name__,str(exc),flush=True)
        return 2
    (args.out_dir/"elg_extra_membership_failure.json").unlink(missing_ok=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
