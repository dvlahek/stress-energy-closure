#!/usr/bin/env python3
"""SHA-first, blinded 9-ID realistic EZmock galaxy×own-random cross-LS audit.

Every one of 72 previously SHA-pinned compressed source files is fully
rehashed BEFORE the first new FITS header or galaxy/random row is opened.
The original fixed mock0001 sample arrays and all forward/reverse pair
histogram SHAs must reproduce before reporting other mock results.
The nine fixed mock IDs are a code-transport pilot, NOT an inferential
18D covariance, physical odd null, detection or observed odd analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from audit_eboss_dr16_ezmock0001_galaxy_raw_sha import (
    ROOT, digest, load, atomic,
)
from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import (
    sample_catalogue, oriented_pair_terms, PILOT as SINGLE_PROTOCOL,
)
from audit_eboss_dr16_rr_pair_closure import mirrored_closure
from check_eboss_cross_ls_synthetic import (
    LABELS, cross_landy_szalay, audit as synthetic_cross_ls,
)
from eboss_dr16_fiducial import (
    PRIMARY_GEOMETRY, SYSTOT_NUMERICAL_ZERO_TOL, comoving_mpc_over_h,
)
from inspect_eboss_dr16_mock_headers import BASE, mock_path

PROTOCOL = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_protocol_2026-09-26.json"
PASS = "EZMOCK_NINE_GALAXY_CROSS_LS_CODE_TRANSPORT_ALGEBRA_ONLY"
STOP = "EZMOCK_NINE_GALAXY_CROSS_LS_CODE_TRANSPORT_INCOMPLETE_STOP"
IDS = (1,125,250,375,500,625,750,875,1000)
CAPS = ("NGC","SGC")
TRACERS = ("eBOSS_LRG","eBOSS_ELG")
ROLES = ("dat","ran")
S_EDGES = np.asarray([20.,40.,60.,80.,100.,120.,140.],dtype="f8")
MU_EDGES = np.linspace(-1. - 1e-7,1. + 1e-7,25,dtype="f8")
MAPPING = {"D1D2":"D1D2","D1R2":"R1D2","R1D2":"D1R2","R1R2":"R1R2"}
G_REPORT = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_report_2026-09-26.json"
R_REPORT = ROOT / "source_data/eboss_dr16_nine_ezmock_matched_random_sha_report_2026-09-26.json"
G_MANIFEST = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_uploaded_manifest_2026-09-26.json"
R_MANIFEST = ROOT / "source_data/eboss_dr16_nine_ezmock_matched_random_sha_uploaded_manifest_2026-09-26.json"
REFERENCE = ROOT / "source_data/eboss_dr16_nine_ezmock_random_reference_sha_from_20260924_artifact.json"
G_PROTOCOL = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_protocol_2026-09-26.json"
R_PROTOCOL = ROOT / "source_data/eboss_dr16_nine_ezmock_matched_random_sha_protocol_2026-09-26.json"
ORIG_RAW = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_raw_bytes_protocol_2026-09-26.json"
ORIG_PAIR = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_report_2026-09-26.json"
ORIG_PAIR_MANIFEST = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_uploaded_manifest_2026-09-26.json"
ALL = tuple((mid,cap,tracer,role) for mid in IDS for cap in CAPS
            for tracer in TRACERS for role in ROLES)


def sha_bytes(path, expected, *, local=None):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError("Frozen exact uploaded source-report SHA mismatch: " + str(path))
    if local is not None and local.read_bytes()!=raw:
        raise ValueError("Local uploaded source checkpoint differs byte-for-byte from archived source report: "+str(local))
    return json.loads(raw)


def key(mid,cap,tracer):
    return f"{mid:04d}/{cap}/{tracer}"


def check_maps(ref,g,r,gmanifest,rmanifest):
    """Reusable strict source-key+SHA consistency gate; synthetic tamper tested."""
    keys=tuple(key(mid,cap,tracer) for mid in IDS for cap in CAPS for tracer in TRACERS)
    if (
        tuple(g["existing_frozen_0001_sources"])!=keys[:4]
        or tuple(g["new_galaxy_sources"])!=keys[4:]
        or tuple(r["verified_random_sources"])!=keys
        or tuple(ref["full_random_gzip_sha256_by_id_cap_tracer"])!=keys
        or rmanifest["complete_random_source_count"]!=36
        or gmanifest["count_total"]!=36
    ):
        raise ValueError("Frozen 72 unique cap/tracer/ID source keys differ")
    for mid,cap,tracer in (entry[:3] for entry in ALL if entry[3]=="dat"):
        k=key(mid,cap,tracer)
        gr=(g["existing_frozen_0001_sources"][k] if mid==1
            else g["new_galaxy_sources"][k])
        gp=(gmanifest["previously_frozen_0001_sha"][k] if mid==1
            else gmanifest["additional_eight_ids_full_compressed_sha"][k])
        rs=r["verified_random_sources"][k]
        rp=rmanifest["verified_random_sources_by_id_cap_tracer"][k]
        reference=ref["full_random_gzip_sha256_by_id_cap_tracer"][k]
        gsha=(gr["full_compressed_sha256_reverified"] if mid==1
              else gr["full_compressed_sha256_first_seen"])
        gpin=(gp["full_compressed_sha256_reverified"] if mid==1
              else gp["full_compressed_sha256"])
        gsize=gr["compressed_bytes"]
        if (
            gsha!=gpin or gsize!=gp["compressed_bytes"]
            or gr["FITS_headers_or_mock_rows_read"] is not False
            or rs["fully_reverified_compressed_sha256"]!=reference
            or rs["prior_20260924_artifact_expected_sha256"]!=reference
            or rp["previously_pinned_sha256"]!=reference
            or rs["compressed_bytes"]!=rp["compressed_bytes"]
            or rs["sha_matches_original_artifact_reference"] is not True
            or rs["gzip_decompressed_or_FITS_rows_read"] is not False
        ):
            raise ValueError("Frozen galaxy/random/reference source SHA, size or scope mismatch: "+k)
        grel=mock_path(tracer,cap,"dat",mid)
        rrel=mock_path(tracer,cap,"ran",mid)
        if (
            gp.get("released_filename",Path(grel).name)!=Path(grel).name
            or gr.get("release_source_url",gr.get("official_source_url"))!=BASE+grel
            or rs["official_source_url"]!=BASE+rrel
            or rs["relative_released_filename"]!=rrel
        ):
            raise ValueError("Source filename/official URL changed: "+k)
    return keys


def preflight(p, *, require_local):
    gproto=load(G_PROTOCOL)
    rproto=load(R_PROTOCOL)
    origproto=load(SINGLE_PROTOCOL)
    rawproto=load(ORIG_RAW)
    gman=load(G_MANIFEST)
    rman=load(R_MANIFEST)
    ref=load(REFERENCE)
    origmanifest=load(ORIG_PAIR_MANIFEST)
    g=sha_bytes(G_REPORT,p["nine_galaxy_input_source_report_sha256"],
                local=(ROOT/p["nine_galaxy_input_source_report_local"]) if require_local else None)
    r=sha_bytes(R_REPORT,p["nine_random_input_source_report_sha256"],
                local=(ROOT/p["nine_random_input_source_report_local"]) if require_local else None)
    old=sha_bytes(ORIG_PAIR,p["original_0001_cross_ls_report_SHA256"])
    expected_old=origmanifest["exact_uploaded_sha256"]
    if (
        expected_old!=p["original_0001_cross_ls_report_SHA256"]
        or p["original_0001_cross_ls_protocol"]!=str(SINGLE_PROTOCOL.relative_to(ROOT))
        or p["nine_galaxy_input_source_report_archive"]!=str(G_REPORT.relative_to(ROOT))
        or p["nine_random_input_source_report_archive"]!=str(R_REPORT.relative_to(ROOT))
        or p["nine_galaxy_input_source_manifest"]!=str(G_MANIFEST.relative_to(ROOT))
        or p["nine_random_input_source_manifest"]!=str(R_MANIFEST.relative_to(ROOT))
        or p["original_20260924_random_SHA_reference"]!=str(REFERENCE.relative_to(ROOT))
        or p["source_input_total_count"]!=72
        or tuple(p["previously_fixed_ids"])!=IDS
        or tuple(gproto["previously_fixed_nine_ids"])!=IDS
        or tuple(rproto["previously_fixed_ids"])!=IDS
        or tuple(r["fixed_mock_ids"])!=IDS
        or tuple(p["caps"])!=CAPS or tuple(p["tracers"])!=TRACERS
        or tuple(p["roles"])!=ROLES
        or p["fixed_candidate_redshift"]!=origproto["exact_fixed_highz_bin"]
        or p["fixed_sampling"]["data_per_tracer_cap_id"]!=origproto["pilot_counts_per_cap_tracer"]["mock_galaxy_D"]
        or p["fixed_sampling"]["random_per_tracer_cap_id"]!=origproto["pilot_counts_per_cap_tracer"]["mock_random_R"]
        or p["fixed_sampling"]["seed_root"]!=origproto["sampling"]["seed_root"]
        or p["fixed_geometry"]["fiducial"]!=PRIMARY_GEOMETRY
        or p["fixed_geometry"]["line_of_sight"]!=origproto["geometry"]["line_of_sight"]
        or p["fixed_geometry"]["theta_min_deg"]!=origproto["geometry"]["theta_min_deg"]
        or p["fixed_geometry"]["separation_edges_mpc_h"]!=S_EDGES.tolist()
        or origproto["geometry"]["s_edges_mpc_h"]!=S_EDGES.tolist()
        or p["fixed_geometry"]["mu_edges"]!="np.linspace(-1 - 1e-7, 1 + 1e-7, 25), 24 signed bins, physical range clipped only if later independently preregistered multipole projection"
        or p["fixed_weights"]!=origproto["weights"]["formula"]
        or p["fixed_numerical_zero_systot_tol"]!=SYSTOT_NUMERICAL_ZERO_TOL
        or g["status"]!="EZMOCK_PREDECLARED_NINE_GALAXY_RAW_SHA_INVENTORY_ONLY"
        or r["status"]!="EZMOCK_PREDECLARED_NINE_MATCHED_RANDOM_FULL_SHA_REVERIFIED_ONLY"
        or g["errors"]!=[] or r["errors"]!=[]
        or r["complete_random_source_count"]!=36
        or r["no_mock_galaxy_random_or_observed_rows_read"] is not True
        or r["observed_odd_data_vector_read"] is not False
        or g["mock_galaxy_or_random_rows_read"] is not False
        or g["observed_odd_data_vector_read"] is not False
        or g["protocol_sha256"]!=gman["source_protocol_sha256"]
        or r["protocol_sha256"]!=rman["source_protocol_sha256"]
        or gman["uploaded_exact_sha256"]!=p["nine_galaxy_input_source_report_sha256"]
        or rman["uploaded_exact_SHA256"]!=p["nine_random_input_source_report_sha256"]
        or gman["uploaded_exact_bytes"]!=G_REPORT.stat().st_size
        or rman["uploaded_exact_bytes"]!=R_REPORT.stat().st_size
        or ref["source_json_sha256"]!=r["prior_20260924_source_json_sha256"]
        or old["status"]!="EZMOCK0001_MATCHED_GALAXY_CROSS_LS_PILOT_ALGEBRA_ONLY"
        or old["errors"]!=[]
        or old["observed_odd_data_vector_read"] is not False
        or p["observed_odd_data_vector_read"] is not False
        or p["observed_galaxy_rows_read"] is not False
        or p["new_science_selection_applied"] is not False
        or p["inference_or_detection_authorized"] is not False
        or rawproto["mock_realization_id"]!=1
        or old["pilot_mock_id"]!=1
    ):
        raise ValueError("Original mock0001 code/72 source protocol or blinded source report changed")
    check_maps(ref,g,r,gman,rman)
    return g,r,gman,rman,ref,old,gproto,rproto,rawproto,origproto


def resolve_72_sources(p, g,r,gman,rman,ref,gproto,rproto,rawproto):
    """Full 72-file complete compressed-byte SHA; do not read ANY FITS yet."""
    paths={}
    total=0
    for mid,cap,tracer,role in ALL:
        k=key(mid,cap,tracer)
        rel=mock_path(tracer,cap,role,mid)
        if role=="dat":
            if mid==1:
                source=ROOT/rawproto["local_quarantine_dir"]/rel
                row=g["existing_frozen_0001_sources"][k]
                expected=row["full_compressed_sha256_reverified"]
            else:
                source=ROOT/gproto["local_quarantine_dir"]/rel
                row=g["new_galaxy_sources"][k]
                expected=row["full_compressed_sha256_first_seen"]
                if not row["local_path"].endswith("/"+gproto["local_quarantine_dir"]+"/"+rel):
                    raise ValueError("New source local path suffix changed: "+k)
            capbytes=gproto["max_compressed_bytes_each_additional_source"]
        else:
            row=r["verified_random_sources"][k]
            source=Path(row["local_verified_source_path"])
            expected=ref["full_random_gzip_sha256_by_id_cap_tracer"][k]
            capbytes=rproto["max_compressed_bytes_each"]
            from audit_eboss_dr16_nine_ezmock_matched_random_sha import cache_candidates
            approved,_=cache_candidates(rproto,mid,cap,tracer)
            if source.resolve() not in {x.resolve() for x in approved}:
                raise ValueError("Random source path outside fixed approved caches: "+k)
        if (
            not source.is_file() or source.is_symlink()
            or not source.resolve().is_relative_to(ROOT.resolve())
            or (role=="dat" and not str(source.resolve()).endswith("/"+rel))
            or (role=="ran" and not str(source.resolve()).endswith("/"+Path(rel).name))
        ):
            raise ValueError("Missing, external or symlink 72-source input: "+k+"/"+role)
        actual,size=digest(source,capbytes)
        if actual!=expected or size!=row["compressed_bytes"]:
            raise ValueError("Full 72-source compressed SHA/size mismatch: "+k+"/"+role)
        if role=="dat":
            frozen=(gman["previously_frozen_0001_sha"][k] if mid==1
                    else gman["additional_eight_ids_full_compressed_sha"][k])
            oldsha=(frozen["full_compressed_sha256_reverified"] if mid==1
                    else frozen["full_compressed_sha256"])
            if actual!=oldsha:
                raise ValueError("Frozen original galaxy source manifest mismatch: "+k)
        else:
            if actual!=rman["verified_random_sources_by_id_cap_tracer"][k]["previously_pinned_sha256"]:
                raise ValueError("Original 2026-09-24 random source reference mismatch: "+k)
        paths[k+"/"+role]=source
        total+=size
        print("A02_REHASH_72_INPUT_SHA_OK",k,role,size,actual,flush=True)
    if len(paths)!=72 or total!=gman["total_compressed_bytes"]+rman["total_full_compressed_random_bytes"]:
        raise ValueError("All 72 separate full source SHA/size references are not valid")
    print("A02_ALL_72_COMPLETE_GZIP_SHA_OK",len(paths),total,flush=True)
    return paths,total


def original_0001_replay(cap,metadata,fmeta,rmeta,cs,original):
    case=next(c for c in original["cases"] if c["cap"]==cap)
    for tag,item in metadata.items():
        expect=case["input_sample_diagnostics"][tag]
        for attr in ("seed","input_FITS_rows",
                     "eligible_highz_rows_after_fixed_weight_gate",
                     "selected_rows","selected_array_SHA256"):
            if item[attr]!=expect[attr]:
                raise ValueError("Frozen original 0001 mock sample changed: "+cap+"/"+tag+"/"+attr)
    for side,new in (("forward_pair_terms",fmeta),("reverse_pair_terms",rmeta)):
        for name in LABELS:
            got=new[name]
            old=case["cross_ls"][side][name]
            for field in ("accepted_pairs","weighted_histogram_SHA256"):
                if got[field]!=old[field]:
                    raise ValueError("Frozen original 0001 true galaxy weighted pair changed: "+cap+"/"+side+"/"+name+"/"+field)
            if not np.isclose(got["independently_normalized_pair_weight"],
                              old["independently_normalized_pair_weight"],
                              rtol=1e-13,atol=0):
                raise ValueError("Frozen original 0001 pair independent weight norm changed")
    if (cs["RR_supported_s_mu_cells"]!=case["cross_ls"]["RR_supported_s_mu_cells"]
        or cs["RR_total_s_mu_cells"]!=144
        or cs["cross_xi_reverse_max_abs_residual"]>=1e-8):
        raise ValueError("Original 0001 fully supported cross-LS reflection failed")


def one_cap(d_l,d_e,r_l,r_e):
    distance=lambda z:comoving_mpc_over_h(z,PRIMARY_GEOMETRY)
    fwd,fnorm,fmeta=oriented_pair_terms(d_l,d_e,r_l,r_e,distance=distance)
    rev,rnorm,rmeta=oriented_pair_terms(d_e,d_l,r_e,r_l,distance=distance)
    mirrors={}
    for term,revterm in MAPPING.items():
        mirror=mirrored_closure(
            fwd[term],rev[revterm],
            {"pair_normalization":fnorm[term],
             "accepted_pairs":fmeta[term]["accepted_pairs"]},
            {"pair_normalization":rnorm[revterm],
             "accepted_pairs":rmeta[revterm]["accepted_pairs"]},MU_EDGES)
        if mirror["closure_passed"] is not True:
            raise ValueError("True mock-galaxy weighted pair reversal failed: "+term)
        mirrors[term]=mirror
    xi,support=cross_landy_szalay(fwd,fnorm)
    xi_rev,support_rev=cross_landy_szalay(rev,rnorm)
    if not np.array_equal(support,support_rev[:,::-1]):
        raise ValueError("Signed-mu cross-LS RR support differs under tracer reversal")
    if not np.any(support):
        raise ValueError("No RR-supported cell for this fixed cap/ID; preserve case as failure")
    if not np.isfinite(xi[support]).all() or not np.isfinite(xi_rev[support_rev]).all():
        raise ValueError("Nonfinite normalized cross xi on RR-supported cells")
    err=float(np.max(np.abs(xi[support]-xi_rev[:,::-1][support])))
    if not np.isfinite(err) or err>=1e-8:
        raise ValueError("True mock-galaxy cross xi tracer reverse mismatch")
    odd_rr={}
    rr,rrrev=fwd["R1R2"],rev["R1R2"]
    for ell in (1,3):
        if ell==1:
            integ=(MU_EDGES[1:]**2-MU_EDGES[:-1]**2)/2.
        else:
            integ=((5./8.)*(MU_EDGES[1:]**4-MU_EDGES[:-1]**4)
                   -(3./4.)*(MU_EDGES[1:]**2-MU_EDGES[:-1]**2))
        avg=integ/np.diff(MU_EDGES)
        a=float((rr@avg).sum()/rr.sum())
        b=float((rrrev@avg).sum()/rrrev.sum())
        residue=abs(a+b)
        if not np.isfinite(residue) or residue>=1e-10:
            raise ValueError("Mock random RR raw odd parity under tracer reversal failed")
        odd_rr[str(ell)]={"raw_RR_odd_reverse_sum_abs":residue}
    missing=np.argwhere(~support).tolist()
    out={
        "status":"oriented_pair_code_transport_closed_on_supported_cells",
        "RR_supported_s_mu_cells":int(np.count_nonzero(support)),
        "RR_total_s_mu_cells":int(support.size),
        "RR_missing_s_mu_indices":missing,
        "full_RR_support":len(missing)==0,
        "forward_pair_terms":fmeta,
        "reverse_pair_terms":rmeta,
        "four_term_reverse_pair_closure":mirrors,
        "cross_xi_reverse_max_abs_residual":err,
        "RR_odd_reverse_parity":odd_rr,
        "forward_xi_grid_SHA256":hashlib.sha256(np.ascontiguousarray(xi).tobytes()).hexdigest(),
        "physical_odd_multipoles_computed":False,
        "physical_mask_or_pair_window_certified":False,
    }
    return out


def source_header_rows(path, original_case, tracer,role):
    if original_case is not None:
        return int(original_case["input_sample_diagnostics"][tracer+"_"+role]["input_FITS_rows"])
    # Only after ALL 72 full compressed source SHA checks. This is
    # source-internal NAXIS2, NOT a prior independent FITS-header count.
    rows=int(fits.getheader(path,ext=1)["NAXIS2"])
    if rows<=0 or rows>4000000:
        raise ValueError("Authenticated FITS source-internal NAXIS2 outside frozen safety bound")
    return rows


def run(p):
    g,r,gman,rman,ref,original,gproto,rproto,rawproto,origproto=preflight(p,require_local=True)
    paths,total=resolve_72_sources(p,g,r,gman,rman,ref,gproto,rproto,rawproto)
    dest=ROOT/p["local_summary_report"]
    protocol_sha=hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()
    if dest.exists():
        out=load(dest)
        if (
            out.get("protocol_sha256")!=protocol_sha
            or out.get("status") not in (STOP,PASS)
            or out.get("original_0001_cross_ls_SHA256")!=p["original_0001_cross_ls_report_SHA256"]
            or out.get("observed_odd_data_vector_read") is not False
            or out.get("new_science_selection_applied") is not False
            or not isinstance(out.get("cases"),dict)
        ):
            raise ValueError("Prior A-02 code-transport checkpoint mismatches frozen protocol")
        if out["status"]==PASS:
            print("A02_COMPLETE_PREVIOUSLY_CHECKPOINTED_REHASH_72_VALID",flush=True)
            return out
        out["errors"]=[]
    else:
        out={
            "status":STOP,
            "protocol_sha256":protocol_sha,
            "source_galaxy_archive_SHA256":p["nine_galaxy_input_source_report_sha256"],
            "source_random_archive_SHA256":p["nine_random_input_source_report_sha256"],
            "original_0001_cross_ls_SHA256":p["original_0001_cross_ls_report_SHA256"],
            "all_72_full_compressed_gzip_SHA_checked_before_any_FITS_rows":True,
            "total_verified_72_full_compressed_bytes":total,
            "mock_ids":list(IDS),"caps":list(CAPS),
            "fixed_D_R_counts_per_tracer_cap_id":p["fixed_sampling"],
            "cases":{},
            "errors":[],
            "mock_galaxy_and_matched_random_rows_read":False,
            "observed_galaxy_rows_read":False,
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
            "mock_covariance_computed":False,
            "physical_LRG_ELG_pair_window_certified":False,
            "not_a_physical_odd_null_or_detection":True,
        }
    atomic(dest,out)
    nfail=0
    for mid in IDS:
        for cap in CAPS:
            casekey=f"{mid:04d}/{cap}"
            existing=out["cases"].get(casekey)
            if existing is not None and existing.get("status")=="complete":
                if mid==1 and existing.get("exact_original_0001_sample_and_pair_SHA_replay") is not True:
                    raise ValueError("Original 0001 replay has not been certified in previous A-02 checkpoint")
                print("A02_REUSED_IMMUTABLE_SUCCESS_CASE",casekey,flush=True)
                continue
            metadata,cat={},{}
            try:
                if mid!=1 and any(out["cases"].get(f"0001/{c}",{}).get("status")!="complete" for c in CAPS):
                    raise ValueError("Original 0001 source/sample/pair replay must close in both caps before new mock cases")
                prior=next((z for z in original["cases"] if z["cap"]==cap),None) if mid==1 else None
                for tracer in TRACERS:
                    k=key(mid,cap,tracer)
                    for role in ROLES:
                        path=paths[k+"/"+role]
                        rows=source_header_rows(path,prior,tracer,role)
                        values,info=sample_catalogue(
                            path,expected_rows=rows,cap=cap,tracer=tracer,
                            role=role,expected_highz=None,p=origproto)
                        cat[tracer,role]=values
                        metadata[tracer+"_"+role]=info
                        print("A02_FIXED_MOCK_INPUT",casekey,tracer,role,
                              info["eligible_highz_rows_after_fixed_weight_gate"],
                              info["selected_rows"],flush=True)
                cls=one_cap(cat["eBOSS_LRG","dat"],cat["eBOSS_ELG","dat"],
                            cat["eBOSS_LRG","ran"],cat["eBOSS_ELG","ran"])
                if mid==1:
                    original_0001_replay(
                        cap,metadata,cls["forward_pair_terms"],
                        cls["reverse_pair_terms"],cls,original)
                current={
                    "status":"complete","id":mid,"cap":cap,
                    "input_sample_diagnostics":metadata,
                    "cross_ls":cls,
                    "exact_original_0001_sample_and_pair_SHA_replay":mid==1,
                    "physical_odd_multipoles_computed":False,
                    "not_an_18D_mock_covariance":True,
                }
                print("A02_MATCHED_NINE_MOCK_CROSS_LS_CAP",casekey,
                      cls["RR_supported_s_mu_cells"],
                      cls["cross_xi_reverse_max_abs_residual"],flush=True)
            except Exception as exc:
                current={"status":"incomplete_stop","id":mid,"cap":cap,
                         "errors":[str(exc)],
                         "inputs_selected_before_failure":list(metadata),
                         "no_post_result_reselection":True}
                out["errors"].append(casekey+": "+str(exc))
                nfail+=1
                print("A02_MOCK_GALAXY_CASE_FAILURE",casekey,str(exc),flush=True)
            if existing is not None and existing.get("status")=="complete" and existing!=current:
                raise ValueError("Previously successful case result changed unexpectedly: "+casekey)
            out["cases"][casekey]=current
            out["mock_galaxy_and_matched_random_rows_read"]=True
            atomic(dest,out)
            if mid==1 and current["status"]!="complete":
                raise ValueError("Original 0001 replay fails: stop before any new mock ID")
    if len(out["cases"])!=18 or any(x.get("status")!="complete" for x in out["cases"].values()) or nfail:
        out["status"]=STOP
    else:
        out["status"]=PASS
    out["completed_cases"]=sum(x["status"]=="complete" for x in out["cases"].values())
    out["failed_cases"]=sum(x["status"]!="complete" for x in out["cases"].values())
    out["errors"]=list(dict.fromkeys(out["errors"]))
    atomic(dest,out)
    return out


def self_test(p):
    g,r,gm,rm,ref,old,gp,rp,rawp,original=preflight(p,require_local=False)
    if (len(g["new_galaxy_sources"])!=32 or len(r["verified_random_sources"])!=36
            or len(old["cases"])!=2 or tuple(p["previously_fixed_ids"])!=IDS):
        raise AssertionError("72 frozen original galaxy/random source identity contract failed")
    ref_copy=json.loads(json.dumps(ref))
    ref_copy["full_random_gzip_sha256_by_id_cap_tracer"]["0125/NGC/eBOSS_LRG"]="0"*64
    try:
        check_maps(ref_copy,g,r,gm,rm)
    except ValueError:
        pass
    else:
        raise AssertionError("Tampered prior pinned random SHA reference was accepted")
    x=synthetic_cross_ls(.05)
    if x["supported_s_mu_cells"]<=0 or not x["forward_reverse_xi_mirror_passed"]:
        raise AssertionError("Existing independently synthetic oriented cross-LS closure failed")
    print("EBOSS_NINE_MOCK_GALAXY_CROSS_LS_CODE_TRANSPORT_SYNTHETIC_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    p=load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    outpath=ROOT/p["local_summary_report"]
    try:
        out=run(p)
    except Exception as exc:
        if outpath.exists():
            out=load(outpath)
            out["status"]=STOP
            out["errors"]=list(dict.fromkeys(out.get("errors",[])+[str(exc)]))
        else:
            out={
                "status":STOP,"errors":[str(exc)],
                "all_72_full_compressed_gzip_SHA_checked_before_any_FITS_rows":False,
                "observed_galaxy_rows_read":False,
                "observed_odd_data_vector_read":False,
                "new_science_selection_applied":False,
                "mock_covariance_computed":False,
            }
        atomic(outpath,out)
    print("EZMOCK_NINE_GALAXY_CROSS_LS_CODE_TRANSPORT",out["status"],flush=True)
    print("REPORT",outpath,flush=True)
    if out.get("errors"):
        print("ERRORS",*out["errors"],sep="\n",flush=True)
        return 2
    print("COMPLETED_CASES",out.get("completed_cases",len(out.get("cases",{}))),flush=True)
    print("OBSERVED_ODD_DATA_READ",out["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
