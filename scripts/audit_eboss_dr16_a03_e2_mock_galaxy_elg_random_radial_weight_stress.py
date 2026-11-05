#!/usr/bin/env python3
"""A-03E2: prospective mock-only ELG RANDOM radial-weight ±5% LS stress.

Before ANY mock FITS header/row: replay archived immutable parent checks
and rehash ALL 72 complete same-ID mock galaxy/random gzip sources.
All 18 preselected ID/cap cases replay original A-02 sample/weighted-pair
fingerprints before modifying weights of the selected ELG RANDOMS only.
This is a descriptive model-mismatch probe; NOT a physical systematic
calibration, mock-null covariance, wake detection or observed odd analysis.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import (
    sample_catalogue, oriented_pair_terms,
)
from audit_eboss_dr16_ezmock0001_galaxy_odd_projection import project
from audit_eboss_dr16_rr_pair_closure import mirrored_closure
from check_eboss_cross_ls_synthetic import LABELS, cross_landy_szalay
from eboss_dr16_fiducial import PRIMARY_GEOMETRY, comoving_mpc_over_h

ROOT=A02.ROOT
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_protocol_2026-09-26.json"
PARENT_REPORT=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
DECISION=ROOT/"source_data/eboss_dr16_a03_empirical_only_user_decision_2026-09-26.json"
E01=ROOT/"source_data/eboss_dr16_a03_empirical_rr_window_injection_report_2026-09-26.json"
PROJECT_SOURCE=ROOT/"scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py"
PROJECT_SOURCE_BLOB="5d1831570639abd4610eec76044ac0c6f04912a5"
PASS="A03E2_NINE_MOCK_GALAXY_ELG_RANDOM_RADIAL_WEIGHT_STRESS_DESCRIPTIVE_ONLY"
STOP="A03E2_NINE_MOCK_GALAXY_ELG_RANDOM_RADIAL_WEIGHT_STRESS_INCOMPLETE_STOP"
SCENARIOS=("baseline_original","plus_5pct_ELG_R_z_ramp","minus_5pct_ELG_R_z_ramp")
ELLS=(0,1,2,3)
MAX_DELTA=0.05


def blob(path):
    raw=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def guard_protocol(p):
    parent=A02.load(A02.PROTOCOL)
    decision=A02.load(DECISION)
    if (
        p["registered_before"]!="First new A-03E2 mock FITS row read or A-03E2 mock-galaxy stress output"
        or p["user_decision"]!="EMPIRICAL_ONLY_DO_NOT_CONTACT_AUTHORS"
        or decision["author_contact_approved"] is not False
        or decision["email_sent"] is not False
        or p["author_contact_permitted"] is not False
        or p["observed_galaxy_rows_read"] is not False
        or p["observed_odd_data_vector_read"] is not False
        or p["main_branch_mutation_permitted"] is not False
        or p["mock_null_significance_computed"] is not False
        or p["physical_selection_systematics_certified"] is not False
        or p["physical_LRG_ELG_pair_window_certified"] is not False
        or p["ids"]!=list(A02.IDS) or p["caps"]!=list(A02.CAPS)
        or p["tracers"]!=list(A02.TRACERS)
        or p["roles"]!=list(A02.ROLES)
        or p["n_cases"]!=18
        or p["scenarios"]!=list(SCENARIOS)
        or p["diagnostics"]["ell"]!=list(ELLS)
        or p["stress_input"]["max_absolute_relative_weight_shift"]!=MAX_DELTA
        or p["stress_input"]["role"]!="eBOSS_ELG/ran only"
        or parent["fixed_candidate_redshift"]!=[0.9,1.0]
        or p["geometry"]!=
           "Reuse A02 exact s=[20,40,60,80,100,120,140] Mpc/h, 24 signed mu edges np.linspace(-1-1e-7,1+1e-7,25), midpoint LOS, theta_min_deg=0.05, published_multitracer fiducial, four separately normalized weighted LS terms."
        or blob(A02.PROTOCOL)!=p["parent_A02_protocol_git_blob_sha1"]
        or blob(ROOT/p["parent_A02_implementation"])!=p["parent_A02_implementation_git_blob_sha1"]
        or blob(ROOT/p["parent_0001_sampling_implementation"])!=p["parent_0001_sampling_implementation_git_blob_sha1"]
        or blob(PROJECT_SOURCE)!=PROJECT_SOURCE_BLOB
        or file_sha(E01)!=p["A03E0E1_original_report_sha256"]
    ):
        raise ValueError("Frozen E2 parent source, user decision, geometry or no-observed contract changed")
    raw=PARENT_REPORT.read_bytes()
    if (len(raw)!=p["parent_A02_report_bytes"]
        or hashlib.sha256(raw).hexdigest()!=p["parent_A02_report_full_sha256"]
        or blob(PARENT_REPORT)!=p["parent_A02_report_git_blob_sha1"]):
        raise ValueError("Exact 151892-byte original A-02 report source changed")
    old=json.loads(raw)
    expected=[f"{mid:04d}/{cap}" for mid in A02.IDS for cap in A02.CAPS]
    if (
        old["status"]!=A02.PASS or old["errors"]!=[]
        or old["completed_cases"]!=18 or old["failed_cases"]!=0
        or list(old["cases"])!=expected
        or old["all_72_full_compressed_gzip_SHA_checked_before_any_FITS_rows"] is not True
        or old["observed_odd_data_vector_read"] is not False
        or any(c["status"]!="complete"
               or c["cross_ls"]["RR_supported_s_mu_cells"]!=144
               or c["physical_odd_multipoles_computed"] is not False
               for c in old["cases"].values())
        or parent["source_input_total_count"]!=72
    ):
        raise ValueError("Incomplete original 18/18 A-02 mock code transport; E2 is blocked")
    return parent,old


def verify_original_sample(metadata, archived):
    """Whole original deterministic sample metadata, not just selected count."""
    if metadata!=archived:
        raise ValueError("Original A-02 selected galaxy/random sample or ELG chunk/seed/weight fingerprint differs")


def verify_original_pairs(forward,reverse,old):
    for label,actual in (("forward_pair_terms",forward),
                         ("reverse_pair_terms",reverse)):
        archived=old["cross_ls"][label]
        for term in LABELS:
            new,orig=actual[term],archived[term]
            for k in ("accepted_pairs","weighted_histogram_SHA256"):
                if new[k]!=orig[k]:
                    raise ValueError("Original A-02 pair histogram fingerprint changed: "+label+"/"+term+"/"+k)
            for k in ("independently_normalized_pair_weight",
                      "total_weighted_pairs_in_fixed_s_mu_bins",
                      "normalization_relative_residual"):
                if not np.isclose(new[k],orig[k],rtol=1e-12,atol=1e-14):
                    raise ValueError("Original A-02 pair weighted norm changed: "+label+"/"+term+"/"+k)


def reweight_elg_random(cat, amplitude):
    if amplitude not in (0.,MAX_DELTA,-MAX_DELTA):
        raise ValueError("Stress amplitude outside independently frozen baseline and ±5%")
    ra,dec,z,weights=cat
    if not all(a.shape==(1200,) and np.isfinite(a).all() for a in cat):
        raise ValueError("ELG RANDOM original 1200-row source sample changed")
    if not np.all((z>=0.9)&(z<1.0)) or np.any(weights<=0):
        raise ValueError("ELG RANDOM highz/positive weights changed")
    u=2.0*(z-0.9)/0.1-1.0
    fac=1.0+amplitude*u
    if (not np.isfinite(fac).all()
        or np.any(fac<1.0-MAX_DELTA-1e-12)
        or np.any(fac>1.0+MAX_DELTA+1e-12)
        or np.any(fac<=0)):
        raise ValueError("Modified ELG RANDOM factor outside [0.95,1.05]")
    changed=weights*fac
    if np.any(changed<=0) or not np.isfinite(changed).all():
        raise ValueError("Modified ELG RANDOM weights became invalid")
    out=(ra,dec,z,np.asarray(changed,dtype="f8"))
    for x,y in zip(out[:3],cat[:3]):
        if not np.array_equal(x,y):
            raise ValueError("ELG RANDOM positions/redshifts altered in a weight-only stress")
    return out,{
        "factor_min":float(np.min(fac)),
        "factor_max":float(np.max(fac)),
        "max_absolute_relative_factor_change":float(np.max(np.abs(fac-1.0))),
        "selected_original_ELG_R_weight_SHA256":hashlib.sha256(
            np.ascontiguousarray(weights).tobytes()).hexdigest(),
        "modified_selected_ELG_R_weight_SHA256":hashlib.sha256(
            np.ascontiguousarray(out[3]).tobytes()).hexdigest(),
        "selected_ELG_R_rows":len(z),
        "ELG_RANDOM_only_weight_changed":True,
    }


def compute_scenario(cats, elg_random):
    dist=lambda z:comoving_mpc_over_h(z,PRIMARY_GEOMETRY)
    dl,de,rl=(cats["eBOSS_LRG","dat"],
              cats["eBOSS_ELG","dat"],
              cats["eBOSS_LRG","ran"])
    f,fn,fmeta=oriented_pair_terms(dl,de,rl,elg_random,distance=dist)
    rev,rn,rmeta=oriented_pair_terms(de,dl,elg_random,rl,distance=dist)
    mirrored={}
    for term,other in A02.MAPPING.items():
        c=mirrored_closure(f[term],rev[other],
            {"pair_normalization":fn[term],
             "accepted_pairs":fmeta[term]["accepted_pairs"]},
            {"pair_normalization":rn[other],
             "accepted_pairs":rmeta[other]["accepted_pairs"]},
            A02.MU_EDGES)
        if c["closure_passed"] is not True:
            raise ValueError("Independent forward/reverse pair mirror failed: "+term)
        mirrored[term]=c
    xi,support=cross_landy_szalay(f,fn)
    rxi,rsupport=cross_landy_szalay(rev,rn)
    if (
        xi.shape!=(6,24) or not np.all(support)
        or not np.array_equal(support,rsupport[:,::-1])
        or not np.isfinite(xi).all() or not np.isfinite(rxi).all()
    ):
        raise ValueError("Stress scenario missing RR support, changed signed-mu support or nonfinite xi; never zero-fill")
    resid=float(np.max(np.abs(xi-rxi[:,::-1])))
    if not np.isfinite(resid) or resid>=1e-8:
        raise ValueError("Stress scenario xi forward/reverse closure failed")
    proj=project(xi,support,A02.MU_EDGES,ELLS)
    record={
        "status":"complete_synthetic_random_weight_stress_scenario",
        "RR_supported_cells":int(support.sum()),
        "max_abs_xi_tracer_reverse_mirror":resid,
        "forward_xi_grid_sha256":hashlib.sha256(np.ascontiguousarray(xi).tobytes()).hexdigest(),
        "forward_pair_terms":fmeta,
        "reverse_pair_terms":rmeta,
        "all_four_weighted_pair_mirror_closure_passed":all(
            v["closure_passed"] for v in mirrored.values()),
        "mock_only_descriptive_xi_multipoles":proj,
        "observed_odd_vector_read":False,
    }
    return record


def check_unchanged_terms(stress,base):
    expected={"forward_pair_terms":("D1D2","R1D2"),
              "reverse_pair_terms":("D1D2","D1R2")}
    changed={"forward_pair_terms":("D1R2","R1R2"),
             "reverse_pair_terms":("R1D2","R1R2")}
    for side,terms in expected.items():
        for term in terms:
            if stress[side][term]!=base[side][term]:
                raise ValueError("ELG RANDOM-only stress modified an unrelated DATA/LRG pair term: "+side+"/"+term)
    for side,terms in changed.items():
        for term in terms:
            if stress[side][term]["accepted_pairs"]!=base[side][term]["accepted_pairs"]:
                raise ValueError("Changing ELG RANDOM weights unexpectedly changed accepted geometric pair count")
    return True


def finish_case(parent,cats,metadata):
    verify_original_sample(metadata,parent["input_sample_diagnostics"])
    baseline,bi=reweight_elg_random(cats["eBOSS_ELG","ran"],0.)
    if not np.array_equal(baseline[3],cats["eBOSS_ELG","ran"][3]):
        raise ValueError("Baseline reweight changed original ELG RANDOM weight bytes")
    b=compute_scenario(cats,baseline)
    verify_original_pairs(b["forward_pair_terms"],b["reverse_pair_terms"],parent)
    if b["forward_xi_grid_sha256"]!=parent["cross_ls"]["forward_xi_grid_SHA256"]:
        raise ValueError("Original A-02 forward xi grid fingerprint changed before stress")
    if b["RR_supported_cells"]!=parent["cross_ls"]["RR_supported_s_mu_cells"]:
        raise ValueError("Original A-02 RR support changed before stress")
    result={
        "status":"complete","id":parent["id"],"cap":parent["cap"],
        "original_A02_sample_and_all_eight_pair_SHA_exact_replay":True,
        "original_A02_forward_xi_SHA_exact_replay":True,
        "original_selected_array_SHA256":{
            k:v["selected_array_SHA256"] for k,v in metadata.items()},
        "original_exact_ELG_chunk_counts":{
            k:v["ELG_exact_chunk_diagnostic"] for k,v in metadata.items()
            if k.startswith("eBOSS_ELG_")},
        "scenarios":{"baseline_original":{**b,"ELG_R_weight_stress":bi}},
        "diagnostic_deltas_from_original":{},
        "no_observed_rows_or_odd_read":True,
        "not_a_physical_null_or_selection_bias_calibration":True,
    }
    for name,amp in (("plus_5pct_ELG_R_z_ramp",MAX_DELTA),
                     ("minus_5pct_ELG_R_z_ramp",-MAX_DELTA)):
        reweighted,diag=reweight_elg_random(cats["eBOSS_ELG","ran"],amp)
        if all(np.array_equal(reweighted[3],x[3]) for x in (cats["eBOSS_ELG","ran"],)):
            raise ValueError("Nontrivial ±5% radial test left all selected random weights unchanged")
        rec=compute_scenario(cats,reweighted)
        check_unchanged_terms(rec,b)
        result["scenarios"][name]={**rec,"ELG_R_weight_stress":diag}
        dif={}
        for ell in ELLS:
            before=np.asarray(b["mock_only_descriptive_xi_multipoles"][str(ell)]["values_by_fixed_s_bin"])
            after=np.asarray(rec["mock_only_descriptive_xi_multipoles"][str(ell)]["values_by_fixed_s_bin"])
            delta=after-before
            if not np.isfinite(delta).all():
                raise ValueError("Nonfinite mock-only descriptive stress multipole shift")
            dif[str(ell)]={
                "delta_by_fixed_s_bin":delta.tolist(),
                "max_abs_delta_across_s_bins":float(np.max(np.abs(delta))),
            }
        result["diagnostic_deltas_from_original"][name]=dif
    for ell in ELLS:
        p=np.asarray(result["diagnostic_deltas_from_original"]["plus_5pct_ELG_R_z_ramp"][str(ell)]["delta_by_fixed_s_bin"])
        m=np.asarray(result["diagnostic_deltas_from_original"]["minus_5pct_ELG_R_z_ramp"][str(ell)]["delta_by_fixed_s_bin"])
        result.setdefault("nonlinear_stress_symmetric_sum_max_abs_by_ell",{})[str(ell)]=float(np.max(np.abs(p+m)))
    return result


def synthetic_negative_controls(p):
    """Pure synthetic controls; no source gzip/FITS/observed galaxy rows."""
    z=np.asarray([0.90,0.925,0.95,0.975,0.999999]+[0.951]*1195,dtype="f8")
    ra=np.arange(1200,dtype="f8")%360
    cat=(ra,np.zeros(1200),z,np.ones(1200,dtype="f8"))
    unchanged,zero=reweight_elg_random(cat,0.)
    plus,diag=reweight_elg_random(cat,MAX_DELTA)
    minus,_=reweight_elg_random(cat,-MAX_DELTA)
    if (not np.array_equal(unchanged[3],cat[3])
        or not np.array_equal(cat[3],np.ones(1200))
        or not np.array_equal(cat[0],plus[0])
        or not np.array_equal(cat[2],minus[2])
        or abs(float(plus[3][0])-0.95)>1e-14
        or abs(float(minus[3][0])-1.05)>1e-14
        or diag["max_absolute_relative_factor_change"]>MAX_DELTA+1e-12):
        raise AssertionError("ELG RANDOM-only synthetic weight factor control failed")
    for invalid in (0.051,-0.051,0.1):
        try:
            reweight_elg_random(cat,invalid)
        except ValueError:
            pass
        else:
            raise AssertionError("Changed unregistered ±5% scenario accepted")
    rr=np.ones((6,24),dtype="f8")*11
    h={"D1D2":rr*1.025,"D1R2":rr.copy(),
       "R1D2":rr.copy(),"R1R2":rr.copy()}
    norms={term:1. for term in LABELS}
    base,supp=cross_landy_szalay(h,norms)
    if not np.all(supp) or not np.allclose(base,0.025,rtol=0,atol=1e-14):
        raise AssertionError("Synthetic constant nonzero xi source model failed")
    uniform=copy.deepcopy(h)
    uniform["D1R2"]*=1.03
    uniform["R1R2"]*=1.03
    unorm=dict(norms,D1R2=1.03,R1R2=1.03)
    same,_=cross_landy_szalay(uniform,unorm)
    wrong,_=cross_landy_szalay(uniform,norms)
    if (not np.allclose(base,same,rtol=0,atol=1e-14)
        or np.allclose(base,wrong,rtol=0,atol=1e-5)):
        raise AssertionError("Independent cross-LS normalization negative control failed")
    mu=(A02.MU_EDGES[1:]+A02.MU_EDGES[:-1])/2
    mixed=copy.deepcopy(h)
    mixed["D1R2"]*=1+0.04*mu[None,:]
    mixed["R1R2"]*=1+0.03*mu[None,:]
    shifted,sup=cross_landy_szalay(mixed,norms)
    response=project(shifted,sup,A02.MU_EDGES,(1,3))
    if np.max(np.abs(response["1"]["values_by_fixed_s_bin"]))<=1e-4:
        raise AssertionError("Synthetic spatially varying mismatched random did not generate diagnostic odd")
    mirror={A02.MAPPING[k]:v[:,::-1].copy() for k,v in mixed.items()}
    rev_norm={A02.MAPPING[k]:v for k,v in norms.items()}
    rev,rs=cross_landy_szalay(mirror,rev_norm)
    if not np.array_equal(sup,rs[:,::-1]) or not np.allclose(shifted,rev[:,::-1],atol=1e-14):
        raise AssertionError("Signed mu reversed algebraic closure should hold even under wrong random")
    tampered=copy.deepcopy(mixed)
    tampered["R1R2"][0,0]=0.
    tsup=cross_landy_szalay(tampered,norms)[1]
    try:
        project(base,tsup,A02.MU_EDGES,(1,3))
    except ValueError:
        pass
    else:
        raise AssertionError("Unsupported RR cell was projected")
    original={"input_sample_diagnostics":{"eBOSS_ELG_ran":{
        "selected_array_SHA256":"a","selected_rows":1200}}}
    same_meta=copy.deepcopy(original["input_sample_diagnostics"])
    same_meta["eBOSS_ELG_ran"]["selected_array_SHA256"]="b"
    try:
        verify_original_sample(same_meta,original["input_sample_diagnostics"])
    except ValueError:
        pass
    else:
        raise AssertionError("Tampered original deterministic sample SHA accepted")
    f={"D1D2":{"weighted_histogram_SHA256":"abc","accepted_pairs":9,
                "independently_normalized_pair_weight":1.,
                "total_weighted_pairs_in_fixed_s_mu_bins":1.,
                "normalization_relative_residual":0.}}
    prior={"cross_ls":{"forward_pair_terms":copy.deepcopy(f),
                       "reverse_pair_terms":copy.deepcopy(f)}}
    for t in LABELS:
        f.setdefault(t,copy.deepcopy(f["D1D2"]))
        prior["cross_ls"]["forward_pair_terms"].setdefault(t,copy.deepcopy(f["D1D2"]))
        prior["cross_ls"]["reverse_pair_terms"].setdefault(t,copy.deepcopy(f["D1D2"]))
    rev=copy.deepcopy(f)
    f["R1R2"]["weighted_histogram_SHA256"]="tampered"
    try:
        verify_original_pairs(f,rev,prior)
    except ValueError:
        pass
    else:
        raise AssertionError("Tampered original pair-histogram SHA accepted")


def self_test(p):
    parent,old=guard_protocol(p)
    # Archive-only parent source tests. Do NOT require local 72 FITS sources.
    A02.preflight(parent,require_local=False)
    synthetic_negative_controls(p)
    print("A03E2_SOURCE_ONLY_AND_SYNTHETIC_NEGATIVE_CONTROLS_OK",flush=True)
    print("PARENT_A02_COMPLETED_CASES",len(old["cases"]),flush=True)
    print("OBSERVED_ODD_DATA_READ",False,flush=True)


def run(p,parent,old):
    # Stop before ANY mock FITS if local original uploaded source reports differ.
    g,r,gm,rm,ref,original,gp,rp,rawp,original_protocol=A02.preflight(
        parent,require_local=True)
    paths,total=A02.resolve_72_sources(parent,g,r,gm,rm,ref,gp,rp,rawp)
    dest=ROOT/p["local_output"]
    phash=file_sha(PROTOCOL)
    if dest.exists():
        out=A02.load(dest)
        if (out.get("protocol_sha256")!=phash
            or out.get("parent_A02_report_sha256")!=p["parent_A02_report_full_sha256"]
            or out.get("observed_odd_data_vector_read") is not False
            or out.get("observed_galaxy_rows_read") is not False
            or out.get("new_science_selection_applied") is not False
            or out.get("status") not in (STOP,PASS)
            or not isinstance(out.get("cases"),dict)
            or not set(out["cases"]).issubset(set(old["cases"]))):
            raise ValueError("Existing A-03E2 checkpoint does not match frozen protocol")
        if out["status"]==PASS:
            if (len(out["cases"])!=18
                or any(c.get("status")!="complete" for c in out["cases"].values())):
                raise ValueError("Previously completed E2 checkpoint has missing cases")
            print("A03E2_ALREADY_COMPLETED_72_SOURCES_REVERIFIED",flush=True)
            return out
        out.setdefault("previous_errors",[]).extend(out.get("errors",[]))
        out["errors"]=[]
    else:
        out={
            "status":STOP,
            "protocol_sha256":phash,
            "parent_A02_report_sha256":p["parent_A02_report_full_sha256"],
            "parent_A03E0E1_report_sha256":p["A03E0E1_original_report_sha256"],
            "all_72_original_gzip_full_sha256_rechecked_before_any_E2_FITS_row":True,
            "total_original_72_gzip_compressed_bytes":total,
            "fixed_ids":list(A02.IDS),"caps":list(A02.CAPS),
            "scenarios":list(SCENARIOS),
            "cases":{},"errors":[],
            "observed_galaxy_rows_read":False,
            "observed_random_rows_read":False,
            "observed_odd_data_vector_read":False,
            "new_science_selection_applied":False,
            "physical_window_certified":False,
            "18D_covariance_computed":False,
            "not_a_detection_or_exclusion":True,
        }
    A02.atomic(dest,out)
    for mid in A02.IDS:
        for cap in A02.CAPS:
            tag=f"{mid:04d}/{cap}"
            prev=out["cases"].get(tag)
            if prev is not None and prev.get("status")=="complete":
                if prev.get("original_A02_sample_and_all_eight_pair_SHA_exact_replay") is not True:
                    raise ValueError("Previous successful E2 checkpoint is not SHA verified")
                print("A03E2_REUSED_PREVIOUSLY_CHECKPOINTED_CASE",tag,flush=True)
                continue
            if mid!=1 and any(out["cases"].get(f"0001/{c}",{}).get("status")!="complete" for c in A02.CAPS):
                raise ValueError("Both original 0001 E2 cases must succeed before further mock IDs")
            frozen=old["cases"][tag]
            cats,metadata={},{}
            try:
                original_0001=next((c for c in original["cases"] if c["cap"]==cap),None) if mid==1 else None
                for tracer in A02.TRACERS:
                    for role in A02.ROLES:
                        path=paths[A02.key(mid,cap,tracer)+"/"+role]
                        rows=A02.source_header_rows(path,original_0001,tracer,role)
                        values,info=sample_catalogue(
                            path,expected_rows=rows,cap=cap,tracer=tracer,
                            role=role,expected_highz=None,p=original_protocol)
                        cats[tracer,role]=values
                        metadata[tracer+"_"+role]=info
                        if info!=frozen["input_sample_diagnostics"][tracer+"_"+role]:
                            raise ValueError("Original A-02 deterministic sample SHA/metadata changed: "+tracer+"/"+role)
                current=finish_case(frozen,cats,metadata)
                print("A03E2_MOCK_SELECTION_STRESS_CASE",tag,
                    "plus_max_abs_delta_ell1",
                    current["diagnostic_deltas_from_original"]["plus_5pct_ELG_R_z_ramp"]["1"]["max_abs_delta_across_s_bins"],
                    "minus_max_abs_delta_ell3",
                    current["diagnostic_deltas_from_original"]["minus_5pct_ELG_R_z_ramp"]["3"]["max_abs_delta_across_s_bins"],
                    flush=True)
            except Exception as exc:
                current={
                    "status":"incomplete_stop","id":mid,"cap":cap,
                    "errors":[str(exc)],
                    "sample_roles_processed_before_failure":list(metadata),
                    "no_resampling_or_ID_replacement":True,
                }
                out["errors"].append(tag+": "+str(exc))
                print("A03E2_CASE_FAILURE",tag,str(exc),flush=True)
            out["cases"][tag]=current
            A02.atomic(dest,out)
            if mid==1 and current["status"]!="complete":
                raise ValueError("Original 0001 E2 sample/pair replay failed: STOP before other mock IDs")
    out["completed_cases"]=sum(c["status"]=="complete" for c in out["cases"].values())
    out["failed_cases"]=sum(c["status"]!="complete" for c in out["cases"].values())
    if (out["completed_cases"]==18 and out["failed_cases"]==0
        and len(out["cases"])==18 and not out["errors"]):
        out["status"]=PASS
    else:
        out["status"]=STOP
    by_cap={}
    for cap in A02.CAPS:
        if not all(out["cases"].get(f"{mid:04d}/{cap}",{}).get("status")=="complete" for mid in A02.IDS):
            continue
        by_cap[cap]={}
        for scenario in SCENARIOS[1:]:
            by_cap[cap][scenario]={}
            for ell in (1,3):
                vals=np.asarray([
                    out["cases"][f"{mid:04d}/{cap}"]["diagnostic_deltas_from_original"][scenario][str(ell)]["max_abs_delta_across_s_bins"]
                    for mid in A02.IDS],dtype="f8")
                by_cap[cap][scenario][str(ell)]={
                    "all_nine_fixed_ID_max_abs_delta":vals.tolist(),
                    "descriptive_min":float(np.min(vals)),
                    "descriptive_median":float(np.median(vals)),
                    "descriptive_max":float(np.max(vals)),
                }
    out["by_cap_descriptive"]=by_cap
    A02.atomic(dest,out)
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    p=A02.load(PROTOCOL)
    if args.self_test:
        self_test(p)
        return 0
    dest=ROOT/p["local_output"]
    try:
        parent,old=guard_protocol(p)
        out=run(p,parent,old)
    except Exception as exc:
        if dest.exists():
            try:
                out=A02.load(dest)
                out["status"]=STOP
                out["errors"]=list(dict.fromkeys(out.get("errors",[])+[str(exc)]))
            except (OSError,ValueError,KeyError,TypeError):
                out={"status":STOP,"errors":[str(exc)]}
        else:
            out={"status":STOP,"errors":[str(exc)]}
        out["observed_galaxy_rows_read"]=False
        out["observed_random_rows_read"]=False
        out["observed_odd_data_vector_read"]=False
        out["new_science_selection_applied"]=False
        out["physical_window_certified"]=False
        out["18D_covariance_computed"]=False
        out["not_a_detection_or_exclusion"]=True
        A02.atomic(dest,out)
    print("A03E2_MOCK_GALAXY_ELG_RANDOM_RADIAL_WEIGHT_STRESS",out["status"],flush=True)
    print("REPORT",dest,flush=True)
    print("COMPLETED_CASES",out.get("completed_cases",sum(
        c.get("status")=="complete" for c in out.get("cases",{}).values())),flush=True)
    if out.get("errors"):
        print("ERRORS",*out["errors"],sep="\n",flush=True)
        return 2
    print("OBSERVED_ODD_DATA_READ",out["observed_odd_data_vector_read"],flush=True)
    return 0 if out["status"]==PASS else 2


if __name__=="__main__":
    raise SystemExit(main())
