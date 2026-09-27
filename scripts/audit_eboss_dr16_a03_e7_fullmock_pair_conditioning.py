#!/usr/bin/env python3
"""A-03E7: source-only diagnosis of the completed full-galaxy mock0001 LS drift.

POST-RESULT mathematical audit. No independent random realization, FITS access,
physical error tolerance, A03 survey-window certification, covariance, or
observed odd vector. Exact Shapley shares describe algebraic DR/RD/RR
contributions to *this* fixed 4800-to-48000 comparison; they are not
stochastic variance or causal attribution.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/"source_data/eboss_dr16_mock0001_full_galaxy_4800_48000_completed_source_only_audit_2026-09-27.json"
ARCHIVE_GIT_BLOB="aa080ec8a3f5faab1731396bed1b048c7dfa991c"
REPORT_SHA256="187e21e206cdb6ded1108887a2b6670113d50b752e438d28212cfc3d187572da"
PARENT_NGC_SHA256="27f31f62663122bbd59b1d0cfa0ef40b784be7f8332d26bcb432e5d29a147cfe"
TERMS=("D1D2","D1R2","R1D2","R1R2")
CHANGED=("D1R2","R1D2","R1R2")
MU_EDGES=np.linspace(-1.-1.e-7,1.+1.e-7,25,dtype="f8")
OUT_DEFAULT=ROOT/"eboss_workspace/a03_full_eligible_mock/a03_e7_fullmock0001_source_only_conditioning.json"

def require(condition,message):
    if not condition:raise ValueError(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def git_blob(raw):
    return hashlib.sha1(("blob "+str(len(raw))).encode()+b"\x00"+raw).hexdigest()

def canonical(x):
    return sha(json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False).encode())

def original_archive():
    raw=ARCHIVE.read_bytes()
    require(git_blob(raw)==ARCHIVE_GIT_BLOB,
            "Changed original completed-fullmock uploaded-audit Git blob")
    a=json.loads(raw)
    require(a["source_report"]["sha256"]==REPORT_SHA256 and
            a["source_report"]["bytes"]==529771 and
            a["completed_mock_caps"]==["0001/NGC","0001/SGC"] and
            a["decision"]["observed_odd_remains_sealed"] is True and
            a["decision"]["high_48000_convergence_not_certified"] is True and
            a["decision"]["absolute_EV_wake_amplitude_not_calibrated"] is True and
            a["source_only_json_audit"]["weighted_pair_histogram_SHA256"]==32,
            "Original archive provenance/STOP contract changed")
    return a

def clipped_odd_weight(ell):
    low=np.maximum(-1.,MU_EDGES[:-1])
    high=np.minimum(1.,MU_EDGES[1:])
    require(np.all(high>low),"Invalid signed-mu physical cell overlap")
    if ell==1:
        F=lambda z: z*z/2.
    elif ell==3:
        F=lambda z: (5./8.)*z**4-(3./4.)*z*z
    else:
        raise ValueError("Only original odd ell1 and ell3 allowed")
    w=(2*ell+1)/2.*(F(high)-F(low))
    require(np.allclose(w,(-1)**ell*w[::-1],atol=1e-13,rtol=0),
            "Original exact odd projection lost parity")
    return w

WEIGHT={ell:clipped_odd_weight(ell) for ell in (1,3)}

def projected(x,ell):
    require(np.asarray(x).shape==(6,24),
            "Wrong original 6x24 estimator shape")
    return np.asarray(x,dtype="f8")@WEIGHT[ell]

def verified_terms(level):
    require(level["RR_supported_cells"]==144 and
            level["no_observed_galaxies_or_odd_read"] is True and
            set(level["forward"])==set(TERMS),
            "Full 144-cell mock/source-only term gate failed")
    t={}
    for term in TERMS:
        rec=level["forward"][term]
        h=np.ascontiguousarray(np.asarray(rec["histogram_6x24"],dtype="f8"))
        norm=float(rec["pair_normalization"])
        require(h.shape==(6,24) and np.isfinite(h).all() and np.all(h>=0) and
                norm>0 and np.isfinite(norm) and
                sha(h.tobytes())==rec["histogram_SHA256"],
                "Original SHA-pinned weighted pair histogram changed: "+term)
        t[term]=h/norm
    return t

def xi(terms):
    require(np.all(terms["R1R2"]>0),"Zero RR denominator, refuse extrapolation")
    return 1.+(terms["D1D2"]-terms["D1R2"]-terms["R1D2"])/terms["R1R2"]

def shapley(a,b):
    """Order-independent exact three-group change decomposition."""
    require(np.array_equal(a["D1D2"],b["D1D2"]),
            "DD must be original identical mock galaxy pairs")
    changed={term:np.zeros((6,24),dtype="f8") for term in CHANGED}
    for permutation in itertools.permutations(CHANGED):
        mixed=dict(a)
        for term in permutation:
            before=xi(mixed)
            mixed[term]=b[term]
            changed[term]+=(xi(mixed)-before)/6.
    delta=xi(b)-xi(a)
    require(np.allclose(sum(changed.values()),delta,atol=2e-13,rtol=0),
            "Exact Shapley cross-LS telescoping failed")
    return changed,delta

def verify_case(case):
    original=case["levels"]["original_4800"]
    dense=case["levels"]["nested_48000"]
    a=verified_terms(original)
    b=verified_terms(dense)
    for stage,normalized in ((original,a),(dense,b)):
        calc=np.ascontiguousarray(xi(normalized))
        actual=np.ascontiguousarray(np.asarray(stage["xi_6x24"],dtype="f8"))
        require(actual.shape==(6,24) and np.isfinite(actual).all() and
                sha(actual.tobytes())==stage["xi_SHA256"] and
                np.max(np.abs(actual-calc))<5e-13,
                "Original weighted four-term LS xi or SHA changed")
        for ell in (1,3):
            record=np.asarray(
                stage["odd_and_even_multipoles"][str(ell)]["values_by_fixed_s_bin"],
                dtype="f8")
            require(record.shape==(6,) and
                    np.max(np.abs(projected(calc,ell)-record))<5e-13,
                    "Original exact clipped-mu dipole/octupole changed")
    return a,b

def diagnose(cap,case):
    a,b=verify_case(case)
    shares,delta=shapley(a,b)
    reference=xi(b)
    old=xi(a)
    num=np.abs(b["D1D2"])+np.abs(b["D1R2"])+np.abs(b["R1D2"])+np.abs(b["R1R2"])
    signed=b["D1D2"]-b["D1R2"]-b["R1D2"]+b["R1R2"]
    condition=np.divide(num,np.abs(signed),out=np.full((6,24),np.inf),
                        where=np.abs(signed)>0)
    require(np.isfinite(condition).all(),"Nonfinite mock local cancellation ratio")
    result={
        "cap":cap,
        "xi_4800_mean_abs":float(np.abs(old).mean()),
        "xi_48000_mean_abs":float(np.abs(reference).mean()),
        "xi_change_mean_abs":float(np.abs(delta).mean()),
        "xi_change_relative_L1_to_48000":
            float(np.abs(delta).sum()/np.abs(reference).sum()),
        "xi_change_max_abs_cell":float(np.abs(delta).max()),
        "LS_local_cancellation_ratio_median":float(np.median(condition)),
        "LS_local_cancellation_ratio_p90":float(np.percentile(condition,90)),
        "LS_local_cancellation_ratio_max":float(np.max(condition)),
        "DD_normalized_histogram_bitwise_invariant":
            bool(np.array_equal(a["D1D2"],b["D1D2"])),
        "normalized_pair_histogram_4800_vs_48000_relative_L1":{
            k:float(np.abs(b[k]-a[k]).sum()/np.abs(b[k]).sum()) for k in CHANGED},
        "Shapley_exact_matrix_change_L1":{k:float(np.abs(v).sum())
                                          for k,v in shares.items()},
        "Shapley_exact_matrix_max_abs_reconstruction_error":
            float(np.abs(sum(shares.values())-delta).max()),
        "odd":{}
    }
    for ell in (1,3):
        decomposition={term:projected(shares[term],ell) for term in CHANGED}
        before=projected(old,ell)
        after=projected(reference,ell)
        diff=after-before
        require(np.max(np.abs(sum(decomposition.values())-diff))<2e-13,
                "Projected exact Shapley closure failed")
        result["odd"][str(ell)]={
            "at_4800":before.tolist(),
            "at_48000":after.tolist(),
            "change_48000_minus_4800":diff.tolist(),
            "max_abs_change":float(np.abs(diff).max()),
            "change_L1":float(np.abs(diff).sum()),
            "term_order_independent_share":{
                term:decomposition[term].tolist() for term in CHANGED},
        }
    require(math.isclose(result["xi_change_relative_L1_to_48000"],
            case["full_galaxy_4800_to_48000_relative_xi_L1"],
            rel_tol=0,abs_tol=5e-13),
            "Original reported 4800-to-48000 xi L1 differs")
    return result

def synthetic_self_test():
    archive=original_archive()
    rng=np.random.default_rng(27092026)
    a={term:rng.uniform(.25,1.25,size=(6,24)) for term in TERMS}
    b={term:rng.uniform(.3,1.3,size=(6,24)) for term in TERMS}
    b["D1D2"]=a["D1D2"].copy()
    parts,d=shapley(a,b)
    require(np.max(np.abs(sum(parts.values())-d))<2e-13,
            "Synthetic distinct DR/RD/RR denominator closure failed")
    for ell in (1,3):
        require(abs(float(np.ones((1,24))@WEIGHT[ell]))<1e-12,
                "Pure constant even field falsely produces odd output")
        require(np.max(np.abs(sum(projected(v,ell) for v in parts.values())
                    -projected(d,ell)))<2e-13,
                "Synthetic projected term decomposition failed")
    altered=dict(b)
    altered["D1D2"]=b["D1D2"].copy()+.01
    try:shapley(a,altered)
    except ValueError:pass
    else:raise AssertionError("Synthetic changed galaxy DD accepted")
    require(archive["cases"]["0001/NGC"]["xi_relative_L1_difference_to_48000"]>1 and
            archive["cases"]["0001/SGC"]["xi_relative_L1_difference_to_48000"]>1,
            "Original measured drift/certification STOP contract changed")
    print("EBOSS_A03_E7_SOURCE_ONLY_SHAPLEY_CANCELLATION_SELF_TEST_OK",
          "ORIGINAL_ARCHIVE_PINNED NO_FITS OBSERVED_ODD_SEALED",flush=True)

def atomic_new(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,
                "Existing source-only E7 diagnostic differs, do not overwrite")
        return
    with tempfile.NamedTemporaryFile(
            mode="wb",dir=path.parent,prefix=".e7_",delete=False) as f:
        tmp=Path(f.name)
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())
    try:
        require(not path.exists(),"Existing E7 output appeared, do not overwrite")
        os.link(tmp,path)
    finally:
        tmp.unlink(missing_ok=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--report",type=Path,default=ROOT/"eboss_workspace/a03_full_eligible_mock/full_eligible_id0001_4800_48000_resume_v2_report.json")
    ap.add_argument("--output",type=Path,default=OUT_DEFAULT)
    args=ap.parse_args()
    synthetic_self_test()
    if args.self_test:
        return 0
    raw=args.report.read_bytes()
    require(len(raw)==529771 and sha(raw)==REPORT_SHA256,
            "Only the completed user-uploaded real V2.1 mock0001 full JSON is allowed")
    j=json.loads(raw)
    require(j["status"]=="EBOSS_FULL_ELIGIBLE_MOCK0001_SGC_RECOVERY_V2_ENGINEERING_ONLY" and
            j["errors"]==[] and j["completed_cases"]==2 and
            j["parent_v1_NGC_canonical_case_sha256"]==PARENT_NGC_SHA256 and
            j["all_72_original_gzip_rehashed_before_any_SGC_FITS"] is True and
            j["observed_odd_data_vector_read"] is False and
            j["observed_galaxy_rows_read"] is False and
            j["physical_window_certified"] is False and
            j["inferential_covariance_computed"] is False and
            j["wake_detection_or_exclusion_computed"] is False and
            set(j["cases"])=={"0001/NGC","0001/SGC"},
            "Completed two-cap mock-only original SHA or observed seal gate failed")
    result={
        "status":"A03_E7_POSTHOC_SOURCE_ONLY_ALGEBRA_COMPLETE_PHYSICAL_A03_STOP",
        "source_completed_report_sha256":REPORT_SHA256,
        "source_archived_uploaded_audit_git_blob_sha1":ARCHIVE_GIT_BLOB,
        "method":"Exactly average six DR,RD,RR replacement orders with original DD invariant; not sampling variance, physical window, or inference.",
        "signed_mu_original_edges":"np.linspace(-1-1e-7,1+1e-7,25) clipped to [-1,1]",
        "one_mock_only":True,
        "caps":{cap:diagnose(cap,j["cases"]["0001/"+cap])
                for cap in ("NGC","SGC")},
        "physical_absolute_EV_template_available":False,
        "validated_eBOSS_A03_selection_window_available":False,
        "independent_eBOSS_A04_covariance_available":False,
        "48000_random_convergence_certified":False,
        "observed_odd_data_vector_read":False,
        "physical_significance_computed":False,
        "interpretation":"Measured amplification of frozen 4800-to-48000 estimator change due to LS cancellation only; STOP on a numeric physical tolerance.",
    }
    atomic_new(args.output,result)
    print("EBOSS_A03_E7_POSTHOC_SOURCE_ONLY_ALGEBRA_COMPLETE_PHYSICAL_A03_STOP",
          "REPORT",args.output,flush=True)
    for cap,c in result["caps"].items():
        print("E7_CAP",cap,"RELATIVE_XI_L1",c["xi_change_relative_L1_to_48000"],
              "MEDIAN_LOCAL_CANCELLATION_RATIO",
              c["LS_local_cancellation_ratio_median"],flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
