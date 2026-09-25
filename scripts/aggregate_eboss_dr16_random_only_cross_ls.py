#!/usr/bin/env python3
"""Aggregate the complete predeclared nine-mock blinded cross-LS random-only gate.

One fully source-pinned random-only cross-LS report for each of nine
realizations and both DR16 caps is required. The repeated observed
random split must have identical sampled source fingerprints within
each cap. The result tests *estimator algebra only*, never a measured
galaxy odd signal or mock-galaxy covariance.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np

from aggregate_eboss_dr16_9mock_window import IDS, CAPS, TRACERS
from audit_eboss_dr16_9mock_fine_weighted_nz import preflight
from audit_eboss_dr16_random_only_cross_ls_pilot import PROTOCOL_PATH
from audit_eboss_dr16_9mock_fine_weighted_nz import json_write_atomic

LABELS = ("D1D2", "D1R2", "R1D2", "R1R2")
EXPECTED_N = len(IDS)*len(CAPS)


def case_gate(case, kind, cap, mid):
    if (case.get("status")!="split_random_cross_ls_algebra_checked"
            or case.get("kind")!=kind or case.get("cap")!=cap
            or case.get("mock_id")!=(mid if kind=="mock" else None)
            or tuple(case.get("pair_terms",[]))!=LABELS
            or case.get("pseudo_xi_not_an_observed_signal") is not True
            or case.get("pseudo_D_and_R_from_same_parent_random_catalogues")
                is not True
            or case.get("measured_galaxy_odd_data_vector_read") is not False
            or case.get("supported_s_mu_cells",0)<1):
        raise ValueError("Random-only cross-LS case missing fixed data-blind gate")
    rr=case.get("forward_reverse_pair_closure",{})
    if set(rr)!=set(LABELS):
        raise ValueError("One or more independently oriented cross-pair terms absent")
    for label in LABELS:
        x=rr[label]
        if (x.get("closure_passed") is not True
                or x.get("mirror_l1_relative") is None
                or x["mirror_l1_relative"]>=1e-10
                or x.get("pair_normalization_relative_residual",1)>=1e-12):
            raise ValueError("Cross-pair swap or weight normalization failed: "+label)
    if not 0<=case.get("reverse_cross_ls_max_abs_residual",1)<1e-8:
        raise ValueError("Cross-LS field mirror residual fails pinned tolerance")
    odd=case.get("rr_raw_odd_parity",{})
    if set(odd)!={"1","3"} or any(
            not 0<=odd[k].get("reflection_sum_abs",1)<1e-10
            for k in ("1","3")):
        raise ValueError("Raw random odd-moment orientation closure failed")


def aggregate(input_dir,ensemble_json):
    proto=json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    if (tuple(proto["mock_ids"])!=IDS or tuple(proto["caps"])!=CAPS
            or tuple(proto["tracers"])!=TRACERS
            or proto.get("observed_odd_data_vector_read") is not False):
        raise ValueError("Prospective nine-mock cross-LS protocol changed")
    sha_mock,sha_obs,_,_=preflight(
        argparse.Namespace(ensemble_json=ensemble_json))
    seen={}
    observed_fingerprints={}
    cases=[]
    for path in sorted(input_dir.rglob("random_only_cross_ls_pilot.json")):
        x=json.loads(path.read_text(encoding="utf-8"))
        mid,cap=x.get("mock_id"),x.get("cap")
        key=(mid,cap)
        if key in seen:
            raise ValueError("Duplicate fixed mock/cap result: "+str(key))
        if (mid not in IDS or cap not in CAPS or
                x.get("status")!="split_random_cross_ls_random_only_pilot_complete"
                or x.get("errors") or x.get("source_ensemble_run")!=36017670812
                or x.get("observed_galaxy_data_read") is not False
                or x.get("mock_galaxy_data_read") is not False
                or x.get("observed_odd_data_vector_read") is not False
                or x.get("not_statistical_inference") is not True
                or len(x.get("cases",[]))!=2):
            raise ValueError("Unexpected, unblinded or incomplete shard: "+str(path))
        by_kind={v["kind"]:v for v in x["cases"]}
        if set(by_kind)!={"observed","mock"}:
            raise ValueError("Expected observed-random and mock-random pseudo-D/R cases")
        for kind in ("observed","mock"):
            case_gate(by_kind[kind],kind,cap,mid)
            for tracer in TRACERS:
                item=x["inputs"][kind][tracer]
                sha=(sha_obs[cap,tracer] if kind=="observed"
                     else sha_mock[mid,cap,tracer])
                if (item.get("source_sha256")!=sha
                        or item.get("disjoint_split") is not True
                        or item.get("pseudo_D_rows")!=600
                        or item.get("pseudo_R_rows")!=600
                        or not isinstance(item.get("sample_arrays_sha256"),str)
                        or len(item["sample_arrays_sha256"])!=64):
                    raise ValueError("Fixed source SHA or 600+600 split changed")
                if kind=="observed":
                    obskey=(cap,tracer)
                    fp=item["sample_arrays_sha256"]
                    if (obskey in observed_fingerprints
                            and observed_fingerprints[obskey]!=fp):
                        raise ValueError("Repeated observed random split fingerprint changed")
                    observed_fingerprints[obskey]=fp
        seen[key]=path
        cases.append({
            "mock_id":mid,"cap":cap,
            "random_only_paired_cases":2,
            "observed_random_reverse_xi_max_abs_residual":
                by_kind["observed"]["reverse_cross_ls_max_abs_residual"],
            "mock_random_reverse_xi_max_abs_residual":
                by_kind["mock"]["reverse_cross_ls_max_abs_residual"],
            "mock_RR_supported_cells":
                by_kind["mock"]["supported_s_mu_cells"],
            "source_sha256":{
                tracer:x["inputs"]["mock"][tracer]["source_sha256"]
                for tracer in TRACERS},
        })
    if (set(seen)!={(mid,cap) for mid in IDS for cap in CAPS}
            or len(seen)!=EXPECTED_N or len(observed_fingerprints)!=4):
        raise ValueError(
            "Incomplete fixed nine-mock, two-cap random-only cross-LS cohort: "
            +str(sorted(seen)))
    cases.sort(key=lambda c:(IDS.index(c["mock_id"]),CAPS.index(c["cap"])))
    summary={}
    for cap in CAPS:
        vals=[c for c in cases if c["cap"]==cap]
        if len(vals)!=len(IDS):
            raise ValueError("Missing predeclared cap subgroup")
        err=np.array([c["mock_random_reverse_xi_max_abs_residual"]
                      for c in vals],dtype="f8")
        summary[cap]={
            "fixed_mock_count":len(vals),
            "mock_random_reverse_xi_max_abs_residual_min":float(err.min()),
            "mock_random_reverse_xi_max_abs_residual_median":float(np.median(err)),
            "mock_random_reverse_xi_max_abs_residual_max":float(err.max()),
            "all_mock_RR_positive_supported_cells":all(
                c["mock_RR_supported_cells"]>0 for c in vals),
        }
    return {
        "status":"nine_mock_random_only_cross_ls_algebra_complete",
        "mock_ids":list(IDS),"caps":list(CAPS),"tracers":list(TRACERS),
        "n_distinct_mock_realizations":len(IDS),
        "n_cap_by_mock_shards":len(cases),
        "observed_random_split_fingerprints_by_cap_tracer":{
            f"{cap}_{tracer}":observed_fingerprints[cap,tracer]
            for cap in CAPS for tracer in TRACERS},
        "cases":cases,"summary":summary,"errors":[],
        "scope":"All 18 fixed random-only split pseudo-D/R cases passed oriented cross-LS algebra; not independent mock-galaxy closure, a physical mask, a valid observational covariance, a physical wake test, or a detection.",
        "observed_galaxy_data_read":False,
        "mock_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
        "physical_mask_certified":False,
        "physical_window_convolution_validated":False,
        "joint_mock_galaxy_covariance_estimated":False,
        "inference_protocol_frozen":False,
    }


def self_test():
    proto=json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    assert tuple(proto["mock_ids"])==IDS
    assert EXPECTED_N==18
    fake={"status":"split_random_cross_ls_algebra_checked",
          "kind":"mock","cap":"NGC","mock_id":1,
          "pair_terms":list(LABELS),
          "pseudo_xi_not_an_observed_signal":True,
          "pseudo_D_and_R_from_same_parent_random_catalogues":True,
          "measured_galaxy_odd_data_vector_read":False,
          "supported_s_mu_cells":144,
          "reverse_cross_ls_max_abs_residual":0.,
          "forward_reverse_pair_closure":{
              k:{"closure_passed":True,"mirror_l1_relative":0.,
                 "pair_normalization_relative_residual":0.}
              for k in LABELS},
          "rr_raw_odd_parity":{
              k:{"reflection_sum_abs":0.} for k in ("1","3")}}
    case_gate(fake,"mock","NGC",1)
    fake["rr_raw_odd_parity"]["1"]["reflection_sum_abs"]=1e-4
    try:
        case_gate(fake,"mock","NGC",1)
    except ValueError:
        pass
    else:
        raise AssertionError("Fake broken random odd-parity case was accepted")
    print("EBOSS_NINEMOCK_RANDOM_ONLY_CROSS_LS_AGGREGATOR_SELF_TEST_OK")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir",type=Path)
    ap.add_argument("--ensemble-json",type=Path)
    ap.add_argument("--out",type=Path,
                    default=Path("eboss_workspace/aggregate_random_cross_ls.json"))
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.input_dir is None or args.ensemble_json is None:
        ap.error("Require both --input-dir and --ensemble-json")
    result=aggregate(args.input_dir,args.ensemble_json)
    json_write_atomic(args.out,result)
    print("EBOSS_NINEMOCK_RANDOM_ONLY_CROSS_LS_AGGREGATE_OK",
          result["n_cap_by_mock_shards"],args.out,flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
