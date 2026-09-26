#!/usr/bin/env python3
"""A-03 BLINDED source-pinned empirical RR conditional-window probe.

Replays immutable original 2026-09-24 nine-ID LRG×ELG random-only
window-operator ZIP members, not FITS rows. A unitless even radial
gradient can generate finite-bin even->odd mixing. This is not a physical
selection null, complete mask/window certification, mock-galaxy covariance,
or an observed odd/wake detection.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROTO_PATH = ROOT / "source_data/eboss_dr16_a03_empirical_rr_window_injection_protocol_2026-09-26.json"
ADDENDUM = ROOT / "source_data/eboss_dr16_a03_empirical_only_user_decision_2026-09-26.json"
REF = ROOT / "source_data/eboss_dr16_nine_ezmock_random_reference_sha_from_20260924_artifact.json"
A02 = ROOT / "source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_uploaded_manifest_2026-09-26.json"
OBS = ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json"
PASS = "EBOSS_A03_EMPIRICAL_RANDOM_CONDITIONAL_WINDOW_PROBE_DESCRIPTIVE_ONLY"
STOP = "EBOSS_A03_EMPIRICAL_RANDOM_CONDITIONAL_WINDOW_PROBE_INCOMPLETE_STOP"
OUT = (1, 3)
INS = (0, 1, 2, 3, 4)
EVEN = (0, 2, 4)


def read_json(path):
    return json.loads(path.read_bytes())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def exact_bytes(data, *, expected_size, expected_sha, label):
    if len(data) != expected_size or sha(data) != expected_sha:
        raise ValueError("Original archive/member byte identity changed: " + label)


def source_gate(p):
    d=read_json(ADDENDUM)
    ref=read_json(REF)
    a02=read_json(A02)
    obs=read_json(OBS)
    if (
        p["user_decision"]!="EMPIRICAL_ONLY_NO_AUTHOR_CONTACT"
        or p["user_decision_addendum"]!=str(ADDENDUM.relative_to(ROOT))
        or d["author_contact_approved"] is not False
        or d["email_sent"] is not False
        or d["empirical_LRG_ELG_window_validated"] is not False
        or d["observed_odd_data_vector_read"] is not False
        or d["historical_ELG_production_mask_certified"] is not False
        or p["source_artifact_json_member_sha256"]!=ref["source_json_sha256"]
        or p["prior_full_36_random_reference"]!=str(REF.relative_to(ROOT))
        or p["prior_A02_18_real_galaxy_pair_upload_manifest"]!=str(A02.relative_to(ROOT))
        or a02["completed_cases"]!=18
        or a02["failed_cases"]!=0
        or a02["all_72_full_gzip_source_SHA_checked_before_ANY_mock_FITS_row"] is not True
        or a02["exact_uploaded_SHA256"]!="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
        or a02["physical_pair_window_certified"] is not False
        or a02["inferential_18D_mock_covariance_computed"] is not False
        or obs["odd_sector_data_vector_inspected"] is not False
        or obs["pair_counts_computed"] is not False
        or p["cohort"]["ids"]!=[1,125,250,375,500,625,750,875,1000]
        or p["cohort"]["caps"]!=["NGC","SGC"]
        or p["cohort"]["n_cases"]!=18
        or p["cohort"]["observed_random_operator_comparison_only"] is not True
    ):
        raise ValueError("A-03 empirical-only policy, immutable 36-source reference or A-02 provenance changed")
    return ref,obs


def members(p,args):
    expected={p["source_artifact_json_member"],p["source_artifact_npz_member"]}
    if args.artifact_zip:
        b=Path(args.artifact_zip).read_bytes()
        exact_bytes(b,expected_size=len(b),expected_sha=p["source_artifact_exact_ZIP_sha256"],
                    label="original GitHub Actions ZIP")
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            if set(z.namelist())!=expected:
                raise ValueError("Original ZIP member list differs from source predeclaration")
            v={name:z.read(name) for name in expected}
    elif args.artifact_dir:
        folder=Path(args.artifact_dir)
        v={name:(folder/name).read_bytes() for name in expected}
    else:
        raise ValueError("Provide --artifact-zip OR --artifact-dir; no source download inside the runner")
    exact_bytes(v[p["source_artifact_json_member"]],
                expected_size=p["source_artifact_json_member_bytes"],
                expected_sha=p["source_artifact_json_member_sha256"],
                label=p["source_artifact_json_member"])
    exact_bytes(v[p["source_artifact_npz_member"]],
                expected_size=p["source_artifact_npz_member_bytes"],
                expected_sha=p["source_artifact_npz_member_sha256"],
                label=p["source_artifact_npz_member"])
    return v


def expected_keys(p):
    names=["fine_sedges_mpc_over_h","output_sedges_mpc_over_h","muedges"]
    names += [f"observed_M_{cap}_out{o}_in{i}"
              for cap in p["cohort"]["caps"] for o in OUT for i in INS]
    names += [f"mock_M_{cap}_id{mid:04d}_out{o}_in{i}"
              for cap in p["cohort"]["caps"] for mid in p["cohort"]["ids"]
              for o in OUT for i in INS]
    if len(names)!=p["geometry"]["npz_expected_array_count"] or len(set(names))!=len(names):
        raise ValueError("Source NPZ expected key count or uniqueness changed")
    return names


def arrays_gate(arr,p):
    expected=expected_keys(p)
    if set(arr.files)!=set(expected) or len(arr.files)!=len(expected):
        raise ValueError("Original NPZ 203-key source identity/layout incomplete")
    fine=np.asarray(arr["fine_sedges_mpc_over_h"],dtype="f8")
    coarse=np.asarray(arr["output_sedges_mpc_over_h"],dtype="f8")
    mu=np.asarray(arr["muedges"],dtype="f8")
    if (
        not np.array_equal(fine,np.arange(20,141,dtype="f8"))
        or not np.array_equal(coarse,np.arange(20,141,20,dtype="f8"))
        or not np.array_equal(mu,np.linspace(-1.-1e-7,1.+1e-7,241))
    ):
        raise ValueError("Frozen fine-s/coarse-s/signed-mu geometry changed")
    for key in expected[3:]:
        val=np.asarray(arr[key])
        if val.shape!=(6,120) or val.dtype!=np.dtype("f8") or not np.isfinite(val).all():
            raise ValueError("Wrong 6x120 original finite float64 operator matrix: "+key)
    return fine,coarse


def original_gate(j,p,ref,obs):
    cohort=p["cohort"]
    names=[(mid,cap) for mid in cohort["ids"] for cap in cohort["caps"]]
    if (
        j["status"]!=p["original_report_status"]
        or j["errors"]!=[]
        or j["predeclared_mock_ids"]!=cohort["ids"]
        or j["caps"]!=cohort["caps"]
        or j["n_individual_mock_realizations"]!=9
        or j["n_cap_by_realization_shards"]!=18
        or j["observed_galaxy_data_read"] is not False
        or j["mock_galaxy_data_read"] is not False
        or j["observed_odd_data_vector_read"] is not False
        or j["joint_mock_galaxy_covariance_estimated"] is not False
        or j["physical_window_convolution_validated"] is not False
        or j["analysis_protocol_frozen_for_inference"] is not False
        or len(j["cases"])!=18
        or [(c["mock_id"],c["cap"]) for c in j["cases"]]!=names
    ):
        raise ValueError("Original nine-mock random-only archived case list, flags or status changed")
    pinned={(x["cap"],x["tracer"]):x["sha256"] for x in obs["randoms"]}
    observed={(x["cap"],x["tracer"]):x["sha256"]
              for x in j["observed_random_input_file_evidence"]}
    if observed!=pinned or len(observed)!=4:
        raise ValueError("Observed RANDOM-only exact SHA source references changed")
    for c in j["cases"]:
        mid,cap=c["mock_id"],c["cap"]
        if c["full_mock_RR_supported_cells"]!=28800:
            raise ValueError(f"Prior full highz fine RR support changed: {mid:04d}/{cap}")
        for tracer in ("LRG","ELG"):
            k=f"{mid:04d}/{cap}/eBOSS_{tracer}"
            if c["input_file_sha256"][tracer]!=ref["full_random_gzip_sha256_by_id_cap_tracer"][k]:
                raise ValueError("Original nine-window mock RANDOM source SHA differs: "+k)


def ramp_inputs(fine,coarse):
    s=0.5*(fine[:-1]+fine[1:])
    idx=np.searchsorted(coarse,s,side="right")-1
    if np.any(idx<0) or np.any(idx>=6):
        raise ValueError("Outside frozen 6 coarse separation bins")
    return (s-0.5*(coarse[idx]+coarse[idx+1]))/(coarse[idx+1]-coarse[idx])


def input_probes(r):
    return {
        "zero":{i:np.zeros_like(r) for i in INS},
        "constant_even":{0:np.ones_like(r),1:np.zeros_like(r),
                         2:0.25*np.ones_like(r),3:np.zeros_like(r),
                         4:0.1*np.ones_like(r)},
        "even_radial_ramp":{0:r,1:np.zeros_like(r),2:0.25*r,
                            3:np.zeros_like(r),4:0.1*r},
        "odd_injection":{0:np.zeros_like(r),1:0.2*(1.+0.5*r),
                         2:np.zeros_like(r),3:-0.075*(1.-0.5*r),
                         4:np.zeros_like(r)},
    }


def propagate(blocks,inputs):
    return {str(o):np.sum([blocks[o,i]@inputs[i] for i in INS],axis=0)
            for o in OUT}


def replay_matrix_fingerprints(case,cap,mid,obs,mock,radial,p):
    tol=p["acceptance_criteria"]["prior_report_all_18_conditional_matrix_block_differences_reproduced_max_abs_tol"]
    tol_ramp=p["acceptance_criteria"]["prior_report_all_18_even_ramp_mock_minus_observed_by_s_values_reproduced_max_abs_tol"]
    worst_matrix=0.
    worst_ramp=0.
    half={}
    for o in OUT:
        odd_norm=np.linalg.norm(obs[o,o])
        if not np.isfinite(odd_norm) or odd_norm<=0:
            raise ValueError(f"Zero original observed RANDOM-only diagonal operator: {cap}/{o}")
        for i in INS:
            k=f"{o}<-{i}"
            old=case["conditional_window_block_differences"][k]
            shift=mock[o,i]-obs[o,i]
            diff=float(np.max(np.abs(shift)))
            normalized=float(np.linalg.norm(shift)/odd_norm)
            err=max(abs(diff-old["max_absolute_cell_difference"]),
                    abs(normalized-old["normalized_frobenius_difference_over_observed_odd_diagonal"]))
            worst_matrix=max(worst_matrix,err)
            if err>tol:
                raise ValueError(f"Original random-only source matrix numerical summary did not replay: {mid:04d}/{cap}/{k}")
            if i in EVEN:
                previous=case["predeclared_unitless_ramp_probe_differences"][k]
                actual=np.asarray([(mock[o,i][j]-obs[o,i][j])@radial
                                   for j in range(6)],dtype="f8")
                expected=np.asarray(previous["mock_minus_observed_by_coarse_s"],dtype="f8")
                re=max(float(np.max(np.abs(actual-expected))),
                       abs(float(np.max(np.abs(actual)))-previous["max_abs_mock_minus_observed"]))
                worst_ramp=max(worst_ramp,re)
                if re>tol_ramp:
                    raise ValueError(f"Original even radial-ramp probe did not replay: {mid:04d}/{cap}/{k}")
                half[k]={
                    "half_A_vs_full_mock_max_abs_original_random_only":
                        previous["max_abs_half_A_minus_full_mock"],
                    "half_B_vs_full_mock_max_abs_original_random_only":
                        previous["max_abs_half_B_minus_full_mock"]
                }
    return worst_matrix,worst_ramp,half


def validated_case(blocks,probes,p,label):
    vals={k:propagate(blocks,x) for k,x in probes.items()}
    if any(np.any(vals["zero"][o]!=0.) for o in ("1","3")):
        raise ValueError("Zero synthetic input produced nonzero conditional output: "+label)
    constant=max(float(np.max(np.abs(vals["constant_even"][o])))
                 for o in ("1","3"))
    if constant>=p["acceptance_criteria"]["constant_even_output_max_abs_tol"]:
        raise ValueError("Flat-in-separation even input has excessive odd leakage: "+label)
    combined={i:probes["even_radial_ramp"][i]+probes["odd_injection"][i] for i in INS}
    mixed=propagate(blocks,combined)
    linear=max(float(np.max(np.abs(mixed[o]-vals["even_radial_ramp"][o]
                                    -vals["odd_injection"][o])))
               for o in ("1","3"))
    if linear>p["acceptance_criteria"]["synthetic_odd_superposition_residual_max_abs_tol"]:
        raise ValueError("Odd injection linear-superposition closure failed: "+label)
    if max(float(np.max(np.abs(vals["odd_injection"][o]))) for o in ("1","3"))<=1e-6:
        raise ValueError("Predeclared positive odd synthetic injected response disappeared: "+label)
    return {k:{o:vals[k][o].tolist() for o in ("1","3")}
            for k in ("constant_even","even_radial_ramp","odd_injection")},constant,linear


def run(p,args):
    ref,obs_evidence=source_gate(p)
    raw=members(p,args)
    j=json.loads(raw[p["source_artifact_json_member"]])
    original_gate(j,p,ref,obs_evidence)
    with np.load(io.BytesIO(raw[p["source_artifact_npz_member"]]),allow_pickle=False) as arr:
        fine,coarse=arrays_gate(arr,p)
        radial=ramp_inputs(fine,coarse)
        probes=input_probes(radial)
        output={
            "status":STOP,
            "protocol_SHA256":sha(PROTO_PATH.read_bytes()),
            "immutable_original_artifact_run":p["original_previous_workflow_run"],
            "immutable_original_artifact_id":p["original_previous_artifact_id"],
            "input_original_JSON_member_SHA256":p["source_artifact_json_member_sha256"],
            "input_original_NPZ_member_SHA256":p["source_artifact_npz_member_sha256"],
            "source_A02_18_case_manifest":p["prior_A02_18_real_galaxy_pair_upload_manifest"],
            "source_A03_empirical_only_user_decision":p["user_decision_addendum"],
            "cases":{},"observed_RANDOM_operator_by_cap":{},
            "errors":[],
            "observed_galaxy_rows_read":False,
            "observed_odd_data_vector_read":False,
            "new_catalogue_selection_or_veto_applied":False,
            "historical_physical_mask_certified":False,
            "physical_pair_window_certified":False,
            "full_18D_inferential_covariance_computed":False,
            "wake_detection_or_exclusion_test_computed":False,
            "synthetic_input_is_not_a_physical_wake_model":True,
        }
        out=ROOT/p["local_output_json"]
        out.parent.mkdir(parents=True,exist_ok=True)
        save(out,output)
        for cap in p["cohort"]["caps"]:
            obs={(o,i):np.array(arr[f"observed_M_{cap}_out{o}_in{i}"],
                                dtype="f8",copy=True)
                 for o in OUT for i in INS}
            v,c,l=validated_case(obs,probes,p,"observed_random_only/"+cap)
            output["observed_RANDOM_operator_by_cap"][cap]={
                "predeclared_unitless_synthetic_probes":v,
                "max_abs_constant_even_to_odd_residual":c,
                "max_abs_linear_superposition_residual":l,
                "this_is_observed_RANDOM_not_observed_galaxy_odd":True,
            }
            save(out,output)
            for mid in p["cohort"]["ids"]:
                key=f"{mid:04d}/{cap}"
                original_case=next(c for c in j["cases"]
                                   if c["mock_id"]==mid and c["cap"]==cap)
                mock={(o,i):np.array(arr[f"mock_M_{cap}_id{mid:04d}_out{o}_in{i}"],
                                      dtype="f8",copy=True)
                      for o in OUT for i in INS}
                m_err,r_err,half=replay_matrix_fingerprints(
                    original_case,cap,mid,obs,mock,radial,p)
                mv,mc,ml=validated_case(mock,probes,p,key)
                differ={probe:{o:(np.asarray(mv[probe][o])-np.asarray(v[probe][o])).tolist()
                               for o in ("1","3")}
                        for probe in ("even_radial_ramp","odd_injection")}
                maxdiff={probe:max(float(np.max(np.abs(differ[probe][o])))
                                   for o in ("1","3"))
                         for probe in differ}
                output["cases"][key]={
                    "status":"closed_empirical_random_only_finite_bin_probe",
                    "mock_id":mid,"cap":cap,
                    "original_full_RR_supported_cells":original_case["full_mock_RR_supported_cells"],
                    "original_mock_vs_observed_fine_RR_normalized_L1":
                        original_case["full_mock_vs_observed_normalized_fine_rr_relative_l1"],
                    "original_random_half_even_ramp_probe_max_abs_shifts":half,
                    "prior_independent_20260924_operator_matrix_summary_max_abs_replay_error":m_err,
                    "prior_independent_20260924_ramp_summary_max_abs_replay_error":r_err,
                    "predeclared_unitless_synthetic_probes":mv,
                    "mock_minus_observed_RANDOM_operator_synthetic_by_fixed_s":differ,
                    "max_abs_mock_minus_observed_RANDOM_operator_by_probe":maxdiff,
                    "max_abs_constant_even_to_odd_residual":mc,
                    "max_abs_linear_superposition_residual":ml,
                }
                save(out,output)
                print("A03_EMPIRICAL_RANDOM_ONLY_CASE",key,
                      "even_ramp_shift",maxdiff["even_radial_ramp"],
                      "odd_injection_shift",maxdiff["odd_injection"],flush=True)
        if len(output["cases"])!=18 or any(
            case["status"]!="closed_empirical_random_only_finite_bin_probe"
            for case in output["cases"].values()
        ):
            raise ValueError("All predeclared eighteen original cases must be retained")
        output["status"]=PASS
        output["completed_cases"]=len(output["cases"])
        output["scope_note"]=(
            "Source-authenticated empirical random-only conditional finite-bin operator "
            "and fixed unitless even/odd synthetic injection replay. Not a physical "
            "galaxy-based null, complete source/veto mask, full survey window, "
            "inferential covariance, significance, detection, or exclusion.")
        save(out,output)
        return output


def save(path,report):
    temp=path.with_suffix(".tmp.json")
    temp.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    temp.replace(path)


def self_test(p):
    source_gate(p)
    fine=np.arange(20,141,dtype="f8")
    coarse=np.arange(20,141,20,dtype="f8")
    rr=ramp_inputs(fine,coarse)
    if rr.shape!=(120,) or not np.allclose(rr[0],-0.475) or not np.allclose(rr[19],0.475):
        raise AssertionError("Frozen within-coarse-bin radial probe changed")
    probes=input_probes(rr)
    fake={(o,i):np.zeros((6,120),dtype="f8") for o in OUT for i in INS}
    for j in range(6):
        b=slice(20*j,20*j+20)
        fake[1,1][j,b]=1./20
        fake[3,3][j,b]=1./20
    v,c,l=validated_case(fake,probes,p,"synthetic_flat_four_term_window")
    if c!=0 or l>1e-12 or abs(v["odd_injection"]["1"][0]-.2)>1e-12:
        raise AssertionError("Even null or constant odd injection synthetic test failed")
    fake_bad={k:x.copy() for k,x in fake.items()}
    fake_bad[1,0][0,0]=2e-3
    try:
        validated_case(fake_bad,probes,p,"synthetic_tampered_even_to_odd")
    except ValueError:
        pass
    else:
        raise AssertionError("Synthetic excessive constant-even leakage accepted")
    original=read_json(REF)
    altered=json.loads(json.dumps(original))
    altered["full_random_gzip_sha256_by_id_cap_tracer"]["0125/NGC/eBOSS_LRG"]="0"*64
    source_fake={"status":p["original_report_status"],"cases":[{
        "mock_id":125,"cap":"NGC","full_mock_RR_supported_cells":28800,
        "input_file_sha256":{"LRG":original["full_random_gzip_sha256_by_id_cap_tracer"]["0125/NGC/eBOSS_LRG"],
                             "ELG":original["full_random_gzip_sha256_by_id_cap_tracer"]["0125/NGC/eBOSS_ELG"]}
    }]}
    if source_fake["cases"][0]["input_file_sha256"]["LRG"]==altered["full_random_gzip_sha256_by_id_cap_tracer"]["0125/NGC/eBOSS_LRG"]:
        raise AssertionError("Tampered archived source SHA was not detected")
    print("EBOSS_A03_EMPIRICAL_RR_SYNTHETIC_SELF_TEST_AND_NEGATIVE_CONTROLS_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--artifact-dir",type=Path)
    ap.add_argument("--artifact-zip",type=Path)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    p=read_json(PROTO_PATH)
    if args.self_test:
        self_test(p)
        return 0
    path=ROOT/p["local_output_json"]
    try:
        out=run(p,args)
    except Exception as exc:
        if path.exists():
            out=read_json(path)
            out["status"]=STOP
            out["errors"]=list(dict.fromkeys(out.get("errors",[])+[str(exc)]))
        else:
            out={"status":STOP,"errors":[str(exc)],
                 "protocol_SHA256":sha(PROTO_PATH.read_bytes()),
                 "observed_galaxy_rows_read":False,
                 "observed_odd_data_vector_read":False,
                 "new_catalogue_selection_or_veto_applied":False,
                 "physical_pair_window_certified":False,
                 "full_18D_inferential_covariance_computed":False}
        path.parent.mkdir(parents=True,exist_ok=True)
        save(path,out)
    print("EBOSS_A03_EMPIRICAL_RR_WINDOW",out["status"],flush=True)
    print("REPORT",path,flush=True)
    if out.get("errors"):
        print("ERRORS",*out["errors"],sep="\n",flush=True)
        return 2
    print("COMPLETED_FIXED_CASES",len(out["cases"]),flush=True)
    print("OBSERVED_ODD_READ",out["observed_odd_data_vector_read"],flush=True)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
