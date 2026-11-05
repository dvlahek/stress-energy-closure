#!/usr/bin/env python3
"""Source-only SHA-pinned A-03E2 18/18 uploaded archive gate; NO FITS reads."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent
REPORT=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json"
MANIFEST=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_uploaded_manifest_2026-09-26.json"
A02=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_protocol_2026-09-26.json"
FIXED_REPORT_SHA="6ef86f904bb0cdc407adabe2e125146ce4d0671ea1ac932e13d11fc481cd596a"
FIXED_REPORT_GIT_BLOB="8596d7322f2ec7ef164e89641958c3912f687722"
FIXED_MANIFEST_GIT_BLOB="7ec2536429604227b2428637ae227b53b7d42cb4"
FIXED_A02_SHA="15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27"
FIXED_PROTOCOL_SHA="eb83baa966e267fdc66fa0b5896aef8df04ff366a1a645ec19d87918809fd80c"
IDS=(1,125,250,375,500,625,750,875,1000)
CAPS=("NGC","SGC")
SCENARIOS=("baseline_original","plus_5pct_ELG_R_z_ramp","minus_5pct_ELG_R_z_ramp")
TERMS=("D1D2","D1R2","R1D2","R1R2")

def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def read(path):return json.loads(path.read_bytes())

def audit(j,p,a,manifest):
    keys=[f"{mid:04d}/{cap}" for mid in IDS for cap in CAPS]
    if (j["status"]!="A03E2_NINE_MOCK_GALAXY_ELG_RANDOM_RADIAL_WEIGHT_STRESS_DESCRIPTIVE_ONLY"
        or j["protocol_sha256"]!=FIXED_PROTOCOL_SHA
        or j["parent_A02_report_sha256"]!=FIXED_A02_SHA
        or j["parent_A03E0E1_report_sha256"]!=p["A03E0E1_original_report_sha256"]
        or j["completed_cases"]!=18 or j["failed_cases"]!=0 or j["errors"]!=[]
        or list(j["cases"])!=keys or j["fixed_ids"]!=list(IDS)
        or j["caps"]!=list(CAPS) or j["scenarios"]!=list(SCENARIOS)
        or j["all_72_original_gzip_full_sha256_rechecked_before_any_E2_FITS_row"] is not True
        or j["total_original_72_gzip_compressed_bytes"]!=1860198719
        or a["completed_cases"]!=18 or a["errors"]!=[]
        or manifest["n_fixed_cases"]!=18 or manifest["n_scenario_cases"]!=54
        or manifest["archive_byte_identical_to_uploaded_source"] is not True
        or any(j.get(flag) is not False for flag in (
            "observed_galaxy_rows_read","observed_random_rows_read",
            "observed_odd_data_vector_read","new_science_selection_applied",
            "physical_window_certified","18D_covariance_computed"))
        or j["not_a_detection_or_exclusion"] is not True):
        raise ValueError("A-03E2 original archive/scope or parent was altered")
    worst=0.
    for key in keys:
        z=j["cases"][key]
        old=a["cases"][key]
        fm=manifest["fixed_case_source_and_output_fingerprints"][key]
        if (z["status"]!="complete"
            or z["original_A02_sample_and_all_eight_pair_SHA_exact_replay"] is not True
            or z["original_A02_forward_xi_SHA_exact_replay"] is not True
            or z["no_observed_rows_or_odd_read"] is not True
            or z["not_a_physical_null_or_selection_bias_calibration"] is not True
            or z["original_selected_array_SHA256"]!=fm["original_A02_selected_array_SHA256_by_tracer_role"]
            or z["original_selected_array_SHA256"]!={k:v["selected_array_SHA256"]
                for k,v in old["input_sample_diagnostics"].items()}):
            raise ValueError("Original selected mock/sample SHA differs: "+key)
        for tr,chunk in z["original_exact_ELG_chunk_counts"].items():
            if chunk!=old["input_sample_diagnostics"][tr]["ELG_exact_chunk_diagnostic"]:
                raise ValueError("Original ELG chunk selection changed: "+key)
        base=z["scenarios"]["baseline_original"]
        if base["forward_xi_grid_sha256"]!=old["cross_ls"]["forward_xi_grid_SHA256"]:
            raise ValueError("Original full 6x24 A-02 xi SHA changed: "+key)
        for side in ("forward_pair_terms","reverse_pair_terms"):
            for term in TERMS:
                r=base[side][term]
                prior=old["cross_ls"][side][term]
                if r["weighted_histogram_SHA256"]!=prior["weighted_histogram_SHA256"] or r["accepted_pairs"]!=prior["accepted_pairs"]:
                    raise ValueError("Original A-02 pair fingerprint changed: "+key+"/"+side+"/"+term)
        for scenario in SCENARIOS:
            s=z["scenarios"][scenario]
            snapshot=fm["scenarios"][scenario]
            stress=s["ELG_R_weight_stress"]
            worst=max(worst,s["max_abs_xi_tracer_reverse_mirror"])
            if (s["RR_supported_cells"]!=144
                or s["status"]!="complete_synthetic_random_weight_stress_scenario"
                or s["all_four_weighted_pair_mirror_closure_passed"] is not True
                or s["observed_odd_vector_read"] is not False
                or s["forward_xi_grid_sha256"]!=snapshot["forward_xi_grid_sha256"]
                or stress["ELG_RANDOM_only_weight_changed"] is not True
                or stress["selected_ELG_R_rows"]!=1200
                or stress["factor_min"]<.95-1e-12 or stress["factor_max"]>1.05+1e-12
                or s["max_abs_xi_tracer_reverse_mirror"]>=1e-8):
                raise ValueError("Failed 144-cell/stress/reversal/source invariant: "+key+"/"+scenario)
            for side,field in (("forward_pair_terms","forward_weighted_pair_histogram_sha256_by_term"),
                               ("reverse_pair_terms","reverse_weighted_pair_histogram_sha256_by_term")):
                if {term:s[side][term]["weighted_histogram_SHA256"] for term in TERMS}!=snapshot[field]:
                    raise ValueError("Modified pair histogram SHA differs: "+key+"/"+scenario+"/"+side)
            if scenario!="baseline_original":
                for side,unaffected in (("forward_pair_terms",("D1D2","R1D2")),
                                        ("reverse_pair_terms",("D1D2","D1R2"))):
                    for term in unaffected:
                        if s[side][term]!=base[side][term]:
                            raise ValueError("Non-ELG-RANDOM pair term changed: "+key)
                for ell in ("0","1","2","3"):
                    before=base["mock_only_descriptive_xi_multipoles"][ell]["values_by_fixed_s_bin"]
                    after=s["mock_only_descriptive_xi_multipoles"][ell]["values_by_fixed_s_bin"]
                    delta=z["diagnostic_deltas_from_original"][scenario][ell]
                    if (len(before)!=6 or len(after)!=6
                        or len(delta["delta_by_fixed_s_bin"])!=6
                        or any(abs((v-u)-d)>1e-12 for v,u,d in zip(after,before,delta["delta_by_fixed_s_bin"]))
                        or abs(max(map(abs,delta["delta_by_fixed_s_bin"]))-delta["max_abs_delta_across_s_bins"])>1e-12):
                        raise ValueError("Fixed six-bin projected delta changed: "+key+"/"+scenario+"/"+ell)
    if abs(worst-manifest["max_abs_xi_reverse_residual"])>1e-20:
        raise ValueError("Manifest tracer-reversal residual changed")
    return worst

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    report=REPORT.read_bytes()
    manifest_raw=MANIFEST.read_bytes()
    if (len(report)!=413189 or sha(report)!=FIXED_REPORT_SHA
        or blob(report)!=FIXED_REPORT_GIT_BLOB
        or blob(manifest_raw)!=FIXED_MANIFEST_GIT_BLOB
        or sha(PROTOCOL.read_bytes())!=FIXED_PROTOCOL_SHA
        or sha(A02.read_bytes())!=FIXED_A02_SHA):
        raise ValueError("Frozen uploaded report/manifest/protocol/parent raw byte SHA changed")
    j=json.loads(report)
    m=json.loads(manifest_raw)
    if (m["exact_original_uploaded_bytes"]!=len(report)
        or m["exact_original_uploaded_SHA256"]!=FIXED_REPORT_SHA
        or m["exact_original_uploaded_git_blob_sha1"]!=FIXED_REPORT_GIT_BLOB
        or m["archived_git_blob_sha1"]!=FIXED_REPORT_GIT_BLOB):
        raise ValueError("Original exact uploaded archive identity mismatch")
    p=read(PROTOCOL)
    a=read(A02)
    worst=audit(j,p,a,m)
    if args.self_test:
        tampered=copy.deepcopy(j)
        tampered["observed_odd_data_vector_read"]=True
        try:audit(tampered,p,a,m)
        except ValueError:pass
        else:raise AssertionError("Tampered observed-odd-access flag accepted")
        tampered=copy.deepcopy(j)
        tampered["cases"]["0001/NGC"]["scenarios"]["baseline_original"]["forward_xi_grid_sha256"]="0"*64
        try:audit(tampered,p,a,m)
        except ValueError:pass
        else:raise AssertionError("Tampered original A-02 xi SHA accepted")
    print("A03E2_IMMUTABLE_UPLOADED_413189_BYTE_REPORT_18_OF_18_ARCHIVE_OK")
    print("SCENARIOS",54,"RR_SUPPORTED_EACH",144)
    print("MAX_XI_REVERSE_RESIDUAL",worst)
    print("OBSERVED_ODD_READ",False)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
