#!/usr/bin/env python3
"""Independent pure-stdlib E17D2b2 published Abacus c000 metadata bridge replay.

Re-parse immutable CLASS source assignments with AST, independently evaluate
the primordial amplitude using decimal exp, verify every original report
comparison, hard science STOP and negative original SHA. No original worker,
external dataset, Abacus ASDF, CLASS, NumPy, merger tree or observed galaxy.
"""
from __future__ import annotations
import argparse
import ast
import decimal
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b2_abacussummit_c000_metadata_bridge_prereg_2026-09-28.json"
PROTO_BLOB="c1079b666fa8b91af42a8191dbddf56966ea3221"
COSMO=ROOT/"code/class_response_optimize.py"
BASE=ROOT/"code/wake_two_tracer_fisher.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17A_PROTO=ROOT/"source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json"
B1D=ROOT/"source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28"
B1=B1D/"e17d2b1_original_joint_conditional_exterior_assembly_rates.json"
B1I=B1D/"e17d2b1_independent_stdlib_exterior_assembly_replay.json"
B1M=B1D/"archive_manifest.json"
B1P=ROOT/"source_data/eboss_dr16_a03_e17d2b1_conditional_wechsler_exterior_halo_assembly_prereg_2026-09-28.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b2_c000_metadata"

def check(ok,why):
    if not ok:raise ValueError("E17D2B2_INDEPENDENT_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def gate():
    check(blob(PRE.read_bytes())==PROTO_BLOB,"prospective metadata Git blob drift")
    p=json.loads(PRE.read_bytes())
    q=p["immutable_original_repo_parents"]
    for path,key in ((COSMO,"original_CLASS_cosmo_source_code_blob"),
                     (BASE,"original_CLASS_parameter_builder_code_blob"),
                     (E17A_PROTO,"original_E17A_protocol_git_blob"),
                     (B1P,"original_E17D2b1_protocol_git_blob"),
                     (B1M,"original_E17D2b1_manifest_git_blob")):
        check(blob(path.read_bytes())==q[key],"original Git blob "+key)
    for path,key in ((E8,"original_E8_F0_Fplus_Fminus_4000q_CSV_sha256"),
                     (E16,"original_E16_576_geometry_sha256"),
                     (E17A,"original_E17A_CLASS_joint_sha256"),
                     (B1,"original_E17D2b1_joint_sha256"),
                     (B1I,"original_E17D2b1_independent_sha256")):
        check(sha(path.read_bytes())==q[key],"original science SHA "+key)
    a=json.loads(E17A.read_bytes());b=json.loads(B1.read_bytes())
    check(a["geometry_count"]==576 and a["eBOSS_observed_odd_read"] is False
          and a["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED"
          and b["observed_odd_SEALED"] is True
          and b["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and all(p["absolute_STOP"].values())
          and all(bool(v) for v in p["preregistered_QA"].values()),
          "frozen original E16/E17A/E17D2b1 or absolute scope STOP")
    return p

def parse_assignments(path):
    tree=ast.parse(path.read_text(encoding="utf-8"))
    return {target.id:node.value
            for node in tree.body
            if isinstance(node,ast.Assign)
            for target in node.targets
            if isinstance(target,ast.Name)}

def literal_number(assign,name):
    check(name in assign,"missing original CLASS source assignment: "+name)
    node=assign[name]
    check(isinstance(node,ast.Constant) and type(node.value) in (float,int),
          "CLASS source numeric literal mutated: "+name)
    return float(node.value)

def independent_original(p):
    ca=parse_assignments(COSMO)
    ba=parse_assignments(BASE)
    amp=ca["A_S"]
    check(isinstance(amp,ast.Call) and len(amp.args)==1
          and isinstance(amp.func,ast.Name) and amp.func.id=="float",
          "original A_s must remain pinned numpy exp float")
    expr=amp.args[0]
    check(isinstance(expr,ast.BinOp) and isinstance(expr.op,ast.Mult)
          and isinstance(expr.left,ast.Call)
          and isinstance(expr.left.func,ast.Attribute)
          and expr.left.func.attr=="exp"
          and isinstance(expr.left.func.value,ast.Name)
          and expr.left.func.value.id=="np"
          and len(expr.left.args)==1
          and isinstance(expr.left.args[0],ast.Constant)
          and expr.left.args[0].value==3.044
          and isinstance(expr.right,ast.Constant)
          and expr.right.value==1e-10,
          "frozen original CLASS A_s formula changed")
    with decimal.localcontext() as ctx:
        ctx.prec=40
        As=float(decimal.Decimal("3.044").exp()*
                 decimal.Decimal("1e-10"))
    params={"h":literal_number(ca,"H0")/100.,
            "omega_b":literal_number(ca,"OMEGA_B"),
            "omega_cdm":literal_number(ca,"OMEGA_CDM"),
            "n_s":literal_number(ba,"NS"),
            "N_ur":literal_number(ca,"N_UR"),
            "N_ncdm":1,
            "A_s":As,
            "tau_reio":literal_number(ca,"TAU_REIO")}
    evidence=p["original_CLASS_literal_frozen_parameter_evidence"]
    check(params["h"]==evidence["H0_km_s_Mpc"]/100.
          and params["n_s"]==evidence["n_s"]
          and params["omega_b"]==evidence["omega_b"]
          and params["omega_cdm"]==evidence["omega_cdm"]
          and params["N_ur"]==evidence["N_ur"]
          and params["tau_reio"]==evidence["tau_reio"]
          and evidence["states"]==["FD","plus","minus"]
          and evidence["has_custom_psd_files"] is True,
          "original CLASS private source/PSD science fact changed")
    return params

def validate_claim(raw,expected_sha,protocol):
    check(sha(raw)==expected_sha,"original E17D2b2 report SHA mismatch")
    d=json.loads(raw)
    check(d["prospective_protocol_git_blob"]==PROTO_BLOB
          and d["nominal_six_matches"]==6
          and d["real_ASDF_header_or_halo_or_merger_tree_loaded"] is False
          and d["actual_M200c_halo_branch_and_HOD_available_for_physical_model"] is False
          and d["independently_simulated_Fplus_Fminus_wake"] is False
          and d["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and d["observed_odd_SEALED"] is True
          and d["new_CLASS_or_FITS_or_mock_or_ASDF_catalogue_download"] is False,
          "original public-docs-only source audit or absolute physical STOP")
    ours=independent_original(protocol)
    official=protocol["AbacusSummit_c000_official_documented_not_file_metadata"]
    err=0.;matches={}
    for k in protocol["numerical_diagnostic_only"]["compare_exact"]:
        a=d["exact_nominal_parameter_matches"][k]
        x=ours[k];y=official[k]
        check(x==y and a["NOMINALLY_MATCHED"] is True,
              "nominal c000 header-derived versus frozen source mismatch "+k)
        err=max(err,abs(a["original_frozen_source"]-x)/max(1.,abs(x)),
                abs(a["official_public_c000_table"]-y)/max(1.,abs(y)))
        matches[k]=x
    mismatch={}
    for k in protocol["numerical_diagnostic_only"]["compare_mismatch"]:
        a=d["nonzero_original_vs_c000_baseline_mismatches"][k]
        x=ours[k];y=official[k]
        check(x!=y and a["NONZERO_MISMATCH"] is True,
              "source As or tau wrongly equal to official c000")
        difference=x-y
        ratio=x/y-1.
        err=max(err,abs(a["original_frozen_source"]-x)/max(1.,abs(x)),
                abs(a["official_public_c000_table"]-y)/max(1.,abs(y)),
                abs(a["difference_original_minus_c000"]-difference)/max(1.,abs(difference)),
                abs(a["fractional_original_over_c000_minus_one"]-ratio))
        mismatch[k]={"difference_original_minus_c000":difference,
                     "fractional_difference":ratio}
    check(d["CompaSO_L1_SO_mean_epoch_header_NOT_M200c"] is True
          and d["c000_secondary_directory_z095_NOT_verified_ASDF_header_z"] is True
          and d["Abacus_c000_only_smooth_neutrino_Nbody_treatment"] is True
          and d["secondary_PID_only_not_complete_halo_field_particle_RV"] is True
          and d["original_class_mnu0p06_FD_Fplus_Fminus_are_NOT_c000_smooth_neutrino_particle_sim"] is True
          and d["c000_cleaned_tree_product_exists_in_public_docs_but_exact_box_files_not_verified"] is True,
          "original incorrectly promotes Abacus documented candidate to actual physical source")
    check(d["official_source_URLs"]==[v["url"] for v in
        protocol["official_public_source_evidence_as_of_date"]],
          "immutable published source provenance URLs")
    check(err<protocol["preregistered_QA"][
        "independent_stdlib_replay_of_literal_compared_values_scaled_gap_max"],
        "independent decimal AST versus original source metadata QA")
    return err,matches,mismatch

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"independent provenance SHA report would change")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b2i_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D2B2_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    opt=a.parse_args()
    p=gate()
    if opt.self_test:
        s=independent_original(p)
        check(s["h"]==.6736 and s["tau_reio"]==.054 and s["A_s"]>2e-9,
              "independent original CLASS AST and Decimal source preflight")
        print("E17D2B2_INDEPENDENT_SHA_AST_DECIMAL_ORIGINAL_SOURCE_PREFLIGHT_PASS",
              "NO_REAL_ASDF_NO_M200C_NO_NEUTRINO_WAKE",flush=True)
        return
    path=opt.output_dir/"e17d2b2_original_abacus_c000_documented_metadata_bridge.json"
    raw=path.read_bytes()
    if opt.negative_control:
        # SHA deliberately changed, not new halo data or catalogue.
        fake=json.loads(raw);fake["nominal_six_matches"]=7
        bad=(json.dumps(fake,indent=2,allow_nan=False)+"\n").encode()
        try:validate_claim(bad,sha(raw),p)
        except ValueError as e:
            check("original E17D2b2 report SHA mismatch" in str(e),
                  "tampered original report rejected at wrong gate")
            print("E17D2B2_INDEPENDENT_SYNTHETIC_TAMPERED_ORIGINAL_SHA_REJECTED",
                  flush=True)
            return
        raise AssertionError("mutated report SHA was accepted")
    err,matched,mismatch=validate_claim(raw,sha(raw),p)
    result={"date":"2026-09-28",
            "status":"E17D2B2_INDEPENDENT_PURE_STDLIB_AST_DECIMAL_ABACUS_C000_DOCS_PROVENANCE_REPLAY_PASS_PHYSICAL_B_BLOCKED",
            "prospective_protocol_git_blob":PROTO_BLOB,
            "original_source_only_report_sha256":sha(raw),
            "replayed_immutable_original_repo_parents":p["immutable_original_repo_parents"],
            "method":"Separate pure stdlib AST read of exact frozen original CLASS source, Decimal exp primordial amplitude, literal comparison with official documented c000 metadata and original report full SHA. No original worker, CLASS, NumPy/SciPy, data provider API, ASDF, halo trees, observed odd.",
            "six_original_vs_public_c000_nominal_matches":matched,
            "recomputed_nonzero_mismatches":mismatch,
            "max_scaled_independent_provenance_arithmetic_gap":err,
            "true_ABACUS_ASDF_header_z_and_halo_mass_unread":True,
            "cleaned_CosmoHalo_M200c_and_tracer_HOD_not_verified":True,
            "original_Fplus_Fminus_nonlinear_neutrino_halo_response_NOT_simulated":True,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,
            "new_CLASS_FITS_mock_ASDF_catalogue_download":False,
            "main_untouched_PR_draft":True}
    save_once(opt.output_dir/"e17d2b2_independent_ast_decimal_c000_metadata_replay.json",result)
    print("E17D2B2_INDEPENDENT_C000_DOCS_ONLY_PROVENANCE_PASS",
          "SCALED_GAP",err,"HALO_M200C_AND_PHYSICAL_B_BLOCKED",flush=True)
if __name__=="__main__":main()
