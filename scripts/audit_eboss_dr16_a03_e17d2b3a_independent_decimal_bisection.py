#!/usr/bin/env python3
"""Independent pure-stdlib Decimal+bisection E17D2b3a M200c counterexample.

No original worker imports; no Abacus ASDF, catalog, halo trees or particles,
no CLASS/NumPy/SciPy, no observed eBOSS/odd. All values are synthetic TOY.
"""
from __future__ import annotations
import argparse
from decimal import Decimal,localcontext
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3a_mass_reference_nonidentifiability_prereg_2026-09-28.json"
PRE_BLOB="81554b0f7cdefa599f302d828ab3529812bb3a59"
B2D=ROOT/"source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28"
B2=B2D/"e17d2b2_original_abacus_c000_documented_metadata_bridge.json"
B2I=B2D/"e17d2b2_independent_ast_decimal_c000_metadata_replay.json"
B2MAN=B2D/"archive_manifest.json"
B2PROTO=ROOT/"source_data/eboss_dr16_a03_e17d2b2_abacussummit_c000_metadata_bridge_prereg_2026-09-28.json"
B2RUNNER=ROOT/"scripts/audit_eboss_dr16_a03_e17d2b2_abacus_c000_metadata_bridge.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
B1=ROOT/"source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28/e17d2b1_original_joint_conditional_exterior_assembly_rates.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3a_l1_vs_200c"

def check(ok,why):
    if not ok:raise ValueError("E17D2B3A_INDEPENDENT_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"prospective E17D2b3a protocol blob")
    p=json.loads(PRE.read_bytes())
    pa=p["original_pinned"]
    for path,key in ((B2MAN,"E17D2b2_archive_manifest_blob"),
                     (B2PROTO,"E17D2b2_protocol_blob"),
                     (B2RUNNER,"E17D2b2_original_runner_blob")):
        check(blob(path.read_bytes())==pa[key],"prior immutable Git source: "+key)
    for path,key in ((E8,"E8_csv_sha256"),(E16,"E16_full_sha256"),
                     (E17A,"E17A_joint_sha256"),(B1,"E17D2b1_joint_sha256"),
                     (B2,"E17D2b2_original_sha256"),
                     (B2I,"E17D2b2_independent_sha256")):
        check(sha(path.read_bytes())==pa[key],"immutable science original: "+key)
    j=json.loads(B2.read_bytes())
    check(j["real_ASDF_header_or_halo_or_merger_tree_loaded"] is False
          and j["CompaSO_L1_SO_mean_epoch_header_NOT_M200c"] is True
          and j["c000_secondary_directory_z095_NOT_verified_ASDF_header_z"] is True
          and j["observed_odd_SEALED"] is True
          and j["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and all(p["absolute_STOP"].values()),
          "prior no-read and physical B hard STOP changed")
    s=p["synthetic_scope"]
    check(s["h_c000_documented"]==.6736
          and s["omega_b_documented"]==.02237
          and s["omega_cdm_documented"]==.12
          and s["omega_ncdm_documented"]==.0006442
          and s["z_illustrative_ONLY_not_real_header"]==.95
          and s["delta_L1_mean_illustrative_ONLY_not_real_header"]==200
          and s["positive_enclosed_mass_exponents"]==[1,2]
          and s["reference_radius_R_L1_normalized"]==1
          and s["reference_mass_M_L1_normalized"]==1,
          "prospective mathematical toy input changed")
    return p

def independently_evaluate():
    # Construct an exact Decimal route, not original worker's binary floats.
    with localcontext() as c:
        c.prec=55
        d=Decimal
        h=d("0.6736")
        wb=d("0.02237")
        wc=d("0.1200")
        wn=d("0.00064420")
        redshift=d("0.95")
        threshold=d("200")
        om=(wb+wc+wn)/(h*h)
        E2=om*(1+redshift)**3+1-om
        omz=om*(1+redshift)**3/E2
        target=d("200")/omz
        check(d(0)<om<d(1) and d(0)<omz<d(1)
              and target>threshold,
              "independent positive flat matter-Lambda synthetic domain")
        ratio=threshold/target
        # Exact algebra from mean density equation.
        roots={1:ratio.sqrt(),2:ratio}
        root_bisection={}
        for pp in (1,2):
            lo=d(0);hi=d(1)
            for _ in range(200):
                m=(lo+hi)/2
                rho_over_mean=threshold*m**(pp-3)
                if rho_over_mean>target:lo=m
                else:hi=m
            root_bisection[pp]=(lo+hi)/2
            check(abs(roots[pp]-root_bisection[pp])<d("1e-48"),
                  "separately solved numerical root differs from closed form")
        rows=[]
        for pp in (1,2):
            x=roots[pp]
            m=x**pp
            reclose=threshold*m/(x**3*target)
            check(abs(reclose-1)<d("1e-48"),
                  "independent strict 200-critical mean density closure")
            rows.append({
                "power":pp,
                "root":float(x),
                "mass_fraction":float(m),
                "bisection_root":float(root_bisection[pp]),
                "exact_root_minus_bisection_abs":float(abs(x-root_bisection[pp])),
                "200critical_reclosure_abs":float(abs(reclose-1))
            })
        result={
            "omega_m0_flat_matterLambda_TOY":float(om),
            "E2_z_flat_matterLambda_TOY":float(E2),
            "omega_m_z_flat_matterLambda_TOY":float(omz),
            "synthetic_200crit_over_total_matter_mean_threshold":float(target),
            "synthetic_L1_SO_threshold_over_total_matter_mean_NOT_ABACUS_HEADER":float(threshold),
            "synthetic_ratio_threshold_200c_to_L1":float(target/threshold)}
        gap=float(abs(roots[1]**1-roots[2]**2))
        return result,rows,gap

def replay(raw,expected,p):
    check(sha(raw)==expected,"original E17D2b3a SHA was mutated")
    original=json.loads(raw)
    check(original["prospective_protocol_git_blob"]==PRE_BLOB
          and original["actual_Abacus_SODensityL1_header_read"] is False
          and original["actual_Abacus_ASDF_redshift_header_read"] is False
          and original["actual_Abacus_halo_particle_or_merger_tree_read"] is False
          and original["any_real_halo_M200c_remeasured"] is False
          and original["observed_odd_SEALED"] is True
          and original["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and original["new_CLASS_ASDF_FITS_halo_catalogue_mock_download"] is False,
          "original conditional-only science scope or no-read STOP drift")
    independent,profiles,difference=independently_evaluate()
    tolerance=p["predeclared_checks"]["independent_stdlib_Decimal_bisection_scaled_tolerance"]
    err=0.
    for key,new in independent.items():
        old=original["source_only_analytic_flat_matterLambda_TOY"][key]
        err=max(err,abs(new-old)/max(1.,abs(new),abs(old)))
    for v,old in zip(profiles,original["two_predeclared_positive_normalized_radial_profiles"]):
        check(v["power"]==old["enclosed_mass_exponent_synthetic_ONLY"],
              "positive mass-profile family changed")
        for ind,original_field in (("root","r200c_over_shared_R_L1_TOY"),
                                   ("mass_fraction","M200c_over_shared_M_L1_TOY"),
                                   ("200critical_reclosure_abs",
                                    "200critical_mean_density_reclosure_relative")):
            err=max(err,abs(v[ind]-old[original_field])/
                    max(1.,abs(v[ind]),abs(old[original_field])))
        check(v["exact_root_minus_bisection_abs"]<tolerance,
              "independent Decimal bisection root fails")
    err=max(err,abs(difference-original[
        "distinct_M200c_mass_fraction_absolute_gap_synthetic_only"]))
    check(len(original["synthetic_header_directory_z_and_L1_mass_negative_controls"])==3
          and original["same_L1_mass_radius_density_and_synthetic_z"] is True
          and original["idealized_SO_crossing_is_stronger_than_real_SO_radius_metadata"] is True
          and difference>.05 and err<tolerance,
          "independent full toy analytic/bisection versus original source QA failed")
    return err,independent,profiles,difference

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,
                           "independent immutable report collision")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,
            prefix=".e17d2b3ai_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D2B3A_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--negative-control",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    opt=ap.parse_args()
    p=gate()
    if opt.self_test:
        values,profiles,difference=independently_evaluate()
        check(len(profiles)==2 and difference>.05,
              "synthetic profile ambiguity must be genuinely nonzero")
        print("E17D2B3A_INDEPENDENT_DECIMAL_BISECTION_ORIGINAL_SHA_PREFLIGHT_PASS",
              "M200C_NORMALIZED_DIFFERENCE_TOY",difference,
              "NO_ACTUAL_HALO_READ",flush=True)
        return
    path=opt.output_dir/"e17d2b3a_original_l1_200critical_nonidentifiability.json"
    raw=path.read_bytes()
    if opt.negative_control:
        invalid=json.loads(raw)
        invalid["actual_Abacus_ASDF_redshift_header_read"]=True
        tampered=(json.dumps(invalid,indent=2,allow_nan=False)+"\n").encode()
        try:replay(tampered,sha(raw),p)
        except ValueError as e:
            check("original E17D2b3a SHA was mutated" in str(e),
                  "tamper original SHA must be rejected before physical content")
            print("E17D2B3A_INDEPENDENT_TAMPERED_ORIGINAL_SHA_REJECTED",flush=True)
            return
        raise AssertionError("mutated original report was accepted")
    err,background,profiles,difference=replay(raw,sha(raw),p)
    result={"date":"2026-09-28",
            "status":"E17D2B3A_INDEPENDENT_DECIMAL_AND_EXACT_BISECTION_TWO_POSITIVE_SAME_L1_DIFFERENT_M200C_PASS_ACTUAL_HALO_UNREAD",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_full_report_sha256":sha(raw),
            "independent_decimal_synthetic_flat_matterLambda":background,
            "independent_two_exact_200critical_roots_and_mass_fractions":profiles,
            "independent_mass_fraction_nonidentifiability_gap_synthetic_only":difference,
            "max_scaled_independent_original_analytic_vs_decimal_gap":err,
            "method":"Pure stdlib Decimal exact algebra and independent 200-step density-threshold root bisection, SHA-gated original science parents, fail closed no real halo data. NO original E17D2b3a runner import, CLASS, NumPy, ASDF, actual Abacus mass/history or observed odd.",
            "real_Abacus_SODensityL1_or_actual_z_header_loaded":False,
            "real_Abacus_M200c_200critical_mass_remeasured":False,
            "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,
            "new_CLASS_ASDF_FITS_halo_catalogue_mock_download":False,
            "main_untouched_PR_draft":True}
    save_once(opt.output_dir/"e17d2b3a_independent_decimal_bisection_200critical_replay.json",
              result)
    print("E17D2B3A_INDEPENDENT_TWO_POSITIVE_SAME_L1_DIFFERENT_M200C_PASS",
          "MASS_GAP",difference,"MAX_SCALED",err,
          "REAL_HALO_M200C_B_BLOCKED",flush=True)
if __name__=="__main__":main()
