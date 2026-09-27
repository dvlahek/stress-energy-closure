#!/usr/bin/env python3
"""A03E8c: source-only conditional k-mode phase and published Eq20 conventions.

POST-E8 retrospective physics audit. Zero FITS, CLASS, random/mock data,
observed odd vector, new physical fit, galaxy xi or S/N. Computes the ratio of
the STATIC Eq20 wake field and the unperturbed same-mode halo force. Neither
this instantaneous phase nor the published mu=1 proxy is an eBOSS prediction.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE=ROOT/"source_data/eboss_dr16_a03_e8_frozen_kinetic_resonance_and_static_SI_green_result_2026-09-27.json"
ARCHIVE_BLOB="b2978e74384356db72dbb3af2e93b55cc3e16b38"
E8_PROTO=ROOT/"source_data/eboss_dr16_a03_e8_frozen_resonant_occupancy_protocol_2026-09-27.json"
E8_PROTO_BLOB="2856c190b58fde9b39e2cb3f0d0006a77fe9b70c"
E8B_PROTO=ROOT/"source_data/eboss_dr16_a03_e8b_static_halo_green_kernel_protocol_2026-09-27.json"
E8B_PROTO_BLOB="7865f059d312a0a97753f7c13b6faa5bb07a0e16"
ERRATUM=ROOT/"source_data/eboss_dr16_a03_e8c_exact_eq20_phase_and_bibliography_erratum_2026-09-27.json"
ERRATUM_BLOB="549dccc2b923800f06d0081741ea7c80f3181e4b"
OUT=ROOT/"eboss_workspace/a03_physics_source/e8c_exact_conditional_mode_phase_source_only.json"

def require(value, message):
    if not value: raise ValueError(message)

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\x00"+raw).hexdigest()

def pinned_json(path,blob):
    raw=path.read_bytes()
    require(git_blob(raw)==blob,"Frozen source Git blob changed: "+str(path))
    return json.loads(raw)

def close(x,y,rtol=3e-13):
    return math.isclose(float(x),float(y),rel_tol=rtol,abs_tol=0.)

def source():
    a=pinned_json(ARCHIVE,ARCHIVE_BLOB)
    p=pinned_json(E8_PROTO,E8_PROTO_BLOB)
    pb=pinned_json(E8B_PROTO,E8B_PROTO_BLOB)
    e=pinned_json(ERRATUM,ERRATUM_BLOB)
    require(a["guards"]["observed_odd_read"] is False and
            a["guards"]["inference_or_detection"] is False and
            a["source_provenance"]["E8b_protocol_git_blob"]==E8B_PROTO_BLOB and
            a["source_provenance"]["frozen_distribution_4000q_csv_sha256"]==
            "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0",
            "Original source archive / observed seal changed")
    require(p["inviolable_guards"]["observed_odd_read"] is False and
            p["paper_source"].get("doi")=="10.1093/mnras/stx539" and
            e["source_paper"]["verified_DOI"]=="10.1093/mnras/stx560" and
            e["source_paper"]["originals_modified"] is False and
            e["scientific_gate"]["absolute_halo_SI_normalization_vs_published_mu_1_forecast_UNRESOLVED"] is True and
            pb["physical_reference"]["physical_not_comoving_k"] is True and
            pb["guards"]["observed_odd_read"] is False,
            "DOI or published exact FD versus mu normalization erratum changed")
    return a,pb,e

def mode(a,pb,*,z,v_parallel_kms,k_h_per_Mpc,k_kind):
    """The physical Eq20 first-equality ratio, with an explicit k convention."""
    require(k_kind in ("physical","comoving"),
            "Must explicitly identify k as physical or comoving")
    require(math.isfinite(z) and z>=0 and math.isfinite(v_parallel_kms)
            and math.isfinite(k_h_per_Mpc) and k_h_per_Mpc>0,
            "Invalid source-only physical input")
    ref=pb["physical_reference"]
    require(ref["one_ncdm_species"]==1 and
            close(ref["mass_eV"],a["fixed_physics"]["mass_eV"]) and
            close(ref["h"],0.6736),
            "Original E8b species/mass/geometry source convention changed")
    a_scale=1./(1.+z)
    kcom_h=k_h_per_Mpc*(a_scale if k_kind=="physical" else 1.)
    kphys_h=k_h_per_Mpc*(1. if k_kind=="physical" else 1.+z)
    kphys=kphys_h*ref["h"]/ref["Mpc_m"]
    kcom=kcom_h*ref["h"]/ref["Mpc_m"]
    Mkg=ref["mass_eV"]*ref["rest_mass_conversion_eVc2_to_kg"]
    v=v_parallel_kms*1e3
    tnu0=0.71611*2.7255*8.617333262e-5
    q=ref["mass_eV"]*abs(v)/(ref["c_m_per_s"]*tnu0*(1.+z))
    require(q>=0 and q<=20,"Resonant q outside source grid")
    fd=1./(1.+math.exp(q))
    # Only the original z=.95, v=200 physical response has pinned F±
    # occupancy. Do not invent interpolated new F± source responses here.
    orig=a["example_conditional"]
    require(close(z,orig["z_midpoint_not_measured_effective_z"],rtol=0.) and
            close(v_parallel_kms,orig["v_parallel_abs_illustrative_kms"],rtol=0.),
            "Only original frozen illustrative z and velocity are source-pinned")
    ratio_plus=orig["Fplus_over_FD"]
    ratio_minus=orig["Fminus_over_FD"]
    K=2.*ref["one_ncdm_species"]*ref["Newton_G_m3_kg_s2"]*Mkg**4*v*fd/(ref["hbar_Js"]**3*kphys)
    phases={"FD":K/kphys,"Fplus":K*ratio_plus/kphys,
            "Fminus":K*ratio_minus/kphys}
    phases["Fplus_minus_Fminus"]=phases["Fplus"]-phases["Fminus"]
    com_formula=2.*ref["one_ncdm_species"]*ref["Newton_G_m3_kg_s2"]*Mkg**4*(a_scale*a_scale)*v*fd/(ref["hbar_Js"]**3*kcom*kcom)
    require(close(phases["FD"],com_formula),
            "Exact physical k and published a^2/k_com^2 algebra disagree")
    require(close(phases["Fplus_minus_Fminus"]/phases["FD"],
                  orig["Fplus_minus_Fminus_over_FD"]),
            "Original fixed matched source source-ratio changed")
    return {"z_illustrative":z,"v_parallel_kms_illustrative":v_parallel_kms,
            "input_k_kind":k_kind,"input_k_h_per_Mpc":k_h_per_Mpc,
            "physical_k_h_per_Mpc":kphys_h,
            "comoving_k_h_per_Mpc":kcom_h,
            "k_phys_SI_per_m":kphys,
            "FD_exact_occupation":fd,
            "alpha_signed_conditional_mode_per_unit_original_halo_field":phases,
            "published_mu_equal_1_over_exact_FD_prefactor":1./fd,
            "original_analytic_exact_eq20_not_published_mu1_forecast":True,
            "time_derivative_phase_computed":False,
            "galaxy_xi_computed":False}

def selftest():
    a,pb,e=source()
    orig=a["example_conditional"]
    fixed=mode(a,pb,z=.95,v_parallel_kms=200.,k_h_per_Mpc=.05,
               k_kind="physical")
    vals=fixed["alpha_signed_conditional_mode_per_unit_original_halo_field"]
    for name,old in orig["absolute_conditional_static_gravitational_green_per_unit_Phi_SI_per_m"].items():
        expected={"FD":"FD","Fplus":"Fplus","Fminus":"Fminus",
                  "Fplus_minus_Fminus":"Fplus_minus_Fminus"}[name]
        require(close(vals[expected]*fixed["k_phys_SI_per_m"],old,rtol=8e-13),
                "Exact source-only E8b SI reference K or alpha changed: "+name)
    flipped=mode(a,pb,z=.95,v_parallel_kms=-200.,k_h_per_Mpc=.05,
                 k_kind="physical")
    for term in vals:
        require(close(flipped["alpha_signed_conditional_mode_per_unit_original_halo_field"][term],
                      -vals[term]),"Velocity reversal lost odd parity: "+term)
    same_numeric_comoving=mode(a,pb,z=.95,v_parallel_kms=200.,
                               k_h_per_Mpc=.05,k_kind="comoving")
    for term in vals:
        require(close(same_numeric_comoving["alpha_signed_conditional_mode_per_unit_original_halo_field"][term],
                      vals[term]/(1.+.95)**2),
                "Comoving/physical k conversion missing exact (1+z)^-2 factor")
    for invalid in ("unlabelled","fiducial",""):
        try:mode(a,pb,z=.95,v_parallel_kms=200.,k_h_per_Mpc=.05,k_kind=invalid)
        except ValueError:pass
        else:raise AssertionError("Unknown physical/comoving k accepted")
    # A printed mu in [0.7,1) cannot equal exact FD occupation <=0.5 for q>=0.
    require(fixed["FD_exact_occupation"]<=.5 and
            fixed["published_mu_equal_1_over_exact_FD_prefactor"]>2.,
            "Printed Eq20 exact occupation / approximate mu mismatch was silently removed")
    print("EBOSS_A03_E8C_FROZEN_EQ20_EXACT_SOURCE_PHASE_SI_CLOSURE_OK",
          "ALPHA_FD",vals["FD"],"DELTA_ALPHA",vals["Fplus_minus_Fminus"],flush=True)
    print("EBOSS_A03_E8C_PHYSICAL_COMOVING_K_PARITY_DOI_MU_GUARDS_OK",
          "MU1_OVER_EXACT_FD",fixed["published_mu_equal_1_over_exact_FD_prefactor"],
          "COMOVING_SAME_NUMERIC_FD",
          same_numeric_comoving["alpha_signed_conditional_mode_per_unit_original_halo_field"]["FD"],flush=True)
    return a,pb,e,fixed,same_numeric_comoving

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output",type=Path,default=OUT)
    args=ap.parse_args()
    a,pb,e,physical,comoving=selftest()
    if args.self_test:return 0
    result={
        "status":"POSTHOC_E8C_EXACT_STATIC_MODE_PHASE_READY_MU_CONVENTION_UNRESOLVED_GALAXY_A03_STOP",
        "scope":"Conditional one-mode signed neutrino/halo gravitational-force ratio computed strictly from frozen exact FD Eq20 branch. Not the galaxy phase time derivative, LRGxELG cross xi, physical window, mock covariance or inferential likelihood.",
        "parent_E8_archive_git_blob":ARCHIVE_BLOB,
        "parent_E8_original_protocol_git_blob":E8_PROTO_BLOB,
        "parent_E8b_original_protocol_git_blob":E8B_PROTO_BLOB,
        "E8c_bibliographic_and_mu_convention_erratum_git_blob":ERRATUM_BLOB,
        "DOI_verified":"10.1093/mnras/stx560",
        "original_E8_protocol_DOI_typo_unchanged":"10.1093/mnras/stx539",
        "published_exact_eq20_mu_approximation_normalization_UNRESOLVED":True,
        "physical_k_original_reference":physical,
        "same_numeric_k_reinterpreted_explicitly_as_comoving":comoving,
        "k_relabeling_changes_alpha_by":1./(1.+.95)**2,
        "source_relative_Fplus_minus_Fminus_over_FD_frozen":
            a["example_conditional"]["Fplus_minus_Fminus_over_FD"],
        "absolute_LRG_ELG_xi1_xi3_computed":False,
        "evolving_CLASS_time_derivative_computed":False,
        "physical_eBOSS_pair_window_calibrated":False,
        "independent_eBOSS_A04_covariance_available":False,
        "observed_galaxies_or_randoms_or_odd_read":False,
        "new_download_or_seeds_or_cuts":False
    }
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    if args.output.exists():
        require(args.output.read_bytes()==raw,
                "Existing output changed, refuse overwrite")
    else:
        with tempfile.NamedTemporaryFile(dir=args.output.parent,
                prefix=".e8c_",delete=False) as f:
            tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(tmp,args.output)
        finally:tmp.unlink(missing_ok=True)
    print("EBOSS_A03_E8C_SOURCE_ONLY_CONDITIONAL_MODE_REPORT_SAVED",
          args.output,"SHA256",hashlib.sha256(raw).hexdigest(),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
