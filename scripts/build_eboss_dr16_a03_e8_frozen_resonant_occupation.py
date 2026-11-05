#!/usr/bin/env python3
"""A03E8: physical-normalization-preserving F+/- resonant-source response.

This computes an *absolute CLASS-convention kinetic distribution* and the
dimensionless FD-relative resonant occupation factor in the linear collisionless
halo wake (Okoli et al., MNRAS 468, 2164, 2017, Eq. 14--19). It does NOT compute
a physical LRGxELG correlation, a CLASS transfer, a pair window or an S/N.
Never fit an amplitude, renormalize a template by its maximum, or open FITS.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import io
import itertools
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
import class_response_optimize as cro

KINETIC=ROOT/"code/class_response_optimize.py"
DIRECTION=ROOT/"code/wake_phase7_template.py"
E7_PARENT=ROOT/"source_data/eboss_dr16_mock0001_full_galaxy_4800_48000_completed_source_only_audit_2026-09-27.json"
E6_PROTO=ROOT/"source_data/eboss_dr16_a03_e6_physical_template_to_empirical_window_bridge_protocol_2026-09-27.json"
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e8_frozen_resonant_occupancy_protocol_2026-09-27.json"
SOURCE_PINS={
    KINETIC:"44956e353dbb238537d60da687d1a2a26fa13ee6",
    DIRECTION:"63b291fdbbf1bf76df82583583c22ddf7bdc68a7",
    E7_PARENT:"aa080ec8a3f5faab1731396bed1b048c7dfa991c",
    E6_PROTO:"4df20872cff5b3d8aa09831ab50b22fa5f8b6581",
    PROTOCOL:"2856c190b58fde9b39e2cb3f0d0006a77fe9b70c",
}
MASS_EV=.06
Z_MATCH=1100.
FRAC=.30
Z_GRID=(.90,.95,1.0)
V_GRID_KMS=(100.,200.,400.)
C_KMS=299792.458
OUTPUT_DEFAULT=ROOT/"eboss_workspace/a03_physics_source/e8_frozen_kinetic_resonance"

def need(cond,msg):
    if not cond:raise ValueError(msg)

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def source_gate():
    for path,expected in SOURCE_PINS.items():
        need(git_blob(path.read_bytes())==expected,
             "Original source/preregistered source-only content changed: "+str(path))
    proto=json.loads(PROTOCOL.read_bytes())
    fp=proto["fixed_physics"]
    need(proto["scope"]=="FIXED_SOURCE_ONLY_DISTRIBUTION_AND_CONDITIONAL_RESONANT_FACTOR_NO_CLASS_OR_FITS_REQUIRED" and
         fp["relic_mass_eV"]==MASS_EV and fp["match_z"]==Z_MATCH and
         fp["pointwise_deformation"]==FRAC and
         fp["illustrative_z"]==list(Z_GRID) and
         fp["illustrative_line_of_sight_velocity_kms"]==list(V_GRID_KMS) and
         proto["inviolable_guards"]["use_FD_relative_fraction_as_absolute_galaxy_xi"] is False and
         proto["inviolable_guards"]["observed_odd_read"] is False and
         proto["inviolable_guards"]["borrow_DESI_BGS_literature_Fisher_calibration_as_eBOSS_absolute"] is False,
         "Frozen E8 source-only physics scope changed")
    parent=json.loads(E7_PARENT.read_bytes())
    need(parent["source_report"]["sha256"]==
         "187e21e206cdb6ded1108887a2b6670113d50b752e438d28212cfc3d187572da" and
         parent["decision"]["observed_odd_remains_sealed"] is True,
         "Completed source-only mock parent or observed seal changed")
    return proto

def fixed_coefficients():
    """Read literal prior DESI coefficients without importing CLASS-dependent code."""
    tree=ast.parse(DIRECTION.read_text(encoding="utf-8"))
    candidates=[]
    for node in tree.body:
        if (isinstance(node,ast.Assign) and
            any(isinstance(t,ast.Name) and t.id=="COEFF" for t in node.targets)):
            need(isinstance(node.value,ast.Call) and len(node.value.args)>=1,
                 "Original frozen direction no longer has array constructor")
            candidates.append(np.asarray(ast.literal_eval(node.value.args[0]),dtype="f8"))
    need(len(candidates)==1 and candidates[0].shape==(7,) and
         np.isfinite(candidates[0]).all(),
         "Cannot reconstruct exact published seven-coefficient direction")
    # Detect a source-definition change even if the parser still finds 7 numbers.
    need(np.array_equal(candidates[0],np.asarray([
         0.10329109892348824,-0.060769580477018574,
         0.16915381451528330,0.48937980690857140,
         0.68612921559544580,0.47909899899373720,
         0.13123736993120932],dtype="f8")),
         "Published frozen seven-coefficient direction was replaced")
    return candidates[0]

def kinetic_source():
    q=np.linspace(0.,20.,4000,dtype="f8")
    f0,weights,basis,null,shapes,ym,M=cro.kinetic_objects(q,MASS_EV,Z_MATCH)
    coef=fixed_coefficients()
    need(shapes.shape==(7,4000) and
         null.shape[1]==7 and np.all(np.isfinite(shapes)),
         "Original kinetic-moment null-space dimension changed")
    raw=coef@shapes
    norm=float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    need(np.isfinite(norm) and norm>0,"Degenerate frozen null direction")
    direction=raw/norm
    fp=f0+FRAC*direction
    fm=f0-FRAC*direction
    mom=lambda f:cro.moments(f,q,weights)
    mp,mm,m0=mom(fp),mom(fm),mom(f0)
    mismatch=np.abs(mp-mm)/np.maximum(.5*(np.abs(mp)+np.abs(mm)),1e-300)
    fracp=np.abs(fp-f0)/f0
    fracm=np.abs(fm-f0)/f0
    need(np.all(fp>0) and np.all(fm>0) and
         np.max(fracp)<=FRAC*(1+1e-12) and
         np.max(fracm)<=FRAC*(1+1e-12) and
         np.isclose(np.max(fracp),FRAC,rtol=0,atol=2e-13) and
         np.isclose(np.max(fracm),FRAC,rtol=0,atol=2e-13) and
         np.max(mismatch)<1e-12 and
         np.allclose(fp+fm,2*f0,atol=1e-18,rtol=2e-13),
         "Physical normalization, positivity or source-matched moment gate failed")
    need(np.max(np.abs((fp+fm)/(2*f0)-1.))<3e-13,
         "Kinetic pair was not built symmetrically around CLASS-normalized FD")
    # The coefficient basis is SVD-derived; pin two broad source sign checks.
    tq=np.asarray([.125,.25],dtype="f8")
    need(np.all(np.interp(tq,q,fp/f0)<1.) and
         np.all(np.interp(tq,q,fm/f0)>1.),
         "Frozen physically selected null direction changed sign")
    return dict(q=q,f0=f0,fp=fp,fm=fm,normalization=norm,
                moments_0=m0,moments_plus=mp,moments_minus=mm,
                moment_mismatch=mismatch,coef=coef)

def conditional_resonant_ratio(q,f0,fp,fm,mass_ev,z,v_parallel_kms):
    """Exact source ratio conditional on SAME halo v,potential and source evolution.

    From the total-derivative boundary term in Okoli+2017 Eq.(18)--(19):
    FD occupancy at p_parallel=m |v_parallel| becomes isotropic F(p_parallel).
    CLASS f0/f+/f- share the same phase-space prefactor, which cancels in ratio.
    No assertion about modified v, cosmological transfer or galaxy bias follows.
    """
    need(mass_ev>0 and z>=0 and np.isfinite(z) and
         np.isfinite(v_parallel_kms) and v_parallel_kms>=0,
         "Invalid physical resonant momentum inputs")
    tnu0=float(cro.T_NCDM*cro.TCMB_K*cro.KB_EV_K)
    qr=mass_ev*v_parallel_kms/(C_KMS*tnu0*(1+z))
    need(np.isfinite(qr) and float(q[0])<=qr<=float(q[-1]),
         "Resonance lies outside frozen original q-grid; refuse extrapolation")
    a=float(np.interp(qr,q,f0))
    p=float(np.interp(qr,q,fp))
    m=float(np.interp(qr,q,fm))
    need(a>0 and p>0 and m>0,"Interpolated resonant occupancy not positive")
    return {
        "z_illustrative_not_measured_effective_z":float(z),
        "v_parallel_abs_kms_illustrative_not_measured_eBOSS":float(v_parallel_kms),
        "q_res":float(qr),
        "f_plus_over_FD":p/a,
        "f_minus_over_FD":m/a,
        "f_plus_minus_f_minus_over_FD":(p-m)/a,
        "pair_mean_over_FD":(p+m)/(2*a),
        "not_absolute_galaxy_xi":True,
    }

def build_report(src):
    rows=[conditional_resonant_ratio(
          src["q"],src["f0"],src["fp"],src["fm"],MASS_EV,z,v)
          for z,v in itertools.product(Z_GRID,V_GRID_KMS)]
    for row in rows:
        need(abs(row["pair_mean_over_FD"]-1.)<3e-13 and
             max(abs(row["f_plus_over_FD"]-1.),
                 abs(row["f_minus_over_FD"]-1.))<=FRAC+1e-12,
             "Resonant FD-relative source ratio violates exact symmetric 30% cap")
    return {
        "status":"E8_SOURCE_PHYSICAL_FPLUS_FMINUS_RESONANT_COEFFICIENT_READY_GALAXY_AMPLITUDE_NOT_IDENTIFIED",
        "scope":"Physical CLASS-normalized isotropic source F(q), conditional FD-relative halo linear-response resonant factor; no CLASS transfer or galaxy observable",
        "source_protocol_git_blob_sha1":SOURCE_PINS[PROTOCOL],
        "frozen_COEFF":src["coef"].tolist(),
        "CLASS_kinetic_distribution_prefactor":"2/(2*pi)^3",
        "mass_eV":MASS_EV,"source_matching_z":Z_MATCH,
        "maximum_fractional_distortion":FRAC,
        "source_null_direction_max_relative_before_unit_pointwise_normalization":float(src["normalization"]),
        "min_Fplus":float(np.min(src["fp"])),
        "min_Fminus":float(np.min(src["fm"])),
        "maximum_relative_moment_mismatch":float(np.max(src["moment_mismatch"])),
        "FD_moments_n_rho_P":src["moments_0"].tolist(),
        "plus_moments_n_rho_P":src["moments_plus"].tolist(),
        "minus_moments_n_rho_P":src["moments_minus"].tolist(),
        "resonant_source_ratios":rows,
        "resonant_ratio_reference":"Conditional linear-response Eq.(14)-(19) in Okoli et al. MNRAS 468 (2017) 2164-2175; ratio only at equal halo potential, relative line-of-sight speed, source mass and phase-space convention.",
        "redshift_midpoint_z_0p95_is_not_eBOSS_measured_z_eff":True,
        "velocity_values_are_not_eBOSS_measured_relative_velocity":True,
        "raw_CLASS_transfer_and_absolute_wake_phase_not_computed":True,
        "eBOSS_LRG_ELG_absolute_xi_ell_not_computed":True,
        "physical_eBOSS_tracer_coupling_not_calibrated":True,
        "A03_pair_redshift_conditional_window_not_certified":True,
        "A04_eBOSS_covariance_not_available":True,
        "observed_galaxy_rows_read":False,
        "observed_random_rows_read":False,
        "observed_odd_data_vector_read":False,
        "science_acceptance_or_detection":False,
        "no_download_or_new_seeds_or_cuts":True,
    }

def synthetic_self_test():
    proto=source_gate()
    src=kinetic_source()
    r0=conditional_resonant_ratio(src["q"],src["f0"],src["fp"],src["fm"],
                                   MASS_EV,.95,0.)
    need(abs(r0["pair_mean_over_FD"]-1.)<2e-13,
         "Synthetic zero-speed relative normalization failed")
    r=conditional_resonant_ratio(src["q"],src["f0"],src["fp"],src["fm"],
                                 MASS_EV,.95,200.)
    flipped=conditional_resonant_ratio(src["q"],src["f0"],src["fm"],src["fp"],
                                       MASS_EV,.95,200.)
    need(abs(r["f_plus_over_FD"]-flipped["f_minus_over_FD"])<2e-13 and
         abs(r["f_minus_over_FD"]-flipped["f_plus_over_FD"])<2e-13 and
         abs(r["f_plus_minus_f_minus_over_FD"]+
             flipped["f_plus_minus_f_minus_over_FD"])<2e-13,
         "Source sign exchange negative control failed")
    try:
        conditional_resonant_ratio(src["q"],src["f0"],src["fp"],src["fm"],
                                   MASS_EV,.95,1e9)
    except ValueError:pass
    else:raise AssertionError("Source extrapolation was silently accepted")
    need(proto["not_deliverables"][0]=="Absolutely calibrated eBOSS ξ1 or ξ3" and
         build_report(src)["eBOSS_LRG_ELG_absolute_xi_ell_not_computed"] is True,
         "Source-only E8 accidentally claims an absolute eBOSS galaxy signal")
    print("E8_ORIGINAL_FROZEN_FPLUS_FMINUS_CLASS_NORMALIZATION_AND_MOMENT_CLOSURE_OK",
          "MAX_REL_MOMENT_MISMATCH",float(np.max(src["moment_mismatch"])),flush=True)
    print("E8_FD_RELATIVE_CONDITIONAL_RESONANT_SIGN_SWAP_AND_NO_EXTRAPOLATION_OK",
          "NO_CLASS NO_FITS NO_OBSERVED",flush=True)
    return src

def atomic_reproducible(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        need(path.read_bytes()==raw,
             "Existing frozen source output differs; preserve original bytes")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e8_",delete=False) as handle:
        tmp=Path(handle.name);handle.write(raw);handle.flush();os.fsync(handle.fileno())
    try:
        need(not path.exists(),"Concurrent frozen source output appeared")
        os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUTPUT_DEFAULT)
    args=a.parse_args()
    src=synthetic_self_test()
    if args.self_test:return 0
    table=np.column_stack([src[k] for k in ("q","f0","fp","fm")])
    f=io.StringIO()
    np.savetxt(f,table,fmt="%.17e",delimiter=",",
               header="q_dimensionless,F0_CLASS_normalized,Fplus_CLASS_normalized,Fminus_CLASS_normalized",
               comments="")
    table_bytes=f.getvalue().encode("utf-8")
    report=build_report(src)
    report["frozen_source_csv_sha256"]=sha(table_bytes)
    report["frozen_source_csv_rows"]=4000
    report_bytes=(json.dumps(report,indent=2,allow_nan=False)+"\n").encode("utf-8")
    out=args.output_dir
    atomic_reproducible(out/"frozen_matched_distributions_4000q.csv",table_bytes)
    atomic_reproducible(out/"conditional_resonant_source_response.json",report_bytes)
    print("E8_SOURCE_ONLY_CONDITIONAL_RESONANT_FD_RELATIVE_RESPONSE_SAVED",
          "FROZEN_CSV_SHA256",report["frozen_source_csv_sha256"],
          "REPORT",out/"conditional_resonant_source_response.json",flush=True)
    for row in report["resonant_source_ratios"]:
        print("E8_ILLUSTRATIVE_CONDITIONAL",row["z_illustrative_not_measured_effective_z"],
              row["v_parallel_abs_kms_illustrative_not_measured_eBOSS"],
              row["q_res"],row["f_plus_minus_f_minus_over_FD"],flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
