#!/usr/bin/env python3
"""A-03E3 preregistered MOCK-ONLY finite-1200 random-split and DD odd-injection QA.

No new files are downloaded. Every one of the 72 independently SHA-pinned
complete realistic EZmock gzip sources is rehashed BEFORE any mock FITS row.
The unchanged 18-case A02 sample/pair and E2 full-density fingerprints must
replay. Random subsets and pair-level +/-0.02 synthetic odd injections are
strictly DESCRIPTIVE, not physical wake injections or inferential covariance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import audit_eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress as E2
from audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot import (
    sample_catalogue, oriented_pair_terms,
)
from audit_eboss_dr16_ezmock0001_galaxy_odd_projection import (
    project, integration_weights,
)
from audit_eboss_dr16_rr_pair_closure import mirrored_closure
from check_eboss_cross_ls_synthetic import cross_landy_szalay, LABELS
from eboss_dr16_fiducial import PRIMARY_GEOMETRY, comoving_mpc_over_h

ROOT=A02.ROOT
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e3_nested_mock_random_split_synthetic_dd_injection_protocol_2026-09-26.json"
E2_ARCHIVE=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json"
E2_MANIFEST=ROOT/"source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_uploaded_manifest_2026-09-26.json"
A02_ARCHIVE=E2.PARENT_REPORT
E2_SOURCE=ROOT/"scripts/audit_eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress.py"
A02_SOURCE=ROOT/"scripts/audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport.py"
SAMPLER=ROOT/"scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py"
ODD_HELPER=ROOT/"scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py"
PASS="A03E3_NINE_MOCK_NESTED_RANDOM_AND_SYNTHETIC_DD_ODD_DESCRIPTIVE_ONLY"
STOP="A03E3_NINE_MOCK_NESTED_RANDOM_AND_SYNTHETIC_DD_ODD_INCOMPLETE_STOP"
LEVELS=("half_A_600","half_B_600","three_quarter_A_900","full_1200")
SCENARIOS=("baseline_original","plus_5pct_ELG_R_z_ramp","minus_5pct_ELG_R_z_ramp")
INJECTIONS=(0.02,-0.02)
ELLS=(0,1,2,3)
SEED_ROOT=20402026
PER_TRACER=("eBOSS_LRG","eBOSS_ELG")

def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def arr_sha(arr):
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()

def catalogue_sha(cat):
    h=hashlib.sha256()
    for v in cat:
        h.update(np.ascontiguousarray(v).tobytes())
    return h.hexdigest()

def exact_gate(p):
    if (p["registration_before"]!="Any first A03E3 mock FITS header or row, any A03E3 random subset or odd injection result"
        or p["user_decision"]!="EMPIRICAL_ONLY_DO_NOT_CONTACT_AUTHORS"
        or p["author_contact_permitted"] is not False
        or p["main_mutation_permitted"] is not False
        or p["observed_odd_data_vector_read"] is not False
        or p["observed_galaxy_rows_read"] is not False
        or p["observed_random_rows_read"] is not False
        or p["physical_empirical_pair_window_certified"] is not False
        or p["inferential_18D_covariance_computed"] is not False
        or p["number_of_cases"]!=18
        or p["fixed_nine_ids"]!=list(A02.IDS)
        or p["caps"]!=list(A02.CAPS)
        or p["tracers"]!=list(A02.TRACERS)
        or p["roles"]!=list(A02.ROLES)
        or p["nested_random_design"]["levels_in_order"]!=list(LEVELS)
        or p["nested_random_design"]["permutation_seed_formula"]!=
           "20402026 + 100000*mock_ID + 1000*cap_index + 100*tracer_index, cap_index NGC=0 SGC=1, tracer_index LRG=0 ELG=1"
        or p["E2_weight_stress_replay"]["scenarios"]!=list(SCENARIOS)
        or p["E2_weight_stress_replay"]["levels_stressed"]!=
           ["half_A_600","half_B_600","full_1200"]
        or p["synthetic_pair_histogram_injection"]["amplitudes"]!=list(INJECTIONS)
        or p["synthetic_pair_histogram_injection"]["numerical_increment_residual_absolute_threshold"]!=1e-11
        or p["parent_A02_implementation_git_blob"]!=E2.blob(A02_SOURCE)
        or p["parent_0001_sampler_git_blob"]!=E2.blob(SAMPLER)
        or p["parent_A03E2_runner_git_blob"]!=E2.blob(E2_SOURCE)
        or p["parent_odd_projection_helper_git_blob"]!=E2.blob(ODD_HELPER)
        or E2.blob(E2_MANIFEST)!=p["parent_A03E2_manifest_git_blob"]):
        raise ValueError("Prospectively frozen E3 formulas, source implementation or no-observed scope changed")
    decision=A02.load(ROOT/p["archived_user_decision"])
    if decision["email_sent"] is not False or decision["author_contact_approved"] is not False:
        raise ValueError("User's frozen empirical-only no-author-contact decision changed")
    e2p=A02.load(E2.PROTOCOL)
    a02p,a02=E2.guard_protocol(e2p)
    e2raw=E2_ARCHIVE.read_bytes()
    if (len(e2raw)!=p["parent_A03E2_report_bytes"]
        or file_sha(E2_ARCHIVE)!=p["parent_A03E2_report_sha256"]
        or E2.blob(E2_ARCHIVE)!=p["parent_A03E2_report_git_blob"]
        or file_sha(A02_ARCHIVE)!=p["parent_A02_full_sha256"]):
        raise ValueError("Exact parent A02/E2 full uploaded report bytes changed")
    e2=json.loads(e2raw)
    manifest=A02.load(E2_MANIFEST)
    if (e2["status"]!=E2.PASS or e2["protocol_sha256"]!=file_sha(E2.PROTOCOL)
        or e2["completed_cases"]!=18 or e2["failed_cases"]!=0
        or e2["errors"]!=[] or e2["observed_odd_data_vector_read"] is not False
        or e2["all_72_original_gzip_full_sha256_rechecked_before_any_E2_FITS_row"] is not True
        or manifest["exact_original_uploaded_SHA256"]!=p["parent_A03E2_report_sha256"]
        or manifest["n_scenario_cases"]!=54
        or list(e2["cases"])!=list(a02["cases"])
        or a02p["source_input_total_count"]!=72):
        raise ValueError("Parent 18x3 original E2 result not closed or immutable")
    return a02p,a02,e2,e2p

def subset_indices(mid,cap,tracer):
    if mid not in A02.IDS or cap not in A02.CAPS or tracer not in PER_TRACER:
        raise ValueError("Unregistered mock ID/cap/tracer in nested subsets")
    seed=SEED_ROOT+100000*mid+1000*A02.CAPS.index(cap)+100*PER_TRACER.index(tracer)
    perm=np.random.default_rng(seed).permutation(1200)
    index={
        "half_A_600":np.sort(perm[:600]),
        "half_B_600":np.sort(perm[600:1200]),
        "three_quarter_A_900":np.sort(perm[:900]),
        "full_1200":np.arange(1200,dtype=np.int64),
    }
    a,b,nine,full=(index[k] for k in LEVELS)
    if (not np.array_equal(np.union1d(a,b),full)
        or np.intersect1d(a,b).size!=0
        or not np.all(np.isin(a,nine))
        or any(len(index[k])!=v for k,v in zip(LEVELS,(600,600,900,1200)))):
        raise ValueError("Preregistered half-complement or nested 900 random subset failed")
    return seed,index

def random_subcatalogue(original,indices):
    if len(original)!=4 or any(v.shape!=(1200,) for v in original):
        raise ValueError("Original A02 R must contain exactly 1200 selected rows per vector")
    if (np.any(indices<0) or np.any(indices>=1200)
        or np.unique(indices).size!=len(indices)
        or not np.array_equal(indices,np.sort(indices))):
        raise ValueError("Random subset is out of range, duplicated or changed original row order")
    return tuple(np.asarray(v[indices],dtype="f8").copy() for v in original)

def elg_random_stress(cat,amplitude):
    if amplitude not in (0.,.05,-.05):
        raise ValueError("Only originally frozen E2 +/-5% ELG RANDOM stress allowed")
    ra,dec,z,w=cat
    if (len(z) not in (600,1200) or any(len(v)!=len(z) for v in cat)
        or not np.all((z>=.9)&(z<1.0)) or not np.all(w>0)):
        raise ValueError("ELG RANDOM subset not one of frozen E3 stress levels")
    u=2*(z-.9)/.1-1.
    factor=1.+amplitude*u
    if (np.any(factor<.95-1e-12) or np.any(factor>1.05+1e-12)
        or np.any(factor<=0) or not np.isfinite(factor).all()):
        raise ValueError("ELG RANDOM E2 artificial stress exceeds fixed +/-5%")
    new=tuple(cat[:3])+(np.asarray(w*factor,dtype="f8"),)
    return new,{"factor_min":float(factor.min()),"factor_max":float(factor.max()),
                "selected_original_weight_SHA256":arr_sha(w),
                "modified_selected_weight_SHA256":arr_sha(new[3]),
                "only_ELG_RANDOM_weight_modified":True}

def calculate(cats,rl,re):
    dist=lambda z:comoving_mpc_over_h(z,PRIMARY_GEOMETRY)
    dl,de=cats["eBOSS_LRG","dat"],cats["eBOSS_ELG","dat"]
    f,fn,fmeta=oriented_pair_terms(dl,de,rl,re,distance=dist)
    rev,rn,rmeta=oriented_pair_terms(de,dl,re,rl,distance=dist)
    mirror={}
    for term,other in A02.MAPPING.items():
        chk=mirrored_closure(
            f[term],rev[other],
            {"pair_normalization":fn[term],"accepted_pairs":fmeta[term]["accepted_pairs"]},
            {"pair_normalization":rn[other],"accepted_pairs":rmeta[other]["accepted_pairs"]},
            A02.MU_EDGES)
        if chk["closure_passed"] is not True:
            raise ValueError("True independent 4-term tracer-pair reversal failed: "+term)
        mirror[term]=chk
    xi,support=cross_landy_szalay(f,fn)
    rxi,rsupport=cross_landy_szalay(rev,rn)
    if (xi.shape!=(6,24) or not np.array_equal(support,rsupport[:,::-1])
        or not np.any(support) or not np.isfinite(xi[support]).all()
        or not np.isfinite(rxi[rsupport]).all()):
        raise ValueError("No positive RR support / nonfinite supported cross xi / reverse mask mismatch")
    residue=float(np.max(np.abs(xi[support]-rxi[:,::-1][support])))
    if not np.isfinite(residue) or residue>=1e-8:
        raise ValueError("Signed mu reverse cross xi mismatch in reduced-density scenario")
    cells=int(support.sum())
    full=cells==144
    support_sha=arr_sha(np.asarray(support,dtype=np.uint8))
    index=np.flatnonzero(support).astype("<i8")
    supported_sha=hashlib.sha256(index.tobytes()+
        np.asarray(xi[support],dtype="<f8").tobytes()).hexdigest()
    rec={
        "status":"complete_full_RR_support" if full else "complete_sparse_RR_support_no_full_multipoles",
        "RR_supported_cells":cells,
        "RR_missing_indices":np.argwhere(~support).tolist(),
        "RR_support_mask_sha256":support_sha,
        "max_abs_xi_reverse_on_supported_cells":residue,
        "supported_cell_index_and_xi_sha256":supported_sha,
        "full_6x24_xi_SHA256":arr_sha(xi) if full else None,
        "forward_pair_terms":fmeta,
        "reverse_pair_terms":rmeta,
        "original_physical_domain_multipoles":project(xi,support,A02.MU_EDGES,ELLS) if full else None,
        "no_unsupported_RR_zero_fill_or_projection":True,
    }
    internal={"hist":f,"norm":fn,"rev_hist":rev,"rev_norm":rn,
              "xi":xi,"support":support,"reverse_xi":rxi}
    return rec,internal

def compare_e2_level(parent,scenario,rec,random_diag):
    old=parent["scenarios"][scenario]
    if rec["full_6x24_xi_SHA256"]!=old["forward_xi_grid_sha256"] or rec["RR_supported_cells"]!=144:
        raise ValueError("Full 1200 R E2 original xi/RR fingerprint changed: "+scenario)
    for side in ("forward_pair_terms","reverse_pair_terms"):
        for term in LABELS:
            new,original=rec[side][term],old[side][term]
            if (new["weighted_histogram_SHA256"]!=original["weighted_histogram_SHA256"]
                or new["accepted_pairs"]!=original["accepted_pairs"]
                or not np.isclose(new["independently_normalized_pair_weight"],
                                  original["independently_normalized_pair_weight"],
                                  rtol=1e-12,atol=1e-12)):
                raise ValueError("Full 1200R original E2 8-pair SHA/count/norm changed: "+scenario+"/"+side+"/"+term)
    if (random_diag["modified_selected_weight_SHA256"]!=
        old["ELG_R_weight_stress"]["modified_selected_ELG_R_weight_SHA256"]):
        raise ValueError("Full 1200 original E2 ELG RANDOM weight source SHA changed: "+scenario)

def synthetic_dd_odd_pair_injection(internal,amplitude):
    if amplitude not in INJECTIONS:
        raise ValueError("Nonpreregistered signed synthetic DD odd amplitude")
    f,n,rv,rn=(internal[k] for k in ("hist","norm","rev_hist","rev_norm"))
    xi,support=internal["xi"],internal["support"]
    mu=(A02.MU_EDGES[1:]+A02.MU_EDGES[:-1])/2.
    fac=1.+amplitude*mu[None,:]
    rev_fac=1.-amplitude*mu[None,:]
    if (np.min(fac)<.98-1e-12 or np.max(fac)>1.02+1e-12
        or np.min(rev_fac)<.98-1e-12 or np.max(rev_fac)>1.02+1e-12):
        raise ValueError("Synthetic DD pair modification outside 0.98..1.02")
    inj={**f,"D1D2":np.asarray(f["D1D2"]*fac,dtype="f8")}
    rev_inj={**rv,"D1D2":np.asarray(rv["D1D2"]*rev_fac,dtype="f8")}
    ix,sp=cross_landy_szalay(inj,n)
    rix,rsp=cross_landy_szalay(rev_inj,rn)
    if not np.array_equal(sp,support) or not np.array_equal(sp,rsp[:,::-1]):
        raise ValueError("Synthetic DD-only injection modified original RR support")
    target=np.full((6,24),np.nan,dtype="f8")
    rr=f["R1R2"]/n["R1R2"]
    target[sp]=amplitude*np.broadcast_to(mu,(6,24))[sp]*(f["D1D2"][sp]/n["D1D2"])/rr[sp]
    inc=ix[sp]-xi[sp]
    residual=float(np.max(np.abs(inc-target[sp])))
    rev_resid=float(np.max(np.abs(ix[sp]-rix[:,::-1][sp])))
    rev_inc=rix[:,::-1][sp]-internal["reverse_xi"][:,::-1][sp]
    rev_analytic_residual=float(np.max(np.abs(rev_inc-target[sp])))
    increment_mirror_residual=float(np.max(np.abs(inc-rev_inc)))
    if (residual>=1e-11 or rev_analytic_residual>=1e-11
        or rev_resid>=1e-8 or increment_mirror_residual>=1e-8):
        raise ValueError("Signed DD synthetic bin-center exact estimator increment/reverse closure failed")
    full=bool(np.all(sp))
    projections=None
    if full:
        base=project(xi,sp,A02.MU_EDGES,ELLS)
        got=project(ix,sp,A02.MU_EDGES,ELLS)
        expected={
            str(ell):(target@integration_weights(A02.MU_EDGES,ell)).tolist()
            for ell in ELLS
        }
        increments={}
        for ell in ELLS:
            label=str(ell)
            diff=np.asarray(got[label]["values_by_fixed_s_bin"])-np.asarray(base[label]["values_by_fixed_s_bin"])
            if not np.allclose(diff,expected[label],rtol=0,atol=1e-11):
                raise ValueError("Full 144-cell synthetic DD increment multipole recovery failed")
            increments[label]={
                "delta_injected_minus_baseline_by_s_bin":diff.tolist(),
                "expected_finite_bin_pair_analytic_increment_by_s_bin":expected[label],
                "max_abs_delta":float(np.max(np.abs(diff))),
            }
        projections=increments
    return {
        "amplitude":amplitude,
        "post_binning_DD_only_modification":True,
        "original_DD_pair_count_and_normalization_unchanged":True,
        "forward_injected_DD_hist_sha256":arr_sha(inj["D1D2"]),
        "reverse_opposite_sign_injected_DD_hist_sha256":arr_sha(rev_inj["D1D2"]),
        "positive_RR_supported_cell_count":int(sp.sum()),
        "max_abs_analytic_xi_increment_residual_on_supported_cells":residual,
        "max_abs_reverse_analytic_xi_increment_residual_on_supported_cells":rev_analytic_residual,
        "max_abs_forward_reverse_injected_increment_mirror_residual_on_supported_cells":increment_mirror_residual,
        "max_abs_forward_reverse_injected_xi_residual_on_supported_cells":rev_resid,
        "max_abs_injected_xi_increment_on_supported_cells":float(np.max(np.abs(inc))),
        "max_abs_expected_increment_on_supported_cells":float(np.max(np.abs(target[sp]))),
        "full_ell0to3_injected_increment_only_if_144_supported":projections,
        "unsupported_RR_cells_have_no_injected_xi_or_full_multipole":not full,
        "not_a_physical_galaxy_wake_injection":True,
    }

def compare_supported(left,right):
    a,sa=left["xi"],left["support"]
    b,sb=right["xi"],right["support"]
    common=sa&sb
    if not np.any(common):
        return {"common_supported_cells":0,"max_abs_xi_difference_on_common_cells":None,
                "median_abs_xi_difference_on_common_cells":None,
                "no_common_support_warning":True}
    x=np.abs(a[common]-b[common])
    return {"common_supported_cells":int(common.sum()),
            "max_abs_xi_difference_on_common_cells":float(x.max()),
            "median_abs_xi_difference_on_common_cells":float(np.median(x)),
            "no_common_support_warning":False}

def synthetic_self_test(p):
    for cap in A02.CAPS:
        for tracer in PER_TRACER:
            _,idx=subset_indices(1,cap,tracer)
            a,b,nine,full=(idx[k] for k in LEVELS)
            if (len(a)!=600 or len(b)!=600 or len(nine)!=900 or len(full)!=1200
                or np.intersect1d(a,b).size or not np.all(np.isin(a,nine))
                or not np.array_equal(np.union1d(a,b),full)):
                raise AssertionError("Fixed nested and disjoint mock random split failed")
    sample=(np.arange(1200,dtype="f8"),)*4
    _,idx=subset_indices(1,"NGC","eBOSS_ELG")
    try:random_subcatalogue(sample,np.asarray([0,0]))
    except ValueError:pass
    else:raise AssertionError("Duplicated random row accepted")
    cat=(np.arange(600,dtype="f8"),np.zeros(600),
         np.linspace(.9,.999,600),np.ones(600))
    plus,diag=elg_random_stress(cat,.05)
    if (not np.array_equal(cat[3],np.ones(600))
        or not np.array_equal(plus[2],cat[2])
        or diag["factor_min"]<.95-1e-12):
        raise AssertionError("ELG RANDOM weight-only synthetic change failed")
    for wrong in (.06,-.06):
        try:elg_random_stress(cat,wrong)
        except ValueError:pass
        else:raise AssertionError("Unregistered E2 random perturbation accepted")
    # Synthetic 6x24 cross-LS with one isolated missing RR cell.
    ones=np.ones((6,24),dtype="f8")*11
    h={term:ones.copy() for term in LABELS}
    h["D1D2"]*=1.04
    norm={term:1. for term in LABELS}
    x,support=cross_landy_szalay(h,norm)
    rev={A02.MAPPING[k]:v[:,::-1].copy() for k,v in h.items()}
    revnorm={A02.MAPPING[k]:v for k,v in norm.items()}
    rx,rs=cross_landy_szalay(rev,revnorm)
    internal={"hist":h,"norm":norm,"rev_hist":rev,"rev_norm":revnorm,
              "xi":x,"support":support,"reverse_xi":rx}
    result=synthetic_dd_odd_pair_injection(internal,.02)
    if (result["max_abs_analytic_xi_increment_residual_on_supported_cells"]>=1e-11
        or result["max_abs_reverse_analytic_xi_increment_residual_on_supported_cells"]>=1e-11
        or result["max_abs_forward_reverse_injected_increment_mirror_residual_on_supported_cells"]>=1e-8
        or result["full_ell0to3_injected_increment_only_if_144_supported"] is None):
        raise AssertionError("Fixed +/-0.02 synthetic DD odd pair injection failed")
    # WRONG reverse sign must break physical signed-tracer mirror.
    mu=(A02.MU_EDGES[1:]+A02.MU_EDGES[:-1])/2.
    wrong_reverse={**rev,"D1D2":rev["D1D2"]*(1.+.02*mu[None,:])}
    wx,_=cross_landy_szalay(wrong_reverse,revnorm)
    correct_forward={**h,"D1D2":h["D1D2"]*(1.+.02*mu[None,:])}
    fx,_=cross_landy_szalay(correct_forward,norm)
    if np.max(np.abs(fx-wx[:,::-1]))<1e-3:
        raise AssertionError("Wrong same-sign reverse odd injection incorrectly accepted")
    h_missing={k:v.copy() for k,v in h.items()}
    h_missing["R1R2"][0,0]=0.
    sx,ss=cross_landy_szalay(h_missing,norm)
    sparse={**internal,"hist":h_missing,"xi":sx,"support":ss}
    rev_missing={k:v.copy() for k,v in rev.items()}
    rev_missing["R1R2"][0,-1]=0.
    rxm,rsm=cross_landy_szalay(rev_missing,revnorm)
    sparse["rev_hist"]=rev_missing
    sparse["reverse_xi"]=rxm
    outcome=synthetic_dd_odd_pair_injection(sparse,.02)
    if (outcome["positive_RR_supported_cell_count"]!=143
        or outcome["full_ell0to3_injected_increment_only_if_144_supported"] is not None):
        raise AssertionError("One missing RR cell did not disable all-six-bin multipole projection")
    # Constant global rescaling of ELG R hist and independent norms cancels.
    scale={k:v.copy() for k,v in h.items()}
    scale["D1R2"]*=1.04
    scale["R1R2"]*=1.04
    scaled_norm={**norm,"D1R2":1.04,"R1R2":1.04}
    clean,_=cross_landy_szalay(scale,scaled_norm)
    wrong,_=cross_landy_szalay(scale,norm)
    if not np.allclose(clean,x,rtol=0,atol=1e-13) or np.allclose(wrong,x,rtol=0,atol=1e-4):
        raise AssertionError("Independently normalized four-term LS negative control failed")

def run_case(mid,cap,cats,metadata,parent_e2):
    orig=parent_e2["scenarios"]["baseline_original"]
    for key,info in metadata.items():
        if info["selected_array_SHA256"]!=parent_e2["original_selected_array_SHA256"][key]:
            raise ValueError("Original E2 immutable selected array SHA changed: "+key)
        if key.startswith("eBOSS_ELG_") and info["ELG_exact_chunk_diagnostic"]!=parent_e2["original_exact_ELG_chunk_counts"][key]:
            raise ValueError("Original E2 ELG chunk metadata changed: "+key)
    records={"id":mid,"cap":cap,"status":"complete",
             "original_e2_selected_array_sha_replayed":True,
             "full_original_e2_all_three_pair_and_xi_SHA_replayed":False,
             "random_permutation_and_levels":{},"levels":{},
             "supported_cell_comparisons":{},
             "no_observed_rows_read":True,
             "no_physical_wake_injection_or_covariance":True}
    bytracer={}
    for tracer in PER_TRACER:
        seed,idx=subset_indices(mid,cap,tracer)
        cat=cats[tracer,"ran"]
        per_level={}
        for lev in LEVELS:
            sub=random_subcatalogue(cat,idx[lev])
            per_level[lev]={"n":len(idx[lev]),
                "selected_original_1200_index_SHA256":arr_sha(idx[lev]),
                "subcatalogue_four_vector_SHA256":catalogue_sha(sub)}
        records["random_permutation_and_levels"][tracer]={
            "fixed_numpy_PCG64_seed":seed,
            "permutation_original_1200_index_SHA256":arr_sha(np.random.default_rng(seed).permutation(1200)),
            "half_A_and_B_disjoint_and_cover_original_1200":True,
            "half_A_nested_in_900":True,
            "per_level":per_level}
        bytracer[tracer]={lev:random_subcatalogue(cat,idx[lev]) for lev in LEVELS}
        if catalogue_sha(bytracer[tracer]["full_1200"])!=metadata[tracer+"_ran"]["selected_array_SHA256"]:
            raise ValueError("Original 1200 RANDOM full selected catalogue SHA changed")
    originals={}
    baseline_internal={}
    for lev in LEVELS:
        rl,re=(bytracer[tr][lev] for tr in PER_TRACER)
        baseline,bint=calculate(cats,rl,re)
        baseline_internal[lev]=bint
        if lev=="full_1200":
            compare_e2_level(parent_e2,"baseline_original",baseline,
                {"modified_selected_weight_SHA256":arr_sha(re[3])})
        injections={str(a):synthetic_dd_odd_pair_injection(bint,a) for a in INJECTIONS}
        level_record={"LRG_R_rows":len(rl[0]),"ELG_R_rows":len(re[0]),
            "scenarios":{"baseline_original":baseline},
            "DD_pair_level_synthetic_odd_injections":injections,
            "stress_deltas_vs_same_level_baseline":{},
            "not_a_physical_wake_or_random_infinite_density_inference":True}
        if lev!="three_quarter_A_900":
            for scenario,amp in (("plus_5pct_ELG_R_z_ramp",.05),
                                 ("minus_5pct_ELG_R_z_ramp",-.05)):
                weight,diag=elg_random_stress(re,amp)
                r,internal=calculate(cats,rl,weight)
                for side,unaffected in (("forward_pair_terms",("D1D2","R1D2")),
                                        ("reverse_pair_terms",("D1D2","D1R2"))):
                    for term in unaffected:
                        if r[side][term]!=baseline[side][term]:
                            raise ValueError("ELG RANDOM stress altered non-ELG pair term")
                for side,affected in (("forward_pair_terms",("D1R2","R1R2")),
                                      ("reverse_pair_terms",("R1D2","R1R2"))):
                    for term in affected:
                        if r[side][term]["accepted_pairs"]!=baseline[side][term]["accepted_pairs"]:
                            raise ValueError("ELG RANDOM weight-only stress changed accepted geometric pair count")
                if lev=="full_1200":
                    compare_e2_level(parent_e2,scenario,r,diag)
                common=compare_supported(internal,bint)
                if not np.array_equal(internal["support"],bint["support"]):
                    raise ValueError("ELG RANDOM weight-only perturbation altered RR support mask")
                delta_by_ell=None
                if r["RR_supported_cells"]==144:
                    delta_by_ell={}
                    for ell in ELLS:
                        label=str(ell)
                        bv=np.asarray(baseline["original_physical_domain_multipoles"][label]["values_by_fixed_s_bin"])
                        av=np.asarray(r["original_physical_domain_multipoles"][label]["values_by_fixed_s_bin"])
                        dd=av-bv
                        delta_by_ell[label]={"delta_by_fixed_s_bin":dd.tolist(),
                            "max_abs_delta_across_s_bins":float(np.max(np.abs(dd)))}
                else:
                    # Never produce six-bin multipoles when ANY physical mu cell lacks RR.
                    if r["original_physical_domain_multipoles"] is not None:
                        raise ValueError("Sparse RR multipoles were silently projected")
                if lev=="full_1200" and delta_by_ell is not None:
                    for ell in ELLS:
                        old=parent_e2["diagnostic_deltas_from_original"][scenario][str(ell)]
                        if not np.allclose(
                            delta_by_ell[str(ell)]["delta_by_fixed_s_bin"],
                            old["delta_by_fixed_s_bin"],rtol=0,atol=1e-12):
                            raise ValueError("Original E2 full-R delta multipoles changed")
                level_record["scenarios"][scenario]={**r,"ELG_R_weight_stress":diag}
                level_record["stress_deltas_vs_same_level_baseline"][scenario]={
                    **common,"full_multipole_delta_if_144_supported":delta_by_ell}
        records["levels"][lev]=level_record
        originals[lev]=baseline
    records["full_original_e2_all_three_pair_and_xi_SHA_replayed"]=True
    for lev in LEVELS:
        if lev=="full_1200":
            continue
        records["supported_cell_comparisons"][lev+"_vs_full_1200"]=compare_supported(
            baseline_internal[lev],baseline_internal["full_1200"])
    records["supported_cell_comparisons"]["half_A_vs_half_B_600"]=compare_supported(
        baseline_internal["half_A_600"],baseline_internal["half_B_600"])
    return records

def run(p,a02p,original_a02,e2):
    pre=A02.preflight(a02p,require_local=True)
    g,r,gman,rman,ref,original,gproto,rproto,rawproto,original_pilot=pre
    paths,total=A02.resolve_72_sources(a02p,g,r,gman,rman,ref,gproto,rproto,rawproto)
    dest=ROOT/p["local_output"]
    phash=file_sha(PROTOCOL)
    if dest.exists():
        out=A02.load(dest)
        if (out.get("protocol_sha256")!=phash
            or out.get("parent_E2_archive_sha256")!=p["parent_A03E2_report_sha256"]
            or out.get("status") not in (STOP,PASS)
            or out.get("observed_odd_data_vector_read") is not False
            or not isinstance(out.get("cases"),dict)
            or not set(out["cases"]).issubset(set(e2["cases"]))):
            raise ValueError("Existing E3 checkpoint does not match exact preregistered protocol")
        if out.get("errors"):
            raise ValueError("Prior E3 checkpoint contains failed case; preserve and report unchanged; no outcome-driven rerun")
        for key,item in out["cases"].items():
            if item.get("status")!="complete" or item.get("full_original_e2_all_three_pair_and_xi_SHA_replayed") is not True:
                raise ValueError("Prior E3 checkpoint contains unsuccessful case; stop without overwriting")
            raw=json.dumps(item,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
            if hashlib.sha256(raw).hexdigest()!=out.get("per_case_checkpoint_SHA256",{}).get(key):
                raise ValueError("Prior successful E3 checkpoint was tampered")
        if out["status"]==PASS:
            if len(out["cases"])!=18 or len(out["per_case_checkpoint_SHA256"])!=18:
                raise ValueError("Prior E3 PASS checkpoint omitted original fixed cases")
            print("A03E3_ALREADY_COMPLETE_AFTER_REHASH_ALL_72",flush=True)
            return out
    else:
        out={
            "status":STOP,"protocol_sha256":phash,
            "parent_E2_archive_sha256":p["parent_A03E2_report_sha256"],
            "parent_A02_archive_sha256":p["parent_A02_full_sha256"],
            "all_72_full_gzip_SHA_rehashed_before_any_E3_FITS_row":True,
            "total_original_72_compressed_bytes":total,
            "fixed_ids":list(A02.IDS),"caps":list(A02.CAPS),
            "random_levels":list(LEVELS),"weight_scenarios":list(SCENARIOS),
            "DD_pair_injection_amplitudes":list(INJECTIONS),
            "cases":{},"per_case_checkpoint_SHA256":{},"errors":[],
            "observed_galaxy_rows_read":False,"observed_random_rows_read":False,
            "observed_odd_data_vector_read":False,"new_science_selection_applied":False,
            "physical_empirical_pair_window_certified":False,
            "physical_wake_injection_recovery_certified":False,
            "inferential_18D_covariance_computed":False,
            "not_a_detection_or_exclusion":True,
        }
    A02.atomic(dest,out)
    for mid in A02.IDS:
        for cap in A02.CAPS:
            casekey=f"{mid:04d}/{cap}"
            if casekey in out["cases"]:
                print("A03E3_REUSED_IMMUTABLE_CHECKPOINT",casekey,flush=True)
                continue
            if mid!=1 and any(out["cases"].get(f"0001/{x}",{}).get("status")!="complete" for x in A02.CAPS):
                raise ValueError("Both original 0001 cases must close before new mock IDs")
            cats,metadata={},{}
            oldcase=original_a02["cases"][casekey]
            e2case=e2["cases"][casekey]
            try:
                prior=next((x for x in original["cases"] if x["cap"]==cap),None) if mid==1 else None
                for tracer in A02.TRACERS:
                    for role in A02.ROLES:
                        path=paths[A02.key(mid,cap,tracer)+"/"+role]
                        rows=A02.source_header_rows(path,prior,tracer,role)
                        cat,meta=sample_catalogue(path,expected_rows=rows,
                            cap=cap,tracer=tracer,role=role,expected_highz=None,p=original_pilot)
                        k=tracer+"_"+role
                        if meta!=oldcase["input_sample_diagnostics"][k]:
                            raise ValueError("A02 exact source/selected sample/ELG chunk digest changed: "+k)
                        cats[tracer,role]=cat
                        metadata[k]=meta
                current=run_case(mid,cap,cats,metadata,e2case)
            except Exception as exc:
                out["cases"][casekey]={
                    "status":"incomplete_stop","id":mid,"cap":cap,
                    "errors":[str(exc)],
                    "roles_verified_before_failure":list(metadata),
                    "no_post_result_reselection":True}
                out["errors"].append(casekey+": "+str(exc))
                A02.atomic(dest,out)
                raise ValueError("E3 fixed mock/cap failed without post-hoc retuning: "+casekey+" "+str(exc)) from exc
            out["cases"][casekey]=current
            raw=json.dumps(current,sort_keys=True,separators=(",",":"),allow_nan=False).encode()
            out["per_case_checkpoint_SHA256"][casekey]=hashlib.sha256(raw).hexdigest()
            A02.atomic(dest,out)
            support={lev:current["levels"][lev]["scenarios"]["baseline_original"]["RR_supported_cells"]
                     for lev in LEVELS}
            print("A03E3_FIXED_MOCK_CASE",casekey,"RR_SUPPORT",support,flush=True)
    out["completed_cases"]=sum(x["status"]=="complete" for x in out["cases"].values())
    out["failed_cases"]=sum(x["status"]!="complete" for x in out["cases"].values())
    out["status"]=PASS if out["completed_cases"]==18 and out["failed_cases"]==0 and not out["errors"] else STOP
    A02.atomic(dest,out)
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    p=A02.load(PROTOCOL)
    dest=ROOT/p["local_output"]
    try:
        a02p,a02,e2,e2p=exact_gate(p)
        if args.self_test:
            A02.preflight(a02p,require_local=False)
            synthetic_self_test(p)
            print("A03E3_SOURCE_ONLY_AND_SYNTHETIC_NESTED_SPLIT_DD_INJECTION_NEGATIVES_OK",flush=True)
            print("PARENT_E2_CASES",len(e2["cases"]),"OBSERVED_ODD_DATA_READ",False,flush=True)
            return 0
        out=run(p,a02p,a02,e2)
    except Exception as exc:
        if args.self_test:
            raise
        if dest.exists():
            try:
                out=A02.load(dest)
                if out.get("status")==STOP and out.get("errors"):
                    # A failed report may already have been written in run(), or
                    # belong to an earlier attempt. Do not mutate either one's bytes.
                    print("A03E3_PRESERVED_EXISTING_FAILED_REPORT",dest,flush=True)
                    print("A03E3_NESTED_MOCK_RANDOM_SPLIT_AND_DD_ODD_INJECTION",STOP,flush=True)
                    print("REPORT",dest,flush=True)
                    print("ERRORS",*out["errors"],sep="\\n",flush=True)
                    print("CURRENT_EXCEPTION",str(exc),flush=True)
                    return 2
                out["status"]=STOP
                out["errors"]=list(dict.fromkeys(out.get("errors",[])+[str(exc)]))
            except (OSError,ValueError,KeyError,TypeError):
                out={"status":STOP,"errors":[str(exc)]}
        else:
            out={"status":STOP,"errors":[str(exc)]}
        out["observed_galaxy_rows_read"]=False
        out["observed_random_rows_read"]=False
        out["observed_odd_data_vector_read"]=False
        out["physical_empirical_pair_window_certified"]=False
        out["physical_wake_injection_recovery_certified"]=False
        out["inferential_18D_covariance_computed"]=False
        out["not_a_detection_or_exclusion"]=True
        A02.atomic(dest,out)
    print("A03E3_NESTED_MOCK_RANDOM_SPLIT_AND_DD_ODD_INJECTION",out["status"],flush=True)
    print("REPORT",dest,flush=True)
    print("COMPLETED_CASES",out.get("completed_cases",len([
        v for v in out.get("cases",{}).values() if v.get("status")=="complete"])),flush=True)
    if out.get("errors"):
        print("ERRORS",*out["errors"],sep="\n",flush=True)
        return 2
    print("OBSERVED_ODD_DATA_READ",out["observed_odd_data_vector_read"],flush=True)
    return 0 if out["status"]==PASS else 2

if __name__=="__main__":
    raise SystemExit(main())
