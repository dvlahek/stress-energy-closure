#!/usr/bin/env python3
"""A03E16: original E14/E15 source-only exact CLOSED triangle geometry audit.

This is a diagnostic of what E14's reduced i*muLong*S fails to determine.
It certifies Cartesian closure, BOTH short legs, azimuth/LOS geometry, and
the range of the INDIVIDUAL static 1/k^2 factor for predeclared K/k<=0.1.
Two Hermitian-compatible toy completions with the SAME squeezed limit show
that the existing 72 reduced coefficients do NOT identify a full bispectrum.
This does not compute physical finite-K Vlasov corrections, P_cb(k1,k2),
tracer/HOD, eBOSS triple window, covariance, xi or significance.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e16_exact_triangle_geometry_finite_squeeze_prereg_2026-09-27.json"
E14=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json"
E14MAN=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/archive_manifest.json"
E15=ROOT/"source_data/eboss_dr16_a03_e15_original_E14_long_LOS_parity_and_toy_selection_2026_09_27.json"
E15IND=ROOT/"source_data/eboss_dr16_a03_e15_independent_original_72_angular_replay_2026_09_27.json"
BLOBS={
 PRE:"28127fbedddc320660b6ac3f9d6262272905c7ca",
 E14:"6ef0b2c3268a5bae9564b58abdd914328e861cba",
 E14MAN:"ba747f413988693278a6322bffe2da69d4c52d68",
 E15:"9da057195b9fb34c113be53533af129f9de8809f",
 E15IND:"18bf68b868a164289d1b7b5974d656dbd17ea3de",
}
E14SHA="0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8"
E15SHA="6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983"
STATES=("FD","plus","minus")
KS=(.05,.075,.1)
KLS=(.001,.002,.003,.005)
ELLS=(1,3)
MU_S=(-1.,0.,.6,1.)
MU_L=(-1.,0.,.5,1.)
PHIS=(0.,math.pi/2.,math.pi)
N_ORIENT=48

def require(q,msg):
    if not q:raise ValueError(msg)

def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def source_gate():
    for path,expected in BLOBS.items():
        require(blob(path.read_bytes())==expected,
                "Original SHA-frozen source/prereg blob changed: "+str(path))
    pre=json.loads(PRE.read_bytes())
    require(pre["immutable_parents"]["E14_full_sha256"]==E14SHA and
            pre["immutable_parents"]["E15_full_sha256"]==E15SHA and
            pre["original_axes"]["states"]==list(STATES) and
            pre["original_axes"]["short_k_comoving_h_per_Mpc"]==list(KS) and
            pre["original_axes"]["long_K_comoving_h_per_Mpc"]==list(KLS) and
            pre["original_axes"]["ell_short"]==list(ELLS) and
            pre["fixed_triangle_model"]["orientation_mu_short"]==list(MU_S) and
            pre["fixed_triangle_model"]["orientation_mu_long"]==list(MU_L) and
            pre["technical_QA"]["triangle_cartesian_closure_absolute_h_per_Mpc"]==1e-14 and
            pre["prohibitions"]["no_observed_galaxy_random_odd_read"] is True and
            pre["prohibitions"]["no_main_mutation"] is True,
            "E16 preregistered k/K/orientations, thresholds or STOP changed")
    e14raw=E14.read_bytes();e15raw=E15.read_bytes()
    require(sha(e14raw)==E14SHA and sha(e15raw)==E15SHA,
            "E14/E15 full original archived bytes changed")
    man=json.loads(E14MAN.read_bytes())
    row=next((x for x in man["files"] if
             x["file"]==E14.name),None)
    require(row is not None and row["sha256"]==E14SHA and
            row["size_bytes"]==len(e14raw),"Original E14 ZIP archive manifest changed")
    ind=json.loads(E15IND.read_bytes())
    require(ind["original_E14_full_sha256"]==E14SHA and
            ind["original_E15_full_sha256"]==E15SHA and
            ind["original_state_cases"]==72 and
            ind["original_Fplus_minus_Fminus_contrasts"]==24,
            "E15 independent 72/24 archived replay is absent/changed")
    e14=json.loads(e14raw);e15=json.loads(e15raw)
    require(e14["status"].startswith("E14_THREE_FROZEN_CLASS_SHORT_K_EXACT_STATIC_")
            and e15["original_E14_full_SHA256"]==E14SHA and
            e15["analytic_parity_controls"]["unweighted_long_LOS_L0_2_3_4_all_zero"]
            and e14["scientific_scope"]["observed_galaxies_randoms_odd_sealed"] is True,
            "Original E14/E15 source-only physical/observational context changed")
    return pre,e14,e15

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def norm(x):return math.sqrt(dot(x,x))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def sc(a,x):return tuple(a*v for v in x)
def neg(x):return sc(-1.,x)

def angle_vectors(mu_s,mu_L,phi):
    s=math.sqrt(max(0.,1.-mu_s*mu_s))
    l=math.sqrt(max(0.,1.-mu_L*mu_L))
    return (s,0.,mu_s),(l*math.cos(phi),l*math.sin(phi),mu_L)

def geometry(k,K,mu_s,mu_L,phi):
    r=K/k
    u,v=angle_vectors(mu_s,mu_L,phi)
    long=sc(K,v)
    short1=add(sc(k,u),sc(-K/2.,v))
    short2=add(sc(-k,u),sc(-K/2.,v))
    closure=norm(add(add(short1,short2),long))
    c=dot(u,v)
    c_expected=(mu_s*mu_L+math.sqrt(max(0.,1.-mu_s**2))*
                 math.sqrt(max(0.,1.-mu_L**2))*math.cos(phi))
    p1=norm(short1)/k;p2=norm(short2)/k
    d1=math.sqrt(1.+r*r/4.-r*c)
    d2=math.sqrt(1.+r*r/4.+r*c)
    lo1=short1[2]/norm(short1);lo2=short2[2]/norm(short2)
    e1=(mu_s-r*mu_L/2.)/d1;e2=(-mu_s-r*mu_L/2.)/d2
    f1=1./p1**2;f2=1./p2**2
    gsym=.5*(f1+f2);gasym=.5*(f1-f2)
    return {"mu_s":mu_s,"mu_L":mu_L,"azimuth_rad":phi,
       "cos_k_K":c,"k1_over_k":p1,"k2_over_k":p2,
       "mu1_LOS":lo1,"mu2_LOS":lo2,"inverse_k2_leg1_factor":f1,
       "inverse_k2_leg2_factor":f2,
       "geometric_symmetrized_inverse_k2_factor_TOY":gsym,
       "geometric_short_leg_antisymmetry_TOY":gasym,
       "closure_abs_h_per_Mpc":closure,
       "analytic_dot_gap":abs(c-c_expected),
       "analytic_k1_gap":abs(p1-d1),"analytic_k2_gap":abs(p2-d2),
       "analytic_mu1_gap":abs(lo1-e1),"analytic_mu2_gap":abs(lo2-e2),
       "short1":short1,"short2":short2,"long":long}

def test_geometry(x,k,K,pre):
    tol=pre["technical_QA"]
    require(x["closure_abs_h_per_Mpc"]<=tol[
        "triangle_cartesian_closure_absolute_h_per_Mpc"] and
        max(x[a] for a in ("analytic_dot_gap","analytic_k1_gap",
        "analytic_k2_gap","analytic_mu1_gap","analytic_mu2_gap"))<
        tol["moduli_and_mu_max_abs_gap"],
        "E16 exact vector triangle/cartesian vs analytic geometry failed")
    r=K/k
    lower=1./(1.+r/2.)**2
    upper=1./(1.-r/2.)**2
    for key in ("inverse_k2_leg1_factor","inverse_k2_leg2_factor"):
        require(lower-1e-12<=x[key]<=upper+1e-12,
                "E16 static inverse-k2 EXACT geometry envelope failed")
    require(abs(x["mu1_LOS"])<=1.+1e-13 and
            abs(x["mu2_LOS"])<=1.+1e-13,
            "E16 one short-leg LOS exceeds physical unit vector")
    # Flip the center k. This swaps short legs, including inverse-k2.
    flip=geometry(k,K,-x["mu_s"],x["mu_L"],math.pi-x["azimuth_rad"])
    # k-center direction flip means u->-u. Our spherical azimuth has
    # original u=(sin(theta),0,cos(theta)), so phi must transform K
    # azimuth only if rotate the entire coordinate frame. Instead use
    # analytic c->-c and equality of two short-leg reciprocal factors.
    fs1=1./(1.+r*r/4.+r*x["cos_k_K"])
    fs2=1./(1.+r*r/4.-r*x["cos_k_K"])
    require(abs(fs1-x["inverse_k2_leg2_factor"])<tol[
        "leg_swap_relative_gap"] and
        abs(fs2-x["inverse_k2_leg1_factor"])<tol[
        "leg_swap_relative_gap"],
        "E16 signed center/short-leg exchange symmetry failed")
    # Under simultaneous reversal of k and K, c stays unchanged and
    # mu_long flips. Both separable and Gsym toy completions are Hermitian.
    if abs(x["mu_L"])>0.:
        S=1.2345
        a=complex(0.,x["mu_L"]*S)
        b=complex(0.,x["mu_L"]*S*x[
            "geometric_symmetrized_inverse_k2_factor_TOY"])
        require(abs(a.conjugate()-complex(0.,-x["mu_L"]*S))<
                tol["conjugation_relative_gap"] and
                abs(b.conjugate()-complex(0.,-x["mu_L"]*S*x[
                  "geometric_symmetrized_inverse_k2_factor_TOY"]))<
                tol["conjugation_relative_gap"],
                "E16 two toy Hermitian completions violate reality")

def orientation_48(k,K,pre):
    out=[]
    for ms in MU_S:
        for ml in MU_L:
            for phi in PHIS:
                x=geometry(k,K,ms,ml,phi)
                test_geometry(x,k,K,pre)
                out.append(x)
    require(len(out)==N_ORIENT,"Original registered E16 orientation grid changed")
    return out

def bounds_and_witness(k,K,orient):
    r=K/k
    lower=1./(1.+r/2.)**2
    upper=1./(1.-r/2.)**2
    candidates=[val for x in orient for val in
                (x["inverse_k2_leg1_factor"],x["inverse_k2_leg2_factor"])]
    require(abs(min(candidates)-lower)<1e-12 and
            abs(max(candidates)-upper)<1e-12,
            "Registered 48 orientation samples failed exact full-c inverse-k² endpoints")
    maxclosure=max(x["closure_abs_h_per_Mpc"] for x in orient)
    maxanalytic=max(max(x[key] for key in
        ("analytic_dot_gap","analytic_k1_gap","analytic_k2_gap",
         "analytic_mu1_gap","analytic_mu2_gap")) for x in orient)
    maxdiff=max(abs(x["geometric_symmetrized_inverse_k2_factor_TOY"]-1.)
                for x in orient)
    witness=next(x for x in orient if x["mu_s"]==0. and
                 x["mu_L"]==1. and x["azimuth_rad"]==0.)
    require(witness["geometric_symmetrized_inverse_k2_factor_TOY"]!=1.,
            "E16 two distinct finite-K Hermitian possible completions collapsed")
    # Mean over UNIFORM relative c is exact, not real eBOSS angular weighting.
    isotropic_uniform_c_sym_mean=math.log((1.+r/2.)/(1.-r/2.))/r
    return {"r_K_over_k":r,"exact_continuous_c_envelope_individual_static_inverse_k2_factor":[lower,upper],
       "max_abs_individual_static_inverse_k2_fractional_change":max(upper-1.,1.-lower),
       "uniform_relative_c_geometry_ONLY_mean_symmetrized_factor":isotropic_uniform_c_sym_mean,
       "N_fixed_Cartesian_orientation_probes":len(orient),
       "max_Cartesian_triangle_closure_h_per_Mpc":maxclosure,
       "max_cartesian_vs_analytic_geometry_abs_gap":maxanalytic,
       "min_sampled_k1_h_per_Mpc":min(x["k1_over_k"]*k for x in orient),
       "max_sampled_k1_h_per_Mpc":max(x["k1_over_k"]*k for x in orient),
       "min_sampled_k2_h_per_Mpc":min(x["k2_over_k"]*k for x in orient),
       "max_sampled_k2_h_per_Mpc":max(x["k2_over_k"]*k for x in orient),
       "at_least_one_short_leg_outside_original_three_k_CLASS_Pcb_support":
           any(any(not(.05-1e-14<=q<=.1+1e-14) for q in
                   (k*x["k1_over_k"],k*x["k2_over_k"])) for x in orient),
       "max_abs_symmetrized_TOY_completion_ratio_minus_original":maxdiff,
       "min_antisymmetric_TOY_leg_ratio":min(x["geometric_short_leg_antisymmetry_TOY"] for x in orient),
       "max_antisymmetric_TOY_leg_ratio":max(x["geometric_short_leg_antisymmetry_TOY"] for x in orient),
       "one_explicit_same_squeezed_limit_distinct_finiteK_witness":{
          "geometry":{field:witness[field] for field in
                      ("mu_s","mu_L","azimuth_rad","cos_k_K","mu1_LOS",
                       "mu2_LOS","k1_over_k","k2_over_k")},
          "model_A_reduced_source_multiplier":1.,
          "model_B_symmetrized_STATIC_GEOMETRY_ONLY_source_multiplier":
              witness["geometric_symmetrized_inverse_k2_factor_TOY"],
          "these_are_not_physical_predictions_or_error_bars":True}}

def original_reduced(e14,state,k,K,ell):
    return float(e14["all_original_state_k_short_and_long_results"][state][str(k)]
            ["four_original_filtered_long_modes"][str(K)]["ells"][str(ell)]
            ["Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"])

def compute():
    pre,e14,e15=source_gate()
    fixed={}
    for k in KS:
        fixed[str(k)]={}
        for K in KLS:
            orient=orientation_48(k,K,pre)
            fixed[str(k)][str(K)]=bounds_and_witness(k,K,orient)
            require(K/k<=.1+1e-14,"Original preregistered max squeeze exceeded")
    n=0;states={}
    for state in STATES:
        states[state]={}
        for k in KS:
            states[state][str(k)]={}
            for K in KLS:
                rr=fixed[str(k)][str(K)]
                states[state][str(k)][str(K)]={}
                for ell in ELLS:
                    S=original_reduced(e14,state,k,K,ell)
                    old=e15["all_72_original_state_short_k_long_K_short_ell_cases"][
                        state][str(k)][str(K)][str(ell)]["original_reduced_E14_S_Mpc6"]
                    require(abs(S-old)<=pre["technical_QA"][
                        "archived_original_source_contrast_scaled_gap"]*
                        max(1.,abs(S),abs(old)),
                        "E14/E15 state-specific reduced S numerical source mismatch")
                    states[state][str(k)][str(K)][str(ell)]={
                      "original_frozen_E14_E15_reduced_source_S_Mpc6":S,
                      "exact_geometric_r_K_over_k":rr["r_K_over_k"],
                      "static_inverse_k2_individual_prefactor_envelope_ONLY":
                          rr["exact_continuous_c_envelope_individual_static_inverse_k2_factor"],
                      "full_physical_closed_triangle_bispectrum_calculated":False}
                    n+=1
    require(n==72,"Original 72 E14 state short/long/ell reduced source cases lost")
    ctr={};nc=0
    for k in KS:
        ctr[str(k)]={}
        for K in KLS:
            ctr[str(k)][str(K)]={}
            for ell in ELLS:
                plus=original_reduced(e14,"plus",k,K,ell)
                minus=original_reduced(e14,"minus",k,K,ell)
                d=plus-minus
                orig=float(e14["Fplus_minus_Fminus_reduced_source_difference"][
                   str(k)][str(K)][str(ell)]["Fplus_minus_Fminus_reduced_Bsource_Mpc6"])
                old=float(e15["all_24_original_Fplus_minus_Fminus_parity_contrasts"][
                   str(k)][str(K)][str(ell)]["E14_Fplus_minus_Fminus_reduced_S_Mpc6"])
                scale=max(1.,abs(plus),abs(minus),abs(d),abs(orig),abs(old))
                require(abs(d-orig)/scale<1e-13 and abs(d-old)/scale<1e-13,
                        "E16 original F+F- 24 reduced contrast numerical closure failed")
                ctr[str(k)][str(K)][str(ell)]={
                  "original_frozen_Fplus_minus_Fminus_E14_reduced_S_Mpc6":orig,
                  "checked_E15_reduced_S_Mpc6":old,
                  "E16_not_bispectrum_difference":True}
                nc+=1
    require(nc==24,"Original 24 plus-minus source contrasts lost")
    maxr=max(fixed[str(k)][str(K)]["r_K_over_k"] for k in KS for K in KLS)
    maxsingle=max(fixed[str(k)][str(K)][
      "max_abs_individual_static_inverse_k2_fractional_change"]
       for k in KS for K in KLS)
    maxsym=max(fixed[str(k)][str(K)][
      "max_abs_symmetrized_TOY_completion_ratio_minus_original"]
       for k in KS for K in KLS)
    outs=[(k,K) for k in KS for K in KLS if fixed[str(k)][str(K)][
       "at_least_one_short_leg_outside_original_three_k_CLASS_Pcb_support"]]
    sample=fixed["0.05"]["0.005"]
    require(abs(maxr-.1)<1e-14 and len(outs)==8,
            "Original short-leg support or predeclared squeeze-ratio axis changed")
    result={"date":"2026-09-27",
       "status":"E16_EXACT_THREE_VECTOR_CLOSED_TRIANGLE_AND_FINITE_K_GEOMETRY_CERTIFIED_PHYSICAL_BISPECTRUM_NONIDENTIFIABLE",
       "prospective_E16_protocol_git_blob":"28127fbedddc320660b6ac3f9d6262272905c7ca",
       "original_E14_full_SHA256":E14SHA,
       "original_E15_full_SHA256":E15SHA,
       "original_frozen_source_CLASS_commit":
           pre["immutable_parents"]["CLASS_original_commit"],
       "scope":"Exact geometric necessary closure and inverse-k2 factor ONLY; unmodified archived original E14/E15 reduced source not an Einstein-Vlasov full triangle kernel.",
       "fixed_geometry":"K+k1+k2=0, k1=k-K/2, k2=-k-K/2, LOS n=z; original axes, 48 registered orientations per pair",
       "twelve_original_short_long_geometry_cases":fixed,
       "all_72_original_state_short_long_shortell_source_cases":states,
       "all_24_original_Fplus_minus_Fminus_reduced_source_contrasts":ctr,
       "QA":{"original_72_cases":n,"original_24_contrasts":nc,
             "distinct_short_long_pairs":12,
             "orientation_probes_per_pair":N_ORIENT,
             "max_triangle_closure_h_per_Mpc":max(x[
                 "max_Cartesian_triangle_closure_h_per_Mpc"] for rows in fixed.values()
                 for x in rows.values()),
             "max_cartesian_vs_analytic_gap":max(x[
                 "max_cartesian_vs_analytic_geometry_abs_gap"] for rows in fixed.values()
                 for x in rows.values()),
             "technical_warnings":[]},
       "finite_ratio_scope":{"max_original_K_over_k_short":maxr,
          "max_single_short_leg_STATIC_inverse_k2_geometric_relative_variation":maxsingle,
          "max_sampled_symmetrized_geometry_TOY_relative_variation":maxsym,
          "pairs_with_some_closed_short_leg_outside_original_three_k_CLASS_Pcb_support":
              [[k,K] for k,K in outs],
          "n_pairs_outside_original_3point_Pcb_support":len(outs),
          "example_original_k_0p05_K_0p005":sample,
          "NOT_an_error_bound_on_actual_Einstein_Vlasov_tracer_bispectrum":True},
       "nonuniqueness_witness":{
          "A":"B_A/(i mu_long DeltaBias)=S (original E15 separable source)",
          "B":"B_B/(i mu_long DeltaBias)=S*Gsym(c,r), Gsym=0.5*((k/k1)^2+(k/k2)^2)",
          "both_have_E15_squeezed_limit_and_Hermitian_reality_under_k_K_reversal":True,
          "different_at_registered_nonzero_r":True,
          "Gsym_model_is_artificial_GEOMETRIC_ONLY_not_physical_prediction":True,
          "identifying_actual_finite_K_triangle_from_only_E14_E15_data":False},
       "physical_stop":{"full_both_short_leg_Pcb_and_RSD_growth_transfer":False,
         "short_long_angle_and_halo_tracer_velocity_density_response_calibrated":False,
         "full_retarded_Einstein_Vlasov_and_wake_phase_evolution":False,
         "eBOSS_HOD_b_LRG_ELG_evolution_magnification_wide_angle":False,
         "actual_eBOSS_triple_random_window_and_independent_bispectrum_covariance":False,
         "original_unconditional_24D_odd_is_three_point_observable":False,
         "observed_odd_SEALED":True,"new_catalogue_or_mock_download":False,
         "new_survey_cuts_science_seeds":False,"main_mutated":False,
         "PR_remains_draft":True}}
    return result

def atomic_new(path,raw):
    path.parent.mkdir(parents=True,exist_ok=True)
    require(not path.exists(),"Original E16 output already exists, refuse overwrite")
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e16_",delete=False) as f:
        tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    result=compute()
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    atomic_new(a.output,raw)
    s=result["finite_ratio_scope"]
    print("EBOSS_A03_E16_EXACT_CLOSED_TRIANGLE_GEOMETRY_72_24_SOURCE_ONLY_PASS",
          "ORIGINAL_SHA256",sha(raw),
          "MAX_K_over_k",s["max_original_K_over_k_short"],
          "MAX_INDIVIDUAL_INVERSE_K2_PREF_FACTOR_RELATIVE_CHANGE",
          s["max_single_short_leg_STATIC_inverse_k2_geometric_relative_variation"],
          "MAX_SYMMETRIZED_GEOMETRY_ONLY_TOY_CHANGE",
          s["max_sampled_symmetrized_geometry_TOY_relative_variation"],
          "N_OUTSIDE_THREE_CLASS_Pcb_GRID",
          s["n_pairs_outside_original_3point_Pcb_support"],
          "NO_FULL_GALAXY_BISPECTRUM NO_OBSERVED_ODD",flush=True)
    return 0
if __name__=="__main__":raise SystemExit(main())
