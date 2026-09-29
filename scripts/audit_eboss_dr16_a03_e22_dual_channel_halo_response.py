#!/usr/bin/env python3
"""E22 restricted linear-drag vs quadratic-mixed-B physical-GATE source-only math.

Never a physical halo response, new transfer, eBOSS prediction, or detection.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PINS={
"protocol":("source_data/eboss_dr16_a03_e22_dual_channel_halo_response_comparison_protocol_2026-09-29.json","5832921806cdbf4dd5ab3fba06e3ce6903d4f025"),
"E17D1_doc":("docs/EBOSS_DR16_A03_E17D1_RESTRICTED_CAUSAL_SOURCE_HISTORY_NONIDENTIFIABILITY_2026-09-27.md","3a536aa384eec1a202e77c532beb3a55c7333cb2"),
"E17D1_result":("source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_original_joint_causal_history_nonidentifiability.json","9f7e96b034779d0e20f4f9db61d5a642fda03975"),
"E20":("source_data/eboss_dr16_a03_e20_archived_exact_moment_and_cross_power_structural_QA_2026-09-29.json","c736ea3bca7bd0560ab07225294e1dd0930e0561"),
"E21":("source_data/eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json","cc3312c5313cce60fad4a88c3a56ef3d4e675e38"),
"E21_doc":("docs/EBOSS_DR16_A03_E21_ORIGINAL_E13_LINEAR_LOS_CHANNEL_AND_NORMALIZATION_2026-09-29.md","f418e2bd1546f01b2b1f140fb5a7b9e1bf427c43")
}
K=("0.001","0.002","0.003","0.005")
def need(valid,why):
    if not valid:
        raise ValueError("E22_SOURCE_ONLY_PHYSICS_SCOPE_STOP: "+why)
def gitblob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def load():
    data={}
    for name,(rel,sha) in PINS.items():
        path=ROOT/rel
        need(path.is_file() and not path.is_symlink(),"missing pinned source "+name)
        raw=path.read_bytes()
        need(gitblob(raw)==sha,"original E22 pinned Git source drift "+name)
        data[name]=raw if name.endswith("_doc") else json.loads(raw)
    p=data["protocol"];e17=data["E17D1_result"];e20=data["E20"];e21=data["E21"]
    need(p["frozen_parent_head"]=="008bdb4ce5954fdcfbb8eef5c1ac2e8b40e7a224" and
         all(p["STOP"].values()) and
         p["registration_honesty"].startswith("Original E10/E13/E17D1"),
         "E22 explicit prospective-new-QA and physical STOP")
    need(e17["observed_odd_SEALED"] and
         e17["real_halo_potential_history_or_HOD_calibrated"] is False and
         e17["full_physical_finiteK_bispectrum"]=="BLOCKED" and
         e20["no_observed_odd_access"] and
         e20["no_new_real_data_or_source_download"] and
         e21["no_observed_eBOSS_odd_access"] and e21["no_user_WSL"] and
         e21["new_physical_halo_galaxy_response"] is False and
         list(e21["original_long_modes"])==list(K),
         "parent physically missing halo history / observed seal")
    need(p["pinned_sources"]["E21_source_only_result"]["blob"]==PINS["E21"][1] and
         p["pinned_sources"]["E20_source_only_result"]["blob"]==PINS["E20"][1],
         "parent pin map drift")
    return data
def close(a,b,what,rel=1e-13,abs_=1e-12):
    need(isinstance(a,(int,float)) and isinstance(b,(int,float)) and
         math.isfinite(a) and math.isfinite(b) and
         math.isclose(a,b,rel_tol=rel,abs_tol=abs_),what)

def gamma_toy(H,gamma,dt):
    """dimensionless constant-Markov drag response u/v; NOT actual E17 halo."""
    need(H>=0 and gamma>=0 and dt>=0 and
         all(math.isfinite(x) for x in (H,gamma,dt)),
         "unphysical synthetic gamma/H/elapsed")
    if gamma==0 or dt==0:return 0.0
    return gamma/(H+gamma)*(-math.expm1(-(H+gamma)*dt))
def contrast(A_L,A_E,c_L,c_E):
    return A_L*c_E-A_E*c_L
def odd_linear(A_L,A_E,c_L,c_E,mu,C):
    return mu*contrast(A_L,A_E,c_L,c_E)*C
def kaiser_density_cross(k,mu,C,H):
    """Conjugate of -i*k*mu*u/H gives +i*k*mu, P_delta,u=+i*mu*C."""
    need(k>=0 and H>0,"synthetic Kaiser k/H invalid")
    return (complex(0,k*mu/H)*complex(0,mu*C))
def odd_quadratic(lambda_L,lambda_E,bL,bE,Q):
    return (lambda_L*bE-bL*lambda_E)*Q.imag
def reject(name,fn):
    try:fn()
    except (ValueError,KeyError,TypeError):
        print("E22_NEGATIVE_REJECT",name,flush=True)
    else:raise AssertionError("accepted negative E22 control "+name)

def limits(d):
    p=d["protocol"];e=d["E21"]
    ratios={}
    for key in K:
        row=e["original_long_modes"][key]
        a=row["states"]["Fplus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        b=row["states"]["Fminus"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        fd=row["states"]["FD"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        need(min(a,b,fd)>0,"original E21 nonphysical C at "+key)
        close(a-b,row["signed_Fplus_minus_Fminus"]["unnormalized_km_s_Mpc3"],
              "original E21 signed v contrast changed "+key)
        avg=(a+b)/2.
        rel=(a-b)/avg
        ratios[key]={"Fplus_over_Fminus":a/b,
                      "signed_contrast_over_Fplus_Fminus_mean":rel,
                      "signed_contrast_over_FD":(a-b)/fd}
    need(len(ratios)==4 and all(ratios[k]["signed_contrast_over_Fplus_Fminus_mean"]>0 for k in K),
         "original E13 physical-unit four-K cohort")
    # Two physical limiting cases from a deliberately SYNTHETIC constant
    # conformal-drag equation u'+(H+gamma)u=gamma*v_rel, u0=0.
    u0=gamma_toy(2.,0.,.5)
    uplus=gamma_toy(2.,1.,.5)
    u_origin=gamma_toy(2.,1.,0.)
    need(u0==0 and u_origin==0 and 0<uplus<1,
         "retarded toy common-force/differential-force limit")
    # If both halo tracers SHARE drag then same u, but different biases
    # and identical nonderivative Doppler selection still allow odd.
    shared_differential_u=uplus-uplus
    common_c=uplus
    shared_wake_odd=odd_linear(2.,3.,common_c,common_c,.4,1.)
    need(shared_differential_u==0 and shared_wake_odd!=0,
         "mistaken universal shared-wake dipole zero")
    # A single L-E imaginary odd constrains ONLY A_L*c_E-A_E*c_L:
    # (c_L,c_E) proportional to (A_L,A_E) is its null direction.
    need(contrast(2.,3.,4.,6.)==0 and
         contrast(2.,3.,1.,1.)==-1,
         "two-tracer response rank/null direction")
    for key in K:
        C=e["original_long_modes"][key]["states"]["FD"]["C_linear_P_delta_vR16_over_i_mu_km_s_Mpc3"]
        lef=odd_linear(2.,3.,1.,1.,.4,C)
        elf=odd_linear(3.,2.,1.,1.,.4,C)
        murev=odd_linear(2.,3.,1.,1.,-.4,C)
        close(lef,-elf,"LE/EL sign "+key)
        close(lef,-murev,"LOS mu sign "+key)
    kin=kaiser_density_cross(.05,.4,1.,2.)
    kinrev=kaiser_density_cross(.05,-.4,1.,2.)
    need(kin.imag==0 and kin.real!=0 and kin==kinrev,
         "standard spatial derivative must give real even mu²")
    # Quadratic Q=integral Gamma*B: Gaussian B=0 => Q=0.
    need(odd_quadratic(1.,1.,2.,3.,0j)==0,
         "false Gaussian quadratic channel")
    # Real isotropic parity-even scalar cross Q cannot give Im, even with B!=0.
    Q_even=complex(.7,0.)
    need(odd_quadratic(1.,1.,2.,3.,Q_even)==0,
         "false nonzero mixed B sufficient for LOS odd")
    # LOS-selected complex Q is SYNTHETIC, NOT original physical B/kernel.
    Q_selected=complex(.7,.3)
    need(odd_quadratic(1.,1.,2.,3.,Q_selected)!=0 and
         odd_quadratic(1.,1.,2.,3.,Q_selected.conjugate())==
         -odd_quadratic(1.,1.,2.,3.,Q_selected),
         "LOS-selected complex quadratic parity/tracer response")
    return ratios,uplus,shared_wake_odd

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,default=None)
    args=ap.parse_args()
    d=load()
    ratios,uplus,common_odd=limits(d)
    bad=copy.deepcopy(d)
    bad["protocol"]["STOP"]["no_unsealed_observed_galaxy_or_odd"]=False
    reject("OBSERVED_SEAL_TAMPER",
           lambda:need(all(bad["protocol"]["STOP"].values()),"observed seal violated"))
    reject("ZERO_DRAG_NONZERO_HALO_RESPONSE",
           lambda:need(gamma_toy(2.,0.,.5)!=0,"no differential force and u0=0"))
    reject("SHARED_HALO_VELOCITY_IMPLIES_UNIVERSAL_ODD_NULL",
           lambda:need(common_odd==0,"unequal biases can retain Doppler odd"))
    reject("FALSE_TRACER_SWAP_SAME_SIGN",
           lambda:need(odd_linear(2.,3.,1.,1.,.4,1.)==
                       odd_linear(3.,2.,1.,1.,.4,1.),
                       "LE and EL must conjugate"))
    reject("NONZERO_B_MUST_HAVE_ODD",
           lambda:need(odd_quadratic(1.,1.,2.,3.,complex(.7,0.))!=0,
                       "isotropic Q real"))
    reject("GAUSSIAN_QUADRATIC_FIRST_ORDER_NONZERO",
           lambda:need(odd_quadratic(1.,1.,2.,3.,0j)!=0,
                       "Gaussian first-order mixed B zero"))
    bad=copy.deepcopy(d)
    bad["E21"]["original_long_modes"]["0.003"]["signed_Fplus_minus_Fminus"]["unnormalized_km_s_Mpc3"]*=-1
    reject("ORIGINAL_E21_SIGNED_STATE_CONTRAST_TAMPER",lambda:limits(bad))
    print("E22_SEVEN_NEGATIVE_CONTROLS_PASS",flush=True)
    for k,v in ratios.items():
        print("E22_FROZEN_E21_PHYSICAL_VELOCITY_RATIO",k,
              format(v["Fplus_over_Fminus"],".14g"),
              format(v["signed_contrast_over_Fplus_Fminus_mean"],".14g"),
              flush=True)
    print("E22_CAUSAL_DRAG_TOY_ZERO_GAMMA_AND_SHARED_FORCE_CASES_PASS",flush=True)
    print("E22_GAUSSIAN_QUADRATIC_NULL_AND_SEPARATE_LOS_PARITY_GATE_PASS",flush=True)
    print("E22_NO_PHYSICAL_HALO_COEFFICIENT_OR_BISPECTRUM_NO_OBSERVED_ODD",flush=True)
    out={
      "stage":d["protocol"]["stage"],
      "status":"ORIGINAL_E13_SOURCE_RATIO_AND_DUAL_CHANNEL_ANALYTIC_LIMITING_CASES_PASS_NO_PHYSICAL_COUPLING",
      "original_source_git_blobs":{k:v[1] for k,v in PINS.items()},
      "ratios_four_original_E13_K_only":ratios,
      "synthetic_dimensionless_constant_drag_H2_gamma1_DeltaTau0p5_u_over_v":uplus,
      "synthetic_shared_wake_can_have_nonzero_two_tracer_Doppler_odd":common_odd!=0,
      "two_tracer_only_one_linear_response_contrast_identifiable":True,
      "gaussian_quadratic_mixed_B_zero_at_first_order_only":True,
      "nonzero_B_alone_not_sufficient_without_LOS_projection":True,
      "Kaiser_spatial_derivative_even_in_mu":True,
      "negative_controls_passed":7,
      "real_E17D1_halo_force_history_known":False,
      "real_linear_halo_gamma_or_galaxy_c_calculated":False,
      "real_physical_mixed_B_or_Gamma_calculated":False,
      "actual_eBOSS_24D_odd_or_covariance_accessed":False,
      "no_new_CLASS_E8_FITS_mock_ASDF_WSL":True}
    if args.out is not None:
        need(args.out.suffix==".json" and not args.out.exists(),"no overwrite/non-JSON")
        with args.out.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
        print("E22_SMALL_SOURCE_ONLY_OUTPUT",args.out,flush=True)
if __name__=="__main__":
    main()
