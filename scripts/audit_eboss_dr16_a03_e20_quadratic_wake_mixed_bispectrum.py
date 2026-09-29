#!/usr/bin/env python3
"""E20 quadratic wake structural TEST ONLY, not a physical EV or eBOSS result."""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PINS={
"protocol":("source_data/eboss_dr16_a03_e20_quadratic_wake_mixed_bispectrum_structural_protocol_2026-09-29.json","300f73be14306ccf45342592e842e30bda63b8c6"),
"E10":("source_data/eboss_dr16_a03_e10_frozen_conditional_angular_wind_parity_summary_2026-09-27.json","a9e0f49a58b9fd5e0ad975d4d14ba6521a648c28"),
"E19":("source_data/eboss_dr16_a03_e19_posthoc_conditional_Fourier_shape_and_marginal_RR_nonidentifiability_2026-09-29.json","fc5044684b84ff3de86c236988617ead79d8e7e6"),
"E6":("source_data/eboss_dr16_a03_e6_physical_template_to_empirical_window_bridge_protocol_2026-09-27.json","4df20872cff5b3d8aa09831ab50b22fa5f8b6581"),
"E17D1":("docs/EBOSS_DR16_A03_E17D1_RESTRICTED_CAUSAL_SOURCE_HISTORY_NONIDENTIFIABILITY_2026-09-27.md","3a536aa384eec1a202e77c532beb3a55c7333cb2"),
}
def need(ok,msg):
    if not ok:raise ValueError("E20_SOURCE_ONLY_STOP: "+msg)
def load():
    raw={}
    for name,(path,sha) in PINS.items():
        file=ROOT/path
        need(file.is_file() and not file.is_symlink(),"missing original "+name)
        r=file.read_bytes()
        need(hashlib.sha1(b"blob "+str(len(r)).encode()+b"\0"+r).hexdigest()==sha,
             "frozen original Git blob changed "+name)
        raw[name]=r
    p,e10,e19,e6=(json.loads(raw[name]) for name in ("protocol","E10","E19","E6"))
    need(p["registration_status"].startswith("PROSPECTIVE_ONLY_FOR_NEW_E20") and
         all(p["STOP"].values()) and
         p["base_audit_head"]=="8adc0881a296a645e3c6190e1641e143113d9072",
         "E20 prospective protocol or STOP changed")
    need(e10["observed_odd_sealed"] and
         e10["scientific_interpretation"]["actual_velocity_distribution_or_velocity_density_selection_correlation_not_computed"] and
         e19["observed_odd_SEALED"] and
         e19["conditional_shape_algebra"]["not_eBOSS_measurable_shape_or_sigma"] and
         e6["math_contract"]["pilot_joint_dimension"]==24 and
         e6["physics_prediction_status"].startswith("NO_EBOSS_ABSOLUTE_PHYSICAL_WAKE") and
         not e6["observed_odd_data_vector_read"],
         "original E6/E10/E19 physical nonclaim drift")
    need(p["source_pins"]["original_E17D1_document"]["git_blob"]==PINS["E17D1"][1],
         "E17D1 missing physical halo-history STOP drift")
    return p
def avg(rows):return sum(rows,F(0))/len(rows)
def moments(rows):
    need(len(rows)==4,"original fixed four-state witness changed")
    mu=[avg([r[j] for r in rows]) for j in range(3)]
    pair=[avg([r[j]*r[k] for r in rows]) for j,k in ((0,1),(0,2),(1,2))]
    third=avg([r[0]*r[1]*r[2] for r in rows])
    return mu,pair,third
def symmetry_case():
    # Jointly sign-symmetric but non-Gaussian discrete analogue:
    # V=X,D1=X,D2=X+Y all pair-correlated; odd third moment = 0.
    rows=[(F(x),F(x),F(x+y)) for x in (-1,1) for y in (-1,1)]
    mu,pair,third=moments(rows)
    need(mu==[0]*3 and pair==[1]*3 and third==0,
         "two-point correlation did not preserve odd-third zero")
    return rows,mu,pair,third
def non_gauss_case():
    # Moment-only witness: NOT a physical 3D bispectrum or halo model.
    # V=X,D1=Y,D2=X*Y all centered and pairwise uncorrelated, third=1.
    rows=[(F(x),F(y),F(x*y)) for x in (-1,1) for y in (-1,1)]
    mu,pair,third=moments(rows)
    need(mu==[0]*3 and pair==[0]*3 and third==1,
         "connected mixed third moment missing")
    return rows,mu,pair,third
def imaginary(bL,bE,lL,lE,Q):
    # P_LE linear order = lL*bE*Q + bL*lE*conjugate(Q).
    actual=(lL*bE*Q+bL*lE*Q.conjugate()).imag
    expected=(lL*bE-bL*lE)*Q.imag
    need(math.isclose(actual,expected,rel_tol=0,abs_tol=1e-14),
         "first-order cross-power tracer contrast broken")
    return actual
def negative(label,fn):
    try:fn()
    except (ValueError,KeyError,TypeError):
        print("E20_SYNTHETIC_NEGATIVE_REJECT",label,flush=True)
    else:raise AssertionError("accepted negative "+label)
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path)
    a=ap.parse_args()
    p=load()
    gr,mu_g,pair_g,third_g=symmetry_case()
    ng,mu_ng,pair_ng,third_ng=non_gauss_case()
    Q=complex(.25,.125) # purely synthetic, not measured spectrum
    LE=imaginary(2,3,1,1,Q)
    EL=imaginary(3,2,1,1,Q)
    eq=imaginary(2,2,1,1,Q)
    diff=imaginary(2,2,1,2,Q)
    need(LE==-EL and eq==0 and diff!=0,
         "real-tracer reversal or equal response bias contrast failed")
    for k,q in ((F(1,20),F(1,50)),(F(1,10),F(-1,50))):
        gamma=1j*(k-q) # local v.grad(delta) derivative prototype only
        gamma_reversed=1j*(-k+q)
        need(gamma_reversed==gamma.conjugate(),"real-field kernel reality")
    bad=copy.deepcopy(p)
    bad["STOP"]["no_observed_eBOSS_galaxies_or_24D_odd"]=False
    negative("OBSERVED_ODD_STOP_TAMPER",
             lambda:need(all(bad["STOP"].values()),"observed input forbidden"))
    bad=gr[:];bad[-1]=bad[0]
    negative("SIGN_SYMMETRY_BROKEN",
             lambda:need(moments(bad)[2]==0,"false zero third moment"))
    bad=ng[:];bad[-1]=(bad[-1][0],bad[-1][1],-bad[-1][2])
    negative("TRIPLE_WITNESS_BROKEN",
             lambda:need(moments(bad)[2]==1,"false nonzero third moment"))
    negative("FALSE_UNIVERSAL_EQUAL_BIAS_NULL",
             lambda:need(diff==0,"different tracer response is allowed"))
    negative("WRONG_TRACER_REVERSAL",
             lambda:need(LE==EL,"swap conjugates cross-power"))
    print("E20_FIVE_SYNTHETIC_NEGATIVES_PASS",flush=True)
    report={
      "stage":p["stage"],
      "status":"EXACT_SOURCE_ONLY_CENTRAL_SYMMETRY_MIXED_THIRD_AND_COMPLEX_CROSS_POWER_PASS",
      "E20_protocol_git_blob":PINS["protocol"][1],
      "E10_original_git_blob":PINS["E10"][1],
      "E19_original_git_blob":PINS["E19"][1],
      "sign_symmetric_correlated_moments":{"means":list(map(str,mu_g)),"cross_pairs":list(map(str,pair_g)),"mixed_third":str(third_g)},
      "non_gaussian_moment_only_witness":{"means":list(map(str,mu_ng)),"cross_pairs":list(map(str,pair_ng)),"mixed_third":str(third_ng),"not_physical_spatial_bispectrum":True},
      "synthetic_cross_power":{"LE_Im":LE,"EL_Im":EL,"equal_bias_equal_response_Im":eq,"equal_bias_different_response_Im":diff},
      "negative_controls_passed":5,"actual_B_vdelta_delta_calculated":False,
      "actual_halo_galaxy_response_calibrated":False,
      "actual_eBOSS_predicted_xi_or_covariance":False,
      "no_real_FITS_ASDF_CLASS_WSL_or_mock":True,"observed_odd_SEALED":True
    }
    if a.out:
        need(a.out.suffix==".json" and not a.out.exists(),"refuse overwrite")
        with a.out.open("x",encoding="utf-8") as f:
            json.dump(report,f,indent=2,allow_nan=False)
            f.write("\n")
        print("E20_SOURCE_ONLY_OUTPUT",a.out,flush=True)
    print("E20_NONZERO_GAUSSIAN_TWO_POINT_DOES_NOT_FORCE_MIXED_THREE_POINT",flush=True)
    print("E20_ZERO_PAIRWISE_DOES_NOT_EXCLUDE_CONNECTED_THREE_POINT",flush=True)
    print("E20_NO_PHYSICAL_BISPECTRUM_NO_OBSERVED_ODD_NO_WSL",flush=True)
if __name__=="__main__":main()
