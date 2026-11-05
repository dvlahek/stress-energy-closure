#!/usr/bin/env python3
import json, pathlib

A_FD = -0.00201850894575251
A_PLUS = -0.00201220917438763
A_MINUS = -0.00202480872683428
HIGH_FRAC = 0.9422589
E21_RESPONSE_ALL4 = 0.000912183978546639
E19_RATIO = 0.02605283750649608

delta = A_MINUS - A_PLUS
rel = delta / A_FD
high = A_FD * HIGH_FRAC
highk_req = abs(delta)/abs(high)
q_born_only = abs(rel)/(1+abs(rel))

def product_bound(S_hat, B_S, chi_hat, B_chi):
    return abs(chi_hat)*B_S + abs(S_hat)*B_chi + B_S*B_chi

checks = 0
for S_hat,B_S,chi_hat,B_chi in [
    (1.2,.03,.8,.02),(-.4,.07,1.1,.01),(2.,0.,-.5,.2),(0.,.1,1.,0.)
]:
    b=product_bound(S_hat,B_S,chi_hat,B_chi)
    for ds in (-B_S,B_S):
        for dc in (-B_chi,B_chi):
            err=abs((chi_hat+dc)*(S_hat+ds)-chi_hat*S_hat)
            assert err <= b + 1e-15
            checks += 1

assert abs(rel-0.00624200971374) < 1e-13; checks += 1
assert 0 < highk_req < .01; checks += 1
assert 0 < q_born_only < .01; checks += 1

out={
  "stage":"E42_LOCAL_WSL_OPERATOR_LEVEL_PHYSICAL_CLOSURE_BUDGET",
  "status":"PURE_STDLIB_OFFLINE_QA_PASS",
  "frozen_E28":{
    "a_FD":A_FD,"a_Fplus":A_PLUS,"a_Fminus":A_MINUS,
    "delta_Fminus_minus_Fplus":delta,
    "abs_contrast_fraction_of_FD":abs(rel),
    "fraction_FD_force_k_gt_1":HIGH_FRAC,
    "highk_state_differential_fraction_sign_certificate":highk_req
  },
  "born_only_illustrative_q_certificate":q_born_only,
  "E41":{
    "E21_all4_response_fraction":E21_RESPONSE_ALL4,
    "E19_ratio_fraction":E19_RATIO
  },
  "operator_certificate":"For each F: B_Y <= |chi_hat| B_S + |S_hat| B_chi + B_S B_chi; for D=Yminus-Yplus, B_D<=B_Yminus+B_Yplus; B_D<|Dhat| preserves sign.",
  "qa_checks":checks,
  "qa_pass":True,
  "no_ASDF_FITS_CLASS_mock_or_observed_odd":True
}

p=pathlib.Path("source_data/e42_local_offline_operator_closure_budget.json")
p.parent.mkdir(exist_ok=True)
p.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

print("E42_WSL_OFFLINE_PASS")
print("QA_CHECKS", checks)
print("DELTA_A_FMINUS_MINUS_FPLUS", delta)
print("RELATIVE_TO_FD", abs(rel))
print("HIGHK_DIFFERENTIAL_CONTROL_PERCENT", 100*highk_req)
print("ILLUSTRATIVE_BORN_Q_IF_SOLE_ERROR", q_born_only)
print("E21_RESPONSE_ALL4_PERCENT", 100*E21_RESPONSE_ALL4)
print("REPORT", p)
