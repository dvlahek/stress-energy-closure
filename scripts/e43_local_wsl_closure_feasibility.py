#!/usr/bin/env python3
import json
from pathlib import Path

P = Path("source_data/e42_local_offline_operator_closure_budget.json")
if not P.exists():
    raise SystemExit("MISSING " + str(P))

x = json.loads(P.read_text(encoding="utf-8"))
assert x["qa_pass"] is True
assert x["status"] == "PURE_STDLIB_OFFLINE_QA_PASS"

delta = abs(x["frozen_E28"]["delta_Fminus_minus_Fplus"])
fd = abs(x["frozen_E28"]["a_FD"])
rel = x["frozen_E28"]["abs_contrast_fraction_of_FD"]

assert abs(delta/fd-rel) < 1e-15

pair_budget = delta/4.0
state_term_budget = delta/8.0

out = {
    "stage": "E43_LOCAL_WSL_CLOSURE_FEASIBILITY",
    "status": "PASS_NO_PHYSICAL_SIGN_CERTIFICATE",
    "absolute_total_sign_budget": delta,
    "percent_of_abs_FD": 100*rel,
    "equal_four_category_pair_budget": pair_budget,
    "equal_eight_term_budget": state_term_budget,
    "current_bounded_terms": {
        "incoming": False,
        "history": False,
        "profile": False,
        "Born": False,
        "tracer_chi": False
    },
    "recommended_next_gate":
        "small same-object multi-epoch dynamical summary; no whole-column ASDF scan",
    "observed_odd_remains_SEALED": True
}

q = Path("source_data/e43_local_wsl_closure_feasibility.json")
q.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")

print("E43_WSL_PASS")
print("TOTAL_SIGN_BUDGET", delta)
print("PERCENT_OF_ABS_FD", 100*rel)
print("PAIRWISE_CATEGORY_BUDGET_EQUAL4", pair_budget)
print("PER_STATE_TERM_BUDGET_EQUAL8", state_term_budget)
print("PHYSICAL_SIGN_CERTIFIED", False)
print("REPORT", q)
