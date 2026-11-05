#!/usr/bin/env python3
import json, math, statistics
from pathlib import Path

P=Path("source_data/e45_archive_only_synthetic_odd_amplitude_recovery.json")
if not P.exists():
    raise SystemExit("MISSING "+str(P))
x=json.loads(P.read_text(encoding="utf-8"))
assert x["status"]=="PASS_DESCRIPTIVE_NOT_PHYSICAL_WAKE_NOT_INFERENCE"
assert x["interpretation_guardrails"]["observed_odd_used"] is False

def stat(cap):
    vals=[float(r["baseline_apparent_amplitude"]) for r in x["rows"] if r["cap"]==cap]
    av=sorted(abs(v) for v in vals)
    n=len(av)
    def th(target):
        need=math.ceil((2*target-1)*n-1e-12)
        return av[need-1], 0.5+0.5*need/n
    return {
        "n":n,
        "sample_sd":statistics.stdev(vals),
        "threshold_83pct":th(.80)[0],
        "actual_discrete_recovery_83pct":th(.80)[1],
        "threshold_94pct":th(.90)[0],
        "actual_discrete_recovery_94pct":th(.90)[1],
        "threshold_9of9":av[-1]
    }

out={
  "stage":"E46_LOCAL_ARCHIVE_ONLY_SENSITIVITY_BRIDGE_GATE",
  "status":"PASS_PHYSICAL_TEMPLATE_STILL_MISSING",
  "NGC":stat("NGC"),
  "SGC":stat("SGC"),
  "required_physical_bridge":
      "Provide y_phys[cap,12] = post-window dimensionless xi odd vector; then A_phys=(q dot y_phys)/(q dot q).",
  "direct_E28_acceleration_to_A_mapping_allowed":False,
  "observed_odd_used":False,
  "FITS_ASDF_opened":False
}
Q=Path("source_data/e46_local_archive_only_sensitivity_bridge_gate.json")
Q.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

print("E46_WSL_ARCHIVE_ONLY_PASS")
for cap in ("NGC","SGC"):
    s=out[cap]
    print(cap,"SD",s["sample_sd"],
          "A83",s["threshold_83pct"],
          "A94",s["threshold_94pct"],
          "A9OF9",s["threshold_9of9"])
print("PHYSICAL_TEMPLATE_AVAILABLE",False)
print("DIRECT_E28_TO_A_ALLOWED",False)
print("NO_FITS_ASDF_OPENED",True)
print("REPORT",Q)
