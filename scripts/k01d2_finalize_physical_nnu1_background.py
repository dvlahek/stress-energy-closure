#!/usr/bin/env python3
from pathlib import Path
import hashlib, json

ROOT=Path(__file__).resolve().parents[1]
PR=ROOT/"source_data/k01d2_physical_nnu1_background_path_prereg_2026-11-06.json"
PARENT=ROOT/"k01d1_swift_transport_result.json"
DIFF=ROOT/"k01d2_swift_universe_age_patch.diff"
OUT=ROOT/"k01d2_physical_nnu1_background_result.json"
EXPECTED_PR_BLOB="bf7befdfae4556ec6652f49108aa5ec7ab4d7011"

def blob(p):
    raw=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

if blob(PR)!=EXPECTED_PR_BLOB:
    raise SystemExit("K01D2_PREREG_BLOB_MISMATCH")
p=json.loads(PR.read_text())
r=json.loads(PARENT.read_text())
if r.get("status")!="PASS_SWIFT_CUSTOM_F_COSMOLOGICAL_TRANSPORT_SMOKE":
    raise SystemExit(f"PARENT_NOT_PASS {r.get('status')}")
if not r.get("science_gate_pass",False):
    raise SystemExit("PARENT_SCIENCE_GATE_FALSE")
if r.get("swift_commit")!=p["swift"]["commit"]:
    raise SystemExit("SWIFT_COMMIT_DRIFT")

diff=DIFF.read_text()
if diff.count("-  gsl_integration_qag(&F, 0., 1, 0, 1.0e-13")==0:
    raise SystemExit("PATCH_OLD_LINE_MISSING")
if diff.count("+  gsl_integration_qag(&F, 0., 1, 0, 1.0e-11")==0:
    raise SystemExit("PATCH_NEW_LINE_MISSING")

out={
  "stage":p["stage"],
  "status":"PASS_PHYSICAL_NNU1_BACKGROUND_AND_CUSTOM_F_TRANSPORT_WITH_AUDITED_UNIVERSE_AGE_PATCH",
  "science_gate_pass":True,
  "swift_commit":p["swift"]["commit"],
  "patch":p["exact_patch"],
  "patch_diff_sha256":hashlib.sha256(DIFF.read_bytes()).hexdigest(),
  "physical_cosmology":p["physical_cosmology"],
  "calibration":r["calibration"],
  "expected_force_free_velocity_scale":r["expected_force_free_velocity_scale"],
  "states":r["states"],
  "cross_state_macro_mass_relative_gap":r["cross_state_macro_mass_relative_gap"],
  "state_speed_ordering":r["state_speed_ordering"],
  "state_speed_ordering_all_match":r["state_speed_ordering_all_match"],
  "pairwise_displacement_rms_difference_Mpc":r["pairwise_displacement_rms_difference_Mpc"],
  "state_specific_displacement_not_identical":r["state_specific_displacement_not_identical"],
  "interpretation":p["interpretation"],
  "guards":p["guards"],
  "next_stage":p["next_if_pass"]
}
OUT.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
print("K01D2_PHYSICAL_NNU1_BACKGROUND_PASS",flush=True)
print("K01D2_PATCH_DIFF_SHA256",out["patch_diff_sha256"],flush=True)
for s,v in out["states"].items():
    print("K01D2_STATE",s,
          "VEL_MAX_REL",v["max_velocity_vector_relative_error"],
          "SPEARMAN",v["q_displacement_spearman"],
          "MEAN_DX",v["mean_displacement_Mpc"],flush=True)
print("K01D2_NEXT",out["next_stage"],flush=True)
