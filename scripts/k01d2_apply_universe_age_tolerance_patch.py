#!/usr/bin/env python3
from pathlib import Path
import hashlib
import difflib

p=Path("external/SWIFT/src/cosmology.c")
raw=p.read_bytes()
sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
expected="e28492319a4689c903444bd65df4116cc6460445"
if sha!=expected:
    raise SystemExit(f"PINNED_COSMOLOGY_BLOB_MISMATCH got={sha} expected={expected}")

old=raw.decode()
target='''  /* Integrate the time \\int_{0}^{1} dt */\n  gsl_integration_qag(&F, 0., 1, 0, 1.0e-13, GSL_workspace_size,\n                      GSL_INTEG_GAUSS61, space, &result, &abserr);\n  c->universe_age_at_present_day = result;\n'''
replacement='''  /* Integrate the time \\int_{0}^{1} dt */\n  gsl_integration_qag(&F, 0., 1, 0, 1.0e-11, GSL_workspace_size,\n                      GSL_INTEG_GAUSS61, space, &result, &abserr);\n  c->universe_age_at_present_day = result;\n'''
if old.count(target)!=1:
    raise SystemExit(f"TARGET_COUNT_FAIL {old.count(target)}")
new=old.replace(target,replacement,1)

# Ensure exactly one textual line changed.
diff=list(difflib.unified_diff(old.splitlines(),new.splitlines(),fromfile="a/src/cosmology.c",tofile="b/src/cosmology.c",lineterm=""))
changed=[x for x in diff if (x.startswith("+") or x.startswith("-")) and not x.startswith("+++") and not x.startswith("---")]
if len(changed)!=2:
    raise SystemExit(f"UNEXPECTED_DIFF_LINES {changed}")
if "1.0e-13" not in changed[0] or "1.0e-11" not in changed[1]:
    raise SystemExit(f"UNEXPECTED_DIFF_CONTENT {changed}")

p.write_text(new)
Path("k01d2_swift_universe_age_patch.diff").write_text("\n".join(diff)+"\n")

print("K01D2_PINNED_COSMOLOGY_BLOB_PASS",sha,flush=True)
print("K01D2_EXACT_ONE_LINE_PATCH_PASS",flush=True)
print("\n".join(diff),flush=True)
