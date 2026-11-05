#!/usr/bin/env python3

from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
M=ROOT/"source_data/e63b_r6_freeze/e63b_r6_freeze_manifest.json"

def need(x,msg):
    if not x:
        raise RuntimeError(msg)

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

need(M.is_file(),"MISSING_FREEZE_MANIFEST")

m=json.loads(M.read_text())

need(
    m["stage"]=="E63B_R6_TO_E63D_FROZEN_LIGHTCONE_TRANSFER_CHAIN",
    "WRONG_FREEZE_STAGE"
)

need(m["decision"]["E63B_outerrim_lightcone_truth"]=="FAIL",
     "E63B_DECISION_DRIFT")

need(m["decision"]["E64_authorized"] is False,
     "E64_AUTHORIZATION_DRIFT")

need(m["decision"]["observed_eBOSS_authorized"] is False,
     "OBSERVED_EBOSS_GUARD_DRIFT")

need(m["decision"]["observed_odd_authorized"] is False,
     "OBSERVED_ODD_GUARD_DRIFT")

for name,x in m["provenance_files"].items():
    p=ROOT/x["path"]
    need(p.is_file(),f"MISSING_PROVENANCE {name}")
    need(sha(p)==x["sha256"],f"PROVENANCE_SHA_DRIFT {name}")

for stage in ("E63B","E63C","E63D"):
    p=ROOT/m[stage]["result_path"]
    need(p.is_file(),f"MISSING_RESULT {stage}")
    need(sha(p)==m[stage]["result_sha256"],
         f"RESULT_SHA_DRIFT {stage}")

cache=ROOT/m["retrieval"]["large_cache_local_path"]

if cache.is_file():
    need(
        sha(cache)==m["retrieval"]["large_cache_sha256"],
        "LOCAL_R6_CACHE_SHA_DRIFT"
    )
    print("LOCAL_R6_CACHE_SHA_PASS")
else:
    print("LOCAL_R6_CACHE_NOT_PRESENT_SKIP_BYTE_VERIFY")

print("E62_PERIODIC_TRUTH",m["decision"]["E62_periodic_nbody_truth"])
print("E63B_LIGHTCONE_TRUTH",m["decision"]["E63B_outerrim_lightcone_truth"])
print("E64_AUTHORIZED",m["decision"]["E64_authorized"])
print("OBSERVED_EBOSS_AUTHORIZED",m["decision"]["observed_eBOSS_authorized"])
print("OBSERVED_ODD_AUTHORIZED",m["decision"]["observed_odd_authorized"])
print("E63B_R6_TO_E63D_FROZEN_CHAIN_VERIFY_PASS")
