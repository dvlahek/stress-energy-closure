#!/usr/bin/env python3
import json, math, argparse
from pathlib import Path
import numpy as np

C=299792.458
SIGMA=135.86206681100845

ap=argparse.ArgumentParser()
ap.add_argument("--input",default="source_data/e47_local_physical_response_basis_to_e45.json")
ap.add_argument("--b-lrg",type=float,default=2.0)
ap.add_argument("--b-elg",type=float,default=1.0)
ap.add_argument("--growth-f",type=float,default=1.0)
args=ap.parse_args()

x=json.loads(Path(args.input).read_text())
assert x["status"]=="PASS_PARAMETERIZED_PHYSICAL_RESPONSE_BASIS_NOT_CALIBRATED_GALAXY_MODEL"
rows=x["rows"]
u1=SIGMA/C

def minD(cap,m):
    rr=[r for r in rows if r["cap"]==cap]
    best=(float("inf"),None)
    for t in np.linspace(-1,1,10001):
        for uL,uE in ((1,t),(-1,t),(t,1),(t,-1)):
            req=[]
            for r in rr:
                a=float(r["alpha_A_per_unit_g0"])
                b=float(r["beta_A_per_unit_g2"])
                k=a*(args.b_lrg*uE-args.b_elg*uL)+b*args.growth_f*(uE-uL)
                req.append(float("inf") if abs(k)<1e-30 else abs(float(r["baseline_A0"]))/abs(k))
            req.sort()
            if req[m-1]<best[0]:
                best=(req[m-1],(uL,uE))
    return best

out={"stage":"E48_LOCAL_PERTURBATIVITY_TEST",
     "b_L":args.b_lrg,"b_E":args.b_elg,"growth_f":args.growth_f,
     "sigma_R16_LOS_km_s":SIGMA,"U_1sigma":u1,"caps":{},
     "observed_odd_used":False}

for cap in ("NGC","SGC"):
    out["caps"][cap]={}
    for m,label in ((6,"83pct"),(8,"94pct"),(9,"100pct")):
        D,direction=minD(cap,m)
        out["caps"][cap][label]={
            "D_required":D,
            "optimal_direction":direction,
            "D_times_U_1sigma":D*u1
        }

Q=Path("source_data/e48_local_perturbativity_test.json")
Q.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print("E48_WSL_PASS")
print("U_1SIGMA",u1)
print("D_ORDER1_LIMIT",1/u1)
print("D_10PCT_LIMIT",.1/u1)
for cap in ("NGC","SGC"):
    print("CAP",cap,out["caps"][cap])
print("OBSERVED_ODD_USED",False)
print("REPORT",Q)
