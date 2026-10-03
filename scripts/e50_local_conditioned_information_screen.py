#!/usr/bin/env python3
import json, math
from pathlib import Path

P=Path("source_data/e49_direct_vTk_conditioned_velocity_tag_feasibility.json")
x=json.loads(P.read_text())
assert x["status"]=="PASS_SOURCE_ONLY_ALIGNMENT_SCREEN"
assert x["guardrails"]["observed_odd_used"] is False

tp=(0.00022278876746842203,0.00014150296987877026)
tm=(0.00031476164003239810,0.00021061456076441418)

def norm(v): return math.sqrt(sum(z*z for z in v))
def dot(a,b): return sum(x*y for x,y in zip(a,b))

np_=norm(tp); nm=norm(tm)
diff=norm((tm[0]-tp[0],tm[1]-tp[1]))/np_
theta=math.acos(dot(tp,tm)/(np_*nm))
s=math.sin(theta)

print("E50_CONDITIONED_INFORMATION_SCREEN_PASS")
print("E49_MIN_ABS_R",x["summary"]["min_abs_r_rel_baryon_across_F_states"])
print("FPM_DIFF_OVER_PLUS_NORM",diff)
print("FPM_ANGLE_RAD",theta)
print("FPM_ANGLE_DEG",theta*180/math.pi)
print("TAGGED_SNR_3SIG_AMPLITUDE_CALIBRATED",3/diff)
print("TAGGED_SNR_3SIG_SHAPE_ONLY",3/s)

for rr in ("1.00","0.90","0.80","0.70","0.50"):
    qs=[]
    for state in ("FD","plus","minus"):
        if rr=="1.00":
            q=abs(x["states"][state]["ideal_baryon_sign_factor"])
        else:
            q=abs(x["states"][state]["reconstruction_degradation"][rr]["gaussian_sign_correlation_factor"])
        qs.append(q)
    q=min(qs)
    print("RREC",rr,"QMIN",q,
          "UNDERLYING_SNR_3SIG_CAL",3/(diff*q),
          "UNDERLYING_SNR_3SIG_SHAPE",3/(s*q))

print("OBSERVED_ODD_USED",False)
