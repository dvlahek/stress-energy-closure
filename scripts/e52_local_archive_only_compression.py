#!/usr/bin/env python3
import json, math
from pathlib import Path
import numpy as np

P=Path("source_data/e51_conditioned_marked_ls_ninemock_result.json")
x=json.loads(P.read_text())
assert x["status"]=="PASS_MOCK_ONLY_CONDITIONED_ESTIMATOR_BACKGROUND_AND_E19_INJECTION_BASIS"
assert x["observed_odd_used"] is False

def summ(a):
    a=np.asarray(a,float); av=np.sort(np.abs(a)); n=len(a)
    def th(t):
        k=math.ceil((2*t-1)*n-1e-12)
        return float(av[k-1])
    return float(np.std(a,ddof=1)),th(.80),th(.90)

print("E52_ARCHIVE_ONLY_PASS")
for cap in ("NGC","SGC"):
    angles=[]
    print("CAP",cap)
    for field in ("A","B","C","D"):
        md=[]; mr=[]
        for c in x["cases"].values():
            if c["cap"]!=cap: continue
            b=np.asarray(c["baseline_12d_by_tag_field"][field],float)
            qp=np.asarray(c["unit_lambda_E19_plus_12d"],float)
            qm=np.asarray(c["unit_lambda_E19_minus_12d"],float)
            d=qm-qp
            rho=float(qp@qm/(qp@qp))
            r=qm-rho*qp
            md.append(float(d@b/(d@d)))
            mr.append(float(r@b/(r@r)))
            if field=="A": angles.append(float(c["plus_minus_postwindow_angle_rad"]))
        sd,th83,th94=summ(md)
        sds,t83s,t94s=summ(mr)
        print("TAG",field,
              "CAL_SD",sd,"CAL_A83",th83,"CAL_A94",th94,
              "SHAPE_SD",sds,"SHAPE_A83",t83s,"SHAPE_A94",t94s)
    print("ANGLE_MIN_MED_MAX",min(angles),float(np.median(angles)),max(angles))
print("NO_NEW_FITS_ASDF",True)
print("OBSERVED_ODD_USED",False)
