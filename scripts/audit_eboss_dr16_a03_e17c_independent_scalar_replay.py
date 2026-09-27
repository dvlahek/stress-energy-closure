#!/usr/bin/env python3
"""Independent scalar E17C frozen-q and E16-two-leg streaming replay.

Standard library only. Reconstructs every geometry and integral from original
E8 CSV; never imports the original E17C implementation, NumPy or CLASS.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17c_restricted_retarded_kinetic_propagator_prereg_2026-09-27.json"
P_BLOB="33defef2020f2931756e3ef7a3959d905d539165"
CSV=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17B=ROOT/"source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_original_joint_native_CLASS_transfer_resolution.json"
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
      "minus":"Fminus_CLASS_normalized"}
KS=(.05,.075,.1)
KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.)
ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi)
S0=(0.,.25,1.)
DEFAULT=ROOT/"eboss_workspace/a03_physics_source/e17c_restricted_retarded_kinetic"

def require(ok,why):
    if not ok:
        raise ValueError("E17C_INDEPENDENT_FAIL_CLOSED: "+why)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def scalar_input():
    require(blob(PRE.read_bytes())==P_BLOB,"prospective E17C protocol Git blob changed")
    p=json.loads(PRE.read_bytes())
    par=p["immutable_parents"]
    for path,key in ((CSV,"E8_original_4000q_CSV_sha256"),
                     (E16,"E16_original_full_sha256"),
                     (E17A,"E17A_original_joint_sha256"),
                     (E17B,"E17B_original_joint_sha256")):
        require(sha(path.read_bytes())==par[key],"immutable original "+key)
    require(p["frozen_geometry"]["states"]==list(STATES)
            and p["frozen_geometry"]["tracer_order"]==["LRG","ELG"]
            and p["derived_restricted_kinetic_problem"][
                "prospectively_fixed_dimensionless_center_phases"]==list(S0)
            and all(p["stop"].values()),"scope changed")
    with CSV.open(newline="",encoding="utf-8") as fh:
        source=list(csv.DictReader(fh))
    require(len(source)==4000,"frozen source sample count")
    q=[float(x["q_dimensionless"]) for x in source]
    f={state:[float(x[col]) for x in source] for state,col in COLS.items()}
    require(all(q[i+1]>q[i] for i in range(len(q)-1))
            and all(v>0 and math.isfinite(v) for seq in f.values() for v in seq),
            "frozen q or source positivity")
    return p,q,f

def scalar_geometry():
    """Independent analytic dot-product geometry, not original Cartesian worker."""
    for k in KS:
        for K in KLS:
            r=K/k
            for ms in MS:
                for ml in ML:
                    for j,phi in enumerate(PH):
                        c=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(
                            max(0.,1.-ml*ml))*math.cos(phi)
                        k1=k*math.sqrt(1.+r*r/4.-r*c)
                        k2=k*math.sqrt(1.+r*r/4.+r*c)
                        yield (k,K,ms,ml,j,phi,k1,k2)

def weights(q,f):
    """Distinct node-weight quadrature, not np.trapz/np.sinc."""
    v=[q[i]*q[i]*f[i] for i in range(len(q))]
    w=[0.]*len(q)
    w[0]=.5*(q[1]-q[0])*v[0]
    w[-1]=.5*(q[-1]-q[-2])*v[-1]
    for i in range(1,len(q)-1):
        w[i]=.5*(q[i+1]-q[i-1])*v[i]
    return w,math.fsum(w)

def integral(q,w,s):
    if s==0.:
        return math.fsum(w)
    return math.fsum(wi*(math.sin(s*qi)/(s*qi) if s*qi else 1.)
                     for qi,wi in zip(q,w))

def compare(original,joint,p,q,f):
    require(joint["original_E16_geometries"]==576
            and joint["three_state_leg_phase_samples"]==10368
            and joint["full_physical_bispectrum"]=="BLOCKED"
            and joint["observed_odd_SEALED"] is True,
            "original joint source-only STOP")
    geometry=list(scalar_geometry())
    require(len(geometry)==576,"analytic original E16 geometry count")
    maximum=0.
    n=0
    digest={}
    for state in STATES:
        src=original/state
        raw=src.read_bytes()
        require(sha(raw)==joint["state_reports"][state]["sha256"],
                "original E17C worker report SHA mismatch: "+state)
        obj=json.loads(raw)
        rows=obj["original_geometries_both_legs_and_phase_characteristic"]
        require(obj["state"]==state and len(rows)==576
                and obj["n_leg_phase_values"]==3456 and obj["observed_odd_read"] is False
                and obj["full_physical_finiteK_halo_tracer_B"]=="BLOCKED",
                "state source-only scope/geometry missing")
        w,den=weights(q,f[state])
        require(den>0. and math.isfinite(den),"independent normalization")
        maximum=max(maximum,abs(den-obj["truncated_q_number_weighted_integral"])/max(den,1.))
        for expected,row in zip(geometry,rows):
            k,K,ms,ml,j,phi,k1,k2=expected
            require((row["k"],row["K"],row["mu_s"],row["mu_L"],row["phi_index"])==
                    (k,K,ms,ml,j),"independent orientation or pair order")
            for leg,mod in enumerate((k1,k2)):
                stored=row["leg_moduli_h_Mpc"][leg]
                maximum=max(maximum,abs(stored-mod)/max(1.,abs(stored),abs(mod)))
                samples=row["per_leg_three_dimensionless_phases"][leg]
                require(len(samples)==3,"three frozen dimensionless phase nodes")
                for i,s0 in enumerate(S0):
                    got=samples[i]
                    s=s0*mod/.05
                    ii=integral(q,w,s)
                    c=ii/den
                    require(got["s_center"]==s0,"dimensionless phase choice changed")
                    maximum=max(maximum,abs(got["s_leg"]-s)/max(1.,abs(s)))
                    maximum=max(maximum,abs(got["I_unnormalized_truncated_q"]-ii)/den)
                    maximum=max(maximum,abs(got["C_normalized"]-c))
                    n+=1
        digest[state]={"original_sha256":sha(raw),"replayed_samples":len(rows)*6}
    require(n==10368,"missing fixed state/leg/phase samples")
    require(maximum < p["QA_before_numerics"][
        "direct_numpy_vs_independent_stdlib_trapezoid_scaled_gap"],
        "source-only independent full scalar mismatch")
    return maximum,n,digest

def new_file(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,"existing audit output changed")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17ci_",delete=False) as fh:
            temp=Path(fh.name)
            fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17C_INDEPENDENT_AUDIT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--output-dir",type=Path,default=DEFAULT)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    args=a.parse_args()
    p,q,f=scalar_input()
    require(len(list(scalar_geometry()))==576,"geometry algebra count")
    if args.self_test:
        w,d=weights(q,f["FD"])
        require(abs(integral(q,w,0.)/d-1.)<1e-14,"unit free flight")
        print("E17C_INDEPENDENT_PARENT_AND_SCALAR_PREFLIGHT_PASS",flush=True)
        return
    path=args.output_dir/"e17c_original_joint_restricted_retarded_streaming.json"
    raw=path.read_bytes()
    joint=json.loads(raw)
    if args.negative_control:
        copy=json.loads(raw)
        copy["state_reports"]["FD"]["sha256"]="0"*64
        try:
            compare(args.output_dir,copy,p,q,f)
        except ValueError as exc:
            require("SHA mismatch" in str(exc),"tamper did not fail at expected gate")
            print("E17C_INDEPENDENT_TAMPER_NEGATIVE_CONTROL_REJECTED",flush=True)
            return
        raise AssertionError("tampered original result was accepted")
    maximum,n,digests=compare(args.output_dir,joint,p,q,f)
    result={"date":"2026-09-27",
            "status":"E17C_INDEPENDENT_STDLIB_ALL_FROZEN_Q_BOTH_LEG_CAUSAL_PROPAGATOR_REPLAY_PASS_FULL_B_BLOCKED",
            "original_joint_SHA256":sha(raw),
            "original_E8_CSV_SHA256":p["immutable_parents"]["E8_original_4000q_CSV_sha256"],
            "original_E16_full_SHA256":p["immutable_parents"]["E16_original_full_sha256"],
            "state_original_sha256_and_counts":digests,
            "replayed_original_E16_geometries":576,
            "replayed_state_leg_phase_values":n,
            "max_scaled_independent_scalar_replay_gap":maximum,
            "independent_method":"Pure standard library CSV and analytic triangle k1,k2; own node-weight scalar trapezoid and sinc, no original runner import/NumPy/CLASS.",
            "free_streaming_only_not_forced_halo_or_retarded_tracer_response":True,
            "full_physical_finiteK_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,
            "new_catalogue_mock_download":False,
            "main_untouched_PR_draft":True}
    new_file(args.output_dir/"e17c_independent_scalar_source_only_replay.json",result)
    print("E17C_INDEPENDENT_10368_SOURCE_ONLY_VALUES_PASS","MAX",maximum,
          "PHYSICAL_HALO_TRACER_B_BLOCKED",flush=True)

if __name__=="__main__":
    main()
