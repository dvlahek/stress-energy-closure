#!/usr/bin/env python3
"""E17D0 independent pure-stdlib piecewise-linear integration-by-parts audit.

Independent derivative-free formula including finite q boundary, original E16
analytic triangle reconstruction; NO original worker import, NumPy or CLASS.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d0_linearized_vlasov_external_potential_source_prereg_2026-09-27.json"
P_BLOB="66df071f0750201e792ed2b3ffc81e28b13dfce6"
CSV=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
ARCH=ROOT/"source_data/eboss_dr16_a03_e17c_archived_CI_2026_09_27"
PARENTS=[
 (ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json","E16_original_full_sha256"),
 (ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json","E17A_original_joint_sha256"),
 (ROOT/"source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_original_joint_native_CLASS_transfer_resolution.json","E17B_original_joint_sha256"),
 (ARCH/"e17c_original_joint_restricted_retarded_streaming.json","E17C_original_joint_sha256"),
 (ARCH/"e17c_independent_scalar_source_only_replay.json","E17C_independent_sha256")]
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
      "minus":"Fminus_CLASS_normalized"}
KS=(.05,.075,.1);KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.);ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi);S0=(0.,.25,1.)
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d0_external_gravity_source"

def check(c,msg):
    if not c:raise ValueError("E17D0_INDEPENDENT_FAIL_CLOSED: "+msg)

def sha(b):return hashlib.sha256(b).hexdigest()

def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def gate():
    check(blob(PRE.read_bytes())==P_BLOB,"E17D0 prospective Git blob")
    p=json.loads(PRE.read_bytes());parents=p["immutable_parents"]
    check(sha(CSV.read_bytes())==parents["E8_original_4000q_CSV_sha256"],
          "frozen E8 original CSV SHA")
    for path,name in PARENTS:
        check(sha(path.read_bytes())==parents[name],"original parent "+name)
    check(blob((ARCH/"archive_manifest.json").read_bytes())==
          parents["E17C_archive_manifest_git_blob"],"E17C archived manifest")
    check(p["frozen"]["states"]==list(STATES)
          and p["frozen"]["tracer_order"]==["LRG","ELG"]
          and p["frozen"]["geometry_count"]==576
          and p["frozen"]["phase_reference_center"]==list(S0)
          and all(p["stop"].values()),"E17D0 scope or physical STOP")
    with CSV.open(newline="",encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh))
    check(len(rows)==4000,"original E8 q row count")
    q=[float(t["q_dimensionless"]) for t in rows]
    f={state:[float(t[COLS[state]]) for t in rows] for state in STATES}
    check(all(q[i]<q[i+1] for i in range(3999))
          and all(math.isfinite(y) and y>0 for arr in f.values() for y in arr),
          "original 4000q positivity and monotonicity")
    return p,q,f

def geometry():
    for k in KS:
        for K in KLS:
            r=K/k
            for ms in MS:
                for ml in ML:
                    for index,phi in enumerate(PH):
                        c=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(
                            max(0.,1.-ml*ml))*math.cos(phi)
                        yield (k,K,ms,ml,index,
                               k*math.sqrt(1.+r*r/4.-r*c),
                               k*math.sqrt(1.+r*r/4.+r*c))

def derivatives(x):
    """Analytic j0'=d(sin x/x)/dx and j0''; stable origin series."""
    if abs(x)<.02:
        x2=x*x
        return (-x/3.+x*x2/30.-x*x2*x2/840.+x*x2*x2*x2/45360.,
                -1./3.+x2/10.-x2*x2/168.+x2*x2*x2/6480.)
    sin=math.sin(x);cos=math.cos(x)
    return ((x*cos-sin)/(x*x),
            -sin/x-2.*cos/(x*x)+2.*sin/(x*x*x))

def prepare(q,f):
    """Independent GL2 node weights; denominator exact for piecewise linear F."""
    intervals=[]
    norm_terms=[]
    sqrt3=math.sqrt(3.)
    for i in range(len(q)-1):
        lo=q[i];hi=q[i+1];dx=hi-lo
        slope=(f[i+1]-f[i])/dx
        center=(lo+hi)/2.
        z=dx/(2.*sqrt3)
        n1=center-z;n2=center+z
        ff1=f[i]+slope*(n1-lo)
        ff2=f[i]+slope*(n2-lo)
        area=dx/2.
        intervals.append((n1,n2,ff1,ff2,area))
        norm_terms.append(area*(n1*n1*ff1+n2*n2*ff2))
    den=math.fsum(norm_terms)
    check(den>0,"finite q number-weighted integral positive")
    return intervals,den,(q[-1]*q[-1]*f[-1],q[0]*q[0]*f[0])

def ibp_shape(prepared,s,q):
    intervals,den,bd=prepared
    upper=bd[0]*derivatives(s*q[-1])[0]
    lower=bd[1]*derivatives(s*q[0])[0]
    def integrand(node,f):
        j,jj=derivatives(s*node)
        return f*(2.*node*j+s*node*node*jj)
    terms=(area*(integrand(n1,f1)+integrand(n2,f2))
           for n1,n2,f1,f2,area in intervals)
    numerator=-(upper-lower)+math.fsum(terms)
    return numerator,numerator/den

def verify(original,p,q,f,joint):
    check(joint["original_E16_geometries"]==576
          and joint["three_state_leg_phase_samples"]==10368
          and joint["full_physical_finiteK_bispectrum"]=="BLOCKED"
          and joint["observed_odd_SEALED"] is True,
          "original joint scope/size")
    geos=list(geometry())
    check(len(geos)==576,"independent analytic E16 geometry")
    max_gap=0.;n=0;reports={}
    for state in STATES:
        filename=original/("e17d0_external_potential_source_"+state+".json")
        raw=filename.read_bytes()
        check(sha(raw)==joint["original_state_reports"][state]["sha256"],
              "immutable original state SHA "+state)
        obj=json.loads(raw)
        check(obj["state"]==state and obj["n_legs_times_fixed_phases"]==3456
              and obj["full_retarded_halo_tracer_bispectrum"]=="BLOCKED"
              and obj["observed_odd_SEALED"] is True,"state or physical stop "+state)
        piece=prepare(q,f[state])
        den=piece[1]
        max_gap=max(max_gap,abs(den-obj["I_F_finite_original_q_support"])/max(1.,den))
        table=obj["original_geometries_both_legs"]
        check(len(table)==576,"original state orientation count")
        cache={}
        for g,row in zip(geos,table):
            k,K,ms,ml,index,m1,m2=g
            check((row["k"],row["K"],row["mu_s"],row["mu_L"],row["phi_index"])==
                  (k,K,ms,ml,index),"fixed original E16 geometry row order")
            for leg,m in enumerate((m1,m2)):
                gotm=row["leg_moduli_h_Mpc"][leg]
                max_gap=max(max_gap,abs(gotm-m)/max(1.,m))
                for j,s0 in enumerate(S0):
                    rec=row["per_leg_three_dimensionless_phases"][leg][j]
                    check(rec["s_center"]==s0,"prospective phase grid")
                    s=s0*m/.05
                    key=repr(s)
                    if key not in cache:
                        cache[key]=ibp_shape(piece,s,q)
                    J,H=cache[key]
                    max_gap=max(max_gap,abs(s-rec["s_leg"])/max(1.,abs(s)),
                                abs(J-rec["unnormalized_external_potential_impulse_shape_J"])/den,
                                abs(H-rec["H_normalized_unit_symbolic_potential_shape"]))
                    n+=1
        se=p["QA_before_new_numerics"]["small_s_test_positive"]
        sml=ibp_shape(piece,se,q)[1]/se
        max_gap=max(max_gap,abs(sml-obj["direct_small_s_H_over_s"])*se)
        reports[state]={"sha256":sha(raw),"samples":3456,
                        "I_F_denominator":den}
        print("E17D0_INDEPENDENT_STATE_REPLAY",state,3456,
              "MAX_SCALED_GAP_SO_FAR",max_gap,flush=True)
    check(n==10368,"original independent 10368 sample closure")
    check(max_gap<p["QA_before_new_numerics"][
          "independent_stdlib_integration_by_parts_scaled_tolerance"],
          "direct piecewise source vs derivative-free independent IBP")
    return max_gap,n,reports

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"existing independent audit changed")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d0i_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D0_INDEPENDENT_SOURCE_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--output-dir",type=Path,default=OUT)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    args=a.parse_args()
    p,q,f=gate()
    check(len(list(geometry()))==576,"independent E16 576 exact triangles")
    if args.self_test:
        check(derivatives(0.)==(0.,-1./3.),"independent j0 analytic origin")
        print("E17D0_INDEPENDENT_PARENT_AND_ANALYTIC_GATE_OK",flush=True)
        return
    src=args.output_dir/"e17d0_original_joint_external_potential_source_only.json"
    raw=src.read_bytes()
    joint=json.loads(raw)
    if args.negative_control:
        tampered=json.loads(raw)
        tampered["original_state_reports"]["FD"]["sha256"]="0"*64
        try:verify(args.output_dir,p,q,f,tampered)
        except ValueError as exc:
            check("immutable original state SHA" in str(exc),
                  "negative control rejected at unexpected gate")
            print("E17D0_TAMPERED_ORIGINAL_SHA_REJECTED",flush=True)
            return
        raise AssertionError("tampered source report accepted")
    gap,n,report=verify(args.output_dir,p,q,f,joint)
    cert={"date":"2026-09-27",
          "status":"E17D0_INDEPENDENT_STDLIB_DERIVATIVE_FREE_FINITE_Q_BOUNDARY_SOURCE_REPLAY_PASS_FULL_PHYSICAL_B_BLOCKED",
          "prospective_protocol_git_blob":P_BLOB,
          "original_joint_sha256":sha(raw),
          "original_E8_CSV_sha256":p["immutable_parents"]["E8_original_4000q_CSV_sha256"],
          "state_sha256_and_replay_counts":report,
          "total_original_E16_geometries":576,"total_state_leg_phase_samples":n,
          "max_scaled_independent_integration_by_parts_gap":gap,
          "method":"Independent pure-stdlib scalar E16 geometry and derivative-free integration by parts with finite original q0,qmax boundary, independently reconstructed GL2 piecewise-linear original F; no original E17D0 code import, CLASS or NumPy.",
          "full_halo_force_history_and_self_gravity_computed":False,
          "full_physical_finiteK_bispectrum":"BLOCKED",
          "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
          "main_untouched_PR_draft":True}
    save_once(args.output_dir/"e17d0_independent_finite_boundary_scalar_replay.json",cert)
    print("E17D0_INDEPENDENT_10368_FINITE_BOUNDARY_SOURCE_PASS",
          "MAX",gap,"FULL_PHYSICAL_B_BLOCKED",flush=True)

if __name__=="__main__":
    main()
