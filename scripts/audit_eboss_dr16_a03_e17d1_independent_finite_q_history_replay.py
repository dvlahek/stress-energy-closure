#!/usr/bin/env python3
"""Independent E17D1 finite-q scalar time-ordered source-history replay.

Only Python standard library. Original E8 CSV and analytic E16 geometry;
IBP of each original piecewise-linear F(q) including true qmax boundary.
No original E17D1 worker import, no NumPy, no CLASS or survey inputs.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e17d1_causal_external_history_nonidentifiability_prereg_2026-09-27.json"
P_BLOB="791502d5385e35b358fb2d2afebf5fac77bd6ad5"
CSV=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
CARCH=ROOT/"source_data/eboss_dr16_a03_e17c_archived_CI_2026_09_27"
DARCH=ROOT/"source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27"
DJOINT=DARCH/"e17d0_original_joint_external_potential_source_only.json"
DIND=DARCH/"e17d0_independent_finite_boundary_scalar_replay.json"
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
      "minus":"Fminus_CLASS_normalized"}
KS=(.05,.075,.1);KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.);ML=(-1.,0.,.5,1.);PH=(0.,math.pi/2.,math.pi)
S0=(0.,.25,1.)
HIST={"early":(0.,.25,4.),"late":(.75,1.,4.)}
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d1_causal_history"

def require(ok,why):
    if not ok:raise ValueError("E17D1_INDEPENDENT_FAIL_CLOSED: "+why)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def gate():
    require(blob(PRE.read_bytes())==P_BLOB,"prospective E17D1 protocol Git blob changed")
    p=json.loads(PRE.read_bytes())
    pp=p["immutable_parent_sha256"]
    require(sha(CSV.read_bytes())==pp["E8_original_4000q_CSV"]
            and sha(E16.read_bytes())==pp["E16_original_full"]
            and sha((CARCH/"e17c_original_joint_restricted_retarded_streaming.json").read_bytes())==
            pp["E17C_original_joint"]
            and sha(DJOINT.read_bytes())==pp["E17D0_original_joint"]
            and sha(DIND.read_bytes())==pp["E17D0_independent"],
            "original E8/E16/E17C/E17D0 immutable SHA")
    require(blob((DARCH/"archive_manifest.json").read_bytes())==
            p["immutable_parent_git_blob"]["E17D0_archive_manifest"],
            "original E17D0 archive manifest")
    require(all(p["absolute_STOP"].values())
            and p["frozen"]["states"]==list(STATES)
            and p["frozen"]["k_short_comoving_h_Mpc"]==list(KS)
            and p["frozen"]["K_long_comoving_h_Mpc"]==list(KLS)
            and p["frozen"]["mu_short"]==list(MS)
            and p["frozen"]["mu_long"]==list(ML)
            and p["frozen"]["center_dimensionless_phase"]==list(S0)
            and p["frozen"]["E16_geometries"]==576,
            "frozen science/observed STOP or geometry changed")
    for name,(a,b,g) in HIST.items():
        x=p["restricted_mathematical_object"]["source_history_"+name]
        require(x["u_support_inclusive"]==[a,b]
                and x["constant_amplitude"]==g and g*(b-a)==1.,
                "prospective unit area symbolic test histories changed")
    with CSV.open(newline="",encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh))
    require(len(rows)==4000,"original 4000q count")
    q=[float(x["q_dimensionless"]) for x in rows]
    f={state:[float(x[COLS[state]]) for x in rows] for state in STATES}
    require(all(q[i]<q[i+1] for i in range(len(q)-1))
            and all(math.isfinite(y) and y>0 for arr in f.values() for y in arr),
            "frozen original F or q grid drift")
    return p,q,f

def analytic_geometry():
    """Independent dot-product E16 magnitudes, no original worker module."""
    for k in KS:
        for K in KLS:
            ratio=K/k
            for ms in MS:
                for ml in ML:
                    for index,phi in enumerate(PH):
                        c=ms*ml+math.sqrt(max(0.,1-ms*ms))*math.sqrt(
                            max(0.,1-ml*ml))*math.cos(phi)
                        yield (k,K,ms,ml,index,
                               k*math.sqrt(1.+ratio*ratio/4.-ratio*c),
                               k*math.sqrt(1.+ratio*ratio/4.+ratio*c))

def sinc_and_prime(x):
    if abs(x)<.02:
        t=x*x
        return (1.-t/6.+t*t/120.-t*t*t/5040.+t*t*t*t/362880.,
                -x/3.+x*t/30.-x*t*t/840.+x*t*t*t/45360.)
    s=math.sin(x); c=math.cos(x)
    return s/x,(x*c-s)/(x*x)

def setup_quadrature(q,f):
    nodes=[];area=[];den=[]
    # Independent scalar GL2 on exact piecewise-linear original source.
    for i in range(len(q)-1):
        lo,hi=q[i],q[i+1]
        step=hi-lo
        s=(f[i+1]-f[i])/step
        mid=(lo+hi)*.5
        delta=step/(2.*math.sqrt(3.))
        for node in (mid-delta,mid+delta):
            fq=f[i]+s*(node-lo)
            a=.5*step
            nodes.append((node,fq,a))
            den.append(a*node*node*fq)
    normalization=math.fsum(den)
    require(normalization>0 and math.isfinite(normalization),
            "independent original finite q denominator")
    return nodes,normalization,q[0]*f[0],q[-1]*f[-1]

def ibp_history(nodes,den,q_lower_F,q_upper_F,sigma,u,history,qmin,qmax):
    a,b,g=HIST[history]
    if sigma==0. or u<=a:return 0.
    bprime=min(b,u)
    if bprime<=a:return 0.
    lo=u-bprime;hi=u-a
    def kernel(q):
        yhi=sigma*q*hi
        ylo=sigma*q*lo
        sh,jh=sinc_and_prime(yhi)
        sl,jl=sinc_and_prime(ylo)
        A=sh-sl
        Ap=sigma*(hi*jh-lo*jl)
        return A+q*Ap
    # Exact finite-original-q integration by parts:
    # R=g/(den*sigma) [ integral F(A+qA') - [q F A]_{qmin}^{qmax} ].
    vals=(weight*f*kernel(q) for q,f,weight in nodes)
    top_hi=sinc_and_prime(sigma*qmax*hi)[0]
    top_lo=sinc_and_prime(sigma*qmax*lo)[0]
    bot_hi=sinc_and_prime(sigma*qmin*hi)[0]
    bot_lo=sinc_and_prime(sigma*qmin*lo)[0]
    boundary=q_upper_F*(top_hi-top_lo)-q_lower_F*(bot_hi-bot_lo)
    return g*(math.fsum(vals)-boundary)/(den*sigma)

def full_replay(original,joint,p,q,fs):
    require(joint["original_E16_geometries"]==576
            and joint["original_state_leg_phase_samples"]==10368
            and joint["two_histories_per_sample"] is True
            and joint["full_physical_finiteK_bispectrum"]=="BLOCKED"
            and joint["observed_odd_SEALED"] is True,
            "original joint source-only scope or count drift")
    geometries=list(analytic_geometry())
    require(len(geometries)==576,"independent original E16 geometry count")
    maximum=0.;n=0;per_state={}
    for state in STATES:
        path=original/("e17d1_causal_history_"+state+".json")
        raw=path.read_bytes()
        require(sha(raw)==joint["original_state_reports"][state]["sha256"],
                "original E17D1 state SHA mismatch: "+state)
        item=json.loads(raw)
        rows=item["original_geometries_both_legs"]
        require(item["state"]==state and len(rows)==576
                and item["original_state_leg_phase_values"]==3456
                and item["full_physical_finiteK_bispectrum"]=="BLOCKED"
                and item["observed_odd_SEALED"] is True,
                "original state geometry/physical STOP")
        nodes,den,qmin_F,qmax_F=setup_quadrature(q,fs[state])
        maximum=max(maximum,abs(den-item["original_I_F"])/max(1.,den))
        cache={}
        def replay(sigma,name,u=1.):
            key=(repr(sigma),name,u)
            if key not in cache:
                cache[key]=ibp_history(nodes,den,qmin_F,qmax_F,
                                       sigma,u,name,q[0],q[-1])
            return cache[key]
        for expect,got in zip(geometries,rows):
            k,K,ms,ml,index,k1,k2=expect
            require((got["k"],got["K"],got["mu_s"],got["mu_L"],
                     got["phi_index"])==(k,K,ms,ml,index),
                    "original geometry labels or sample order")
            for leg,mod in enumerate((k1,k2)):
                max_mod=got["leg_moduli_h_Mpc"][leg]
                maximum=max(maximum,abs(max_mod-mod)/max(1.,abs(mod)))
                for phase,s0 in enumerate(S0):
                    row=got["per_leg_three_dimensionless_phases"][leg][phase]
                    require(row["s_center"]==s0,"prospectively frozen phases")
                    sigma=s0*mod/.05
                    maximum=max(maximum,abs(sigma-row["sigma_leg"])/
                                max(1.,abs(sigma)))
                    for name in ("early","late"):
                        computed=replay(sigma,name)
                        stored=row["R_"+name+"_u1_unit_area"]
                        maximum=max(maximum,abs(computed-stored)/
                                    max(1.,abs(computed),abs(stored)))
                    pre=replay(sigma,"late",.5)
                    maximum=max(maximum,abs(pre-row[
                        "R_late_u_half_causal_control"]))
                    require(pre==0.,"future forcing leaks into u=.5 retarded output")
                    n+=1
        # Compare an independently computed original fixed witness, not a fitted point.
        if state=="FD":
            k,K,ms,ml,index,mod,_=geometries[0]
            sigma=.25*mod/.05
            witness=abs(replay(sigma,"early")-replay(sigma,"late"))
            maximum=max(maximum,abs(witness-joint["fixed_FD_history_witness_abs"]))
            require(witness>p["pre_registered_QA"][
                "FD_first_leg_fixed_history_difference_min_abs"],
                "fixed frozen FD counterexample does not survive independent replay")
        per_state[state]={"original_sha256":sha(raw),
                          "replayed_original_state_leg_phase_samples":len(rows)*6,
                          "finite_q_normalization":den}
        print("E17D1_INDEPENDENT_FULL_STATE_REPLAY",state,3456,
              "MAX_SCALED_GAP",maximum,flush=True)
    require(n==10368,"original source-only state leg phase sample count")
    require(maximum<p["pre_registered_QA"][
            "independent_full_stdlib_finite_q_IBP_scaled_gap_max"],
            "full independent derivative-free finite-q source audit mismatch")
    return maximum,n,per_state

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"independent certificate immutability")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d1i_",delete=False) as fh:
            temp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17D1_INDEPENDENT_JSON",path,"SHA256",sha(raw),flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--output-dir",type=Path,default=OUT)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--negative-control",action="store_true")
    args=a.parse_args()
    p,q,f=gate()
    require(len(list(analytic_geometry()))==576,"independent E16 576 geometry gate")
    require(sinc_and_prime(0.)==(1.,0.),"analytic sinc origin")
    if args.self_test:
        print("E17D1_INDEPENDENT_ORIGINAL_SHA_AND_ANALYTIC_GEOMETRY_GATE_PASS",
              "FROZEN_Q",len(q),flush=True)
        return
    path=args.output_dir/"e17d1_original_joint_causal_history_nonidentifiability.json"
    raw=path.read_bytes()
    joint=json.loads(raw)
    if args.negative_control:
        bad=json.loads(raw)
        bad["original_state_reports"]["FD"]["sha256"]="0"*64
        try:full_replay(args.output_dir,bad,p,q,f)
        except ValueError as exc:
            require("original E17D1 state SHA mismatch" in str(exc),
                    "negative control rejected at unexpected gate")
            print("E17D1_INDEPENDENT_TAMPERED_ORIGINAL_SHA_REJECTED",flush=True)
            return
        raise AssertionError("tampered original JSON SHA was accepted")
    maximum,n,per_state=full_replay(args.output_dir,joint,p,q,f)
    output={"date":"2026-09-27",
            "status":"E17D1_INDEPENDENT_STDLIB_FULL_FINITE_Q_CAUSAL_HISTORY_REPLAY_PASS_PHYSICAL_HALO_B_BLOCKED",
            "prospective_E17D1_protocol_git_blob":P_BLOB,
            "original_joint_sha256":sha(raw),
            "original_E8_CSV_sha256":p["immutable_parent_sha256"][
                "E8_original_4000q_CSV"],
            "state_original_SHA_and_samples":per_state,
            "original_E16_geometries":576,
            "replayed_state_leg_phase_samples":n,
            "replayed_two_unit_area_symbolic_histories_per_sample":True,
            "maximum_scaled_independent_finite_q_IBP_history_gap":maximum,
            "independent_method":"Pure stdlib, independently reconstructed E16 analytic leg moduli and original 4000q piecewise-linear F; derivative-free q integration by parts with full original qmin/qmax boundary, exact top-hat time integral; no original E17D1 runner import, NumPy, CLASS or eBOSS observed inputs.",
            "physical_halo_time_history_supplied":False,
            "full_physical_finiteK_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(args.output_dir/"e17d1_independent_full_finite_q_history_replay.json",output)
    print("E17D1_INDEPENDENT_10368_TWICE_HISTORY_SOURCE_ONLY_PASS",
          "MAX",maximum,"PHYSICAL_HALO_B_BLOCKED",flush=True)

if __name__=="__main__":
    main()
