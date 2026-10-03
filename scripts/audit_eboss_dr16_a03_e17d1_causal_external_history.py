#!/usr/bin/env python3
"""E17D1 restricted retarded external-force history test, NOT a halo B.

For each unchanged E16 triangle and E8 F(q), convolve the original frozen
E17D0 unit-symbolic-potential susceptibility with two preregistered unit-area
external test histories. Time u and phases are DIMENSIONLESS TEST VARIABLES,
not physical halo formation, cosmic conformal time, LRG or ELG properties.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import sys
import tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e17d0_external_potential_vlasov_source as e17d0

PRE=ROOT/"source_data/eboss_dr16_a03_e17d1_causal_external_history_nonidentifiability_prereg_2026-09-27.json"
PRE_BLOB="791502d5385e35b358fb2d2afebf5fac77bd6ad5"
E17D0_RUNNER_BLOB="f8ed0cfcea70251c005a1dd1bb7be403984e9854"
E17D0_ARCH=ROOT/"source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27"
E17D0_JOINT=E17D0_ARCH/"e17d0_original_joint_external_potential_source_only.json"
E17D0_IND=E17D0_ARCH/"e17d0_independent_finite_boundary_scalar_replay.json"
E17D0_MAN=E17D0_ARCH/"archive_manifest.json"
STATES=("FD","plus","minus")
S0=(0.,.25,1.)
HISTORIES={"early":(0.,.25,4.),"late":(.75,1.,4.)}
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d1_causal_history"

def require(c,msg):
    if not c:raise ValueError("E17D1_FAIL_CLOSED: "+msg)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"existing original report changed")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d1_",delete=False) as fh:
            temp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(temp,path)
        finally:temp.unlink(missing_ok=True)
    print("E17D1_ORIGINAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"E17D1 prospective protocol Git blob changed")
    p=json.loads(PRE.read_bytes())
    require(blob(Path(e17d0.__file__).read_bytes())==E17D0_RUNNER_BLOB,
            "original E17D0 physical runner changed")
    pp=p["immutable_parent_git_blob"]
    require(pp["E17D0_original_runner"]==E17D0_RUNNER_BLOB
            and blob(E17D0_MAN.read_bytes())==pp["E17D0_archive_manifest"]
            and blob(e17d0.PRE.read_bytes())==pp["E17D0_protocol"],
            "E17D0 immutable original runner/protocol/manifest blob changed")
    q=p["immutable_parent_sha256"]
    for path,key in ((e17d0.CSV,"E8_original_4000q_CSV"),
                     (e17d0.E16,"E16_original_full"),
                     (e17d0.E17C,"E17C_original_joint"),
                     (E17D0_JOINT,"E17D0_original_joint"),
                     (E17D0_IND,"E17D0_independent")):
        require(sha(path.read_bytes())==q[key],"E17D1 immutable parent SHA: "+key)
    require(all(p["absolute_STOP"].values())
            and p["frozen"]["states"]==list(STATES)
            and p["frozen"]["tracer_order"]==["LRG","ELG"]
            and p["frozen"]["E16_geometries"]==576
            and p["frozen"]["center_dimensionless_phase"]==list(S0)
            and p["frozen"]["mu_short"]==list(e17d0.MS)
            and p["frozen"]["mu_long"]==list(e17d0.ML)
            and p["frozen"]["k_short_comoving_h_Mpc"]==list(e17d0.KS)
            and p["frozen"]["K_long_comoving_h_Mpc"]==list(e17d0.KLS),
            "prospective frozen geometry or hard science STOP")
    for name,(a,b,amp) in HISTORIES.items():
        x=p["restricted_mathematical_object"]["source_history_"+name]
        require(x["u_support_inclusive"]==[a,b] and x["constant_amplitude"]==amp
                and amp*(b-a)==1.,"preregistered equal-unit-area source histories changed")
    require(p["restricted_mathematical_object"]["main_report_observation_u"]==1
            and p["restricted_mathematical_object"]["dimensionless_observation_u"]==[.5,1]
            and p["pre_registered_QA"]["no_posthoc_tolerance_retuning"],
            "fixed dimensionless observation or original QA drift")
    old=json.loads(E17D0_JOINT.read_bytes())
    require(old["three_state_leg_phase_samples"]==10368
            and old["full_physical_finiteK_bispectrum"]=="BLOCKED"
            and old["observed_odd_SEALED"] is True,
            "E17D0 original source-only science scope drift")
    return p

def stable_sinc(x):
    return np.sinc(x/np.pi)

def window_response(prepared,sigma,u,window):
    a,b,amplitude=window
    if sigma==0. or u<=a:return 0.
    bprime=min(b,u)
    if bprime<=a:return 0.
    lag_lo=u-bprime
    lag_hi=u-a
    q1,q2,slope,area,den,_=prepared
    A1=stable_sinc(sigma*q1*lag_hi)-stable_sinc(sigma*q1*lag_lo)
    A2=stable_sinc(sigma*q2*lag_hi)-stable_sinc(sigma*q2*lag_lo)
    numerator=float(np.sum(area*slope*(q1*A1+q2*A2)))
    return -amplitude*numerator/(den*sigma)

def independent_time_gauss_probe(prepared,sigma,window,u=1.,n=48):
    a,b,amplitude=window
    if sigma==0. or u<=a:return 0.
    b=min(b,u)
    if b<=a:return 0.
    nodes,weights=np.polynomial.legendre.leggauss(n)
    center=(a+b)/2.
    half=(b-a)/2.
    result=0.
    for x,w in zip(nodes,weights):
        eval_s=sigma*(u-(center+half*x))
        result+=w*e17d0.source(prepared,eval_s)[1]
    return amplitude*half*result

def run_state(state,out):
    p=gate()
    require(state in STATES,"unexpected neutrino state")
    p0,q,fs,sym=e17d0.gate()
    prepared=e17d0.prepare(q,fs[state])
    rows=e17d0.geometry()
    require(len(rows)==576,"original 576 E16 triangles")
    qv=p["pre_registered_QA"]
    cache={}
    def response(sigma,window,u=1.):
        key=(repr(sigma),window,repr(u))
        if key not in cache:
            cache[key]=window_response(prepared,sigma,u,HISTORIES[window])
        return cache[key]
    qa={"max_s0_zero_abs":0.,"max_signed_phase_oddness_abs":0.,
        "max_late_u_half_causal_abs":0.,"max_short_leg_exchange_abs":0.,
        "max_original_triangle_closure_h_Mpc":0.,
        "maximum_time_Gauss48_spot_scaled_gap":0.}
    output=[]
    for row in rows:
        leg_vals=[]
        qa["max_original_triangle_closure_h_Mpc"]=max(
            qa["max_original_triangle_closure_h_Mpc"],row["closure_h_Mpc"])
        for leg in row["leg_moduli_h_Mpc"]:
            ph=[]
            for s0 in S0:
                sigma=s0*leg/.05
                early=response(sigma,"early")
                late=response(sigma,"late")
                late_half=response(sigma,"late",.5)
                qa["max_late_u_half_causal_abs"]=max(
                    qa["max_late_u_half_causal_abs"],abs(late_half))
                if s0==0.:
                    qa["max_s0_zero_abs"]=max(
                        qa["max_s0_zero_abs"],abs(early),abs(late))
                for name,actual in (("early",early),("late",late)):
                    qa["max_signed_phase_oddness_abs"]=max(
                        qa["max_signed_phase_oddness_abs"],
                        abs(actual+response(-sigma,name)))
                ph.append({"s_center":s0,"sigma_leg":sigma,
                           "R_early_u1_unit_area":early,
                           "R_late_u1_unit_area":late,
                           "R_late_u_half_causal_control":late_half})
            leg_vals.append(ph)
        row["per_leg_three_dimensionless_phases"]=leg_vals
        output.append(row)
    require(qa["max_s0_zero_abs"]<qv["sigma_zero_both_history_abs_max"]
            and qa["max_signed_phase_oddness_abs"]<
            qv["signed_sigma_oddness_abs_max"]
            and qa["max_late_u_half_causal_abs"]<
            qv["late_history_causal_zero_at_u_half_abs_max"],
            "source-only causal control or odd response violated")
    require(qa["max_original_triangle_closure_h_Mpc"]<
            qv["original_short_leg_geometry_closure_abs_h_Mpc_max"],
            "original E16 closed triangle control")
    mp={(r["k"],r["K"],r["mu_s"],r["mu_L"],r["phi_index"]):r
        for r in output}
    for row in output:
        other=mp.get((row["k"],row["K"],-row["mu_s"],row["mu_L"],
                      2-row["phi_index"]))
        if other is None:
            ms,ml,phi=row["mu_s"],row["mu_L"],row["phi_rad"]
            cos=ms*ml+math.sqrt(max(0.,1-ms*ms))*math.sqrt(
                max(0.,1-ml*ml))*math.cos(phi)
            ratio=row["K"]/row["k"]
            reflected=(row["k"]*math.sqrt(1.+ratio*ratio/4.+ratio*cos),
                       row["k"]*math.sqrt(1.+ratio*ratio/4.-ratio*cos))
            for j in range(2):
                qa["max_short_leg_exchange_abs"]=max(
                    qa["max_short_leg_exchange_abs"],
                    abs(row["leg_moduli_h_Mpc"][j]-reflected[1-j]))
        else:
            for j in range(2):
                qa["max_short_leg_exchange_abs"]=max(
                    qa["max_short_leg_exchange_abs"],
                    abs(row["leg_moduli_h_Mpc"][j]-other["leg_moduli_h_Mpc"][1-j]))
                for i in range(3):
                    a=row["per_leg_three_dimensionless_phases"][j][i]
                    b=other["per_leg_three_dimensionless_phases"][1-j][i]
                    for name in ("early","late"):
                        key="R_"+name+"_u1_unit_area"
                        qa["max_short_leg_exchange_abs"]=max(
                            qa["max_short_leg_exchange_abs"],abs(a[key]-b[key]))
    require(qa["max_short_leg_exchange_abs"]<qv["short_leg_exchange_abs_max"],
            "E16 exact short-leg reversal/analytic unsampled orientation")
    for row in (output[0],output[-1]):
        for j,mod in enumerate(row["leg_moduli_h_Mpc"]):
            for s0 in (.25,1.):
                sigma=s0*mod/.05
                for name,win in HISTORIES.items():
                    result=response(sigma,name)
                    numerical=independent_time_gauss_probe(prepared,sigma,win)
                    gap=abs(result-numerical)/max(1.,abs(result),abs(numerical))
                    qa["maximum_time_Gauss48_spot_scaled_gap"]=max(
                        qa["maximum_time_Gauss48_spot_scaled_gap"],gap)
    require(qa["maximum_time_Gauss48_spot_scaled_gap"]<
            qv["time_quadrature_spot_relative_gap_max"],
            "original exact-time primitive vs independent Gauss48 Duhamel convolution")
    witness=output[0]["per_leg_three_dimensionless_phases"][0][1]
    difference=abs(witness["R_early_u1_unit_area"]-witness[
                   "R_late_u1_unit_area"])
    if state=="FD":
        require(difference>qv["FD_first_leg_fixed_history_difference_min_abs"],
                "preregistered fixed FD restricted-model source-history witness absent")
    data={"date":"2026-09-27",
          "status":"E17D1_ORIGINAL_CAUSAL_EQUAL_AREA_EXTERNAL_HISTORY_RESTRICTED_MODEL_ONLY_FULL_B_BLOCKED",
          "state":state,"prospective_protocol_git_blob":PRE_BLOB,
          "original_E17D0_joint_sha256":p["immutable_parent_sha256"][
              "E17D0_original_joint"],
          "same_original_E8_E16_E17C_E17D0_without_physical_halo_force":True,
          "original_q_support":[float(q[0]),float(q[-1])],
          "original_piecewise_linear_q_intervals":len(q)-1,
          "original_I_F":prepared[4],
          "unit_area_early_late_histories":True,
          "histories_are_symbolic_dimensionless_external_forcing_NOT_halo_potentials":True,
          "u_observation_dimensionless":1.,
          "QA":qa,
          "fixed_first_geometry_first_leg_phase_quarter_difference_abs":difference,
          "original_geometries_both_legs":output,
          "original_state_leg_phase_values":len(output)*2*3,
          "physically_specified_real_halo_source_or_Frw_cosmic_history":False,
          "full_physical_finiteK_bispectrum":"BLOCKED",
          "observed_odd_SEALED":True,"no_new_CLASS_FITS_mock":True,
          "main_untouched_PR_draft":True}
    save_once(out/("e17d1_causal_history_"+state+".json"),data)
    print("E17D1_ORIGINAL_STATE",state,"SAMPLES",data["original_state_leg_phase_values"],
          "WITNESS",difference,"TIME_GAUSS48_GAP",
          qa["maximum_time_Gauss48_spot_scaled_gap"],
          "FULL_PHYSICAL_B_BLOCKED",flush=True)

def aggregate(out):
    p=gate()
    reports={}
    for state in STATES:
        path=out/("e17d1_causal_history_"+state+".json")
        raw=path.read_bytes()
        o=json.loads(raw)
        require(o["state"]==state and o["original_state_leg_phase_values"]==3456
                and len(o["original_geometries_both_legs"])==576
                and o["full_physical_finiteK_bispectrum"]=="BLOCKED"
                and o["observed_odd_SEALED"] is True,
                "original E17D1 state report or science STOP")
        reports[state]={"sha256":sha(raw),"size_bytes":len(raw),
                        "QA":o["QA"],
                        "fixed_first_leg_history_witness":o[
                            "fixed_first_geometry_first_leg_phase_quarter_difference_abs"]}
    result={"date":"2026-09-27",
            "status":"E17D1_RESTRICTED_CAUSAL_UNIT_AREA_HISTORY_DEPENDENCE_CERTIFIED_FULL_PHYSICAL_HALO_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_parent_sha256":p["immutable_parent_sha256"],
            "original_state_reports":reports,
            "original_E16_geometries":576,
            "original_state_leg_phase_samples":10368,
            "two_histories_per_sample":True,
            "fixed_FD_history_witness_abs":reports["FD"][
                "fixed_first_leg_history_witness"],
            "identifiability_statement":"Under original fixed-epoch linearized Vlasov response and frozen E8–E17D0 source, two specified equal-area positive symbolic forcing histories give distinct retarded responses. Therefore these original sources alone do not determine a unique forced response without a source time history; this is NOT a proof about fully specified self-consistent physical halo models.",
            "real_halo_potential_history_or_HOD_calibrated":False,
            "full_physical_finiteK_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"no_new_CLASS_FITS_mock":True,
            "main_untouched_PR_draft":True}
    save_once(out/"e17d1_original_joint_causal_history_nonidentifiability.json",result)
    print("E17D1_ORIGINAL_SOURCE_ONLY_CAUSAL_HISTORY_WITNESS_PASS",
          result["fixed_FD_history_witness_abs"],"FULL_PHYSICAL_B_BLOCKED",flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--state",choices=STATES)
    g.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    args=a.parse_args()
    if args.self_test:
        p=gate()
        _,q,fs,sym=e17d0.gate()
        require(len(e17d0.geometry())==576 and sym<
                p["pre_registered_QA"][
                    "original_Fplus_Fminus_pointwise_average_scaled_gap_max"],
                "original physical parent/geometry or E8 pair average")
        print("E17D1_ORIGINAL_PROSPECTIVE_PARENT_AND_SOURCE_GATE_PASS",
              "FROZEN_Q",len(q),"PAIR",sym,"HALO_B_BLOCKED",flush=True)
    elif args.state:run_state(args.state,args.output_dir)
    else:aggregate(args.output_dir)

if __name__=="__main__":
    main()
