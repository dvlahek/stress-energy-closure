#!/usr/bin/env python3
"""E17D0 frozen external-gravity susceptibility, NOT retarded physical halo/tracer B.

Fixed-epoch, leading nonrelativistic, linearized Newtonian Vlasov impulse
per unit *symbolic* external potential. Frozen E8 F(q) piecewise-linear slopes,
Gauss-Legendre 2 quadrature per original interval, exact original E16 legs.
No halo potential/time history, GR completion, second-order force, HOD or data.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
import tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d0_linearized_vlasov_external_potential_source_prereg_2026-09-27.json"
PRE_BLOB="66df071f0750201e792ed2b3ffc81e28b13dfce6"
ARCH=ROOT/"source_data/eboss_dr16_a03_e17c_archived_CI_2026_09_27"
CSV=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17B=ROOT/"source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_original_joint_native_CLASS_transfer_resolution.json"
E17C=ARCH/"e17c_original_joint_restricted_retarded_streaming.json"
E17CI=ARCH/"e17c_independent_scalar_source_only_replay.json"
MAN=ARCH/"archive_manifest.json"
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
      "minus":"Fminus_CLASS_normalized"}
KS=(.05,.075,.1); KLS=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.); ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi); S0=(0.,.25,1.)
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d0_external_gravity_source"

def require(ok,why):
    if not ok:
        raise ValueError("E17D0_FAIL_CLOSED: "+why)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def save_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,"output differs; immutable archive: "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d0_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D0_ORIGINAL_SOURCE_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective original protocol Git blob")
    p=json.loads(PRE.read_bytes());parent=p["immutable_parents"]
    paths=((CSV,"E8_original_4000q_CSV_sha256"),
           (E16,"E16_original_full_sha256"),
           (E17A,"E17A_original_joint_sha256"),
           (E17B,"E17B_original_joint_sha256"),
           (E17C,"E17C_original_joint_sha256"),
           (E17CI,"E17C_independent_sha256"))
    for path,key in paths:
        require(sha(path.read_bytes())==parent[key],
                "original E8–E17C SHA parent changed: "+key)
    require(blob(MAN.read_bytes())==parent["E17C_archive_manifest_git_blob"],
            "original E17C archive manifest changed")
    old=json.loads(E17C.read_bytes())
    require(old["original_E16_geometries"]==576
            and old["three_state_leg_phase_samples"]==10368
            and old["full_physical_bispectrum"]=="BLOCKED"
            and old["observed_odd_SEALED"] is True,
            "E17C source-only parent scope changed")
    require(json.loads(E16.read_bytes())["QA"]["original_72_cases"]==72
            and json.loads(E17B.read_bytes())["physical_finite_K_bispectrum"]=="BLOCKED",
            "original fixed source or physical STOP changed")
    ff=p["frozen"]
    require(ff["states"]==list(STATES) and ff["tracer_order"]==["LRG","ELG"]
            and ff["short_k_comoving_h_Mpc"]==list(KS)
            and ff["long_K_comoving_h_Mpc"]==list(KLS)
            and ff["mu_short"]==list(MS) and ff["mu_long"]==list(ML)
            and ff["geometry_count"]==576
            and ff["phase_reference_center"]==list(S0)
            and ff["E8_q_samples"]==4000
            and all(p["stop"].values()),
            "registered E17D0 scope or STOP changed")
    a=np.genfromtxt(CSV,delimiter=",",names=True)
    require(a.shape==(4000,) and set(COLS.values()).issubset(a.dtype.names or ()),
            "original 4000q schema")
    q=np.asarray(a["q_dimensionless"],float)
    f={s:np.asarray(a[col],float) for s,col in COLS.items()}
    require(np.all(np.isfinite(q)) and np.all(np.diff(q)>0)
            and all(np.all(v>0) and np.all(np.isfinite(v)) for v in f.values()),
            "original q monotonicity or F positivity")
    symmetry=float(np.max(abs((f["plus"]+f["minus"])*.5-f["FD"]))/
                   np.max(abs(f["FD"])))
    require(symmetry<p["QA_before_new_numerics"][
            "Fplus_minus_pointwise_average_scaled_tolerance"],
            "frozen original pair average altered")
    return p,q,f,symmetry

def geometry():
    rows=[]
    for k in KS:
        for K in KLS:
            for ms in MS:
                for ml in ML:
                    for idx,phi in enumerate(PH):
                        u=(math.sqrt(max(0.,1.-ms*ms)),0.,ms)
                        v=(math.sqrt(max(0.,1.-ml*ml))*math.cos(phi),
                           math.sqrt(max(0.,1.-ml*ml))*math.sin(phi),ml)
                        lv=tuple(K*x for x in v)
                        x=tuple(k*u[i]-.5*lv[i] for i in range(3))
                        y=tuple(-k*u[i]-.5*lv[i] for i in range(3))
                        m1=math.sqrt(sum(xx*xx for xx in x))
                        m2=math.sqrt(sum(xx*xx for xx in y))
                        closure=math.sqrt(sum((x[i]+y[i]+lv[i])**2 for i in range(3)))
                        require(closure<1e-13,"E16 original triangle closure")
                        rows.append({"k":k,"K":K,"mu_s":ms,"mu_L":ml,
                                     "phi_index":idx,"phi_rad":phi,
                                     "leg_moduli_h_Mpc":[m1,m2],
                                     "closure_h_Mpc":closure})
    require(len(rows)==576,"E16 original 48×12 orientation drift")
    return rows

def jp(x):
    a=np.asarray(x,float)
    out=np.empty_like(a)
    small=np.abs(a)<.02
    t=a[small]
    out[small]=-t/3.+t**3/30.-t**5/840.+t**7/45360.
    t=a[~small]
    out[~small]=(t*np.cos(t)-np.sin(t))/(t*t)
    return out

def prepare(q,f):
    lo=q[:-1]; hi=q[1:]; step=hi-lo
    slope=(f[1:]-f[:-1])/step
    mid=(lo+hi)/2.
    d=step/(2.*math.sqrt(3.))
    qn1=mid-d;qn2=mid+d
    fl1=f[:-1]+slope*(qn1-lo)
    fl2=f[:-1]+slope*(qn2-lo)
    area=step/2.
    den=float(np.sum(area*(qn1**2*fl1+qn2**2*fl2)))
    require(den>0 and math.isfinite(den),"original q finite-support I_F")
    boundary=float(q[-1]**3*f[-1]-q[0]**3*f[0])
    return (qn1,qn2,slope,area,den,boundary)

def source(prepared,s):
    qn1,qn2,slope,area,den,_=prepared
    j=float(-np.sum(area*slope*(qn1**2*jp(s*qn1)+qn2**2*jp(s*qn2))))
    return j,j/den

def run_state(state,out):
    p,q,fs,sym=gate()
    require(state in STATES,"only original FD/plus/minus accepted")
    f=fs[state]
    prep=prepare(q,f)
    den=prep[4];boundary=prep[5]
    cache={}
    def calc(s):
        key=repr(s)
        if key not in cache:cache[key]=source(prep,s)
        return cache[key]
    records=[]
    qa={"max_zero_lag_abs":0.,"max_odd_s_gap":0.,
        "max_triangle_closure_h_Mpc":0.,"max_short_leg_exchange_abs":0.,
        "max_squeezed_H_abs_gap":0.}
    for row in geometry():
        leg_values=[]
        qa["max_triangle_closure_h_Mpc"]=max(
            qa["max_triangle_closure_h_Mpc"],row["closure_h_Mpc"])
        for mod in row["leg_moduli_h_Mpc"]:
            one=[]
            for s0 in S0:
                s=s0*mod/.05
                j,H=calc(s)
                jneg,Hneg=calc(-s)
                qa["max_odd_s_gap"]=max(qa["max_odd_s_gap"],abs(H+Hneg))
                if s0==0.:
                    qa["max_zero_lag_abs"]=max(qa["max_zero_lag_abs"],abs(H))
                one.append({"s_center":s0,"s_leg":s,
                            "unnormalized_external_potential_impulse_shape_J":j,
                            "H_normalized_unit_symbolic_potential_shape":H})
            leg_values.append(one)
        row["per_leg_three_dimensionless_phases"]=leg_values
        records.append(row)
    require(qa["max_zero_lag_abs"]<p["QA_before_new_numerics"]["H_zero_lag_abs"]
            and qa["max_odd_s_gap"]<p["QA_before_new_numerics"]["odd_s_reversal_abs"],
            "unit impulse zero-lag or odd time phase failed")
    mapping={(r["k"],r["K"],r["mu_s"],r["mu_L"],r["phi_index"]):r for r in records}
    for row in records:
        other=mapping.get((row["k"],row["K"],-row["mu_s"],row["mu_L"],
                           2-row["phi_index"]))
        if other is None:
            c=(row["mu_s"]*row["mu_L"]+
               math.sqrt(max(0.,1.-row["mu_s"]**2))*
               math.sqrt(max(0.,1.-row["mu_L"]**2))*math.cos(row["phi_rad"]))
            r=row["K"]/row["k"]
            reflected=(row["k"]*math.sqrt(1.+r*r/4.+r*c),
                       row["k"]*math.sqrt(1.+r*r/4.-r*c))
            for i in range(2):
                qa["max_short_leg_exchange_abs"]=max(
                    qa["max_short_leg_exchange_abs"],
                    abs(row["leg_moduli_h_Mpc"][i]-reflected[1-i]))
        else:
            for i in range(2):
                qa["max_short_leg_exchange_abs"]=max(
                    qa["max_short_leg_exchange_abs"],
                    abs(row["leg_moduli_h_Mpc"][i]-
                        other["leg_moduli_h_Mpc"][1-i]))
                for n in range(3):
                    qa["max_short_leg_exchange_abs"]=max(
                        qa["max_short_leg_exchange_abs"],
                        abs(row["per_leg_three_dimensionless_phases"][i][n][
                            "H_normalized_unit_symbolic_potential_shape"]-
                            other["per_leg_three_dimensionless_phases"][1-i][n][
                            "H_normalized_unit_symbolic_potential_shape"]))
    require(qa["max_short_leg_exchange_abs"]<
            p["QA_before_new_numerics"]["short_leg_exchange_abs"],
            "short-leg exchange with original E16 frozen orientation")
    eps=p["QA_before_new_numerics"]["squeezed_test_epsilon_K_over_k"]
    for k in KS:
        for s0 in S0:
            center=calc(s0*k/.05)[1]
            for mult in (-.5,+.5):
                x=s0*(k*(1.+mult*eps))/.05
                qa["max_squeezed_H_abs_gap"]=max(
                    qa["max_squeezed_H_abs_gap"],abs(calc(x)[1]-center))
    require(qa["max_squeezed_H_abs_gap"]<
            p["QA_before_new_numerics"]["squeezed_H_abs_gap"],
            "mathematical center K/k->0 test")
    se=p["QA_before_new_numerics"]["small_s_test_positive"]
    hs=calc(se)[1]
    predicted=-1.+boundary/(3.*den)
    small_gap=abs(hs/se-predicted)/max(1.,abs(predicted))
    qa["small_s_limit_relative_gap"]=small_gap
    require(small_gap<p["QA_before_new_numerics"][
            "small_s_limit_relative_tolerance"],"finite q-boundary small-s identity")
    outj={"date":"2026-09-27",
          "status":"E17D0_ORIGINAL_FROZEN_EXTERNAL_POTENTIAL_VLASOV_IMPULSE_SHAPE_ONLY_FULL_B_BLOCKED",
          "state":state,"prospective_protocol_git_blob":PRE_BLOB,
          "original_frozen_E8_CSV_sha256":p["immutable_parents"]["E8_original_4000q_CSV_sha256"],
          "E17C_parent_joint_sha256":p["immutable_parents"]["E17C_original_joint_sha256"],
          "q_support_original":[float(q[0]),float(q[-1])],
          "piecewise_linear_finite_q_intervals":len(q)-1,
          "source_pointwise_pair_average_scaled_gap":sym,
          "I_F_finite_original_q_support":den,
          "finite_support_q3F_upper_minus_lower":boundary,
          "analytic_small_s_H_over_s_limit_with_boundary":predicted,
          "direct_small_s_H_over_s":hs/se,
          "QA":qa,"original_geometries_both_legs":records,
          "n_legs_times_fixed_phases":len(records)*2*3,
          "actual_external_halo_potential_history_or_FRW_time_supplied":False,
          "full_retarded_halo_tracer_bispectrum":"BLOCKED",
          "observed_odd_SEALED":True,"new_CLASS_or_FITS_or_mock":False,
          "main_untouched_PR_draft":True}
    save_once(out/("e17d0_external_potential_source_"+state+".json"),outj)
    print("E17D0_ORIGINAL_STATE",state,"N",outj["n_legs_times_fixed_phases"],
          "SMALL_S_GAP",small_gap,"SQUEEZED",qa["max_squeezed_H_abs_gap"],
          "FULL_PHYSICAL_B_BLOCKED",flush=True)

def aggregate(out):
    p,q,fs,sym=gate()
    obj={}
    for state in STATES:
        file=out/("e17d0_external_potential_source_"+state+".json")
        raw=file.read_bytes()
        row=json.loads(raw)
        require(row["state"]==state and row["n_legs_times_fixed_phases"]==3456
                and row["full_retarded_halo_tracer_bispectrum"]=="BLOCKED"
                and row["observed_odd_SEALED"] is True,
                "state report or physical scope drift")
        obj[state]={"sha256":sha(raw),"bytes":len(raw),"QA":row["QA"]}
    reports={s:json.loads((out/("e17d0_external_potential_source_"+s+".json")).read_bytes())
             for s in STATES}
    maxsym=0.
    for rows in zip(*(reports[s]["original_geometries_both_legs"] for s in STATES)):
        fd,plus,minus=rows
        for leg in range(2):
            for ph in range(3):
                def J(row):
                    return row["per_leg_three_dimensionless_phases"][leg][ph][
                        "unnormalized_external_potential_impulse_shape_J"]
                u=.5*(J(plus)+J(minus))
                maxsym=max(maxsym,abs(u-J(fd))/max(1.,abs(J(fd)),abs(u)))
    require(maxsym<p["QA_before_new_numerics"][
        "original_pair_average_unnormalized_source_scaled_tolerance"],
        "original pair linear force numerator fails frozen F average")
    result={"date":"2026-09-27",
            "status":"E17D0_THREE_STATE_UNIT_EXTERNAL_POTENTIAL_VLASOV_SOURCE_SUSCEPTIBILITY_COMPLETE_PHYSICAL_FINITEK_B_BLOCKED",
            "prospective_protocol_git_blob":PRE_BLOB,
            "original_parents_sha256":p["immutable_parents"],
            "original_E16_geometries":576,"three_state_leg_phase_samples":10368,
            "original_state_reports":obj,
            "max_Fplus_minus_unnormalized_gravitational_source_pair_average_scaled_gap":maxsym,
            "physically_specified_halo_force_potential_history":False,
            "evolving_FRW_self_consistent_Vlasov_solution":False,
            "galactic_short_long_tracer_coupling_calibrated":False,
            "full_physical_finiteK_bispectrum":"BLOCKED",
            "observed_odd_SEALED":True,"new_CLASS_FITS_mock":False,
            "main_untouched_PR_draft":True}
    save_once(out/"e17d0_original_joint_external_potential_source_only.json",result)
    print("E17D0_THREE_STATE_EXTERNAL_SOURCE_ONLY_COMPLETE",maxsym,
          "FULL_PHYSICAL_B_BLOCKED",flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--state",choices=STATES)
    g.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    x=a.parse_args()
    if x.self_test:
        p,q,fs,sym=gate()
        require(len(geometry())==576,"original E16 triangles")
        print("E17D0_FROZEN_ORIGINAL_AND_PROTOCOL_SHA_PREFLIGHT_OK",
              len(q),sym,flush=True)
    elif x.state:run_state(x.state,x.output_dir)
    else:aggregate(x.output_dir)

if __name__=="__main__":
    main()
