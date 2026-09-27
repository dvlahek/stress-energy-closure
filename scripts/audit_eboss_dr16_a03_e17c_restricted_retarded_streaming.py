#!/usr/bin/env python3
"""E17C: frozen-source, restricted collisionless retarded streaming factor.

This is ONLY the homogeneous force-free propagator between external impulses.
The halo force, self-consistent Vlasov source, tracer coupling and physical
finite-K galaxy bispectrum are deliberately NOT supplied or inferred.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRE = ROOT / "source_data/eboss_dr16_a03_e17c_restricted_retarded_kinetic_propagator_prereg_2026-09-27.json"
PRE_BLOB = "33defef2020f2931756e3ef7a3959d905d539165"
CSV = ROOT / "source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16 = ROOT / "source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E16IND = ROOT / "source_data/eboss_dr16_a03_e16_independent_original_scalar_geometry_replay_2026_09_27.json"
E17A = ROOT / "source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
E17B = ROOT / "source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_original_joint_native_CLASS_transfer_resolution.json"
A_MAN = ROOT / "source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/archive_manifest.json"
B_MAN = ROOT / "source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/archive_manifest.json"
STATES = ("FD", "plus", "minus")
COLS = {"FD":"F0_CLASS_normalized", "plus":"Fplus_CLASS_normalized",
        "minus":"Fminus_CLASS_normalized"}
KS = (.05,.075,.1)
KLS = (.001,.002,.003,.005)
MUS = (-1.,0.,.6,1.)
MUL = (-1.,0.,.5,1.)
PHIS = (0.,math.pi/2.,math.pi)
PHASES = (0.,.25,1.)
OUT = ROOT / "eboss_workspace/a03_physics_source/e17c_restricted_retarded_kinetic"

def check(ok, message):
    if not ok:
        raise ValueError("E17C_FAIL_CLOSED: " + message)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def blob(b):
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()

def pack(obj):
    return (json.dumps(obj, indent=2, allow_nan=False) + "\n").encode()

def save_once(path, obj):
    raw = pack(obj)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        check(path.read_bytes() == raw, "existing result differs: " + str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".e17c_", delete=False) as fh:
            tmp = Path(fh.name)
            fh.write(raw)
            fh.flush()
            os.fsync(fh.fileno())
        try:
            os.link(tmp, path)
        finally:
            tmp.unlink(missing_ok=True)
    print("E17C_ORIGINAL_SOURCE_ONLY_JSON", path, "SHA256", sha(raw), flush=True)

def gate():
    check(blob(PRE.read_bytes()) == PRE_BLOB, "prospective protocol changed")
    pre = json.loads(PRE.read_bytes())
    par = pre["immutable_parents"]
    for path, key in ((CSV,"E8_original_4000q_CSV_sha256"),
                      (E16,"E16_original_full_sha256"),
                      (E16IND,"E16_independent_sha256"),
                      (E17A,"E17A_original_joint_sha256"),
                      (E17B,"E17B_original_joint_sha256")):
        check(sha(path.read_bytes()) == par[key], "immutable original parent SHA: " + key)
    check(blob(A_MAN.read_bytes()) == par["E17A_manifest_git_blob"]
          and blob(B_MAN.read_bytes()) == par["E17B_manifest_git_blob"],
          "archived original manifest Git identity")
    e16 = json.loads(E16.read_bytes())
    e17a = json.loads(E17A.read_bytes())
    e17b = json.loads(E17B.read_bytes())
    check(e16["QA"]["original_72_cases"] == 72
          and e16["QA"]["original_24_contrasts"] == 24
          and e16["QA"]["orientation_probes_per_pair"] == 48
          and e17a["geometry_count"] == 576
          and e17a["eBOSS_observed_odd_read"] is False
          and e17b["physical_finite_K_bispectrum"] == "BLOCKED"
          and e17b["observed_odd_read"] is False,
          "original E16/E17A/E17B scope altered")
    fr = pre["frozen_geometry"]
    check(fr["states"] == list(STATES) and fr["tracer_order"] == ["LRG","ELG"]
          and fr["short_k_comoving_h_Mpc"] == list(KS)
          and fr["long_K_comoving_h_Mpc"] == list(KLS)
          and fr["mu_short"] == list(MUS) and fr["mu_long"] == list(MUL)
          and fr["geometries"] == 576
          and fr["z_CLASS_nodes"] == [.945,.95,.955]
          and pre["derived_restricted_kinetic_problem"][
              "prospectively_fixed_dimensionless_center_phases"] == list(PHASES)
          and all(pre["stop"].values()),
          "frozen scientific or STOP scope changed")
    arr = np.genfromtxt(CSV, delimiter=",", names=True)
    check(arr.shape == (4000,) and
          set(COLS.values()).issubset(set(arr.dtype.names or ())) and
          "q_dimensionless" in arr.dtype.names,
          "frozen 4000q CSV schema")
    q = np.asarray(arr["q_dimensionless"], dtype=float)
    f = {state: np.asarray(arr[col], dtype=float) for state,col in COLS.items()}
    check(np.all(np.isfinite(q)) and np.all(np.diff(q)>0)
          and all(np.all(np.isfinite(v)) and np.all(v>0) for v in f.values()),
          "q grid, positivity or finite kinetic inputs")
    sym = float(np.max(np.abs(.5*(f["plus"]+f["minus"])-f["FD"]))/
                np.max(np.abs(f["FD"])))
    check(sym < pre["QA_before_numerics"]["pointwise_Fplus_Fminus_average_FD_scaled_tolerance"],
          "frozen source pointwise pair averaging changed")
    return pre, q, f, sym

def geometry_rows():
    rows = []
    for k in KS:
        for K in KLS:
            for ms in MUS:
                for ml in MUL:
                    for idx,phi in enumerate(PHIS):
                        u = (math.sqrt(max(0.,1.-ms*ms)),0.,ms)
                        v = (math.sqrt(max(0.,1.-ml*ml))*math.cos(phi),
                             math.sqrt(max(0.,1.-ml*ml))*math.sin(phi),ml)
                        long = tuple(K*t for t in v)
                        leg1 = tuple(k*u[j]-.5*long[j] for j in range(3))
                        leg2 = tuple(-k*u[j]-.5*long[j] for j in range(3))
                        mod = [math.sqrt(sum(x*x for x in leg)) for leg in (leg1,leg2)]
                        closure = math.sqrt(sum((leg1[j]+leg2[j]+long[j])**2 for j in range(3)))
                        check(closure <= 1e-13,"E16 fixed triangle closure")
                        rows.append({"k":k,"K":K,"mu_s":ms,"mu_L":ml,
                                     "phi_index":idx,"phi_rad":phi,
                                     "leg_moduli_h_Mpc":mod,
                                     "closure_h_Mpc":closure})
    check(len(rows)==576,"original 576 E16 orientations")
    return rows

def characteristic(q, f, s, den):
    z = np.sinc((s*q)/math.pi)
    integral = float(np.trapz(q*q*f*z, q))
    return integral, integral/den

def state_calc(state, out):
    check(state in STATES,"unknown original E8 neutrino state")
    pre, q, fs, sym = gate()
    f = fs[state]
    den = float(np.trapz(q*q*f,q))
    check(math.isfinite(den) and den>0, "invalid original number-weighted normalization")
    records = []
    max_zero = 0.
    max_reversal = 0.
    bound = pre["QA_before_numerics"]["positive_weight_bound_C_abs"]
    for g in geometry_rows():
        measures = []
        for ki in g["leg_moduli_h_Mpc"]:
            leg = []
            for s0 in PHASES:
                s = s0*ki/.05
                ii, cc = characteristic(q, f, s, den)
                _, rev = characteristic(q, f, -s, den)
                check(math.isfinite(cc) and abs(cc) <= bound,
                      "positive-distribution Fourier characteristic bound failed")
                max_reversal = max(max_reversal, abs(cc-rev))
                if s0 == 0.:
                    max_zero = max(max_zero, abs(cc-1.))
                leg.append({"s_center":s0,"s_leg":s,
                            "I_unnormalized_truncated_q":ii, "C_normalized":cc})
            measures.append(leg)
        g["per_leg_three_dimensionless_phases"] = measures
        records.append(g)
    check(max_zero < pre["QA_before_numerics"]["all_C_zero_equals_one_abs"]
          and max_reversal < pre["QA_before_numerics"][
              "simultaneous_fourier_reversal_even_C_relative_gap"],
          "free-streaming zero-lag or real Fourier reversal")
    table = {(g["k"],g["K"],g["mu_s"],g["mu_L"],g["phi_index"]):g for g in records}
    max_swap=0.
    for g in records:
        other = table[(g["k"],g["K"],-g["mu_s"],g["mu_L"],2-g["phi_index"])]
        for j in range(2):
            max_swap=max(max_swap,abs(g["leg_moduli_h_Mpc"][j]-
                                       other["leg_moduli_h_Mpc"][1-j]))
            for i in range(len(PHASES)):
                max_swap=max(max_swap,abs(
                    g["per_leg_three_dimensionless_phases"][j][i]["C_normalized"]-
                    other["per_leg_three_dimensionless_phases"][1-j][i]["C_normalized"]))
    check(max_swap < pre["QA_before_numerics"][
          "leg_exchange_at_fixed_K_and_reversed_center_relative_gap"],
          "fixed E16 center reversal and both-leg exchange")
    result = {"date":"2026-09-27",
              "status":"E17C_RESTRICTED_ORIGINAL_FROZEN_SOURCE_FREE_STREAMING_PROPAGATOR_ONLY_PHYSICAL_HALO_B_BLOCKED",
              "state":state,"prospective_E17C_protocol_git_blob":PRE_BLOB,
              "original_E8_CSV_sha256":pre["immutable_parents"]["E8_original_4000q_CSV_sha256"],
              "source_frozen_pointwise_pair_average_scaled_gap":sym,
              "original_q_support_dimensionless":[float(q[0]),float(q[-1])],
              "original_q_samples":len(q),
              "truncated_q_number_weighted_integral":den,
              "n_original_E16_geometry":len(records),
              "n_leg_phase_values":len(records)*2*len(PHASES),
              "dimensionless_phase_grid_at_k0_not_a_cosmological_time":list(PHASES),
              "QA":{"zero_lag_max_abs":max_zero,"fourier_reversal_max_abs":max_reversal,
                    "leg_swap_max_abs":max_swap,
                    "max_triangle_closure_h_Mpc":max(g["closure_h_Mpc"] for g in records)},
              "original_geometries_both_legs_and_phase_characteristic":records,
              "interpretation":"Truncated original-q, isotropic fixed-epoch leading nonrelativistic field-free retarded transport between impulses ONLY; no force, halo formation, Einstein-Vlasov self-gravity, CLASS time propagation or physical observable.",
              "observed_odd_read":False,"new_catalogue_or_mock_download":False,
              "full_physical_finiteK_halo_tracer_B":"BLOCKED",
              "main_untouched_PR_draft":True}
    save_once(out/("e17c_restricted_streaming_"+state+".json"),result)
    print("E17C_STATE",state,"SOURCE_ONLY",len(records),"LEGSxPHASE",
          result["n_leg_phase_values"],"ZERO",max_zero,"SWAP",max_swap,
          "PHYSICAL_B_BLOCKED",flush=True)

def aggregate(out):
    pre, q, fs, sym = gate()
    reports = {}
    for state in STATES:
        path=out/("e17c_restricted_streaming_"+state+".json")
        raw=path.read_bytes()
        x=json.loads(raw)
        check(x["state"]==state and x["n_original_E16_geometry"]==576
              and x["n_leg_phase_values"]==3456 and x["observed_odd_read"] is False
              and x["full_physical_finiteK_halo_tracer_B"]=="BLOCKED",
              "original three-state result missing or science STOP drift")
        reports[state]={"sha256":sha(raw),"size_bytes":len(raw),"QA":x["QA"]}
    joint={"date":"2026-09-27",
           "status":"E17C_FROZEN_THREE_STATE_RESTRICTED_CAUSAL_STREAMING_COMPONENT_COMPLETE_FULL_PHYSICAL_B_BLOCKED",
           "prospective_protocol_git_blob":PRE_BLOB,
           "parents":pre["immutable_parents"],"state_reports":reports,
           "original_E16_geometries":576,"three_state_leg_phase_samples":10368,
           "actual_evolving_FRW_time_or_external_halo_potential_specified":False,
           "actual_finiteK_long_short_halo_tracer_coupling_calculated":False,
           "original_E14_E15_reduced_72_24_unchanged":True,
           "full_physical_bispectrum":"BLOCKED",
           "observed_odd_SEALED":True,"new_survey_download":False,
           "main_untouched_PR_draft":True}
    save_once(out/"e17c_original_joint_restricted_retarded_streaming.json",joint)
    print("E17C_RESTRICTED_CAUSAL_STREAMING_SOURCE_ONLY_COMPLETE",
          "FULL_PHYSICAL_B_BLOCKED",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--self-test",action="store_true")
    group.add_argument("--state",choices=STATES)
    group.add_argument("--aggregate",action="store_true")
    ap.add_argument("--output-dir",type=Path,default=OUT)
    args=ap.parse_args()
    if args.self_test:
        pre,q,fs,sym=gate()
        check(len(geometry_rows())==576,"frozen geometry")
        print("E17C_PROSPECTIVE_PARENT_AND_PHYSICS_STOP_GATE_OK",
              "ORIGINAL_Q_N",len(q),"PAIR_SYM",sym,flush=True)
    elif args.state:
        state_calc(args.state,args.output_dir)
    else:
        aggregate(args.output_dir)

if __name__=="__main__":
    main()
