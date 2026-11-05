#!/usr/bin/env python3
"""E10: original-source, conditional LRGxELG Fourier-odd wind parity gate.

Uses only the SHA-pinned original frozen E8 CSV + E9 3-state CLASS report and
the original post-E9 preregistration. It computes a CONDITIONAL, per-unit
(b_LRG-b_ELG)*P_cb angular response at ONE original comoving k, not ξ_l(s).
For a symmetric +/- LOS wind sample independent of tracer selection the
intrinsic ensemble-mean odd multipoles cancel. Even-to-odd RR leakage is separate.
Do not use this program to claim absolute eBOSS ξ or survey significance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss, legval

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_eboss_dr16_a03_e8_frozen_resonant_occupation as e8

PROTOCOL = ROOT / "source_data/eboss_dr16_a03_e10_wind_parity_angular_odd_bridge_protocol_2026-09-27.json"
PROTOCOL_BLOB = "d9de48012103eaa07d531e11bca9cef285e62545"
E9_SUMMARY = ROOT / "source_data/eboss_dr16_a03_e9_CLASS_quasistatic_phase_independently_audited_summary_2026-09-27.json"
E9_SUMMARY_BLOB = "7596540697124e3c014369b6d1c051d6145d62dc"
ORIG_CSV_SHA = "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
ORIG_E9_SHA = "450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890"
STATES = ("FD", "plus", "minus")
ZCENTER = .95
HSTEP = .005
VELSTEP = .004
KCOM = .05


def require(ok, message):
    if not ok:
        raise ValueError(message)


def data_gate(csv_path, e9_path):
    e8.source_gate()
    for path, sha1 in ((PROTOCOL, PROTOCOL_BLOB), (E9_SUMMARY, E9_SUMMARY_BLOB)):
        require(e8.git_blob(path.read_bytes()) == sha1,
                "Immutable original/preregistered source Git blob changed: " + str(path))
    proto = json.loads(PROTOCOL.read_bytes())
    require(proto["fixed_angular_test"]["k_com_h_per_Mpc"] == KCOM
            and proto["fixed_angular_test"]["z_center"] == ZCENTER
            and proto["fixed_angular_test"]["phase_redshift_stencil_h"] == HSTEP
            and proto["fixed_angular_test"]["velocity_density_relative_z_step"] == VELSTEP
            and proto["fixed_angular_test"]["mu_quadrature_nodes"] == [256, 512]
            and proto["guards"]["observed_odd_sealed"] is True
            and proto["guards"]["not_abs_eBOSS_xi"] is True,
            "E10 preregistration/observational STOP was changed")
    summary = json.loads(E9_SUMMARY.read_bytes())
    require(summary["provenance"]["full_CI_report_sha256"] == ORIG_E9_SHA,
            "Original E9 archived summary changed its reported source SHA")
    require(e8.sha(csv_path.read_bytes()) == ORIG_CSV_SHA and
            e8.sha(e9_path.read_bytes()) == ORIG_E9_SHA,
            "Exact frozen E8 CSV or original E9 CLASS output SHA mismatch")
    e9 = json.loads(e9_path.read_bytes())
    require(e9["status"].startswith("E9_CLASS_WIND_QUASISTATIC_PHASE_DERIVATIVE_READY_")
            and e9["CLASS_commit_expected"] == e8.cro.CLASS_COMMIT
            and e9["source_only_E8_original_CSV_sha256"] == ORIG_CSV_SHA
            and e9["observed_odd_data_vector_read"] is False
            and e9["absolute_LRG_ELG_xi1_xi3_or_signal_to_noise_computed"] is False
            and len(e9["engineering_QA_warnings_not_survey_acceptance"]) == 0,
            "Parent original E9 physical assumptions or sealed data status changed")
    a = np.genfromtxt(csv_path, delimiter=",", names=True)
    require(a.dtype.names == ("q_dimensionless", "F0_CLASS_normalized",
                            "Fplus_CLASS_normalized", "Fminus_CLASS_normalized")
            and a.shape == (4000,), "Exact frozen E8 source grid changed")
    q = np.asarray(a["q_dimensionless"], float)
    f = {s: np.asarray(a[field], float) for s, field in
         (("FD", "F0_CLASS_normalized"),
          ("plus", "Fplus_CLASS_normalized"),
          ("minus", "Fminus_CLASS_normalized"))}
    require(np.isfinite(q).all() and np.all(np.diff(q) > 0)
            and all(np.isfinite(v).all() and np.all(v > 0) for v in f.values()),
            "Frozen E8 q-grid or physical PSD invalid")
    cases = {(float(row["wind_velocity_derivative_relative_z_step"]),
              float(row["z"])): row for row in
             e9["all_own_and_common_FD_wind_cases"]}
    require(len(cases) == 14 and
            all((VELSTEP, z) in cases for z in (ZCENTER-HSTEP, ZCENTER, ZCENTER+HSTEP)),
            "Original E9 2x7 original source cases are missing")
    return proto, q, f, cases


def fd(q):
    return 1. / (1. + np.exp(np.asarray(q)))


def occupation(qsource, fs, state, qarg):
    qarg = np.asarray(qarg, float)
    require(np.all(qarg >= qsource[0]) and np.all(qarg <= qsource[-1]),
            "Resonant momentum outside frozen original q grid")
    f = fd(qarg)
    if state == "FD":
        return f
    return f * np.interp(qarg, qsource, fs[state]) / np.interp(qarg, qsource, fs["FD"])


def angular_phase(q, fs, case, state, mu, wind_sign):
    """Published exact Eq20 occupancy with wind explicitly parallel to LOS."""
    require(state in STATES and wind_sign in (-1, 1), "Unknown frozen state or sign")
    mu = np.asarray(mu, float)
    z = float(case["z"])
    obj = case["state"][state]
    v = float(obj["own_CLASS_v_R16_1sigma_los_kms"])
    alpha1 = float(obj["alpha_with_own_CLASS_velocity"])
    require(v > 0 and alpha1 > 0, "Original positive conditional E9 wind/phase changed")
    q1 = e8.MASS_EV * v / (
        e8.C_KMS * e8.cro.T_NCDM * e8.cro.TCMB_K * e8.cro.KB_EV_K * (1. + z))
    occ1 = float(occupation(q, fs, state, np.asarray(q1)))
    angular_occ = occupation(q, fs, state, q1 * np.abs(mu))
    # alpha_aligned is original E9 alpha at mu=khat dot vhat=+1. Under
    # uniform +LOS vector v, angular k-dot-v and q resonance both vary.
    out = wind_sign * alpha1 * mu * angular_occ / occ1
    require(np.isfinite(out).all(), "Nonfinite exact-occupancy angular response")
    return out


def projected(proto, q, fs, cases, nmu, state, wind_sign):
    nodes, wt = leggauss(nmu)
    zm = cases[(VELSTEP, ZCENTER-HSTEP)]
    zp = cases[(VELSTEP, ZCENTER+HSTEP)]
    pminus = angular_phase(q, fs, zm, state, nodes, wind_sign)
    pplus = angular_phase(q, fs, zp, state, nodes, wind_sign)
    dalpha = -(1. + ZCENTER) * (pplus - pminus) / (2. * HSTEP)
    theta = nodes**2 * dalpha
    pells = {}
    for ell in (1, 3):
        pol = np.zeros(ell + 1); pol[-1] = 1.
        pells[str(ell)] = float((2*ell+1)/2. *
                                np.sum(wt * legval(nodes, pol) * theta))
    parity_error = float(np.max(np.abs(theta + theta[::-1])))
    require(parity_error < 2e-15, "LOS mu reflection lost odd parity")
    return pells, theta, nodes


def self_tests(proto, q, fs, cases):
    center = cases[(VELSTEP, ZCENTER)]
    mu = np.array([-1., -.7, 0., .25, .9, 1.])
    for s in STATES:
        aligned = float(center["state"][s]["alpha_with_own_CLASS_velocity"])
        one = float(angular_phase(q, fs, center, s, np.asarray(1.), 1))
        require(math.isclose(one, aligned, rel_tol=3e-14),
                "E10 angular extension fails original E9 aligned-mode closure")
        pos = angular_phase(q, fs, center, s, mu, 1)
        neg = angular_phase(q, fs, center, s, mu, -1)
        require(np.allclose(pos, -neg, atol=2e-18, rtol=1e-13)
                and pos[2] == 0,
                "Exact wind reversal/zero-angular-projection control failed")
        for v in (-1., -.5, 0., .5, 1.):
            x = float(angular_phase(q, fs, center, s, np.asarray(v), 1))
            y = float(angular_phase(q, fs, center, s, np.asarray(-v), 1))
            require(abs(x+y) < 5e-17, "Angular Fourier-mu oddness failed")
    # A8: Im C_LE / P = (bL-bE) mu^2 (dphi/dln a). Equal biases
    # and complex-conjugating tracer order are exact algebraic nulls.
    def imcross(x, bl, be, p=1.):
        return (bl-be)*p*x
    require(imcross(1.1, 1., 1.) == 0
            and imcross(1.1, 2., 1.) == -imcross(1.1, 1., 2.),
            "Two-tracer RSD equal-bias/null and conjugate orientation failed")
    # Exact constant-occupation angular control: theta(mu)=mu^3.
    muq, ww = leggauss(256)
    theta = muq**3
    p1 = 1.5 * float(np.sum(ww*muq*theta))
    p3 = 3.5 * float(np.sum(ww*legval(muq, [0,0,0,1])*theta))
    require(abs(p1-.6) < 1e-13 and abs(p3-.4)<1e-13,
            "Constant-occupation analytic mu^3=(3/5)P1+(2/5)P3 fails")
    print("E10_FROZEN_E9_AND_ORIGINAL_E8_SOURCE_SHA_ANGULAR_EQA8_NULL_CONTROLS_OK",
          "CONDITIONAL_ONLY", flush=True)


def compute(csv_path, e9_path):
    proto, q, fs, cases = data_gate(csv_path, e9_path)
    self_tests(proto, q, fs, cases)
    byn={}
    for nmu in (256, 512):
        pos={}; neg={}; combined={}
        for s in STATES:
            plus, theta_plus, nodes = projected(proto, q, fs, cases, nmu, s, 1)
            minus, theta_minus, nodes2 = projected(proto, q, fs, cases, nmu, s, -1)
            require(np.array_equal(nodes,nodes2) and
                    np.max(np.abs(theta_plus+theta_minus))<1e-15,
                    "Equal-weight signed wind reflection failed at every mu")
            pos[s]=plus;neg[s]=minus
            combined[s]={
               ell:float(.5*(plus[ell]+minus[ell]))
               for ell in ("1","3")}
            require(max(abs(combined[s]["1"]),abs(combined[s]["3"]))<2e-15,
                    "Spurious unconditional intrinsic odd survived +/- symmetric wind")
        byn[str(nmu)]={"conditional_positive_LOS":pos,
                        "conditional_negative_LOS":neg,
                        "symmetric_equal_weight_ensemble_intrinsic_odd_mean":combined,
                        "conditional_positive_plus_minus_difference":{
                           ell:pos["plus"][ell]-pos["minus"][ell] for ell in ("1","3")}}
    qa={}
    for state in STATES:
        for ell in ("1","3"):
            a=byn["256"]["conditional_positive_LOS"][state][ell]
            b=byn["512"]["conditional_positive_LOS"][state][ell]
            gap=abs(a-b)/max(abs(a),abs(b),1e-14)
            qa[state+"_ell"+ell]=gap
    for ell in ("1","3"):
        a=byn["256"]["conditional_positive_plus_minus_difference"][ell]
        b=byn["512"]["conditional_positive_plus_minus_difference"][ell]
        qa["delta_plus_minus_ell"+ell]=abs(a-b)/max(abs(a),abs(b),1e-14)
    require(max(qa.values()) < proto["fixed_angular_test"][
            "quad_256_vs_512_relative_QA_tolerance"],
            "Registered 256/512 angular quadrature QA failed (no retuning)")
    result={
        "date":"2026-09-27",
        "status":"E10_CONDITIONAL_ANGULAR_FOURIER_SOURCE_TESTS_PASS_UNCONDITIONAL_EBOSS_ODD_MEAN_NOT_ESTABLISHED",
        "protocol_git_blob_sha1":PROTOCOL_BLOB,
        "original_E8_4000q_CSV_sha256":ORIG_CSV_SHA,
        "original_E9_three_CLASS_JSON_sha256":ORIG_E9_SHA,
        "original_CLASS_commit":e8.cro.CLASS_COMMIT,
        "definition":"Original Eq A8 C_LE=<delta_LRG^s delta_ELG^s*>; unit-DeltaBias/unit-Pcb conditional T(mu)=mu² dphi(kcom,mu)/dln a, phi from exact FD Eq20 occupancy and original E9 +1sigma LOS wind, explicitly aligned with LOS",
        "frozen_example_not_measured_eBOSS":{"k_com_h_per_Mpc":KCOM,
                "z":ZCENTER,"mu_quad_nodes":[256,512],
                "z_stencil_step":HSTEP,
                "density_velocity_step_relative_1plusz":VELSTEP},
        "conditional_Fourier_source_projected_ells_not_xi":byn,
        "quadrature_256_512_relative_gaps":qa,
        "technical_QA_warnings":[],
        "physical_ensemble_statement":"Equal +/- LOS wind at same magnitude, independent of tracer selection => intrinsic mean conditional response zero for each ell. This is a necessary null for this simplified symmetric model, NOT a theorem for nonlinear correlated real halo/tracer selection.",
        "physical_ensemble_second_moment_qualifier":"Equal +/- fixed |v_LOS| two-point ensemble RMS equals absolute positive conditional projected coefficient; NOT the true velocity-distributed variance or survey S/N.",
        "tracer_response":"Im P_LRG_ELG(k,mu)=(b_LRG-b_ELG)P_cb(k,z)T(mu) under original Eq A8 assumptions; actual biases, P_cb covariance, halo velocity-density cross response and redshift-pair window NOT determined here.",
        "cannot_directly_propagate_E9_positive_LOS_into_unconditional_eBOSS_24D":True,
        "A03_physical_pair_z_RR_window_certified":False,
        "A04_independent_eBOSS_covariance_available":False,
        "observed_galaxies_or_randoms_or_odd_read":False,
        "new_catalogue_or_mock_download":False,
        "new_science_seeds_or_cuts":False,
        "full_retarded_Einstein_Vlasov_halo_response_computed":False,
        "absolute_eBOSS_xi1_xi3_or_detection_calculated":False
    }
    return result


def atomic(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == raw, "Original output exists with different bytes; fail closed")
        return
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e10_",delete=False) as f:
        tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
    try:os.link(tmp,path)
    finally:tmp.unlink(missing_ok=True)


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--e8-csv",type=Path,required=True)
    a.add_argument("--e9-json",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    z=a.parse_args()
    result=compute(z.e8_csv,z.e9_json)
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    atomic(z.output,raw)
    print("EBOSS_A03_E10_ORIGINAL_CLASS_SOURCE_ONLY_ANGULAR_PARITY_GATE_PASS",
          "RESULT_SHA256",hashlib.sha256(raw).hexdigest(),
          "MAX_QUAD_REL_GAP",max(result["quadrature_256_512_relative_gaps"].values()),
          "MEAN_SYMMETRIC_WIND_INTRINSIC_ODD_ZERO NO_SURVEY_XI NO_OBSERVED",
          flush=True)
    for state in STATES:
        print("E10_POSITIVE_CONDITIONAL_ELL1_ELL3",state,
              json.dumps(result["conditional_Fourier_source_projected_ells_not_xi"]["512"]
                          ["conditional_positive_LOS"][state],sort_keys=True),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
