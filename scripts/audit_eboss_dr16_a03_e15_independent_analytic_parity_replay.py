#!/usr/bin/env python3
"""Retrospective independent E15 analytic replay from immutable archived E14/E15.

No import of original E15 executable, no numpy, no CLASS, no FITS,
no observed odd. Uses exact Legendre orthogonality and synthetic
W_even/odd polynomial integrals. This is an independent ALGEBRA check
of the restricted E14 i*muLong*S model, NOT independent physical dynamics.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
E14=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json"
E15=ROOT/"source_data/eboss_dr16_a03_e15_original_E14_long_LOS_parity_and_toy_selection_2026_09_27.json"
E14_SHA="0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8"
E15_SHA="6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983"
STATES=("FD","plus","minus")
KS=("0.05","0.075","0.1")
LONG=("0.001","0.002","0.003","0.005")
SHORT_ELL=("1","3")
EPS=0.1
OUT=ROOT/"eboss_workspace/a03_physics_source/e15_independent_analytic_replay.json"

def sha(raw):return hashlib.sha256(raw).hexdigest()
def need(test,msg):
    if not test:raise ValueError(msg)
def scaled_gap(got,want,scale):
    return abs(float(got)-float(want))/max(1.,abs(float(scale)))

def main():
    e14raw=E14.read_bytes()
    e15raw=E15.read_bytes()
    need(sha(e14raw)==E14_SHA and sha(e15raw)==E15_SHA,
         "Original archived E14/E15 SHA gate failed")
    e14=json.loads(e14raw)
    e15=json.loads(e15raw)
    need(e14["original_E8_CSV_sha256"]==
         "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
         and e15["original_E14_full_SHA256"]==E14_SHA
         and e15["E15_prospective_protocol_git_blob"]==
         "71a8cb29314036f5899af5fd265d3a8779e46be5"
         and not e15["technical_QA_warnings"]
         and e15["physical_limits"]["original_observed_odd_sealed"],
         "E15 physics or original pre-registration provenance changed")
    ncases=0
    ncontrast=0
    max_scaled=0.
    worst=None
    for state in STATES:
        for k in KS:
            for K in LONG:
                for ell in SHORT_ELL:
                    source=e14["all_original_state_k_short_and_long_results"][state][k][
                        "four_original_filtered_long_modes"][K]["ells"][ell][
                        "Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"]
                    obs=e15["all_72_original_state_short_k_long_K_short_ell_cases"][state][k][K][ell]
                    scale=source
                    checks=[("frozen E14 real S",obs["original_reduced_E14_S_Mpc6"],source),
                        ("K positive imaginary",obs["actual_complex_B_at_muLong_plus1_unitBias"]["imag"],source),
                        ("K negative imaginary",obs["actual_complex_B_at_muLong_minus1_unitBias"]["imag"],-source)]
                    for n in ("32","64"):
                        projected=obs["quadrature_imaginary_coefficient"][n]
                        for longell in range(5):
                            checks.append((f"Legendre n{n} L{longell}",
                             projected["imag_long_multipoles_over_unit_bias_Mpc6"][str(longell)],
                             source if longell==1 else 0.))
                        mon=projected["synthetic_weighted_imag_monopole_over_unit_bias_Mpc6"]
                        dip=projected["synthetic_weighted_imag_dipole_over_unit_bias_Mpc6"]
                        checks.extend([
                            (f"even toy monopole n{n}",mon["even_W_1_plus_0p1_P2"],0.),
                            (f"odd toy monopole n{n}",mon["odd_W_1_plus_0p1_P1"],EPS*source/3.),
                            (f"even toy dipole n{n}",dip["even_W_1_plus_0p1_P2"],source*(1.+2.*EPS/5.)),
                            (f"odd toy dipole n{n}",dip["odd_W_1_plus_0p1_P1"],source)])
                    for label,got,want in checks:
                        gap=scaled_gap(got,want,scale)
                        if gap>max_scaled:max_scaled=gap;worst={"state":state,"k":k,"K":K,"ell":ell,"check":label,"gap":gap}
                        need(gap<1e-12,"Independent E15 analytic failure: "+str((state,k,K,ell,label,gap)))
                    ncases+=1
    for k in KS:
        for K in LONG:
            for ell in SHORT_ELL:
                Splus=e14["all_original_state_k_short_and_long_results"]["plus"][k][
                  "four_original_filtered_long_modes"][K]["ells"][ell][
                  "Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"]
                Sminus=e14["all_original_state_k_short_and_long_results"]["minus"][k][
                  "four_original_filtered_long_modes"][K]["ells"][ell][
                  "Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"]
                expected=Splus-Sminus
                archived=e14["Fplus_minus_Fminus_reduced_source_difference"][k][K][ell][
                    "Fplus_minus_Fminus_reduced_Bsource_Mpc6"]
                e15row=e15["all_24_original_Fplus_minus_Fminus_parity_contrasts"][k][K][ell]
                checks=[("E14 plus minus",archived,expected),
                  ("E15 plus minus",e15row["source_reconstructed_from_state_values_Mpc6"],expected),
                  ("E15 long dipole",e15row["imag_long_dipole_Mpc6"],expected),
                  ("E15 long monopole",e15row["imag_long_monopole_unweighted_Mpc6"],0.),
                  ("E15 odd toy leakage",e15row["synthetic_odd_selection_imag_monopole_Mpc6"],EPS*expected/3.)]
                for label,got,want in checks:
                    gap=scaled_gap(got,want,max(abs(Splus),abs(Sminus)))
                    max_scaled=max(max_scaled,gap)
                    need(gap<1e-12,"Independent E15 contrast failure: "+str((k,K,ell,label,gap)))
                ncontrast+=1
    need(ncases==72 and ncontrast==24,"Independent original E15 case counts wrong")
    out={"date":"2026-09-27",
      "status":"RETROSPECTIVE_INDEPENDENT_E15_ORIGINAL_72_ANGULAR_PARITY_AND_24_CONTRAST_ANALYTIC_REPLAY_PASS",
      "original_E14_full_sha256":E14_SHA,
      "original_E15_full_sha256":E15_SHA,
      "independence":"Pure Python standard-library exact Legendre/selection polynomial replay from original archived E14 and E15 values. E15 script not imported. Same RESTRICTED E14 source model, NOT independent physics, numerical CLASS run, real galaxy triangle or real survey-window validation.",
      "analytical_identities":{
        "B_unit_DeltaBias":"i muLong S",
        "long_L1_imag":"S",
        "long_L0_L2_L3_L4_imag":"0",
        "even_synthetic_selection_W":"1+0.1 P2, monopole 0 and dipole 1.04 S",
        "odd_synthetic_selection_W":"1+0.1 P1, monopole S/30 and dipole S",
        "long_K_reverse":"B(-K)=conjugate(B(K)), imaginary sign reversed"},
      "original_state_cases":ncases,
      "original_Fplus_minus_Fminus_contrasts":ncontrast,
      "max_scaled_absolute_replay_gap":max_scaled,
      "worst_case":worst,
      "physical_eBOSS_bispectrum_or_physical_window_calculated":False,
      "original_observed_odd_data_sealed":True,
      "new_CLASS_FITS_mock_download_or_science_seed_cut":False,
      "main_untouched":True}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(out,indent=2,allow_nan=False)+"\n").encode()
    if OUT.exists():need(OUT.read_bytes()==raw,"Original independent E15 output altered")
    else:OUT.write_bytes(raw)
    print("E15_INDEPENDENT_ORIGINAL_72_PLUS_24_ANALYTIC_REPLAY_PASS",
          "MAX_SCALED_GAP",max_scaled,"OUTPUT_SHA256",sha(raw),
          "NO_REAL_EBOSS_NO_OBSERVED",flush=True)

if __name__=="__main__":main()
