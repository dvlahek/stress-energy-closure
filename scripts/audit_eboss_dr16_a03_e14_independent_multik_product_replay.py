#!/usr/bin/env python3
"""Independent retrospective archived E14 multi-k original source replay.

No CLASS, E14 executable, FITS, original observed data or eBOSS covariance.
Recompute all 3 source x 3 kshort x 4 Klong x 2 ell reduced source
coefficients DIRECTLY from E13 original independently audited source and
E14 three original state-specific CLASS Pcb(k). Original E14 joint file is
used only as the held-out numerical comparison, not as a source input.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
E13=ROOT/"source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json"
E14=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27"
E14_MANIFEST=E14/"archive_manifest.json"
E13_SHA="537365b2e758ec215810bb687557771440851cf2869d02578569b1d47d4ec98d"
E14_SHA="0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8"
E14_MANIFEST_BLOB="ba747f413988693278a6322bffe2da69d4c52d68"
STATES=("FD","plus","minus")
KSHORT=(.05,.075,.1)
LONGK=(.001,.002,.003,.005)
ELLS=(1,3)
SOURCE_SHAS={"FD":"639eb2346f074bdcdf87aecaeee1f0e3f497ee065272b1373c8a804cb32582f4",
             "plus":"01dd19c93551ad846505980d673c454ad6dc36d9d0f03c7dda79f9c9388c9a30",
             "minus":"6d7c5b657aef2c9d7d05d7e028dd83b053649123192866312c0835aa81e67d32"}
OUT=ROOT/"eboss_workspace/a03_physics_source/e14_independent_multik_replay.json"

def need(x,msg):
    if not x:raise ValueError(msg)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def check_source():
    need(sha(E13.read_bytes())==E13_SHA,"Original E13 full report SHA changed")
    need(blob(E14_MANIFEST.read_bytes())==E14_MANIFEST_BLOB,
         "E14 archive manifest original blob changed")
    man=json.loads(E14_MANIFEST.read_bytes())
    record={x["file"]:x for x in man["files"]}
    joint=E14/"e14_frozen_multik_direct_rank_unitbias_source_only.json"
    need(sha(joint.read_bytes())==E14_SHA and record[joint.name]["sha256"]==E14_SHA,
         "Original E14 joint source report SHA changed")
    e14=json.loads(joint.read_bytes())
    e13=json.loads(E13.read_bytes())
    need(e13["status"].startswith("E13_RANK_MATCHED_DIRECT_VTK_GAUSSIAN")
         and e14["status"].startswith("E14_THREE_FROZEN_CLASS_SHORT_K")
         and e14["scientific_scope"]["observed_galaxies_randoms_odd_sealed"] is True,
         "Original E13/E14 scope changed or observed data unsealed")
    states={}
    for state in STATES:
        path=E14/("e14_short_Pcb_"+state+".json")
        need(record[path.name]["sha256"]==SOURCE_SHAS[state]
             and sha(path.read_bytes())==SOURCE_SHAS[state],
             "Original E14 CLASS state SHA changed")
        obj=json.loads(path.read_bytes())
        need(obj["state"]==state and obj["observed_galaxy_random_or_odd_read"] is False
             and obj["short_k_comoving_h_per_Mpc"]==list(KSHORT),
             "Original fixed original-F CLASS state/survey seal changed")
        states[state]=obj
    return e13,e14,states

def main():
    e13,e14,states=check_source()
    max_abs=0.
    max_rel=0.
    count=0
    delta_count=0
    result={}
    for state in STATES:
        e=e13["conditional_Bsource_unit_DeltaBias_per_E8_state"][state]
        s=states[state]["CLASS_Pcb_short_per_state_Mpc3"]
        a=e14["all_original_state_k_short_and_long_results"][state]
        records={}
        for k in KSHORT:
            kr={}
            pk=float(s[str(k)])
            scale=(.05/k)**2
            need(math.isclose(pk,a[str(k)]["Pcb_short_CLASS_Mpc3"],
                              rel_tol=0.,abs_tol=1e-10),
                 "E14 joint report does not reproduce original state Pcb")
            for K in LONGK:
                lr={}
                source_long=e["long_K_fixed"][str(K)][
                    "P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"]
                hold=a[str(k)]["four_original_filtered_long_modes"][str(K)]
                need(math.isclose(source_long,
                                  hold["original_E13_filtered_long_P_delta_cb_r_over_i_mu_Mpc3"],
                                  rel_tol=0.,abs_tol=1e-9),
                     "Original E13 R16 long cross changed in E14")
                for ell in ELLS:
                    gamma=e["Gamma_direct_ell"+str(ell)]*scale
                    expected=gamma*pk*source_long
                    reported=hold["ells"][str(ell)][
                        "Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"]
                    absdiff=abs(expected-reported)
                    reldiff=absdiff/max(abs(expected),abs(reported),1e-14)
                    max_abs=max(max_abs,absdiff);max_rel=max(max_rel,reldiff)
                    need(reldiff<4e-14,
                         "Independent original-source E14 product numerical mismatch")
                    lr[str(ell)]={"independent_Mpc6":expected,
                                  "original_E14_Mpc6":reported,
                                  "relative_gap":reldiff}
                    count+=1
                kr[str(K)]=lr
            records[str(k)]=kr
        result[state]=records
    contrasts={}
    for k in KSHORT:
        for K in LONGK:
            for ell in ELLS:
                a=result["plus"][str(k)][str(K)][str(ell)]["independent_Mpc6"]
                b=result["minus"][str(k)][str(K)][str(ell)]["independent_Mpc6"]
                reported=e14["Fplus_minus_Fminus_reduced_source_difference"][str(k)][str(K)][str(ell)][
                    "Fplus_minus_Fminus_reduced_Bsource_Mpc6"]
                gap=abs((a-b)-reported)/max(abs(a-b),abs(reported),1e-14)
                need(gap<4e-14,"E14 Fplus-Fminus state-only contrast replay failed")
                contrasts[f"k_{k}_K_{K}_ell_{ell}"]={"independent":a-b,
                     "original_E14":reported,"relative_gap":gap}
                delta_count+=1
    need(count==72 and delta_count==24,"Not all pre-registered E14 combinations were replayed")
    out={"date":"2026-09-27",
         "status":"E14_INDEPENDENT_3x3x4x2_FULL_SOURCE_ONLY_REPLAY_PASS",
         "audit_type":"RETROSPECTIVE_INDEPENDENT_FROM_E14_EXECUTABLE_ORIGINAL_ARCHIVED_INPUTS",
         "original_E13_sha256":E13_SHA,"original_E14_sha256":E14_SHA,
         "original_E14_archive_manifest_git_blob":E14_MANIFEST_BLOB,
         "E14_per_state_original_SHA256":SOURCE_SHAS,
         "numerical_case_count":count,
         "Fplus_minus_Fminus_contrast_case_count":delta_count,
         "max_absolute_independent_source_replay_gap_Mpc6":max_abs,
         "max_relative_independent_source_replay_gap":max_rel,
         "source_formula":"Archived E13 independent direct-vTk Gaussian Gamma_ell(k0)*(.05/kshort)^2 times ORIGINAL E14 frozen state Pcb(kshort) times original E13 W_R16-filtered long P_delta_cb,r(K). No E14 executable imported.",
         "full_numeric_state_k_short_Klong_ell_replay":result,
         "full_Fplus_minus_Fminus_replay":contrasts,
         "physical_STOP":{"not_actual_eBOSS_galaxy_bispectrum":True,
             "original_24D_unconditional_odd_not_converted_to_3point":True,
             "real_tracer_bias_HOD_evolution_selection_and_GR_uncalibrated":True,
             "three_short_modes_not_full_squeezed_triangle_observable":True,
             "survey_triple_window_and_independent_eBOSS_covariance_absent":True,
             "no_CLASS_or_FITS_or_mock_download":True,
             "observed_odd_sealed":True,
             "main_untouched":True}}
    OUT.parent.mkdir(parents=True,exist_ok=True)
    raw=(json.dumps(out,indent=2,allow_nan=False)+"\n").encode()
    if OUT.exists():need(OUT.read_bytes()==raw,"Existing independent E14 audit differs")
    else:OUT.write_bytes(raw)
    print("EBOSS_A03_E14_INDEPENDENT_ALL_72_PRODUCTS_AND_24_CONTRASTS_SOURCE_ONLY_REPLAY_PASS",
          "MAX_RELATIVE_GAP",max_rel,
          "AUDIT_SHA256",sha(raw),
          "NO_REAL_EBOSS_BISPECTRUM_NO_OBSERVED",flush=True)

if __name__=="__main__":
    main()
