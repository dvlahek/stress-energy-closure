#!/usr/bin/env python3
"""Independent, standard-library E16 scalar triangle and E14/E15 72/24 replay.

Does NOT import original E16 implementation, NumPy, CLASS, mocks or FITS.
Independently derives k1,k2, LOS projection, exact inverse-k2 envelopes, toy
non-unique finite-K witnesses and checks original 72/24 frozen source values.
This is independent ALGEBRA, not independent physical galaxy simulation.
"""
from __future__ import annotations
import hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
E14=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json"
E15=ROOT/"source_data/eboss_dr16_a03_e15_original_E14_long_LOS_parity_and_toy_selection_2026_09_27.json"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
P=ROOT/"source_data/eboss_dr16_a03_e16_exact_triangle_geometry_finite_squeeze_prereg_2026-09-27.json"
SHA={
 E14:"0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8",
 E15:"6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983",
 E16:"5d81f95a5571175f06eb032b865285ba34804701693f76af392be7fc31fb06db",
}
BLOBS={
 E14:"6ef0b2c3268a5bae9564b58abdd914328e861cba",
 E15:"9da057195b9fb34c113be53533af129f9de8809f",
 E16:"fcea9f3d73c39282d4833b993bfc5ca8e2bb7ff5",
 P:"28127fbedddc320660b6ac3f9d6262272905c7ca",
}
KS=(.05,.075,.1)
KL=(.001,.002,.003,.005)
MS=(-1.,0.,.6,1.)
ML=(-1.,0.,.5,1.)
PH=(0.,math.pi/2.,math.pi)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def digest(x):return hashlib.sha256(x).hexdigest()
def blob(x):return hashlib.sha1(b"blob "+str(len(x)).encode()+b"\0"+x).hexdigest()
def gap(a,b):return abs(a-b)/max(1.,abs(a),abs(b))
def comp(a,b,key):
    g=gap(a,b)
    need(g<1e-12,"E16 independent scalar audit mismatch: "+str(key))
    return g

def main():
    for path in BLOBS:
        raw=path.read_bytes()
        need(blob(raw)==BLOBS[path],"Original E16/E14/E15/protocol Git blob changed")
        if path in SHA:
            need(digest(raw)==SHA[path],"Original permanent full result SHA256 changed")
    a=json.loads(E14.read_bytes());b=json.loads(E15.read_bytes())
    q=json.loads(E16.read_bytes());p=json.loads(P.read_bytes())
    need(q["status"].startswith("E16_EXACT_THREE_VECTOR_CLOSED_TRIANGLE_")
         and q["original_E14_full_SHA256"]==SHA[E14]
         and q["original_E15_full_SHA256"]==SHA[E15]
         and q["QA"]["original_72_cases"]==72
         and q["QA"]["original_24_contrasts"]==24
         and q["physical_stop"]["observed_odd_SEALED"] is True
         and q["physical_stop"]["main_mutated"] is False
         and p["prohibitions"]["no_full_bispectrum_detection_significance_or_physical_acceptance"],
         "E16 physical STOP/original axes changed")
    worst=0.;checks=0;cases=[]
    outside=[]
    for k in KS:
        for K in KL:
            r=K/k
            s=q["twelve_original_short_long_geometry_cases"][str(k)][str(K)]
            ori=[]
            for ms in MS:
                for ml in ML:
                    for phi in PH:
                        c=(ms*ml+math.sqrt(max(0.,1.-ms*ms))*
                           math.sqrt(max(0.,1.-ml*ml))*math.cos(phi))
                        dplus=1.+r*r/4.-r*c
                        dminus=1.+r*r/4.+r*c
                        # This audit uses only scalar norm formulas, NO
                        # Cartesian vector construction from original E16.
                        q1=k*math.sqrt(dplus)
                        q2=k*math.sqrt(dminus)
                        m1=(ms-r*ml/2.)/math.sqrt(dplus)
                        m2=(-ms-r*ml/2.)/math.sqrt(dminus)
                        need(abs(m1)<=1.+1e-12 and abs(m2)<=1.+1e-12,
                             "Independent E16 scalar unit-LOS consistency failure")
                        f1=1./dplus;f2=1./dminus
                        g=.5*(f1+f2);odd=.5*(f1-f2)
                        ori.append((q1,q2,f1,f2,g,odd))
            need(len(ori)==48,"Independent original 48 orientations changed")
            allf=[f for row in ori for f in row[2:4]]
            limitlo=1./(1.+r/2.)**2
            limithi=1./(1.-r/2.)**2
            rawvals={
              "r_K_over_k":r,
              "max_abs_individual_static_inverse_k2_fractional_change":
                    max(limithi-1.,1.-limitlo),
              "uniform_relative_c_geometry_ONLY_mean_symmetrized_factor":
                    math.log((1.+r/2.)/(1.-r/2.))/r,
              "min_sampled_k1_h_per_Mpc":min(o[0] for o in ori),
              "max_sampled_k1_h_per_Mpc":max(o[0] for o in ori),
              "min_sampled_k2_h_per_Mpc":min(o[1] for o in ori),
              "max_sampled_k2_h_per_Mpc":max(o[1] for o in ori),
              "max_abs_symmetrized_TOY_completion_ratio_minus_original":
                    max(abs(o[4]-1.) for o in ori),
              "min_antisymmetric_TOY_leg_ratio":min(o[5] for o in ori),
              "max_antisymmetric_TOY_leg_ratio":max(o[5] for o in ori),
            }
            for key,val in rawvals.items():
                worst=max(worst,comp(float(s[key]),val,(k,K,key)));checks+=1
            for got,expected in zip(s[
                 "exact_continuous_c_envelope_individual_static_inverse_k2_factor"],
                  (limitlo,limithi)):
                worst=max(worst,comp(float(got),expected,("envelope",k,K)))
                checks+=1
            need(abs(min(allf)-limitlo)<1e-12 and
                 abs(max(allf)-limithi)<1e-12,
                 "Original 48 sample grid does not realize exact endpoint inverse-k2 bounds")
            flag=any(z<.05-1e-14 or z>.1+1e-14 for row in ori for z in row[:2])
            need(flag==s[
                 "at_least_one_short_leg_outside_original_three_k_CLASS_Pcb_support"],
                 "E16 original CLASS Pcb sampling envelope independent replay mismatch")
            if flag:outside.append([k,K])
            # In exact aligned-long LOS and transverse short k, both closed
            # short legs have modulus k sqrt(1+r²/4). Gsym=1/(1+r²/4)
            # and is a distinct finite-K completion of original S.
            wit=s["one_explicit_same_squeezed_limit_distinct_finiteK_witness"]
            worst=max(worst,comp(float(wit[
               "model_B_symmetrized_STATIC_GEOMETRY_ONLY_source_multiplier"]),
               1./(1.+r*r/4.),("nonunique witness",k,K)))
            cases.append({"k_short_h_Mpc":k,"K_long_h_Mpc":K,
                          "ratio":r,"fractional_single_leg_static_inverse_k2_envelope":
                            [limitlo,limithi],"sample_48":True,
                          "outside_original_three_CLASS_k_support":flag})
    need(len(outside)==8 and q["finite_ratio_scope"][
        "n_pairs_outside_original_3point_Pcb_support"]==8,
        "Independent source support coverage != 8 of 12")
    statecount=0;dcount=0
    for state in ("FD","plus","minus"):
        for k in KS:
            for K in KL:
                for ell in (1,3):
                    key=(state,str(k),str(K),str(ell))
                    S=float(a["all_original_state_k_short_and_long_results"][state][str(k)]
                        ["four_original_filtered_long_modes"][str(K)]["ells"][str(ell)]
                        ["Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"])
                    original=float(b["all_72_original_state_short_k_long_K_short_ell_cases"][
                        state][str(k)][str(K)][str(ell)]["original_reduced_E14_S_Mpc6"])
                    replay=float(q["all_72_original_state_short_long_shortell_source_cases"][
                        state][str(k)][str(K)][str(ell)][
                        "original_frozen_E14_E15_reduced_source_S_Mpc6"])
                    worst=max(worst,comp(S,original,key),comp(replay,original,key))
                    statecount+=1
    for k in KS:
        for K in KL:
            for ell in (1,3):
                key=(str(k),str(K),str(ell))
                S=lambda state:float(a["all_original_state_k_short_and_long_results"][
                     state][str(k)]["four_original_filtered_long_modes"][str(K)]["ells"][
                     str(ell)]["Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"])
                d=S("plus")-S("minus")
                e15=float(b["all_24_original_Fplus_minus_Fminus_parity_contrasts"][
                    str(k)][str(K)][str(ell)]["E14_Fplus_minus_Fminus_reduced_S_Mpc6"])
                e16=float(q["all_24_original_Fplus_minus_Fminus_reduced_source_contrasts"][
                    str(k)][str(K)][str(ell)][
                    "original_frozen_Fplus_minus_Fminus_E14_reduced_S_Mpc6"])
                worst=max(worst,comp(d,e15,key),comp(e16,e15,key));dcount+=1
    need(statecount==72 and dcount==24,
        "Original state/contrast case counts changed")
    result={"date":"2026-09-27",
        "status":"E16_INDEPENDENT_STANDARD_LIBRARY_ORIGINAL_E14_E15_E16_EXACT_SCALAR_GEOMETRY_AND_72_24_REPLAY_PASS",
        "original_full_SHA256":{"E14":SHA[E14],"E15":SHA[E15],"E16":SHA[E16]},
        "method":"Completely separate source-only scalar-triangle and SHA-gated original 72/24 replay, no E16 original code import, NumPy, CLASS, FITS or survey data.",
        "original_short_long_pairs":12,"independent_scalar_orientation_probes_per_pair":48,
        "original_state_source_cases":statecount,
        "original_plus_minus_contrasts":dcount,
        "independent_numeric_checks":checks,
        "max_scaled_replay_gap":worst,
        "pairs_with_unavailable_original_Pcb_short_leg_support":outside,
        "each_case":cases,
        "geometrical_fractional_inverse_k2_prefactor_is_NOT_actual_bispectrum_uncertainty":True,
        "finite_K_physical_three_point_source_calibrated":False,
        "observed_odd_sealed":True,"new_mock_catalogue_DOWNLOAD":False,
        "main_untouched":True}
    data=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    dest=ROOT/"eboss_workspace/a03_physics_source/e16_independent_original_scalar_geometry_replay.json"
    dest.parent.mkdir(parents=True,exist_ok=True)
    need(not dest.exists(),"Independent original E16 result exists; no overwrite")
    dest.write_bytes(data)
    print("EBOSS_A03_E16_INDEPENDENT_SCALAR_TRIANGLE_72_24_REPLAY_PASS",
          "ORIGINAL_E16_SHA256",SHA[E16],
          "REPORT_SHA256",digest(data),
          "MAX_REPLAY_GAP",worst,
          "N_OUTSIDE_SHORT_CLASS_PCB_SUPPORT",len(outside),
          "NO_ACTUAL_BISPECTRUM NO_OBSERVED_ODD",flush=True)
if __name__=="__main__":main()
