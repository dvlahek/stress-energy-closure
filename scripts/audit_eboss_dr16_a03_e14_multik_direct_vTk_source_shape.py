#!/usr/bin/env python3
"""A03E14: preregistered multi-kshort source-only E13 direct-rank shape.

Three separately evolved original F0/F+/- pinned CLASS P_cb(kshort,z=.95).
Only the restricted E13 quasistatic exact-Eq20 phase allows Gamma(kshort)
=Gamma(k0)*[k0/kshort]**2; all E13 direct-vTk Gaussian rank and R16
long-density-rank cross transfer are IDENTICAL original archived reports.
Unit DeltaBias fixed-mode source, not eBOSS galaxy bispectrum or 24D xi.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import os
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/"code"),str(ROOT/"scripts")]
PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e14_pre_registered_multik_direct_vTk_source_shape_2026-09-27.json"
P_BLOB="688c56b5486395911ff39ff6f862b405c1f446d7"
ARCH811=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27"
ARCH13=ROOT/"source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27"
INDEPENDENT=ROOT/"source_data/eboss_dr16_a03_e13_independent_original_source_replay_2026_09_27.json"
KINETIC=ROOT/"code/class_response_optimize.py"
BASE=ROOT/"code/wake_two_tracer_fisher.py"
E13_CODE=ROOT/"scripts/audit_eboss_dr16_a03_e13_rank_matched_direct_vTk_squeezed_source.py"
EXPECTED_BLOBS={
    KINETIC:"44956e353dbb238537d60da687d1a2a26fa13ee6",
    BASE:"473e6d9b5941d8f67d1ecc1f5f170124d9ef8ea2",
    E13_CODE:"b031c7d9bb95fd55eb4cffac19a46d0674325916",
    INDEPENDENT:"18b5cc834db70fb2619b619ea2026e84eefd5cb2",
    ARCH811/"archive_manifest.json":"83748b452f26b31da34a51c2c5de509d6b7ae3af",
    ARCH13/"archive_manifest.json":"67ea3407393ddd2b8fd679c493ef36550ad15642",
}
ORIGINAL_E8_SHA="bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
ORIGINAL_E13_SHA="537365b2e758ec215810bb687557771440851cf2869d02578569b1d47d4ec98d"
STATE_SHA={"FD":"f32527307024480e031e3d7952caf681b577f528d70b86fae1efa4ee7828bfde",
           "plus":"40d0da5e772ee99dce809a48179064a964ae20f46876dacde8d529158a8c7dd2",
           "minus":"bf94d3b701867492f91768b182839b149de68cce147de2d5b21a813c2d907a7a"}
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized","minus":"Fminus_CLASS_normalized"}
KSHORT=(.05,.075,.1)
LONGK=(.001,.002,.003,.005)
Z=.95
K0=.05
CLASS_SHA="e85808324f51fc694d12e3ed7439552a3c3f9540"
OUT=ROOT/"eboss_workspace/a03_physics_source/e14_multik_direct_rank_source"

def require(ok,msg):
    if not ok:raise ValueError(msg)

def sha(raw):return hashlib.sha256(raw).hexdigest()

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def stable(obj):
    return (json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()

def save_new(path,obj):
    raw=stable(obj)
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        require(path.read_bytes()==raw,"An E14 output already exists with DIFFERENT bytes")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e14_",delete=False) as fh:
            tmp=Path(fh.name);fh.write(raw);fh.flush();os.fsync(fh.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E14_SOURCE_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    require(blob(PROTOCOL.read_bytes())==P_BLOB,"E14 registered scope changed")
    for path,digest in EXPECTED_BLOBS.items():
        require(blob(path.read_bytes())==digest,
                "Original source/version/archive parent changed: "+str(path))
    p=json.loads(PROTOCOL.read_bytes())
    req=p["physics_and_grid"]
    require(req["states"]==list(STATES)
            and req["k_short_comoving_h_per_Mpc"]==list(KSHORT)
            and req["K_long_comoving_h_per_Mpc"]==list(LONGK)
            and req["z_math_midpoint"]==Z
            and req["exact_Eq20_FD_only"] is True
            and p["limits"]["observed_odd_sealed"] is True
            and p["limits"]["main_untouched"] is True,
            "E14 prereg physics, z/k, or observed-data STOP changed")
    manifest=json.loads((ARCH13/"archive_manifest.json").read_bytes())
    mfiles={entry["file"]:entry["sha256"] for entry in manifest["files"]}
    ep=ARCH13/"e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json"
    require(sha(ep.read_bytes())==ORIGINAL_E13_SHA
            and mfiles[ep.name]==ORIGINAL_E13_SHA,
            "Original fully archived E13 report changed")
    e13=json.loads(ep.read_bytes())
    require(e13["status"].startswith("E13_RANK_MATCHED_DIRECT_VTK_GAUSSIAN_WICK")
            and e13["original_E11_Gamma_IS_NOT_USED_IN_FINAL_Bsource"] is True
            and e13["observed_galaxy_random_or_odd_read"] is False
            and e13["original_short_k_COMOVING_h_Mpc"]==K0
            and e13["original_long_K_COMOVING_h_Mpc"]==list(LONGK),
            "Original E13 restricted source-only physics changed")
    source=ARCH811/"E8/frozen_matched_distributions_4000q.csv"
    require(sha(source.read_bytes())==ORIGINAL_E8_SHA,
            "Original E8 4000q source changed")
    arr=np.genfromtxt(source,delimiter=",",names=True)
    require(arr.shape==(4000,) and set(COLS.values()).issubset(arr.dtype.names or ())
            and np.all(np.diff(arr["q_dimensionless"])>0),
            "Original E8 source array schema changed")
    ind=json.loads(INDEPENDENT.read_bytes())
    require(ind["original_E13_SHA256"]==ORIGINAL_E13_SHA
            and ind["technical_replay_full_source_no_new_CLASS_or_FITS"] is True
            and ind["observed_odd_sealed"] is True,
            "Original E13 independent certificate changed")
    perstate={}
    for state in STATES:
        f=ARCH13/("e13_direct_wind_"+state+".json")
        require(mfiles.get(f.name)==STATE_SHA[state] and
                sha(f.read_bytes())==STATE_SHA[state],
                "Archived original direct CLASS state source SHA changed")
        obj=json.loads(f.read_bytes())
        require(obj["state"]==state and obj["CLASS_commit"]==CLASS_SHA
                and obj["observed_galaxy_random_or_odd_read"] is False,
                "Original E13 CLASS state/survey stop changed")
        perstate[state]=obj
    return p,arr,e13,perstate

def self_test():
    p,arr,e13,old=gate()
    for state in STATES:
        row=e13["conditional_Bsource_unit_DeltaBias_per_E8_state"][state]
        require(math.isclose(row["Pcb_short_Mpc3"],
            old[state]["CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95"],
            rel_tol=0.,abs_tol=2e-10),
            "Original E13 joint/individual Pcb k0 mismatch")
        for longk in LONGK:
            orig=row["long_K_fixed"][str(longk)]
            w=old[state]["fresh_direct_filtered_Plong_over_i_mu_fixed_long_K"][str(longk)]
            require(math.isclose(orig["P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"],
                                 w["P_delta_cb_rDIRECT_R16_filtered_over_i_mu_Mpc3"],
                                 rel_tol=1e-14),
                    "Original E13 filtered P_Dr state/longK mismatch")
            for ell in (1,3):
                ans=row["Gamma_direct_ell"+str(ell)]*row["Pcb_short_Mpc3"]*orig[
                    "P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"]
                got=orig["ell"+str(ell)+"_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"]
                require(math.isclose(ans,got,rel_tol=2e-14),
                        "Original E13 source product does not close at k0")
    print("E14_PRE_CLASS_E13_ALL_FROZEN_PCB_GAMMA_AND_R16_LONG_PRODUCT_CLOSURES_PASS",flush=True)
    return p,arr,e13,old

def run_state(state,out):
    require(state in STATES,"Only the frozen source triplet is allowed")
    p,arr,ref,old=self_test()
    import class_response_optimize as cro
    import wake_two_tracer_fisher as base
    from classy import Class
    require(cro.CLASS_COMMIT==CLASS_SHA and base.NS==.9649
            and math.isclose(cro.H0,67.36,rel_tol=0.,abs_tol=1e-12),
            "Pinned original CLASS settings changed")
    out.mkdir(parents=True,exist_ok=True)
    q=np.asarray(arr["q_dimensionless"],float)
    f=np.asarray(arr[COLS[state]],float)
    text="\n".join(f"{x:.14e} {y:.14e}" for x,y in zip(q,f))+"\n"
    psd=out/("e14_original_frozen_"+state+".dat")
    if psd.exists():require(psd.read_text()==text,"Original frozen PSD file already changed")
    else:psd.write_text(text)
    require(sha(psd.read_bytes())==old[state]["reconstructed_PSD_from_original_archived_E8_SHA256"],
            "E14 exact original E13-state PSD byte SHA mismatch")
    base.ZBINS=np.asarray([.9,1.],float)
    params=base.class_params(psd,.06)
    params["output"]="mPk,dTk,vTk"
    require(params["gauge"]=="newtonian" and params["z_max_pk"]>=1.05
            and params["P_k_max_h/Mpc"]>=max(KSHORT)
            and math.isclose(params["H0"],67.36,rel_tol=0.,abs_tol=1e-12),
            "E13 original CLASS configuration changed")
    cosmo=Class();cosmo.set(params);cosmo.compute()
    try:
        h=float(cosmo.h())
        require(math.isclose(h,.6736,rel_tol=0,abs_tol=1e-13),
                "Original E13 frozen source H0/h changed")
        powers={str(k):float(cosmo.pk_cb_lin(k*h,Z)) for k in KSHORT}
        require(all(math.isfinite(v) and v>0 for v in powers.values()),
                "Nonpositive or nonfinite frozen state CLASS Pcb(k)")
        anchor=old[state]["CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95"]
        gap=abs(powers[str(K0)]-anchor)/anchor
        require(gap<p["technical_QA_not_physical_budget"][
                "original_E13_k0_Pcb_and_full_product_relative_tolerance"],
                "Fresh CLASS frozen Pcb(k0) differs from archived E13 original beyond prereg QA")
        result={"status":"E14_ORIGINAL_FROZEN_STATE_CLASS_MULTIK_PCB_SOURCE_ONLY_COMPLETE",
                "state":state,"original_CLASS_commit":CLASS_SHA,
                "original_E8_CSV_SHA256":ORIGINAL_E8_SHA,
                "original_E13_full_JSON_SHA256":ORIGINAL_E13_SHA,
                "E14_preregistration_git_blob":P_BLOB,
                "source_exact_E13_PSD_SHA256":sha(psd.read_bytes()),
                "CLASS_gauge":"newtonian","CLASS_output":"mPk,dTk,vTk",
                "redshift_z_not_eBOSS_effective":Z,
                "h":h,"short_k_comoving_h_per_Mpc":list(KSHORT),
                "CLASS_Pcb_short_per_state_Mpc3":powers,
                "original_E13_k0_Pcb_relative_QA_gap":gap,
                "observed_galaxy_random_or_odd_read":False,
                "new_catalogues_or_mocks_or_science_cuts":False}
        save_new(out/("e14_short_Pcb_"+state+".json"),result)
        print("E14_CLASS_ORIGINAL_STATE_PCB",state,json.dumps(powers,sort_keys=True),
              "ANCHOR_REL_GAP",gap,flush=True)
    finally:cosmo.struct_cleanup();cosmo.empty()
    return 0

def aggregate(out):
    p,arr,e13,old=self_test()
    states={}
    for state in STATES:
        src=out/("e14_short_Pcb_"+state+".json")
        j=json.loads(src.read_bytes())
        require(j["status"]=="E14_ORIGINAL_FROZEN_STATE_CLASS_MULTIK_PCB_SOURCE_ONLY_COMPLETE"
                and j["state"]==state
                and j["original_E13_full_JSON_SHA256"]==ORIGINAL_E13_SHA
                and list(map(float,j["CLASS_Pcb_short_per_state_Mpc3"]))==list(KSHORT)
                and j["observed_galaxy_random_or_odd_read"] is False,
                "One of three prereg original E14 CLASS state reports is missing or invalid")
        states[state]=j
    rows={}
    maxanchor=0.
    for state in STATES:
        parent=e13["conditional_Bsource_unit_DeltaBias_per_E8_state"][state]
        source=states[state]
        krows={}
        for k in KSHORT:
            scale=(K0/k)**2
            pk=source["CLASS_Pcb_short_per_state_Mpc3"][str(k)]
            kval={}
            for ell in (1,3):
                kval[str(ell)]={"Gamma_direct_ell_conditional_at_k":parent["Gamma_direct_ell"+str(ell)]*scale,
                                 "original_E13_exact_k0_Gamma":parent["Gamma_direct_ell"+str(ell)]}
            longrows={}
            for K in LONGK:
                original=parent["long_K_fixed"][str(K)]
                pd=original["P_delta_cb_rdirect_R16_filtered_over_i_mu_Mpc3"]
                ellrows={}
                for ell in (1,3):
                    gamma=kval[str(ell)]["Gamma_direct_ell_conditional_at_k"]
                    computed=pk*gamma*pd
                    ellrows[str(ell)]={"Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6":computed,
                                      "Gamma_direct_rank":gamma}
                    if k==K0:
                        ref=original["ell"+str(ell)+"_Bsource_reduced_over_i_muLong_unit_DeltaBias_Mpc6"]
                        gap=abs(computed-ref)/max(abs(ref),1e-20)
                        maxanchor=max(maxanchor,gap)
                        require(gap<p["technical_QA_not_physical_budget"][
                            "original_E13_k0_Pcb_and_full_product_relative_tolerance"],
                            "New CLASS multik k0 joint product differs from archived E13 original beyond prereg QA")
                longrows[str(K)]={"original_E13_filtered_long_P_delta_cb_r_over_i_mu_Mpc3":pd,
                                  "K_over_k_short":K/k,
                                  "ells":ellrows}
                require(K/k<=p["physics_and_grid"]["max_K_long_over_k_short"]+1e-14,
                        "A predeclared triangle is outside restricted K/k<=.1 boundary")
            krows[str(k)]={"Pcb_short_CLASS_Mpc3":pk,"relative_quasistatic_Gamma_k_scaling_from_k0":scale,
                            "conditional_ell_Gamma":kval,"four_original_filtered_long_modes":longrows}
        rows[state]=krows
    delta={}
    for k in KSHORT:
        kk={}
        for K in LONGK:
            ll={}
            for ell in (1,3):
                path=lambda s:rows[s][str(k)]["four_original_filtered_long_modes"][str(K)]["ells"][str(ell)][
                    "Bsource_reduced_over_i_mu_long_unit_DeltaBias_Mpc6"]
                a=path("plus");d=path("minus")
                ll[str(ell)]={"Fplus_minus_Fminus_reduced_Bsource_Mpc6":a-d,
                    "separate_frozen_Pcb_and_direct_wind_and_long_cross_for_each_state":True}
                if k==K0:
                    ref=e13["conditional_Fplus_minus_Fminus_Bsource_unit_DeltaBias"][str(K)][
                        "delta_ell"+str(ell)+"_unit_DeltaBias_Bsource_over_i_muLong_Mpc6"]
                    gap=abs(a-d-ref)/max(abs(ref),1e-20)
                    maxanchor=max(maxanchor,gap)
                    require(gap<p["technical_QA_not_physical_budget"][
                        "original_E13_k0_Pcb_and_full_product_relative_tolerance"],
                        "Frozen original Fplus-minus contrast k0 E13 anchor mismatch")
            kk[str(K)]=ll
        delta[str(k)]=kk
    outj={"date":"2026-09-27","status":"E14_THREE_FROZEN_CLASS_SHORT_K_EXACT_STATIC_E13_SOURCE_SHAPE_PASS_NOT_EBOSS_BISPECTRUM",
       "preregistered_protocol_git_blob":P_BLOB,
       "original_E8_CSV_sha256":ORIGINAL_E8_SHA,
       "original_E13_full_JSON_sha256":ORIGINAL_E13_SHA,
       "original_E13_independent_source_certificate_git_blob":
           EXPECTED_BLOBS[INDEPENDENT],
       "CLASS_commit":CLASS_SHA,
       "original_three_E14_state_JSON_sha256":{
            s:sha((out/("e14_short_Pcb_"+s+".json")).read_bytes()) for s in STATES},
       "fixed_z_math_midpoint":Z,"short_k_comoving_h_Mpc":list(KSHORT),
       "four_original_long_K_comoving_h_Mpc":list(LONGK),
       "source_model":"Exact original E13 restricted Eq20 direct-vTk R16 Gaussian-rank Gamma scales as kshort^-2 at same z and original fixed source; original E13 state-specific R16-smoothed long cross P_Dr reused EXACTLY; new CLASS Pcb_short(k) separately for three original sources. Coefficient Bsource/(i mu_long DeltaBias), Mpc^6.",
       "all_original_state_k_short_and_long_results":rows,
       "Fplus_minus_Fminus_reduced_source_difference":delta,
       "original_E13_fixed_k0_max_relative_anchor_gap":maxanchor,
       "technical_QA_warnings":[],
       "scientific_scope":{"only_three_short_modes_K_over_k_le_0p1_not_continuous_triangles":True,
          "quasistatic_kminus2_assumption_not_full_Vlasov":True,
          "full_halo_bias_HOD_and_GR_selection_uncalibrated":True,
          "physical_eBOSS_triple_window_and_independent_covariance_absent":True,
          "existing_unconditional_24D_odd_unchanged":True,
          "observed_galaxies_randoms_odd_sealed":True,
          "no_new_catalogues_mocks_or_survey_cuts":True,
          "main_untouched":True,"PR_draft":True}}
    save_new(out/"e14_frozen_multik_direct_rank_unitbias_source_only.json",outj)
    print("E14_FROZEN_THREE_CLASS_MULTIK_DIRECT_RANK_FIXED_K_SOURCE_ONLY_PASS",
          "REPORT_SHA256",sha(stable(outj)),
          "MAX_ORIGINAL_E13_ANCHOR_REL_GAP",maxanchor,
          "NO_REAL_EBOSS_BISPECTRUM_OR_OBSERVED",flush=True)
    for k in KSHORT:
        print("E14_DELTA_FPLUS_MINUS_FMINUS_FIXED_SHORT_K",str(k),
              json.dumps(delta[str(k)],sort_keys=True),flush=True)
    return 0

def main():
    a=argparse.ArgumentParser(description=__doc__)
    x=a.add_mutually_exclusive_group(required=True)
    x.add_argument("--self-test",action="store_true")
    x.add_argument("--state",choices=STATES)
    x.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT)
    arg=a.parse_args()
    if arg.self_test:self_test();return 0
    if arg.state:return run_state(arg.state,arg.output_dir)
    return aggregate(arg.output_dir)

if __name__=="__main__":
    raise SystemExit(main())
