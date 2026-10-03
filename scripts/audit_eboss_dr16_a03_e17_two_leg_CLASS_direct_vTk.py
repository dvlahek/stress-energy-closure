#!/usr/bin/env python3
"""E17A original frozen E16 576 triangles: both-leg CLASS linear Pcb and direct-vTk transfers.
This is NOT a finite-K halo/tracer bispectrum. All original E8-E16 physical sources are immutable.
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
sys.path[:0]=[str(ROOT/"code"),str(ROOT/"scripts")]
import audit_eboss_dr16_a03_e14_multik_direct_vTk_source_shape as e14
import audit_eboss_dr16_a03_e16_exact_closed_triangle_geometry as e16

PRE=ROOT/"source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json"
PRE_BLOB="f4bd3ea2a0519c4091ea94d0600ac62e57b3a2d4"
E16RAW=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E16IND=ROOT/"source_data/eboss_dr16_a03_e16_independent_original_scalar_geometry_replay_2026_09_27.json"
E14DIR=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27"
E13DIR=ROOT/"source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27"
STATES=("FD","plus","minus")
KS=(.05,.075,.1); KLS=(.001,.002,.003,.005)
MUS=(-1.,0.,.6,1.); MUL=(-1.,0.,.5,1.); PHIS=(0.,math.pi/2.,math.pi)
ZS=(.945,.95,.955)
CLASS_SHA="e85808324f51fc694d12e3ed7439552a3c3f9540"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17_two_leg_class"

def check(c,msg):
    if not c: raise ValueError("E17_FAIL_CLOSED: "+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def pack(o):return (json.dumps(o,indent=2,allow_nan=False)+"\n").encode()
def save_once(path,o):
    raw=pack(o);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"Existing output has different bytes; no overwrite: "+str(path))
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17_ORIGINAL_JSON",path,"SHA256",sha(raw),flush=True)

def gate():
    check(blob(PRE.read_bytes())==PRE_BLOB,"prospective E17 protocol blob changed")
    p=json.loads(PRE.read_bytes());parents=p["immutable_parents"]
    check(p["status"]=="E17_PROSPECTIVELY_REGISTERED_BEFORE_ANY_NEW_E17_CLASS_EXECUTION"
          and p["physics_STOP"]["observed_odd_SEALED"]
          and p["physics_STOP"]["no_main_change_or_PR_merge"],"E17 science/observed stop drift")
    check(parents["E14_original_full_sha256"]==e16.E14SHA
          and parents["E15_original_full_sha256"]==e16.E15SHA
          and parents["original_CLASS_git_commit"]==CLASS_SHA
          and parents["E14_runner_git_blob"]==blob(Path(e14.__file__).read_bytes())
          and parents["E16_protocol_git_blob"]==blob(e16.PRE.read_bytes()),"original runners or CLASS parents changed")
    check(sha(E16RAW.read_bytes())==parents["E16_original_full_sha256"]
          and sha(E16IND.read_bytes())==parents["E16_independent_sha256"],"original E16 full/independent SHA missing")
    r=json.loads(E16RAW.read_bytes());i=json.loads(E16IND.read_bytes())
    check(r["QA"]["original_72_cases"]==72 and r["QA"]["original_24_contrasts"]==24
          and r["QA"]["orientation_probes_per_pair"]==48
          and i["original_short_long_pairs"]==12
          and i["independent_scalar_orientation_probes_per_pair"]==48
          and i["observed_odd_sealed"] is True,"E16 original 72/24/576 archival identity changed")
    e16.source_gate()
    old=e14.self_test()
    for state in STATES:
        row=(E14DIR/("e14_short_Pcb_"+state+".json"))
        check(sha(row.read_bytes())==parents["E14_original_per_state_sha256"][state],
              "E14 per-state original three-mode anchor missing: "+state)
    f=p["fixed_physics_and_grid"]
    check(f["states"]==list(STATES) and f["original_tracer_order"]==["LRG","ELG"]
          and f["original_short_k_comoving_h_per_Mpc"]==list(KS)
          and f["original_long_K_comoving_h_per_Mpc"]==list(KLS)
          and f["z_math_nodes_not_eBOSS_z_eff"]==list(ZS)
          and f["orientations_per_kK"]==48
          and f["total_closed_geometries"]==576
          and f["CLASS_gauge"]=="newtonian","original E17 scope changed")
    return p,old

def geometry():
    """E16 frozen 4x4x3 orientation order and exact closed Cartesian legs."""
    rows=[]
    for k in KS:
        for K in KLS:
            for ms in MUS:
                for ml in MUL:
                    for ph in PHIS:
                        u=(math.sqrt(max(0.,1.-ms*ms)),0.,ms)
                        v=(math.sqrt(max(0.,1.-ml*ml))*math.cos(ph),
                           math.sqrt(max(0.,1.-ml*ml))*math.sin(ph),ml)
                        Kvec=tuple(K*x for x in v)
                        k1=tuple(k*u[j]-.5*Kvec[j] for j in range(3))
                        k2=tuple(-k*u[j]-.5*Kvec[j] for j in range(3))
                        m1=math.sqrt(sum(x*x for x in k1))
                        m2=math.sqrt(sum(x*x for x in k2))
                        close=math.sqrt(sum((k1[j]+k2[j]+Kvec[j])**2 for j in range(3)))
                        check(close<=1e-13,"closed triangle failure")
                        c=sum(u[j]*v[j] for j in range(3))
                        r=K/k
                        check(abs(m1/k-math.sqrt(1+r*r/4-r*c))<1e-13
                              and abs(m2/k-math.sqrt(1+r*r/4+r*c))<1e-13,
                              "independent exact leg moduli mismatch")
                        rows.append({"k":k,"K":K,"mu_s":ms,"mu_L":ml,"phi":ph,"c":c,
                            "K_vec":Kvec,"k1_vec":k1,"k2_vec":k2,"k1_h_Mpc":m1,
                            "k2_h_Mpc":m2,"mu1":k1[2]/m1,"mu2":k2[2]/m2,
                            "triangle_closure_h_Mpc":close})
    check(len(rows)==576,"E16 576 frozen orientations lost")
    return rows

def run_state(state,out):
    check(state in STATES,"unknown state");p,old=gate()
    import class_response_optimize as cro
    import wake_two_tracer_fisher as base
    import audit_eboss_dr16_a03_e12_direct_class_vTk_long_cross as e12
    import build_eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase as e9
    from classy import Class
    frozen_occupation,_,Eq20_physical=e9.source_only_self_test()
    check(cro.CLASS_COMMIT==CLASS_SHA and base.NS==.9649 and
          math.isclose(cro.H0,67.36,abs_tol=1e-12),"frozen CLASS parameters changed")
    arr=old[1];original_e13=old[3][state]
    q=np.asarray(arr["q_dimensionless"],float)
    f=np.asarray(arr[e14.COLS[state]],float)
    check(np.isfinite(f).all() and (f>0).all(),"frozen PSD positivity lost")
    out.mkdir(parents=True,exist_ok=True)
    psd=out/("e17_frozen_E13_"+state+".dat")
    txt="\n".join("%.14e %.14e"%(a,b) for a,b in zip(q,f))+"\n"
    if psd.exists():check(psd.read_text()==txt,"frozen source PSD changed")
    else:psd.write_text(txt)
    check(sha(psd.read_bytes())==original_e13["reconstructed_PSD_from_original_archived_E8_SHA256"],
          "E17 frozen E13 PSD exact SHA not reproduced")
    base.ZBINS=np.asarray([.9,1.],float)
    params=base.class_params(psd,.06);params["output"]="mPk,dTk,vTk"
    check(params["gauge"]=="newtonian" and params["P_k_max_h/Mpc"]>=.105
          and params["z_max_pk"]>=.955 and params["use_ncdm_psd_files"]==1,
          "original CLASS transfer/Pcb settings incompatible")
    cosmo=Class();cosmo.set(params);cosmo.compute()
    try:
        h=float(cosmo.h());check(abs(h-.6736)<1e-12,"frozen h changed")
        rows=geometry();legs=np.asarray([a[x] for a in rows for x in ("k1_h_Mpc","k2_h_Mpc")],float)
        check(legs.min()>=.0475-1e-12 and legs.max()<=.1025+1e-12,"new short k grid detected")
        transfer={};kin0=None
        for z in ZS:
            tk=cosmo.get_transfer(z=z,output_format="class")
            check({"k (h/Mpc)","t_ncdm[0]","t_cdm","d_b","d_cdm"}.issubset(tk),
                  "direct CLASS vTk/cb columns absent")
            kin=np.asarray(tk["k (h/Mpc)"],float)
            check(np.isfinite(kin).all() and np.all(np.diff(kin)>0)
                  and kin[0]<legs.min() and kin[-1]>legs.max(),"CLASS transfer support insufficient: no extrapolation")
            if kin0 is None:kin0=kin
            else:check(np.array_equal(kin0,kin),"CLASS transfer k grid changed across z")
            a=np.searchsorted(kin,legs,side="right")
            check((a>0).all() and (a<len(kin)).all(),"one or more short legs not bracketed by CLASS transfer")
            theta=np.interp(legs,kin,np.asarray(tk["t_ncdm[0]"],float)-np.asarray(tk["t_cdm"],float))
            dcb=np.interp(legs,kin,(cro.OMEGA_B*np.asarray(tk["d_b"],float)+
                    cro.OMEGA_CDM*np.asarray(tk["d_cdm"],float))/(cro.OMEGA_B+cro.OMEGA_CDM))
            pc=np.asarray([float(cosmo.pk_cb_lin(float(k)*h,z)) for k in legs])
            check(np.isfinite(pc).all() and (pc>0).all() and np.isfinite(theta).all()
                  and np.isfinite(dcb).all(),"CLASS numerical short-leg transfer or Pcb invalid")
            v=-e12.C*theta/(legs*h)*np.sqrt(e12.primordial_delta2(legs*h,cro.A_S,base.NS))
            transfer[str(z)]=(pc,theta,dcb,v,a)
        powers={str(k):float(cosmo.pk_cb_lin(k*h,.95)) for k in KS}
        anchors=json.loads((E14DIR/("e14_short_Pcb_"+state+".json")).read_bytes())["CLASS_Pcb_short_per_state_Mpc3"]
        anchor_gaps={str(k):abs(powers[str(k)]-anchors[str(k)])/anchors[str(k)] for k in KS}
        check(max(anchor_gaps.values())<p["numerical_QA_registered_before_CLASS"]["archived_E14_three_short_k_Pcb_relative_anchor_max"],
              "fresh CLASS Pcb fails original E14 3k anchors")
        origpk=original_e13["CLASS_pk_cb_short_Mpc3_at_kcom_0p05_z0p95"]
        check(abs(powers["0.05"]-origpk)/origpk<1e-8,"fresh CLASS E13 center anchor changed")
        maxsqueeze=0.
        for row in rows:
            k=row["k"];K=row["K"];c=row["c"]
            for sgn in (-1,1):
                kk=k*math.sqrt(1.+(K/k*1e-5)**2/4.+sgn*(K/k*1e-5)*c)
                testpk=float(cosmo.pk_cb_lin(kk*h,.95))
                pk0=float(cosmo.pk_cb_lin(k*h,.95))
                maxsqueeze=max(maxsqueeze,abs(testpk/pk0-1.))
        check(maxsqueeze<p["numerical_QA_registered_before_CLASS"]["continuous_squeezed_Pcb_limit_relative_diagnostic_at_epsilon_1e_minus_5"],
              "linear Pcb no continuous squeezed-limit diagnostic")
        for n,row in enumerate(rows):
            row["z_nodes"]={}
            for z in ZS:
                pc,t,dc,v,a=transfer[str(z)]
                row["z_nodes"][str(z)]={
                    "leg1":{"P_cb_CLASS_Mpc3":float(pc[2*n]),"theta_rel_DIRECT_vTk_per_R_Mpc_inv":float(t[2*n]),
                            "delta_cb_Newtonian_per_R":float(dc[2*n]),"vTk_E12_convention_single_mode_kms_NOT_R16_sigma":float(v[2*n]),
                            "CLASS_transfer_bracket_indices":[int(a[2*n]-1),int(a[2*n])]},
                    "leg2":{"P_cb_CLASS_Mpc3":float(pc[2*n+1]),"theta_rel_DIRECT_vTk_per_R_Mpc_inv":float(t[2*n+1]),
                            "delta_cb_Newtonian_per_R":float(dc[2*n+1]),"vTk_E12_convention_single_mode_kms_NOT_R16_sigma":float(v[2*n+1]),
                            "CLASS_transfer_bracket_indices":[int(a[2*n+1]-1),int(a[2*n+1])]}}
            # MODEL ONLY: carry the original direct-R16 rank into exact Eq20.
            # Short-leg CLASS theta is NOT substituted for this conditioned rank.
            for z in ZS:
                oldwind=original_e13["direct_R16_sigma_LOS_kms_three_fixed_z"][str(z)]
                for lab,kh in (("leg1",row["k1_h_Mpc"]),("leg2",row["k2_h_Mpc"])):
                    val=e9.alpha_given_v(frozen_occupation,state,z,oldwind,kh,Eq20_physical)
                    check(math.isfinite(val) and val>0,"Eq20 source occupancy/phase invalid")
                    row["z_nodes"][str(z)][lab]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]=val
            for lab in ("leg1","leg2"):
                am=row["z_nodes"]["0.945"][lab]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                ap=row["z_nodes"]["0.955"][lab]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                row.setdefault("Eq20_dphase_dln_a_MODEL_ONLY",{})[lab]=-(1.+.95)*(ap-am)/(.955-.945)
            for z in ZS:
                a1=row["z_nodes"][str(z)]["leg1"]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                a2=row["z_nodes"][str(z)]["leg2"]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                check(abs(a1*row["k1_h_Mpc"]**2-a2*row["k2_h_Mpc"]**2)/
                      max(a1*row["k1_h_Mpc"]**2,a2*row["k2_h_Mpc"]**2,1e-30)<1e-12,
                      "Eq20 k^-2 both-leg static ratio lost")
            row["central_original_Pcb_z095_Mpc3"]=powers[str(row["k"])]
            row["static_Eq20_individual_k_inverse_square_geometric_ratio_MODEL_ONLY"]=[
                (row["k"]/row["k1_h_Mpc"])**2,(row["k"]/row["k2_h_Mpc"])**2]
        ans={"date":"2026-09-27","status":"E17A_ORIGINAL_BOTH_LEG_CLASS_PCB_DIRECT_VTK_LINEAR_TRANSFER_COMPLETE_PHYSICAL_B_BLOCKED",
            "state":state,"preregistered_E17_git_blob":PRE_BLOB,
            "parent_original_E16_SHA256":p["immutable_parents"]["E16_original_full_sha256"],
            "parent_original_E14_state_SHA256":p["immutable_parents"]["E14_original_per_state_sha256"][state],
            "original_CLASS_commit":CLASS_SHA,"original_E13_PSD_SHA256":sha(psd.read_bytes()),
            "CLASS_gauge":"newtonian","CLASS_output":"mPk,dTk,vTk","h":h,
            "CLASS_transfer_grid":{"N":len(kin0),"min_k_h_Mpc":float(kin0[0]),"max_k_h_Mpc":float(kin0[-1]),
                                   "interpolation":"linear within directly verified CLASS grid, no extrapolation"},
            "old_E13_DIRECT_R16_rank_sigma_LOS_kms_by_z":original_e13["direct_R16_sigma_LOS_kms_three_fixed_z"],
            "original_E14_center_short_k_Pcb_anchor_replay_relative_gap":anchor_gaps,
            "max_Pcb_squeezed_epsilon1e_minus5_relative_gap":maxsqueeze,
            "rows_original_576_geometry_three_z_two_legs":rows,
            "CLASS_calculated":["P_cb(k1,k2,z) using direct CLASS pk_cb_lin","theta_ncdm-theta_cdm DIRECT CLASS vTk at each actual |ki| using within-grid transfer interpolation","Newtonian delta_cb at each actual |ki| and fixed three times"],
            "Eq20_MODEL_ONLY":["original exact frozen occupation at positive direct E13 R16 LOS rank on both actual legs at three fixed z; static phase and centered d/dln(a), NOT retarded dynamics or bispectrum"],
            "BLOCKED":["retarded nonlinear Vlasov halo response","finite-K halo/tracer three-point coupling","LRG/ELG bias/HOD/evolution/magnification/GR and triple survey window","independent eBOSS 3pt covariance"],
            "observed_odd_read":False,"new_catalogue_mock_download_science_seed_cut":False,"physical_eBOSS_B_xi_SNR_significance":None}
        save_once(out/("e17_original_both_legs_"+state+".json"),ans)
        print("E17_STATE_CLASS_PASS",state,"GRID_N",len(kin0),"ANCHOR_GAP_MAX",max(anchor_gaps.values()),
              "SQUEEZED_PCB_EPS1e-5_MAX",maxsqueeze,flush=True)
    finally:cosmo.struct_cleanup();cosmo.empty()

def aggregate(out):
    p,_=gate();reports={}
    for state in STATES:
        path=out/("e17_original_both_legs_"+state+".json")
        raw=path.read_bytes();a=json.loads(raw)
        check(a["state"]==state and len(a["rows_original_576_geometry_three_z_two_legs"])==576
              and a["observed_odd_read"] is False and a["physical_eBOSS_B_xi_SNR_significance"] is None,
              "one original per-state CLASS report invalid")
        reports[state]={"sha256":sha(raw),"bytes":len(raw),"max_E14_Pcb_anchor_relative_gap":
              max(a["original_E14_center_short_k_Pcb_anchor_replay_relative_gap"].values()),
              "max_squeezed_Pcb_diagnostic":a["max_Pcb_squeezed_epsilon1e_minus5_relative_gap"]}
    ans={"date":"2026-09-27","status":"E17A_ORIGINAL_THREE_STATE_576_EXACT_TRIANGLES_BOTH_CLASS_LEGS_LINEAR_TRANSFER_DONE_PHYSICAL_B_BLOCKED",
         "prospective_E17_protocol_git_blob":PRE_BLOB,
         "immutable_parent_SHA256":p["immutable_parents"],"state_reports":reports,
         "geometry_count":576,"state_geometry_count":1728,"redshift_nodes":list(ZS),
         "original_E14_E15_72_source_24_contrasts_UNMODIFIED":True,
         "both_short_leg_direct_CLASS_Pcb_and_vTk_linear_transfer_CALCULATED":True,
         "static_Eq20_distinct_from_CLASS_linear_transfer":True,
         "retarded_halo_tracer_physical_bispectrum":"BLOCKED",
         "tracer_swap_formal_antisymmetric_Delta_b_only_not_physical_amplitude":"Delta_b(LRG,ELG)=-Delta_b(ELG,LRG)",
         "Hermitian_reality_for_real_linear_CLASS_transfer_only":"Pcb(-ki)=Pcb(ki); i mu_long S has original E15 parity; full B cannot be asserted",
         "eBOSS_observed_odd_read":False,"eBOSS_24D_two_point_is_not_three_point":True,
         "new_catalogue_or_mock_download":False,"main_untouched":True,"PR_remains_draft":True}
    save_once(out/"e17_original_joint_three_state_CLASS_both_short_legs_source_only.json",ans)
    print("E17_ORIGINAL_JOINT_LINEAR_CLASS_DONE_NO_PHYSICAL_B",flush=True)

def main():
    a=argparse.ArgumentParser();g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true");g.add_argument("--state",choices=STATES);g.add_argument("--aggregate",action="store_true")
    a.add_argument("--output-dir",type=Path,default=OUT);o=a.parse_args()
    if o.self_test:
        gate();check(len(geometry())==576,"geometry preflight");print("E17_PRECLASS_SHA_576_GEOMETRY_GATE_PASS",flush=True)
    elif o.state:run_state(o.state,o.output_dir)
    else:aggregate(o.output_dir)
if __name__=="__main__":main()
