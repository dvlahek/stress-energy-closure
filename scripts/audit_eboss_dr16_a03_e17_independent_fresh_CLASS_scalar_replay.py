#!/usr/bin/env python3
"""Independent E17A re-execution: scalar E16 geometry and fresh pinned CLASS all-leg replay.
Does not import E17 original implementation or consume its intermediate arrays.
Uses original E8 source bytes, the shared frozen cosmology backend, and direct CLASS calls.
No full physical three-point inference.
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
PRE=ROOT/"source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json"
PROTO_BLOB="f4bd3ea2a0519c4091ea94d0600ac62e57b3a2d4"
CSV=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E13D=ROOT/"source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27"
E14D=ROOT/"source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E16I=ROOT/"source_data/eboss_dr16_a03_e16_independent_original_scalar_geometry_replay_2026_09_27.json"
STATES=("FD","plus","minus")
COLS={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized","minus":"Fminus_CLASS_normalized"}
KS=(.05,.075,.1);KLS=(.001,.002,.003,.005)
MUS=(-1.,0.,.6,1.);MUL=(-1.,0.,.5,1.);PHIS=(0.,math.pi/2.,math.pi)
ZS=(.945,.95,.955)
C=299792.458

def check(a,msg):
    if not a:raise ValueError("E17_INDEPENDENT_FAIL_CLOSED: "+msg)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def relative(a,b):return abs(a-b)/max(abs(a),abs(b),1e-13)
def manual_transfer(ks,values,x):
    j=int(np.searchsorted(ks,x,side="right"))
    check(0<j<len(ks),"independent transfer extrapolation rejected")
    lo=float(ks[j-1]);hi=float(ks[j])
    return float(values[j-1]+(x-lo)/(hi-lo)*(values[j]-values[j-1])),(j-1,j)
def save_once(path,j):
    raw=(json.dumps(j,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"independent output changed; refuse overwrite")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17ind_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17_INDEPENDENT_CERTIFICATE",path,"SHA256",sha(raw),flush=True)

def run(original_dir,output):
    p=json.loads(PRE.read_bytes());pins=p["immutable_parents"]
    check(blob(PRE.read_bytes())==PROTO_BLOB and
          sha(CSV.read_bytes())==pins["E8_original_4000q_CSV_sha256"] and
          sha(E16.read_bytes())==pins["E16_original_full_sha256"] and
          sha(E16I.read_bytes())==pins["E16_independent_sha256"] and
          p["physics_STOP"]["observed_odd_SEALED"],"registered source SHA or observed STOP")
    jointraw=(original_dir/"e17_original_joint_three_state_CLASS_both_short_legs_source_only.json").read_bytes()
    joint=json.loads(jointraw)
    check(joint["prospective_E17_protocol_git_blob"]==PROTO_BLOB and
          joint["geometry_count"]==576 and joint["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED",
          "joint original/source STOP corrupted")
    import class_response_optimize as cro
    import wake_two_tracer_fisher as base
    from classy import Class
    check(cro.CLASS_COMMIT==pins["original_CLASS_git_commit"] and base.NS==.9649,
          "independent CLASS code revision changed")
    data=np.genfromtxt(CSV,delimiter=",",names=True)
    check(data.shape==(4000,) and np.all(np.diff(data["q_dimensionless"])>0),
          "independent original E8 CSV invalid")
    h=.6736
    qa={"replayed_state_count":0,"geometry_count_per_state":576,"long_short_pairs":12,
        "orientation_replays":0,"CLASS_Pcb_evaluations":0,"direct_vTk_theta_transfer_evaluations":0,
        "max_geometric_relative_gap":0.,"max_CLASS_Pcb_relative_gap":0.,
        "max_DIRECT_theta_relative_gap":0.,"max_cb_density_relative_gap":0.,
        "max_single_mode_velocity_relative_gap":0.,"max_static_Eq20_ratio_relative_gap":0.,
        "max_static_Eq20_time_stencil_relative_gap":0.,"max_anchor_E14_Pcb_gap":0.,
        "max_epsilon1e_minus5_squeezed_relative_gap":0.,"max_triangle_closure_h_Mpc":0.,
        "all_three_original_redshift_nodes":list(ZS)}
    original_shas={}
    for state in STATES:
        raw=(original_dir/("e17_original_both_legs_"+state+".json")).read_bytes()
        check(sha(raw)==joint["state_reports"][state]["sha256"],"original state changed after aggregate")
        original_shas[state]=sha(raw)
        a=json.loads(raw)
        check(a["state"]==state and a["original_CLASS_commit"]==cro.CLASS_COMMIT
              and len(a["rows_original_576_geometry_three_z_two_legs"])==576
              and a["physical_eBOSS_B_xi_SNR_significance"] is None
              and a["observed_odd_read"] is False,"state result or observed stop invalid")
        ref=json.loads((E13D/("e13_direct_wind_"+state+".json")).read_bytes())
        q=np.asarray(data["q_dimensionless"],float)
        f=np.asarray(data[COLS[state]],float)
        check(np.isfinite(f).all() and (f>0).all(),"independent frozen F positivity")
        outdir=output.parent;outdir.mkdir(parents=True,exist_ok=True)
        psd=outdir/("independent_E13_frozen_"+state+".dat")
        psdtext="\n".join("%.14e %.14e"%(x,y) for x,y in zip(q,f))+"\n"
        if psd.exists():check(psd.read_text()==psdtext,"independent PSD collision")
        else:psd.write_text(psdtext)
        check(sha(psd.read_bytes())==ref["reconstructed_PSD_from_original_archived_E8_SHA256"],
              "independent original E13 PSD SHA failed")
        base.ZBINS=np.asarray([.9,1.],float)
        params=base.class_params(psd,.06);params["output"]="mPk,dTk,vTk"
        check(params["gauge"]=="newtonian" and params["use_ncdm_psd_files"]==1 and
              params["P_k_max_h/Mpc"]>=.105,"independent CLASS cosmology drift")
        cos=Class();cos.set(params);cos.compute()
        try:
            check(abs(float(cos.h())-h)<1e-12,"independent h changed")
            tks={z:cos.get_transfer(z=z,output_format="class") for z in ZS}
            kin=np.asarray(tks[.95]["k (h/Mpc)"],float)
            check(np.all(np.diff(kin)>0) and kin[0]<.0475 and kin[-1]>.1025,
                  "independent short transfer k support lost")
            for z in ZS:check(np.array_equal(kin,np.asarray(tks[z]["k (h/Mpc)"],float)),
                               "independent CLASS transfer z sampling changed")
            e14=json.loads((E14D/("e14_short_Pcb_"+state+".json")).read_bytes())
            check(sha((E14D/("e14_short_Pcb_"+state+".json")).read_bytes())==
                  pins["E14_original_per_state_sha256"][state],"independent E14 anchor source changed")
            for k in KS:
                exact=float(cos.pk_cb_lin(k*h,.95))
                stored=e14["CLASS_Pcb_short_per_state_Mpc3"][str(k)]
                qa["max_anchor_E14_Pcb_gap"]=max(qa["max_anchor_E14_Pcb_gap"],relative(exact,stored))
                check(relative(exact,stored)<1e-8,"independent CLASS E14 anchor inconsistent")
            ix=0
            for k in KS:
                for K in KLS:
                    for ms in MUS:
                        for ml in MUL:
                            for ph in PHIS:
                                row=a["rows_original_576_geometry_three_z_two_legs"][ix];ix+=1
                                c=ms*ml+math.sqrt(max(0.,1.-ms*ms))*math.sqrt(max(0.,1.-ml*ml))*math.cos(ph)
                                r=K/k
                                moduli=[k*math.sqrt(1.+r*r/4.-r*c),k*math.sqrt(1.+r*r/4.+r*c)]
                                uk=(math.sqrt(max(0.,1.-ms*ms)),0.,ms)
                                ul=(math.sqrt(max(0.,1.-ml*ml))*math.cos(ph),
                                    math.sqrt(max(0.,1.-ml*ml))*math.sin(ph),ml)
                                k1=tuple(k*uk[j]-.5*K*ul[j] for j in range(3))
                                k2=tuple(-k*uk[j]-.5*K*ul[j] for j in range(3))
                                closure=math.sqrt(sum((k1[j]+k2[j]+K*ul[j])**2 for j in range(3)))
                                qa["max_triangle_closure_h_Mpc"]=max(qa["max_triangle_closure_h_Mpc"],closure)
                                check((row["k"],row["K"],row["mu_s"],row["mu_L"],row["phi"])==
                                      (k,K,ms,ml,ph),"changed original 48-orientation ordering")
                                for i,(v,m) in enumerate(zip((k1,k2),moduli),1):
                                    saved=row["k"+str(i)+"_h_Mpc"]
                                    qa["max_geometric_relative_gap"]=max(qa["max_geometric_relative_gap"],relative(saved,m))
                                    check(relative(saved,m)<1e-12,"independent geometry short leg mismatch")
                                    check(all(relative(t,u)<1e-12 for t,u in zip(v,row["k"+str(i)+"_vec"])),
                                          "independent vector short leg mismatch")
                                    check(abs(m-math.sqrt(sum(t*t for t in v)))<1e-13,
                                          "scalar and cartesian geometry disagree")
                                    check(abs(row["mu"+str(i)]-v[2]/m)<1e-12,
                                          "short LOS projection mismatch")
                                    for z in ZS:
                                        t=tks[z]
                                        pp=float(cos.pk_cb_lin(m*h,z))
                                        le=row["z_nodes"][str(z)]["leg"+str(i)]
                                        qa["CLASS_Pcb_evaluations"]+=1
                                        qa["max_CLASS_Pcb_relative_gap"]=max(qa["max_CLASS_Pcb_relative_gap"],
                                                                            relative(pp,le["P_cb_CLASS_Mpc3"]))
                                        theta_n,br=manual_transfer(kin,np.asarray(t["t_ncdm[0]"],float),m)
                                        theta_c,_=manual_transfer(kin,np.asarray(t["t_cdm"],float),m)
                                        theta=theta_n-theta_c
                                        db,_=manual_transfer(kin,np.asarray(t["d_b"],float),m)
                                        dc,_=manual_transfer(kin,np.asarray(t["d_cdm"],float),m)
                                        dcb=(cro.OMEGA_B*db+cro.OMEGA_CDM*dc)/(cro.OMEGA_B+cro.OMEGA_CDM)
                                        vel=-C*theta/(m*h)*math.sqrt(cro.A_S*(m*h/.05)**(base.NS-1.))
                                        qa["direct_vTk_theta_transfer_evaluations"]+=1
                                        for key,got,refkey in (
                                            ("max_DIRECT_theta_relative_gap",theta,"theta_rel_DIRECT_vTk_per_R_Mpc_inv"),
                                            ("max_cb_density_relative_gap",dcb,"delta_cb_Newtonian_per_R"),
                                            ("max_single_mode_velocity_relative_gap",vel,"vTk_E12_convention_single_mode_kms_NOT_R16_sigma")):
                                            qa[key]=max(qa[key],relative(got,le[refkey]))
                                        check(list(br)==le["CLASS_transfer_bracket_indices"],
                                              "independent transfer support bracket mismatch")
                                        alpha=le["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                                        check(math.isfinite(alpha) and alpha>0,"Eq20 model not positive")
                                for z in ZS:
                                    aa=row["z_nodes"][str(z)]["leg1"]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                                    bb=row["z_nodes"][str(z)]["leg2"]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                                    qa["max_static_Eq20_ratio_relative_gap"]=max(
                                        qa["max_static_Eq20_ratio_relative_gap"],
                                        relative(aa*moduli[0]**2,bb*moduli[1]**2))
                                for i in (1,2):
                                    lg="leg"+str(i)
                                    am=row["z_nodes"]["0.945"][lg]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                                    ap=row["z_nodes"]["0.955"][lg]["Eq20_positive_R16_rank_static_phase_MODEL_ONLY"]
                                    d=-(1.+.95)*(ap-am)/(.955-.945)
                                    qa["max_static_Eq20_time_stencil_relative_gap"]=max(
                                        qa["max_static_Eq20_time_stencil_relative_gap"],
                                        relative(d,row["Eq20_dphase_dln_a_MODEL_ONLY"][lg]))
                                # A true scalar linear CLASS spectrum is even under all-vector reversal.
                                # This checks real-field geometry only, NOT a complex physical B.
                                check(closure<1e-13 and
                                      relative(math.sqrt(sum(x*x for x in k1)),
                                               math.sqrt(sum((-x)**2 for x in k1)))<1e-14,
                                      "Hermitian real-field geometry failed")
                                # Independent vanishing-K Pcb diagnostic, not a finite-K B limit.
                                for sgn in (-1.,1.):
                                    ks=k*math.sqrt(1.+(r*1e-5)**2/4.+sgn*r*1e-5*c)
                                    pt=float(cos.pk_cb_lin(ks*h,.95))
                                    pk=float(cos.pk_cb_lin(k*h,.95))
                                    qa["max_epsilon1e_minus5_squeezed_relative_gap"]=max(
                                        qa["max_epsilon1e_minus5_squeezed_relative_gap"],relative(pt,pk))
                                qa["orientation_replays"]+=1
            check(ix==576,"independent all 576 orientations incomplete")
        finally:cos.struct_cleanup();cos.empty()
        qa["replayed_state_count"]+=1
    for key in ("max_geometric_relative_gap","max_CLASS_Pcb_relative_gap",
                "max_DIRECT_theta_relative_gap","max_cb_density_relative_gap",
                "max_single_mode_velocity_relative_gap","max_static_Eq20_ratio_relative_gap",
                "max_static_Eq20_time_stencil_relative_gap"):
        check(qa[key]<1e-8,"independent replay failed predeclared relative QA: "+key)
    check(qa["replayed_state_count"]==3 and qa["orientation_replays"]==1728
          and qa["CLASS_Pcb_evaluations"]==10368 and qa["direct_vTk_theta_transfer_evaluations"]==10368
          and qa["max_triangle_closure_h_Mpc"]<1e-13 and
          qa["max_anchor_E14_Pcb_gap"]<1e-8 and
          qa["max_epsilon1e_minus5_squeezed_relative_gap"]<1e-4,
          "independent full matrix / original E14 anchor / squeezed limit failed")
    ans={"date":"2026-09-27",
         "status":"E17A_INDEPENDENT_FRESH_PINNED_CLASS_ALL_10368_LEGS_TRANSFER_AND_SCALAR_GEOMETRY_REPLAY_PASS",
         "method":"Separate Python runner, no import of original E17 implementation; independently generate every original E16 scalar and Cartesian leg; re-run original three frozen CLASS states; separate linear transfer bracketing; exact E13 source PSD SHA. Same CLASS physics backend is NOT independent Einstein-Vlasov halo evolution.",
         "prospective_E17_protocol_git_blob":PROTO_BLOB,
         "original_joint_SHA256":sha(jointraw),"original_state_SHA256":original_shas,
         "numerical_QA":qa,"original_E14_E15_72_24_unchanged":True,
         "reversed_tracer_order_formal_factor_only":"(b_LRG-b_ELG)=-(b_ELG-b_LRG), no physical triple model",
         "Hermitian_scope":"Verified real CLASS scalar spectra and geometry under vector inversion, NOT physical full bispectrum",
         "retarded_Vlasov_halo_tracer_bispectrum":"BLOCKED",
         "no_observed_eBOSS_odd_read":True,"no_new_catalogue_mock_download":True,
         "no_main_change":True,"PR_remains_draft":True}
    save_once(output,ans)
    print("E17_INDEPENDENT_FRESH_THREE_STATE_CLASS_REPLAY_PASS",json.dumps(qa,sort_keys=True),
          "NO_PHYSICAL_B",flush=True)

if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--original-dir",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    x=a.parse_args();run(x.original_dir,x.output)
