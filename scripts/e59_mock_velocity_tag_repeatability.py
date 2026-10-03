#!/usr/bin/env python3
"""
E59: mock-only density-to-velocity-tag repeatability screen.

E58 closed the conditioned pair-estimator random-integration gate. E59 is the
first survey-derived velocity-tag step, but remains MOCK ONLY.

For each fixed EZmock ID and cap, E59:
  * re-verifies all 72 original compressed sources before FITS rows;
  * loads all frozen source-eligible LRG and ELG mock galaxies;
  * uses the already-frozen E58 R1 48k random subset as a fixed selection field;
  * constructs deterministic science-like probe midpoints from full mock
    LRGxELG DD pairs in 20--140 Mpc/h;
  * reconstructs a linear velocity-direction proxy separately from LRG and ELG
    density contrast using the exact R=16 Mpc/h real-space top-hat-smoothed
    1/r^2 kernel;
  * quantifies galaxy split-half repeatability and 192/256 Mpc/h kernel-cutoff
    stability;
  * uses mocks only to select (or reject) a tracer reconstruction for a later
    truth-labelled velocity calibration.

The same fixed random field is subtracted from both data halves. Therefore
split-half repeatability isolates galaxy sampling noise conditional on that
selection realization; it is NOT a random-mask convergence test and is NOT a
correlation with true baryon/halo velocity.

No observed galaxy rows, observed odd vector, covariance inverse, p-value,
detection, or physical absolute EV amplitude is used.
"""
from __future__ import annotations

import argparse, hashlib, json, math, os, sys, tempfile
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
import run_eboss_dr16_full_eligible_mock0001_sparse as V1
import e55_local_conditioned_full_eligible_ninemock_two48k as E55
from audit_eboss_dr16_rr_kdtree_pilot import cartesian

E58=ROOT/"source_data/e58_conditioned_full_eligible_ninemock_seven48k_result.json"
OUT=ROOT/"source_data/e59_mock_velocity_tag_repeatability.json"

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES

R_SMOOTH=16.0
RMAX_PRIMARY=256.0
RMAX_CHECK=192.0
PAIR_SMIN=20.0
PAIR_SMAX=140.0
NPROBE=512
N_RANDOM_MASK=48000

# E58 exact random-replica seed contract. E59 uses R1 only as a fixed selection
# realization during split-half galaxy repeatability; it does not reinterpret
# R1 as a converged pair-estimator random integral.
E58_SEED_ROOT=202609305400

# Prospectively linked to the E49 r_rec=0.70 screening point under the
# equal-independent-noise split-half diagnostic:
# r_half = r_rec^2/(2-r_rec^2) for r_rec=0.70.
RHALF_TARGET=(0.70**2)/(2.0-0.70**2)
SIGN_AGREE_TARGET=0.60
CUTOFF_CORR_TARGET=0.90
CUTOFF_SIGN_TARGET=0.90

SEED_ROOT=202609590000

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e59_",delete=False) as f:
        q=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(q,path)
    finally: q.unlink(missing_ok=True)

def arr_sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def cat_xyz(cat):
    return cartesian(cat,lambda z:V1.comoving_mpc_over_h(z,V1.PRIMARY_GEOMETRY))

def split_labels(n,seed):
    p=np.random.default_rng(seed).permutation(n)
    lab=np.empty(n,np.int8); lab[p[:n//2]]=0; lab[p[n//2:]]=1
    return lab

def local_seed(mid,cap,tracer,role,kind):
    return (SEED_ROOT + 1_000_000*IDS.index(mid) + 100_000*CAPS.index(cap)
            + 10_000*TRACERS.index(tracer) + 1_000*ROLES.index(role)
            + {"split":101,"probe":211}[kind])

def e58_r1_seed(mid,cap,tracer):
    return (E58_SEED_ROOT + 1_000_000*IDS.index(mid) + 10_000
            + 1_000*CAPS.index(cap) + 100*TRACERS.index(tracer))

def sample_fixed_random_mask(pool,mid,cap,tracer):
    need(len(pool[0])>=N_RANDOM_MASK,"Eligible random pool < 48k")
    pos=np.sort(np.random.default_rng(e58_r1_seed(mid,cap,tracer)).choice(
        len(pool[0]),size=N_RANDOM_MASK,replace=False))
    cat=tuple(np.asarray(v[pos],dtype="f8").copy() for v in pool)
    return cat,pos

def top_hat_velocity_kernel(delta,R):
    """Real-space kernel of a top-hat-smoothed linear velocity/gravity field.

    K(r)=r/R^3 for r<R and K(r)=r/r^3 outside. Overall positive factors cancel
    from the E59 correlations/signs.
    """
    r2=np.einsum("ij,ij->i",delta,delta)
    r=np.sqrt(r2)
    fac=np.empty_like(r)
    inside=r<R
    fac[inside]=1.0/(R**3)
    outside=~inside
    fac[outside]=1.0/np.maximum(r[outside]**3,1e-300)
    return delta*fac[:,None],r

def weighted_split_field(tree,xyz,w,labels,probes,rmax,chunk=64):
    out={key:np.zeros((len(probes),3),float)
         for key in ("full_primary","A_primary","B_primary",
                     "full_check","A_check","B_check")}
    for a0 in range(0,len(probes),chunk):
        a1=min(a0+chunk,len(probes))
        nbs=tree.query_ball_point(probes[a0:a1],rmax,workers=2)
        for local,js in enumerate(nbs):
            if not len(js): continue
            j=np.asarray(js,dtype=np.int64)
            delta=xyz[j]-probes[a0+local]
            K,r=top_hat_velocity_kernel(delta,R_SMOOTH)
            wk=K*w[j,None]
            out["full_primary"][a0+local]=wk.sum(axis=0)
            for lab,name in ((0,"A_primary"),(1,"B_primary")):
                m=labels[j]==lab
                if np.any(m): out[name][a0+local]=wk[m].sum(axis=0)
            mc=r<=RMAX_CHECK
            if np.any(mc):
                wkc=wk[mc]; jc=j[mc]
                out["full_check"][a0+local]=wkc.sum(axis=0)
                for lab,name in ((0,"A_check"),(1,"B_check")):
                    m=labels[jc]==lab
                    if np.any(m): out[name][a0+local]=wkc[m].sum(axis=0)
    return out

def weighted_full_field(tree,xyz,w,probes,rmax,chunk=64):
    out={key:np.zeros((len(probes),3),float) for key in ("primary","check")}
    for a0 in range(0,len(probes),chunk):
        a1=min(a0+chunk,len(probes))
        nbs=tree.query_ball_point(probes[a0:a1],rmax,workers=2)
        for local,js in enumerate(nbs):
            if not len(js): continue
            j=np.asarray(js,dtype=np.int64)
            delta=xyz[j]-probes[a0+local]
            K,r=top_hat_velocity_kernel(delta,R_SMOOTH)
            wk=K*w[j,None]
            out["primary"][a0+local]=wk.sum(axis=0)
            mc=r<=RMAX_CHECK
            if np.any(mc): out["check"][a0+local]=wk[mc].sum(axis=0)
    return out

def contrast_field(data,ranmask,probes,mid,cap,tracer,chunk):
    xd,ud,wd=cat_xyz(data); xr,ur,wr=cat_xyz(ranmask)
    ld=split_labels(len(xd),local_seed(mid,cap,tracer,"dat","split"))
    td=cKDTree(xd,leafsize=32); tr=cKDTree(xr,leafsize=32)
    D=weighted_split_field(td,xd,wd,ld,probes,RMAX_PRIMARY,chunk)
    R=weighted_full_field(tr,xr,wr,probes,RMAX_PRIMARY,chunk)

    sums={
      "D_full":float(np.sum(wd)),
      "D_A":float(np.sum(wd[ld==0])),"D_B":float(np.sum(wd[ld==1])),
      "R_full":float(np.sum(wr)),
    }
    need(min(sums.values())>0,"Nonpositive E59 catalogue weight sum")

    out={}
    for rad in ("primary","check"):
        rr=R[rad]/sums["R_full"]
        out["full_"+rad]=D["full_"+rad]/sums["D_full"]-rr
        out["A_"+rad]=D["A_"+rad]/sums["D_A"]-rr
        out["B_"+rad]=D["B_"+rad]/sums["D_B"]-rr
    return out,{
      "data_rows":len(xd),"random_mask_rows":len(xr),
      "data_split_rows":[int(np.sum(ld==0)),int(np.sum(ld==1))],
      "data_weight_sums":[sums["D_A"],sums["D_B"]],
      "random_mask_weight_sum":sums["R_full"],
      "random_mask_policy":"fixed E58 R1 48k; common to both data halves",
    }

def pair_probes(lrg,elg,mid,cap,nprobe):
    xl,ul,wl=cat_xyz(lrg); xe,ue,we=cat_xyz(elg)
    tree=cKDTree(xe,leafsize=32)
    rng=np.random.default_rng(local_seed(mid,cap,"eBOSS_LRG","dat","probe"))
    order=rng.permutation(len(xl))
    coslim=float(np.cos(np.deg2rad(.05)))
    probes=[]; seps=[]
    for i in order:
        js=tree.query_ball_point(xl[i],PAIR_SMAX)
        if not len(js): continue
        j=np.asarray(js,dtype=np.int64)
        delta=xe[j]-xl[i]
        s=np.linalg.norm(delta,axis=1)
        cost=ue[j]@ul[i]
        good=(s>=PAIR_SMIN)&(s<=PAIR_SMAX)&(cost<=coslim)&(cost>-1.)
        j=j[good]; s=s[good]
        if not len(j): continue
        q=int(rng.integers(len(j)))
        probes.append(.5*(xl[i]+xe[j[q]])); seps.append(float(s[q]))
        if len(probes)>=nprobe: break
    need(len(probes)==nprobe,
         f"E59 insufficient deterministic DD midpoint probes: {len(probes)}<{nprobe}")
    p=np.asarray(probes,float)
    return p,{
      "n":len(p),"probe_cartesian_SHA256":arr_sha(p),
      "separation_mean":float(np.mean(seps)),
      "separation_min":float(np.min(seps)),"separation_max":float(np.max(seps)),
    }

def pearson(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    a=a-np.mean(a); b=b-np.mean(b)
    den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(a@b/den) if den>0 else 0.0

def sign_stats(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    nz=(a!=0)&(b!=0)
    need(np.count_nonzero(nz)>=0.95*len(a),"Too many zero reconstructed LOS values")
    sf=float(np.mean(np.sign(a[nz])*np.sign(b[nz])))
    return {"sign_factor":sf,"sign_agreement":0.5*(1+sf),
            "n_nonzero":int(np.count_nonzero(nz))}

def los(v,probes):
    n=probes/np.linalg.norm(probes,axis=1)[:,None]
    return np.einsum("ij,ij->i",v,n)

def metrics(field,probes):
    A=los(field["A_primary"],probes)
    B=los(field["B_primary"],probes)
    F=los(field["full_primary"],probes)
    C=los(field["full_check"],probes)
    rhalf=pearson(A,B)
    s=sign_stats(A,B)
    rcut=pearson(F,C); scut=sign_stats(F,C)
    reliability=max(0.0,min(1.0,2*rhalf/(1+rhalf))) if rhalf>-0.999999 else 0.0
    implied=math.sqrt(reliability)
    return {
      "split_half_pearson":rhalf,
      "split_half_sign_factor":s["sign_factor"],
      "split_half_sign_agreement":s["sign_agreement"],
      "equal_independent_noise_implied_full_vs_latent_r_DIAGNOSTIC_ONLY":implied,
      "cutoff_192_vs_256_pearson":rcut,
      "cutoff_192_vs_256_sign_agreement":scut["sign_agreement"],
      "full_primary_los_rms":float(np.sqrt(np.mean(F*F))),
      "full_check_los_rms":float(np.sqrt(np.mean(C*C))),
    }

def describe(vals):
    a=np.asarray(vals,float)
    need(np.isfinite(a).all(),"Nonfinite E59 aggregate metric")
    return {"n":len(a),"mean":float(np.mean(a)),"median":float(np.median(a)),
            "min":float(np.min(a)),"max":float(np.max(a)),
            "sample_sd":float(np.std(a,ddof=1))}

def aggregate(state):
    out={}
    for cap in CAPS:
        out[cap]={}
        cases=[state["cases"][f"{mid:04d}/{cap}"] for mid in IDS]
        for tr in TRACERS:
            rr=[c["tracers"][tr]["metrics"]["split_half_pearson"] for c in cases]
            ss=[c["tracers"][tr]["metrics"]["split_half_sign_agreement"] for c in cases]
            cc=[c["tracers"][tr]["metrics"]["cutoff_192_vs_256_pearson"] for c in cases]
            cs=[c["tracers"][tr]["metrics"]["cutoff_192_vs_256_sign_agreement"] for c in cases]
            ir=[c["tracers"][tr]["metrics"]["equal_independent_noise_implied_full_vs_latent_r_DIAGNOSTIC_ONLY"]
                for c in cases]
            out[cap][tr]={
              "split_half_pearson":describe(rr),
              "split_half_sign_agreement":describe(ss),
              "cutoff_192_vs_256_pearson":describe(cc),
              "cutoff_192_vs_256_sign_agreement":describe(cs),
              "implied_full_vs_latent_r_DIAGNOSTIC_ONLY":describe(ir),
            }
        cross=[c["LRG_vs_ELG_full_primary"] for c in cases]
        out[cap]["LRG_vs_ELG_full_primary_pearson"]=describe(cross)

    candidates={}
    for tr in TRACERS:
        min_r=min(out[c][tr]["split_half_pearson"]["median"] for c in CAPS)
        min_sign=min(out[c][tr]["split_half_sign_agreement"]["median"] for c in CAPS)
        min_cut_r=min(out[c][tr]["cutoff_192_vs_256_pearson"]["median"] for c in CAPS)
        min_cut_sign=min(out[c][tr]["cutoff_192_vs_256_sign_agreement"]["median"] for c in CAPS)
        passed=(min_r>=RHALF_TARGET and min_sign>=SIGN_AGREE_TARGET and
                min_cut_r>=CUTOFF_CORR_TARGET and min_cut_sign>=CUTOFF_SIGN_TARGET)
        candidates[tr]={
          "minimum_cap_median_split_half_pearson":min_r,
          "minimum_cap_median_split_half_sign_agreement":min_sign,
          "minimum_cap_median_cutoff_pearson":min_cut_r,
          "minimum_cap_median_cutoff_sign_agreement":min_cut_sign,
          "repeatability_screen_pass":bool(passed),
        }
    passed=[t for t in TRACERS if candidates[t]["repeatability_screen_pass"]]
    selected=(max(passed,key=lambda t:candidates[t]["minimum_cap_median_split_half_pearson"])
              if passed else None)
    return out,candidates,selected

def self_test():
    d=np.array([[0.,0.,0.],[8.,0.,0.],[16.,0.,0.],[32.,0.,0.]])
    K,r=top_hat_velocity_kernel(d,R_SMOOTH)
    need(np.allclose(K[0],0),"E59 kernel r=0 fail")
    need(np.isclose(K[1,0],8/R_SMOOTH**3),"E59 inside top-hat kernel fail")
    need(np.isclose(K[2,0],16/16**3),"E59 boundary kernel fail")
    need(np.isclose(K[3,0],32/32**3),"E59 outside kernel fail")
    need(abs(RHALF_TARGET-0.3245033112582781)<1e-14,"E59 rhalf target drift")
    a=split_labels(100,12); b=split_labels(100,12)
    need(np.array_equal(a,b) and np.sum(a==0)==50,"E59 deterministic split fail")
    need(e58_r1_seed(1,"NGC","eBOSS_LRG")==202609315400,
         "E59 no longer reuses E58 R1 NGC LRG seed")
    print("E59_SYNTHETIC_VELOCITY_RECONSTRUCTION_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--probes",type=int,default=NPROBE)
    ap.add_argument("--chunk",type=int,default=64)
    args=ap.parse_args()
    self_test()
    if args.self_test: return
    need(args.run,"Use --run for E59 mock-only reconstruction")
    need(args.probes==NPROBE,"E59 probe count is frozen")
    need(16<=args.chunk<=256,"Unsafe E59 chunk")
    need(E58.is_file(),"Missing local E58 result")
    e58=json.loads(E58.read_text())
    need(e58["status"]=="PASS_NINEMOCK_FULL_ELIGIBLE_SEVEN48K_BACKGROUND_QUANTIFIED"
         and e58["decision"]["sevenrep_random_mean_subdominant_all_caps"] is True,
         "E58 random-integration gate is not closed")
    need(e58["observed_odd_used"] is False and e58["observed_galaxy_rows_used"] is False,
         "E58 observation guardrail changed")

    p,a02,e4=V1.source_only_gate()
    protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        protocol,require_local=True)
    paths,total=A02.resolve_72_sources(protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=a02
    need(archive["completed_cases"]==18 and archive["observed_odd_data_vector_read"] is False,
         "A02 parent changed")

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state["stage"]=="E59_MOCK_VELOCITY_TAG_REPEATABILITY",
             "Existing E59 checkpoint incompatible")
    else:
        state={
          "stage":"E59_MOCK_VELOCITY_TAG_REPEATABILITY",
          "date":"2026-09-30","status":"INCOMPLETE",
          "parent_E58_gate_closed":True,
          "mock_ids":list(IDS),"caps":list(CAPS),"tracers":list(TRACERS),
          "physics":{
            "smoothing":"exact real-space top-hat-smoothed linear velocity kernel",
            "R_Mpc_over_h":R_SMOOTH,
            "primary_kernel_cutoff_Mpc_over_h":RMAX_PRIMARY,
            "cutoff_stability_check_Mpc_over_h":RMAX_CHECK,
            "probe_midpoint_s_range_Mpc_over_h":[PAIR_SMIN,PAIR_SMAX],
            "probes_per_case":NPROBE,
          },
          "random_selection_policy":{
            "rows_per_tracer_cap":N_RANDOM_MASK,
            "subset":"exact frozen E58 R1 without-replacement seed",
            "same_selection_field_subtracted_from_both_data_halves":True,
            "purpose":"isolate galaxy split-half repeatability conditional on fixed numerical mask"
          },
          "screen":{
            "split_half_pearson_min":RHALF_TARGET,
            "split_half_sign_agreement_min":SIGN_AGREE_TARGET,
            "cutoff_pearson_min":CUTOFF_CORR_TARGET,
            "cutoff_sign_agreement_min":CUTOFF_SIGN_TARGET,
            "selection_rule":"among passing tracer reconstructions choose largest minimum-cap median split-half Pearson",
          },
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,"observed_odd_used":False,
          "true_velocity_available":False,
          "physical_velocity_truth_correlation_calibrated":False,
          "absolute_EV_amplitude_calibrated":False,
          "covariance_inverse_used":False,"pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        case=state["cases"].setdefault(ck,{"status":"INCOMPLETE","id":mid,"cap":cap})
        if case.get("status")=="complete":
            print("E59_REUSE_CASE",ck,flush=True); continue

        archived=archive["cases"][ck]
        prior=next((z for z in old0001["cases"] if z["cap"]==cap),None) if mid==1 else None
        full_data={}; random_pool={}; replay={}
        for tr in TRACERS:
          k=A02.key(mid,cap,tr)
          for role in ROLES:
            src=paths[k+"/"+role]
            rows=A02.source_header_rows(src,prior,tr,role)
            fixed,full,info=E55.load_full_verified(
                src,mid=mid,cap=cap,tracer=tr,role=role,expected_rows=rows,
                archived_info=archived["input_sample_diagnostics"][tr+"_"+role],
                origproto=origproto)
            (full_data if role=="dat" else random_pool)[tr]=full
            replay[tr+"_"+role]={
              "eligible_rows":len(full[0]),
              "A02_selected_array_SHA256":info["selected_array_SHA256"]
            }

        random_mask={}; random_meta={}
        for tr in TRACERS:
            random_mask[tr],pos=sample_fixed_random_mask(random_pool[tr],mid,cap,tr)
            random_meta[tr]={
              "eligible_pool_rows":len(random_pool[tr][0]),
              "selected_rows":N_RANDOM_MASK,
              "E58_R1_seed":e58_r1_seed(mid,cap,tr),
              "selected_positions_SHA256":arr_sha(pos),
              "selected_catalogue_SHA256":E55.catalogue_sha(random_mask[tr]),
            }

        probes,pmeta=pair_probes(full_data["eBOSS_LRG"],full_data["eBOSS_ELG"],
                                 mid,cap,NPROBE)
        tresults={}
        for tr in TRACERS:
            field,meta=contrast_field(full_data[tr],random_mask[tr],probes,mid,cap,tr,args.chunk)
            tresults[tr]={"metrics":metrics(field,probes),"input":meta}
            print("E59_TRACER_PASS",ck,tr,
                  "RHALF",tresults[tr]["metrics"]["split_half_pearson"],
                  "SIGN",tresults[tr]["metrics"]["split_half_sign_agreement"],
                  "RCUT",tresults[tr]["metrics"]["cutoff_192_vs_256_pearson"],flush=True)
            tresults[tr]["full_primary_los"]=los(field["full_primary"],probes).tolist()

        case.update({
          "source_replay":replay,"fixed_random_mask":random_meta,
          "probe_midpoints":pmeta,"tracers":tresults,
          "LRG_vs_ELG_full_primary":pearson(
              np.asarray(tresults["eBOSS_LRG"]["full_primary_los"]),
              np.asarray(tresults["eBOSS_ELG"]["full_primary_los"])),
          "status":"complete"
        })
        atomic(OUT,state)
        print("E59_CASE_PASS",ck,"CROSS",case["LRG_vs_ELG_full_primary"],flush=True)

    need(len(state["cases"])==18 and all(c["status"]=="complete" for c in state["cases"].values()),
         "Not all E59 cases complete")
    agg,candidates,selected=aggregate(state)
    state["aggregate"]=agg
    state["candidate_reconstructions"]=candidates
    state["decision"]={
      "selected_tracer_reconstruction_for_truth_calibration":selected,
      "repeatability_screen_any_tracer_pass":bool(selected is not None),
      "interpretation":(
        "MOCK_INTERNAL_RECONSTRUCTION_REPEATABILITY_PASS; NEXT_REQUIRE_TRUE_VELOCITY_CALIBRATION"
        if selected is not None else
        "MOCK_INTERNAL_RECONSTRUCTION_REPEATABILITY_FAIL; DO_NOT_USE_OBSERVED_TAG"
      ),
      "important_limit":(
        "Split-half repeatability is conditional on one fixed E58 R1 random selection field "
        "and is not correlation with true baryon/halo velocity. Shared survey geometry and "
        "large-scale modes can make it optimistic. A truth-labelled N-body/lightcone "
        "calibration, including reconstruction-mask robustness, is required before observations."
      )
    }
    state["status"]="PASS_MOCK_VELOCITY_TAG_REPEATABILITY_QUANTIFIED"
    atomic(OUT,state)

    print("E59_WSL_MOCK_VELOCITY_TAG_REPEATABILITY_PASS",flush=True)
    for cap in CAPS:
      for tr in TRACERS:
        z=state["aggregate"][cap][tr]
        print("CAP",cap,"TRACER",tr,
              "RHALF_MED",z["split_half_pearson"]["median"],
              "SIGN_MED",z["split_half_sign_agreement"]["median"],
              "RCUT_MED",z["cutoff_192_vs_256_pearson"]["median"],
              "IMPLIED_RREC_DIAGNOSTIC",z["implied_full_vs_latent_r_DIAGNOSTIC_ONLY"]["median"],
              flush=True)
    print("CANDIDATES",json.dumps(candidates,sort_keys=True),flush=True)
    print("SELECTED_FOR_TRUTH_CALIBRATION",selected,flush=True)
    print("DECISION",state["decision"]["interpretation"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
