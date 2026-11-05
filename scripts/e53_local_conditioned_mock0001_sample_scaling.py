#!/usr/bin/env python3
"""
E53: sample-size scaling of the E51 conditioned marked estimator on mock0001.

Stage 1 (default --run):
  - exact full-eligible mock galaxies for ID0001, NGC+SGC;
  - exact original nested 4800 randoms per tracer/cap;
  - same four frozen E51 technical sign fields;
  - same E19 F+/F- technical DD injection basis;
  - checkpoint after every pair term.

Stage 2 (--run --include-48000):
  - additionally reuses the exact already-defined nested 48000 randoms;
  - checkpointed separately;
  - compares 4800 -> 48000 finite-random stability.

This is a mock-only engineering/scientific sensitivity test. No observed
galaxies, no observed odd vector, no physical velocity reconstruction,
no covariance inverse, no p-values, no retuning of fields after E51.
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
from audit_eboss_dr16_rr_kdtree_pilot import cartesian
from eboss_dr16_fiducial import PRIMARY_GEOMETRY

E51=ROOT/"source_data/e51_conditioned_marked_ls_ninemock_result.json"
OUT=ROOT/"source_data/e53_conditioned_mock0001_sample_scaling.json"

CAPS=("NGC","SGC")
TRACERS=V1.TRACERS
ROLES=V1.ROLES
TERMS=("D1D2","D1R2","R1D2","R1R2")
S_EDGES=V1.S_EDGES
MU_EDGES=V1.MU_EDGES
LEVELS=("original_4800","nested_48000")

TPLUS=(0.00022278876746842203,0.00014150296987877026)
TMINUS=(0.00031476164003239810,0.00021061456076441418)

FIELDS={
 "A":((1.,2.,3.),900.,0.37,0.55,(-2.,1.,1.),600.,1.11),
 "B":((2.,-1.,2.),750.,0.83,0.55,(1.,3.,-2.),500.,2.07),
 "C":((-1.,2.,2.),1100.,1.43,0.55,(3.,1.,1.),650.,0.19),
 "D":((1.,-3.,2.),850.,2.41,0.55,(-2.,-1.,3.),550.,1.73),
}
FIELD_NAMES=tuple(FIELDS)

def need(c,msg):
    if not c: raise RuntimeError(msg)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e53_",delete=False) as f:
        t=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(t,path)
    finally: t.unlink(missing_ok=True)

def sha_array(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def unit(v):
    a=np.asarray(v,float); return a/np.linalg.norm(a)

def tag_values(mid):
    out={}
    for name,(v1,L1,p1,a2,v2,L2,p2) in FIELDS.items():
        n1,n2=unit(v1),unit(v2)
        u=(np.sin(2*np.pi*(mid@n1)/L1+p1)
           +a2*np.sin(2*np.pi*(mid@n2)/L2+p2))
        out[name]=np.where(u>=0.,1.,-1.)
    return out

def P1(mu): return mu
def P3(mu): return 0.5*(5*mu**3-3*mu)

_g=np.linspace(-1.,1.,200001)
_COMMON=max(abs(TPLUS[0]*P1(_g)+TPLUS[1]*P3(_g)))

def h_e19(mu,state):
    T=TPLUS if state=="plus" else TMINUS
    return (T[0]*P1(mu)+T[1]*P3(mu))/_COMMON

def clipped_odd_weight(ell):
    lo=np.maximum(-1.,MU_EDGES[:-1]); hi=np.minimum(1.,MU_EDGES[1:])
    if ell==1: F=lambda z:z*z/2.
    elif ell==3: F=lambda z:(5./8.)*z**4-(3./4.)*z*z
    else: raise ValueError
    return (2*ell+1)/2.*(F(hi)-F(lo))
ODD={ell:clipped_odd_weight(ell) for ell in (1,3)}

def project12(xi):
    out=[]
    for ell in (1,3): out.extend((xi@ODD[ell]).tolist())
    return np.asarray(out,float)

def sparse_marked(first,second,distance,*,chunk,max_candidates):
    x1,u1,w1=cartesian(first,distance)
    x2,u2,w2=cartesian(second,distance)
    tree=cKDTree(x2,leafsize=32)
    radius=float(S_EDGES[-1])
    counts=tree.query_ball_point(x1,radius,return_length=True,workers=2)
    candidates=int(np.sum(counts,dtype=np.int64))
    need(candidates<=max_candidates,
         f"E53 candidate-pair budget exceeded: {candidates}>{max_candidates}")

    shape=(len(S_EDGES)-1,len(MU_EDGES)-1)
    ordinary=np.zeros(shape,float)
    tagged={f:np.zeros(shape,float) for f in FIELD_NAMES}
    inj={s:np.zeros(shape,float) for s in ("plus","minus")}
    accepted=0
    coslim=float(np.cos(np.deg2rad(.05)))

    for a0 in range(0,len(x1),chunk):
        a1=min(a0+chunk,len(x1))
        nbs=tree.query_ball_point(x1[a0:a1],radius,workers=2)
        cnt=np.fromiter((len(v) for v in nbs),dtype=np.int64,count=a1-a0)
        if not np.any(cnt): continue
        ii=np.repeat(np.arange(a0,a1,dtype=np.int64),cnt)
        jj=np.concatenate([np.asarray(v,dtype=np.int64) for v in nbs if len(v)])
        delta=x2[jj]-x1[ii]
        s=np.linalg.norm(delta,axis=1)
        mid2=x2[jj]+x1[ii]
        midnorm=np.linalg.norm(mid2,axis=1)
        den=s*midnorm
        good=den>0
        mu=np.zeros(len(s),float)
        np.divide(np.einsum("ij,ij->i",delta,mid2),den,out=mu,where=good)
        cost=np.einsum("ij,ij->i",u1[ii],u2[jj])
        si=np.searchsorted(S_EDGES,s,side="right")-1
        mi=np.searchsorted(MU_EDGES,mu,side="right")-1
        ok=(good&(si>=0)&(si<shape[0])&(mi>=0)&(mi<shape[1])
            &(cost<=coslim)&(cost>-1.))
        if not np.any(ok): continue

        accepted += int(np.count_nonzero(ok))
        flat=si[ok]*shape[1]+mi[ok]
        pw=w1[ii[ok]]*w2[jj[ok]]
        ordinary += np.bincount(flat,weights=pw,minlength=shape[0]*shape[1]).reshape(shape)

        # Cartesian midpoint itself, not midpoint2, for exact E51 tag definition.
        mid=0.5*mid2[ok]
        tv=tag_values(mid)
        for f in FIELD_NAMES:
            tagged[f]+=np.bincount(flat,weights=pw*tv[f],
                                   minlength=shape[0]*shape[1]).reshape(shape)
        muv=mu[ok]
        for state in ("plus","minus"):
            inj[state]+=np.bincount(flat,weights=pw*h_e19(muv,state),
                                    minlength=shape[0]*shape[1]).reshape(shape)

    norm=float(np.sum(w1,dtype="f8")*np.sum(w2,dtype="f8"))
    need(norm>0 and np.isfinite(norm),"Invalid pair norm")
    return {
      "ordinary":ordinary,
      "tagged":tagged,
      "inj":inj,
      "meta":{
        "candidate_neighbour_pairs":candidates,
        "accepted_pairs":accepted,
        "pair_normalization":norm,
        "ordinary_SHA256":sha_array(ordinary)
      }
    }

def encode_term(rec):
    return {
      "ordinary":rec["ordinary"].tolist(),
      "tagged":{f:rec["tagged"][f].tolist() for f in FIELD_NAMES},
      "inj":{s:rec["inj"][s].tolist() for s in ("plus","minus")},
      "meta":rec["meta"],
    }

def decode_term(rec):
    return {
      "ordinary":np.asarray(rec["ordinary"],float),
      "tagged":{f:np.asarray(rec["tagged"][f],float) for f in FIELD_NAMES},
      "inj":{s:np.asarray(rec["inj"][s],float) for s in ("plus","minus")},
      "meta":rec["meta"],
    }

def summarize_level(termmap):
    H={k:decode_term(v) for k,v in termmap.items()}
    norms={k:H[k]["meta"]["pair_normalization"] for k in TERMS}
    rr=H["R1R2"]["ordinary"]/norms["R1R2"]
    need(np.all(rr>0),"E53 ordinary RR lacks full 144-cell support")

    baseline={}
    for f in FIELD_NAMES:
        num=(H["D1D2"]["tagged"][f]/norms["D1D2"]
             -H["D1R2"]["tagged"][f]/norms["D1R2"]
             -H["R1D2"]["tagged"][f]/norms["R1D2"]
             +H["R1R2"]["tagged"][f]/norms["R1R2"])
        baseline[f]=project12(num/rr)

    q={}
    for state in ("plus","minus"):
        q[state]=project12(
            (H["D1D2"]["inj"][state]/norms["D1D2"])/rr
        )

    qp,qm=q["plus"],q["minus"]
    d=qm-qp
    rho=float(qp@qm/(qp@qp))
    shape=qm-rho*qp
    angle=float(math.acos(np.clip(float(qp@qm)/
                  math.sqrt(float(qp@qp)*float(qm@qm)),-1,1)))

    out={
      "full_RR_support":True,
      "q_plus_12d":qp.tolist(),
      "q_minus_12d":qm.tolist(),
      "postwindow_angle_rad":angle,
      "difference_norm_over_plus_norm":float(np.linalg.norm(d)/np.linalg.norm(qp)),
      "rho_minus_on_plus":rho,
      "shape_residual_norm_over_plus_norm":float(np.linalg.norm(shape)/np.linalg.norm(qp)),
      "by_tag_field":{}
    }
    for f,b in baseline.items():
        out["by_tag_field"][f]={
          "baseline_12d":b.tolist(),
          "plus_apparent_lambda":float(qp@b/(qp@qp)),
          "minus_apparent_lambda":float(qm@b/(qm@qm)),
          "calibrated_Fminus_minus_Fplus_coordinate":float(d@b/(d@d)),
          "shape_only_coordinate":float(shape@b/(shape@shape)),
          "baseline_L2":float(np.linalg.norm(b)),
        }
    return out

def e51_small_case(cap):
    e=json.loads(E51.read_text())
    c=e["cases"][f"0001/{cap}"]
    qp=np.asarray(c["unit_lambda_E19_plus_12d"],float)
    qm=np.asarray(c["unit_lambda_E19_minus_12d"],float)
    d=qm-qp
    rho=float(qp@qm/(qp@qp))
    shape=qm-rho*qp
    out={"postwindow_angle_rad":float(c["plus_minus_postwindow_angle_rad"]),
         "by_tag_field":{}}
    for f,bv in c["baseline_12d_by_tag_field"].items():
        b=np.asarray(bv,float)
        out["by_tag_field"][f]={
          "plus_apparent_lambda":float(qp@b/(qp@qp)),
          "minus_apparent_lambda":float(qm@b/(qm@qm)),
          "calibrated_Fminus_minus_Fplus_coordinate":float(d@b/(d@d)),
          "shape_only_coordinate":float(shape@b/(shape@shape)),
          "baseline_L2":float(np.linalg.norm(b)),
        }
    return out

def scaling_summary(small,big):
    out={}
    for f in FIELD_NAMES:
        a=small["by_tag_field"][f]; b=big["by_tag_field"][f]
        out[f]={}
        for key in ("plus_apparent_lambda","minus_apparent_lambda",
                    "calibrated_Fminus_minus_Fplus_coordinate",
                    "shape_only_coordinate","baseline_L2"):
            av=float(a[key]); bv=float(b[key])
            out[f][key]={
              "small_600D1200R":av,
              "new":bv,
              "abs_ratio_new_over_small":(
                  abs(bv)/abs(av) if abs(av)>1e-15 else None
              ),
              "absolute_change":bv-av,
            }
    return out

def synthetic_self_test():
    rng=np.random.default_rng(53009)
    c1=(rng.uniform(120,121,31),rng.uniform(10,11,31),
        rng.uniform(.91,.99,31),rng.uniform(.7,1.5,31))
    c2=(rng.uniform(120,121,33),rng.uniform(10,11,33),
        rng.uniform(.91,.99,33),rng.uniform(.7,1.5,33))
    dist=lambda z:np.asarray(z,float)*2800.
    rec=sparse_marked(c1,c2,dist,chunk=16,max_candidates=100000)
    need(rec["ordinary"].shape==(6,24) and rec["meta"]["pair_normalization"]>0,
         "E53 synthetic sparse marked pair test failed")
    print("E53_SYNTHETIC_SPARSE_MARKED_SELF_TEST_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--run",action="store_true")
    ap.add_argument("--include-48000",action="store_true")
    ap.add_argument("--chunk",type=int,default=64)
    ap.add_argument("--max-candidate-pairs-per-term",type=int,default=150000000)
    args=ap.parse_args()

    synthetic_self_test()
    if args.self_test: return
    need(args.run,"Use --run for mock-only source reads")
    need(16<=args.chunk<=256,"Unsafe chunk")
    need(args.max_candidate_pairs_per_term>=1_000_000,"Candidate budget too small")
    need(E51.exists(),"Missing E51 parent result")
    e51=json.loads(E51.read_text())
    need(e51["status"]=="PASS_MOCK_ONLY_CONDITIONED_ESTIMATOR_BACKGROUND_AND_E19_INJECTION_BASIS"
         and e51["observed_odd_used"] is False,
         "E51 parent changed")
    # Freeze exact E51 technical field definitions.
    for f in FIELD_NAMES:
        need(e51["external_tag_fields"][f]["parameters"]==json.loads(json.dumps(FIELDS[f])),
             "E51 tag field definition changed: "+f)

    p,a02,e4=V1.source_only_gate()

    # Full 72-source provenance before first FITS row.
    original_protocol=A02.load(A02.PROTOCOL)
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,p0=A02.preflight(
        original_protocol,require_local=True)
    paths,total=A02.resolve_72_sources(
        original_protocol,g,r,gman,rman,ref,gproto,rproto,rawproto)

    if OUT.exists():
        state=json.loads(OUT.read_text())
        need(state.get("stage")=="E53_CONDITIONED_MOCK0001_SAMPLE_SCALING"
             and state.get("observed_odd_used") is False,
             "Existing E53 checkpoint incompatible")
    else:
        state={
          "stage":"E53_CONDITIONED_MOCK0001_SAMPLE_SCALING",
          "status":"INCOMPLETE",
          "mock_id":1,
          "caps":list(CAPS),
          "levels_requested":["original_4800"],
          "same_E51_fields":True,
          "all_72_full_gzip_rehashed_before_any_FITS":True,
          "source_72_total_compressed_bytes":total,
          "cases":{},
          "observed_galaxy_rows_used":False,
          "observed_odd_used":False,
          "physical_velocity_reconstruction_used":False,
          "pvalue_or_detection_sigma":False,
        }
        atomic(OUT,state)

    requested=["original_4800"]+(["nested_48000"] if args.include_48000 else [])
    state["levels_requested"]=requested
    atomic(OUT,state)

    for cap in CAPS:
        ck=f"0001/{cap}"
        case=state["cases"].setdefault(ck,{
          "status":"INCOMPLETE","cap":cap,
          "small_E51_600D1200R":e51_small_case(cap),
          "source_full_eligible_diagnostics":{},
          "levels":{}
        })

        metadata=a02["cases"][ck]["input_sample_diagnostics"]
        prior=e4["cases"][ck]
        cats={}
        info={}
        for tracer in TRACERS:
            for role in ROLES:
                label=tracer+"_"+role
                src=paths[A02.key(1,cap,tracer)+"/"+role]
                cats[(tracer,role)],info[label]=V1.load_original_full(
                    src,cap,tracer,role,metadata[label],prior,p0)
                print("E53_SOURCE_OK",ck,label,
                      "FULL_ELIGIBLE",info[label]["full_eligible_rows"],flush=True)
        case["source_full_eligible_diagnostics"]=info
        dist=V1.cached_exact_distance(cats)

        # Mandatory old 600D/1200R sparse/dense replay before new full-D marked counts.
        bench=V1.validate_sparse_against_original_dense(cats,a02["cases"][ck],dist)
        case["original_600D1200R_sparse_dense_replay"]=bench
        atomic(OUT,state)
        print("E53_PARENT_PAIR_REPLAY_OK",ck,flush=True)

        dl=cats["eBOSS_LRG","dat"]["full"]
        de=cats["eBOSS_ELG","dat"]["full"]

        for level in requested:
            lev=case["levels"].setdefault(level,{"terms":{},"status":"INCOMPLETE"})
            rl=cats["eBOSS_LRG","ran"][level]
            re=cats["eBOSS_ELG","ran"][level]
            samples={
              "D1D2":(dl,de),
              "D1R2":(dl,re),
              "R1D2":(rl,de),
              "R1R2":(rl,re)
            }
            for term in TERMS:
                if term in lev["terms"]:
                    print("E53_REUSE_TERM",ck,level,term,flush=True)
                    continue
                rec=sparse_marked(*samples[term],dist,chunk=args.chunk,
                                  max_candidates=args.max_candidate_pairs_per_term)
                # Exact E4 4800 RR replay when applicable.
                if level=="original_4800" and term=="R1R2":
                    oldrr=prior["levels"]["nested_4800"]["baseline_original"]["forward_pair_terms"]["R1R2"]
                    need(rec["meta"]["accepted_pairs"]==oldrr["accepted_pairs"],
                         "E53 original 4800 RR pair count changed "+ck)
                    need(abs(float(rec["ordinary"].sum())/
                             oldrr["total_weighted_pairs_in_fixed_s_mu_bins"]-1)<1e-11,
                         "E53 original 4800 RR weighted sum changed "+ck)
                lev["terms"][term]=encode_term(rec)
                atomic(OUT,state)
                print("E53_TERM_PASS",ck,level,term,
                      "CAND",rec["meta"]["candidate_neighbour_pairs"],
                      "ACC",rec["meta"]["accepted_pairs"],flush=True)

            lev["summary"]=summarize_level(lev["terms"])
            lev["status"]="complete"
            if level=="original_4800":
                case["scaling_600D1200R_to_fullD4800R"]=scaling_summary(
                    case["small_E51_600D1200R"],lev["summary"])
            if level=="nested_48000" and "original_4800" in case["levels"] and \
                    case["levels"]["original_4800"].get("status")=="complete":
                case["scaling_fullD4800R_to_fullD48000R"]=scaling_summary(
                    case["levels"]["original_4800"]["summary"],lev["summary"])
            atomic(OUT,state)
            print("E53_LEVEL_PASS",ck,level,
                  "ANGLE",lev["summary"]["postwindow_angle_rad"],flush=True)

        case["status"]=(
          "complete_4800_48000" if "nested_48000" in requested
          and case["levels"].get("nested_48000",{}).get("status")=="complete"
          else "complete_4800"
        )
        atomic(OUT,state)

    state["status"]=(
      "PASS_FULLD_4800_AND_48000_CONDITIONED_SAMPLE_SCALING"
      if args.include_48000 else
      "PASS_FULLD_4800_CONDITIONED_SAMPLE_SCALING_STAGE1"
    )
    state["interpretation_guardrails"]={
      "single_mock_realization_only":True,
      "cannot_measure_mock_to_mock_scatter":True,
      "technical_tag_fields_not_physical_velocity":True,
      "E19_injection_not_absolute_eBOSS_prediction":True,
      "observed_data_unsealed":False,
      "purpose":"determine whether E51 background is strongly reduced by full eligible galaxies and denser frozen randoms before expanding conditioned estimator to all nine mocks"
    }
    atomic(OUT,state)

    print("E53_WSL_PASS",state["status"],flush=True)
    for cap in CAPS:
        c=state["cases"][f"0001/{cap}"]
        for level in requested:
            s=c["levels"][level]["summary"]
            print("CAP",cap,"LEVEL",level,
                  "ANGLE",s["postwindow_angle_rad"],
                  "DIFF_OVER_PLUS",s["difference_norm_over_plus_norm"])
            for f in FIELD_NAMES:
                z=s["by_tag_field"][f]
                print("TAG",f,
                      "CAL_DIFF",z["calibrated_Fminus_minus_Fplus_coordinate"],
                      "PLUS",z["plus_apparent_lambda"],
                      "MINUS",z["minus_apparent_lambda"],
                      "L2",z["baseline_L2"])
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
