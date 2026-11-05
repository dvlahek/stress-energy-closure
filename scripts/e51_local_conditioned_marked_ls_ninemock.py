#!/usr/bin/env python3
"""
E51 mock-only velocity-conditioned marked Landy-Szalay feasibility test.

Scientific object
-----------------
For an external pair-midpoint sign tag tau(x_mid) in {-1,+1}, define the
tagged/marked cross-Landy-Szalay field

  xi_tau = [DD_tau/NDD - DR_tau/NDR - RD_tau/NRD + RR_tau/NRR]
           / [RR_1/NRR],

where each tagged pair count is sum(w_i w_j tau_mid) in the fixed (s,mu)
cell and RR_1 is the ordinary positive RR histogram.

This estimates a sign-weighted conditional cross-correlation without splitting
the random catalogue into +/- branches. The same external field is evaluated
for DD, DR, RD and RR, so survey geometry/selection from the tag is subtracted
by the four-term estimator. Under tracer reversal the midpoint tag is unchanged
while mu changes sign.

Scope
-----
- exact already-pinned 9 eBOSS EZmock realizations;
- exact original 600D/1200R deterministic samples;
- four fixed synthetic smooth external sign fields, chosen before this run;
- E19 F+ and F- angular injections, pair-level DD only, as TECHNICAL probes;
- NO observed galaxies, NO observed odd vector, NO physical velocity reconstruction,
  NO covariance inverse, NO p-values.

The synthetic external tag fields are NOT claimed to be neutrino-CDM velocities.
They test estimator background/noise and algebra. E49 separately supplies the
physical sign-retention factor expected from a baryon velocity proxy.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys, tempfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))

import audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport as A02
from audit_eboss_dr16_rr_pair_closure import cartesian
from eboss_dr16_fiducial import PRIMARY_GEOMETRY, comoving_mpc_over_h

IDS=A02.IDS
CAPS=A02.CAPS
TRACERS=A02.TRACERS
ROLES=A02.ROLES
S_EDGES=A02.S_EDGES
MU_EDGES=A02.MU_EDGES
TERMS=("D1D2","D1R2","R1D2","R1R2")

OLD=ROOT/"source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json"
OUT=ROOT/"source_data/e51_conditioned_marked_ls_ninemock_result.json"

TPLUS=(0.00022278876746842203,0.00014150296987877026)
TMINUS=(0.00031476164003239810,0.00021061456076441418)

# Four fixed smooth null-tag fields. These are technical geometry probes only.
# Each field: (n1, L1[Mpc/h], phi1, amp2, n2, L2[Mpc/h], phi2).
FIELDS={
 "A":((1.,2.,3.),900.,0.37,0.55,(-2.,1.,1.),600.,1.11),
 "B":((2.,-1.,2.),750.,0.83,0.55,(1.,3.,-2.),500.,2.07),
 "C":((-1.,2.,2.),1100.,1.43,0.55,(3.,1.,1.),650.,0.19),
 "D":((1.,-3.,2.),850.,2.41,0.55,(-2.,-1.,3.),550.,1.73),
}
FIELD_NAMES=tuple(FIELDS)

def need(c,msg):
    if not c: raise RuntimeError(msg)

def sha_array(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()

def unit(v):
    a=np.asarray(v,float); return a/np.linalg.norm(a)

def tag_values(midpoints):
    out={}
    for name,(v1,L1,p1,a2,v2,L2,p2) in FIELDS.items():
        n1=unit(v1); n2=unit(v2)
        u=(np.sin(2*np.pi*(midpoints@n1)/L1+p1)
           + a2*np.sin(2*np.pi*(midpoints@n2)/L2+p2))
        # measure-zero exact zero assigned +1 deterministically.
        out[name]=np.where(u>=0.,1.,-1.)
    return out

def P1(mu): return mu
def P3(mu): return 0.5*(5*mu**3-3*mu)

# Common normalization preserves the physical F-/F+ amplitude difference.
_grid=np.linspace(-1.,1.,200001)
_COMMON=max(abs(TPLUS[0]*P1(_grid)+TPLUS[1]*P3(_grid)))
need(_COMMON>0,"Invalid E19 common normalization")

def h_e19(mu,state):
    T=TPLUS if state=="plus" else TMINUS
    return (T[0]*P1(mu)+T[1]*P3(mu))/_COMMON

def clipped_odd_weight(ell):
    low=np.maximum(-1.,MU_EDGES[:-1])
    high=np.minimum(1.,MU_EDGES[1:])
    if ell==1:
        F=lambda z:z*z/2.
    elif ell==3:
        F=lambda z:(5./8.)*z**4-(3./4.)*z*z
    else:
        raise ValueError
    w=(2*ell+1)/2.*(F(high)-F(low))
    need(np.allclose(w,(-1)**ell*w[::-1],atol=1e-13,rtol=0),
         "Odd projection weights lost parity")
    return w
ODD_W={ell:clipped_odd_weight(ell) for ell in (1,3)}

def project12(xi):
    out=[]
    for ell in (1,3):
        out.extend((xi@ODD_W[ell]).tolist())
    return np.asarray(out,float)

def pair_histograms(first,second,distance,block=128):
    """
    Return ordinary, four tagged, and E19 DD-injection-basis histograms.
    Injection bases are returned for every term for code simplicity, but are
    used only for D1D2 in the marked-LS response.
    """
    pos1,dir1,w1=cartesian(first,distance)
    pos2,dir2,w2=cartesian(second,distance)
    r1sq=np.einsum("ij,ij->i",pos1,pos1)
    r2sq=np.einsum("ij,ij->i",pos2,pos2)
    nsep,nmu=len(S_EDGES)-1,len(MU_EDGES)-1
    ordinary=np.zeros((nsep,nmu),float)
    tagged={f:np.zeros((nsep,nmu),float) for f in FIELD_NAMES}
    inj={s:np.zeros((nsep,nmu),float) for s in ("plus","minus")}
    theta_cos_limit=np.cos(np.deg2rad(.05))
    accepted_total=0

    for i0 in range(0,len(pos1),block):
        i1=min(len(pos1),i0+block)
        a=pos1[i0:i1]
        dot=a@pos2.T
        s2=r1sq[i0:i1,None]+r2sq[None,:]-2*dot
        m2=r1sq[i0:i1,None]+r2sq[None,:]+2*dot
        np.maximum(s2,0,out=s2); np.maximum(m2,0,out=m2)
        denom=np.sqrt(s2*m2)
        mu=np.zeros_like(denom)
        good=denom>0
        np.divide(r2sq[None,:]-r1sq[i0:i1,None],denom,out=mu,where=good)
        cost=dir1[i0:i1]@dir2.T
        si=np.searchsorted(S_EDGES,np.sqrt(s2),side="right")-1
        mi=np.searchsorted(MU_EDGES,mu,side="right")-1
        ok=(good&(si>=0)&(si<nsep)&(mi>=0)&(mi<nmu)
            &(cost<=theta_cos_limit)&(cost>-1))
        if not np.any(ok): continue
        accepted_total += int(np.count_nonzero(ok))
        flat=si[ok]*nmu+mi[ok]
        pw=(w1[i0:i1,None]*w2[None,:])[ok]
        ordinary += np.bincount(flat,weights=pw,minlength=nsep*nmu).reshape(nsep,nmu)

        mid=0.5*(a[:,None,:]+pos2[None,:,:])
        midok=mid[ok]
        tags=tag_values(midok)
        for f in FIELD_NAMES:
            tagged[f]+=np.bincount(flat,weights=pw*tags[f],
                                   minlength=nsep*nmu).reshape(nsep,nmu)
        muv=mu[ok]
        for state in ("plus","minus"):
            inj[state]+=np.bincount(flat,weights=pw*h_e19(muv,state),
                                    minlength=nsep*nmu).reshape(nsep,nmu)

    norm=float(np.sum(w1,dtype="f8")*np.sum(w2,dtype="f8"))
    need(norm>0 and np.isfinite(norm),"Bad pair normalization")
    return ordinary,tagged,inj,{
        "accepted_pairs":accepted_total,
        "pair_normalization":norm,
        "ordinary_SHA256":sha_array(ordinary),
    }

def marked_xi(hist,norms):
    rr=hist["R1R2"]["ordinary"]/norms["R1R2"]
    need(np.all(rr>0),"E51 requires the archived 144/144 ordinary RR support")
    out={}
    for field in FIELD_NAMES:
        num=(hist["D1D2"]["tagged"][field]/norms["D1D2"]
             -hist["D1R2"]["tagged"][field]/norms["D1R2"]
             -hist["R1D2"]["tagged"][field]/norms["R1D2"]
             +hist["R1R2"]["tagged"][field]/norms["R1R2"])
        out[field]=num/rr
    return out

def injection_response(hist,norms,state):
    # For DD -> DD*(1+lambda*tau*h_F), the tagged DD numerator gains
    # lambda * sum(w*h_F) because tau^2=1.
    rr=hist["R1R2"]["ordinary"]/norms["R1R2"]
    q=(hist["D1D2"]["inj"][state]/norms["D1D2"])/rr
    return q

def synthetic_self_test():
    # Pure algebra: constant geometry with exact signed cancellation and injection.
    rng=np.random.default_rng(5102026)
    shape=(6,24)
    rr=rng.uniform(.4,1.2,size=shape)
    n={k:1. for k in TERMS}
    # identical tagged four terms => zero marked LS
    t=rng.normal(size=shape)
    h={k:{"ordinary":rr.copy(),"tagged":{f:t.copy() for f in FIELD_NAMES},
          "inj":{"plus":np.zeros(shape),"minus":np.zeros(shape)}} for k in TERMS}
    z=marked_xi(h,n)
    need(max(np.max(np.abs(v)) for v in z.values())<1e-14,
         "Synthetic marked-LS null failed")
    # Tracer reversal parity of the fixed external midpoint tag is encoded by
    # mu mirror only; odd projection must change sign.
    x=rng.normal(size=shape)
    for ell in (1,3):
        a=x@ODD_W[ell]; b=x[:,::-1]@ODD_W[ell]
        need(np.max(np.abs(a+b))<2e-13,"Synthetic odd reversal failed")
    print("E51_SYNTHETIC_MARKED_LS_SELF_TEST_PASS",flush=True)

def atomic(path,obj):
    raw=(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e51_",delete=False) as f:
        tmp=Path(f.name); f.write(raw); f.flush(); os.fsync(f.fileno())
    try: os.replace(tmp,path)
    finally: tmp.unlink(missing_ok=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--block",type=int,default=128)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    need(args.block>=16 and args.block<=512,"Unsafe block size")
    synthetic_self_test()
    if args.self_test: return

    p=json.loads(A02.PROTOCOL.read_text())
    g,r,gman,rman,ref,old0001,gproto,rproto,rawproto,origproto=A02.preflight(
        p,require_local=True)
    paths,total=A02.resolve_72_sources(p,g,r,gman,rman,ref,gproto,rproto,rawproto)
    archive=json.loads(OLD.read_text())
    need(archive["status"]==A02.PASS and archive["completed_cases"]==18
         and archive["observed_odd_data_vector_read"] is False,
         "Archived A02 nine-mock parent changed")

    if OUT.exists():
        result=json.loads(OUT.read_text())
        need(result.get("stage")=="E51_CONDITIONED_MARKED_LS_NINEMOCK"
             and result.get("observed_odd_used") is False,
             "Existing E51 checkpoint incompatible")
    else:
        result={
          "stage":"E51_CONDITIONED_MARKED_LS_NINEMOCK",
          "status":"INCOMPLETE",
          "mock_ids":list(IDS),"caps":list(CAPS),
          "external_tag_fields":{k:{
             "definition":"two fixed smooth sinusoidal Cartesian modes; technical null only",
             "parameters":list(v)} for k,v in FIELDS.items()},
          "estimator":"[DD_tau/NDD-DR_tau/NDR-RD_tau/NRD+RR_tau/NRR]/[RR/NRR]",
          "E19_common_angular_normalization":_COMMON,
          "cases":{},
          "observed_galaxy_rows_used":False,
          "observed_odd_used":False,
          "physical_velocity_reconstruction_used":False,
          "inferential_covariance":False,
        }
        atomic(OUT,result)

    distance=lambda z:comoving_mpc_over_h(z,PRIMARY_GEOMETRY)

    for mid in IDS:
      for cap in CAPS:
        ck=f"{mid:04d}/{cap}"
        if result["cases"].get(ck,{}).get("status")=="complete":
            print("E51_REUSE",ck,flush=True); continue
        parent=archive["cases"][ck]
        need(parent["status"]=="complete" and parent["cross_ls"]["full_RR_support"] is True,
             "A02 parent case incomplete "+ck)

        cats={}; metadata={}
        prior=next((z for z in old0001["cases"] if z["cap"]==cap),None) if mid==1 else None
        for tracer in TRACERS:
          kk=A02.key(mid,cap,tracer)
          for role in ROLES:
            path=paths[kk+"/"+role]
            rows=A02.source_header_rows(path,prior,tracer,role)
            cat,info=A02.sample_catalogue(path,expected_rows=rows,cap=cap,tracer=tracer,
                                          role=role,expected_highz=None,p=origproto)
            expected=parent["input_sample_diagnostics"][tracer+"_"+role]
            need(info["selected_array_SHA256"]==expected["selected_array_SHA256"],
                 "A02 selected catalogue SHA drift "+ck+"/"+tracer+"/"+role)
            cats[tracer,role]=cat; metadata[tracer+"_"+role]=info

        samples={
          "D1D2":(cats["eBOSS_LRG","dat"],cats["eBOSS_ELG","dat"]),
          "D1R2":(cats["eBOSS_LRG","dat"],cats["eBOSS_ELG","ran"]),
          "R1D2":(cats["eBOSS_LRG","ran"],cats["eBOSS_ELG","dat"]),
          "R1R2":(cats["eBOSS_LRG","ran"],cats["eBOSS_ELG","ran"]),
        }
        H={}; norms={}
        for term,(c1,c2) in samples.items():
            ordinary,tagged,inj,meta=pair_histograms(c1,c2,distance,args.block)
            expected=parent["cross_ls"]["forward_pair_terms"][term]
            need(meta["ordinary_SHA256"]==expected["weighted_histogram_SHA256"],
                 "E51 ordinary pair SHA fails exact A02 replay "+ck+"/"+term)
            need(meta["accepted_pairs"]==expected["accepted_pairs"],
                 "E51 accepted pair count fails A02 replay "+ck+"/"+term)
            need(abs(meta["pair_normalization"]/
                     expected["independently_normalized_pair_weight"]-1)<1e-12,
                 "E51 pair norm fails A02 replay "+ck+"/"+term)
            H[term]={"ordinary":ordinary,"tagged":tagged,"inj":inj}
            norms[term]=meta["pair_normalization"]

        base=marked_xi(H,norms)
        qplus=project12(injection_response(H,norms,"plus"))
        qminus=project12(injection_response(H,norms,"minus"))
        need(float(qplus@qplus)>0 and float(qminus@qminus)>0,"Zero injection response")

        cased={
          "status":"complete","id":mid,"cap":cap,
          "baseline_12d_by_tag_field":{f:project12(base[f]).tolist() for f in FIELD_NAMES},
          "unit_lambda_E19_plus_12d":qplus.tolist(),
          "unit_lambda_E19_minus_12d":qminus.tolist(),
          "plus_minus_postwindow_angle_rad":
             float(math.acos(np.clip(float(qplus@qminus)/
                    math.sqrt(float(qplus@qplus)*float(qminus@qminus)),-1,1))),
          "A02_selected_catalogues_exact_replay":True,
          "A02_four_forward_pair_histograms_exact_replay":True,
        }

        # Actual new tag-specific reversal closure only for ID0001, both caps.
        if mid==1:
            revsamples={
              "D1D2":(cats["eBOSS_ELG","dat"],cats["eBOSS_LRG","dat"]),
              "D1R2":(cats["eBOSS_ELG","dat"],cats["eBOSS_LRG","ran"]),
              "R1D2":(cats["eBOSS_ELG","ran"],cats["eBOSS_LRG","dat"]),
              "R1R2":(cats["eBOSS_ELG","ran"],cats["eBOSS_LRG","ran"]),
            }
            rev={}
            for term,(c1,c2) in revsamples.items():
                ordinary,tagged,inj,meta=pair_histograms(c1,c2,distance,args.block)
                rev[term]={"ordinary":ordinary,"tagged":tagged,"inj":inj}
            # Mapping under tracer reversal.
            mapping={"D1D2":"D1D2","D1R2":"R1D2","R1D2":"D1R2","R1R2":"R1R2"}
            mx=0.
            for term,rt in mapping.items():
              for f in FIELD_NAMES:
                mx=max(mx,float(np.max(np.abs(
                    H[term]["tagged"][f]-rev[rt]["tagged"][f][:,::-1]))))
            need(mx<5e-12,"Tagged midpoint reversal closure failed "+ck)
            cased["tagged_midpoint_reversal_max_abs"]=mx

        result["cases"][ck]=cased
        atomic(OUT,result)
        print("E51_CASE_PASS",ck,"ANGLE",cased["plus_minus_postwindow_angle_rad"],flush=True)

    need(len(result["cases"])==18 and all(v["status"]=="complete" for v in result["cases"].values()),
         "Not all E51 cases complete")

    # Descriptive, no covariance inverse: apparent lambda background per
    # tag field and template, and post-window F+/F- angle.
    summaries={}
    for cap in CAPS:
      summaries[cap]={}
      for field in FIELD_NAMES:
        rr=[v for v in result["cases"].values() if v["cap"]==cap]
        rec={}
        for state,key in (("plus","unit_lambda_E19_plus_12d"),
                          ("minus","unit_lambda_E19_minus_12d")):
            vals=[]
            for c in rr:
                b=np.asarray(c["baseline_12d_by_tag_field"][field],float)
                q=np.asarray(c[key],float)
                vals.append(float(q@b/(q@q)))
            rec[state]={
              "baseline_apparent_lambda":vals,
              "mean":float(np.mean(vals)),
              "sample_sd":float(np.std(vals,ddof=1)),
              "median":float(np.median(vals)),
              "MAD":float(np.median(np.abs(vals-np.median(vals)))),
            }
        summaries[cap][field]=rec
      angles=[v["plus_minus_postwindow_angle_rad"] for v in result["cases"].values()
              if v["cap"]==cap]
      summaries[cap]["postwindow_angle_rad"]={
        "min":float(min(angles)),"median":float(np.median(angles)),"max":float(max(angles))}
    result["descriptive_summary"]=summaries
    result["status"]="PASS_MOCK_ONLY_CONDITIONED_ESTIMATOR_BACKGROUND_AND_E19_INJECTION_BASIS"
    result["interpretation_guardrails"]={
      "tag_fields_are_physical_velocity_reconstructions":False,
      "E19_injection_is_absolute_eBOSS_prediction":False,
      "nine_mocks_are_sufficient_for_12d_covariance_inverse":False,
      "detection_significance_or_pvalue":False,
      "purpose":"estimator null/background scatter and technical E19 injection response before any real velocity-tag field"
    }
    atomic(OUT,result)
    print("E51_WSL_CONDITIONED_MOCK_ESTIMATOR_PASS",flush=True)
    for cap in CAPS:
      print("CAP",cap,"ANGLE",summaries[cap]["postwindow_angle_rad"],flush=True)
      for f in FIELD_NAMES:
        print("TAG",f,
              "PLUS_SD",summaries[cap][f]["plus"]["sample_sd"],
              "MINUS_SD",summaries[cap][f]["minus"]["sample_sd"],flush=True)
    print("OBSERVED_ODD_USED",False,flush=True)
    print("REPORT",OUT,flush=True)

if __name__=="__main__":
    main()
