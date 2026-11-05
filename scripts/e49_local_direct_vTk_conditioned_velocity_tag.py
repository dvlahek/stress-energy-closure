#!/usr/bin/env python3
"""
E49: direct-CLASS velocity-tag feasibility for the conditioned Einstein-Vlasov estimator.

This is a SOURCE-ONLY screening calculation. It asks whether an independently
observable baryon/electron bulk-velocity tag could preserve the sign of the
neutrino-CDM relative velocity that controls the E10/E19 conditional wake.

For each frozen F0/F+/F- state at z=0.95 and R=16 Mpc/h:
  1. rerun CLASS with direct vTk,
  2. reproduce the frozen E21 R16 LOS sigma,
  3. compute r(v_nu-v_cdm, v_b) and r(v_nu-v_cdm, v_cdm),
  4. convert Gaussian field correlation r into ideal sign-tag fidelity
         E[sign X sign Y] = 2/pi asin(r),
  5. optionally degrade by representative reconstruction correlations r_rec.

NO eBOSS FITS, NO ASDF, NO observed galaxy rows, NO observed odd vector.
"""
from pathlib import Path
import sys, json, math, hashlib
import numpy as np
from classy import Class

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))

import class_response_optimize as cro
import wake_two_tracer_fisher as base

E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E8_SHA="bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
CLASS_SHA="e85808324f51fc694d12e3ed7439552a3c3f9540"
OUT=ROOT/"source_data/e49_direct_vTk_conditioned_velocity_tag_feasibility.json"

Z=.95
MASS=.06
R16=16.0
C=299792.458
KPIV=.05
KVEL=np.geomspace(1e-4,.15,180,dtype=float)
CUTOFF=.1
RREC=(1.0,.9,.8,.7,.5)

STATES=("FD","plus","minus")
COLS={
    "FD":"F0_CLASS_normalized",
    "plus":"Fplus_CLASS_normalized",
    "minus":"Fminus_CLASS_normalized",
}
E21_SIGMA={
    "FD":135.67798143801986,
    "plus":135.86206681100845,
    "minus":135.48537162977118,
}

def need(c,msg):
    if not c:
        raise RuntimeError(msg)

def sha(b):
    return hashlib.sha256(b).hexdigest()

def tophat(x):
    x=np.asarray(x,float)
    out=np.ones_like(x)
    m=np.abs(x)>1e-5
    xm=x[m]
    out[m]=3.0*(np.sin(xm)-xm*np.cos(xm))/xm**3
    return out

def cumulative_logk(y,k):
    out=np.zeros_like(y)
    out[1:]=np.cumsum(.5*(y[1:]+y[:-1])*np.diff(np.log(k)))
    return out

def integ_to_cutoff(y):
    return float(np.interp(CUTOFF,KVEL,cumulative_logk(y,KVEL)))

def filtered_stats(a,b):
    w=tophat(KVEL*R16)
    aa=max(integ_to_cutoff(a*a*w*w),0.0)
    bb=max(integ_to_cutoff(b*b*w*w),0.0)
    ab=integ_to_cutoff(a*b*w*w)
    den=math.sqrt(max(aa*bb,1e-300))
    r=max(-1.0,min(1.0,ab/den))
    return {
        "r":r,
        "abs_r":abs(r),
        "sigma_a_3d_proxy_kms":math.sqrt(aa),
        "sigma_b_3d_proxy_kms":math.sqrt(bb),
        "sigma_a_LOS_kms":math.sqrt(aa/3.0),
        "sigma_b_LOS_kms":math.sqrt(bb/3.0),
        "regression_a_on_b":ab/max(bb,1e-300),
    }

def sign_factor(r):
    # Joint-zero-mean-Gaussian identity.
    r=max(-1.0,min(1.0,float(r)))
    return 2.0/math.pi*math.asin(r)

def write_psd(path,q,f):
    path.parent.mkdir(parents=True,exist_ok=True)
    txt="\n".join(f"{qi:.14e} {fi:.14e}" for qi,fi in zip(q,f))+"\n"
    path.write_text(txt,encoding="utf-8")

def velocity_fields(q,f,state):
    work=ROOT/"eboss_workspace/a03_physics_source/e49_direct_velocity_tag"
    psd=work/f"frozen_{state}.dat"
    write_psd(psd,q,f)

    base.ZBINS=np.array([.9,1.],dtype=float)
    pars=base.class_params(psd,MASS)
    pars["output"]="mPk,dTk,vTk"
    need(pars["gauge"]=="newtonian","Unexpected CLASS gauge")
    need(cro.CLASS_COMMIT==CLASS_SHA,"CLASS source pin changed")

    c=Class(); c.set(pars); c.compute()
    try:
        h=float(c.h())
        t=c.get_transfer(z=Z,output_format="class")
        keys=set(t.keys())
        required={"k (h/Mpc)","t_ncdm[0]","t_cdm"}
        need(required.issubset(keys),"Missing required direct vTk columns: "+repr(sorted(required-keys)))
        # Prefer direct baryon theta. Fail closed instead of silently swapping
        # to a density-derivative proxy.
        need("t_b" in keys,
             "CLASS vTk output has no t_b. STOP: do not replace with density derivative in E49.")
        kin=np.asarray(t["k (h/Mpc)"],float)
        need(kin[0]<=KVEL[0] and kin[-1]>=KVEL[-1],"CLASS k grid does not cover E49 KVEL")

        def interp(key):
            return np.interp(KVEL,kin,np.asarray(t[key],float))

        tn=interp("t_ncdm[0]")
        tc=interp("t_cdm")
        tb=interp("t_b")

        k_mpc=KVEL*h
        prim=np.sqrt(float(cro.A_S)*(k_mpc/KPIV)**(float(base.NS)-1.0))
        common=-C*prim/k_mpc

        vnu=common*tn
        vcdm=common*tc
        vb=common*tb
        vrel=vnu-vcdm

        return h,vrel,vb,vcdm,sorted(keys)
    finally:
        c.struct_cleanup(); c.empty()

def main():
    need(E8.is_file(),"MISSING "+str(E8))
    need(sha(E8.read_bytes())==E8_SHA,"Frozen E8 CSV SHA changed")

    a=np.genfromtxt(E8,delimiter=",",names=True)
    q=np.asarray(a["q_dimensionless"],float)
    need(len(q)==4000,"Frozen E8 q grid changed")

    results={}
    for state in STATES:
        f=np.asarray(a[COLS[state]],float)
        need(np.min(f)>0 and np.isfinite(f).all(),"Invalid frozen "+state)
        h,vrel,vb,vcdm,keys=velocity_fields(q,f,state)

        rel_b=filtered_stats(vrel,vb)
        rel_c=filtered_stats(vrel,vcdm)

        sigma=rel_b["sigma_a_LOS_kms"]
        gap=abs(sigma-E21_SIGMA[state])/E21_SIGMA[state]
        need(gap<2e-5,
             f"{state} direct-vTk sigma does not replay E21: rel gap {gap}")

        bsign=sign_factor(rel_b["r"])
        csign=sign_factor(rel_c["r"])
        degraded={}
        for rr in RREC:
            reff=rel_b["r"]*rr
            degraded[f"{rr:.2f}"]={
                "assumed_reconstruction_corr_with_true_baryon_velocity":rr,
                "effective_r_rel_tag":reff,
                "gaussian_sign_correlation_factor":sign_factor(reff),
                "ideal_correct_sign_probability":
                    0.5*(1.0+abs(sign_factor(reff))),
            }

        results[state]={
            "h":h,
            "E21_sigma_LOS_target_kms":E21_SIGMA[state],
            "direct_sigma_LOS_replay_kms":sigma,
            "E21_sigma_relative_gap":gap,
            "relative_vs_baryon":rel_b,
            "relative_vs_cdm":rel_c,
            "ideal_baryon_sign_factor":bsign,
            "ideal_baryon_correct_sign_probability":0.5*(1+abs(bsign)),
            "ideal_cdm_sign_factor":csign,
            "ideal_cdm_correct_sign_probability":0.5*(1+abs(csign)),
            "reconstruction_degradation":degraded,
            "CLASS_transfer_keys":keys,
        }

    abs_rb=[results[s]["relative_vs_baryon"]["abs_r"] for s in STATES]
    abs_sign=[abs(results[s]["ideal_baryon_sign_factor"]) for s in STATES]
    out={
        "stage":"E49_DIRECT_VTK_CONDITIONED_VELOCITY_TAG_FEASIBILITY",
        "date":"2026-09-30",
        "status":"PASS_SOURCE_ONLY_ALIGNMENT_SCREEN",
        "physics":{
            "z_mathematical":Z,
            "R_Mpc_over_h":R16,
            "k_h_range":[float(KVEL[0]),float(CUTOFF)],
            "mass_eV":MASS,
            "tag":"direct CLASS baryon bulk velocity; prospective proxy for kSZ/velocity reconstruction",
            "target":"direct CLASS neutrino-CDM relative velocity",
            "gaussian_sign_identity":"E[sign X sign Y]=(2/pi) asin(r)",
        },
        "states":results,
        "summary":{
            "min_abs_r_rel_baryon_across_F_states":min(abs_rb),
            "median_abs_r_rel_baryon_across_F_states":float(np.median(abs_rb)),
            "max_abs_r_rel_baryon_across_F_states":max(abs_rb),
            "min_abs_ideal_sign_factor_across_F_states":min(abs_sign),
            "median_abs_ideal_sign_factor_across_F_states":float(np.median(abs_sign)),
            "interpretation_gate":{
                "strong":"min |r_rel,b| >= 0.7 -> external baryon velocity is a promising direction tag; proceed to conditioned estimator mock design",
                "intermediate":"0.4 <= min |r_rel,b| < 0.7 -> conditioned route possible but tag noise is a major loss",
                "weak":"min |r_rel,b| < 0.4 -> baryon/kSZ sign tag is a poor proxy; do not build survey estimator around it"
            }
        },
        "guardrails":{
            "observed_odd_used":False,
            "observed_galaxy_rows_used":False,
            "FITS_used":False,
            "ASDF_used":False,
            "survey_velocity_reconstruction_performed":False,
            "kSZ_data_used":False,
            "physical_conditioned_eBOSS_covariance":False,
        }
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

    print("E49_WSL_DIRECT_VTK_ALIGNMENT_PASS")
    for state in STATES:
        z=results[state]
        print("STATE",state)
        print("E21_SIGMA_REL_GAP",z["E21_sigma_relative_gap"])
        print("R_REL_BARYON",z["relative_vs_baryon"]["r"])
        print("SIGN_FACTOR_BARYON",z["ideal_baryon_sign_factor"])
        print("IDEAL_CORRECT_SIGN_PROB",z["ideal_baryon_correct_sign_probability"])
        for rr in ("0.80","0.70","0.50"):
            d=z["reconstruction_degradation"][rr]
            print("RREC",rr,
                  "R_EFF",d["effective_r_rel_tag"],
                  "SIGN_FACTOR",d["gaussian_sign_correlation_factor"],
                  "P_CORRECT",d["ideal_correct_sign_probability"])
    print("NO_FITS_ASDF_OBSERVED_ODD",True)
    print("REPORT",OUT)

if __name__=="__main__":
    main()
