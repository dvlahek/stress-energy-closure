#!/usr/bin/env python3
"""Blind ACT DR6 lensing fit of the pre-defined fixed-omega0 hidden mode.

This uses the public ACT DR6 lensing likelihood in lens-only mode, hence the
provided CMB-marginalized covariance and no primary-CMB likelihood corrections.
The hidden direction, mass and abundance control were fixed before inspecting
ACT data.  No ACT-driven optimization of the kinetic state is performed.

The CLASS text-output pp column is [L(L+1)]^2 C_L^{phiphi}/(2 pi).  ACT expects
C_L^{kappa kappa}; therefore C_L^{kk} = (2 pi / 4) pp.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

TRIM_LMAX = 2499  # source CLASS spectra cover ell=2..2500
CLASS_PP_TO_KK = 2.0 * math.pi / 4.0
VARIANTS = ("act_baseline", "actplanck_baseline")


def read_table(path: Path):
    header=[]
    with path.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                break
    arr=np.loadtxt(path)
    if arr.ndim==1:
        arr=arr[None,:]
    import re
    labels={}
    for line in header:
        for num,name in re.findall(r"(\d+):([^\s]+)",line):
            labels[name.strip().lower()]=int(num)-1
    return arr,labels


def unique(root: Path, name: str) -> Path:
    hits=list(root.rglob(name))
    if len(hits)!=1:
        raise RuntimeError(f"Expected one {name} below {root}, found {hits}")
    return hits[0]


def load_clkk(root: Path, stem: str):
    p=unique(root,f"{stem}_00_cl.dat")
    a,labels=read_table(p)
    j=None
    for key in ("pp","phiphi","phi-phi","phi_phi"):
        if key in labels:
            j=labels[key]; break
    if j is None:
        if a.shape[1] >= 6:
            j=5
        else:
            raise RuntimeError(f"No pp column in {p}; labels={labels}")
    ell=a[:,0].astype(int)
    pp=a[:,j]
    if ell[0] > 2 or ell[-1] < TRIM_LMAX+1:
        raise RuntimeError(f"Need ell through {TRIM_LMAX+1}; got {ell[0]}..{ell[-1]}")
    return ell, pp*CLASS_PP_TO_KK


def binned_theory(alike, data, ell, clkk):
    s=alike.standardize(ell,clkk,TRIM_LMAX)
    b=data["binmat_act"] @ s
    if data["include_planck"]:
        b=np.append(b,data["binmat_planck"] @ s)
    return b


def chi2(data_vec, theory, cinv):
    r=data_vec-theory
    return float(r @ cinv @ r)


def profile_linear_nuisance(data_vec, theory, cols, cinv):
    if not cols:
        return chi2(data_vec,theory,cinv), []
    B=np.column_stack(cols)
    r=data_vec-theory
    N=B.T @ cinv @ B
    g=B.T @ cinv @ r
    beta=np.linalg.pinv(N,rcond=1e-12) @ g
    rr=r-B @ beta
    return float(rr @ cinv @ rr), beta.tolist()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--spectra-root",type=Path,required=True)
    ap.add_argument("--act-like-root",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)

    import sys
    sys.path.insert(0,str(args.act_like_root.resolve()))
    import act_dr6_lenslike as alike

    spectra={}
    for stem in ("fiducial","winner_plus","winner_minus"):
        spectra[stem]=load_clkk(args.spectra_root,stem)
    ell=spectra["fiducial"][0]
    for e,_ in spectra.values():
        if not np.array_equal(ell,e):
            raise RuntimeError("ell grids differ")

    out={
        "experiment":"ACT DR6 real-data blind hidden-mode lensing fit",
        "source_class_spectra_run_id":34677322082,
        "hidden_template":{
            "mass_eV":0.60,
            "control":"fixed present-day relic omega_ncdm",
            "alpha_definition":"alpha=+1/-1 are the pre-defined +/-30% winner; alpha=0 is the midpoint",
            "physical_validated_range":[-1.0,1.0],
            "no_ACT_data_driven_mass_or_direction_optimization":True,
        },
        "act_likelihood":{
            "mode":"lens_only",
            "like_corrections":False,
            "trim_lmax":TRIM_LMAX,
            "class_pp_to_clkk_factor":CLASS_PP_TO_KK,
        },
        "variants":{},
    }

    for variant in VARIANTS:
        data=alike.load_data(variant,lens_only=True,like_corrections=False,
                             trim_lmax=TRIM_LMAX)
        y=np.asarray(data["data_binned_clkk"],float)
        cinv=np.asarray(data["cinv"],float)
        theory={k:binned_theory(alike,data,*v) for k,v in spectra.items()}
        t0=theory["fiducial"]
        tp=theory["winner_plus"]
        tm=theory["winner_minus"]
        a=0.5*(tp-tm)
        even=0.5*(tp+tm)-t0
        fisher=float(a @ cinv @ a)
        sigma_alpha=float(1.0/math.sqrt(fisher)) if fisher>0 else float("inf")
        r=y-t0
        alpha_hat=float((a @ cinv @ r)/fisher) if fisher>0 else float("nan")
        alpha_phys=float(np.clip(alpha_hat,-1.0,1.0))
        tphys=t0+alpha_phys*a
        chi={
            "minus1_exact":chi2(y,tm,cinv),
            "zero_exact":chi2(y,t0,cinv),
            "plus1_exact":chi2(y,tp,cinv),
            "best_linear_physical":chi2(y,tphys,cinv),
        }
        cbest=min(chi["minus1_exact"],chi["zero_exact"],chi["plus1_exact"],chi["best_linear_physical"])

        # Conservative shape-only check: freely profile a single overall
        # convergence-spectrum amplitude around the midpoint theory.
        prof={}
        for label,t in (("minus1",tm),("zero",t0),("plus1",tp)):
            c,b=profile_linear_nuisance(y,t,[t0],cinv)
            prof[label]={"chi2":c,"amplitude_additive_coefficient":b[0]}
        pbest=min(v["chi2"] for v in prof.values())
        for v in prof.values():
            v["delta_chi2_from_best_three_point_profile"]=v["chi2"]-pbest

        endpoint_sep=float(math.sqrt(max((tp-tm) @ cinv @ (tp-tm),0.0)))
        even_norm=float(math.sqrt(max(even @ cinv @ even,0.0)))
        out["variants"][variant]={
            "n_bandpowers":int(y.size),
            "fixed_background":{
                "alpha_hat_linear_template":alpha_hat,
                "sigma_alpha_linear_template":sigma_alpha,
                "alpha_over_sigma":alpha_hat/sigma_alpha if np.isfinite(sigma_alpha) else None,
                "alpha_best_in_physical_range":alpha_phys,
                "chi2":chi,
                "delta_chi2":{k:v-cbest for k,v in chi.items()},
                "endpoint_separation_SN":endpoint_sep,
                "midpoint_even_nonlinearity_norm":even_norm,
            },
            "free_overall_lensing_amplitude_profile":prof,
            "scope":(
                "Actual lensing bandpowers with the official lens-only CMB-marginalized covariance. "
                "The fixed-background result is not a full cosmological-parameter marginalization; "
                "the free-amplitude profile is a conservative shape-only diagnostic."
            ),
        }

    p=args.outdir/"act_dr6_hidden_mode_fit.json"
    p.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":
    main()
