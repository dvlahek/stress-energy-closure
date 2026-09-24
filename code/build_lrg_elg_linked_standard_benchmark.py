#!/usr/bin/env python3
"""Build an externally linked DESI LRG->ELG standard odd-dipole benchmark.

This is a geometry benchmark, not a final DR1 likelihood model.

Inputs:
  * pre-window basis from build_lrg_elg_physical_odd_basis.py
  * Bonvin et al. 2023 Eq. (8) for the linked nu1 + wide-angle dipole
  * Friedman-Shaw et al. JCAP 03 (2025) 059, 0.8<z<1.0 DESI ELG/LRG
    bias inputs for North and South, including the ELG Doppler bias b_D.

Friedman-Shaw convention adds b_D to the line-of-sight velocity response.
Bonvin Eq. (8) uses f_evo inside alpha with the opposite sign in that
response, so we map
    f_evo,eff = f_evo - b_D.
This convention is written explicitly to avoid silently folding selection
physics into an unconstrained nuisance amplitude.

The North/South signals are combined with the 5000/9000 deg^2 survey-area
weights used in Friedman-Shaw et al. for their DESI forecast.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import build_lrg_elg_physical_odd_basis as phys


FS2025 = {
    "north": {
        "area_deg2": 5000.0,
        "b_lrg": 2.0,
        "b_elg": 1.3,
        "s_lrg": 1.0,
        "s_elg": 0.44,
        "fevo_lrg": 8.1,
        "fevo_elg": -1.2,
        "bd_lrg": 0.0,
        "bd_elg": -4.7,
    },
    "south": {
        "area_deg2": 9000.0,
        "b_lrg": 2.0,
        "b_elg": 1.3,
        "s_lrg": 1.0,
        "s_elg": 0.38,
        "fevo_lrg": 10.2,
        "fevo_elg": -2.1,
        "bd_lrg": 0.0,
        "bd_elg": -5.9,
    },
}


def cosine(a, b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    den=float(np.sqrt(np.dot(a,a)*np.dot(b,b)))
    return float(np.dot(a,b)/den) if den>0 else None


def metric_cosine(a,b,P):
    a=np.asarray(a,float); b=np.asarray(b,float)
    aa=float(a@P@a); bb=float(b@P@b)
    return float((a@P@b)/np.sqrt(aa*bb)) if aa>0 and bb>0 else None


def linked_signal(row, f, geom, nu1_raw, wa_raw, include_bd=True):
    fe_lrg=float(row["fevo_lrg"])
    fe_elg=float(row["fevo_elg"])
    if include_bd:
        fe_lrg -= float(row["bd_lrg"])
        fe_elg -= float(row["bd_elg"])
    Arel,Awa=phys.physical_coefficients(
        float(row["b_lrg"]), float(row["b_elg"]),
        float(row["s_lrg"]), float(row["s_elg"]),
        fe_lrg, fe_elg, float(f), geom
    )
    sig=(geom["Hconf_h_over_Mpc"]/geom["H0_h_over_Mpc"])*Arel*nu1_raw + Awa*wa_raw
    return np.asarray(sig,float), {
        "A_relativistic":float(Arel),
        "A_wide_angle":float(Awa),
        "fevo_eff_lrg":float(fe_lrg),
        "fevo_eff_elg":float(fe_elg),
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--basis",required=True,
                    help="lrg_elg_physical_odd_basis.csv")
    ap.add_argument("--basis-summary",required=True,
                    help="basis_summary.json from the same build")
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--covariance",
                    help="optional broad-bin mock covariance npz for C^-1 metric cosines")
    args=ap.parse_args()

    a=np.genfromtxt(args.basis,delimiter=",",names=True)
    S=json.loads(Path(args.basis_summary).read_text())
    s=np.asarray(a["s_Mpc_over_h"],float)
    wake=np.asarray(a["wake_shape"],float)

    scales=S["normalization_scales"]
    nu1_raw=np.asarray(a["relativistic_nu1_shape"],float)*float(scales["relativistic_nu1"])
    wa_raw=np.asarray(a["wide_angle_shape"],float)*float(scales["wide_angle"])
    f=float(S["growth_f_reference"])
    geom=S["geometry"]

    signals={}
    meta={}
    for region,row in FS2025.items():
        sig,coeff=linked_signal(row,f,geom,nu1_raw,wa_raw,include_bd=True)
        sig0,coeff0=linked_signal(row,f,geom,nu1_raw,wa_raw,include_bd=False)
        signals[region]=sig
        signals[region+"_no_bd"]=sig0
        meta[region]={"inputs":row,"with_bd":coeff,"without_bd":coeff0}

    Atot=sum(FS2025[r]["area_deg2"] for r in ("north","south"))
    wN=FS2025["north"]["area_deg2"]/Atot
    wS=FS2025["south"]["area_deg2"]/Atot
    linked=wN*signals["north"]+wS*signals["south"]
    linked0=wN*signals["north_no_bd"]+wS*signals["south_no_bd"]

    def norm(x):
        y,scale=phys.normalized(x)
        return y,float(scale)

    nN,sN=norm(signals["north"])
    nS,sS=norm(signals["south"])
    nA,sA=norm(linked)
    nA0,sA0=norm(linked0)

    cols=np.column_stack([s,wake,nN,nS,nA,nA0])
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    csv=out/"lrg_elg_linked_standard_benchmark.csv"
    np.savetxt(
        csv,cols,delimiter=",",
        header=(
            "s_Mpc_over_h,wake_shape,"
            "standard_fs2025_north_shape,standard_fs2025_south_shape,"
            "standard_fs2025_area_shape,standard_fs2025_area_no_bd_shape"
        ),comments=""
    )

    cos={
        "unweighted":{
            "wake_vs_north":cosine(wake,nN),
            "wake_vs_south":cosine(wake,nS),
            "wake_vs_area":cosine(wake,nA),
            "wake_vs_area_no_bd":cosine(wake,nA0),
            "area_with_vs_without_bd":cosine(nA,nA0),
        }
    }

    covinfo=None
    if args.covariance:
        B=np.load(args.covariance,allow_pickle=False)
        C=np.asarray(B["cov_xi1"],float)
        nmock=int(np.asarray(B["n_mocks"]).item())
        p=C.shape[0]
        if p!=len(wake):
            raise RuntimeError(f"covariance dimension {p} != template length {len(wake)}")
        alpha=float((nmock-p-2)/(nmock-1)) if nmock>p+2 else None
        if alpha is None:
            raise RuntimeError("not enough mocks for Hartlap precision")
        P=alpha*np.linalg.pinv(C,rcond=1e-12)
        cos["covariance_metric"]={
            "hartlap_factor":alpha,
            "wake_vs_north":metric_cosine(wake,nN,P),
            "wake_vs_south":metric_cosine(wake,nS,P),
            "wake_vs_area":metric_cosine(wake,nA,P),
            "wake_vs_area_no_bd":metric_cosine(wake,nA0,P),
            "area_with_vs_without_bd":metric_cosine(nA,nA0,P),
        }
        covinfo={"path":args.covariance,"n_mocks":nmock}

    summary={
        "scope":"External linked-standard geometry benchmark for DESI LRG->ELG",
        "applicability":"Friedman-Shaw et al. inputs are averaged over 0.8<z<1.0; use as a benchmark near z~0.9, not as the final 0.8<z<1.1 DR1 model.",
        "source_model":"Bonvin et al. 2023 Eq. (8), midpoint-LOS dipole",
        "source_biases":"Friedman-Shaw et al., JCAP 03 (2025) 059, DESI ELG/LRG 0.8<z<1.0 North/South averages",
        "doppler_bias_mapping":"f_evo_eff = f_evo - b_D",
        "area_weights":{"north":wN,"south":wS},
        "growth_f":f,
        "geometry":geom,
        "region_coefficients":meta,
        "normalization_scales":{
            "north":sN,"south":sS,"area":sA,"area_no_bd":sA0
        },
        "template_cosines":cos,
        "covariance":covinfo,
        "output_csv":str(csv),
        "analysis_scope":(
            "We combine standard radial bases using externally specified tracer coefficients. "
            "The survey-specific model requires the corresponding tracer calibration and window response."
        ),
    }
    (out/"linked_standard_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
