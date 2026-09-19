#!/usr/bin/env python3
"""Build one linked z-resolved standard odd template for DESI LRG->ELG.

The radial basis is fixed by build_lrg_elg_zresolved_physical_basis.py.
Tracer coefficients are linked through Bonvin et al. 2023 Eq. (8), rather than
letting nu1 and wide-angle amplitudes float independently.

Default benchmark inputs:
  * b_LRG=2.0, b_ELG=1.3
  * s_LRG=1.0, s_ELG=area-weighted 0.40143
  * b_D,ELG=area-weighted -5.47143 for z_eff<1 and 0 above z=1
These come from Friedman-Shaw et al. JCAP 03 (2025) 059 for DESI 0.8<z<1.
Evolution bias is read from a DR1-specific calibration JSON produced by
estimate_lrg_elg_dr1_evolution_bias.py.

This remains a benchmark because the published b,s,b_D inputs are not a
full DR1 per-redshift calibration and no DESI window convolution is applied.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import build_lrg_elg_physical_odd_basis as phys


B_LRG=2.0
B_ELG=1.3
S_LRG=1.0
S_ELG=(5000.0*0.44+9000.0*0.38)/14000.0
BD_ELG_08_10=(5000.0*(-4.7)+9000.0*(-5.9))/14000.0


def cosine(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    den=float(np.sqrt(np.dot(a,a)*np.dot(b,b)))
    return float(np.dot(a,b)/den) if den>0 else None


def lookup_fevo(summary,tracer,lo,hi):
    for row in summary[tracer]:
        if abs(float(row["zlo"])-lo)<1e-9 and abs(float(row["zhi"])-hi)<1e-9:
            return float(row["fevo_benchmark_median"]), float(row["fevo_fit_half_range"])
    raise KeyError(f"missing {tracer} evolution bias for {lo}-{hi}")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--basis",required=True)
    ap.add_argument("--basis-summary",required=True)
    ap.add_argument("--evolution-summary",required=True)
    ap.add_argument("--outdir",required=True)
    ap.add_argument("--b-lrg",type=float,default=B_LRG)
    ap.add_argument("--b-elg",type=float,default=B_ELG)
    ap.add_argument("--s-lrg",type=float,default=S_LRG)
    ap.add_argument("--s-elg",type=float,default=S_ELG)
    ap.add_argument("--bd-elg-lowz",type=float,default=BD_ELG_08_10)
    ap.add_argument("--bd-transition-z",type=float,default=1.0)
    args=ap.parse_args()

    a=np.genfromtxt(args.basis,delimiter=",",names=True)
    bs=json.loads(Path(args.basis_summary).read_text())
    es=json.loads(Path(args.evolution_summary).read_text())

    scale_nu=float(bs["global_normalization_scales"]["relativistic_nu1"])
    scale_wa=float(bs["global_normalization_scales"]["wide_angle"])
    nu=np.asarray(a["relativistic_nu1_shape"],float)*scale_nu
    wa=np.asarray(a["wide_angle_shape"],float)*scale_wa
    wake=np.asarray(a["wake_shape"],float)

    zlo=np.asarray(a["zlo"],float); zhi=np.asarray(a["zhi"],float)
    ze=np.asarray(a["z_effective"],float); s=np.asarray(a["s_Mpc_over_h"],float)

    linked=np.zeros(len(a),float)
    linked_no_bd=np.zeros(len(a),float)
    rows_meta=[]
    for meta in bs["z_bins"]:
        lo=float(meta["zlo"]); hi=float(meta["zhi"]); z=float(meta["z_effective"])
        m=(np.abs(zlo-lo)<1e-9)&(np.abs(zhi-hi)<1e-9)
        if not np.any(m): raise RuntimeError(f"basis lacks z bin {lo}-{hi}")
        f=float(meta["growth_f_reference"])
        geom=meta["geometry"]
        fe_lrg,err_lrg=lookup_fevo(es,"LRG",lo,hi)
        fe_elg,err_elg=lookup_fevo(es,"ELG",lo,hi)

        bd_elg=float(args.bd_elg_lowz) if z < float(args.bd_transition_z) else 0.0
        # Friedman-Shaw velocity-response convention maps to Bonvin as
        # f_evo,eff = f_evo - b_D.
        fe_lrg_eff=fe_lrg
        fe_elg_eff=fe_elg-bd_elg

        Arel,Awa=phys.physical_coefficients(
            args.b_lrg,args.b_elg,args.s_lrg,args.s_elg,
            fe_lrg_eff,fe_elg_eff,f,geom
        )
        Arel0,Awa0=phys.physical_coefficients(
            args.b_lrg,args.b_elg,args.s_lrg,args.s_elg,
            fe_lrg,fe_elg,f,geom
        )
        linked[m]=(geom["Hconf_h_over_Mpc"]/geom["H0_h_over_Mpc"])*Arel*nu[m]+Awa*wa[m]
        linked_no_bd[m]=(geom["Hconf_h_over_Mpc"]/geom["H0_h_over_Mpc"])*Arel0*nu[m]+Awa0*wa[m]
        rows_meta.append({
            "zlo":lo,"zhi":hi,"z_effective":z,
            "b_lrg":args.b_lrg,"b_elg":args.b_elg,
            "s_lrg":args.s_lrg,"s_elg":args.s_elg,
            "fevo_lrg":fe_lrg,"fevo_elg":fe_elg,
            "fevo_half_range_lrg":err_lrg,"fevo_half_range_elg":err_elg,
            "bD_elg":bd_elg,
            "fevo_eff_elg":fe_elg_eff,
            "growth_f":f,
            "A_relativistic":Arel,"A_wide_angle":Awa,
            "A_relativistic_no_bd":Arel0,"A_wide_angle_no_bd":Awa0,
        })

    linked_n,scale=phys.normalized(linked)
    linked0_n,scale0=phys.normalized(linked_no_bd)

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    np.savetxt(
        out/"lrg_elg_zresolved_linked_standard.csv",
        np.column_stack([zlo,zhi,ze,s,wake,linked_n,linked0_n]),
        delimiter=",",
        header="zlo,zhi,z_effective,s_Mpc_over_h,wake_shape,standard_linked_shape,standard_linked_no_bd_shape",
        comments=""
    )

    summary={
        "scope":"DR1-specific z-resolved linked-standard benchmark",
        "source_model":"Bonvin et al. 2023 Eq. (8)",
        "external_bias_inputs":{
            "b_lrg":args.b_lrg,"b_elg":args.b_elg,
            "s_lrg":args.s_lrg,"s_elg":args.s_elg,
            "bD_elg_below_transition":args.bd_elg_lowz,
            "bD_transition_z":args.bd_transition_z,
        },
        "evolution_bias_source":args.evolution_summary,
        "per_bin":rows_meta,
        "global_scales":{"with_bd":scale,"without_bd":scale0},
        "unweighted_cosines":{
            "wake_vs_linked":cosine(wake,linked_n),
            "wake_vs_linked_no_bd":cosine(wake,linked0_n),
            "linked_with_vs_without_bd":cosine(linked_n,linked0_n),
        },
        "guardrail":(
            "This is a linked benchmark, not the final DR1 standard-odd prediction. "
            "The evolution bias is DR1-specific, but b,s and b_D use external 0.8<z<1 averages; "
            "the z>1 Doppler-bias benchmark is set to zero. Do not tune these inputs on the odd data."
        ),
    }
    (out/"linked_standard_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
