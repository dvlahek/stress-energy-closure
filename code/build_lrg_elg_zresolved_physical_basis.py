#!/usr/bin/env python3
"""Build globally normalized z-resolved wake and standard odd radial bases.

The measurement fixes the three z bins and their data effective-redshift proxies.
The hidden-state wake is evaluated independently at each z and then normalized
once over the full 18D vector, preserving its predicted relative redshift
evolution.  The Bonvin nu1 and leading wide-angle bases are treated similarly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_phase7_template as wpt
import build_lrg_elg_physical_odd_basis as phys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--measurement", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--fine-step", type=float, default=0.25)
    args = ap.parse_args()

    meas = np.genfromtxt(args.measurement, delimiter=",", names=True)
    zlo = np.asarray(meas["zlo"], float)
    zhi = np.asarray(meas["zhi"], float)
    ze = np.asarray(meas["z_effective"], float)
    sdata = np.asarray(meas["s_Mpc_over_h"], float)

    bins = []
    for lo, hi in dict.fromkeys((float(a),float(b)) for a,b in zip(zlo,zhi)):
        m = (zlo == lo) & (zhi == hi)
        ss = sdata[m]
        if len(ss) < 2:
            raise RuntimeError("Need at least two separation bins per z slice")
        step = float(np.median(np.diff(ss)))
        edges = np.r_[ss[0]-0.5*step, ss+0.5*step]
        bins.append((lo,hi,float(np.mean(ze[m])),ss,edges))

    zgrid = np.asarray([b[2] for b in bins], float)
    base.ZBINS = zgrid

    q = np.linspace(0.0,20.0,4000)
    f0, weights, basis, Nnull, shapes, y, M = cro.kinetic_objects(q,args.mass,args.z_match)
    if shapes.shape[0] != len(wpt.COEFF):
        raise RuntimeError("production kinetic basis changed")
    raw = wpt.COEFF @ shapes
    shape = raw / float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
    fp = f0 + args.frac*shape
    fm = f0 - args.frac*shape

    out = Path(args.outdir)
    out.mkdir(parents=True,exist_ok=True)
    p0=out/"fd_reference.dat"; pp=out/"hidden_plus.dat"; pm=out/"hidden_minus.dat"
    base.write_psd(p0,q,f0); base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
    st0=base.build_state(p0,args.mass); stp=base.build_state(pp,args.mass); stm=base.build_state(pm,args.mass)

    rows=[]
    raw_w=[]; raw_nu=[]; raw_wa=[]
    metadata=[]
    for lo,hi,z,ss,edges in bins:
        info=st0["z"][z]
        k=np.asarray(base.KOBS,float); pk=np.asarray(info["pk"],float)
        geom=phys.class_geometry(p0,args.mass,z)
        hid=wpt.signed_hidden_response(q,f0,fp,fm,st0,stp,stm,args.mass,z)

        sfine=np.arange(edges[0],edges[-1]+0.5*args.fine_step,args.fine_step)
        wf=wpt.hankel_dipole(k,pk,hid,sfine)
        nf=phys.trapz_hankel(k,pk,k*geom["H0_h_over_Mpc"],1,sfine)
        muf=phys.trapz_hankel(k,pk,k*k,2,sfine)
        waf=(sfine/geom["r_Mpc_over_h"])*muf

        wb=phys.volume_bin_average(sfine,wf,edges)
        nb=phys.volume_bin_average(sfine,nf,edges)
        wab=phys.volume_bin_average(sfine,waf,edges)
        raw_w.extend(wb.tolist()); raw_nu.extend(nb.tolist()); raw_wa.extend(wab.tolist())
        metadata.append({"zlo":lo,"zhi":hi,"z_effective":z,"geometry":geom})

    raw_w=np.asarray(raw_w,float); raw_nu=np.asarray(raw_nu,float); raw_wa=np.asarray(raw_wa,float)
    wn,sw=phys.normalized(raw_w); nn,sn=phys.normalized(raw_nu); wan,swa=phys.normalized(raw_wa)

    i=0
    for lo,hi,z,ss,edges in bins:
        for s in ss:
            rows.append((lo,hi,z,float(s),wn[i],nn[i],wan[i])); i+=1

    np.savetxt(
        out/"lrg_elg_zresolved_physical_basis.csv",
        np.asarray(rows,float),
        delimiter=",",
        header="zlo,zhi,z_effective,s_Mpc_over_h,wake_shape,relativistic_nu1_shape,wide_angle_shape",
        comments="",
    )

    mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
    mis=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
    def cosine(a,b):
        den=float(np.sqrt(np.dot(a,a)*np.dot(b,b)))
        return float(np.dot(a,b)/den) if den>0 else None

    summary={
        "scope":"Globally normalized 18D z-resolved physical basis for DESI DR1 LRGxELG",
        "z_bins":metadata,
        "global_normalization_scales":{"wake":sw,"relativistic_nu1":sn,"wide_angle":swa},
        "flattened_unweighted_template_cosines":{
            "wake_vs_relativistic_nu1":cosine(wn,nn),
            "wake_vs_wide_angle":cosine(wn,wan),
            "relativistic_nu1_vs_wide_angle":cosine(nn,wan),
        },
        "max_relative_moment_mismatch":float(np.max(mis)),
        "guardrail":(
            "One normalization is used over the full (z,s) vector, so relative redshift evolution "
            "is preserved. The nu1 and wide-angle columns are radial-basis diagnostics until tracer "
            "bias, magnification bias, evolution bias and the DESI window are linked physically."
        ),
    }
    (out/"basis_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
