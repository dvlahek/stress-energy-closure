#!/usr/bin/env python3
"""Estimate DR1-specific LRG/ELG evolution-bias benchmarks from observed N(z).

This is an input-calibration diagnostic for the linked odd-dipole model, not a
final selection-function measurement.  It uses WEIGHT (not FKP) and removes the
geometric shell-volume factor before differentiating the comoving number
density.  To avoid selecting a smoothing choice from the odd data, the reported
benchmark is the median over four predeclared smooth fits:
  dz = 0.01, 0.02 crossed with polynomial degree = 2, 3.

Definition:
    f_evo = -(1+z) d ln nbar / dz
which matches Bonvin et al. Eq. (5) and the b_e convention used by
Friedman-Shaw et al. 2025.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

import desi_dr1_lrg_elg_exact as exact


def read_zw(paths, zmin, zmax):
    zs=[]; ws=[]
    for filename in paths:
        with fits.open(filename,memmap=True) as hdul:
            d=hdul[1].data
            n=len(d)
            z=exact.col(d,"Z").astype("f8",copy=False)
            w=exact.col(d,"WEIGHT",default=np.ones(n)).astype("f8",copy=False)
            good=np.isfinite(z)&np.isfinite(w)&(w>0)&(z>=zmin)&(z<zmax)
            zs.append(np.asarray(z[good],float))
            ws.append(np.asarray(w[good],float))
    return np.concatenate(zs),np.concatenate(ws)


def shell_volume_per_sr(zedges):
    # Mpc/h distances; an overall constant survey area cancels from d ln nbar/dz.
    r=np.asarray(exact.COSMO.comoving_distance(zedges).value*exact.H,float)
    return (r[1:]**3-r[:-1]**3)/3.0


def one_fit(z,w,lo,hi,zeff,dz,degree):
    edges=np.arange(lo,hi+0.5*dz,dz)
    if edges[-1] < hi-1e-10:
        edges=np.r_[edges,hi]
    else:
        edges[-1]=hi
    counts,_=np.histogram(z,bins=edges,weights=w)
    vol=shell_volume_per_sr(edges)
    cen=0.5*(edges[:-1]+edges[1:])
    nbar=counts/np.maximum(vol,1e-300)
    good=(counts>0)&np.isfinite(nbar)&(nbar>0)
    if np.count_nonzero(good) < degree+2:
        raise RuntimeError(f"insufficient populated bins for z={lo}-{hi}, dz={dz}, deg={degree}")
    x=cen[good]-zeff
    y=np.log(nbar[good])
    p=np.polyfit(x,y,degree)
    dp=np.polyder(p)
    slope=float(np.polyval(dp,0.0))
    fevo=float(-(1.0+zeff)*slope)
    yhat=np.polyval(p,x)
    rms=float(np.sqrt(np.mean((y-yhat)**2)))
    return {
        "dz":float(dz),"degree":int(degree),"fevo":fevo,
        "dln_nbar_dz":slope,"log_nbar_fit_rms":rms,
        "fine_bin_count":int(np.count_nonzero(good)),
    }


def tracer_summary(z,w,bins):
    out=[]
    for b in bins:
        vals=[]
        for dz in (0.01,0.02):
            for degree in (2,3):
                vals.append(one_fit(z,w,b["zlo"],b["zhi"],b["z_effective"],dz,degree))
        fs=np.asarray([v["fevo"] for v in vals],float)
        out.append({
            **b,
            "fevo_benchmark_median":float(np.median(fs)),
            "fevo_fit_min":float(np.min(fs)),
            "fevo_fit_max":float(np.max(fs)),
            "fevo_fit_half_range":float(0.5*(np.max(fs)-np.min(fs))),
            "fits":vals,
            "weighted_object_count":float(np.sum(w[(z>=b["zlo"])&(z<b["zhi"])])),
            "raw_object_count":int(np.count_nonzero((z>=b["zlo"])&(z<b["zhi"]))),
        })
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--lrg-data",nargs="+",required=True)
    ap.add_argument("--elg-data",nargs="+",required=True)
    ap.add_argument("--measurement",required=True,
                    help="z-resolved measurement CSV defining frozen bins and z_eff")
    ap.add_argument("--out",required=True)
    args=ap.parse_args()

    a=np.genfromtxt(args.measurement,delimiter=",",names=True)
    bins=[]
    seen=set()
    for lo,hi,ze in zip(a["zlo"],a["zhi"],a["z_effective"]):
        key=(float(lo),float(hi))
        if key in seen: continue
        seen.add(key)
        m=(np.asarray(a["zlo"],float)==float(lo))&(np.asarray(a["zhi"],float)==float(hi))
        bins.append({"zlo":float(lo),"zhi":float(hi),"z_effective":float(np.mean(np.asarray(a["z_effective"],float)[m]))})
    zmin=min(b["zlo"] for b in bins); zmax=max(b["zhi"] for b in bins)

    zl,wl=read_zw(args.lrg_data,zmin,zmax)
    ze,we=read_zw(args.elg_data,zmin,zmax)

    summary={
        "scope":"DR1-specific evolution-bias benchmark from observed weighted N(z)",
        "definition":"f_evo = -(1+z) d ln(nbar_comoving)/dz",
        "weights":"WEIGHT only; FKP deliberately excluded from selection-function calibration",
        "geometry":"weighted dN/dz divided by fiducial comoving shell volume per steradian; survey area cancels in logarithmic derivative",
        "predeclared_fit_ensemble":"dz in {0.01,0.02} x polynomial degree in {2,3}; benchmark is median, half-range is smoothing sensitivity diagnostic",
        "LRG":tracer_summary(zl,wl,bins),
        "ELG":tracer_summary(ze,we,bins),
        "guardrail":(
            "Observed N(z) contains cosmic variance and residual selection effects. "
            "Use these values as DR1-specific benchmark inputs, not as final calibrated evolution biases. "
            "Final inference should validate against the official DESI selection function / random construction."
        ),
    }
    p=Path(args.out); p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
