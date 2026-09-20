#!/usr/bin/env python3
"""Compare legacy Astropy and DESI fiducial redshift-to-distance mappings."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import desi_dr1_lrg_elg_exact as exact

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--zmin",type=float,default=0.8)
    ap.add_argument("--zmax",type=float,default=1.1)
    ap.add_argument("--nz",type=int,default=61)
    ap.add_argument("--out",default="desi_distance_cosmology_audit.json")
    args=ap.parse_args()
    z=np.linspace(args.zmin,args.zmax,args.nz)
    exact.set_distance_cosmology("legacy_astropy")
    old=exact.distance_mpc_over_h(z)
    exact.set_distance_cosmology("desi")
    new=exact.distance_mpc_over_h(z)
    rel=(new-old)/old
    stride=max(1,len(z)//6)
    result={
        "scope":"Redshift-to-distance audit for exact DESI LRGxELG estimator",
        "z_range":[float(args.zmin),float(args.zmax)],
        "n_grid":int(args.nz),
        "legacy":"Astropy FlatLambdaCDM(H0=67.4, Om0=0.315, Tcmb0=2.7255) times h=0.674",
        "production":"cosmoprimo.fiducial.DESI().comoving_radial_distance(z), Mpc/h",
        "max_abs_distance_shift_Mpc_over_h":float(np.max(np.abs(new-old))),
        "max_abs_fractional_shift":float(np.max(np.abs(rel))),
        "rms_fractional_shift":float(np.sqrt(np.mean(rel**2))),
        "samples":[
            {"z":float(zz),"legacy_Mpc_over_h":float(oo),"desi_Mpc_over_h":float(nn),"fractional_shift":float(rr)}
            for zz,oo,nn,rr in zip(z[::stride],old[::stride],new[::stride],rel[::stride])
        ],
        "guardrail":"Even a small mapping change moves pairs across fixed s bins. Final data and mock covariance must therefore use the same production mapping."
    }
    Path(args.out).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
