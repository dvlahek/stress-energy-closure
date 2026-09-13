#!/usr/bin/env python3
"""Survey-specific hidden-state neutrino-wake Fisher for a published DESI BGS split.

This forecast uses the flux-limited 50/50 bright/faint DESI-BGS mock catalogue
of Bonvin et al., MNRAS 525, 4611 (2023), Table 4.  Their catalogue was built
specifically to study a parity-odd cross-correlation dipole and provides, in
each redshift bin, the BGS-footprint galaxy counts and measured tracer biases.

We keep the kinetic/wake modelling and conservative odd-nuisance projection of
wake_two_tracer_fisher.py, but replace the idealized constant tracer densities,
biases and equal bin volumes by the published survey values.  This therefore
answers a narrower and more defensible question: can two stress-energy-matched
neutrino distributions be distinguished by a parity-odd observable in a tracer
configuration already demonstrated for DESI BGS?

Published flux-limited 50/50 BGS split used here (Bonvin et al. Table 4):
 z=0.25: N_F=1,178,444, N_B=1,177,602, b_F=1.16,  b_B=1.45
 z=0.35: N_F=  474,971, N_B=  472,874, b_F=1.380, b_B=1.89
 z=0.45: N_F=   78,120, N_B=   78,349, b_F=1.38,  b_B=2.10
The bins have width Delta z=0.1 and represent a DESI-BGS-sized footprint.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path
import numpy as np
from classy import Class
import class_response_optimize as cro
import wake_two_tracer_fisher as base

C_KMS = base.C_KMS
KOBS = base.KOBS
MU, MUW = base.MU, base.MUW
R_FILTER = base.R_FILTER

# Bonvin et al. 2023, flux-limited Case 3, Table 4.
SURVEY = [
    {"z":0.25,"zlo":0.20,"zhi":0.30,"N_faint":1178444,"N_bright":1177602,"b_faint":1.16,"b_bright":1.45},
    {"z":0.35,"zlo":0.30,"zhi":0.40,"N_faint":474971,"N_bright":472874,"b_faint":1.380,"b_bright":1.89},
    {"z":0.45,"zlo":0.40,"zhi":0.50,"N_faint":78120,"N_bright":78349,"b_faint":1.38,"b_bright":2.10},
]
AREA_DEG2 = 14000.0
FSKY = AREA_DEG2 / (4.0*np.pi*(180.0/np.pi)**2)


def shell_volume_from_class(cosmo, zlo, zhi):
    """Comoving shell volume in (Gpc/h)^3 over 14,000 deg^2."""
    h=float(cosmo.h())
    # CLASS angular distance: comoving distance = (1+z) d_A in Mpc.
    rlo=(1.0+zlo)*float(cosmo.angular_distance(zlo))*h
    rhi=(1.0+zhi)*float(cosmo.angular_distance(zhi))*h
    return FSKY*(4.0*np.pi/3.0)*(rhi**3-rlo**3)/1.0e9


def mode_noise_weight(info, volume_gpc3, bl, bf, nl, nf):
    pk, fg = info["pk"], info["f"]
    dk=np.gradient(KOBS)
    rows=[]
    for ik,k in enumerate(KOBS):
        for im,mu in enumerate(MU):
            al=bl+fg[ik]*mu**2
            af=bf+fg[ik]*mu**2
            Cll=al**2*pk[ik]+1.0/nl
            Cff=af**2*pk[ik]+1.0/nf
            Cfl=al*af*pk[ik]
            det=max(Cll*Cff-Cfl*Cfl,1e-300)
            pref2=pk[ik]**2*(bl-bf)**2*mu**4/(0.5*det)
            mode=volume_gpc3*1e9/(4.0*np.pi**2)*k**2*dk[ik]*MUW[im]
            rows.append((ik,im,mode*pref2))
    return rows


def desired_fd_sn2(mass,z,Ntot,delta_b):
    # Okoli et al. Eq. 35 literature calibration used in our baseline forecast.
    nth=1.7e7*(mass/0.05)**(-6.0)*(28.5**z)/(delta_b**2)
    return 9.0*Ntot/nth


def prepare_survey(state0):
    result=[]
    for s in SURVEY:
        z=float(s["z"])
        # Volume comes from the same CLASS fiducial used for the wake transfer.
        # build_state has already closed CLASS, so reconstruct shell volume with
        # a short fiducial CLASS instance here.
        result.append(dict(s))
    p=Class()
    # Fiducial geometry only. ncdm details have negligible effect on shell volume,
    # but use the same background parameters as the kinetic run.
    p.set({"H0":cro.H0,"omega_b":cro.OMEGA_B,"omega_cdm":cro.OMEGA_CDM,
           "A_s":cro.A_S,"n_s":base.NS,"tau_reio":cro.TAU_REIO,
           "N_ur":cro.N_UR,"N_ncdm":0})
    p.compute()
    try:
        for r in result:
            V=shell_volume_from_class(p,r["zlo"],r["zhi"])
            r["volume_hminus3_Gpc3"]=float(V)
            r["n_faint_h3_Mpc3"]=float(r["N_faint"]/(V*1e9))
            r["n_bright_h3_Mpc3"]=float(r["N_bright"]/(V*1e9))
            r["delta_b"]=float(r["b_bright"]-r["b_faint"])
    finally:
        p.struct_cleanup(); p.empty()
    return result


def calibrated_linear_fisher(q,f0,shapes,state0,mass,frac,survey,nuisance_level):
    nd=shapes.shape[0]; nz=len(survey)
    nn=nz*(1 if nuisance_level=="amp" else 4)
    F=np.zeros((nd,nd)); G=np.zeros((nd,nn)); N=np.zeros((nn,nn)); perz=[]
    for iz,s in enumerate(survey):
        z=float(s["z"]); info=state0["z"][z]
        bl,bf=float(s["b_bright"]),float(s["b_faint"])
        nl,nf=float(s["n_bright_h3_Mpc3"]),float(s["n_faint_h3_Mpc3"])
        V=float(s["volume_hminus3_Gpc3"])
        rows=mode_noise_weight(info,V,bl,bf,nl,nf)
        cache={}; raw_base=0.0
        for ik,im,w in rows:
            if ik not in cache:
                cache[ik]=base.linear_basis_moments(q,f0,shapes,state0,mass,z,ik,frac)
            m00,m0,mm=cache[ik]; raw_base += w*m00
        target=desired_fd_sn2(mass,z,s["N_faint"]+s["N_bright"],bl-bf)
        cal=target/max(raw_base,1e-300)
        perz.append({"z":z,"volume_hminus3_Gpc3":V,"n_faint":nf,"n_bright":nl,
                     "b_faint":bf,"b_bright":bl,"delta_b":bl-bf,
                     "N_total":int(s["N_faint"]+s["N_bright"]),
                     "FD_SN":math.sqrt(target),"calibration":cal,
                     "sigma_v_R16_kmax_kms":float(info["sigma_kms"][-1])})
        for ik,im,w0 in rows:
            w=w0*cal; m00,m0,mm=cache[ik]
            ns=base.nuisance_shapes(iz,KOBS[ik],MU[im],nz,nuisance_level)
            F += w*mm
            G += w*np.outer(m0,ns)
            N += w*m00*np.outer(ns,ns)
    P=F-G@np.linalg.pinv(N,rcond=1e-11)@G.T
    return P,F,perz


def nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,mass,survey,nuisance_level):
    nz=len(survey); nn=nz*(1 if nuisance_level=="amp" else 4)
    Fdd=0.0; g=np.zeros(nn); N=np.zeros((nn,nn)); raw=0.0; perz=[]
    for iz,s in enumerate(survey):
        z=float(s["z"]); info=state0["z"][z]
        bl,bf=float(s["b_bright"]),float(s["b_faint"])
        nl,nf=float(s["n_bright_h3_Mpc3"]),float(s["n_faint_h3_Mpc3"])
        V=float(s["volume_hminus3_Gpc3"])
        rows=mode_noise_weight(info,V,bl,bf,nl,nf)
        cache={}; raw_base=0.0
        for ik,im,w in rows:
            if ik not in cache:
                cache[ik]=base.stochastic_moments(q,f0,fp,fm,state0,statep,statem,mass,z,ik)
            m00,mdd,m0d=cache[ik]; raw_base += w*m00
        target=desired_fd_sn2(mass,z,s["N_faint"]+s["N_bright"],bl-bf)
        cal=target/max(raw_base,1e-300); zun=0.0
        for ik,im,w0 in rows:
            w=w0*cal; m00,mdd,m0d=cache[ik]
            ns=base.nuisance_shapes(iz,KOBS[ik],MU[im],nz,nuisance_level)
            Fdd += w*mdd; raw += w*mdd; zun += w*mdd
            g += w*m0d*ns
            N += w*m00*np.outer(ns,ns)
        perz.append({"z":z,"unprojected_hidden_SN":math.sqrt(max(zun,0.0))})
    proj=max(Fdd-float(g@np.linalg.pinv(N,rcond=1e-11)@g),0.0)
    return math.sqrt(max(raw,0.0)),math.sqrt(proj),perz


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",default="wake_desi_bonvin_output")
    ap.add_argument("--mass",type=float,default=0.06)
    ap.add_argument("--z-match",type=float,default=1100.0)
    ap.add_argument("--frac",type=float,default=0.30)
    ap.add_argument("--samples",type=int,default=50000)
    ap.add_argument("--validate",type=int,default=5)
    ap.add_argument("--seed",type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    base.ZBINS=np.array([s["z"] for s in SURVEY],float)
    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/"fd_reference.dat"; base.write_psd(p0,q,f0)
    state0=base.build_state(p0,args.mass)
    survey=prepare_survey(state0)

    Bfull,Braw,perz=calibrated_linear_fisher(q,f0,shapes,state0,args.mass,args.frac,survey,"full")
    selected=base.candidate_pool(Bfull,shapes,f0,args.frac,args.samples,args.seed,args.validate)
    validations=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f"cand_{rank:02d}_plus.dat"; pm=out/f"cand_{rank:02d}_minus.dat"
        base.write_psd(pp,q,fp); base.write_psd(pm,q,fm)
        statep=base.build_state(pp,args.mass); statem=base.build_state(pm,args.mass)
        sn_un,sn_full,pz=nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,args.mass,survey,"full")
        _,sn_amp,_=nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,args.mass,survey,"amp")
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={"rank":rank,"predicted_full_projected_SN":math.sqrt(max(pred2,0.0)),
             "validated_unprojected_hidden_SN":sn_un,
             "validated_after_wake_amplitude_projection_SN":sn_amp,
             "validated_after_full_odd_nuisance_projection_SN":sn_full,
             "max_relative_moment_mismatch":float(mismatch.max()),
             "per_redshift":pz,"coefficients":coeff.tolist()}
        validations.append(row); print(json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    if best is None: raise RuntimeError("no validated candidates")
    sn,row,fp,fm,shape=best
    np.savetxt(out/"best_pair.csv",np.column_stack([q,f0,fp,fm,shape]),delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape",comments="")
    summary={"class_commit":cro.CLASS_COMMIT,"mass_eV":args.mass,"z_match":args.z_match,
             "fractional_distortion_cap":args.frac,
             "survey_source":"Bonvin et al., MNRAS 525, 4611 (2023), Table 4, flux-limited 50/50 DESI BGS split",
             "area_deg2":AREA_DEG2,"survey_bins":survey,"FD_calibration":perz,
             "nuisance_model_full":"independent wake amplitude, k^-1 odd template, k odd broadband, and k^2 odd broadband in each redshift bin",
             "linear_best_predicted_SN":math.sqrt(max(selected[0][0],0.0)),
             "nonlinear_candidates_validated":len(validations),"best":row,
             "all_validations":validations,
             "interpretation":"Published DESI-BGS tracer counts and measured mock biases replace all idealized tracer densities/biases. Absolute wake normalization remains literature-calibrated to Okoli et al. Eq.35. This is a survey-specific forecast diagnostic, not a DESI collaboration likelihood."}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("SURVEY\n"+json.dumps(survey,indent=2)); print("BEST\n"+json.dumps(row,indent=2))

if __name__=="__main__": main()
