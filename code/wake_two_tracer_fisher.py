#!/usr/bin/env python3
"""Low-z two-tracer parity-odd Fisher forecast for hidden neutrino phase-space states.

The calculation follows the covariance structure of Okoli et al. 2017,
Eq. (34)/Appendix A, but compares two stress-energy-matched distributions
F_+ and F_- at fixed mass.  The relative neutrino-CDM velocity is reconstructed
from CLASS density-transfer derivatives, avoiding ambiguity in velocity-transfer
normalization.  Absolute per-redshift normalization is calibrated to the
published Okoli et al. Eq. (35) FD forecast, while the k-, mu-, redshift-,
and hidden-state dependence is recomputed here.

The matched null space is optimized after projection over a standard wake
amplitude and conservative odd broadband templates in every redshift bin.
Top candidates are then validated with full CLASS transfer functions for F_+
and F_- at the final 30% pointwise cap.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path
import numpy as np
from classy import Class
import class_response_optimize as cro

C_KMS = 299792.458
NS = 0.9649
KPIV_MPC = 0.05
TNU0_EV = cro.T_NCDM * cro.TCMB_K * cro.KB_EV_K
ZBINS = np.array([0.05, 0.20, 0.40])
KOBS = np.geomspace(0.003, 0.10, 48)       # h/Mpc
KVEL = np.geomspace(1.0e-4, 0.15, 180)     # h/Mpc
MU, MUW = np.polynomial.legendre.leggauss(16)
BL, BF = 2.0, 1.0
NL = NF = 0.02                              # h^3/Mpc^3
R_FILTER = 16.0                             # Mpc/h
GHX0, GHW0 = np.polynomial.hermite.hermgauss(48)
GHX = np.sqrt(2.0) * GHX0
GHW = GHW0 / np.sqrt(np.pi)


def write_psd(path: Path, q, f):
    if np.min(f) <= 0:
        raise RuntimeError(f"non-positive PSD in {path}: {np.min(f)}")
    np.savetxt(path, np.column_stack([q, f]), fmt="%.14e")


def class_params(psd: Path, mass: float):
    return {
        "output": "mPk,dTk",
        "modes": "s",
        "gauge": "newtonian",
        "H0": cro.H0,
        "omega_b": cro.OMEGA_B,
        "omega_cdm": cro.OMEGA_CDM,
        "A_s": cro.A_S,
        "n_s": NS,
        "tau_reio": cro.TAU_REIO,
        "N_ur": cro.N_UR,
        "N_ncdm": 1,
        "use_ncdm_psd_files": 1,
        "ncdm_psd_filenames": str(psd.resolve()),
        "m_ncdm": mass,
        "T_ncdm": cro.T_NCDM,
        "deg_ncdm": 1.0,
        "P_k_max_h/Mpc": 0.25,
        "z_max_pk": 0.7,
    }


def interp_transfer(cosmo, z, key, kh):
    t = cosmo.get_transfer(z=float(z), output_format="class")
    kin = np.asarray(t["k (h/Mpc)"], float)
    val = np.asarray(t[key], float)
    return np.interp(kh, kin, val)


def delta_rel(cosmo, z, kh):
    return interp_transfer(cosmo, z, "d_ncdm[0]", kh) - interp_transfer(cosmo, z, "d_cdm", kh)


def delta_v_kms(cosmo, z, kh):
    dz = 0.004 * (1.0 + z)
    zm = max(0.0, z - dz)
    zp = z + dz
    if zm == 0.0 and z < dz:
        d0 = delta_rel(cosmo, z, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - d0) / (zp - z)
    else:
        dm = delta_rel(cosmo, zm, kh)
        dp = delta_rel(cosmo, zp, kh)
        d_dz = (dp - dm) / (zp - zm)
    h = float(cosmo.h())
    k_mpc = kh * h
    H_mpc = float(cosmo.Hubble(float(z)))
    primordial = np.sqrt(cro.A_S * (k_mpc / KPIV_MPC) ** (NS - 1.0))
    # a*dot(delta)/k = -H*d(delta)/dz/k because a(1+z)=1.
    return np.abs(H_mpc * d_dz / k_mpc) * primordial * C_KMS


def tophat(x):
    x = np.asarray(x)
    out = np.ones_like(x)
    m = np.abs(x) > 1e-5
    xm = x[m]
    out[m] = 3.0 * (np.sin(xm) - xm * np.cos(xm)) / xm**3
    return out


def cumulative_trapz_logk(y, k):
    lk = np.log(k)
    out = np.zeros_like(y)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(lk))
    return out


def velocity_summary(cosmo, z):
    dv = delta_v_kms(cosmo, z, KVEL)
    dz = 0.012 * (1.0 + z)
    dvm = delta_v_kms(cosmo, max(0.0, z - dz), KVEL)
    dvp = delta_v_kms(cosmo, z + dz, KVEL)
    ddv_dz = (dvp - dvm) / ((z + dz) - max(0.0, z - dz))
    a = 1.0 / (1.0 + z)
    w = tophat(KVEL * R_FILTER)
    sig2_cum = cumulative_trapz_logk(dv**2 * w**2, KVEL)
    combo = 3.0 * a**3 * dv - a**2 * ddv_dz
    b2_cum = cumulative_trapz_logk(combo**2 * w**2, KVEL)
    sigma = np.sqrt(np.maximum(np.interp(KOBS, KVEL, sig2_cum), 0.0))
    # Eq. A12 shape up to constants common to F0,F+,F- and calibrated below.
    bshape = np.sqrt(np.maximum(np.interp(KOBS, KVEL, b2_cum), 0.0)) / np.maximum(KOBS, 1e-30)
    return sigma, bshape


def pk_cb_h(cosmo, z):
    h = float(cosmo.h())
    return np.array([cosmo.pk_cb_lin(float(k*h), float(z)) * h**3 for k in KOBS])


def growth_f(cosmo, z):
    dz = 0.008 * (1.0 + z)
    p0 = pk_cb_h(cosmo, z)
    if z > dz:
        pm, pp = pk_cb_h(cosmo, z-dz), pk_cb_h(cosmo, z+dz)
        dlnp = (np.log(pp)-np.log(pm))/(2*dz)
    else:
        pp = pk_cb_h(cosmo, z+dz)
        dlnp = (np.log(pp)-np.log(p0))/dz
    return -0.5*(1.0+z)*dlnp


def build_state(psd: Path, mass: float):
    c = Class(); c.set(class_params(psd, mass)); c.compute()
    try:
        zs = {}
        for z in ZBINS:
            sigma, bshape = velocity_summary(c, float(z))
            zs[float(z)] = {
                "sigma_kms": sigma,
                "bshape": bshape,
                "pk": pk_cb_h(c, float(z)),
                "f": growth_f(c, float(z)),
                "H": float(c.Hubble(float(z))),
            }
        h = float(c.h())
    finally:
        c.struct_cleanup(); c.empty()
    return {"h": h, "z": zs}


def mode_noise_weight(zinfo, volume_gpc3):
    """Eq. A10 covariance shape integrated over k and LOS mu, up to z-bin calibration."""
    pk, fg = zinfo["pk"], zinfo["f"]
    dk = np.gradient(KOBS)
    rows=[]
    for ik,k in enumerate(KOBS):
        for im,mu in enumerate(MU):
            al = BL + fg[ik]*mu**2
            af = BF + fg[ik]*mu**2
            Cll = al**2*pk[ik] + 1.0/NL
            Cff = af**2*pk[ik] + 1.0/NF
            Cfl = al*af*pk[ik]
            det = max(Cll*Cff-Cfl*Cfl, 1e-300)
            # Eq. A8/A10 signal factor. H and constants cancel in per-bin calibration.
            signal_pref2 = pk[ik]**2 * (BL-BF)**2 * mu**4 / (0.5*det)
            mode = volume_gpc3*1e9/(4.0*np.pi**2) * k**2 * dk[ik] * MUW[im]
            rows.append((ik,im,mode*signal_pref2))
    return rows


def occupancy_values(qgrid, f, sigma_kms, mass, z):
    # GHX is a unit normal draw of the coherent LOS velocity.
    qres = mass * np.abs(GHX) * sigma_kms / (C_KMS*TNU0_EV*(1.0+z))
    return np.interp(qres, qgrid, f, left=f[0], right=f[-1])


def stochastic_moments(qgrid, f0, fp, fm, state0, statep, statem, mass, z, ik):
    s0 = state0["z"][z]["sigma_kms"][ik]
    sp = statep["z"][z]["sigma_kms"][ik]
    sm = statem["z"][z]["sigma_kms"][ik]
    b0 = state0["z"][z]["bshape"][ik]
    bp = statep["z"][z]["bshape"][ik]
    bm = statem["z"][z]["bshape"][ik]
    F0 = occupancy_values(qgrid,f0,s0,mass,z)
    Fp = occupancy_values(qgrid,fp,sp,mass,z)
    Fm = occupancy_values(qgrid,fm,sm,mass,z)
    x=GHX
    y0=b0*x*F0; yp=bp*x*Fp; ym=bm*x*Fm; yd=yp-ym
    return (
        float(np.sum(GHW*y0*y0)),
        float(np.sum(GHW*yd*yd)),
        float(np.sum(GHW*y0*yd)),
    )


def linear_basis_moments(qgrid, f0, shapes, state0, mass, z, ik, frac):
    sigma=state0["z"][z]["sigma_kms"][ik]
    b=state0["z"][z]["bshape"][ik]
    qres=mass*np.abs(GHX)*sigma/(C_KMS*TNU0_EV*(1.0+z))
    F0=np.interp(qres,qgrid,f0,left=f0[0],right=f0[-1])
    S=np.array([np.interp(qres,qgrid,s,left=s[0],right=s[-1]) for s in shapes])
    x=GHX
    base=b*x*F0
    resp=2.0*frac*b*(S*x[None,:])
    m00=float(np.sum(GHW*base*base))
    m0=np.array([np.sum(GHW*base*r) for r in resp])
    mm=(resp*GHW[None,:])@resp.T
    return m00,m0,mm


def nuisance_shapes(z_index, k, mu, n_z, level="full"):
    vals=[]
    for iz in range(n_z):
        active=1.0 if iz==z_index else 0.0
        vals.append(active)  # free standard-wake amplitude in each z bin
        if level=="full":
            vals.append(active*(0.01/k))     # conservative gravitational-redshift-like odd template
            vals.append(active*(k/0.01))    # smooth odd broadband
            vals.append(active*(k/0.01)**2) # smooth curvature
    return np.asarray(vals,float)


def desired_fd_sn2(mass,z,volume_bin):
    ngal=(NL+NF)*volume_bin*1e9
    nth=1.7e7*(mass/0.05)**(-6.0)*(28.5**z)/((BL-BF)**2)
    return 9.0*ngal/nth


def calibrated_linear_fisher(qgrid,f0,shapes,state0,mass,frac,volume_total,nuisance_level):
    nd=shapes.shape[0]
    nz=len(ZBINS)
    nn=nz*(1 if nuisance_level=="amp" else 4)
    F=np.zeros((nd,nd)); G=np.zeros((nd,nn)); N=np.zeros((nn,nn))
    perz=[]
    volume_bin=volume_total/nz
    for iz,zv in enumerate(ZBINS):
        z=float(zv); info=state0["z"][z]
        rows=mode_noise_weight(info,volume_bin)
        cache={}
        raw_base=0.0
        for ik,im,w in rows:
            if ik not in cache:
                cache[ik]=linear_basis_moments(qgrid,f0,shapes,state0,mass,z,ik,frac)
            m00,m0,mm=cache[ik]
            raw_base += w*m00
        target=desired_fd_sn2(mass,z,volume_bin)
        cal=target/max(raw_base,1e-300)
        perz.append({"z":z,"FD_SN":math.sqrt(target),"calibration":cal,
                     "sigma_v_R16_kmax_kms":float(info["sigma_kms"][-1])})
        for ik,im,w0 in rows:
            w=w0*cal; m00,m0,mm=cache[ik]
            ns=nuisance_shapes(iz,KOBS[ik],MU[im],nz,nuisance_level)
            F += w*mm
            G += w*np.outer(m0,ns)
            N += w*m00*np.outer(ns,ns)
    Ninv=np.linalg.pinv(N,rcond=1e-11)
    P=F-G@Ninv@G.T
    return P,F,perz


def candidate_pool(Bproj, shapes, f0, frac, samples, seed, keep):
    nd=shapes.shape[0]
    evals,evecs=np.linalg.eigh(0.5*(Bproj+Bproj.T))
    cand=[evecs[:,j] for j in np.argsort(evals)[::-1]]
    for j in range(nd):
        e=np.zeros(nd); e[j]=1.0; cand.extend([e,-e])
    rng=np.random.default_rng(seed)
    rr=rng.normal(size=(samples,nd)); rr/=np.linalg.norm(rr,axis=1,keepdims=True)
    cand.extend(rr)
    scored=[]
    for c in cand:
        raw=c@shapes
        maxrel=float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
        if not np.isfinite(maxrel) or maxrel<=0: continue
        norm=1.0/maxrel
        score2=float(norm*norm*(c@Bproj@c))
        scored.append((score2,c.copy(),norm,raw*norm))
    scored.sort(key=lambda x:x[0],reverse=True)
    sel=[]
    for item in scored:
        u=item[1]/np.linalg.norm(item[1])
        if all(abs(float(u@(v[1]/np.linalg.norm(v[1]))))<0.9995 for v in sel):
            sel.append(item)
        if len(sel)>=keep: break
    return sel


def nonlinear_fisher(qgrid,f0,fp,fm,state0,statep,statem,mass,volume_total,nuisance_level):
    nz=len(ZBINS); nn=nz*(1 if nuisance_level=="amp" else 4)
    Fdd=0.0; g=np.zeros(nn); N=np.zeros((nn,nn)); raw_unproj=0.0
    volume_bin=volume_total/nz
    perz=[]
    for iz,zv in enumerate(ZBINS):
        z=float(zv); info=state0["z"][z]
        rows=mode_noise_weight(info,volume_bin)
        cache={}; raw_base=0.0
        for ik,im,w in rows:
            if ik not in cache:
                cache[ik]=stochastic_moments(qgrid,f0,fp,fm,state0,statep,statem,mass,z,ik)
            m00,mdd,m0d=cache[ik]
            raw_base += w*m00
        target=desired_fd_sn2(mass,z,volume_bin)
        cal=target/max(raw_base,1e-300)
        z_un=0.0
        for ik,im,w0 in rows:
            w=w0*cal; m00,mdd,m0d=cache[ik]
            ns=nuisance_shapes(iz,KOBS[ik],MU[im],nz,nuisance_level)
            Fdd += w*mdd; raw_unproj += w*mdd; z_un += w*mdd
            g += w*m0d*ns
            N += w*m00*np.outer(ns,ns)
        perz.append({"z":z,"unprojected_hidden_SN":math.sqrt(max(z_un,0.0))})
    proj=max(Fdd-float(g@np.linalg.pinv(N,rcond=1e-11)@g),0.0)
    return math.sqrt(max(raw_unproj,0.0)),math.sqrt(proj),perz


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--outdir",default="wake_two_tracer_output")
    ap.add_argument("--mass",type=float,default=0.06)
    ap.add_argument("--z-match",type=float,default=1100.0)
    ap.add_argument("--frac",type=float,default=0.30)
    ap.add_argument("--volume",type=float,default=1.0,help="total low-z volume in (Gpc/h)^3")
    ap.add_argument("--samples",type=int,default=50000)
    ap.add_argument("--validate",type=int,default=5)
    ap.add_argument("--seed",type=int,default=20260913)
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)

    q=np.linspace(0.0,20.0,4000)
    f0,weights,basis,Nnull,shapes,y,M=cro.kinetic_objects(q,args.mass,args.z_match)
    p0=out/"fd_reference.dat"; write_psd(p0,q,f0)
    state0=build_state(p0,args.mass)

    Bfull,Braw,perz=calibrated_linear_fisher(q,f0,shapes,state0,args.mass,args.frac,args.volume,"full")
    Bamp,_,_=calibrated_linear_fisher(q,f0,shapes,state0,args.mass,args.frac,args.volume,"amp")
    selected=candidate_pool(Bfull,shapes,f0,args.frac,args.samples,args.seed,args.validate)
    validations=[]; best=None
    for rank,(pred2,coeff,norm,shape) in enumerate(selected):
        fp=f0+args.frac*shape; fm=f0-args.frac*shape
        pp=out/f"cand_{rank:02d}_plus.dat"; pm=out/f"cand_{rank:02d}_minus.dat"
        write_psd(pp,q,fp); write_psd(pm,q,fm)
        statep=build_state(pp,args.mass); statem=build_state(pm,args.mass)
        sn_un,sn_full,pz=nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,args.mass,args.volume,"full")
        _,sn_amp,_=nonlinear_fisher(q,f0,fp,fm,state0,statep,statem,args.mass,args.volume,"amp")
        mp,mm=cro.moments(fp,q,weights),cro.moments(fm,q,weights)
        mismatch=np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)),1e-300)
        row={
            "rank":rank,
            "predicted_full_projected_SN":math.sqrt(max(pred2,0.0)),
            "validated_unprojected_hidden_SN":sn_un,
            "validated_after_wake_amplitude_projection_SN":sn_amp,
            "validated_after_full_odd_nuisance_projection_SN":sn_full,
            "max_relative_moment_mismatch":float(mismatch.max()),
            "per_redshift":pz,
            "coefficients":coeff.tolist(),
        }
        validations.append(row); print(json.dumps(row))
        if best is None or sn_full>best[0]: best=(sn_full,row,fp,fm,shape)
    if best is None: raise RuntimeError("no validated candidates")
    sn,row,fp,fm,shape=best
    np.savetxt(out/"best_pair.csv",np.column_stack([q,f0,fp,fm,shape]),delimiter=",",
               header="q,f_FD,f_plus,f_minus,normalized_delta_shape",comments="")
    summary={
        "class_commit":cro.CLASS_COMMIT,
        "mass_eV":args.mass,
        "z_match":args.z_match,
        "mass_over_Tnu_at_match":float(y),
        "fractional_distortion_cap":args.frac,
        "survey":{
            "total_volume_hminus3_Gpc3":args.volume,
            "redshift_bins":ZBINS.tolist(),
            "n_l_h3_Mpc3":NL,"n_f_h3_Mpc3":NF,
            "b_l":BL,"b_f":BF,"delta_b":BL-BF,
            "k_h_Mpc":[float(KOBS.min()),float(KOBS.max())],
            "R_filter_hminus1_Mpc":R_FILTER,
        },
        "FD_calibration":perz,
        "nuisance_model_full":"independent wake amplitude, k^-1 odd template, k odd broadband, and k^2 odd broadband in each redshift bin",
        "linear_best_predicted_SN":math.sqrt(max(selected[0][0],0.0)),
        "nonlinear_candidates_validated":len(validations),
        "best":row,
        "all_validations":validations,
        "interpretation":"Literature-calibrated low-z two-tracer Fisher using Okoli et al. 2017 Eq.34/A10/A13 covariance and Eq.35 absolute FD normalization. CLASS recomputes relative-velocity transfer for F+/F-. Conservative odd templates are projected per redshift bin. This is a forecast diagnostic, not a survey collaboration likelihood.",
    }
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("BEST\n"+json.dumps(row,indent=2))

if __name__=="__main__": main()
