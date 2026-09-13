#!/usr/bin/env python3
"""DESI-scale full-shape + bispectrum diagnostic for matched kinetic states.

This is a conservative Fisher diagnostic, not the DESI likelihood. It combines
linear Kaiser P_s(k,mu) with a tree-level real-space galaxy bispectrum over a
set of triangle shapes. The matched-state response is projected over standard
cosmological parameters and, independently in each redshift bin, b1, b2, bs2,
EFT-like k^2 power-spectrum shapes and stochastic P/B amplitudes. Gaussian
power and bispectrum covariances are used and P-B cross-covariance is neglected.
Top null-space directions are then validated with full CLASS evaluations at the
final 30 percent pointwise cap.
"""
from __future__ import annotations
import argparse, json, itertools, sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import class_response_optimize as cro
try:
    from classy import Class
except Exception as exc:
    raise RuntimeError("classy required") from exc

NS=0.9649
ZB=np.array([0.5,0.8,1.1,1.4])
BIAS=1.0+0.84*ZB
NBAR=np.array([3.0e-4,3.0e-4,2.5e-4,2.0e-4])
KP=np.geomspace(0.02,0.20,28)
MU,MUW=np.polynomial.legendre.leggauss(10)
KTRI=np.array([0.03,0.05,0.08,0.11,0.15,0.19])
DKTRI=0.015


def triangles():
    out=[]
    for a,b,c in itertools.combinations_with_replacement(KTRI,3):
        k1,k2,k3=sorted((a,b,c),reverse=True)
        if k2+k3>k1+1e-12: out.append((k1,k2,k3))
    return out
TRI=triangles()


def write_psd(path,q,f):
    if np.min(f)<=0: raise RuntimeError("nonpositive PSD")
    np.savetxt(path,np.column_stack([q,f]),fmt="%.14e")


def params(psd,mass,over=None):
    p={"output":"mPk","modes":"s","H0":cro.H0,"omega_b":cro.OMEGA_B,"omega_cdm":cro.OMEGA_CDM,
       "A_s":cro.A_S,"n_s":NS,"tau_reio":cro.TAU_REIO,"N_ur":cro.N_UR,"N_ncdm":1,
       "use_ncdm_psd_files":1,"ncdm_psd_filenames":str(Path(psd).resolve()),"m_ncdm":mass,
       "T_ncdm":cro.T_NCDM,"deg_ncdm":1.0,"P_k_max_h/Mpc":0.6,"z_max_pk":1.6}
    if over: p.update(over)
    return p


def pkgrid(c,z,k):
    h=float(c.h()); return np.array([c.pk_cb_lin(float(x*h),float(z)) for x in k])

def growth(c,z,k):
    p0=pkgrid(c,z,k); dz=0.01*(1+z)
    pm=pkgrid(c,max(0,z-dz),k); pp=pkgrid(c,z+dz,k)
    if z>dz: d=(np.log(pp)-np.log(pm))/(2*dz)
    else: d=(np.log(pp)-np.log(p0))/dz
    return p0,-0.5*(1+z)*d

def F2(k1,k2,k3):
    mu=(k3*k3-k1*k1-k2*k2)/(2*k1*k2)
    return 5/7+0.5*mu*(k1/k2+k2/k1)+2/7*mu*mu

def S2(k1,k2,k3):
    mu=(k3*k3-k1*k1-k2*k2)/(2*k1*k2)
    return mu*mu-1/3

def bg_tree(k1,k2,k3,p1,p2,p3,b1,b2,bs2):
    def term(ka,kb,kc,pa,pb):
        z2=b1*F2(ka,kb,kc)+0.5*b2+0.5*bs2*S2(ka,kb,kc)
        return 2*b1*b1*z2*pa*pb
    return term(k1,k2,k3,p1,p2)+term(k2,k3,k1,p2,p3)+term(k3,k1,k2,p3,p1)

def symfac(k1,k2,k3):
    if np.isclose(k1,k2) and np.isclose(k2,k3): return 6.0
    if np.isclose(k1,k2) or np.isclose(k2,k3) or np.isclose(k1,k3): return 2.0
    return 1.0


def blocks_from_cosmo(c,bias_override=None):
    vec=[]; weights=[]; meta=[]
    vbin=50.0e9/len(ZB)
    dkp=np.gradient(KP)
    for iz,(z,b0,nbar) in enumerate(zip(ZB,BIAS,NBAR)):
        b1,b2,bs2=(float(b0),0.3*(b0-1),-4/7*(b0-1))
        if bias_override and iz in bias_override:
            bo=bias_override[iz]; b1=bo.get("b1",b1); b2=bo.get("b2",b2); bs2=bo.get("bs2",bs2)
        p,f=growth(c,float(z),KP)
        ps=p[:,None]*(b1+f[:,None]*MU[None,:]**2)**2
        ptot=ps+1/nbar
        wp=(vbin/(8*np.pi**2))*(KP[:,None]**2)*dkp[:,None]*MUW[None,:]/np.maximum(ptot**2,1e-300)
        vec.extend(ps.ravel()); weights.extend(wp.ravel()); meta.extend([("P",iz)]*ps.size)
        # Bispectrum monopole proxy in real space. P at triangle sides is interpolated from cb spectrum.
        for k1,k2,k3 in TRI:
            pp=np.interp([k1,k2,k3],KP,p)
            B=bg_tree(k1,k2,k3,*pp,b1,b2,bs2)
            pg=(b1*b1*pp)+1/nbar
            VB=8*np.pi**2*k1*k2*k3*(DKTRI**3)
            var=symfac(k1,k2,k3)*(2*np.pi)**3*np.prod(pg)/(vbin*VB)
            vec.append(B); weights.append(1/max(var,1e-300)); meta.append(("B",iz))
    return np.asarray(vec),np.asarray(weights),meta


def observable(psd,mass,over=None,bias_override=None):
    c=Class(); c.set(params(psd,mass,over)); c.compute()
    try: return blocks_from_cosmo(c,bias_override)
    finally: c.struct_cleanup(); c.empty()

def projector(N,w):
    A=np.sqrt(w)[:,None]*N; q,r=np.linalg.qr(A,mode="reduced")
    d=np.abs(np.diag(r)); keep=d>1e-10*max(1,float(d.max()))
    return q[:,keep]

def sn(delta,w,Q):
    x=np.sqrt(w)*delta; f=float(x@x); y=x-Q@(Q.T@x); p=float(y@y)
    return np.sqrt(max(f,0)),np.sqrt(max(p,0)),p/f if f>0 else 0


def nuisance(psd,mass,ref,weights,meta):
    cols=[]; names=[]
    specs=[("H0","H0",0.25,cro.H0),("omega_b","omega_b",1e-4,cro.OMEGA_B),("omega_cdm","omega_cdm",4e-4,cro.OMEGA_CDM),("lnA_s","A_s",0.01,cro.A_S),("n_s","n_s",0.003,NS)]
    for lab,key,step,cen in specs:
        if lab=="lnA_s": op=observable(psd,mass,{key:cen*np.exp(step)})[0]; om=observable(psd,mass,{key:cen*np.exp(-step)})[0]
        else: op=observable(psd,mass,{key:cen+step})[0]; om=observable(psd,mass,{key:cen-step})[0]
        cols.append((op-om)/(2*step)); names.append(lab)
    # Per-bin galaxy bias derivatives.
    for iz,z in enumerate(ZB):
        baseb=float(BIAS[iz]); vals={"b1":baseb,"b2":0.3*(baseb-1),"bs2":-4/7*(baseb-1)}
        for par,step in [("b1",0.01),("b2",0.01),("bs2",0.01)]:
            bp={iz:dict(vals)}; bm={iz:dict(vals)}; bp[iz][par]+=step; bm[iz][par]-=step
            op=observable(psd,mass,bias_override=bp)[0]; om=observable(psd,mass,bias_override=bm)[0]
            cols.append((op-om)/(2*step)); names.append(f"{par}_z{z:g}")
        # Conservative shape/stochastic directions constructed directly in data space.
        for tag in ["eft0","eft2","Pshot","Bshot0","BshotP"]:
            d=np.zeros_like(ref); ip=0
            # reconstruct index by metadata and local sequence
            pcount=0; bcount=0
            for j,(kind,jz) in enumerate(meta):
                if jz!=iz: continue
                if kind=="P":
                    kk=KP[pcount//len(MU)]; mu=MU[pcount%len(MU)]
                    if tag=="eft0": d[j]=kk*kk*ref[j]
                    elif tag=="eft2": d[j]=kk*kk*mu*mu*ref[j]
                    elif tag=="Pshot": d[j]=1.0
                    pcount+=1
                else:
                    k1,k2,k3=TRI[bcount]
                    if tag=="Bshot0": d[j]=1.0
                    elif tag=="BshotP": d[j]=(k1+k2+k3) # generic smooth stochastic shape
                    bcount+=1
            if np.any(d): cols.append(d); names.append(f"{tag}_z{z:g}")
    return np.column_stack(cols),names


def candidate_pool(B,shapes,f0,frac,samples,seed,keep):
    ev,ec=np.linalg.eigh(B); cs=[ec[:,j] for j in np.argsort(ev)[::-1]]
    nd=shapes.shape[0]
    for j in range(nd): e=np.zeros(nd); e[j]=1; cs.extend([e,-e])
    rng=np.random.default_rng(seed); rr=rng.normal(size=(samples,nd)); rr/=np.linalg.norm(rr,axis=1,keepdims=True); cs.extend(rr)
    scored=[]
    for c in cs:
        raw=c@shapes; mr=float(np.max(np.abs(raw)/np.maximum(f0,1e-300)))
        if not np.isfinite(mr) or mr<=0: continue
        n=1/mr; amp=2*frac*n; s2=float(amp*amp*(c@B@c)); scored.append((s2,c.copy(),raw*n))
    scored.sort(key=lambda x:x[0],reverse=True); out=[]
    for it in scored:
        c=it[1]/np.linalg.norm(it[1])
        if all(abs(float(c@(q[1]/np.linalg.norm(q[1]))))<0.9995 for q in out): out.append(it)
        if len(out)>=keep: break
    return out


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--outdir",default="fullshape_bispectrum_output"); ap.add_argument("--mass",type=float,default=0.60); ap.add_argument("--z-match",type=float,default=1100); ap.add_argument("--probe-frac",type=float,default=.01); ap.add_argument("--final-frac",type=float,default=.30); ap.add_argument("--samples",type=int,default=12000); ap.add_argument("--keep",type=int,default=14); ap.add_argument("--seed",type=int,default=20260913); a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    q=np.linspace(0,20,4000); f0,wm,basis,N,shapes,y,M=cro.kinetic_objects(q,a.mass,a.z_match); f0p=out/"fd.dat"; write_psd(f0p,q,f0)
    ref,w,meta=observable(f0p,a.mass); nuis,names=nuisance(f0p,a.mass,ref,w,meta); Q=projector(nuis,w)
    R=[]
    for j,s in enumerate(shapes):
        fp=f0+a.probe_frac*s; fm=f0-a.probe_frac*s; pp=out/f"probe{j}_p.dat"; pm=out/f"probe{j}_m.dat"; write_psd(pp,q,fp); write_psd(pm,q,fm)
        op=observable(pp,a.mass)[0]; om=observable(pm,a.mass)[0]; R.append((op-om)/(2*a.probe_frac))
    R=np.column_stack(R); RW=np.sqrt(w)[:,None]*R; RP=RW-Q@(Q.T@RW); B=RP.T@RP
    cand=candidate_pool(B,shapes,f0,a.final_frac,a.samples,a.seed,a.keep); rows=[]; best=None
    maskP=np.array([x[0]=="P" for x in meta]); maskB=~maskP
    for rank,(pred2,c,shape) in enumerate(cand):
        fp=f0+a.final_frac*shape; fm=f0-a.final_frac*shape; pp=out/f"cand{rank}_p.dat"; pm=out/f"cand{rank}_m.dat"; write_psd(pp,q,fp); write_psd(pm,q,fm)
        op=observable(pp,a.mass)[0]; om=observable(pm,a.mass)[0]; d=op-om; fixed,proj,ret=sn(d,w,Q)
        # component S/N after projecting the same nuisance columns restricted to each block
        NP=nuis[maskP]; NB=nuis[maskB]; QP=projector(NP,w[maskP]); QB=projector(NB,w[maskB]); fpS,ppS,rpS=sn(d[maskP],w[maskP],QP); fbS,pbS,rbS=sn(d[maskB],w[maskB],QB)
        row={"rank":rank,"predicted_projected_SN":float(np.sqrt(max(pred2,0))),"validated_fixed_SN":fixed,"validated_projected_SN":proj,"retained_fraction":ret,"P_projected_SN":ppS,"B_projected_SN":pbS}
        rows.append(row); print(json.dumps(row))
        if best is None or proj>best[0]: best=(proj,row,shape,fp,fm)
    proj,brow,shape,fp,fm=best; mp,mm=cro.moments(fp,q,wm),cro.moments(fm,q,wm); mis=np.abs(mp-mm)/np.maximum(.5*(np.abs(mp)+np.abs(mm)),1e-300)
    np.savetxt(out/"best_pair.csv",np.column_stack([q,f0,fp,fm,shape]),delimiter=",",header="q,f0,fplus,fminus,shape",comments="")
    summary={"class_commit":cro.CLASS_COMMIT,"mass_eV":a.mass,"z_match":a.z_match,"redshift_bins":ZB.tolist(),"k_power_range":[float(KP.min()),float(KP.max())],"triangle_count":len(TRI),"total_effective_volume_hminus3_Gpc3":50.0,"nuisance_parameters":names,"final_fractional_distortion_cap":a.final_frac,"max_relative_moment_mismatch":float(mis.max()),"validated_candidates":len(rows),"best":brow,"all_candidates":rows,"interpretation":"DESI-scale Gaussian P+B diagnostic with tree-level galaxy bispectrum and conservative nuisance projection. Not the DESI likelihood; P-B cross-covariance and nonlinear EFT loop terms are omitted."}
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n"); print("BEST",json.dumps(brow,indent=2))
if __name__=="__main__": main()
