#!/usr/bin/env python3
# Exact FLRW particle-mass sweep with fresh same-(n,rho,P) pair at each mass.

import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

C=1/(2*np.pi)**3
def integ(y,x):
    if hasattr(np,"trapezoid"): return np.trapezoid(y,x)
    return np.trapz(y,x)
def gauss(q,s): return np.exp(-0.5*(q/s)**2)

def moms(q,F,m,a=1.):
    E=np.sqrt(q*q+(a*m)**2)
    n=C*4*np.pi/a**3*integ(q*q*F,q)
    rho=C*4*np.pi/a**4*integ(q*q*E*F,q)
    P=C*4*np.pi/(3*a**4)*integ(np.divide(q**4,E,out=np.zeros_like(q),where=E>0)*F,q)
    Q=C*4*np.pi*integ(np.divide(q**4,E**3,out=np.zeros_like(q),where=E>0)*F,q)
    return n,rho,P,Q

def pair(q,m):
    sig=np.array([0.35,0.75,1.4,2.6])
    B=np.array([gauss(q,s) for s in sig])
    V=np.array([[moms(q,b,m)[i] for b in B] for i in range(3)])
    scales=np.maximum(np.linalg.norm(V,axis=1),1e-300)
    A=V/scales[:,None]
    _,_,VT=np.linalg.svd(A)
    c=VT[-1]; c/=np.max(np.abs(c))
    base=np.ones(4)/4; mix=0.18
    cp=base+mix*c; cm=base-mix*c
    Fp=cp@B; Fm=cm@B
    rp=moms(q,Fp,m)[1]; rm=moms(q,Fm,m)[1]
    scale=0.5*(rp+rm); Fp/=scale; Fm/=scale
    return Fp,Fm

def evolve(q,F,m,tmax=4.):
    def rhs(t,y):
        a=y[0]; rho=moms(q,F,m,a)[1]
        return [a*np.sqrt(rho)]
    s=solve_ivp(rhs,(0,tmax),[1.0],rtol=2e-10,atol=2e-12,max_step=0.01,dense_output=True)
    t=np.linspace(0,tmax,401)
    return t,s.sol(t)[0]

def main():
    out=Path("ev_np_campaign_output"); out.mkdir(exist_ok=True)
    q=np.linspace(1e-6,15,9000)
    masses=[0.0,0.05,0.1,0.2,0.5,1.0,2.0,5.0]
    rows=[]
    for m in masses:
        Fp,Fm=pair(q,m)
        np_,rp,pp,qp=moms(q,Fp,m); nm,rm,pm,qm=moms(q,Fm,m)
        H=np.sqrt(0.5*(rp+rm))
        dQ=qp-qm
        d3=0.5*m*m*H*dQ
        t,ap=evolve(q,Fp,m); _,am=evolve(q,Fm,m)
        rows.append([m,np_-nm,rp-rm,pp-pm,dQ,d3,np.max(np.abs(ap-am))])
        print("m=",m," delta_n,rho,P,Q,a3,maxda =",rows[-1][1:])
    arr=np.array(rows,float)
    np.savetxt(out/"mass_sweep.csv",arr,delimiter=",",
               header="m,delta_n,delta_rho,delta_P,delta_Q,delta_a3,max_delta_a",comments="")
    plt.figure(figsize=(7,4.5))
    plt.plot(arr[:,0],np.abs(arr[:,5]),marker="o")
    plt.xlabel("particle mass m"); plt.ylabel("|Delta a'''(0)|")
    plt.title("Massless universality to massive kinetic distinguishability")
    plt.tight_layout(); plt.savefig(out/"mass_sweep_delta_a3.png",dpi=180); plt.close()
    plt.figure(figsize=(7,4.5))
    plt.plot(arr[:,0],np.abs(arr[:,6]),marker="o")
    plt.xlabel("particle mass m"); plt.ylabel("max_t |Delta a(t)|")
    plt.title("Exact FLRW geometry separation versus particle mass")
    plt.tight_layout(); plt.savefig(out/"mass_sweep_geometry.png",dpi=180); plt.close()
    with open(out/"mass_sweep_summary.json","w") as f: json.dump({"rows":rows},f,indent=2)

if __name__=="__main__": main()
