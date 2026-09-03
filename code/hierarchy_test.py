#!/usr/bin/env python3
# Finite gravitational-source jet hierarchy test.

import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def integ(y,x):
    if hasattr(np,"trapezoid"): return np.trapezoid(y,x)
    return np.trapz(y,x)

def poly_next(c):
    a=np.zeros(len(c)+1)
    a[:len(c)] += -4*c
    a[1:] += c
    if len(c)>1:
        d=np.arange(1,len(c))*c[1:]
        a[1:1+len(d)] += 2*d
        a[2:2+len(d)] += -2*d
    return a

def polynomials(nmax):
    P=[np.array([1.0])]
    for _ in range(nmax): P.append(poly_next(P[-1]))
    return P

def peval(c,y): return np.polynomial.polynomial.polyval(y,c)

def bump(q,c,w):
    x=(q-c)/w
    z=np.zeros_like(q)
    inside=np.abs(x)<1
    z[inside]=np.exp(-1.0/(1.0-x[inside]**2))
    return z

def main():
    out=Path("ev_np_campaign_output"); out.mkdir(exist_ok=True)
    q=np.linspace(0.02,6.0,20000)
    E=np.sqrt(q*q+1); y=1/(q*q+1)
    w0=q*q*E
    Ps=polynomials(9)
    rows=[]; details={}
    for r in range(0,7):
        J=r+2
        centers=np.linspace(0.35,4.8,J)
        spacing=(centers[1]-centers[0]) if J>1 else 0.5
        width=min(0.18,0.35*spacing)
        B=np.array([bump(q,c,width) for c in centers])
        L=np.array([[integ(w0*peval(Ps[n],y)*b,q) for b in B] for n in range(r+2)])
        A=L[:r+1]
        scales=np.maximum(np.linalg.norm(A,axis=1),1e-300)
        _,S,VT=np.linalg.svd(A/scales[:,None])
        c=VT[-1]; c/=np.max(np.abs(c))
        base=np.ones(J)/J
        mix=0.35/J
        cp=base+mix*c; cm=base-mix*c
        vals_p=L@cp; vals_m=L@cm
        dif=vals_p-vals_m
        denom=np.maximum(np.maximum(np.abs(vals_p[:r+1]),np.abs(vals_m[:r+1])),1e-30)
        matched=np.max(np.abs(dif[:r+1])/denom)
        next_rel=abs(dif[r+1])/max(abs(vals_p[r+1]),abs(vals_m[r+1]),1e-30)
        rows.append([r,matched,next_rel,abs(dif[r+1]),S[-1]])
        details[str(r)]={"matched_max_relative":float(matched),
                         "next_relative_difference":float(next_rel),
                         "next_absolute_difference":float(abs(dif[r+1])),
                         "coeff_plus":cp.tolist(),"coeff_minus":cm.tolist()}
        print("r=",r," matched=",matched," next_rel=",next_rel)
    arr=np.array(rows,float)
    np.savetxt(out/"hierarchy_test.csv",arr,delimiter=",",
               header="r,max_relative_matched_jet_error,next_relative_difference,next_abs_difference,min_singular_value",comments="")
    plt.figure(figsize=(7,4.5))
    plt.semilogy(arr[:,0],np.maximum(arr[:,1],1e-18),marker="o",label="matched jets error")
    plt.semilogy(arr[:,0],np.maximum(arr[:,2],1e-18),marker="o",label="next unmatched jet")
    plt.xlabel("highest matched source-jet order r"); plt.ylabel("relative magnitude")
    plt.title("Finite jet matching leaves a next kinetic direction")
    plt.legend(); plt.tight_layout(); plt.savefig(out/"hierarchy_test.png",dpi=180); plt.close()
    with open(out/"hierarchy_summary.json","w") as f:
        json.dump({"polynomials":[p.tolist() for p in Ps],"details":details},f,indent=2)

if __name__=="__main__": main()
