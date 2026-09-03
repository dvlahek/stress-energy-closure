#!/usr/bin/env python3
# Direct phase-space Vlasov transport versus eliminated retarded memory kernel.

import argparse, json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib.pyplot as plt

C=1.0/(2*np.pi)**3

def integ(y,x):
    if hasattr(np,"trapezoid"): return np.trapezoid(y,x,axis=0)
    return np.trapz(y,x,axis=0)

def phi(p,s): return np.exp(-0.5*(p/s)**2)
def dphi(p,s): return -(p/s**2)*phi(p,s)

def build_pair(p,m=1.0):
    sig=np.array([0.3,0.8,2.0])
    E=np.sqrt(p*p+m*m)
    B=np.array([phi(p,s) for s in sig])
    dB=np.array([dphi(p,s) for s in sig])
    rho=np.array([C*4*np.pi*integ(p*p*E*b,p) for b in B])
    P=np.array([C*4*np.pi/3*integ(p**4/E*b,p) for b in B])
    psi=B/rho[:,None]; dpsi=dB/rho[:,None]
    w=P/rho
    c=np.array([w[1]-w[2],w[2]-w[0],w[0]-w[1]])
    c/=np.max(np.abs(c))
    base=np.ones(3)/3
    cp=base+0.30*c; cm=base-0.30*c
    return cp@psi,cm@psi,cp@dpsi,cm@dpsi

def kernel(p,dF,taus,mu,wmu,k,m):
    E=np.sqrt(p*p+m*m); vp=p/E
    radial=p**5/E*dF
    angw=wmu*(1-mu*mu)**2
    out=np.empty(len(taus))
    pref=C*np.pi/4
    for j,tau in enumerate(taus):
        phase=np.cos(k*vp[:,None]*mu[None,:]*tau)
        A=phase@angw
        out[j]=pref*integ(radial*A,p)
    return out

def pi_from_g(p,g,mu,wmu,m):
    E=np.sqrt(p*p+m*m)
    angw=wmu*(1-mu*mu)**2
    angular=np.real(g)@angw
    return C*np.pi/4*integ(p**4/E*angular,p)

def direct_evolve(p,dF,mu,wmu,k,m,gamma,dt,tmax):
    E=np.sqrt(p*p+m*m); vp=p/E
    omega=k*vp[:,None]*mu[None,:]
    src=(p*dF)[:,None]
    n=int(round(tmax/dt))+1
    t=np.arange(n)*dt
    h=np.zeros(n); vel=np.zeros(n); piv=np.zeros(n)
    h[0]=1.0
    g=np.zeros((len(p),len(mu)),dtype=np.complex128)

    def rhs(hh,vv,gg):
        pi=pi_from_g(p,gg,mu,wmu,m)
        return vv, -k*k*hh+gamma*pi, -1j*omega*gg+src*vv, pi

    for i in range(n-1):
        k1=rhs(h[i],vel[i],g)
        k2=rhs(h[i]+0.5*dt*k1[0],vel[i]+0.5*dt*k1[1],g+0.5*dt*k1[2])
        k3=rhs(h[i]+0.5*dt*k2[0],vel[i]+0.5*dt*k2[1],g+0.5*dt*k2[2])
        k4=rhs(h[i]+dt*k3[0],vel[i]+dt*k3[1],g+dt*k3[2])
        h[i+1]=h[i]+dt*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        vel[i+1]=vel[i]+dt*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        g=g+dt*(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        piv[i]=k1[3]
    piv[-1]=pi_from_g(p,g,mu,wmu,m)
    return t,h,vel,piv

def causal_convolution(K,v,dt):
    n=len(v); out=np.zeros(n)
    for i in range(1,n):
        s=0.5*K[i]*v[0]+0.5*K[0]*v[i]
        if i>1: s+=np.dot(K[1:i][::-1],v[1:i])
        out[i]=dt*s
    return out

def run(label,np_,nmu,dt,out):
    m=1.; k=1.; gamma=0.5; tmax=12.
    p=np.linspace(0,10,np_)
    mu,wmu=leggauss(nmu)
    Fp,Fm,dFp,dFm=build_pair(p,m)
    result={}
    for name,dF in [("plus",dFp),("minus",dFm)]:
        t,h,v,pi_direct=direct_evolve(p,dF,mu,wmu,k,m,gamma,dt,tmax)
        K=kernel(p,dF,t,mu,wmu,k,m)
        pi_mem=causal_convolution(K,v,dt)
        sl=slice(1,None)
        rel=np.sqrt(integ((pi_direct[sl]-pi_mem[sl])**2,t[sl])/
                    max(integ(pi_direct[sl]**2,t[sl]),1e-300))
        maxerr=np.max(np.abs(pi_direct-pi_mem))
        result[name]={"relative_L2_pi_error":float(rel),"max_abs_pi_error":float(maxerr)}
        np.savetxt(out/f"direct_memory_{name}_{label}.csv",
                   np.column_stack([t,h,v,pi_direct,pi_mem,pi_direct-pi_mem,K]),
                   delimiter=",",header="t,h,v,pi_direct,pi_memory,dpi,K",comments="")
        plt.figure(figsize=(7.0,4.5))
        plt.plot(t,pi_direct,label="direct Vlasov")
        plt.plot(t,pi_mem,label="memory kernel")
        plt.xlabel("t"); plt.ylabel("TT anisotropic stress")
        plt.title(f"Direct versus memory response ({name}, {label})")
        plt.legend(); plt.tight_layout()
        plt.savefig(out/f"direct_memory_{name}_{label}.png",dpi=180); plt.close()
    result["label"]=label; result["np"]=np_; result["nmu"]=nmu; result["dt"]=dt
    return result

def main():
    ap=argparse.ArgumentParser()
    grp=ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--quick",action="store_true"); grp.add_argument("--full",action="store_true")
    args=ap.parse_args()
    out=Path("ev_np_campaign_output"); out.mkdir(exist_ok=True)
    cfg=[("quick",70,32,0.03)] if args.quick else [
        ("coarse",100,48,0.02),("fine",150,72,0.01)]
    rr=[run(*c,out) for c in cfg]
    for r in rr: print(r)
    with open(out/"direct_vs_memory_summary.json","w") as f: json.dump(rr,f,indent=2)
    print("Direct-vs-memory finished.")

if __name__=="__main__": main()
