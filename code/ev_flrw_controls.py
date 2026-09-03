#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

C = 1.0/(2.0*np.pi)**3

def integ(y,x):
    if hasattr(np,"trapezoid"):
        return np.trapezoid(y,x)
    return np.trapz(y,x)

def gauss(q,s):
    return np.exp(-0.5*(q/s)**2)

def moments(q,F,m,a):
    E=np.sqrt(q*q+(a*m)**2)
    n=C*4*np.pi/a**3*integ(q*q*F,q)
    rho=C*4*np.pi/a**4*integ(q*q*E*F,q)
    pweight=np.divide(q**4,E,out=np.zeros_like(q),where=E>0)
    P=C*4*np.pi/(3*a**4)*integ(pweight*F,q)
    if m>0:
        qweight=np.divide(q**4,E**3,out=np.zeros_like(q),where=E>0)
        Q=C*4*np.pi*integ(qweight*F,q)
    else:
        Q=0.0
    return n,rho,P,Q

def evolve(q,F,m,a0,tmax,lam,rtol,atol,max_step):
    def rhs(t,y):
        a=float(y[0]); rho=moments(q,F,m,a)[1]
        return [a*np.sqrt(lam*rho)]
    sol=solve_ivp(rhs,(0,tmax),[a0],dense_output=True,rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

def massless_control(q,a0,tmax,lam,rtol,atol,max_step,nout):
    F1=gauss(q,0.45)
    F2=0.65*gauss(q,0.22)+0.35*gauss(q,1.7)
    F1/=moments(q,F1,0.0,a0)[1]
    F2/=moments(q,F2,0.0,a0)[1]
    n1,rho1,P1,_=moments(q,F1,0.0,a0)
    n2,rho2,P2,_=moments(q,F2,0.0,a0)
    s1=evolve(q,F1,0.0,a0,tmax,lam,rtol,atol,max_step)
    s2=evolve(q,F2,0.0,a0,tmax,lam,rtol,atol,max_step)
    t=np.linspace(0,tmax,nout); a1=s1.sol(t)[0]; a2=s2.sol(t)[0]
    res={
        "rho1":rho1,"rho2":rho2,"P1":P1,"P2":P2,"n1":n1,"n2":n2,
        "max_abs_delta_a":float(np.max(np.abs(a1-a2))),
        "max_abs_delta_F":float(np.max(np.abs(F1-F2))),
        "same_rho":bool(abs(rho1-rho2)<2e-9),
        "radiation_eos_1":bool(abs(P1-rho1/3)<2e-9),
        "radiation_eos_2":bool(abs(P2-rho2/3)<2e-9),
        "different_profiles":bool(np.max(np.abs(F1-F2))>1e-3),
        "identical_geometry":bool(np.max(np.abs(a1-a2))<2e-9),
    }
    return res,t,a1,a2

def massive_strong(q,a0,tmax,lam,rtol,atol,max_step,nout):
    m=1.0
    sigmas=np.array([0.56948686,1.02498214,1.39514640,1.53581232])
    phis=np.array([gauss(q,s) for s in sigmas])
    V=np.array([moments(q,phi,m,a0) for phi in phis]).T
    A=V[:3,:]
    U,S,VT=np.linalg.svd(A)
    c=VT[-1,:]; c/=np.max(np.abs(c))
    base=np.ones(4)/4; mix=0.24
    cp=base+mix*c; cm=base-mix*c
    if min(cp.min(),cm.min())<=0: raise RuntimeError("positivity failed")
    Fp=cp@phis; Fm=cm@phis
    r0=0.5*(moments(q,Fp,m,a0)[1]+moments(q,Fm,m,a0)[1])
    Fp/=r0; Fm/=r0
    np_,rp,pp,qp=moments(q,Fp,m,a0)
    nm,rm,pm,qm=moments(q,Fm,m,a0)
    Hp=np.sqrt(lam*rp); Hm=np.sqrt(lam*rm)
    dPp=-4*pp/a0-m*m*qp/(3*a0**3)
    dPm=-4*pm/a0-m*m*qm/(3*a0**3)
    addot_p=-(lam/2)*a0*(rp+3*pp); addot_m=-(lam/2)*a0*(rm+3*pm)
    rhodot_p=-3*Hp*(rp+pp); rhodot_m=-3*Hm*(rm+pm)
    Pdot_p=dPp*a0*Hp; Pdot_m=dPm*a0*Hm
    adddot_p=-(lam/2)*(a0*Hp*(rp+3*pp)+a0*(rhodot_p+3*Pdot_p))
    adddot_m=-(lam/2)*(a0*Hm*(rm+3*pm)+a0*(rhodot_m+3*Pdot_m))
    predicted=(lam/2)*(a0*Hp)*(m*m/a0**2)*(qp-qm)
    sp=evolve(q,Fp,m,a0,tmax,lam,rtol,atol,max_step)
    sm=evolve(q,Fm,m,a0,tmax,lam,rtol,atol,max_step)
    t=np.linspace(0,tmax,nout); ap=sp.sol(t)[0]; am=sm.sol(t)[0]
    res={
        "sigmas":sigmas.tolist(),"singular_values_constraints":S.tolist(),"null_direction":c.tolist(),
        "coeff_plus":cp.tolist(),"coeff_minus":cm.tolist(),
        "n_plus":np_,"n_minus":nm,"rho_plus":rp,"rho_minus":rm,"P_plus":pp,"P_minus":pm,
        "Q_plus":qp,"Q_minus":qm,"delta_n":np_-nm,"delta_rho":rp-rm,"delta_P":pp-pm,"delta_Q":qp-qm,
        "delta_addot":addot_p-addot_m,"delta_adddot":adddot_p-adddot_m,
        "delta_adddot_predicted":predicted,"adddot_identity_error":abs((adddot_p-adddot_m)-predicted),
        "max_abs_delta_a":float(np.max(np.abs(ap-am))),"final_delta_a":float(ap[-1]-am[-1]),
        "positive_coefficients":bool(min(cp.min(),cm.min())>0),
        "same_n":bool(abs(np_-nm)<2e-9),"same_rho":bool(abs(rp-rm)<2e-9),"same_P":bool(abs(pp-pm)<2e-9),
        "same_initial_acceleration":bool(abs(addot_p-addot_m)<2e-9),"different_Q":bool(abs(qp-qm)>1e-5),
        "different_third_derivative":bool(abs(adddot_p-adddot_m)>1e-5),
        "third_derivative_identity":bool(abs((adddot_p-adddot_m)-predicted)<2e-9),
        "different_geometry":bool(np.max(np.abs(ap-am))>1e-6),
    }
    return res,t,ap,am

def main():
    ap=argparse.ArgumentParser(); mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--quick",action="store_true"); mode.add_argument("--full",action="store_true")
    ap.add_argument("--out",default="ev_flrw_controls_output"); args=ap.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    a0=1.0; tmax=4.0; lam=1.0
    configs=[(7000,0.015,"quick")] if args.quick else [(9000,0.008,"full_coarse"),(16000,0.004,"full_fine")]
    allres=[]
    for nq,max_step,label in configs:
        q=np.linspace(0,15,nq)
        ml,t,a1,a2=massless_control(q,a0,tmax,lam,2e-10,2e-12,max_step,801)
        ms,t,ap_,am_=massive_strong(q,a0,tmax,lam,2e-10,2e-12,max_step,801)
        ml_pass=all(v for k,v in ml.items() if isinstance(v,(bool,np.bool_)))
        ms_pass=all(v for k,v in ms.items() if isinstance(v,(bool,np.bool_)))
        result={"label":label,"nq":nq,"max_step":max_step,"massless":ml,"massive_same_N_and_T":ms,
                "massless_PASS":bool(ml_pass),"massive_PASS":bool(ms_pass),"PASS":bool(ml_pass and ms_pass)}
        allres.append(result)
        np.savetxt(out/f"massless_timeseries_{label}.csv",np.column_stack([t,a1,a2,a1-a2]),delimiter=",",header="t,a1,a2,delta_a",comments="")
        np.savetxt(out/f"massive_timeseries_{label}.csv",np.column_stack([t,ap_,am_,ap_-am_]),delimiter=",",header="t,a_plus,a_minus,delta_a",comments="")
        plt.figure(figsize=(7.2,4.6)); plt.plot(t,a1-a2); plt.xlabel("t"); plt.ylabel("delta a"); plt.title("Massless control: identical FLRW geometry"); plt.tight_layout(); plt.savefig(out/f"massless_delta_a_{label}.png",dpi=180); plt.close()
        plt.figure(figsize=(7.2,4.6)); plt.plot(t,ap_-am_); plt.xlabel("t"); plt.ylabel("delta a"); plt.title("Massive: same N and T, different future geometry"); plt.tight_layout(); plt.savefig(out/f"massive_same_NT_delta_a_{label}.png",dpi=180); plt.close()
        print("="*78); print(label); print("MASSLESS max |a1-a2| =",ml["max_abs_delta_a"],"PASS =",ml_pass)
        print("MASSIVE |dn| =",abs(ms["delta_n"]),"|drho| =",abs(ms["delta_rho"]),"|dP| =",abs(ms["delta_P"]))
        print("MASSIVE |dQ| =",abs(ms["delta_Q"]),"|da'''| =",abs(ms["delta_adddot"]),"max |da| =",ms["max_abs_delta_a"],"PASS =",ms_pass)
        print("OVERALL PASS =",result["PASS"])
    summary={"results":allres}
    if len(allres)==2:
        c,f=allres; summary["convergence"]={
            "massless_max_delta_a_change":abs(c["massless"]["max_abs_delta_a"]-f["massless"]["max_abs_delta_a"]),
            "massive_delta_Q_change":abs(c["massive_same_N_and_T"]["delta_Q"]-f["massive_same_N_and_T"]["delta_Q"]),
            "massive_delta_adddot_change":abs(c["massive_same_N_and_T"]["delta_adddot"]-f["massive_same_N_and_T"]["delta_adddot"]),
            "massive_max_delta_a_change":abs(c["massive_same_N_and_T"]["max_abs_delta_a"]-f["massive_same_N_and_T"]["max_abs_delta_a"])}
    with open(out/"summary.json","w",encoding="utf-8") as fh: json.dump(summary,fh,indent=2)
    with open(out/"certificate.txt","w",encoding="utf-8") as fh:
        for r in allres:
            fh.write(f"{r['label']}: {'PASS' if r['PASS'] else 'FAIL'}\n")
            fh.write(f"  massless_PASS: {r['massless_PASS']}\n  massive_same_N_and_T_PASS: {r['massive_PASS']}\n")
            fh.write(f"  massless_max_delta_a: {r['massless']['max_abs_delta_a']:.12e}\n")
            m=r['massive_same_N_and_T']
            for key in ['delta_n','delta_rho','delta_P','delta_Q','delta_adddot','max_abs_delta_a']:
                fh.write(f"  massive_{key}: {m[key]:.12e}\n")
            fh.write("\n")
    print("Output:",out.resolve()); print("Zip the whole output folder and send it back.")

if __name__=="__main__": main()
