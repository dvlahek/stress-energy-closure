#!/usr/bin/env python3
"""Independent two-tracer Fisher sanity check.

We calculate the two-tracer limit independently of the wake forecast modules and verify
from first principles that the covariance-matrix expression used by the
multi-tracer code reduces, for two tracers, to the closed analytic formula

  Tr(C^-1 D C^-1 D) = 2 P^2 mu^4 (b1-b2)^2 / det(C),

with C_ij = P a_i a_j + delta_ij/n_i and
D_12 = i P mu^2 (b1-b2), D_21 = D_12*.

We report the analytic identity and a deterministic numerical comparison to the covariance-matrix expression.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np


def direct_matrix(P, mu, b1, b2, f, n1, n2):
    a1=b1+f*mu*mu; a2=b2+f*mu*mu
    C=np.array([[P*a1*a1+1.0/n1, P*a1*a2],
                [P*a1*a2, P*a2*a2+1.0/n2]],float)
    D=np.array([[0.0,1j*P*mu*mu*(b1-b2)],
                [-1j*P*mu*mu*(b1-b2),0.0]],complex)
    Ci=np.linalg.inv(C)
    return float(np.trace(Ci@D@Ci@D).real), float(np.linalg.det(C))


def closed_form(P, mu, b1, b2, f, n1, n2):
    a1=b1+f*mu*mu; a2=b2+f*mu*mu
    det=(P*a1*a1+1.0/n1)*(P*a2*a2+1.0/n2)-(P*a1*a2)**2
    val=2.0*P*P*mu**4*(b1-b2)**2/det
    return float(val), float(det)


def symbolic_check():
    import sympy as sp
    P,u,b1,b2,f,n1,n2=sp.symbols('P u b1 b2 f n1 n2', positive=True, finite=True)
    a1=b1+f*u**2; a2=b2+f*u**2
    C=sp.Matrix([[P*a1**2+1/n1,P*a1*a2],[P*a1*a2,P*a2**2+1/n2]])
    D=sp.Matrix([[0,sp.I*P*u**2*(b1-b2)],[-sp.I*P*u**2*(b1-b2),0]])
    lhs=sp.factor(sp.trace(C.inv()*D*C.inv()*D))
    rhs=sp.factor(2*P**2*u**4*(b1-b2)**2/C.det())
    return bool(sp.simplify(lhs-rhs)==0), str(lhs), str(rhs)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--outdir',default='fisher_independent_check')
    ap.add_argument('--samples',type=int,default=10000); ap.add_argument('--seed',type=int,default=20260913)
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    symbolic_ok,lhs,rhs=symbolic_check()
    rng=np.random.default_rng(args.seed); rel=[]; abs_err=[]
    worst=None
    for _ in range(args.samples):
        P=10**rng.uniform(2.0,5.0); mu=rng.uniform(0.02,1.0)
        b1=rng.uniform(0.7,4.5); b2=rng.uniform(0.7,4.5)
        if abs(b1-b2)<0.05: b2+=0.1
        f=rng.uniform(0.3,1.1); n1=10**rng.uniform(-5.0,-1.2); n2=10**rng.uniform(-5.0,-1.2)
        a,det1=direct_matrix(P,mu,b1,b2,f,n1,n2); c,det2=closed_form(P,mu,b1,b2,f,n1,n2)
        e=abs(a-c); r=e/max(abs(c),1e-300); rel.append(r); abs_err.append(e)
        if worst is None or r>worst['relative_error']:
            worst={'relative_error':r,'direct':a,'closed':c,'P':P,'mu':mu,'b1':b1,'b2':b2,'f':f,'n1':n1,'n2':n2,'det_direct':det1,'det_closed':det2}
    summary={'symbolic_identity_exact':symbolic_ok,'symbolic_lhs':lhs,'symbolic_rhs':rhs,
             'samples':args.samples,'seed':args.seed,'max_relative_error':float(max(rel)),
             'median_relative_error':float(np.median(rel)),'max_absolute_error':float(max(abs_err)),
             'worst_case':worst}
    if not symbolic_ok or summary['max_relative_error']>1e-9:
        raise RuntimeError(json.dumps(summary,indent=2))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
