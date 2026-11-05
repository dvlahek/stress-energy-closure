# E37 — sharp two-snapshot causal-history bound once a physical time-derivative envelope exists

**Offline analytic result; no new Abacus/CLASS/FITS data.** E34 proved that two endpoint halo profiles do not determine a retarded source history. E37 asks what additional input would make two snapshots quantitatively useful.

Let `f(t)` be any scalar source factor on `[0,T]` with fixed endpoints `f0,f1` and a physically justified Lipschitz bound `|df/dt|<=L`. Write `d=(f1-f0)/T` and the endpoint linear interpolation `ell(t)=f0+d t`. Necessarily `|d|<=L`. Endpoint cones give

`f-ell <= min((L-d)t,(L+d)(T-t))`

and

`ell-f <= min((L+d)t,(L-d)(T-t))`.

Therefore

`sup_t |f(t)-ell(t)| <= (L^2-d^2)T/(2L)`

for `L>0`. The bound is sharp. For any retarded weight `W(t)`,

`|I-I_lin| <= integral |W(t)| B(t) dt <= Bmax integral |W(t)| dt`.

For an evolving Fourier halo profile this may be applied independently at each `k` with `f=u_src(k,t)` only after an external physical argument supplies `L_k=sup|partial_t u_src|`.

Toy QA: `T=1`, `f0=1`, `f1=1.2`, `L=.8`, positive retarded weight `exp[-2(T-t)]`. Sharp global uncertainty = 0.375. All 409 checks passed.

**Project implication:** a spline through two snapshots is not a physical error bound. E37 identifies the missing dynamical quantity but does not determine it, certify E28 high-k physics or alter the sealed eBOSS observable.
