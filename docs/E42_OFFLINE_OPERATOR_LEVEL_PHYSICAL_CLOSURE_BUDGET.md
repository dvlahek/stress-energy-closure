# E42 — operator-level physical closure budget

E42 combines E35-E41 into one fail-closed criterion. It does not add a new halo model or survey cut.

For each neutrino state F in {plus, minus}, write:

S_F = S_hat_F + dS_init,F + dS_hist,F + dS_prof,F + dS_Born,F

At galaxy level:

Y_F = (chi_hat_F + dchi_F) * S_F

If B_S,F = B_init,F + B_hist,F + B_prof,F + B_Born,F and |dchi_F| <= B_chi,F, then

|Y_F-Y_hat_F| <= |chi_hat_F| B_S,F + |S_hat_F| B_chi,F + B_chi,F B_S,F.

For the contrast D = Y_minus - Y_plus,

|D-D_hat| <= B_Y,minus + B_Y,plus.

Therefore B_Y,minus + B_Y,plus < |D_hat| is a sufficient sign-preservation certificate.

## Frozen E28 numerical gate

For low-mass alpha=.4 E28R1:

- a_Fminus-a_Fplus = -1.25995524466503e-05 (km/s)^2/Mpc
- |delta|/|a_FD| = 0.00624200971 = 0.6242 percent
- 94.22589 percent of the registered FD finite-band force is reported from k>1/Mpc.

Thus at direct halo-acceleration level, the sum of all state-differential physical closure errors must be below 1.25995524466503e-05 (km/s)^2/Mpc to guarantee the frozen contrast sign.

If all unknown differential error were concentrated in the reported high-k contribution, the sufficient limit is 0.662 percent of the high-k FD contribution.

## What each term requires

Incoming kinetic state: E27 set it to zero. E36 showed that matching low moments does not identify hidden phase-space structure. A usable bound must control the propagated homogeneous Vlasov term.

Halo history: E37 supplies the endpoint+Lipschitz formula once a physical L_k = sup|partial_t u(k,t)| is known. Two snapshots alone do not supply L_k.

Finite halo size/profile: E35/E38 supply a conservative second-moment bound. Physical radial moments through the relevant history are still missing.

Born truncation: if a perturbative step operator obeys ||f_(n+1)|| <= q ||f_n||, q<1, the omitted remainder beyond first order is <= q/(1-q) times the first-order norm. We have no physical q bound. If illustratively all of the E28 relative F-state budget were assigned to this one term, q would need to be below 0.00620329. This is not a measured Born parameter.

Tracer response: this is a separate structural gate. E40 proves that unconstrained state-dependent chi_F(K) can absorb the source-state difference. E41 shows that even the optimistic common-amplitude E21 source geometry would need response control below about 0.0912 percent to keep all four long-K component intervals disjoint in that simple source-level model.

## Decision

1. close or bound halo-acceleration source terms: init + history + profile + Born;
2. independently calibrate high-z LRG/ELG response chi_F;
3. only then propagate the physical template through A03;
4. only after that build A04 covariance and consider unblinding.

At present the physical closure inequality cannot be certified. That is a scientific result about identifiability, not a numerical failure.
