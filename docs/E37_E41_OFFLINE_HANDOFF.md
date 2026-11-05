# EinsteinVlasovNP — E37–E40 offline handoff (30. 9. 2026)

## New results

**E37:** two snapshots become quantitatively useful only after an independent time-derivative bound. For an `L`-Lipschitz history with endpoint slope `d`, the sharp maximum deviation from endpoint linear interpolation is `(L^2-d^2)T/(2L)`.

**E38:** E35 and E37 combine into one sufficient absolute error certificate:
`|A-A_point| <= integral |W|[B_time + delta_test +(1-s)delta0+s delta1]dt`.
Needed inputs are physical second-moment bounds and a physical `L_k` or stronger dynamical envelope.

**E39:** frozen E28R1 low-mass alpha=.4 has signed Fminus-Fplus `-1.259955e-05 (km/s)^2/Mpc`, 0.6242% of FD, while 94.22589% of finite-band FD force is reported from `k>1/Mpc`. A conservative sufficient sign-robustness condition would require state-differential high-k model error below about 0.662% of that high-k FD contribution.

**E40:** with free state-specific tracer response `chi_F(K)`, the E21 linear channel is exactly non-identifying for F. Under the optimistic common-scalar-response restriction, frozen E21 F+/F- four-K source vectors have only 0.0969% amplitude-free norm separation (`theta=9.686e-04 rad`). Frozen E19 conditional T1/T3 source angle is about 24.6x larger, but it is a different conditioned estimand.

## Consequence

The blocker is no longer 'we need more metadata'. A physical forward model needs:
1. a dynamical halo+environment history supplying `L_k` and radial moments through time;
2. an independently specified incoming neutrino phase-space state or bound on its homogeneous Vlasov term;
3. external high-z LRG/ELG response constraints strong enough that `chi_F` cannot absorb the source-state difference;
4. short-scale/Born physical error below the state-differential level relevant to any F+/F- claim.

Until then E27/E28 remain conditional finite-band benchmarks. A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED, observed eBOSS 24D odd SEALED. No GitHub/main/PR mutation is performed by this package.


## E41 addition

A simple symmetric multiplicative response-uncertainty model quantifies the conditioning problem. Frozen E21 F+/F− long-K components require `eta` below 0.0912% to keep **all four** component intervals disjoint. The frozen E19 conditional T3/T1 ratio remains disjoint up to about 2.605% ratio uncertainty, roughly 28.6 times looser by this narrow source-level metric. This is not a covariance or sensitivity forecast; it reinforces that the original unconditioned two-point route needs extremely strong external tracer-response calibration, while a wind-conditioned angular fingerprint is better conditioned but is a different experiment.
