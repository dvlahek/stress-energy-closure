# E46 — can the physical Einstein–Vlasov result be placed on the E45 sensitivity curve?

## Short answer

**Not yet as a number.** E45's amplitude `A` is the coefficient of the archived synthetic pair-level DD modulation `1 + A mu_bin`. It is not the same object as the E28 halo acceleration and it is not automatically the E19/E21 Fourier source amplitude.

The exact comparison is nevertheless now defined.

For each cap/case let `q_c` be the archived 12D unit synthetic response in the six `s` bins for `ell=1,3`. If a physical model produces a **dimensionless, post-window eBOSS odd vector** `y_phys,c` in exactly the same ordering, then

`A_phys,c = (q_c dot y_phys,c)/(q_c dot q_c)`.

That `A_phys` can be placed directly on the E45 practical-sensitivity curve.

## What E45 says

The archived nine-mock baseline scatter is:

- NGC SD = 0.094948401
- SGC SD = 0.064236974.

For symmetric +/- sign recovery on these same nine mocks, the discrete empirical amplitude thresholds are approximately:

- ~83.3% recovery: `A > 0.040599` NGC and `A > 0.046332` SGC.
- ~94.4% recovery: `A > 0.151322` NGC and `A > 0.090913` SGC.
- 9/9 sign recovery for both signs requires `A` above the largest absolute baseline amplitude: 0.171306 NGC and 0.124060 SGC.

These are descriptive thresholds, not detection significances.

## Why E28 cannot simply be divided by 0.02

E28 gives a conditional halo-acceleration contrast

`a_Fminus-a_Fplus = -1.25995524466503e-05 (km/s)^2/Mpc`.

E25 shows that a halo acceleration is not itself the observed galaxy-number-count odd coefficient. In the geodesic full-metric limit, the wake must not be double-counted as an extra Euler residual. The actual local galaxy response needs independently specified `beta`, `epsilon`, selection response and the ordinary number-count coefficient `D_a`.

If one *purely diagnostically* wrote `A = kappa_a * delta_a`, then reaching one baseline SD would require

- `kappa_a ~ 7.54e+03 Mpc/(km/s)^2` in NGC,
- `kappa_a ~ 5.1e+03 Mpc/(km/s)^2` in SGC.

No existing result derives such a `kappa_a`; therefore these are **required gains**, not predictions.

## Why E19/E21 also do not close the bridge

E19 is conditional Fourier source geometry at one `k=0.05 h/Mpc` and `z=.95`. It is explicitly not a physical eBOSS `xi_l(s)`.

E21 supplies a source `C_v(K)` but leaves the high-z LRG/ELG response coefficient unmeasured.

E6 defines the correct mathematical operator once a physical pre-window `xi_L(s,z)` exists, but the current wake builder normalizes the wake shape to unit maximum and explicitly sets `absolute_wake_prediction=False`.

Therefore the next useful physics task is **not another sensitivity test**. It is to build exactly one prospective absolute `y_phys` under a clearly labeled physical response scenario, with no tuning to observed odd data. Only then can E45 tell us if that scenario is practically measurable.
