# E43 — closure feasibility matrix: what can actually be certified now?

E42 passed locally in WSL. E43 uses only that uploaded result and the already established E35–E41 logic.

The frozen low-mass, alpha=.4 E28R1 Fminus-Fplus acceleration contrast is

`1.25995524466503e-05 (km/s)^2/Mpc`

in absolute value, equal to **0.6242% of |a_FD|**.

That number is the entire available sufficient sign-preservation budget. Any physical closure error budget must fit below it.

## What is currently bounded?

None of the five required physical blocks is yet numerically certified strongly enough:

- **Incoming kinetic state:** unbounded. E27's zero incoming wake is a benchmark, not a reconstructed initial condition.
- **Halo history:** E37 supplies a theorem, but no physical `L_k` is known.
- **Finite-size/profile:** E35/E38 supply a theorem, but no physical radial-moment history is known.
- **Born remainder:** no physical operator-norm or second-order force bound exists yet.
- **Tracer response chi:** E40 shows exact reparameterization freedom when chi is unconstrained; E41 shows the source geometry is poorly conditioned.

Therefore the current physical E28 F-state sign is **not certifiable**, despite the successful numerical quadrature.

## How small is the available budget?

For orientation only, if the four halo-level source uncertainties (incoming, history, profile, Born) were assigned equal shares, each *pairwise Fminus+Fplus category* would get only

`3.149888e-06 (km/s)^2/Mpc`

or **0.1561% of |a_FD|**.

If one divided equally across four categories and both states, each state/category term would get only

`1.574944e-06 (km/s)^2/Mpc`

or **0.0780% of |a_FD|**.

These are not physical priors. They show how tight the combined closure problem is.

## What is the first WSL task worth doing?

Not another broad ASDF column read.

A useful next local-data gate should target a **small same-object multi-epoch dynamical summary** with:

1. documented linkage/tree identity,
2. a consistent mass definition across epochs,
3. radial second moments or a sufficiently rich profile summary,
4. centre and bulk velocity,
5. preferably at least three epochs.

That is the minimum kind of information that could begin to bound `L_k` or the profile-history term. Two endpoints alone cannot do so.

Until such a bound exists, more randoms, more covariance work or observed-odd unblinding does not solve the physical closure problem.
