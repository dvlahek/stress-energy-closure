# E39 — E28 numerical convergence is much stronger than its physical short-scale certification

This uses only frozen E28R1 values; no new force integral was run.

For the low-mass, `alpha=.4` benchmark:

- `a_FD = -0.00201850894575251 (km/s)^2/Mpc`
- `a_Fminus-a_Fplus = -1.25995524466503e-05 (km/s)^2/Mpc`
- signed F-state contrast = 0.624201% of `|a_FD|`
- reported finite-band FD contribution from `k>1 Mpc^-1` = 94.22589%
- complementary `k<=1` share in that same finite-band bookkeeping = 5.77411%.

The E28R1 numerical 24-vs-48 grid gap is only 0.419452% of the signed Fminus-Fplus difference, an absolute acceleration difference of about 5.2849e-08. This certifies quadrature of the chosen model, not the model itself.

A conservative sufficient condition is informative. If the entire unknown state-differential high-k modeling error were bounded by `epsilon*|a_FD,high|`, preserving the frozen finite-band sign is guaranteed only if

`epsilon < 0.00662452`

or about **0.662% of the high-k FD contribution**.

This is sufficient, not necessary. Correlated physical errors can cancel. The point is that the currently resolved numerical sign would require sub-percent physical state-differential short-scale control before becoming a robust physical statement.
