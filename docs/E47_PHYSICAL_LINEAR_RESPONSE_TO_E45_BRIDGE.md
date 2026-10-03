# E47 — physical linear-response basis for the E45 practical-sensitivity metric

E47 is the first bridge that can put a physically dimensioned Einstein–Vlasov source difference into the same **dimensionless eBOSS odd-vector space** used by E45 without pretending that the unknown LRG/ELG response is already calibrated.

The restricted E25 local linear model is written with two dimensionless combinations

`g0 = b_L d_E - b_E d_L`

and

`g2 = f (d_E-d_L)`,

where `d_a` is the E25 local response to the dimensionless neutrino–CDM LOS velocity. Under the deliberately restricted assumption that the tracer response is the same functional response for F+ and F−, the source-state difference is

`Delta Im P_LE = mu [g0 + g2 mu^2] Delta C_U(k)`,

with `Delta C_U=(C_v,Fminus-C_v,Fplus)/c`.

Using `mu^3 = 3/5 P1 + 2/5 P3`, the two physical pre-window basis multipoles are

- unit `g0`: `xi1 = -H1[Delta C_U]`, `xi3=0`
- unit `g2`: `xi1 = -(3/5)H1[Delta C_U]`, `xi3=(2/5)H3[Delta C_U]`

for the Fourier/correlation convention used by the project.

The local script:

1. reads the exact archived E8 4000q F+ and F− distributions;
2. reruns CLASS source-only direct `vTk` and verifies the four E21 source points before any projection;
3. loads only the exact 250104-byte E0/E1 random-window NPZ;
4. projects the two unit physical bases through each fixed mock/cap conditional RR operator;
5. extracts each E3/E45 synthetic matched-response vector `q`;
6. returns, for every fixed case,

`A_E45 = alpha_case*g0 + beta_case*g2`;

7. reports the descriptive magnitude of `g0` or `g2` required to cross the E45 83% and 94% symmetric sign-recovery levels when either parameter is considered alone.

This is **not** a fitted physical LRG/ELG model. The useful result is the transfer function from physical response coefficients to the already measured practical-sensitivity axis. If the required `g0`/`g2` values are implausibly large compared with an independently derived galaxy response, the original unconditioned 2pt route is practically weak. If they are naturally small enough, the route merits full physical closure and larger mock covariance.

No observed odd vector, observed galaxy rows, FITS, ASDF or inference is used.
