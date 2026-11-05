# E40 — how much F-state information survives once the unknown tracer response is separated from the source

For a linear source channel

`y_j = chi_F(K_j) C_v,F(K_j)`,

if `chi_F(K_j)` is an unconstrained state-specific response, **F is exactly non-identifiable**: for any target vector `y`, choose `chi_F,j=y_j/C_F,j`. Precision in `y` cannot recover the state without an external response constraint.

Under the much stronger toy restriction that the same single scalar response amplitude multiplies all four frozen E21 long modes, the `F+` and `F-` vectors are not exactly collinear. Their angle is

`theta_E21 = 0.000968554696135 rad`

and after optimizing one scalar amplitude the residual norm is

`sin(theta_E21) = 0.000968554544588`,

about **0.0969% of the source-vector norm**. Componentwise F+/F- differences are 0.712%, 0.481%, 0.188%, 0.182% from `K=.001` to `.005 h/Mpc`.

For comparison only, the frozen E19 conditional `(T1,T3)` vectors at one `k=.05 h/Mpc` have angle `0.02384118 rad`, an amplitude-free separation about **24.6 times larger** than the E21 four-K source geometry.

This does not make E19 an eBOSS observable. E19 is velocity-conditioned; E21 is a different linear density-velocity source. Neither has the physical high-z LRG/ELG response, survey window or A04 covariance.
