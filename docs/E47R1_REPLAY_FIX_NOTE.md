# E47R1 — correction of the E21 replay gate

The first E47 local run stopped at the E21 source replay with a 5.6577% Fplus gap.

This was traced to an implementation error in E47, not a failure of the frozen E21 result. The first E47 script formed the nonlinear quantity `C_v(k)` on the native CLASS transfer grid and then interpolated that already-multiplied product to the four E21 long modes. The original E12/E13 construction instead interpolates the primitive transfer fields (`delta_cb` and `theta_ncdm-theta_cdm`) to the requested `K` first and only then forms `P_R delta_cb c theta/k`.

At the very low E21 long modes, interpolation and nonlinear multiplication do not commute well enough for the strict replay test, producing the observed ~5.7% discrepancy.

E47R1 restores the original operation order and also restores the mandatory E13 `W_R16(K)` top-hat factor. The exact replay target remains unchanged and the tolerance is not relaxed.

If E47R1 passes the frozen four-K E21 replay, only then is the dense source evaluated and projected through the E0/E1 window. No observed data are used.
