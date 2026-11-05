# E30/E32 — offline derivation and arithmetic QA protocol

Date: 2026-09-29. This is an explicitly **post-existing-results mathematical follow-up**. It is not a prospective experiment on new physical data, not an observed result, and does not revise frozen E8–E29 protocols. The toy arithmetic below is chosen to test the derived identities and is never reported as a simulated or measured galaxy/halo signal.

## Inputs and scope
- Read-only existing project documents: E10 (conditional wind parity), E22 (linear vs quadratic galaxy response), E24 (linear two-tracer sign degeneracy), E25 (total-metric Euler bookkeeping), E26 (one-adiabatic-mode rank), E28 (UV-dominated finite-band benchmark), E29 (evolving source/test profiles) and A-04 (eBOSS 24D readiness). Existing branch is `audit/eboss-elg-bit8-ra-orientation-20260925`, draft PR #1; main unchanged.
- Only user-uploaded small E29 JSON may be read for metadata integrity (optional). No ASDF, PID arrays, FITS, observed galaxies, observed odd, mock retrieval, network calls, simulation re-runs, new science cuts, alternate F states or tuned physical cutoffs.
- E30 derives the exact **conditional** sign-cancellation criterion, and distinguishes it from nonzero Gaussian two-field density–velocity covariance and from Gaussian-leading zero of the mixed three-field correlator. It exhibits the unknown high-z tracer-response contrast and the exact 24D vs nine-mock rank obstacle.
- E32 constructs two positive radial shell distributions with exactly the same 10 published Abacus L2 percentile radii and `r100`, but different spherical Fourier transforms. It derives a conservative percentile-only Fourier enclosure from `sinc(x)=integral_0^1 cos(tx) dt` and `|sinc'(x)|<=1/2`. All radii and wavenumbers in the toy examples are dimensionless.

## Exact QA / fail-closed controls
1. E30 enumerates finite rational ensembles: symmetric marginal wind with density-correlated wind; linear two-field expectation nonzero; mixed three-field expectation zero under central symmetry; conditional sign-cancellation with balanced signs; explicit counterexample to marginal balance alone; signed selection state-specific asymmetry; two-tracer contrast, null, reversal, and nine-mock centered-rank 8.
2. E32 checks mass normalization, all 10 exact cumulative quantiles, that transforms differ at stated positive toy wave numbers, every computed profile lies inside the proven Lipschitz envelope, and that an intentionally perturbed mass distribution fails the matching-quantile check.
3. Any failed assertion stops with nonzero exit and produces no success report. Save separate original script and JSON outputs. For final user-facing statements, distinguish proven algebra, illustrative toy numbers, source-derived prior results, and remaining physical unknowns.

## Physical STOP
No actual high-z `c_{a,F}`, `beta`, `epsilon`, `s_a^w`, state-dependent signed-selection law, `M200c(a)`, physical `u_src/test`, second-order Vlasov force, UV closure, eBOSS `xi_1,xi_3` or A-04 covariance is inferred here. Opaženi eBOSS 24D odd ostaje SEALED; A-03 PHYSICAL_UNCERTIFIED, A-04 BLOCKED. Do not change GitHub, `main` or draft PR in this offline task.
