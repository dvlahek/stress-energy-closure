# Final manuscript validation status — 2026-09-13

This file is the final manuscript-facing numerical record for the Nature Physics-oriented Einstein–Vlasov state-completeness study. It supersedes older observational numbers in `FINAL_VALIDATION_STATUS_2026-09-12.md`, `MANUSCRIPT_PLAN.md`, and the README where they conflict.

The numerical campaign is now **locked**. No additional large CLASS campaign is justified unless manuscript review identifies a concrete implementation or physics issue.

## Scientific statement supported by the numerical campaign

The calculations support a distinction between three levels:

1. **instantaneous source equivalence** — distinct kinetic states can have the same matched low-order source moments;
2. **dynamical distinguishability / identifiability** — their later gravitational response can differ, and the massive ideal response is injective under the analytic assumptions proved in the manuscript;
3. **practical observability** — the response can be strongly compressed into directions degenerate with ordinary cosmological parameters.

The cosmological calculations are a controlled realization of the analytic result. They are not needed for the proof of injectivity or nonlinear source-jet separation.

## Frozen physical control

All final cosmological results use the controlled massive-relic benchmark

- `m_ncdm = 0.60 eV`;
- `z_match = 1100`;
- matched `(n, rho, P)`;
- fixed present-day `omega_ncdm = 0.0010752048549100347`;
- fixed `omega_cdm = 0.12`;
- fixed early total neutrino `N_eff = 3.046`;
- `deg_ncdm = 0.166669632040279`;
- `N_ur = 2.877130328816789`;
- smooth positive hidden deformations with a `+/-30%` pointwise cap;
- CLASS commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.

This control fixes both the present-day relic abundance and the baseline early radiation normalization. The remaining response is therefore associated with the intermediate massive transition and momentum-dependent kinetic evolution, not a trivial change in the present-day relic density or early radiation budget.

## High-precision controlled cosmological response

Workflow run: `34694221873`  
Artifact: `scalar-fixed-omega0-neff-global-highprec`  
Artifact SHA256: `abafddd4e64b92a7d62213c7ceb907ac74a6b79781fd9ee669d143ac2e1649f8`

For the final controlled `m=0.60 eV` scalar winner:

- joint TT+TE+EE ideal full-sky Gaussian-CV `S/N = 0.2502411773074913`;
- TT auto `S/N = 0.08907212006036914`;
- EE auto `S/N = 0.17711004466391997`;
- phi-phi auto `S/N = 0.1594145402996132`;
- max TT relative difference `0.024757307691490288%` at `ell=298`;
- max EE relative difference `0.03410878691272394%` at `ell=478`;
- max phi-phi relative difference `0.043427956493385274%` at `L=134`;
- max `P(k,z=0)` relative difference `0.057089264396426184%` at `k=0.032405016784 h/Mpc`.

The response is real but already sub-unity in ideal fixed-cosmology CMB distinguishability.

## Seven-dimensional hidden-subspace response spectrum

Workflow run: `34711814836`  
Artifact: `hidden-mode-7d-observable-spectrum-m060`  
Artifact ID: `10303790497`  
Artifact SHA256: `5b55738509a12c203e079f5330eb2b323ee0db5156dc020231524259422bef5a`

Using the stated number-density-weighted RMS fractional-deformation metric:

- TT+TE+EE: first hidden mode carries `83.6018546733%`; first two carry `99.2625186987%`;
- phi-phi: first mode carries `99.1516090273%`; first two carry `99.8771567326%`.

Thus the accessible response is effectively low-dimensional within the tested seven-dimensional stress-energy-null subspace. This statement is conditional on the chosen kinetic metric and local linear-response construction.

## Frozen-direction amplitude convergence

Workflow run: `34715371731`  
Artifact: `hidden-mode-amplitude-010-intermediate-convergence`  
Artifact ID: `10306188124`  
Artifact SHA256: `8253948133b1d3e892ab9b2bbdf951685d71fa64791abb78cda33389b6a72572`

The final high-precision hidden direction was rescaled without re-optimization to 10%, 20%, and 30% pointwise caps.

Distribution L2 separation over the Fermi–Dirac baseline:

- 10%: `0.14646581909543194`;
- 20%: `0.2929316381908638`;
- 30%: `0.4393974572862958`.

Joint TT+TE+EE ideal-CV response:

- 10%: `S/N = 0.08354300312105441`;
- 20%: `S/N = 0.1673035704362991`;
- 30%: `S/N = 0.2502411773074913`.

Matter-power maximum relative differences:

- 10%: `0.01902879682745334%`;
- 20%: `0.03805943408124809%`;
- 30%: `0.057089264396426184%`.

The largest CV-weighted departure from linear response over the tested range is

`0.007758579356208998 = 0.7758579356%`,

well below the prespecified 5% tolerance. The earlier anomalous 10% TTEE result was a numerical precision floor and disappears after tightening the CLASS precision settings.

**Manuscript statement:** the 30% benchmark is not a strong-deformation artifact; the controlled response is linear to sub-percent accuracy over the tested 10–30% range.

## Survey-aware primary-CMB observability

Workflow run: `34717947306`  
Artifact: `final-control-survey-forecast`  
Artifact ID: `10306331974`  
Artifact SHA256: `234e5b7b0cdbde9441a9f2b1e310752ebea5c71baea04b7be0afe79e2a963d44`

The forecast uses Gaussian TT/TE/EE covariance and a local linearized standard-cosmology Fisher projection.

### Planck-2018-like

- fixed-cosmology joint TTEE `S/N = 0.11058711796762878`;
- marginalized joint TTEE `S/N = 0.018090244932299306`;
- retained `Delta chi^2` fraction `0.026759612121736765`;
- absorbed fraction `0.9732403878782633`.

### Simons Observatory LAT baseline

- fixed `S/N = 0.27493459523702884`;
- marginalized `S/N = 0.04700669691164685`;
- retained `Delta chi^2` fraction `0.02923214527388671`;
- absorbed fraction `0.9707678547261133`.

### CMB-S4 reference

- fixed `S/N = 0.3351343134960461`;
- marginalized `S/N = 0.06664372383927371`;
- retained `Delta chi^2` fraction `0.03954401110761918`;
- absorbed fraction `0.9604559888923808`.

The actual pinned CLASS lensed-spectrum artifact was separately checked to be dimensionless `[ell(ell+1)/(2pi)] C_ell`; therefore the workflow conversion from microkelvin-arcmin map depth to dimensionless `Delta T/T` noise is correct.

**Conclusion:** the primary-CMB response survives the physical controls but is overwhelmingly degenerate with standard cosmological directions after local marginalization.

## ACT DR6 real-data lensing check

Workflow run: `34694365955`  
Artifact: `act-dr6-realdata-hidden-mode-fit`.

For the predefined physical hidden-mode amplitude `alpha in [-1,1]`:

- ACT-only endpoint `Delta chi^2 ~ 1.03e-4`;
- ACT+Planck endpoint `Delta chi^2 ~ 2.28e-4`.

The unconstrained linear best fit lies at an amplitude of order `10^7` and is outside the validated physical class, so its formal significance has no detection interpretation.

**Conclusion:** current lensing data cannot distinguish the tested hidden mode.

## Idealized three-dimensional LSS observability

Original CLASS screen: run `34722070280`. The expensive 29-case CLASS grid completed successfully; the original workflow failed only in post-processing because valid CLASS spectra have tiny endpoint differences in their finite `k` support.

Finite-support recovery: run `34733700067` — **SUCCESS**.  
Final artifact: `lss-observability-final-recovered`  
Artifact ID: `10311340767`  
Artifact SHA256: `d57c15f6324c47a45c809b9ae63d39f996afa594deb6a12c526680dafc7db071`.

The recovery uses the intersection of finite CLASS support across all spectra entering each Fisher block and performs no extrapolation.

The hidden direction was optimized in the seven-dimensional matched-moment null space for marginalized tomographic matter-power distinguishability, with fresh CLASS validation of the optimized plus/minus pair.

Optimized pair checks:

- relative number-density mismatch `3.579686991212226e-16`;
- relative energy-density mismatch `0`;
- relative pressure mismatch `0`;
- pointwise deformation cap `+/-30%`;
- distribution L2 separation over FD `0.43939745728629526`.

For `k_max = 0.3 h/Mpc` and total effective volume `50 (Gpc/h)^3`:

- optimized fixed-cosmology `S/N = 0.5801527803698333`;
- six-direction numerical Fisher projection `S/N = 0.06521388853657625`;
- retained `Delta chi^2` fraction `0.012635587450189119`;
- absorbed fraction `0.9873644125498109`.

For `100 (Gpc/h)^3`, the corresponding marginalized value is only `S/N = 0.09222636562351344`.

The optimized direction gives essentially no gain over the frozen CMB-selected direction, showing that the weak post-marginalization signal is not caused by an accidentally poor hidden-state choice within the tested local null space.

## LSS Fisher adversarial audit

Workflow run: `34741320456` — **SUCCESS**.  
Artifact: `lss-projection-adversarial-audit`  
Artifact ID: `10312453719`  
Artifact SHA256: `c738417e5bd3507dd9bbf9f7575705f329a0723c0941ff71214ec3df83f4c037`.

This audit performs **no new CLASS calculations**. It reprojects the already validated spectra on exactly the same common finite `k` support while changing nuisance subsets and pseudoinverse cutoffs.

The matter-power field has no physically useful reionization-optical-depth degree of freedom, so the manuscript-facing LSS marginalization should use the five physically relevant local matter-power directions

`H0, omega_b, omega_cdm, lnAs, n_s`,

with `tau_reio` retained only as a numerical robustness check.

At `k_max=0.3 h/Mpc`, `V=50 (Gpc/h)^3`, for the validated optimized pair:

- fixed `S/N = 0.5801527803698042`;
- five-parameter, no-tau marginalized `S/N = 0.06523236275135469`;
- retained `Delta chi^2` fraction `0.012642747447111056`;
- absorbed fraction `0.987357252552889`.

The old six-direction value is `0.0652138885344226`, a negligible difference.

Pseudoinverse sensitivity for the physically relevant five-parameter projection:

- `rcond=1e-8`: `S/N = 0.06633718905835097`;
- `rcond=1e-10`: `S/N = 0.06523236275135469`;
- `rcond=1e-12`: `S/N = 0.06523236275135469`;
- `rcond=1e-14`: `S/N = 0.06523236275135469`.

Thus the conclusion is insensitive to the pseudoinverse cutoff. The audit also reproduces the six-direction reference value to an absolute error of `2.15e-12`.

Amplitude and tilt alone already absorb about `98.6%` of the optimized signal norm at `k_max=0.3`, showing that the observational suppression is primarily a broad-shape degeneracy, not an artifact of the near-null `tau_reio` derivative.

**Final LSS conclusion:** within the tested smooth, moderate-amplitude isotropic relic class, three-dimensional linear matter clustering contains a real hidden-state response, but standard cosmological shape directions absorb almost all of its ideal fixed-cosmology distinguishability.

## Final observational interpretation

The final numerical chain is therefore

`controlled hidden-state response -> converged amplitude scaling -> survey projection -> 3D LSS projection -> adversarial nuisance audit`.

Both primary CMB and idealized linear three-dimensional matter clustering support the same distinction:

> dynamical gravitational information can exist without being stably or practically accessible in the chosen observable sector.

This is an observability-bound result, not a failed detection claim.

## Numerical stop rule — CLOSED

The prespecified amplitude-convergence test passed. The corrected primary-CMB forecast completed. The LSS screen and finite-support recovery completed. The final LSS nuisance/pseudoinverse adversarial audit passed.

**No further large numerical campaign should be started for this manuscript unless a concrete problem is identified during manuscript or referee review.**

Remaining work is manuscript construction, theorem wording audit, figure/table selection, and repository/source-data packaging.
