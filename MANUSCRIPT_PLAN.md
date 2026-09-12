# Manuscript and final-validation plan

This document is the working baseline for the final Nature Physics-oriented version of **Gravitational response distinguishes hidden kinetic states**. It records the interpretation that should be preserved while the remaining validation jobs are completed.

## Current physical baseline

The central theoretical result remains unchanged:

- distinct kinetic states can share the same instantaneous particle current and stress-energy tensor;
- for massive isotropic collisionless matter, the continuous causal TT response kernel is injective with respect to the radial kinetic distribution under the stated assumptions;
- the isotropic massless limit loses radial momentum information because all momenta propagate at the same speed;
- finite source jets can be matched to arbitrary finite order while the next kinetic direction, and therefore the later metric evolution, remains distinct.

The cosmological CLASS campaign now separates the hidden-state response from background and abundance effects.

### Original optimized massive-relic benchmark

For the original thermal normalization, the high-precision global winner occurs at `m = 0.60 eV` and gives

- ideal full-sky joint TT+TE+EE `S/N = 1.5193052363`;
- ideal full-sky phi-phi `S/N = 0.9470909494`;
- maximum `P(k,z=0)` difference `0.3391697439%`.

This benchmark changes both kinematics and the relic gravitational weight, so it must not be interpreted as a pure mass effect.

### Fixed-total-matter control

Holding the present-day total matter density fixed leaves the result essentially unchanged:

- high-precision joint TT+TE+EE `S/N = 1.5315759525` at `m = 0.60 eV`.

Therefore the large response is not explained by the trivial change in the total present-day matter density.

### Fixed relic energy density at the matching epoch

Holding the Fermi-Dirac relic energy density fixed at `z_match = 1100` reduces, but does not remove, the response:

- high-precision joint TT+TE+EE `S/N = 1.0735578430` at `m = 0.60 eV`.

This control equalizes the hidden component's gravitational weight at the matching epoch, but it does not keep its present-day abundance fixed.

### Strict fixed-present-day relic abundance control

The strict control keeps

`omega_ncdm(z=0) = 0.0010752048549100347`

fixed across the mass sweep while also keeping `omega_cdm = 0.12`. The high-precision global winner is again `m = 0.60 eV`, now along the lensing-optimized hidden direction:

- joint TT+TE+EE ideal full-sky `S/N = 0.2793350592`;
- TT auto `S/N = 0.1033954143`;
- EE auto `S/N = 0.2043851583`;
- phi-phi auto `S/N = 0.1748187606`;
- maximum TT difference `0.0260067%`;
- maximum EE difference `0.0540628%`;
- maximum phi-phi difference `0.0442833%`;
- maximum `P(k,z=0)` difference `0.0577727%`.

Thus about `18.4%` of the original ideal joint CMB signal remains after the present-day relic abundance is fixed exactly.

The intermediate strict-control sweep is

| mass [eV] | optimized ideal S/N |
|---:|---:|
| 0.03 | 0.07407165 |
| 0.06 | 0.07149882 |
| 0.10 | 0.07236767 |
| 0.18 | 0.08632430 |
| 0.30 | 0.09853171 |
| 0.60 | 0.12107608 |

The light/relativistic regime is therefore approximately flat, followed by a systematic opening of the response toward the more massive regime. Do **not** state a general theorem that response strength is monotonic in mass. The theorem establishes massive injectivity versus the isotropic massless degeneracy, not monotonic sensitivity.

## Working physical interpretation

The cosmological enhancement contains at least two contributions:

1. **gravitational weighting/background evolution**, which explains most of the large original `m = 0.60 eV` amplitude;
2. **mass-dependent kinetic response**, which remains nonzero after the present-day relic abundance and total matter budget are fixed.

The preferred manuscript statement is:

> The large cosmological enhancement is partly gravitational-weight driven, but a distinct mass-dependent kinetic response remains after fixing the relic abundance.

A stronger causal statement should be made only after the full 7D hidden-subspace spectrum is examined.

## Remaining validation sequence

### 1. Survey-aware Planck / SO / CMB-S4 forecast

Run ID: `34677322082`.

The forecast includes finite sky fraction, beam, white instrumental noise, experiment-specific multipole cuts, joint TT/TE/EE Gaussian covariance, and linearized marginalization over

`H0, omega_b, omega_cdm, ln(A_s), n_s, tau`.

The main reported quantities are

- fixed-cosmology survey-aware `S/N`;
- LCDM-marginalized survey-aware `S/N`.

This forecast determines the strength of the observational paragraph. It does not determine the validity of the identifiability result.

### 2. Full 7D hidden-subspace spectrum

Use all seven stress-energy-null directions already computed at each mass. Construct the covariance-whitened local response Gram matrices and their singular spectra.

Primary diagnostic:

`singular_values(m/T_nu(z_match))` for joint TT+TE+EE and phi-phi.

The aim is to distinguish a structural opening of the inverse problem from a single optimizer-selected direction. The most useful result would be a systematic change in the leading singular modes or conditioning as the system leaves the relativistic regime.

### 3. Amplitude-linearity control

For the same `m = 0.60 eV` strict fixed-omega0 winning direction, evaluate pointwise caps of

`10%, 20%, 30%`.

Check if the observable response and ideal S/N scale approximately linearly with amplitude. This excludes the possibility that the result is driven by a special large-deformation corner.

### 4. Stop criterion

After the survey forecast, hidden-subspace spectrum, and amplitude control are complete, do not start further large numerical campaigns unless one of these tests exposes a specific unresolved problem.

The current high-precision CLASS checks already establish numerical stability of the principal benchmark and controls.

## Final four-panel characterization figure

The final characterization workflow is stored in

`.github/workflows/class_response_mass_matrix.yml`

and is now repurposed as a manual final-validation workflow.

The intended main/Extended Data figure is:

### Panel (a): background and abundance controls

Optimized ideal response versus relic mass for

- original thermal normalization;
- fixed total matter;
- fixed relic rho at `z=1100`;
- strict fixed `omega_ncdm(z=0)`.

Purpose: separate the contribution of background/relic weighting from the residual mass-dependent response.

### Panel (b): hidden-subspace singular spectrum

Leading covariance-whitened singular modes versus

`m/T_nu(z_match)`.

Purpose: test if response accessibility is a property of the full stress-energy-null subspace rather than one optimized direction.

### Panel (c): amplitude scaling

High-precision response for 10%, 20%, and 30% caps along the same strict-control winning direction.

Purpose: test local linearity and deformation robustness.

### Panel (d): experimental accessibility

Planck-like, SO LAT, and CMB-S4 reference forecasts, each shown before and after LCDM marginalization.

Purpose: distinguish mathematical identifiability, ideal distinguishability, and survey accessibility.

If the marginalized survey signal is very small, panel (d) should move to Extended Data and the hidden-subspace spectrum should receive greater emphasis in the main text. If SO or CMB-S4 retains a meaningful marginalized signal, panel (d) can remain in the main characterization figure.

## Manuscript integration

### Abstract

Keep the theorem as the main claim. Add at most one sentence stating that controlled cosmological calculations retain a residual massive-state response after relic-abundance matching. Mention survey accessibility only if the marginalized forecast is genuinely informative.

### Introduction

Explicitly separate three facts:

1. a kinetic distribution contains more information than its low moments;
2. distinct kinetic states can share the same stress-energy tensor;
3. the new question is if gravitational dynamics itself can become injective with respect to the hidden kinetic information.

Add the closest prior art on Einstein-Vlasov degeneracy and nonthermal-relic higher moments. Use a careful `to our knowledge` first claim for the inverse gravitational response result, not for the general non-uniqueness of moment descriptions.

Preferred novelty statement:

> To our knowledge, no previous work has established that the gravitational response itself can become injective with respect to kinetic information absent from the instantaneous particle current and stress-energy tensor.

### Linear theory

Preserve the current injectivity theorem. Strengthen the physical bridge through

`v(p) = p/sqrt(p^2+m^2)`.

In the massless limit `v=1`, radial momentum information collapses. For massive particles, momentum-dependent velocities create distinct temporal/scale signatures.

### Nonlinear FLRW result

Keep the finite-source-jet theorem as an independent nonlinear result. It demonstrates that the central issue is not an artifact of linear response theory.

### Cosmological realization

Replace any simple statement that larger mass causes a larger signal with the controlled decomposition:

- original thermal sweep;
- fixed-total-matter control;
- fixed-rho-at-match control;
- strict fixed-present-day-relic-abundance control.

The main conclusion should be that abundance/background evolution amplifies the response strongly, while a smaller but distinct massive-regime kinetic contribution remains.

### Observability section

Maintain the hierarchy

`identifiability != stable inversion != observational detectability`.

Report ideal cosmic-variance numbers separately from survey-aware and marginalized numbers. Do not describe an ideal or Fisher-level result as a detection.

### Discussion

The preferred conceptual statement is:

> The stress-energy tensor is a complete instantaneous Einstein source, but it is not a complete dynamical state of kinetic matter. Information removed by moment compression need not remain gravitationally hidden, because it can reappear in the later response.

This should lead naturally to the interpretation as gravitational spectroscopy/tomography of kinetic state, without saying that spacetime itself stores a memory of hidden matter information.

## Claim guardrails

Do not claim:

- that stress-energy non-uniqueness is new;
- that different velocity distributions producing different observables is new;
- that response strength must increase monotonically with particle mass;
- that `m=0.60 eV` is a realistic standard-neutrino cosmology;
- that ideal or survey Fisher `S/N` is an experimental detection;
- that the current smooth 10-function, 30%-cap class exhausts all possible hidden-state spectra.

The `m=0.60 eV` case is a controlled massive-relic benchmark used to expose the mechanism.

## Reproducibility anchors

- CLASS commit: `e85808324f51fc694d12e3ed7439552a3c3f9540`
- original optimized mass sweep: run `34611544016`
- fixed-total / fixed-rho-at-match controls: run `34642785031`
- strict fixed-omega0 control: run `34674788325`
- survey-aware likelihood forecast: run `34677322082`
- strict fixed-omega0 high-precision artifact: `scalar-fixed-relic-omega0-global-highprec`
- strict fixed-omega0 high-precision artifact digest: `sha256:73088d44c0206231687f0bc871ad9a05f0e55a80689b11ba08cf7f3e14fd0867`
