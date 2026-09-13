# Final manuscript integration plan

This document is the working baseline for the final Nature Physics-oriented Einstein–Vlasov manuscript. The numerical campaign is complete; manuscript construction and proof wording are now the only planned scientific tasks.

The exact final numerical anchors are recorded in [`FINAL_VALIDATION_STATUS_2026-09-13.md`](FINAL_VALIDATION_STATUS_2026-09-13.md).

## Central claim

The paper is not based on the generic observation that a kinetic distribution contains more information than a few moments. The sharper result is:

> **Instantaneous stress-energy is not a complete gravitationally relevant description of kinetic matter.**

The logical structure is

`exact source-level indistinguishability -> future gravitational separation -> massive-response injectivity -> nonlinear persistence -> controlled cosmological realization -> observability boundary`.

In equations, the central distinction is

`F1 != F2`, with the same instantaneous particle current and stress-energy tensor,

but different later gravitational response.

The manuscript must distinguish throughout:

`instantaneous sourcing != dynamical identifiability != stable inversion != practical observability`.

Do not use language implying that spacetime itself stores the hidden matter information. The additional state remains in the kinetic matter distribution.

## Analytic result chain

### 1. Exact source degeneracy

Construct distinct isotropic collisionless distributions with matched particle current and stress-energy tensor at the comparison time.

This is the starting counterexample, not the novelty claim by itself.

### 2. Massive linear-response injectivity

For fixed known `m>0`, nonzero wave number, response normalization and the stated weighted-decay class, equality of the ideal continuous causal TT response kernel on a non-empty time interval implies equality of the radial distribution.

Physical mechanism:

`v(p)=p/sqrt(p^2+m^2)`

is one-to-one, so different momentum regions acquire different dynamical phases.

Reviewer-proof wording requirements:

- state explicitly that particle mass, wave number and response normalization/coupling are fixed and known;
- in the Green-function corollary, state a known nonzero coupling and equality on a domain where the causal Laplace transforms exist;
- keep the weighted-decay/boundary assumptions next to the theorem statement.

### 3. Massless boundary

For isotropic massless matter `v=1`, radial momentum information collapses to the appropriate integrated moment in the TT response. This gives a clean structural boundary for the massive injectivity result.

Do not claim monotonic response strength in mass.

### 4. Compact forward map and unstable inverse

On every fixed finite momentum/response window the restricted forward map is compact. Since it is injective on an infinite-dimensional domain, the inverse on its range is unbounded.

This is the rigorous reason that exact identifiability does not imply robust finite-data recovery.

### 5. Nonlinear FLRW separation

Two distributions can have matched `(n,rho,P)` and therefore the same initial scale factor data through acceleration, while a different higher kinetic functional enters a subsequent derivative and separates the later geometry.

The current formulae for `dP/da`, `Delta dot P`, and `Delta a^(3)` have passed the adversarial algebraic audit.

### 6. Finite local source-jet theorem

For every finite order, distinct isotropic kinetic states can be constructed with the same finite local FLRW stress-energy jet while the next independent kinetic direction differs.

Use the phrase **“complete local FLRW stress-energy jet”** or **“within the isotropic FLRW symmetry class”**. Do not make an unrestricted statement about arbitrary tensor jets in general spacetimes.

## Cosmological realization

The publication-facing realization uses a controlled massive relic, not the older unconstrained mass sweep.

Frozen physical control:

- `m_ncdm=0.60 eV`;
- `z_match=1100`;
- matched `(n,rho,P)`;
- fixed present-day `omega_ncdm`;
- fixed `omega_cdm`;
- fixed early total `N_eff=3.046`;
- smooth positive deformations with a `+/-30%` cap;
- pinned CLASS commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.

The main fixed-cosmology anchor is

- joint ideal TT+TE+EE `S/N = 0.2502411773`;
- phi-phi `S/N = 0.1594145403`;
- max `P(k,z=0)` difference `0.0570892644%`.

The 10/20/30% frozen-direction test shows a maximum CV-weighted nonlinearity of `0.776%`; therefore the 30% result is not a strong-deformation artifact.

The older thermal-normalization, fixed-total-matter and fixed-rho-at-match sweeps are useful provenance/controls but should not drive the main narrative. They can be summarized in Methods, Extended Data or Supplementary Information if space permits.

## Observability result

The observational part is a boundary result, not a detection attempt.

### Primary CMB

Final local standard-cosmology marginalized values:

- Planck-like: `S/N = 0.01809`;
- SO LAT: `S/N = 0.04701`;
- CMB-S4: `S/N = 0.06664`.

About 96–97% of the fixed-cosmology `Delta chi^2` is absorbed by standard cosmological directions.

The ACT DR6 lensing check gives negligible endpoint `Delta chi^2` within the physical hidden-mode amplitude range.

### Three-dimensional linear matter clustering

At `k_max=0.3 h/Mpc`, total effective volume `50 (Gpc/h)^3`:

- fixed-cosmology optimized-pair `S/N = 0.58015`;
- five-parameter matter-power marginalized `S/N = 0.06523`;
- retained `Delta chi^2` fraction `1.264%`;
- absorbed fraction `98.736%`.

The five matter-power nuisance directions are

`H0, omega_b, omega_cdm, lnAs, n_s`.

`tau_reio` should not be presented as a physical matter-power nuisance direction. It was retained only in the numerical robustness audit. Removing it changes the answer negligibly. The five-parameter result is stable from pseudoinverse `rcond=1e-10` through `1e-14`; even the more aggressive `1e-8` cutoff changes `S/N` only to about `0.06634`.

Amplitude and tilt alone absorb about 98.6% of the optimized `k_max=0.3` signal norm, so the suppression is mainly a broad-shape degeneracy rather than a numerical near-null-direction effect.

Direct hidden-null-space optimization gives essentially no improvement over the frozen direction. Therefore the weak post-projection signal is not an accidentally poor choice of the tested hidden deformation.

### Interpretation

The preferred observational statement is:

> The hidden kinetic response survives all physical and numerical controls, but in the tested smooth isotropic relic class it lies predominantly along observational directions that are degenerate with standard cosmology.

The stronger conceptual conclusion is:

> Dynamical gravitational identifiability does not imply practical observability.

## Nature Physics framing

### Editorial danger

The most likely desk-rejection reading is:

> “A distribution function contains higher moments than the stress-energy tensor, so of course later observables can differ.”

The first page must explicitly explain why this is not the result.

The novelty is the combined chain:

- exact equality as instantaneous gravitational sources;
- later gravitational distinction;
- an injective ideal gravitational response map in the massive case;
- a sharp isotropic massless boundary;
- nonlinear finite-source-jet separation;
- a controlled cosmological realization;
- a quantitative boundary between identifiability and observational access.

### Preferred title direction

Working title:

**Kinetic matter carries gravitational information beyond instantaneous stress-energy**

Alternative if a slightly more conservative formulation is needed:

**Gravitational response distinguishes stress-energy-degenerate kinetic states**

### Abstract structure

1. Einstein equations use instantaneous stress-energy as the local source.
2. Ask if it is also a complete state for predicting the later gravitational influence of kinetic matter.
3. Give the source-degenerate/future-distinct result.
4. State massive injectivity and massless loss.
5. State nonlinear finite-jet persistence.
6. Give controlled cosmological realization.
7. End with the observational boundary: response exists, but most accessible CMB/LSS distinguishability is absorbed by standard cosmological directions.

Avoid beginning the abstract with neutrinos, CLASS or nonthermal relic spectra.

### Introduction structure

Paragraph 1: instantaneous Einstein source versus complete dynamical state.

Paragraph 2: explain that “distribution contains higher moments” is not the result; formulate exact source equivalence versus future gravitational distinction.

Paragraphs 3–4: massive injectivity and massless boundary.

Paragraph 5: nonlinear FLRW/finite-jet persistence.

Paragraph 6: identifiability versus stable inversion versus observability.

Paragraph 7: controlled cosmological realization and final observability result.

The first 1–1.5 pages should remain problem-first and largely independent of detailed cosmological numerics.

## Main-text result order

Recommended order:

1. **Stress-energy-degenerate kinetic states** — concise counterexample and statement of the problem.
2. **Dynamical recovery through massive response** — injectivity theorem and physical mechanism.
3. **Massless boundary and unstable finite-window inversion** — establish what injectivity does and does not imply.
4. **Nonlinear persistence beyond the instantaneous source** — FLRW separation plus finite local source-jet theorem.
5. **Controlled cosmological realization** — final fixed-abundance/fixed-early-radiation benchmark and amplitude convergence.
6. **Observability boundary** — primary CMB, idealized 3D LSS, nuisance projection and ACT check.
7. **Discussion** — state completeness, not hidden gravitational degrees of freedom.

## Figure hierarchy

Do not force all validation material into one four-panel main figure.

Recommended main figures:

### Figure 1 — concept and analytic mechanism

- two distributions with identical instantaneous source moments;
- different massive response kernels;
- schematic source projection versus dynamical response;
- optional massless-collapse inset.

### Figure 2 — nonlinear/state-completeness result

- matched FLRW initial source/geometry;
- subsequent geometric separation;
- finite source-jet matching illustration or hierarchy summary.

### Figure 3 — cosmological realization and observability boundary

- final controlled observable response;
- amplitude convergence;
- CMB fixed versus marginalized accessibility;
- LSS fixed versus marginalized accessibility.

Seven-dimensional singular-spectrum details, older abundance controls, full survey tables and nuisance/pseudoinverse audits are better suited to Extended Data/Supplementary Information.

## Claim guardrails

Do not claim:

- that non-uniqueness of moment descriptions is new;
- that different velocity distributions producing different observables is new;
- that the stress-energy tensor is incomplete as the instantaneous Einstein source;
- that spacetime itself stores hidden matter memory;
- that response strength increases monotonically with mass;
- that `m=0.60 eV` is a standard-neutrino best fit or realistic neutrino-mass claim;
- that the ideal/Fisher signals are experimental detections;
- that the tested smooth ten-function, 30%-cap relic class exhausts kinetic-state space;
- that the LSS forecast is a realistic galaxy-survey likelihood.

Preferred exact wording:

> The stress-energy tensor is the complete instantaneous Einstein source, but it is not a complete dynamical state variable for kinetic matter.

and

> Two kinetic states can be gravitationally indistinguishable now and gravitationally distinguishable later.

## Theorem wording audit status

A theorem-by-theorem adversarial pass found no result-breaking gap in:

- massive response injectivity;
- weighted-decay implication for the required moments;
- causal Green-function transfer corollary;
- compactness/unbounded inverse;
- FLRW pressure derivative and `a^(3)` separation;
- independence of the higher kinetic functional;
- finite local source-jet construction;
- isotropic massless collapse.

Remaining proof work is presentation-level:

- make all fixed/known response parameters explicit;
- qualify source-jet completeness by the FLRW symmetry class;
- make the causal/Laplace domain in the Green-function corollary explicit;
- keep theorem assumptions adjacent to claims in the main text or Methods.

## Numerical stop rule — CLOSED

The final amplitude test, corrected CMB survey forecast, LSS screen/recovery and adversarial nuisance/pseudoinverse audit have all completed successfully.

Do not initiate another broad CMB, lensing or LSS campaign for the present manuscript without a specific identified deficiency.

Remaining tasks:

1. replace the outdated cosmology section of the current TeX manuscript;
2. rewrite title, abstract and first 1–1.5 pages of the Introduction in the final framing;
3. integrate theorem wording fixes;
4. choose main versus Extended Data figures/tables;
5. perform final prose, citation and reproducibility audit;
6. build the submission package.

## Reproducibility anchors

- CLASS commit: `e85808324f51fc694d12e3ed7439552a3c3f9540`
- final controlled high-precision response: run `34694221873`
- hidden-subspace spectrum: run `34711814836`
- amplitude convergence: run `34715371731`
- final CMB survey forecast: run `34717947306`
- ACT DR6 real-data check: run `34694365955`
- original LSS CLASS screen: run `34722070280`
- finite-support LSS recovery: run `34733700067`
- LSS projection adversarial audit: run `34741320456`

The exact artifact IDs, digests and final values are in `FINAL_VALIDATION_STATUS_2026-09-13.md`.
