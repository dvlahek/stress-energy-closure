# Final manuscript validation status — 2026-09-12

This file records the manuscript-facing numerical anchors after the fixed-present-day relic abundance plus fixed-early-Neff control. It supersedes older observational numbers in `MANUSCRIPT_PLAN.md` where they conflict.

## Final controlled high-precision anchor

Source workflow run: `34694221873`  
Artifact: `scalar-fixed-omega0-neff-global-highprec`  
Artifact SHA256: `abafddd4e64b92a7d62213c7ceb907ac74a6b79781fd9ee669d143ac2e1649f8`

Final winner:

- mass: `0.60 eV`;
- hidden direction: scalar-optimized winner;
- `omega_cdm = 0.12`;
- fixed present-day `omega_ncdm = 0.0010752048549100347`;
- `deg_ncdm = 0.166669632040279`;
- `N_ur = 2.877130328816789`;
- controlled early total neutrino `N_eff = 3.046`;
- CLASS commit: `e85808324f51fc694d12e3ed7439552a3c3f9540`.

High-precision response:

- joint TT+TE+EE ideal full-sky Gaussian CV `S/N = 0.2502411773074913`;
- TT auto `S/N = 0.08907212006036914`;
- EE auto `S/N = 0.17711004466391997`;
- phi-phi auto `S/N = 0.1594145402996132`;
- maximum TT relative difference `0.024757307691490288%` at `ell=298`;
- maximum EE relative difference `0.03410878691272394%` at `ell=478`;
- maximum phi-phi relative difference `0.043427956493385274%` at `L=134`;
- maximum `P(k,z=0)` difference `0.057089264396426184%` at `k=0.032405016784 h/Mpc`.

Manuscript interpretation: the hidden-state response is real and numerically stable after both present-day abundance and early radiation are fixed, but the ideal distinguishability remains below unity. This is the final primary cosmological anchor unless a validation step reveals an implementation error.

## Seven-dimensional hidden-subspace observable spectrum

Workflow run: `34711814836`  
Artifact: `hidden-mode-7d-observable-spectrum-m060`  
Artifact ID: `10303790497`  
Artifact SHA256: `5b55738509a12c203e079f5330eb2b323ee0db5156dc020231524259422bef5a`

Method: generalized eigenproblem `G v = lambda H v`, where `G` is the ideal response Gram matrix and `H` is the number-density-weighted RMS fractional-deformation metric

`int q^2 f0 (delta f/f0)^2 dq / int q^2 f0 dq`.

The response has numerical rank seven but is effectively low-dimensional:

- TT+TE+EE: first mode carries `83.6018546733%` of the metric-defined information; first two carry `99.2625186987%`;
- phi-phi: first mode carries `99.1516090273%`; first two carry `99.8771567326%`.

This supports the statement that the accessible gravitational response is concentrated in one or two kinetic combinations and is not an artifact of inspecting a single optimizer-selected direction. The information fractions are conditional on the stated kinetic metric and the local linear-response approximation.

## ACT DR6 real-data check

Workflow run: `34694365955`  
Artifact: `act-dr6-realdata-hidden-mode-fit`.

Within the physical predefined hidden-mode range `alpha in [-1,1]`:

- ACT-only endpoint difference: `Delta chi^2 ~ 1.03e-4`;
- ACT+Planck endpoint difference: `Delta chi^2 ~ 2.28e-4`;
- profiling an overall lensing amplitude makes the shape distinction still smaller.

The unconstrained linear best fit lies at `alpha ~ 10^7`; the associated formal `~27 sigma` ratio is not a detection because it is far outside the validated physical amplitude range. The real-data conclusion is that current ACT/ACT+Planck lensing cannot distinguish this tested hidden kinetic mode.

## Amplitude-linearity control

Workflow run: `34711963013`.

The frozen final `m=0.60 eV` direction is being tested at 10%, 20%, and 30% pointwise deformation caps without re-optimization. The existing validated 30% high-precision spectra are reused; only the four 10%/20% endpoint spectra are recomputed.

Input pair checks already pass at machine precision:

- 10% maximum `(n,rho,P)` pair mismatch: `2.085155943954705e-16`;
- 20% maximum `(n,rho,P)` pair mismatch: `1.193228997070742e-16`;
- 30% maximum `(n,rho,P)` pair mismatch: `4.772915988282968e-16`;
- distribution L2 separations: `0.1464658191`, `0.2929316382`, `0.4393974573`, exactly proportional to amplitude.

Final manuscript use depends on the completed CLASS response check. The prespecified pass condition in the workflow is <=5% CV-weighted deviation from linear scaling for both TT+TE+EE and phi-phi.

## Final survey-aware forecast

Workflow run: `34712036907`.

This run is explicitly gated on successful completion of the amplitude-linearity workflow. It uses the final fixed-omega0 plus fixed-early-Neff winner and corrects both issues in the older survey pipeline:

1. map depths in microkelvin-arcmin are converted to dimensionless `DeltaT/T` noise before combination with CLASS spectra;
2. CLASS inputs use the final-control `N_ur = 2.877130328816789` instead of the historical `2.0328`.

The forecast includes Planck-like, Simons Observatory LAT baseline, and CMB-S4 reference configurations, with fixed-cosmology and six-parameter linearized LCDM-marginalized TT/TE/EE results. It remains a Gaussian/Fisher forecast, not an official experimental likelihood or a detection claim.

## Superseded Planck nonlinear run

Run `34693285951` was cancelled during the first nonlinear bounded Planck likelihood optimization. It used an older control and a computationally impractical design (`lmax=2500`, roughly 20-minute CLASS evaluations, `maxfev=120`). It is not a manuscript blocker and should not be rerun in this form. A nonlinear Planck profile can be added only if specifically required by review.

## Stop rule

After the amplitude-linearity and corrected final survey forecast complete successfully, do not start another large numerical campaign unless either result exposes a concrete unresolved problem. At that point the numerical campaign is considered locked for manuscript integration.
