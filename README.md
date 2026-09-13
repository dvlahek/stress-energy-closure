# Gravitational response distinguishes hidden kinetic states

Code and source data accompanying the Einstein–Vlasov kinetic-state completeness study.

The central question is not simply that a kinetic distribution contains more information than its low-order moments. The paper asks a sharper gravitational question: can two collisionless kinetic states be exactly indistinguishable as instantaneous Einstein sources and still become gravitationally distinguishable later?

The answer is yes. Distinct kinetic states can share the same particle current and stress-energy tensor at an initial time while retaining different momentum-space structure. For massive isotropic collisionless matter, at fixed known particle mass and response normalization, the ideal continuous causal transverse-traceless response kernel at any fixed nonzero wave number uniquely identifies the radial kinetic distribution in the stated weighted-decay class. In the isotropic massless limit, all momenta propagate at the same speed and radial-profile injectivity is lost.

The nonlinear FLRW construction gives the complementary result: matched instantaneous source data, and more generally matched finite local FLRW source jets, do not exhaust the kinetic information relevant to later geometry. The manuscript therefore separates

`instantaneous source equivalence != dynamical identifiability != stable inversion != practical observability`.

On every fixed finite momentum/response window the restricted gravitational forward map is compact, so its inverse is unbounded even though the ideal massive response is injective.

## Final numerical status

The numerical campaign is **locked**. The authoritative manuscript-facing record is [`FINAL_VALIDATION_STATUS_2026-09-13.md`](FINAL_VALIDATION_STATUS_2026-09-13.md).

All final cosmological results use the controlled `m_ncdm=0.60 eV` benchmark with

- matched `(n,rho,P)` at `z_match=1100`;
- fixed `omega_ncdm(z=0)=0.0010752048549100347`;
- fixed `omega_cdm=0.12`;
- fixed early total neutrino `N_eff=3.046`;
- smooth positive hidden distributions with a `+/-30%` pointwise cap;
- CLASS pinned to commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.

The final high-precision controlled response gives

- ideal full-sky joint TT+TE+EE `S/N = 0.2502411773`;
- phi-phi `S/N = 0.1594145403`;
- maximum `P(k,z=0)` relative difference `0.0570892644%`.

A frozen-direction 10/20/30% convergence test gives a maximum CV-weighted departure from linear scaling of only `0.776%`, below the prespecified 5% tolerance.

Survey-aware primary-CMB forecasts show strong standard-cosmology degeneracy:

| Forecast | fixed TTEE S/N | marginalized TTEE S/N | absorbed Delta chi^2 |
|---|---:|---:|---:|
| Planck-like | 0.11059 | 0.01809 | 97.32% |
| SO LAT baseline | 0.27493 | 0.04701 | 97.08% |
| CMB-S4 reference | 0.33513 | 0.06664 | 96.05% |

The idealized tomographic three-dimensional linear-matter screen reaches `S/N = 0.58015` at fixed cosmology for `k_max=0.3 h/Mpc` and total effective volume `50 (Gpc/h)^3`, but only `S/N = 0.06523` after projection over the five physically relevant local matter-power directions `(H0, omega_b, omega_cdm, lnAs, n_s)`. About `98.74%` of the fixed-cosmology `Delta chi^2` is absorbed.

A dedicated post-processing adversarial audit confirms that this LSS suppression is insensitive to removing the near-null `tau_reio` direction and to the pseudoinverse cutoff. No additional CLASS calculation was needed for that audit.

The observational conclusion is therefore not a detection claim. The controlled cosmological realization demonstrates a real response, while the CMB and linear-LSS calculations identify a strong boundary between dynamical information and practical recovery.

## Main scripts

- `code/ev_flrw_controls.py` — exact flat-FLRW Einstein–Vlasov controls, including the massive same-`N^mu`/same-`T_munu` construction and massless profile-universality check.
- `code/injectivity_tomography.py` — stress-energy-matched response construction and regularized recovery of the hidden radial distribution.
- `code/prediction_transfer.py` — model-level amplitude, phase and complex metric-response separation.
- `code/noise_sweep.py` — deterministic finite-data stability diagnostics.
- `code/mass_sweep.py` — matched matter pairs across particle mass and FLRW geometry separation.
- `code/hierarchy_test.py` — finite local source-jet matching with compact-support bump functions.
- `code/direct_vs_memory.py` — direct phase-space evolution compared with the eliminated retarded-response representation.
- `code/class_scalar_response_optimize.py` — scalar CMB/lensing response optimization in the exact stress-energy null space.
- `code/class_scalar_fixed_omega0_neff_run.py` — final fixed-present-day-abundance plus fixed-early-`N_eff` control.
- `code/class_cmb_survey_likelihood_forecast.py` and final-control wrappers — Planck-like, SO and CMB-S4 Gaussian/Fisher forecasts.
- `code/class_lss_observability_campaign.py` — tomographic idealized 3D matter-power observability campaign.
- `code/class_lss_observability_campaign_finite_support.py` — common-finite-support recovery, with no extrapolation.
- `code/lss_projection_adversarial_audit.py` — nuisance-subset and pseudoinverse sensitivity audit using already validated spectra.

## Publication-facing validation runs

Key frozen runs are:

- `34694221873` — final fixed-`omega_ncdm(z=0)` + fixed-early-`N_eff` high-precision control;
- `34711814836` — seven-dimensional hidden-subspace observable spectrum;
- `34715371731` — final 10/20/30% amplitude convergence;
- `34717947306` — final survey-aware primary-CMB forecast;
- `34694365955` — ACT DR6 real-data lensing check;
- `34722070280` — original 29-case LSS CLASS screen;
- `34733700067` — successful finite-support LSS recovery and fresh optimized-pair validation;
- `34741320456` — final LSS Fisher adversarial audit.

The exact artifact names, digests and manuscript-facing values are recorded in [`FINAL_VALIDATION_STATUS_2026-09-13.md`](FINAL_VALIDATION_STATUS_2026-09-13.md).

## Reproducing the non-CLASS calculations

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

The CLASS calculations require a local CLASS build or can be reproduced through the supplied GitHub Actions workflows. Fixed random seeds are used only for synthetic-noise realizations and response-optimization search heuristics.

## What the numerical calculations establish

1. Smooth isotropic massive Vlasov states can share the same initial particle current and stress-energy tensor while producing different later FLRW geometry.
2. Massive states with the same matched low-order source moments have distinct causal TT response kernels.
3. Finite noisy inversion illustrates instability; exact injectivity, compactness and the unbounded inverse are analytic results.
4. In the isotropic massless limit, radial-profile dependence collapses as predicted analytically.
5. Arbitrarily many local FLRW source derivatives can be matched while a subsequent kinetic direction remains distinct.
6. Direct phase-space evolution and the reduced retarded-response representation agree under refinement.
7. Controlled CLASS calculations retain a small but converged hidden-state response after both present-day relic abundance and early radiation normalization are fixed.
8. Primary CMB and idealized linear 3D matter clustering strongly suppress practical distinguishability after standard-cosmology projection.

The numerical calculations illustrate and validate the analytic arguments. The source-level degeneracy, massive-response injectivity, massless boundary, nonlinear FLRW separation and inverse ill-posedness do not depend on a cosmological detection.

## Repository structure

```text
code/                    deterministic calculation and audit scripts
observational_forecast/  CLASS provenance and precision settings
source_data/             compact source-data summaries
.github/workflows/       deterministic validation workflows
MANUSCRIPT_PLAN.md        final manuscript integration plan
FINAL_VALIDATION_STATUS_2026-09-13.md  locked numerical anchors
```

## Citation

If you use this code, please cite the associated manuscript. Bibliographic information will be updated when the paper is published.

## License

The code and accompanying source data are released under the MIT License.
