# Gravitational response distinguishes hidden kinetic states

Code and source data accompanying the manuscript **“Gravitational response distinguishes hidden kinetic states”**.

The paper studies a distinction between the instantaneous gravitational source and the information contained in dynamical response. Two collisionless kinetic states can have the same particle current and stress-energy tensor while retaining different momentum-space structure. For massive isotropic collisionless matter, at fixed known particle mass and response normalization, the ideal continuous causal transverse-traceless response kernel at any fixed nonzero wave number uniquely identifies the radial kinetic distribution in a weighted-decay class that includes compactly supported and exponentially decaying profiles.

The injectivity argument is an instance of a broader phase-mixing identifiability principle: for responses of the form

`K_F(tau) = integral w(p) F'(p) A(k v(p) tau) dp`,

kinetic information is retained when the known momentum-to-velocity map is one-to-one and the analytic response has non-vanishing even Taylor coefficients, together with the stated weighted-decay and boundary conditions. The massive relativistic velocity map satisfies these conditions. In the isotropic massless limit all momenta propagate at the same speed, the radial profile collapses to the energy-density moment in the response, and injectivity fails.

On every fixed finite momentum and response window the restricted gravitational forward map is compact, so the inverse is necessarily unbounded in its H^1 domain norm even though it is unique. The manuscript therefore distinguishes exact identifiability from stable finite-data recovery. A direct model-level consequence is that stress-energy-matched massive states can produce different frequency-dependent amplitude and phase responses to the same gravitational perturbation.

The repository contains the deterministic calculations used in the manuscript, together with the CLASS campaigns used to separate kinetic response from background and relic-abundance effects.

The main analytic/numerical scripts include:

- `ev_flrw_controls.py` — exact flat-FLRW Einstein–Vlasov controls, including the massless profile-universality test and the massive same-`N^mu`/same-`T_munu` construction.
- `injectivity_tomography.py` — stress-energy-matched response construction and regularized reconstruction of the hidden radial distribution.
- `prediction_transfer.py` — dimensionless model-level amplitude, phase and complex metric-response separation for the stress-energy-matched pair.
- `noise_sweep.py` — finite-data stability test over several response-noise levels and deterministic noise realizations.
- `mass_sweep.py` — matched matter pairs across particle mass and the corresponding FLRW geometry separation.
- `hierarchy_test.py` — finite source-jet matching with compact-support bump functions.
- `direct_vs_memory.py` — direct phase-space Vlasov evolution compared with the reduced retarded-memory representation.
- `class_scalar_response_optimize.py` — scalar CMB/lensing response optimization in the exact stress-energy null space.
- `class_scalar_control_run.py` and `class_scalar_fixed_omega0_run.py` — controlled mass sweeps with fixed background/relic quantities.
- `class_cmb_survey_likelihood_forecast.py` — survey-aware Planck-like, SO and CMB-S4 Gaussian/Fisher forecast with linearized LCDM marginalization.

The CLASS workflows pin the public solver to commit `e85808324f51fc694d12e3ed7439552a3c3f9540`. The publication-facing workflows are:

- `.github/workflows/class_scalar_mass_sweep.yml` — optimized scalar/lensing mass sweep.
- `.github/workflows/class_scalar_control_sweeps.yml` — fixed-total-matter and fixed-relic-rho-at-match controls.
- `.github/workflows/class_scalar_fixed_omega0_control.yml` — strict fixed-present-day relic-abundance control.
- `.github/workflows/class_cmb_survey_likelihood_forecast.yml` — survey-aware Planck/SO/CMB-S4 forecast.
- `.github/workflows/class_response_mass_matrix.yml` — manual final characterization workflow that assembles the 7D hidden-subspace spectrum, 10/20/30% amplitude control and the final four-panel characterization figure.

The current interpretation and manuscript integration plan are recorded in [`MANUSCRIPT_PLAN.md`](MANUSCRIPT_PLAN.md).

## Reproducing the calculations

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

Run the non-CLASS calculations from the repository root:

```bash
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

The CLASS calculations require a local CLASS build or can be reproduced through the supplied GitHub Actions workflows. The calculations are deterministic. Fixed random seeds are used only for the reported synthetic noise realizations and response-optimization search heuristics.

## What is checked numerically

1. Two smooth isotropic massive Vlasov states can have the same initial particle current and stress-energy tensor while producing different exact FLRW metric evolutions.
2. Two massive states with the same `n`, `rho` and `P` have distinct causal TT response kernels.
3. The hidden radial distributions can be reconstructed from finite noisy response data with non-negative regularized inversion.
4. A separate noise sweep shows the practical loss of reconstruction accuracy as response noise increases. This illustrates stability only; injectivity, compactness and the unbounded inverse are proved analytically in the manuscript.
5. In the isotropic massless limit, profile dependence collapses and the matched FLRW geometries coincide to numerical precision.
6. An arbitrary finite number of local gravitational source jets can be matched while the next jet remains distinct.
7. Direct phase-space evolution and the eliminated retarded-memory description converge to the same transverse-traceless response.
8. Scalar CMB, lensing and matter-power calculations retain sensitivity to stress-energy-matched nonthermal relic states.
9. The original optimized `m=0.60 eV` benchmark reaches high-precision ideal full-sky joint TT+TE+EE `S/N = 1.5193`, but this configuration changes both kinematics and the relic gravitational weight.
10. Fixing the total present-day matter density leaves the high-mass response essentially unchanged (`S/N = 1.5316`), excluding a trivial total-matter-density explanation.
11. Fixing the relic energy density at `z=1100` reduces the high-mass response to `S/N = 1.0736`, showing that relic weighting matters but does not remove the effect.
12. The strict control with fixed `omega_ncdm(z=0)` and fixed `omega_cdm` gives a high-precision `m=0.60 eV` joint TT+TE+EE `S/N = 0.2793` and phi-phi `S/N = 0.1748`. About 18.4% of the original ideal CMB signal remains after the present-day relic abundance is fixed exactly.
13. In that strict control, the light/relativistic masses `0.03–0.10 eV` form an approximate response plateau, followed by a systematic increase through `0.18`, `0.30` and `0.60 eV`. This is treated as an empirical massive-regime trend, not as a monotonic-mass theorem.
14. The survey-aware likelihood workflow tests how much of the strict-control signal survives finite sky coverage, beams, detector noise and six-parameter LCDM marginalization.

The numerical calculations illustrate and validate the analytic arguments. The FLRW existence, response injectivity and inverse ill-posedness results do not depend on the numerical examples.

## Repository structure

```text
code/                    deterministic calculation scripts
observational_forecast/  CLASS provenance and precision settings
source_data/             compact source-data summaries
.github/workflows/       deterministic CLASS reproduction workflows
MANUSCRIPT_PLAN.md        working physical interpretation and final manuscript plan
```

Full figure-level source data are supplied with the manuscript submission package and are also produced by the final characterization workflow.

## Citation

If you use this code, please cite the associated manuscript. Bibliographic information will be updated when the paper is published.

## License

The code and accompanying source data are released under the MIT License.
