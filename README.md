# Gravitational response recovers kinetic information beyond stress-energy

Code, source-data summaries and reproducibility workflows accompanying the manuscript **“Gravitational response recovers kinetic information beyond stress-energy”**.

The project asks what information gravity retains after collisionless matter is compressed to its instantaneous particle current and stress-energy tensor. The main analytic results separate four levels of information: instantaneous source equivalence, dynamical inequivalence, ideal causal identifiability and finite-data recoverability. For massive isotropic collisionless matter, the complete causal transverse-traceless response at fixed known nonzero mass and wave number uniquely determines the radial kinetic distribution in the stated weighted-decay class. The isotropic massless limit removes this radial encoding exactly. On finite response windows the forward map remains injective but compact, so inversion is unstable.

The numerical calculations validate the exact Einstein–Vlasov constructions, the response theory and the observational controls. The observational forecasts are sensitivity studies, not detections.

## Reproducing the core calculations

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

Core calculations that do not require CLASS can be run from the repository root:

```bash
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

The CLASS-based calculations use the public solver pinned to commit `e85808324f51fc694d12e3ed7439552a3c3f9540`. They can be run locally with a compatible CLASS installation or through the publication workflows in `.github/workflows/`.

## Publication calculations

The repository keeps only the production-facing observational chain on `main`.

- `code/class_observational_forecast.py` reproduces the controlled CLASS tensor realization.
- `code/class_response_optimize.py` performs the smooth matched-moment response optimization and mass sweep.
- `code/rsd_hidden_state_forecast.py` defines the even-parity RSD observable and nuisance projection.
- `code/rsd_hidden_state_forecast_nonlinear.py` performs the final nonlinear-in-distribution RSD validation.
- `code/wake_two_tracer_fisher.py` defines the parity-odd wake response and baseline Fisher machinery.
- `code/wake_desi_bonvin_fisher.py` supplies the relativistic odd-sector survey model and physical calibration.
- `code/wake_desi_multitracer_fisher.py` implements the full covariance-consistent multi-tracer Fisher calculation.
- `code/wake_desi_multitracer_fullgrid.py` scans the validated HOD threshold grid and tests tracer-resolution saturation.
- `code/wake_desi_robustness.py` repeats the final mass, deformation-cap and scale-cut controls.
- `code/wake_desi_survey_design.py` varies survey area and number density after fixing the final nine-tracer hidden state and physical calibration at the baseline survey.

Production GitHub Actions are:

- `class_observational_forecast.yml`
- `class_response_mass_matrix.yml`
- `class_response_resolution_convergence.yml`
- `rsd_hidden_state_forecast.yml`
- `wake_desi_multitracer_fullgrid.yml`
- `wake_desi_robustness.yml`
- `wake_desi_survey_design.yml`

## Final observational controls

The even-parity cosmological calculations show that most of the hidden-state response lies inside ordinary cosmological nuisance directions. The nonlinear RSD control reaches a maximum projected value of `S/N = 0.28792` in the reported 0.60 eV test.

The parity-odd wake is more informative because its leading response samples the relic distribution at resonant momentum. For the final DESI-BGS-like calculation at `m_nu = 0.06 eV`, a 30% pointwise deformation cap and nine disjoint tracer populations give:

| quantity | S/N |
| --- | ---: |
| before odd-sector projection | 3.19098 |
| after wake-amplitude projection | 0.89809 |
| after the full odd nuisance projection | 0.8605447 |
| linear prediction for the final projected result | 0.8607548 |

The maximum matched-moment mismatch is below `6.0e-16`. Increasing tracer resolution from five to seven to nine populations gives projected values `0.8530`, `0.8589` and `0.8605`, showing saturation. At fixed number density, exact area scaling of the final baseline places `S/N = 1` at approximately `1.8905e4 deg^2`.

The reported robustness sweep gives:

| control | values | projected S/N |
| --- | --- | --- |
| pointwise cap | 10%, 20%, 30% | 0.28687, 0.57371, 0.86054 |
| `k_max` [h/Mpc] | 0.05, 0.075, 0.10 | 0.34327, 0.58626, 0.86054 |
| relic mass [eV] | 0.05, 0.06, 0.08, 0.10 | 0.42511, 0.86054, 2.57537, 5.79156 |

For the mass sweep, the physical wake calibration is recomputed for each mass. Those absolute values therefore inherit the normalization assumptions of the adopted literature-calibrated forecast model. They should not be read as a present observational detection claim.

## Repository structure

```text
code/                    deterministic analysis and validation scripts
observational_forecast/  CLASS provenance and diagnostic notes
source_data/             compact result summaries used for reproducibility
.github/workflows/       publication-facing deterministic workflows
```

Development and diagnostic branches preserve the exploratory history. The `main` branch is the compact publication-facing state of the project.

## Determinism and provenance

The calculations are deterministic. Fixed random seeds are used for synthetic-noise realizations and for candidate searches. Production wake and RSD workflows use seed `20260913`. CLASS is fetched at the pinned commit above and the resolved solver commit is written into every workflow artifact.

## Citation

If you use this code, please cite the associated manuscript. Bibliographic information will be updated after publication.

## License

The code and accompanying source-data summaries are released under the MIT License.
