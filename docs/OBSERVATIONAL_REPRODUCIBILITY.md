# Observational reproducibility

This document records the frozen observational analysis, validation layers and rerun information used for the reported DESI DR1 results.

## Analysis status

The observational validation suite is complete.

- DESI DR1 BGS full-sample luminosity-rank inference: complete.
- Gfinder mass-proxy robustness test: complete.
- EZmock 30-realization geometry/covariance placebo ensemble: complete.
- AbacusSummit physical luminosity-ranked validation: 25/25 realizations complete.
- production seed plus three independent stochastic closure seeds: complete.
- odd `ell=3` octupole control: complete.

No observational result is presented as a neutrino-wake detection.

## Frozen scientific choices

The hidden-state template is generated at `m_nu = 0.06 eV`, `z_match = 1100`, with a 30% pointwise deformation cap. CLASS-based calculations use pinned `class_public` commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.

The real-data odd vector uses public DESI DR1 BGS clustering catalogs over `0.10 < z < 0.40`, five luminosity-ranked tracer populations, separation bins from `20` to `140 h^-1 Mpc`, angular exclusion `theta < 0.05 deg`, 48 uniformly sampled eligible neighbours per anchor with inverse-probability pair weights, 30 jackknife regions and 32 frozen-geometry null permutations.

The production pair-Monte-Carlo seed is `20260913`. Independent closure seeds are `20260917`, `20260929` and `20261007`.

## Primary observational inference

The primary coefficient is the conservative full-sample five-tracer luminosity-rank DESI DR1 result in `source_data/phase7_desi_fullsample_summary.json`:

`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`), with empirical two-sided permutation `p = 0.39394`.

This coefficient is primary because it is the direct frozen real-data analysis. Gfinder, EZmock, AbacusSummit and `ell=3` are independent validation/control layers, not alternative primary estimates.

The analysis was not fully blinded because real-data outputs were inspected during pipeline development. Final tracer definitions, separation/redshift bins, nuisance model, production seed and independent closure seeds were fixed before the final closure tests.

## Real-data estimator

Main script:

`code/desi_dr1_phase7_multitracer_fullsample.py`

Retained outputs:

- `source_data/phase7_desi_fullsample_summary.json`
- `source_data/wake_phase7_multitracer_real_vector.csv`

The five luminosity ranks are defined independently within narrow redshift cells and NGC/SGC. The 18-component odd data vector is fitted with the frozen wake template and a conservative per-redshift odd nuisance basis.

## Gfinder mass-proxy robustness

Main script: `code/desi_phase7_gfinder_massproxy.py`

Output: `source_data/phase7_gfinder_massproxy_summary.json`

The conservative result is `A_wake = +0.0687839 +/- 0.0827608` (`0.831 sigma`). This layer tests tracer-proxy dependence only.

## EZmock geometry/covariance placebo

Main scripts:

- `code/desi_phase7_ezmock_placebo_realization.py`
- `code/desi_phase7_ezmock_aggregate.py`

Local runner: `scripts/run_ezmock_placebo_local.sh`

The retained ensemble contains 30 homogeneous realization-specific EZmocks, IDs `3-32`. Each mock uses its own released NGC/SGC clustering random catalogs. The released files used here do not provide the luminosity quantity needed for the physical luminosity-ranked tracer construction, so this layer is used for geometry, pair compression, covariance conditioning and false-positive validation.

Final aggregate:

- covariance rank `18/18`;
- OAS condition number `10.678`;
- real-vector cross-check `A_wake = -0.114812 +/- 0.081995`;
- empirical two-sided placebo `p = 0.09677`;
- unit-injection recovery mean `1.0000`.

Compact summaries are under `source_data/phase7_ezmock_local/`.

## AbacusSummit physical-mock validation

Main scripts:

- `code/desi_phase7_mock_realization.py`
- `code/desi_phase7_mock_aggregate.py`

GitHub Actions workflow: `.github/workflows/desi_phase7_abacus_window_covariance.yml`

The final ensemble contains 25 completed realizations, mock IDs `0-24`, for an 18-component data vector. The sample covariance is full rank.

OAS cross-check: `A_wake = -0.105616 +/- 0.059371` (`-1.779 sigma`), empirical two-sided mock `p = 0.07692`.

Raw-sample/Hartlap control: `A_wake = -0.061003 +/- 0.084503` (`-0.722 sigma`).

Because the ensemble contains only 25 mocks for 18 vector components, the covariance-estimator dependence is treated as a finite-mock diagnostic. Abacus validates the physical luminosity-ranked survey window, covariance response and injection recovery; it does not replace the primary DESI coefficient.

The compressed-vector forward-window amplitude sweep over `A = -2, -1, -0.5, 0, 0.5, 1, 2` recovers unit slope to numerical precision. This is a linear-response calibration test, not a separate catalog-level signal injection.

Outputs:

- `source_data/phase7_abacus_final_summary.json`
- `source_data/phase7_abacus_injection_amplitude_sweep.json`

## Stochastic seed closure

The production seed and three independent closure seeds all retain the conservative result within `|z| < 1.02`. Across the four seeds, `mean(A_wake) = -0.0647682`, with seed-to-seed standard deviation `0.0149695`.

Output: `source_data/phase7_seed_control_summary.json`

The conservative nuisance-projected fit is the production inference. The minimal two-template fit is retained only as a diagnostic because it is seed-sensitive.

## Odd `ell=3` control

Main script: `code/desi_dr1_phase7_octupole_control.py`

Local runner: `scripts/run_phase7_octupole_control_local.sh`

The `ell=3` observable uses the same DESI sample and identical sampled pair geometry as the dipole. No wake template is fitted. The dipole is reproduced exactly by the same code path.

The octupole result is

- raw jackknife covariance rank `18/18`;
- global empirical permutation `p = 0.81818`;
- maximum single-bin diagnostic `|z| = 1.3221`.

For inversion the code uses `C_reg = C_JK + r I` with

`r = max(1e-7 lambda_max(C_JK), 1e-6 median(diag(C_JK)), 1e-14)`.

The reported `kappa = 440.8145` is `cond(C_reg)`, the condition number of the ridge-regularized covariance used in the Mahalanobis inverse. It is not the raw jackknife covariance condition number.

Outputs are archived in `source_data/phase7_octupole_control/`.

## Machine-readable hierarchy

`source_data/phase7_validation_manifest.json` is the top-level machine-readable declaration of the primary result and validation roles.

## Local checkout

```bash
git clone https://github.com/dvlahek/stress-energy-closure.git
cd stress-energy-closure
python -m pip install -r requirements.txt
```

Core theory calculations can then be run with the commands in the root `README.md`. Observational workflows requiring public survey catalogs or CLASS use the external inputs specified by the corresponding scripts/workflows; large external catalogs are intentionally not stored in this repository.
