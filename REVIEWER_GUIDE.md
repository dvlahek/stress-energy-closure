# Reviewer guide

This file is the shortest route through the publication-facing repository.

## 1. Core claims

| Manuscript result | Main code | Publication-facing source data |
| --- | --- | --- |
| Same instantaneous source, different gravitational future | `code/ev_flrw_controls.py` | `source_data/Fig1_prediction_summary.csv` |
| Causal response distinguishes matched kinetic states | `code/injectivity_tomography.py` | `source_data/Fig2_identifiability_summary.csv` |
| Finite-window inversion is ill-conditioned | `code/prediction_transfer.py`, `code/noise_sweep.py` | `source_data/Fig3_tomography_summary.csv`, `source_data/Fig3_noise_sweep.csv` |
| Massive/massless identifiability boundary | `code/mass_sweep.py` | `source_data/Fig4_mass_sweep.csv` |
| Finite source-jet hierarchy | `code/hierarchy_test.py` | `source_data/Fig5_hierarchy.csv` |
| Direct versus memory-form validation | `code/direct_vs_memory.py` | `source_data/ExtendedData_direct_memory_convergence.csv` |
| Optimized tensor CMB response | `code/class_response_optimize.py` | `source_data/CLASS_response_optimization_mass_sweep.csv` |
| Parity-odd wake forecast | `code/wake_desi_multitracer_fisher.py`, `code/wake_desi_robustness.py` | `source_data/wake_final_forecast_summary.json`, `source_data/wake_robustness_summary.csv` |
| Primary DESI DR1 odd-sector inference | `code/desi_dr1_phase7_multitracer_fullsample.py` | `source_data/phase7_desi_fullsample_summary.json`, `source_data/wake_phase7_multitracer_real_vector.csv` |
| Mass-proxy robustness | `code/desi_phase7_gfinder_massproxy.py` | `source_data/phase7_gfinder_massproxy_summary.json` |
| EZmock geometry/covariance placebo | `code/desi_phase7_ezmock_placebo_realization.py`, `code/desi_phase7_ezmock_aggregate.py` | `source_data/phase7_ezmock_local/aggregate_validation_summary.json` |
| Abacus physical-mock validation | `code/desi_phase7_mock_realization.py`, `code/desi_phase7_mock_aggregate.py` | `source_data/phase7_abacus_final_summary.json` |
| Odd `ell=3` control | `code/desi_dr1_phase7_octupole_control.py` | `source_data/phase7_octupole_control/summary_octupole_control_compact.json` |

## 2. Observational inference hierarchy

There is exactly one publication-facing headline observational coefficient:

`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`), empirical two-sided permutation
`p = 0.39394`.

It is the conservative full-sample five-tracer luminosity-rank DESI DR1 fit. The Gfinder, EZmock,
AbacusSummit, stochastic-seed and `ell=3` results test tracer definition, survey geometry/covariance,
physical mocks, pair-Monte-Carlo stability and higher odd multipoles, respectively. They do not replace
the primary coefficient.

The machine-readable version of this hierarchy is
`source_data/phase7_validation_manifest.json`.

## 3. Final mock status

The physical AbacusSummit validation contains **25 completed realizations**, mock IDs `0-24`, for an
18-component data vector. The OAS covariance cross-check gives `-1.779 sigma`, while the raw-sample/
Hartlap control gives `-0.722 sigma`; therefore this finite-mock layer is used as validation, not as the
headline likelihood.

The EZmock layer contains 30 homogeneous realization-specific placebo mocks (IDs `3-32`) and is used
for survey-geometry/covariance/false-positive validation because the released files used here do not
provide the luminosity information needed for the physical tracer split.

## 4. Reproduction

Core calculations can be rerun after

```bash
python -m pip install -r requirements.txt
```

The exact commands and fixed observational choices are in `docs/PHASE7_REPRODUCIBILITY.md`.
CLASS-based calculations use pinned CLASS commit
`e85808324f51fc694d12e3ed7439552a3c3f9540`.
