# Reproducibility guide

This file maps the main reported results to the scripts and compact source-data products used to reproduce them.

## Result-to-code map

| Result | Main code | Source data |
| --- | --- | --- |
| Same instantaneous source, different gravitational future | `code/ev_flrw_controls.py` | `source_data/Fig1_prediction_summary.csv` |
| Causal response distinguishes matched kinetic states | `code/injectivity_tomography.py` | `source_data/Fig2_identifiability_summary.csv` |
| Finite-window inversion is ill-conditioned | `code/prediction_transfer.py`, `code/noise_sweep.py` | `source_data/Fig3_tomography_summary.csv`, `source_data/Fig3_noise_sweep.csv` |
| Massive/massless identifiability boundary | `code/mass_sweep.py` | `source_data/Fig4_mass_sweep.csv` |
| Finite source-jet hierarchy | `code/hierarchy_test.py` | `source_data/Fig5_hierarchy.csv` |
| Direct versus memory-form validation | `code/direct_vs_memory.py` | `source_data/ExtendedData_direct_memory_convergence.csv` |
| Optimized tensor CMB response | `code/class_response_optimize.py` | `source_data/CLASS_response_optimization_mass_sweep.csv` |
| Linear relative-velocity retention control | `code/hidden_channel_theta_kprofile.py`, `code/class_state_worker_velocity.py` | `source_data/hidden_channel_retention_direct.json` |
| Parity-odd wake forecast | `code/wake_desi_multitracer_fisher.py`, `code/wake_desi_robustness.py` | `source_data/wake_final_forecast_summary.json`, `source_data/wake_robustness_summary.csv` |
| Primary DESI DR1 odd-sector inference | `code/desi_dr1_phase7_multitracer_fullsample.py` | `source_data/phase7_desi_fullsample_summary.json`, `source_data/wake_phase7_multitracer_real_vector.csv` |
| Mass-proxy robustness | `code/desi_phase7_gfinder_massproxy.py` | `source_data/phase7_gfinder_massproxy_summary.json` |
| EZmock geometry/covariance placebo | `code/desi_phase7_ezmock_placebo_realization.py`, `code/desi_phase7_ezmock_aggregate.py` | `source_data/phase7_ezmock_local/aggregate_validation_summary.json` |
| Abacus physical-mock validation | `code/desi_phase7_mock_realization.py`, `code/desi_phase7_mock_aggregate.py` | `source_data/phase7_abacus_final_summary.json` |
| Odd `ell=3` control | `code/desi_dr1_phase7_octupole_control.py` | `source_data/phase7_octupole_control/summary_octupole_control_compact.json` |

## Core local calculations

Install the Python dependencies with

```bash
python -m pip install -r requirements.txt
```

Then run the core theory and inversion calculations with

```bash
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

CLASS-based calculations use pinned `class_public` commit
`e85808324f51fc694d12e3ed7439552a3c3f9540`.

## Observational inference hierarchy

The primary DESI DR1 coefficient is

`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`), empirical two-sided permutation `p = 0.39394`.

Gfinder, EZmock, AbacusSummit, stochastic-seed and `ell=3` results test tracer definition, survey geometry/covariance, physical mocks, pair-sampling stability and higher odd multipoles. They are validation layers, not alternative primary estimates.

The machine-readable hierarchy is stored in `source_data/phase7_validation_manifest.json`.

## Retention-control status

The precision-grade transfer-level control is the reference/CREF deformation. At `nk=96`, the direct `theta_(nu-cdm) P_cb` half-pair RMS changes from `5.2816e-4` to `5.3228e-4` under the moderate precision test, while the wake response is `0.230844`. The corresponding integrated wake-to-linear contrast changes from `437.07` to `433.69`.

A response-independent orthogonal direction gives a qualitative cross-check with integrated contrasts `94.0` and `79.4` at standard and moderate precision. The response-optimized development direction is excluded from quantitative results because it is precision sensitive.

## Observational reruns

The frozen observational definitions, validation layers and local checkout instructions are in `docs/OBSERVATIONAL_REPRODUCIBILITY.md`.

Large external survey catalogs are not stored in this repository. Scripts and workflows document the required public inputs and retained compact outputs.
