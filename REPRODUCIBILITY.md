# Reproducibility guide

This file maps the main reported results to the scripts and compact source-data products used to reproduce them.

## Result-to-code map

| Result | Main code | Source data |
| --- | --- | --- |
| Same instantaneous source, different gravitational future | `code/ev_flrw_controls.py` | `source_data/Fig1_prediction_summary.csv` |
| Causal response distinguishes matched kinetic states | `code/injectivity_tomography.py` | `source_data/Fig2_identifiability_summary.csv` |
| Finite-window inversion is ill-conditioned | `code/prediction_transfer.py`, `code/noise_sweep.py` | `source_data/Fig3_tomography_summary.csv`, `source_data/Fig3_noise_sweep.csv` |
| Massive/massless identifiability boundary | `code/mass_sweep.py` | `source_data/Fig4_mass_sweep.csv` |
| Several-species aggregate-response ambiguity | analytic result; no numerical dependency | not applicable |
| Finite source-jet hierarchy | `code/hierarchy_test.py` | `source_data/Fig5_hierarchy.csv` |
| Direct versus memory-form validation | `code/direct_vs_memory.py` | `source_data/ExtendedData_direct_memory_convergence.csv` |
| Optimized tensor CMB response | `code/class_response_optimize.py` | `source_data/CLASS_response_optimization_mass_sweep.csv` |
| Linear relative-velocity retention control | `code/hidden_channel_theta_kprofile.py`, `code/class_state_worker_velocity.py` | `source_data/hidden_channel_retention_direct.json` |
| kSZ-tagged wake screening control | `code/ksz_wake_screening.py`, `code/ksz_relative_velocity_alignment.py` | `source_data/ksz_wake_screening_summary.json` |
| Parity-odd wake forecast | `code/wake_desi_multitracer_fisher.py`, `code/wake_desi_robustness.py` | `source_data/wake_final_forecast_summary.json`, `source_data/wake_robustness_summary.csv` |
| Primary DESI DR1 odd-sector inference | `code/desi_dr1_phase7_multitracer_fullsample.py` | `source_data/phase7_desi_fullsample_perm256_summary.json`, `source_data/wake_phase7_multitracer_real_vector_perm256.csv` |
| Secondary exact DESI DR1 LRG×ELG odd-sector consistency test | `code/desi_dr1_lrg_elg_exact_zresolved.py`, `code/fit_lrg_elg_exact_zresolved.py`, `code/audit_lrg_elg_zresolved_matched.py`, `code/build_lrg_elg_zresolved_rr_window.py`, `code/build_lrg_elg_zresolved_windowed_forward.py` | `source_data/lrg_elg_exact_final_manifest_2026-09-23.json`, `source_data/lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json` |
| Mass-proxy robustness | `code/desi_phase7_gfinder_massproxy.py` | `source_data/phase7_gfinder_massproxy_summary.json` |
| EZmock geometry/covariance placebo | `code/desi_phase7_ezmock_placebo_realization.py`, `code/desi_phase7_ezmock_aggregate.py` | `source_data/phase7_ezmock_local/aggregate_validation_summary.json` |
| Abacus physical-mock validation | `code/desi_phase7_mock_realization.py`, `code/desi_phase7_mock_aggregate.py` | `source_data/phase7_abacus_final_summary.json` |
| Odd `ell=3` control | `code/desi_dr1_phase7_octupole_control.py` | `source_data/phase7_octupole_control/summary_octupole_control_compact.json` |

The several-species result is a closed analytic statement: the complete tensor response determines the summed speed-space measure, while the species decomposition remains non-unique without additional information. It therefore has no separate numerical reproduction step.

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

The final reported DESI DR1 null calibration uses 256 frozen-geometry luminosity-mark permutations and gives

`A_wake = -0.0768409 +/- 0.0851660` (`-0.902 sigma`), empirical two-sided permutation `p = 0.37354`.

To reproduce the refined null calibration on Linux, run

```bash
bash scripts/run_desi_perm256_full_linux.sh
```

Gfinder, EZmock, AbacusSummit, stochastic-seed and `ell=3` results test tracer definition, survey geometry/covariance, physical mocks, pair-sampling stability and higher odd multipoles. They are validation layers, not alternative primary estimates.

The machine-readable hierarchy is stored in `source_data/phase7_validation_manifest.json`. Earlier development realizations are retained under `archive/` and are not part of the standard reproduction path.

### Secondary exact LRG×ELG branch

A separate exact-pair DESI DR1 LRG×ELG analysis is retained as a secondary consistency branch. It uses an exact cross-Landy–Szalay estimator, a frozen 18-dimensional z-resolved dipole, 120 identically processed r4 mocks, and an explicit RR-window forward model including linear Kaiser even-to-odd leakage.

Its **final window-convolved** matched-filter result is `Z=2.21019`. Leave-one-out finite-mock calibration gives `p+1=0.06612` (1.84 sigma two-sided equivalent) from the absolute matched-filter score and `p+1=0.04959` (1.96 sigma) from the Sellentin–Heavens likelihood-ratio tail. The result is interpreted as a stable approximately two-sigma hint, not as a detection.

This branch does **not** replace the primary phase7 DESI coefficient. Its final analysis definition, validation chain and reported result are recorded in `source_data/lrg_elg_exact_final_manifest_2026-09-23.json` and `docs/DESI_EXACT_LRG_ELG_FINAL_STATUS.md`. Development-stage checkpoints are retained under `archive/desi_exact/`.

## Retention-control status

The precision-grade transfer-level control is the reference/CREF deformation. At `nk=96`, the direct `theta_(nu-cdm) P_cb` half-pair RMS changes from `5.2816e-4` to `5.3228e-4` under the moderate precision test, while the wake response is `0.230844`. The corresponding integrated wake-to-linear contrast changes from `437.07` to `433.69`.

A response-independent orthogonal direction gives a qualitative cross-check with integrated contrasts `94.0` and `79.4` at standard and moderate precision. The response-optimized development direction is excluded from quantitative results because it is precision sensitive.

A separate kSZ-tagged screening control adds one line-of-sight velocity weight to the wake kernel. At `z=0.3` and `sigma_v=200 km/s`, the untagged and tagged pair fractions are `0.461688` and `0.499525`, an `8.20%` gain. Across the tested grid, the gain remains below `13.16%`. This is a kernel-level screening control, not an absolute kSZ survey forecast.

## Observational reruns

The frozen observational definitions, validation layers and local checkout instructions are in `docs/OBSERVATIONAL_REPRODUCIBILITY.md`.

Large external survey catalogs are not stored in this repository. Scripts and workflows document the required public inputs and retained compact outputs.
