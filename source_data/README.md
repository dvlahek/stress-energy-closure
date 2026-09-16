# Source data

This directory contains compact numerical outputs and provenance for the study **“Gravitational response retains kinetic information beyond stress-energy.”**

## Core theory and forecast outputs

- `Fig1_prediction_summary.csv` — stress-energy-matched nonlinear response example.
- `Fig2_identifiability_summary.csv` — matched-moment and response-separation checks.
- `Fig3_tomography_summary.csv`, `Fig3_noise_sweep.csv` — finite-window reconstruction and noise sensitivity.
- `Fig4_mass_sweep.csv` — massive/massless control.
- `Fig5_hierarchy.csv` — finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — direct-versus-memory convergence.
- `CLASS_forecast_validation_summary.json` — CLASS benchmark/validation.
- `CLASS_response_optimization_mass_sweep.csv` — optimized 0.03–0.60 eV tensor-response sweep.
- `hidden_channel_retention_direct.json` — direct `theta_(nu-cdm) P_cb` versus wake retention and precision controls.
- `wake_final_forecast_summary.json`, `wake_robustness_summary.csv` — parity-odd wake forecast and robustness controls.
- `wake_independent_validation_summary.json` — independent Fisher and related validation checks.

## DESI DR1 observational outputs

The primary observational inference is `phase7_desi_fullsample_summary.json`, the conservative full-sample five-tracer luminosity-rank DESI DR1 result:

`A_wake = -0.0739012 +/- 0.0851661`, empirical two-sided permutation `p = 0.39394`.

Associated files:

- `wake_phase7_multitracer_real_vector.csv` — frozen 18-component full-sample DESI odd vector.
- `wake_phase7_desi_dr1_zresolved_summary.json` — redshift-resolved baseline validation.
- `phase7_gfinder_massproxy_summary.json` — independent physical mass-proxy robustness result.
- `phase7_seed_control_summary.json` — production seed plus three independent stochastic closure seeds.
- `phase7_abacus_final_summary.json` — final 25-realization AbacusSummit physical-mock validation.
- `phase7_abacus_injection_amplitude_sweep.json` — compressed-vector forward-window amplitude calibration.
- `phase7_validation_manifest.json` — machine-readable primary/validation hierarchy.
- `phase7_ezmock_local/aggregate_validation_summary.json` — compact 30-realization EZmock placebo summary and hashes.
- `phase7_ezmock_local/aggregate/` — final EZmock aggregate covariance/vector/window products.
- `phase7_ezmock_local/realizations/mock_03` through `mock_32` — retained homogeneous EZmock realization outputs.
- `phase7_octupole_control/` — `ell=3` odd-multipole control outputs.

## Interpretation guardrails

- CMB, RSD and wake forecast files are sensitivity calculations, not detections.
- The full-sample luminosity-rank conservative DESI fit is the primary observational coefficient.
- Gfinder is a tracer-proxy robustness control.
- EZmock is a survey-geometry/covariance/placebo layer because the released files used here do not provide the luminosity field required for the physical split.
- AbacusSummit is the physical luminosity-ranked high-fidelity mock validation layer; all 25 final realizations completed.
- The `ell=3` octupole is an orthogonal null control and is not a separate wake measurement.
- No observational result is presented as a neutrino-wake detection.

See `../REPRODUCIBILITY.md` and `../docs/OBSERVATIONAL_REPRODUCIBILITY.md` for the result-to-code map and rerun definitions.
