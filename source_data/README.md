# Source data

This directory contains the compact numerical products used by the manuscript and its reported validation analyses.

Development-stage checkpoints and superseded realizations are stored under `../archive/` and are not part of the standard manuscript reproduction path.

## Theory and forecast products

- `Fig1_prediction_summary.csv` — matched-source nonlinear response example.
- `Fig2_identifiability_summary.csv` — response-based identifiability test.
- `Fig3_tomography_summary.csv`, `Fig3_noise_sweep.csv` — finite-window reconstruction and noise sensitivity.
- `Fig4_mass_sweep.csv` — massive/massless control.
- `Fig5_hierarchy.csv` — finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — direct-versus-memory convergence.
- `CLASS_forecast_validation_summary.json` — CLASS validation summary.
- `CLASS_response_optimization_mass_sweep.csv` — optimized response sweep over relic mass.
- `hidden_channel_retention_direct.json` — direct transfer-level retention control.
- `ksz_wake_screening_summary.json` — kSZ-tagged screening control.
- `wake_final_forecast_summary.json`, `wake_robustness_summary.csv` — parity-odd wake forecast and robustness.
- `wake_independent_validation_summary.json` — independent forecast validation.

## Primary DESI DR1 analysis

The primary reported observational coefficient is the frozen five-tracer luminosity-rank odd-sector analysis:

- `phase7_desi_fullsample_perm256_summary.json`
- `wake_phase7_multitracer_real_vector_perm256.csv`
- `phase7_validation_manifest.json`

Reported result:

[
A_{\rm wake}=-0.0768409\pm0.0851660,
]

with empirical two-sided permutation (p=0.37354).

Supporting validation products retained here because they are part of the manuscript validation chain include:

- `phase7_gfinder_massproxy_summary.json`
- `phase7_seed_control_summary.json`
- `phase7_abacus_final_summary.json`
- `phase7_abacus_injection_amplitude_sweep.json`
- `phase7_ezmock_local/`
- `phase7_octupole_control/`
- `wake_phase7_desi_dr1_zresolved_summary.json`

The earlier 32-permutation realization is preserved under `../archive/phase7/intermediate/`.

## Secondary exact LRG×ELG consistency analysis

The exact-pair branch is reported as an independent consistency analysis, not as the primary DESI coefficient.

Final paper-facing products:

- `lrg_elg_exact_final_manifest_2026-09-23.json` — machine-readable final analysis definition and result hierarchy.
- `lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json` — compact final inference summary.
- `lrg_elg_windowed_forward_r4/` — survey-window matrix, final templates, forward-model summary and full 120-mock leave-one-out audit.
- `lrg_elg_r4_inference_inputs/` — frozen data vector and theory inputs.
- `lrg_elg_covariance_r4_120_checkpoint_2026-09-23.json` — final covariance diagnostics.

Validation products retained in the paper-facing tree:

- `lrg_elg_bruteforce_pair_closure_2026-09-21.json` — independent NumPy versus pycorr/Corrfunc estimator closure.
- `lrg_elg_mubin_convergence_r4.json` — angular-discretization convergence.
- `lrg_elg_zresolved_random_density_audit_10pair.json` — random-catalog density validation.

Final result:

[
A_{\rm wake}=0.00327466\pm0.00148162,
qquad
Z_{\rm nominal}=2.21019.
]

The 120-mock leave-one-out calibration gives (1.84\sigma) from the absolute matched-filter score and (1.96\sigma) from the Sellentin–Heavens likelihood-ratio tail. The result is described as an approximately (2\sigma) hint, not as a detection.

Earlier 40-mock, regional, pre-window and intermediate checkpoints are preserved under `../archive/desi_exact/`.

## Reproducibility

See `../REPRODUCIBILITY.md` for the result-to-code map and `../docs/OBSERVATIONAL_REPRODUCIBILITY.md` for the observational analysis definitions.
