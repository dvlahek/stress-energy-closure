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
- `ksz_wake_screening_summary.json` — kSZ-tagged wake screening control.
- `wake_final_forecast_summary.json`, `wake_robustness_summary.csv` — parity-odd wake forecast and robustness controls.
- `wake_independent_validation_summary.json` — independent Fisher and related validation checks.

## DESI DR1 observational outputs

The final reported null calibration is stored in `phase7_desi_fullsample_perm256_summary.json` and uses 256 frozen-geometry luminosity-mark permutations:

`A_wake = -0.0768409 +/- 0.0851660`, empirical two-sided permutation `p = 0.37354`.

The corresponding 18-component null-corrected vector is `wake_phase7_multitracer_real_vector_perm256.csv`.

For auditability, the originally frozen 32-permutation realization is retained as `phase7_desi_fullsample_summary.json` with vector `wake_phase7_multitracer_real_vector.csv`. It gave `A_wake = -0.0739012 +/- 0.0851661` and `p = 0.39394`. The 256-permutation refinement changes only the number of null shuffles. Tracer definitions, bins, pair geometry, covariance, nuisance model, templates and production seed are unchanged.

Associated validation files:

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

- `lrg_elg_exact_broad_40mock_checkpoint.json` — historical exact broad-bin 40-mock covariance, omnibus, wake-only and linked-standard checkpoint.
- `lrg_elg_zresolved_random_density_audit_10pair.json` — paired nested r1/r4 random-catalog density audit that froze production at `nrandom=4`.
- `lrg_elg_exact_development_checkpoint_2026-09-21.json` — consolidated exact LRGxELG development timeline covering broad 40-mock, z-resolved 40-mock, finite-mock calibration, DESI-fiducial upgrade, random-density decision, regional status and 200-mock guardrails.
- `lrg_elg_ngc_sgc_regional_checkpoint_2026-09-21.json` — cap-specific NGC/SGC robustness checkpoint including nominal fits and finite-mock interpretation.
- `lrg_elg_regional_cap_specific_nuisance_2026-09-21.json` — NGC/SGC-specific magnification, Doppler and evolution-bias nuisance inputs used in the regional linked-standard templates.
- `lrg_elg_ngc_mock0008_health_2026-09-21.json` — technical audit showing the single extreme NGC mock is valid and retained.
- `lrg_elg_mubin_convergence_r4.json` — 120/240/480 angular-discretization convergence audit; production remains frozen at 240 mu bins.

- `lrg_elg_bruteforce_pair_closure_2026-09-21.json` — deterministic independent NumPy-vs-pycorr/Corrfunc closure of pair counts, cross-Landy-Szalay `xi(s,mu)` and odd multipoles; passes at floating-point precision.

## Interpretation guardrails

- CMB, RSD and wake forecast files are sensitivity calculations, not detections.
- The full-sample luminosity-rank conservative DESI fit is the primary observational coefficient.
- The 256-permutation result is a null-calibration refinement of the frozen production estimator, not a new optimization.
- Gfinder is a tracer-proxy robustness control.
- EZmock is a survey-geometry/covariance/placebo layer because the released files used here do not provide the luminosity field required for the physical split.
- AbacusSummit is the physical luminosity-ranked high-fidelity mock validation layer; all 25 final realizations completed.
- The `ell=3` octupole is an orthogonal null control and is not a separate wake measurement.
- No observational result is presented as a neutrino-wake detection.

See `../REPRODUCIBILITY.md` and `../docs/OBSERVATIONAL_REPRODUCIBILITY.md` for the result-to-code map and rerun definitions.
