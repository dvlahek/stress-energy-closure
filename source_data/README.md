# Source data

This directory contains compact source-data summaries and validation snapshots for
**“Gravitational response recovers kinetic information beyond stress-energy.”**

Core theory/forecast source files retain their historical filenames for compatibility:
- `Fig1_prediction_summary.csv` — stress-energy-matched response example.
- `Fig2_identifiability_summary.csv` — matched-moment and response-separation checks.
- `Fig3_tomography_summary.csv`, `Fig3_noise_sweep.csv` — finite-window reconstruction.
- `Fig4_mass_sweep.csv` — massive/massless control.
- `Fig5_hierarchy.csv` — finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — direct-versus-memory convergence.
- `CLASS_forecast_validation_summary.json` — CLASS benchmark/validation.
- `CLASS_response_optimization_mass_sweep.csv` — optimized 0.03–0.60 eV mass sweep.
- `wake_final_forecast_summary.json` and `wake_robustness_summary.csv` — final wake forecast and sensitivity controls.
- `wake_independent_validation_summary.json` — independent Fisher, published benchmark, proxy-scatter and public null controls.

## Phase-7 observational validation

The publication-facing **primary observational inference** is
`phase7_desi_fullsample_summary.json`: the full-sample five-tracer luminosity-rank conservative
DESI DR1 result. The mock-calibrated coefficients below are validation cross-checks and do not
replace this headline coefficient.

Files:
- `wake_phase7_desi_dr1_zresolved_summary.json` — redshift-resolved baseline.
- `wake_phase7_multitracer_real_vector.csv` — frozen 18-component full-sample DESI vector.
- `phase7_desi_fullsample_summary.json` — **primary** five-tracer luminosity-rank real-data result with run provenance.
- `phase7_gfinder_massproxy_summary.json` — physical mass-proxy tracer-definition robustness result with exact reproduction hashes.
- `phase7_seed_control_summary.json` — production seed plus three predeclared independent pair-MC seeds.
- `phase7_abacus_final_summary.json` — 25/25 Abacus high-fidelity physical-mock covariance/window cross-check; not the headline likelihood.
- `phase7_abacus_injection_amplitude_sweep.json` — compressed-vector forward-window amplitude sweep.
- `phase7_validation_manifest.json` — top-level role/provenance/status manifest and primary-inference declaration.
- `phase7_ezmock_local/aggregate_validation_summary.json` — compact 30-realization EZmock PASS and SHA256 manifest.
- `phase7_ezmock_local/aggregate/summary_ezmock_placebo_covariance.json` — raw final EZmock aggregate JSON supplied by the production run.
- `phase7_ezmock_local/realizations/mock_03` through `mock_32` — compact production EZmock vector/window/summary outputs.
- `phase7_ezmock_local/mock_02_*` — smoke/transport validation only; excluded from production covariance.

Important interpretation guardrails:
- CMB/RSD/wake forecast files are sensitivity calculations, not detections.
- The full-sample luminosity-rank conservative DESI fit is the sole headline observational coefficient.
- Gfinder is a tracer-proxy robustness check.
- EZmock is a random-rank geometry/covariance/systematics placebo ensemble because the released files used here lack luminosity fields needed for the physical split.
- Abacus is the physical luminosity-ranked high-fidelity mock layer, but with 25 mocks for an 18-dimensional vector its OAS/Hartlap covariance sensitivity keeps it in a cross-check role.
- The analysis was not fully blinded; final closure choices and independent seed tests were fixed before closure, but no formal blinded-measurement claim is made.
- No observational result is presented as a neutrino-wake detection.

See `docs/PHASE7_REPRODUCIBILITY.md` and `docs/PHASE7_HISTORY.md`.
