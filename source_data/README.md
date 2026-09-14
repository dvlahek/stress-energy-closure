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

Phase-7 observational validation:
- `wake_phase7_desi_dr1_zresolved_summary.json` — redshift-resolved baseline.
- `wake_phase7_multitracer_real_vector.csv` — frozen 18-component full-sample DESI vector.
- `phase7_desi_fullsample_summary.json` — compact five-tracer luminosity-rank result with run provenance.
- `phase7_gfinder_massproxy_summary.json` — compact physical mass-proxy robustness result with exact reproduction hashes.
- `phase7_validation_manifest.json` — top-level provenance/status manifest.
- `phase7_ezmock_local/aggregate_validation_summary.json` — compact 30-realization EZmock PASS and SHA256 manifest.
- `phase7_ezmock_local/aggregate/summary_ezmock_placebo_covariance.json` — raw final EZmock aggregate JSON supplied by the production run.
- `phase7_ezmock_local/mock_02_*` — smoke/transport validation only; excluded from production covariance.

Important interpretation guardrails:
- CMB/RSD/wake forecast files are sensitivity calculations, not detections.
- Real DESI results are null/estimator validation results.
- EZmock is a random-rank geometry/covariance/systematics placebo ensemble because the released files used here lack luminosity fields needed for the physical split.
- Abacus is the physical luminosity-ranked high-fidelity mock layer.

See `docs/PHASE7_REPRODUCIBILITY.md` and `docs/PHASE7_HISTORY.md`.
