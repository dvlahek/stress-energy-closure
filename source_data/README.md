# Source data

This directory contains compact source-data summaries used for the manuscript **“Gravitational response recovers kinetic information beyond stress-energy”** and its numerical controls.

The retained legacy filenames preserve compatibility with the original calculation scripts. Their correspondence to the current manuscript is:

- `Fig1_prediction_summary.csv` — model-level complex-response, amplitude and phase separations for the stress-energy-matched prediction example.
- `Fig2_identifiability_summary.csv` — matched-moment accuracy, kernel separation and metric-transfer separation.
- `Fig3_tomography_summary.csv` and `Fig3_noise_sweep.csv` — finite-window synthetic reconstruction and noise sweep.
- `Fig4_mass_sweep.csv` — massless/massive control.
- `Fig5_hierarchy.csv` — finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — direct-versus-memory convergence summary.
- `CLASS_forecast_validation_summary.json` — Planck-anchored CLASS benchmark and numerical validation.
- `CLASS_response_optimization_mass_sweep.csv` — response-optimized mass sweep over 0.03–0.60 eV.
- `CLASS_response_optimization_resolution_convergence.json` — resolution-convergence control for the best tested optimized case.
- `wake_final_forecast_summary.json` — publication-facing nine-tracer wake and RSD summary, including the fixed-density area scaling.
- `wake_robustness_summary.csv` — final wake mass, pointwise-cap and scale-cut robustness values.
- `wake_independent_validation_summary.json` — independent Fisher identity check, published threshold-shape benchmark, halo-mass proxy-scatter sweep and public DESI DR1 parity-odd null test.
- `phase7_ezmock_local/` — compact locally executed EZmock random-rank placebo outputs and provenance. Mock 2 is retained only as a smoke/transport validation; the homogeneous fixed-full-random production ensemble begins at mock 3.

The CLASS and survey calculations are controlled forecasts and sensitivity studies, not observational detections. The DESI DR1 results are estimator/likelihood validation and null tests, not a claimed wake detection. The EZmock layer is explicitly a geometry/covariance/systematics placebo control because the released mock files used here do not provide the luminosity fields required for a luminosity-matched split. The Abacus layer remains the physical luminosity-ranked mock validation. Full figure-level source data are supplied with the manuscript submission package. The scripts in `code/` and retained publication workflows in `.github/workflows/` regenerate the corresponding spectra and diagnostics.
