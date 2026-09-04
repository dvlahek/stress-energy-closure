# Source data

This directory contains compact source-data summaries used for the manuscript **“Gravitational response distinguishes hidden kinetic states”** and its numerical controls.

The retained legacy filenames preserve compatibility with the original calculation scripts. Their correspondence to the current manuscript is:

- `Fig1_prediction_summary.csv` — current Fig. 1c: peak complex-response, amplitude and phase separations for the stress-energy-matched prediction example.
- `Fig2_identifiability_summary.csv` — current Fig. 1a,b and response summary: matched-moment accuracy, kernel separation and metric-transfer separation.
- `Fig3_tomography_summary.csv` and `Fig3_noise_sweep.csv` — current Fig. 2 synthetic reconstruction and noise sweep.
- `Fig4_mass_sweep.csv` — current Extended Data massless/massive control.
- `Fig5_hierarchy.csv` — current Extended Data finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — current Extended Data direct-versus-memory convergence summary.
- `CLASS_forecast_validation_summary.json` — Planck-anchored CLASS benchmark, stable multipole range, moment matching and run-to-run validation.
- `CLASS_response_optimization_mass_sweep.csv` — response-optimized mass sweep over 0.03–0.60 eV.
- `CLASS_response_optimization_resolution_convergence.json` — resolution-convergence control for the best tested optimized case.

The CLASS calculations are controlled forecasts, not observational likelihood analyses or detectability claims. Full figure-level source data are supplied with the manuscript, while the workflows in `code/` and `.github/workflows/` deterministically regenerate the corresponding spectra and diagnostic outputs.
