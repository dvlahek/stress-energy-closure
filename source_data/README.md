# Source data

This directory contains compact source data used for the manuscript figures and numerical controls for **“Gravitational response distinguishes hidden kinetic states”**.

The filenames were retained to preserve compatibility with the calculation scripts. Their correspondence to the current manuscript is:

- `Fig1_prediction_summary.csv` — current Fig. 1c: peak complex-response, amplitude and phase separations for the stress-energy-matched prediction example.
- `Fig2_identifiability_summary.csv` — current Fig. 1a,b and response summary: matched-moment accuracy, kernel separation and metric-transfer separation for the response-identifiability example.
- `Fig3_tomography_summary.csv` — current Fig. 2: representative synthetic reconstruction errors, regularization strengths and inverse conditioning.
- `Fig3_noise_sweep.csv` — current Fig. 2d: reconstruction errors for the finite-data response-noise sweep.
- `Fig4_mass_sweep.csv` — current Extended Data massless/massive control: mass-dependent FLRW separation measures.
- `Fig5_hierarchy.csv` — current Extended Data finite source-jet hierarchy.
- `ExtendedData_direct_memory_convergence.csv` — current Extended Data direct-versus-memory convergence summary.

The prediction values are dimensionless proof-of-principle outputs for the validation coupling used in the manuscript and are not observational sensitivity forecasts. Full frequency-series, time-series and profile data are generated deterministically by the scripts in `code/` and are supplied as Source Data with the manuscript.
