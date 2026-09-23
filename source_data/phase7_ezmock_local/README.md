# Phase-7 EZmock placebo ensemble

This directory contains the publication-facing outputs of the DESI DR1 BGS EZmock
geometry/covariance/placebo validation.

The released EZmock BGS products used here do not expose the luminosity field required for the
physical five-tracer luminosity split. The tracer assignment is therefore an equal-count rank within
the same narrow-redshift and NGC/SGC cells used by the real-data analysis. This ensemble tests survey
geometry, pair compression, covariance conditioning and false positives; it is not a luminosity-matched
physical mock likelihood.

## Final ensemble

The retained homogeneous ensemble contains 30 realizations, mock IDs `3-32`. Every realization uses
its own released NGC/SGC clustering random catalogs.

Final aggregate properties:

- vector dimension: `18`;
- sample covariance rank: `18/18`;
- OAS shrinkage: `0.3551`;
- OAS condition number: `10.68`;
- raw sample-covariance condition number: `325.90`;
- Hartlap factor: `0.3448`.

The mock-calibrated real-vector OAS cross-check gives

`A_wake = -0.1148 +/- 0.0820` (`-1.40 sigma`),

with empirical two-sided placebo `p = 0.0968`. The raw-sample/Hartlap control gives
`-1.03 sigma`.

Unit-injection recovery is unbiased to numerical precision: mean recovered amplitude `1.0000`,
standard deviation `0.0615`. The forward window is stable at the percent level across the retained
realizations.

## Publication-facing files

- `aggregate_validation_summary.json` — compact final validation summary and SHA256 hashes.
- `aggregate/summary_ezmock_placebo_covariance.json` — full aggregate numerical summary.
- `aggregate/ezmock_placebo_covariance_oas.csv` — OAS covariance.
- `aggregate/ezmock_placebo_covariance_sample.csv` — raw sample covariance.
- `aggregate/ezmock_placebo_vectors.csv` — retained realization vectors.
- `aggregate/ezmock_placebo_window_templates.csv` — retained forward-window templates.
- Per-realization outputs are preserved under `../../archive/phase7/ezmock_realizations/` and are not required for the manuscript-facing aggregate.

The complete ensemble can be regenerated with `scripts/run_ezmock_placebo_local.sh`; the retained aggregate and its role in the analysis are documented in `../../docs/OBSERVATIONAL_REPRODUCIBILITY.md`.
