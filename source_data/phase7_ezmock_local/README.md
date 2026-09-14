# Phase-7 local EZmock placebo ensemble

This directory archives compact outputs from the local DESI DR1 BGS EZmock placebo-covariance pipeline.

The public DR1 EZmock BGS products used here do not expose luminosity columns. The tracer split is therefore an equal-count rank constructed from the released `RAN_NUM_0_1` field within the same narrow-redshift and NGC/SGC cells used by the real-data luminosity proxy. These files are a survey-geometry/covariance/systematics control, not a luminosity-matched mock likelihood and not a halo-mass wake constraint. The Abacus layer remains the physical luminosity-ranked mock validation.

## Random-catalog requirement

Every production realization must use its **own released EZmock NGC/SGC clustering random catalogs**. An attempted optimization that reused the validated mock-2 random pair across other realizations was rejected after repeated catastrophic estimator failures (including mocks 6, 18, 33 and 50). Passing fixed-random examples are not retained as production results because the construction is not homogeneous across the released mock selection/fiber-assignment realizations.

The earlier compact shared-angular/redshift-remapping route was also rejected because it produced order-unity estimator artifacts and distorted the forward window. Both shared-random shortcuts have been removed from the production pipeline.

## Validated 30-mock aggregate

Production mocks 3--32 form the retained homogeneous ensemble. All 30 use realization-specific released clustering random catalogs. The 18-dimensional sample covariance has full rank. OAS shrinkage is the primary covariance estimate: shrinkage = 0.3551 and condition number = 10.68, compared with condition number 325.90 for the raw sample covariance. The Hartlap factor for the finite sample covariance is 0.3448.

The mock-calibrated real-data conservative OAS fit gives `A_wake = -0.1148 +/- 0.0820` (`-1.40 sigma`), with chi-square 6.31 for 8 dof and p = 0.612. The two-sided empirical placebo p-value is 0.0968 (2 of 30 placebo amplitudes have larger absolute amplitude; finite-sample corrected p = 3/31). The Hartlap/sample-covariance control gives `A_wake = -0.0969 +/- 0.0938` (`-1.03 sigma`). No statistically significant wake-like excess is present.

Unit-injection recovery is unbiased: mean recovered amplitude = 1.0000, standard deviation = 0.0615, median = 1.0110, with 80% of injections inside the nominal OAS one-sigma width. The forward window is stable across realizations; the median fractional scatter is 0.88% for the wake template and 0.85% for the Doppler template.

`aggregate_validation_summary.json` is the compact publication-facing archive of these results and contains SHA256 hashes for the complete local aggregate outputs.

## Archived realization

- `mock_02_*` — validated local smoke realization with 2x selected random density and realization-specific released EZmock random catalogs. It is retained only as a transport/estimator validation and is excluded from the final covariance ensemble.

## Local-file SHA256 provenance

- `mock_02_summary.json`: `d453bcad27d4385cefcd5f7db3bde01ce82ff106373fa694c328aee9d11bf68c`
- `mock_02_vector.csv`: `68ddd23f03a82c77c8e9ffd2316c30dfac68819f1ea718f7986d8d719d2bf75f`
- `mock_02_window.csv`: `9ff6e47cb8351bc45fcc89b4f6a35137ef5341ace8cd9afc38b747fa646f0fc3`

The complete 30-mock aggregate remains reproducible from the production scripts. Its publication-facing summary is retained here because the EZmock layer is a validation/control and not an absolute luminosity-matched likelihood.
