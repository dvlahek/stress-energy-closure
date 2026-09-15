# Phase-7 octupole control source data

This directory contains the publication-facing outputs of the fixed DESI DR1 odd-octupole (`ell=3`)
control performed on the frozen five-tracer luminosity-rank sample.

The primary octupole statistic is the global permutation-calibrated Mahalanobis distance. No physical
neutrino-wake template is fitted to `ell=3`.

Final result: `CONSISTENT_WITH_NULL`, with empirical permutation `p = 0.8181818` and maximum
single-bin diagnostic `|z| = 1.3221`.

The dipole reproduction gate passed exactly against the frozen Phase-7 vector: maximum absolute
difference `0`, RMS difference `0`, correlation `1`.

Publication-facing files:

- `summary_octupole_control_compact.json`;
- `data_vector_octupole_control.csv`;
- `jackknife_covariance_octupole.csv`;
- `REPOSITORY_COMMIT.txt`;
- `STOCHASTIC_SEED.txt`;
- `SHA256SUMS.txt`.

`SHA256SUMS.txt` also records hashes for the locally generated dipole covariance and dipole/octupole
permutation matrices. These matrices are reproducible from the recorded code commit and seed and are
not required as primary source-data tables.

Run provenance:

- code commit used for the measurement: `e24251e68d2f5c669341e0d0af458452af3f7a73`;
- stochastic seed: `20260913`.

This is an orthogonal control only. It is not a second neutrino-wake measurement and not an
alternative headline result.
