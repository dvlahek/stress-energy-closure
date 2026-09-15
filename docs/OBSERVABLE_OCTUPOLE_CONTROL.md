# Phase-7 odd-octupole control

This document records the final `ell=3` odd-multipole control used in the manuscript reproducibility
snapshot. It does not modify the primary DESI dipole inference.

## Observable

The odd dipole uses

\[
(2\ell+1)P_\ell(\mu)=3\mu,\qquad \ell=1.
\]

The control additionally measures

\[
7P_3(\mu)=\frac{7}{2}(5\mu^3-3\mu),\qquad \ell=3,
\]

from the same DESI DR1 BGS galaxies, luminosity-rank marks, sampled pair geometry, redshift bins,
separation bins, random catalogues and jackknife regions.

No neutrino-wake template is fitted to the octupole. It is an orthogonal odd-sector null/consistency
observable and is not a candidate headline coefficient.

## Fixed test

The primary statistic is the global Mahalanobis distance of the permutation-null-corrected
18-component octupole vector, calibrated by the same 32 within-stratum mark permutations used for the
dipole.

The fixed decision rule is:

- empirical `p >= 0.05`: consistent with the null;
- empirical `p < 0.05`: investigate estimator/survey systematics before physical interpretation.

The maximum single-bin `|z|` is diagnostic only.

## Dipole reproduction gate

The same run remeasures the dipole using the identical pair sample and compares the null-corrected
vector against `source_data/wake_phase7_multitracer_real_vector.csv`.

The reproduction gate passes exactly:

- maximum absolute difference: `0.0`;
- RMS difference: `0.0`;
- correlation: `1.0`.

Thus the additional `ell=3` code path leaves the frozen dipole observable unchanged.

## Covariance convention

The octupole jackknife covariance `C_JK` is formed from 30 delete-one regions and has raw rank
`18/18`. For inversion the code uses

\[
C_{\rm reg}=C_{\rm JK}+rI,
\]

with

\[
r=\max\!\left(10^{-7}\lambda_{\max}(C_{\rm JK}),
10^{-6}\,\mathrm{median}[\mathrm{diag}(C_{\rm JK})],10^{-14}\right).
\]

The reported condition number

`kappa = 440.8145`

is `cond(C_reg)`, the condition number of the ridge-regularized covariance used in the Mahalanobis
inverse. It is not the condition number of the raw jackknife covariance. The raw covariance is full
rank; its condition number is not separately reported by this run. The octupole ridge is
`3.373602582674635e-10`.

## Final result

- raw jackknife covariance rank: `18/18`;
- ridge-regularized covariance condition number: `440.815`;
- global Mahalanobis statistic: `15.6936`;
- empirical permutation p-value: `0.8181818`;
- asymptotic chi-square diagnostic p-value: `0.6139214`;
- maximum single-bin diagnostic: `|z| = 1.3221`.

The largest single-bin fluctuation occurs at `0.2 < z < 0.3` and `s = 110 h^-1 Mpc`, with
`xi3_null_corrected = -0.0132833 +/- 0.0100474` (`z = -1.3221`).

**STATUS: CONSISTENT_WITH_NULL / PASS_NULL_CONTROL**

This is an angular-control result only. The sole publication-facing headline observational
coefficient remains the conservative full-sample luminosity-rank DESI dipole result,
`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`).

## Archived outputs

Publication-facing outputs are under `source_data/phase7_octupole_control/`:

- `summary_octupole_control_compact.json`;
- `data_vector_octupole_control.csv`;
- `jackknife_covariance_octupole.csv`;
- `REPOSITORY_COMMIT.txt`;
- `STOCHASTIC_SEED.txt`;
- `SHA256SUMS.txt`.

`SHA256SUMS.txt` also records hashes of the locally generated dipole covariance and permutation
matrices so a rerun can be verified without storing those larger diagnostic matrices in the source-data
snapshot.

## Local rerun

```bash
git fetch origin
git switch nature-physics-submission
git pull --ff-only origin nature-physics-submission
bash scripts/run_phase7_octupole_control_local.sh
```
