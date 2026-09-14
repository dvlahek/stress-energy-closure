# Phase-7 odd-octupole control

This branch adds one exploratory observable control after the frozen Phase-7 dipole analysis. It does not modify the primary publication-facing DESI inference.

## Observable

The existing odd dipole uses

\[
(2\ell+1)P_\ell(\mu)=3\mu,\qquad \ell=1.
\]

The control additionally measures

\[
7P_3(\mu)=\frac{7}{2}(5\mu^3-3\mu),\qquad \ell=3,
\]

from the same DESI DR1 BGS galaxies, luminosity-rank marks, sampled pair geometry, redshift bins, separation bins, random catalogues and jackknife regions.

No neutrino-wake template is fit to the octupole. The octupole is only an orthogonal odd-sector null/consistency observable and is not a fifth candidate headline coefficient.

## Predeclared primary test

The primary statistic is the global Mahalanobis distance of the permutation-null-corrected 18-component octupole vector, calibrated by the same 32 within-stratum mark permutations used for the dipole.

Decision rule fixed before inspecting the octupole:

- empirical p >= 0.05: consistent with the null;
- empirical p < 0.05: investigate estimator/survey systematics before any physical interpretation.

The maximum single-bin |z| is diagnostic only.

## Dipole reproduction gate

The run remeasures the dipole using the identical pair sample and compares the null-corrected vector against `source_data/wake_phase7_multitracer_real_vector.csv`.

The reproduction gate passed exactly:

- maximum absolute difference: `0.0`
- RMS difference: `0.0`
- correlation: `1.0`
- dipole global empirical permutation p-value: `0.7878788`

Therefore the added octupole code path does not alter the frozen Phase-7 dipole observable.

## Covariance and condition-number definition

The octupole jackknife covariance is first formed as the standard delete-one covariance `C_JK` from 30 regions. Its reported raw rank is `18/18`. For inversion, the code adds a small diagonal ridge

\[
r = \max\left(10^{-7}\lambda_{\max}(C_{\rm JK}),\;10^{-6}\,\mathrm{median}[\mathrm{diag}(C_{\rm JK})],\;10^{-14}\right)
\]

and uses `C_reg = C_JK + r I`.

The quoted octupole condition number

`kappa = 440.8145`

is **the condition number of the ridge-regularized covariance `C_reg` used in the Mahalanobis inverse**, not the condition number of the raw jackknife covariance. The raw covariance is full rank, but its condition number was not separately reported by this run. The ridge added for the octupole is `3.373602582674635e-10`.

For comparison, the exactly reproduced dipole follows the same convention and has regularized condition number `932.3921`.

## Final octupole result

The predeclared octupole null test gives:

- raw jackknife covariance rank: `18/18`
- ridge-regularized covariance condition number: `440.815`
- global Mahalanobis statistic: `15.6936`
- empirical permutation p-value: `0.8181818`
- asymptotic chi-square diagnostic p-value: `0.6139214`
- maximum single-bin diagnostic: `|z| = 1.3221`

The largest single-bin fluctuation occurs at `0.2 < z < 0.3` and `s = 110 h^-1 Mpc`, with
`xi3_null_corrected = -0.0132833 +/- 0.0100474` (`z = -1.3221`).

The predeclared decision rule is therefore satisfied:

**STATUS: CONSISTENT_WITH_NULL / PASS_NULL_CONTROL**

This does not constitute a second neutrino-wake measurement. It is an orthogonal angular-control result showing no significant higher odd-multipole structure in the same frozen tracer sample and pair geometry.

## Observational inference hierarchy

The sole publication-facing headline observational coefficient remains the conservative full-sample five-tracer luminosity-rank DESI DR1 dipole fit, `A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`), because it is the direct real-data analysis fixed before the later mock-covariance closure; Gfinder, EZmock, Abacus and the present `ell=3` octupole are validation/control layers, not alternative headline estimates.

## Archived outputs

Compact publication-facing outputs are archived under `source_data/phase7_octupole_control/`:

- `summary_octupole_control_compact.json`
- `data_vector_octupole_control.csv`
- `jackknife_covariance_octupole.csv`
- `REPOSITORY_COMMIT.txt`
- `STOCHASTIC_SEED.txt`
- `SHA256SUMS.txt`

`SHA256SUMS.txt` also records the hashes of the locally generated dipole covariance and dipole/octupole permutation matrices, so an exact rerun can be verified without treating those large diagnostic matrices as publication-facing source data.

The run used repository commit `e24251e68d2f5c669341e0d0af458452af3f7a73` and stochastic seed `20260913`.

## Local rerun

```bash
git fetch origin
git switch observable-octupole-control
git pull --ff-only origin observable-octupole-control
bash scripts/run_phase7_octupole_control_local.sh
```

Outputs are written to `phase7_octupole_control/` by default.
