# Phase-7 odd-octupole control

This branch adds one exploratory observable control after the frozen Phase-7 dipole analysis. It does not modify the publication-reproducibility snapshot.

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

No neutrino-wake template is fit to the octupole. The octupole is only an orthogonal odd-sector null/consistency observable.

## Predeclared primary test

The primary statistic is the global Mahalanobis distance of the permutation-null-corrected 18-component octupole vector, calibrated by the same 32 within-stratum mark permutations used for the dipole.

Decision rule fixed before inspecting the octupole:

- empirical p >= 0.05: consistent with the null;
- empirical p < 0.05: investigate estimator/survey systematics before any physical interpretation.

The maximum single-bin |z| is diagnostic only.

## Dipole reproduction gate

The run also remeasures the dipole using the identical pair sample and compares the null-corrected vector against `source_data/wake_phase7_multitracer_real_vector.csv`. This verifies that the new code path reproduces the frozen Phase-7 observable before interpreting the octupole.

## Local run

```bash
git fetch origin
git switch observable-octupole-control
git pull --ff-only origin observable-octupole-control
bash scripts/run_phase7_octupole_control_local.sh
```

Outputs are written to `phase7_octupole_control/` by default. The final compact result is `summary_octupole_control.json`.
