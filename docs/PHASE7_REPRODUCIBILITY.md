# Phase-7 observational reproducibility

This document records the production-facing observational validation chain for the manuscript
**“Gravitational response recovers kinetic information beyond stress-energy.”**

## Frozen scientific choices

The hidden-state template is generated at `m_nu = 0.06 eV`, `z_match = 1100`, and a 30% pointwise deformation cap.
All CLASS-based calculations use `class_public` commit
`e85808324f51fc694d12e3ed7439552a3c3f9540`.

The real-data odd vector uses DESI DR1 BGS clustering catalogs over `0.10 < z < 0.40`,
five tracer populations, separation bins from `20` to `140 h^-1 Mpc`, angular exclusion
`theta < 0.05 deg`, 48 sampled eligible neighbours per anchor, 30 jackknife regions,
and 32 frozen-geometry null permutations.

The production seed is `20260913`. The stochastic closure workflow uses three predeclared
independent seeds: `20260917`, `20260929`, and `20261007`.

## Real DESI DR1 analyses

### Full-sample luminosity-rank test

Run `code/desi_dr1_phase7_multitracer_fullsample.py` using the clustering-ready NGC/SGC
BGS data and released random catalogs. The compact result is in
`source_data/phase7_desi_fullsample_summary.json`; the 18-component real-data vector is
`source_data/wake_phase7_multitracer_real_vector.csv`.

The conservative fit is a null result:
`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`).

### Gfinder physical mass-proxy robustness test

Run `code/desi_phase7_gfinder_massproxy.py` through
`.github/workflows/desi_phase7_gfinder_massproxy.yml`.
The exact reproduced data-vector SHA256 is
`a64e522738d5b42ccc1d8cad1cb51767fea3a8d1d019b0d0bf2b0431e49c7ab9`.

The conservative fit is:
`A_wake = +0.0687839 +/- 0.0827608` (`0.831 sigma`).
The minimal ~2.2 sigma coefficient is not stable under the conservative odd-sector nuisance model
and is not treated as evidence for a wake.

## EZmock placebo covariance/systematics ensemble

Released DR1 EZmock BGS files used here do not contain the luminosity fields required for the
physical luminosity split. EZmock is therefore used only as an equal-count random-rank placebo
ensemble testing survey geometry, pair compression, covariance conditioning and false positives.

Every production realization must use its own released clustering random catalogs. Shared-random
experiments are excluded.

Local production command:

```bash
cd ~/stress-energy-closure
ARIA_CONNECTIONS=16 bash scripts/run_ezmock_placebo_local.sh 3 32
```

This produces 30 homogeneous realizations and then runs
`code/desi_phase7_ezmock_aggregate.py`.

Primary covariance is OAS shrinkage. The frozen aggregate gives:
- rank `18/18`
- OAS condition number `10.678`
- `A_wake = -0.114812 +/- 0.081995` (`-1.40 sigma`)
- empirical two-sided placebo `p = 0.09677`
- unit-injection recovery mean `0.99999999999994`

The full raw aggregate JSON is
`source_data/phase7_ezmock_local/aggregate/summary_ezmock_placebo_covariance.json`.
The compact publication-facing check, including SHA256 values for the covariance/vector files,
is `source_data/phase7_ezmock_local/aggregate_validation_summary.json`.

## Abacus high-fidelity validation

AbacusSummit is the physical high-fidelity luminosity-ranked mock layer. The workflow is
`.github/workflows/desi_phase7_abacus_window_covariance.yml`.

At the 2026-09-14 reproducibility snapshot, 24/25 science realizations had completed; the only
unresolved realization was a transport-only download retry for mock 11. The final aggregate is
intentionally not frozen in `phase7_validation_manifest.json` until the closure run completes.

## Production versus exploratory outputs

Production:
- real DESI full-sample luminosity-rank analysis
- Gfinder physical mass-proxy robustness test
- EZmock realizations 3–32 with realization-specific released random catalogs
- Abacus complete cut-sky mocks after final aggregate closure

Non-production:
- EZmock mocks 1–2 as early smoke/transport checks
- fixed/shared-random EZmock experiments
- any output failing the explicit estimator sanity gates
- any significance-selected tracer or mass threshold

No observational result in this repository is presented as a neutrino-wake detection.
