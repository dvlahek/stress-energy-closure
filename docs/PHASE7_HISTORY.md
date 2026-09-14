# Phase-7 analysis history

This file records the methodological decisions that define the final observational validation chain.
It is intentionally short and focuses on changes that affect scientific interpretation or reproducibility.

## 1. Public DESI estimator smoke test

A first parity-odd test on the public Gfinder catalog gave a clean null (`chi2 = 1.299`, 4 dof,
`p = 0.862`). This established that the basic odd estimator was not pathologically biased, but it
was not treated as a DESI LSS wake likelihood.

## 2. Redshift-resolved DESI DR1 baseline

The first Phase-7 clustering analysis used 23,998 data galaxies and 47,996 randoms.
The conservative per-redshift odd-sector nuisance model gave
`A_wake = 0.03184 +/- 0.22249`, a clean null.

## 3. Full-sample five-tracer luminosity ranking

The analysis was upgraded to all 300,043 selected BGS galaxies and 600,086 randoms.
Five weighted luminosity ranks were defined within narrow redshift and Galactic-cap cells.
The conservative result became
`A_wake = -0.07390 +/- 0.08517`, improving the uncertainty by about 2.6 relative to the
subsampled redshift-resolved baseline.

## 4. Independent physical mass-proxy control

DESI BGS galaxies were matched to the Gfinder `GRP_LOGM` proxy, yielding 271,243 matched objects.
Five weighted mass-proxy tracers gave
`A_wake = +0.06878 +/- 0.08276` under the conservative model.
The luminosity and mass-proxy tests have comparable uncertainty but opposite central signs,
supporting a null interpretation.

A completely separate reproduction produced an identical data vector with SHA256
`a64e522738d5b42ccc1d8cad1cb51767fea3a8d1d019b0d0bf2b0431e49c7ab9`.

## 5. EZmock covariance development

The released EZmock BGS files were verified not to provide the luminosity quantities required for
the physical tracer split. EZmock was therefore redefined as a random-rank placebo/systematics layer.

Two exploratory random-catalog strategies were rejected:
1. synthetic/shared angular random geometry
2. reusing one realization's released random catalogs across different mocks

Both can create unstable or pathological odd vectors. Production was changed to one released
NGC/SGC random-catalog pair per EZmock realization.

The final homogeneous ensemble contains mock IDs 3–32. Its OAS covariance is full rank and
well-conditioned (`kappa = 10.68`). The conservative real-data coefficient is
`-0.11481 +/- 0.08200`; empirical placebo `p = 0.0968`. Unit-signal injections recover a mean
amplitude equal to one to numerical precision.

## 6. Abacus and stochastic closure

AbacusSummit was introduced as the physical high-fidelity luminosity-ranked mock layer.
Transport failures are kept distinct from science failures. At the snapshot recorded on 2026-09-14,
24 of 25 mock science jobs had passed and the remaining job was a download-only retry.

A three-seed closure was predeclared to verify that pair-Monte-Carlo sampling does not control the
full-sample conclusion. No seed is selected post hoc.

The production story is therefore:
real DESI null + physical mass-proxy null + EZmock placebo/covariance PASS + Abacus physical
mock validation, with explicit separation between estimator validation and a physical detection claim.
