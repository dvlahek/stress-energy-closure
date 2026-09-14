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

This full-sample luminosity-rank fit is retained as the publication-facing **primary observational
inference**. It is the direct real-data analysis that predates the later external mock-covariance
closure and therefore avoids choosing a headline coefficient after comparing mock covariance
estimators. This hierarchy is documented after closure and is not claimed as a formal preregistration.

## 4. Independent physical mass-proxy control

DESI BGS galaxies were matched to the Gfinder `GRP_LOGM` proxy, yielding 271,243 matched objects.
Five weighted mass-proxy tracers gave
`A_wake = +0.06878 +/- 0.08276` under the conservative model.
The luminosity and mass-proxy tests have comparable uncertainty but opposite central signs,
supporting a null interpretation.

A completely separate reproduction produced an identical data vector with SHA256
`a64e522738d5b42ccc1d8cad1cb51767fea3a8d1d019b0d0bf2b0431e49c7ab9`.

The mass-proxy fit is treated as a tracer-definition robustness check, not an alternative headline
coefficient.

## 5. EZmock covariance development

The released EZmock BGS files were verified not to provide the luminosity quantities required for
the physical tracer split. EZmock was therefore redefined as a random-rank placebo/systematics layer.

Two exploratory random-catalog strategies were rejected:
1. synthetic/shared angular random geometry
2. reusing one realization's released random catalogs across different mocks

Both can create unstable or pathological odd vectors. Production was changed to one released
NGC/SGC random-catalog pair per EZmock realization.

The final homogeneous ensemble contains mock IDs 3–32. Its OAS covariance is full rank and
well-conditioned (`kappa = 10.68`). The real-vector covariance cross-check is
`-0.11481 +/- 0.08200`; empirical placebo `p = 0.0968`. Unit-signal injections recover a mean
amplitude equal to one to numerical precision.

Because the released EZmock files do not support the physical luminosity-ranked split, this layer is
used for survey-geometry, covariance-conditioning and false-positive validation, not for the headline
DESI coefficient.

## 6. Abacus and stochastic closure

AbacusSummit was introduced as the physical high-fidelity luminosity-ranked mock layer.
Transport failures were kept distinct from science failures. Final run `34780454743` completed
successfully with all 25 mock IDs `0-24` and produced the final aggregate.

The 18-dimensional Abacus sample covariance is full rank. OAS shrinkage gives `kappa = 29.44` and
the conservative real-vector cross-check
`A_wake = -0.10562 +/- 0.05937` (`-1.779 sigma`), with empirical mock `p = 0.0769`.
The raw sample-covariance/Hartlap control gives
`A_wake = -0.06100 +/- 0.08450` (`-0.722 sigma`).

The Abacus ensemble contains only 25 realizations for 18 vector components (`25/18 = 1.39`). The
OAS/Hartlap difference is therefore treated as a finite-mock covariance diagnostic. Abacus is a
high-fidelity physical-mock cross-check of survey window, covariance behavior, null amplitudes and
injection recovery, not the primary headline likelihood.

A three-seed closure was predeclared to verify that pair-Monte-Carlo sampling does not control the
full-sample conclusion. No seed was selected post hoc. Production seed `20260913` plus independent
seeds `20260917`, `20260929`, and `20261007` all give conservative `|z| < 1.02`; the four-seed
amplitude mean is `-0.06477` with seed scatter `0.01497`.

## 7. Documentation correction and final inference hierarchy

An earlier revision of `docs/PHASE7_REPRODUCIBILITY.md` incorrectly retained the pre-closure
statement that only 24/25 Abacus realizations had completed even though the workflow had subsequently
closed successfully. This was a stale documentation state, not a science-based exclusion. The Git
history preserves the original wording. The corrected final state is 25/25 PASS with the aggregate
frozen in the validation manifest.

The final publication-facing hierarchy is:

1. **Primary real-data inference:** full-sample luminosity-rank DESI conservative fit
   (`-0.07390 +/- 0.08517`, `-0.868 sigma`).
2. **Tracer-proxy robustness:** Gfinder mass-proxy fit (`+0.06878 +/- 0.08276`).
3. **Geometry/systematics placebo:** 30 EZmocks with realization-specific released randoms.
4. **Physical-mock covariance/window validation:** 25/25 AbacusSummit realizations.
5. **Stochastic closure:** production seed plus three predeclared independent seeds.

The analysis was not fully blinded. Real-data outputs were inspected during pipeline development.
Final tracer definitions, redshift/separation bins, nuisance structure, production seed and independent
closure seeds were fixed before the final closure tests.

## 8. Predeclared ell=3 odd-multipole control

After the Phase-7 dipole inference was frozen, an orthogonal ell=3 control was predeclared and run on
the same 300,043-galaxy sample and identical sampled pair geometry. No wake template was fit to ell=3.
The dipole was reproduced exactly (`max |Delta xi1| = 0`, correlation `= 1`).

The octupole gave global empirical permutation `p = 0.81818` and maximum single-bin `|z| = 1.3221`,
so the predeclared null criterion passed.

The ell=3 raw jackknife covariance is full rank (`18/18`). The quoted `kappa = 440.8145` is the
condition number of the ridge-regularized covariance `C_reg = C_JK + r I` used in the Mahalanobis
inverse, not the raw jackknife covariance. The raw condition number was not separately reported.

This octupole result is an orthogonal control only. It is not a fifth candidate for the headline
observational result.

**Single headline rule:** the sole publication-facing headline observational coefficient remains
`A_wake = -0.07390 +/- 0.08517` (`-0.868 sigma`) from the full-sample luminosity-rank DESI fit because
that direct real-data analysis was fixed before later mock-covariance and orthogonal-control closure.

No observational layer is presented as a neutrino-wake detection.
