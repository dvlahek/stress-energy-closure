# Phase-7 observational reproducibility

This document records the production-facing observational validation chain for the manuscript
**“Gravitational response recovers kinetic information beyond stress-energy.”**

## Closure status and documentation correction

Phase-7 observational validation is closed. All 25 AbacusSummit science realizations completed,
the final Abacus aggregate was produced successfully, and the three predeclared stochastic seed
controls completed successfully.

An earlier revision of this document retained the pre-closure statement that only 24/25 Abacus
realizations had completed. That wording was stale relative to the completed workflow and was
inconsistent with `source_data/phase7_validation_manifest.json`. The Git history preserves that
revision. The current document supersedes it: the final state is **25/25 Abacus realizations PASS**
and the aggregate is frozen in the publication-facing source-data snapshot. No scientific result
was changed by this documentation correction.

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

## Primary observational inference

The publication-facing **primary real-data inference** is the conservative full-sample
five-tracer luminosity-rank DESI DR1 fit in
`source_data/phase7_desi_fullsample_summary.json`:

`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`),
with conservative empirical two-sided permutation `p = 0.39394`.

**This is the sole headline observational coefficient because it is the direct real-data analysis fixed before the later mock-covariance closure; Gfinder, EZmock, Abacus, and the ell=3 octupole are validation/control layers rather than alternative headline estimates.**

This hierarchy avoids selecting a headline coefficient after comparing several covariance estimators.
It is documented after closure and is **not claimed as a formal preregistration**.

The analysis was not fully blinded. Real-data outputs were inspected during pipeline development.
However, the final tracer definition, separation/redshift bins, nuisance model, production pair-MC
seed, and independent closure seeds were fixed before the final closure tests. No result is described
as a blinded measurement.

## Real DESI DR1 analyses

### Full-sample luminosity-rank test — primary

Run `code/desi_dr1_phase7_multitracer_fullsample.py` using the clustering-ready NGC/SGC
BGS data and released random catalogs. The compact result is in
`source_data/phase7_desi_fullsample_summary.json`; the 18-component real-data vector is
`source_data/wake_phase7_multitracer_real_vector.csv`.

The conservative fit is the primary null result:
`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`).

### Gfinder physical mass-proxy robustness test — proxy cross-check

Run `code/desi_phase7_gfinder_massproxy.py` through
`.github/workflows/desi_phase7_gfinder_massproxy.yml`.
The exact reproduced data-vector SHA256 is
`a64e522738d5b42ccc1d8cad1cb51767fea3a8d1d019b0d0bf2b0431e49c7ab9`.

The conservative fit is:
`A_wake = +0.0687839 +/- 0.0827608` (`0.831 sigma`).
The minimal ~2.2 sigma coefficient is not stable under the conservative odd-sector nuisance model
and is not treated as evidence for a wake. The mass-proxy result is a tracer-definition robustness
check, not an alternative headline coefficient.

## EZmock placebo covariance/systematics ensemble — geometry/systematics cross-check

Released DR1 EZmock BGS files used here do not contain the luminosity fields required for the
physical luminosity split. EZmock is therefore used only as an equal-count random-rank placebo
ensemble testing survey geometry, pair compression, covariance conditioning and false positives.

Every production realization uses its own released clustering random catalogs. Shared-random
experiments are excluded.

The final 30-realization ensemble gives:
- rank `18/18`
- OAS condition number `10.678`
- real-vector cross-check `A_wake = -0.114812 +/- 0.081995` (`-1.40 sigma`)
- empirical two-sided placebo `p = 0.09677`
- unit-injection recovery mean `0.99999999999994`

Because the released EZmock catalogs do not support the physical luminosity-ranked tracer
construction, this coefficient is not used as the publication headline inference.

## Abacus high-fidelity validation — physical-mock cross-check

AbacusSummit is the physical high-fidelity luminosity-ranked mock layer. Final run `34780454743`
completed successfully with all **25/25** mock IDs `0-24` and a full-rank 18-dimensional sample
covariance.

The OAS covariance has shrinkage `0.22540` and condition number `29.437`. Using this covariance,
the conservative real-vector cross-check is

`A_wake = -0.105616 +/- 0.059371` (`-1.779 sigma`),
with empirical two-sided mock `p = 0.07692`.

This is **not the primary observational coefficient**. The ensemble contains only 25 realizations
for an 18-dimensional vector (`25/18 = 1.39`), so covariance-estimator uncertainty is non-negligible.
The raw sample-covariance/Hartlap control gives

`A_wake = -0.061003 +/- 0.084503` (`-0.722 sigma`).

The OAS/Hartlap difference is treated as a finite-mock covariance diagnostic. Abacus is retained as
an end-to-end physical-mock validation of survey window, null behavior and signal recovery, not as
the headline likelihood.

## Stochastic seed closure

Three independent seeds (`20260917`, `20260929`, `20261007`) were declared before the closure run.
Together with production seed `20260913`, the conservative luminosity-rank coefficient remains
within `|z| < 1.02`. The mean amplitude over the four seeds is `-0.06477` with seed-to-seed standard
deviation `0.01497`. No seed is selected post hoc.

The minimal two-template coefficient is seed-sensitive and is retained only as a diagnostic. The
conservative nuisance-projected fit is the production inference.

## Orthogonal ell=3 odd-multipole control

The predeclared ell=3 control uses the same frozen DESI sample and sampled pair geometry as the
dipole and fits no wake template. Its global permutation-calibrated result is `p = 0.81818`, with
maximum single-bin `|z| = 1.3221`, so it is consistent with the null.

The raw ell=3 jackknife covariance is full rank (`18/18`). For inversion the code uses
`C_reg = C_JK + r I`, with
`r = max(1e-7 lambda_max(C_JK), 1e-6 median(diag(C_JK)), 1e-14)`.
The reported `kappa = 440.8145` is `cond(C_reg)`, the condition number of this ridge-regularized
covariance used in the Mahalanobis inverse. It is **not** the condition number of the raw jackknife
covariance; the raw condition number was not separately reported.

This ell=3 result is an orthogonal control, not a fifth candidate headline measurement.

## Production versus control outputs

Production headline:
- full-sample luminosity-rank DESI conservative fit — sole primary observational inference

Validation/control layers:
- Gfinder physical mass-proxy robustness test
- EZmock 30-realization geometry/covariance/placebo ensemble
- Abacus 25/25 cut-sky physical mocks and final aggregate
- three-seed stochastic closure
- predeclared ell=3 odd-octupole null control

Excluded exploratory outputs:
- EZmock mocks 1–2 as early smoke/transport checks
- fixed/shared-random EZmock experiments
- any output failing explicit estimator sanity gates
- any significance-selected tracer or mass threshold

No observational result in this repository is presented as a neutrino-wake detection.
