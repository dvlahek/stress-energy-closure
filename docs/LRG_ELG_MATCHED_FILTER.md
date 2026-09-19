# DESI DR1 LRG x ELG exact matched-filter stage

This stage starts only after the exact Corrfunc data vector is stable against
random-catalog density. It does not optimize scale cuts, redshift cuts or
template choices on the observed data.

## 1. Run the identical exact estimator on mocks

Each mock realization must be processed with the same settings as the frozen
data vector:

- LRG x ELG orientation,
- 0.8 < z < 1.1,
- separation edges 20,40,...,140 Mpc/h,
- 240 mu bins,
- midpoint line of sight,
- theta >= 0.05 degree,
- identical object and FKP weight convention,
- identical survey window and random treatment.

Each output directory must contain:

```
lrg_elg_exact_odd_multipoles.csv
```

No data-dependent scale selection is allowed between the data and mock runs.

## 2. Build the covariance

For mock directories such as `mocks/exact_0000`, `mocks/exact_0001`, ...:

```bash
python code/build_lrg_elg_mock_covariance.py \
  --glob 'mocks/exact_*' \
  --out lrg_elg_mocks/covariance_exact.npz
```

The bundle contains:

- dipole covariance,
- octupole covariance,
- dipole-octupole cross covariance,
- the full joint covariance,
- mock mean vectors,
- finite-mock Hartlap factors and covariance diagnostics.

The covariance is valid only for an identically processed data vector.

## 3. Freeze the wake and nuisance templates

Build the wake shape only after the production effective redshift definition
has been finalized:

```bash
python code/build_lrg_elg_wake_template.py \
  --outdir lrg_elg_template \
  --z <PAIR_WEIGHTED_ZEFF> \
  --mass 0.06 --z-match 1100 --frac 0.30 \
  --s-min 20 --s-max 140 --s-step 20
```

The matched-filter code records SHA256 hashes of the template files. Template
columns therefore have an explicit frozen identity in every result.

The current wake-template builder also provides `doppler_shape`. That column
is useful for a plumbing/provisional matched-filter check. It is not by itself
the final standard odd-sector nuisance model. Final inference must use frozen
templates for the relevant relativistic/Doppler, wide-angle, evolution and
magnification contributions.

## 4. Exact null tests and matched filter

With the r4 data vector:

```bash
python code/fit_lrg_elg_exact_matched.py \
  --measurement lrg_elg_exact_r4/lrg_elg_exact_odd_multipoles.csv \
  --covariance lrg_elg_mocks/covariance_exact.npz \
  --wake-template lrg_elg_template/lrg_elg_shape_template.csv \
  --nuisance-cols doppler_shape \
  --outdir lrg_elg_exact_matched
```

The output `matched_filter_summary.json` reports:

1. dipole null chi-square,
2. octupole null chi-square,
3. joint dipole + octupole null chi-square using the full cross covariance,
4. the profiled fixed-shape wake amplitude and uncertainty,
5. Delta chi-square for adding the wake template,
6. the signed matched-filter Z,
7. the Hartlap correction used for the inverse mock covariance.

The matched-filter significance is a one-template test only when the wake
shape and nuisance basis were fixed independently of the observed data.

## 5. Production interpretation

Do not increase significance by selecting only bins that look favorable in
the observed vector. If a scale window, redshift weighting or tracer split is
optimized, the identical optimization must be repeated on every mock and its
global look-elsewhere probability must be reported.

The preferred sequence is:

1. exact random-density closure,
2. estimator freeze,
3. mock covariance,
4. ordinary dipole/octupole/joint null tests,
5. frozen matched filter,
6. full standard-odd nuisance model,
7. only then any hidden-state interpretation.
