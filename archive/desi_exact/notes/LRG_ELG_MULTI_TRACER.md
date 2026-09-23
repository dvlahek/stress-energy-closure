# DESI DR1 genuine multi-tracer test: LRG x ELG

The existing workflow `.github/workflows/desi_dr1_phase7_multitracer_fullsample.yml`
is **not** a genuine cross-population analysis. It uses five luminosity quantiles
inside BGS.

The clean public-DR1 genuine tracer overlap is LRG x ELG at roughly
`0.8 < z < 1.1`. BGS lives at much lower redshift, so it should not be forced
into the same cross-pair estimator.

## Files

- `code/desi_dr1_lrg_elg_odd.py`
- `code/build_lrg_elg_wake_template.py`
- `code/fit_lrg_elg_wake.py`

Use clustering-ready DESI DR1 v1.5 files for:

- `LRG_NGC_clustering.dat.fits`
- `LRG_SGC_clustering.dat.fits`
- `ELG_LOPnotqso_NGC_clustering.dat.fits`
- `ELG_LOPnotqso_SGC_clustering.dat.fits`

and the corresponding `_0_clustering.ran.fits` random catalogs.

## 1. Smoke test

```bash
python code/desi_dr1_lrg_elg_odd.py \
  --lrg-data LRG_NGC_clustering.dat.fits LRG_SGC_clustering.dat.fits \
  --elg-data ELG_LOPnotqso_NGC_clustering.dat.fits ELG_LOPnotqso_SGC_clustering.dat.fits \
  --lrg-random LRG_NGC_0_clustering.ran.fits LRG_SGC_0_clustering.ran.fits \
  --elg-random ELG_LOPnotqso_NGC_0_clustering.ran.fits ELG_LOPnotqso_SGC_0_clustering.ran.fits \
  --outdir lrg_elg_smoke --zmin 0.80 --zmax 1.10 \
  --max-data-per-tracer 200000 --neighbors-per-anchor 32 \
  --jackknife 20 --seed 20260918
```

## 2. Production data vector

Run the same command with:

```text
--outdir lrg_elg_full
--max-data-per-tracer 0
--neighbors-per-anchor 48
--jackknife 30
```

Then read `z_effective_pair_weighted` from
`lrg_elg_full/summary_lrg_elg.json`.

## 3. Build the frozen hidden-state shape

For example, if the measured effective redshift is 0.93:

```bash
python code/build_lrg_elg_wake_template.py \
  --outdir lrg_elg_template \
  --z 0.93 --mass 0.06 --z-match 1100 --frac 0.30 \
  --s-min 20 --s-max 140 --s-step 20
```

This is intentionally a **shape-only** template. It does not import the low-z
BGS HOD amplitude calibration. The free fitted coefficient absorbs the overall
LRG-ELG normalization.

## 4. Fit

```bash
python code/fit_lrg_elg_wake.py \
  --measurement-dir lrg_elg_full \
  --template lrg_elg_template/lrg_elg_shape_template.csv \
  --outdir lrg_elg_fit
```

Return these outputs for interpretation:

- `lrg_elg_full/summary_lrg_elg.json`
- `lrg_elg_full/lrg_elg_odd_multipoles.csv`
- `lrg_elg_fit/fit_summary.json`

## Before promotion into the manuscript

Run:

1. three independent pair-sampling seeds,
2. an LRG->ELG versus ELG->LRG sign-reversal check,
3. the odd-octupole null control,
4. a mock/window closure test.

The conservative wake coefficient in
`fit_conservative_odd_plus_wake` is the quantity to inspect only after these
checks.
