# DESI DR1 LRG x ELG production odd-multipole analysis

This branch replaces the exploratory nearest-neighbour Monte-Carlo pair
compression with exact Corrfunc pair counts through `pycorr`.

The exploratory run remains useful only as a plumbing/sensitivity check.
It must not be used for the final significance of the LRG x ELG odd sector.

## Scientific status

The current exploratory full-sample run showed:

- raw dipole consistent with zero,
- non-zero octupole structure at the exploratory-covariance level,
- no wake claim.

The production sequence is deliberately staged so the measured odd sector is
not interpreted until pair counting, random noise, covariance and standard
odd physics are controlled.

## 1. Pull and activate the environment

```bash
git pull origin main
source .venv/bin/activate
```

## 2. Install the exact pair-counting backend

`pycorr` currently uses the DESI Corrfunc branch as its CPU pair-counting
engine.

```bash
USE_GPU=0 python -m pip install 'git+https://github.com/cosmodesi/pycorr#egg=pycorr[corrfunc]'
```

Check:

```bash
python - <<'PY'
import pycorr, Corrfunc
print("pycorr", getattr(pycorr, "__version__", "unknown"))
print("Corrfunc", getattr(Corrfunc, "__version__", "unknown"))
PY
```

## 3. Download several DESI random realizations

The previous exploratory analysis used one random realization and then
subsampled it. Production analysis should instead test convergence with
complete independent DESI random realizations.

Start with four:

```bash
bash scripts/download_desi_lrg_elg_randoms.sh 4
```

This verifies/downloads random realizations 0,1,2,3 for LRG and
ELG_LOPnotqso, NGC and SGC.

## 4. Exact one-random baseline

Use the full data catalogs and complete random realization 0:

```bash
python code/desi_dr1_lrg_elg_exact.py \
  --lrg-data \
    data/desi_dr1_lrg_elg/LRG_NGC_clustering.dat.fits \
    data/desi_dr1_lrg_elg/LRG_SGC_clustering.dat.fits \
  --elg-data \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_clustering.dat.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_clustering.dat.fits \
  --lrg-random \
    data/desi_dr1_lrg_elg/LRG_NGC_0_clustering.ran.fits \
    data/desi_dr1_lrg_elg/LRG_SGC_0_clustering.ran.fits \
  --elg-random \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_0_clustering.ran.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_0_clustering.ran.fits \
  --outdir lrg_elg_exact_r1 \
  --zmin 0.80 --zmax 1.10 \
  --mu-bins 240 \
  --theta-min-deg 0.05
```

The script also runs the reversed ELG->LRG estimator. For a correct oriented
odd estimator, both odd multipoles must reverse sign. The summary records the
forward-plus-reverse residual.

## 5. Random-density convergence

Repeat with complete random realizations 0--1 and then 0--3.

For the four-random run, pass all four NGC/SGC files to each random argument.
Example:

```bash
python code/desi_dr1_lrg_elg_exact.py \
  --lrg-data \
    data/desi_dr1_lrg_elg/LRG_NGC_clustering.dat.fits \
    data/desi_dr1_lrg_elg/LRG_SGC_clustering.dat.fits \
  --elg-data \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_clustering.dat.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_clustering.dat.fits \
  --lrg-random \
    data/desi_dr1_lrg_elg/LRG_NGC_0_clustering.ran.fits data/desi_dr1_lrg_elg/LRG_SGC_0_clustering.ran.fits \
    data/desi_dr1_lrg_elg/LRG_NGC_1_clustering.ran.fits data/desi_dr1_lrg_elg/LRG_SGC_1_clustering.ran.fits \
    data/desi_dr1_lrg_elg/LRG_NGC_2_clustering.ran.fits data/desi_dr1_lrg_elg/LRG_SGC_2_clustering.ran.fits \
    data/desi_dr1_lrg_elg/LRG_NGC_3_clustering.ran.fits data/desi_dr1_lrg_elg/LRG_SGC_3_clustering.ran.fits \
  --elg-random \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_0_clustering.ran.fits data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_0_clustering.ran.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_1_clustering.ran.fits data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_1_clustering.ran.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_2_clustering.ran.fits data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_2_clustering.ran.fits \
    data/desi_dr1_lrg_elg/ELG_LOPnotqso_NGC_3_clustering.ran.fits data/desi_dr1_lrg_elg/ELG_LOPnotqso_SGC_3_clustering.ran.fits \
  --outdir lrg_elg_exact_r4 \
  --zmin 0.80 --zmax 1.10 \
  --mu-bins 240 \
  --theta-min-deg 0.05
```

Compare runs with:

```bash
python code/compare_lrg_elg_random_convergence.py \
  --runs lrg_elg_exact_r1 lrg_elg_exact_r2 lrg_elg_exact_r4 \
  --out lrg_elg_random_convergence.json
```

Do not declare convergence from a percentage threshold alone. The final
criterion is that the random-density change be negligible compared with the
statistical uncertainty from the mock covariance.

## 6. Outputs to return before moving to inference

Return:

- `lrg_elg_exact_r1/summary_exact.json`
- `lrg_elg_exact_r1/lrg_elg_exact_odd_multipoles.csv`
- `lrg_elg_exact_r4/summary_exact.json`
- `lrg_elg_exact_r4/lrg_elg_exact_odd_multipoles.csv`
- `lrg_elg_random_convergence.json`

At this stage we compare exact pairs against the exploratory MC vector and
decide how many random realizations are required.

## 7. Covariance and physical model -- next production stage

No significance should be quoted from the exact data vector alone.

The next stage is:

1. run the identical exact estimator on the public DESI DR1 mock ensemble,
2. use the mock covariance as the primary covariance,
3. retain compact-sky jackknife only as a secondary check,
4. fit a frozen standard odd-sector model including relativistic/Doppler,
   evolution and wide-angle contributions,
5. add the hidden-state wake shape as a separate frozen template,
6. require tracer-order sign reversal and odd-multipole consistency before
   any physical claim.

The LRG x ELG analysis tests the same multi-tracer physical principle as the
low-z BGS forecast. It is not the same survey configuration and must not be
described as a direct realization of the BGS forecast.
