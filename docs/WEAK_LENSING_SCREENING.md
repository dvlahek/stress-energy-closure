# Weak-lensing wake screening

This checkpoint is intentionally **not** a manuscript result. It is a cheap
kinetic-response test used to decide if a full neutrino-wake weak-lensing
forecast is worth building.

The published wake-lensing proposal of Zhu et al. (PRL 116, 141301; arXiv:1412.1660)
uses the velocity-aligned dipole of the convergence field. The present script
keeps only the hidden-state dependence inherited from the same nonrelativistic
Landau-pole wake factor used elsewhere in this repository,

`Y_F ~ |u| F(m |u| / [Tnu0 (1+z)])`.

All halo, lensing-geometry and survey-noise factors are common to `F_FD`, `F_+`
and `F_-` in this screening layer. The useful output is therefore the pair
separation relative to the fiducial FD wake amplitude,

`R_pair = |A(F_+) - A(F_-)| / A(F_FD)`.

If a future/full calculation gives a fiducial wake-lensing detection significance
`S0`, the corresponding hidden-state separation is approximately
`S_hidden = S0 * R_pair` under this common-noise scaling.

## Run

From the repository root:

```bash
python code/wake_lensing_screening.py --outdir weak_lensing_screening
```

This requires only the normal repository NumPy dependencies and runs in seconds.
It reconstructs the stable reference wake direction at `m=0.06 eV`, `z_match=1100`
and the 30% pointwise cap.

If you have the **validated** `best_pair.csv` from the final wake run, prefer it:

```bash
python code/wake_lensing_screening.py \
  --pair-csv /path/to/best_pair.csv \
  --outdir weak_lensing_screening_validated
```

Optionally attach an assumed fiducial wake-lensing S/N only as a rescaling check:

```bash
python code/wake_lensing_screening.py \
  --pair-csv /path/to/best_pair.csv \
  --baseline-fd-sn 3 \
  --outdir weak_lensing_screening_validated
```

## Outputs

- `screening_grid.csv`: redshift/velocity grid with coherent and RMS kinetic
  response ratios.
- `summary.json`: best robust response and the fiducial wake-lensing S/N required
  for 1-sigma and 3-sigma hidden-state separation.
- `pair_used.csv`: exact distributions used by the screening calculation.

`robust_pair_fraction` is the minimum of the coherent-stack and RMS response
fractions. It is deliberately conservative with respect to these two kinetic
summaries.

## Decision rule

Do **not** add weak lensing to the paper from this script alone. The only
question here is if the hidden-state deformation survives the wake-to-lensing
projection strongly enough to justify a full calculation.

If the fraction is large and stable over physically relevant relative velocities,
the next step is a full velocity-aligned convergence-dipole Fisher calculation
with the lensing kernel, halo population, source redshift distribution, shape
noise and covariance. If it is small, stop and leave weak lensing out of the
manuscript.
