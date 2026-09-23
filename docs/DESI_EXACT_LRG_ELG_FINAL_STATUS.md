# DESI DR1 exact LRG×ELG odd-sector consistency analysis

Date: 2026-09-23  
Status: **final manuscript analysis**

## Role in the manuscript

This exact-pair LRG×ELG branch is a **secondary observational consistency analysis**. It does not replace the repository's primary DESI DR1 luminosity-rank `phase7` inference.

The purpose of this branch is narrower: test a pre-specified parity-odd wake shape with an exact cross-Landy–Szalay estimator, an identically processed mock ensemble, and an explicit survey-window forward model.

No result in this branch is presented as a detection.

## Frozen estimator

- tracers: DESI DR1 LRG × ELG
- orientation: LRG→ELG
- estimator: cross Landy–Szalay in `(s,mu)`
- line of sight: midpoint
- primary multipole: `ell=1`
- control multipole: `ell=3`
- redshift bins: `0.80–0.90`, `0.90–1.00`, `1.00–1.10`
- separation bins: `20–40`, `40–60`, ..., `120–140 Mpc/h`
- angular bins: 240
- angular exclusion: `theta >= 0.05 deg`
- random density: four DESI random realizations per cap/tracer
- weights: `WEIGHT * WEIGHT_FKP`
- cosmology for distances: `cosmoprimo.fiducial.TabulatedDESI`
- primary vector: 18-dimensional dipole

The reported inference uses 120 identically processed mock realizations.

## Validation chain

1. **Estimator validation**
   - independent NumPy brute-force pair closure agrees with pycorr/Corrfunc to floating-point precision;
   - 240 mu bins pass the 120/240/480 convergence audit;
   - the random-density audit freezes production at `nrandom=4`.

2. **120-mock covariance**
   - 18D dipole covariance is positive definite;
   - condition number: `34.07385`;
   - Hartlap factor: `0.840336`.

3. **Null and nuisance validation**
   - zero-null: `chi2=12.66197` for 18 dof, `p=0.81125`;
   - linked-standard nuisance-only residual: `p=0.77184`.

4. **Survey-window forward model**
   - fine DESI `R_LRG R_ELG(s,mu)` counts measured on the frozen r4 random set;
   - wake pre-window/windowed cosine: `0.99998447`;
   - even-to-odd linear Kaiser leakage RMS: `7.38%` of the physical odd-standard RMS;
   - the explicit RR window therefore changes the template modestly but non-negligibly.

## Final result

Using the **window-convolved wake template** and the **window-convolved total standard nuisance template**:

- wake amplitude:
  `0.00327466 +/- 0.00148162`
- nominal matched-filter significance:
  `Z = 2.21019`
- `Delta chi2 = 4.88493`
- nominal two-sided Gaussian `p = 0.02709`

Finite-mock calibration with 120 leave-one-out mock scores gives:

- matched-filter `|Z|` tail:
  - 7 mocks at least as extreme as the data
  - `p+1 = 0.06612`
  - two-sided equivalent: **1.84 sigma**
- Sellentin–Heavens likelihood-ratio tail:
  - 5 mocks at least as extreme as the data
  - `p+1 = 0.04959`
  - two-sided equivalent: **1.96 sigma**

Explicit survey-window convolution changes the nominal matched-filter value only marginally, and the finite-mock tail counts are unchanged.

## Interpretation

> A pre-specified matched filter applied to the DESI DR1 LRG–ELG odd sector yields a nominal 2.21-sigma excess. Leave-one-out finite-mock calibration gives 1.84 sigma using the matched-filter statistic and 1.96 sigma using the Sellentin–Heavens likelihood-ratio tail. We therefore treat the result as a stable approximately two-sigma hint, not as a detection. Explicit survey-window convolution leaves the inference essentially unchanged.

This is the analysis definition used for the reported manuscript result. Earlier regional, pre-window and finite-ensemble development checkpoints are retained under `archive/desi_exact/`.

## Reproduction map

Core scripts:

- `code/desi_dr1_lrg_elg_exact_zresolved.py`
- `code/build_lrg_elg_zresolved_mock_covariance.py`
- `code/fit_lrg_elg_exact_zresolved.py`
- `code/audit_lrg_elg_zresolved_matched.py`
- `code/build_lrg_elg_zresolved_rr_window.py`
- `code/build_lrg_elg_zresolved_windowed_forward.py`

Relevant runners:

- `scripts/run_desi_dr1_lrg_elg_exact_zresolved_data.sh`
- `scripts/process_desi_ezmock_lrg_elg_zresolved_streaming.sh`
- `scripts/run_desi_r4_200_chunk.sh`
- `scripts/run_desi_lrg_elg_rr_window_r4.sh`
- `scripts/run_desi_lrg_elg_windowed_forward.sh`

Final machine-readable manifest:

- `source_data/lrg_elg_exact_final_manifest_2026-09-23.json`

Final compact inference checkpoint:

- `source_data/lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json`
