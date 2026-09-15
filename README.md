# Gravitational response recovers kinetic information beyond stress-energy

Code, source-data summaries and reproducibility workflows accompanying the manuscript
**“Gravitational response recovers kinetic information beyond stress-energy.”**

This repository contains the publication-facing analysis and the material required to reproduce the
reported calculations and observational validation.

## Reviewer quick start

The central result separates three questions:

1. can two collisionless states have the same instantaneous gravitational source?
2. can their later gravitational responses differ and identify the hidden kinetic state?
3. how much of that information survives projection into realistic observables?

For massive isotropic collisionless matter, the complete ideal causal transverse-traceless response at
fixed known nonzero mass and wave number identifies the radial kinetic distribution in the stated class.
The isotropic massless limit removes this radial encoding. On finite response windows the forward map
remains injective but compact, so inversion is unstable.

The publication-facing DESI DR1 observational result is a null test. The **sole headline coefficient** is

`A_wake = -0.0739012 +/- 0.0851661` (`-0.868 sigma`),

with empirical two-sided permutation `p = 0.39394`.

Gfinder, EZmock, AbacusSummit, stochastic-seed and `ell=3` octupole results are validation/control
layers and are not alternative headline estimates.

Start with:

- `REVIEWER_GUIDE.md` — claim-to-code/source-data map.
- `docs/PHASE7_REPRODUCIBILITY.md` — final observational definitions and rerun instructions.
- `docs/HIDDEN_CHANNEL_RETENTION.md` — direct same-pair wake versus linear-transfer retention control.
- `source_data/README.md` — publication-facing numerical outputs.
- `source_data/phase7_validation_manifest.json` — machine-readable observational hierarchy.
- `HISTORY.md` — development history for the recent observable-retention tests.

## Reproducing the core theory calculations

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

CLASS-based calculations use `class_public` commit
`e85808324f51fc694d12e3ed7439552a3c3f9540`.

## Observable projection

Within the tested smooth ten-function matched-moment class, the optimized ideal full-sky tensor
B-mode response remains below `S/N = 0.04` for all six tested relic masses from 0.03 to 0.60 eV.
The nonlinear even-parity RSD control reaches projected `S/N = 0.28792` in the reported reference
test. The final DESI-BGS-like nine-tracer parity-odd wake forecast at `m_nu = 0.06 eV` and a 30%
pointwise deformation cap gives projected `S/N = 0.8605447` (`0.8607548` in the corresponding linear
calculation).

A separate transfer-level control compares source-matched kinetic pairs under the direct CLASS
neutrino-CDM velocity-divergence kernel and the resonant wake kernel at `m_nu = 0.06 eV` and `z = 0.3`.
For the reference/CREF deformation with `nk=96`, the density-weighted direct
`theta_(nu-cdm) P_cb` half-pair RMS is `5.2816e-4` at standard precision and `5.3228e-4` after
tightening the integration tolerances, while the corresponding resonant wake half-pair response is
`0.230844`. The integrated wake-to-linear contrast therefore remains about `4.3e2` and changes by less
than one percent under this precision test. The independently reconstructed density-derived
`v_(nu-cdm) P_cb` proxy agrees closely with the direct velocity-transfer result in both precision
settings.

A deterministic response-independent orthogonal null direction provides a separate qualitative
robustness control. Its integrated wake-to-`theta P_cb` contrast is `94.0` at standard precision and
`79.4` at moderate precision. The corresponding linear RMS changes by about `18.5%`, driven mainly by
the high-`k` third, so this second direction is not used as a precision-grade headline coefficient.
It nevertheless remains smooth, has no zero crossings in either precision setting, and preserves a
large wake/linear separation without response optimization.

A response-optimized direction is documented only as development history: its standard-precision
linear response collapses under tighter tolerances and is excluded from quantitative claims.

These response comparisons are transfer-level diagnostics, not absolute bispectrum amplitudes, kSZ/RSD
survey forecasts, or detections. See `docs/HIDDEN_CHANNEL_RETENTION.md`,
`source_data/hidden_channel_retention_direct.json`, and `HISTORY.md`.

## DESI DR1 validation

Primary real-data inference:

- full-sample five-tracer luminosity-rank DESI DR1:
  `A_wake = -0.0739012 +/- 0.0851661`, empirical two-sided `p = 0.39394`.

Independent controls:

- Gfinder mass-proxy tracer-definition check: `+0.06878 +/- 0.08276`;
- 30-realization EZmock geometry/covariance placebo: empirical two-sided `p = 0.0968`;
- 25-realization AbacusSummit physical luminosity-ranked mock validation: OAS cross-check
  `-1.779 sigma`, with raw-sample/Hartlap control `-0.722 sigma`;
- production seed plus three independent fixed closure seeds: conservative `|z| < 1.02`;
- odd `ell=3` control: global empirical permutation `p = 0.81818`, with exact dipole reproduction.

All **25 AbacusSummit realizations** used in the final physical-mock ensemble completed successfully.

No observational result in this repository is presented as a neutrino-wake detection.

## Repository structure

```text
code/                    analysis and validation scripts
scripts/                 local publication rerun helpers
observational_forecast/  CLASS provenance and diagnostic notes
source_data/             publication-facing numerical snapshots and provenance
docs/                    final reproducibility notes
.github/workflows/       publication-facing workflows
```

The default `main` branch is the canonical publication-facing snapshot. The matching
`nature-physics-submission` branch is kept as an explicit submission reference.

## Citation and license

If you use this code, please cite the associated manuscript. Bibliographic information will be updated
after publication. Code and source-data summaries are released under the MIT License.
