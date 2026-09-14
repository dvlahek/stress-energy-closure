# Gravitational response recovers kinetic information beyond stress-energy

Code, source-data summaries and reproducibility workflows accompanying the manuscript
**“Gravitational response recovers kinetic information beyond stress-energy”**.

The project asks what information gravity retains after collisionless matter is compressed to its
instantaneous particle current and stress-energy tensor. The analytic results separate source
equivalence, dynamical inequivalence, ideal causal identifiability and finite-data recoverability.
For massive isotropic collisionless matter, the complete ideal causal transverse-traceless response
at fixed known nonzero mass and wave number identifies the radial kinetic distribution in the stated
class. The isotropic massless limit removes this radial encoding exactly. Finite response windows
remain ill-conditioned even when the ideal map is injective.

## Reproducing the core calculations

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

## Final forecast controls

For the final DESI-BGS-like nine-tracer wake calculation at `m_nu = 0.06 eV` and a 30% pointwise
deformation cap, the fully projected result is `S/N = 0.8605447`; the corresponding linear prediction
is `0.8607548`. The controlled CLASS/CMB response optimization remains below `S/N = 0.04` for all
six tested masses from 0.03 to 0.60 eV within the tested smooth ten-function class.

Independent controls include a full-covariance Fisher identity check, a published wake benchmark,
and halo-mass proxy scatter. At 0.2 dex proxy scatter the projected wake value changes only from
`0.86054` to `0.83810`.

## Phase-7 DESI DR1 validation

The observational layer is a validation/null test, not a claimed wake detection.

- Redshift-resolved conservative baseline: `A_wake = 0.03184 +/- 0.22249`.
- Full-sample five-tracer luminosity rank (300,043 galaxies):
  `A_wake = -0.07390 +/- 0.08517`.
- Gfinder physical mass-proxy robustness (271,243 matched galaxies):
  `A_wake = +0.06878 +/- 0.08276`.
- 30-realization EZmock random-rank placebo covariance:
  `A_wake = -0.11481 +/- 0.08200`, empirical two-sided `p = 0.0968`.
  OAS covariance is rank 18/18 with condition number 10.68; unit-injection recovery has mean 1.0.

The EZmock result is explicitly a geometry/covariance/systematics placebo control, not a
luminosity-matched physical covariance. AbacusSummit is the physical high-fidelity luminosity-ranked
mock layer.

The top-level snapshot is `source_data/phase7_validation_manifest.json`.
Exact production choices and rejected exploratory paths are documented in
`docs/PHASE7_REPRODUCIBILITY.md` and `docs/PHASE7_HISTORY.md`.

## Publication calculations

Production-facing scripts include:
- `code/class_observational_forecast.py`
- `code/class_response_optimize.py`
- `code/rsd_hidden_state_forecast.py`
- `code/rsd_hidden_state_forecast_nonlinear.py`
- `code/wake_desi_multitracer_fisher.py`
- `code/wake_desi_robustness.py`
- `code/wake_fisher_independent_check.py`
- `code/wake_published_benchmark.py`
- `code/wake_mass_proxy_scatter.py`
- `code/desi_dr1_phase7_multitracer_fullsample.py`
- `code/desi_phase7_gfinder_massproxy.py`
- `code/desi_phase7_ezmock_placebo_realization.py`
- `code/desi_phase7_ezmock_aggregate.py`
- `code/desi_phase7_mock_realization.py`
- `code/desi_phase7_mock_aggregate.py`

GitHub Actions under `.github/workflows/` reproduce the CLASS, DESI, Gfinder, Abacus and stochastic
closure calculations. The local EZmock production runner is `scripts/run_ezmock_placebo_local.sh`.

## Repository structure

```text
code/                    analysis and validation scripts
scripts/                 local production helpers
observational_forecast/  CLASS provenance and diagnostic notes
source_data/             compact numerical snapshots and provenance
docs/                    Phase-7 reproducibility and analysis history
.github/workflows/       publication-facing workflows
```

Exploratory history is preserved in Git. Production status should be taken from the explicit
source-data manifests and, after release, the cited release/tag rather than inferred from branch names.

## Citation and license

If you use this code, please cite the associated manuscript. Bibliographic information will be
updated after publication. Code and source-data summaries are released under the MIT License.
