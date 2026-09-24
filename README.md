# Kinetic information beyond stress-energy in gravitational response

This repository provides the code, numerical source data and validation results for our Einstein–Vlasov study of gravitational response and kinetic-state identifiability.

Two collisionless distributions can agree in their instantaneous stress-energy tensor while differing in their momentum-space structure. We construct source-matched states and calculate their subsequent gravitational response. We then examine which distinctions survive finite observation windows and the projection into cosmological observables.

## Results

The calculations establish response-level distinctions between the source-matched massive states considered here. Finite-window reconstruction is more poorly conditioned than the ideal causal-response problem. We also evaluate the massless isotropic limit and the dependence on the kinetic hierarchy.

The observational analysis uses DESI DR1 in two complementary tests. The five-tracer luminosity-rank analysis is our primary reported DESI measurement. Its conservative wake coefficient is $A_{\rm wake}=-0.07684\pm0.08517$, with empirical two-sided permutation $p=0.3735$ from 256 permutations.

The separate exact LRG–ELG dipole analysis uses 120 mock realizations and an explicit random-pair survey-window model. It gives a nominal matched-filter value of $2.21\sigma$. Empirical finite-mock calibration gives $1.84\sigma$ from the absolute matched-filter statistic and $1.96\sigma$ from the Sellentin–Heavens likelihood-ratio statistic. We treat this result as a tentative consistency indication, not a detection.

## Reproducing the calculations

Install the Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

The central theory calculations can be run with:

```bash
python code/ev_flrw_controls.py --full
python code/injectivity_tomography.py
python code/prediction_transfer.py
python code/noise_sweep.py
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

The CLASS calculations use `class_public` commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.

The observational estimator definitions, required public inputs and result-to-code mapping are documented in `REPRODUCIBILITY.md` and `docs/OBSERVATIONAL_REPRODUCIBILITY.md`. The exact LRG–ELG analysis is described separately in `docs/DESI_EXACT_LRG_ELG_FINAL_STATUS.md`. The eBOSS DR16 [prospective replication protocol](docs/EBOSS_DR16_PROSPECTIVE_REPLICATION.md) records catalogue and mock metadata only; no eBOSS odd-sector measurement is reported. Numerical products are indexed in `source_data/README.md`.

## Directory structure

- `code/`: numerical models, estimators and validation.
- `scripts/`: local survey-data and mock workflows.
- `source_data/`: numerical products underlying the reported results.
- `docs/`: analysis definitions and reproducibility notes.
- `observational_forecast/`: forecast inputs and supporting calculations.

Earlier implementations and intermediate results are retained on the [research archive branch](https://github.com/dvlahek/stress-energy-closure/tree/archive/research-development-2026-09-24). The archive is separate from the manuscript reproduction path.

## Citation and license

Citation metadata is provided in `CITATION.cff`. The code is distributed under the MIT License.
