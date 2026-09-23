# Gravitational response retains kinetic information beyond stress-energy

This repository contains the analysis code and compact source data for the study **“Gravitational response retains kinetic information beyond stress-energy.”**

The central question is simple: two collisionless states can have the same instantaneous stress-energy tensor and therefore the same gravitational source at one time, while retaining different kinetic structure. We test if that hidden structure can be recovered from the subsequent gravitational response and how much of the information survives projection into realistic observables.

## Main results

The theoretical analysis establishes that, for massive isotropic collisionless matter in the stated class, the complete ideal causal response contains information that is not determined by the instantaneous stress-energy tensor alone. The corresponding finite-window inverse problem remains identifiable but becomes ill-conditioned. The massless isotropic limit removes the same radial kinetic encoding.

The observational part contains two DESI DR1 analyses with different roles:

- **Primary DESI analysis:** a frozen five-tracer luminosity-rank odd-sector estimator with 256 permutation realizations. The conservative wake coefficient is
  [
  A_{\rm wake}=-0.07684\pm0.08517,
  ]
  with empirical two-sided permutation (p=0.3735).
- **Secondary exact LRG×ELG consistency analysis:** an exact cross-Landy–Szalay dipole estimator with 120 mock realizations and explicit survey-window convolution. The final nominal matched-filter value is (2.21\sigma). Finite-mock calibration gives (1.84\sigma) from the absolute matched-filter score and (1.96\sigma) from the Sellentin–Heavens likelihood-ratio tail.

The second result is therefore treated as an approximately (2\sigma) consistency hint, not as a detection. No observational result in this repository is presented as evidence for a neutrino-wake detection.

## Reproducibility

Install the Python dependencies with

```bash
python -m pip install -r requirements.txt
```

The central theory calculations can be reproduced with

```bash
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

For the observational analyses and the exact result-to-code mapping, see:

- `REPRODUCIBILITY.md`
- `docs/OBSERVATIONAL_REPRODUCIBILITY.md`
- `docs/DESI_EXACT_LRG_ELG_FINAL_STATUS.md`
- `source_data/README.md`

Machine-readable analysis summaries are provided in:

- `source_data/phase7_validation_manifest.json`
- `source_data/lrg_elg_exact_final_manifest_2026-09-23.json`

## Repository structure

```text
code/                    analysis and validation code
scripts/                 reproducible local workflows
source_data/             final numerical products used by the manuscript
docs/                    methods and reproducibility notes
observational_forecast/  forecast inputs and supporting notes
archive/                 superseded and development-stage material
```

The `archive/` directory is retained for provenance. Its contents are not part of the reported manuscript inference and are not required for the standard reproduction path.

## Citation and license

If you use this code or source data, please cite the associated study. Bibliographic information can be added to `CITATION.cff` when available.

Code and compact source-data products are released under the MIT License.
