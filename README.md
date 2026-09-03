# Stress-energy does not close gravitational dynamics

Code and source data accompanying the manuscript **“Stress-energy does not close gravitational dynamics”**.

The calculations in this repository illustrate a simple distinction in collisionless gravity: the stress-energy tensor is the local source in Einstein's equations, but it does not generally contain enough information to replace the underlying kinetic state in the subsequent dynamics.

The repository contains four reproducible calculations used in the manuscript:

- `ev_flrw_controls.py` — exact flat-FLRW Einstein–Vlasov controls, including the massless profile-universality test and the massive same-`N^mu`/same-`T_munu` counterexample.
- `mass_sweep.py` — reconstruction of matched matter pairs across particle mass and the corresponding FLRW geometry separation.
- `hierarchy_test.py` — finite source-jet matching with compact-support bump functions.
- `direct_vs_memory.py` — direct phase-space Vlasov evolution compared with the reduced retarded-memory representation.
- `run_all.py` — convenience script for the three numerical campaigns used for Figs. 2–4.

## Reproducing the calculations

Python 3.10 or newer is recommended.

```bash
python -m pip install -r requirements.txt
```

The individual calculations can then be run from the repository root:

```bash
python code/ev_flrw_controls.py --full
python code/mass_sweep.py
python code/hierarchy_test.py
python code/direct_vs_memory.py --full
```

The scripts are deterministic; no random seed is required. Numerical outputs are written to local output folders created by the scripts. The CSV files in `source_data/` contain the data used for the main manuscript figures.

## Main numerical checks

The publication calculations verify the following points.

1. Two smooth isotropic massive Vlasov states can have the same initial particle current and stress-energy tensor while producing different exact FLRW metric evolutions.
2. In the isotropic massless limit, profile dependence collapses and the matched geometries coincide to numerical precision.
3. An arbitrary finite number of local gravitational source jets can be matched while the next jet remains distinct.
4. Direct phase-space evolution and the eliminated retarded-memory description converge to the same transverse-traceless response.

These calculations support the analytic results in the manuscript; they are not used as a substitute for the proofs.

## Repository structure

```text
code/          deterministic calculation scripts
source_data/   source data for the main figures
```

## Citation

If you use this code, please cite the associated manuscript. Bibliographic information will be updated here when the paper is published.

## License

The code and accompanying source data are released under the MIT License.
