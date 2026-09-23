# Development archive

This directory preserves intermediate calculations and superseded analysis products for provenance.

Files under `archive/` are **not part of the reported manuscript inference** and are not required for the standard reproduction path. They are retained to document the development history of the analysis.

The paper-facing repository is defined by the root `README.md`, `REPRODUCIBILITY.md`, `docs/OBSERVATIONAL_REPRODUCIBILITY.md`, and the final products under `source_data/`.

- `desi_exact/intermediate/` contains earlier DESI LRG×ELG checkpoints, regional diagnostics, pre-window fits, and analysis-freeze records.
- `desi_exact/notes/` contains working notes superseded by the final exact-pair documentation.
- `desi_exact/code/` and `desi_exact/scripts/` contain superseded broad-bin, regional and development-only implementations.
- `phase7/intermediate/` contains the earlier 32-permutation realization superseded by the final 256-permutation result.
- `phase7/ezmock_realizations/` contains per-realization mock outputs; the paper-facing aggregate remains under `source_data/phase7_ezmock_local/aggregate/`.
- `theory_development/code/` contains diagnostic scripts used during method development but not required by the final manuscript pipeline.
- `development/HISTORY.md` preserves the chronological development log.

Git history remains the authoritative record of all changes.
