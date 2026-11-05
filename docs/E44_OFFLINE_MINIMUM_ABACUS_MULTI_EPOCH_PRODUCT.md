# E44 — minimum Abacus product that could actually improve the physical closure

E43 showed that more endpoint metadata is not enough. E44 inspected the official AbacusSummit data-product/access documentation and the official `abacusutils` cleaned-catalog loader/cleaning source.

## Source-verified result

The best next data object is **not another raw `halo_info` column**.

AbacusSummit secondary snapshots are explicitly designed to support merger-tree construction from persistent particle IDs. The cleaned halo catalog machinery additionally exposes main-progenitor history fields:

- `N_mainprog`
- `vcirc_max_L2com_mainprog`
- `sigmav3d_L2com_mainprog`
- `haloindex_mainprog`
- `v_L2com_mainprog`

and the cleaned header carries `TimeSliceRedshiftsPrev` / `NumTimeSliceRedshiftsPrev`.

This is substantially closer to the E37/E38 need than the two raw endpoint snapshots because it can provide **multi-epoch same-tree dynamical proxies**.

## What these fields can and cannot do

`N_mainprog` can trace the main-progenitor L1 particle-count history. It remains an L1 count history, not `M200c(t)`.

`vcirc_max_L2com_mainprog` and `sigmav3d_L2com_mainprog` can detect substantial internal dynamical evolution across prior timeslices. They do not uniquely determine the Fourier mass profile `u(k,t)`.

`haloindex_mainprog` and `v_L2com_mainprog` provide immediate-tree identity and previous-snapshot velocity information. They do not supply the full environment/wake history.

Therefore this product may let us **reject an approximately frozen halo-history model** or construct an empirical dynamical timescale proxy. By itself it still cannot certify the E42 history/profile/Born terms.

## Safe acquisition order

1. **Inventory only** whether cleaned/tree auxiliary files already exist locally. No ASDF arrays opened.
2. If exactly one relevant cleaned file exists, perform a **YAML/header-only** probe to confirm `TimeSliceRedshiftsPrev` and descriptor shapes.
3. Only after file size and block structure are known, design a bounded targeted extraction.
4. Do not call `CompaSOHaloCatalog(..., fields='all')` on the local base catalog: the official loader notes main-progenitor loading can be slow, and previous WSL crashes show we must avoid blind decompression.

If the product is absent locally, the next decision should be an exact Globus path/file-size inventory before any transfer. The official cleaned-halo release is self-contained, but bulk redshift aggregates are not an appropriate first download for this project.
