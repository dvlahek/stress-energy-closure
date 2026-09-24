#!/usr/bin/env python3
"""Audit the public eBOSS DR16 LRG and ELG random-catalogue selection.

This is a pre-unblinding selection diagnostic. It downloads the four released
random FITS catalogues, verifies their hashes/columns and counts their
redshift-dependent sky occupancy. It never forms pairs or reads an odd vector.
A gridded random intersection is not an exact survey mask or window function.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from astropy.io import fits

from inspect_eboss_dr16_selection import BASE, WEIGHTS, Z_EDGES, download

RANDOMS = {
    ("LRG", "NGC"): ("eBOSS_LRG_clustering_random-NGC-vDR16.fits", 5460719),
    ("LRG", "SGC"): ("eBOSS_LRG_clustering_random-SGC-vDR16.fits", 3453453),
    ("ELG", "NGC"): ("eBOSS_ELG_clustering_random-NGC-vDR16.fits", 3728363),
    ("ELG", "SGC"): ("eBOSS_ELG_clustering_random-SGC-vDR16.fits", 3609460),
}
CANDIDATE_LO, CANDIDATE_HI = 0.6, 1.0
# Equal-area coordinate uses sin(dec). Pixel areas are approximately uniform
# in the (RA, sin(dec)) grid, but pixel boundaries do not define the survey mask.
RA_STEP = 0.5
SIN_DEC_STEP = 0.01
N_RA = int(round(360 / RA_STEP))
N_SIN = int(round(2 / SIN_DEC_STEP)) + 1
N_PIX = N_RA * N_SIN
MIN_RANDOMS_PER_PIXEL = 25


def grid_indices(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ix = np.floor(ra / RA_STEP).astype(np.int64)
    iy = np.floor((np.sin(np.deg2rad(dec)) + 1.0) / SIN_DEC_STEP).astype(np.int64)
    return np.clip(ix, 0, N_RA - 1) * N_SIN + np.clip(iy, 0, N_SIN - 1)


def examine(path: Path, tracer: str, cap: str, expected: int, sha: str,
            size: int, chunk_rows: int) -> tuple[dict, np.ndarray, list[np.ndarray]]:
    cell_counts = np.zeros(N_PIX, dtype=np.int64)
    per_bin_cells = [np.zeros(N_PIX, dtype=np.int64) for _ in range(len(Z_EDGES) - 1)]
    zcounts = np.zeros(len(Z_EDGES) - 1, dtype=np.int64)
    wdiag = {name: {"nonfinite": 0, "nonpositive": 0,
                    "minimum_finite": None, "maximum_finite": None}
             for name in WEIGHTS}
    invalid_coordinates = 0
    candidate_count = 0
    valid_count = 0

    with fits.open(path, memmap=True) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError(f"Expected one binary table in {path.name}")
        table = tables[0]
        names = set(table.columns.names)
        needed = {"RA", "DEC", "Z", *WEIGHTS}
        if not needed.issubset(names):
            raise ValueError(f"Missing random columns: {path.name}: {sorted(needed - names)}")
        total = int(table.header["NAXIS2"])
        if total != expected:
            raise ValueError(f"Random FITS row count mismatch: {path.name}: {total} != {expected}")
        for start in range(0, total, chunk_rows):
            stop = min(total, start + chunk_rows)
            batch = table.data[start:stop]
            ra = np.asarray(batch["RA"], dtype=np.float64)
            dec = np.asarray(batch["DEC"], dtype=np.float64)
            z = np.asarray(batch["Z"], dtype=np.float64)
            valid = (np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                     & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
            invalid_coordinates += int(np.count_nonzero(~valid))
            valid_count += int(np.count_nonzero(valid))
            for name in WEIGHTS:
                w = np.asarray(batch[name], dtype=np.float64)
                finite = np.isfinite(w)
                d = wdiag[name]
                d["nonfinite"] += int(np.count_nonzero(~finite))
                d["nonpositive"] += int(np.count_nonzero(finite & (w <= 0)))
                if finite.any():
                    wmin = float(np.min(w[finite]))
                    wmax = float(np.max(w[finite]))
                    d["minimum_finite"] = (wmin if d["minimum_finite"] is None
                                            else min(d["minimum_finite"], wmin))
                    d["maximum_finite"] = (wmax if d["maximum_finite"] is None
                                            else max(d["maximum_finite"], wmax))
            candidate = valid & (z >= CANDIDATE_LO) & (z < CANDIDATE_HI)
            candidate_count += int(np.count_nonzero(candidate))
            if candidate.any():
                idx = grid_indices(ra[candidate], dec[candidate])
                cell_counts += np.bincount(idx, minlength=N_PIX)
            for i, (lo, hi) in enumerate(zip(Z_EDGES[:-1], Z_EDGES[1:])):
                selected = valid & (z >= lo) & (z < hi)
                zcounts[i] += int(np.count_nonzero(selected))
                if selected.any():
                    idx = grid_indices(ra[selected], dec[selected])
                    per_bin_cells[i] += np.bincount(idx, minlength=N_PIX)

    return {
        "tracer": tracer, "cap": cap, "filename": path.name,
        "file_sha256": sha, "downloaded_bytes": size, "header_rows": expected,
        "valid_coordinate_z_rows": valid_count,
        "invalid_coordinate_z_rows": invalid_coordinates,
        "candidate_z_interval": [CANDIDATE_LO, CANDIDATE_HI],
        "candidate_random_rows": candidate_count,
        "redshift_edges": list(Z_EDGES), "redshift_bin_counts": zcounts.tolist(),
        "weight_column_diagnostics": wdiag,
        "candidate_occupied_pixels_ge1": int(np.count_nonzero(cell_counts)),
        "candidate_supported_pixels_ge_threshold": int(np.count_nonzero(
            cell_counts >= MIN_RANDOMS_PER_PIXEL)),
    }, cell_counts, per_bin_cells


def overlap(a: np.ndarray, b: np.ndarray, threshold: int) -> dict:
    sa = a >= threshold
    sb = b >= threshold
    both = sa & sb
    n_a, n_b, n_both = map(lambda x: int(np.count_nonzero(x)), (sa, sb, both))
    union = int(np.count_nonzero(sa | sb))
    return {
        "minimum_randoms_per_tracer_per_pixel": threshold,
        "lrg_supported_pixels": n_a, "elg_supported_pixels": n_b,
        "joint_supported_pixels": n_both, "union_supported_pixels": union,
        "joint_fraction_of_lrg_pixels": n_both / n_a if n_a else None,
        "joint_fraction_of_elg_pixels": n_both / n_b if n_b else None,
        "pixel_area_square_deg_approx": RA_STEP * SIN_DEC_STEP * (180 / np.pi) ** 2,
        "joint_pixel_area_square_deg_approx": (
            n_both * RA_STEP * SIN_DEC_STEP * (180 / np.pi) ** 2),
        "interpretation": (
            "Diagnostic intersection of tracer random supports on the stated grid; "
            "not an exact joint mask, pair-window estimate or effective survey area."
        ),
    }


def self_test() -> None:
    ra = np.array([0.0, 0.49, 0.5, 359.99])
    dec = np.array([0.0, 0.0, 0.0, 90.0])
    indices = grid_indices(ra, dec)
    assert len(set(indices.tolist())) == 3
    assert np.all((indices >= 0) & (indices < N_PIX))
    a = np.array([0, 25, 30, 0], dtype=np.int64)
    b = np.array([0, 26, 0, 25], dtype=np.int64)
    x = overlap(a, b, 25)
    assert x["lrg_supported_pixels"] == 2
    assert x["elg_supported_pixels"] == 2
    assert x["joint_supported_pixels"] == 1
    assert x["joint_fraction_of_lrg_pixels"] == 0.5
    print("Random sky-grid and overlap self-tests passed")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", default="eboss_workspace/joint_random_selection_audit.json")
    p.add_argument("--cache-dir", default="eboss_workspace/public_random_catalogues")
    p.add_argument("--max-file-mib", type=int, default=1024)
    p.add_argument("--chunk-rows", type=int, default=200000)
    p.add_argument("--timeout", type=float, default=120)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.max_file_mib < 1 or args.chunk_rows < 1000 or args.timeout <= 0:
        p.error("Invalid file-size, chunk-size or timeout limits")

    records = []
    grids: dict[tuple[str, str], np.ndarray] = {}
    redshift_grids: dict[tuple[str, str], list[np.ndarray]] = {}
    errors = []
    cache = Path(args.cache_dir)
    for (tracer, cap), (filename, nrows) in RANDOMS.items():
        try:
            path, sha, nbytes = download(
                filename, cache, args.max_file_mib * 1024 * 1024, args.timeout)
            record, grid, per_z = examine(
                path, tracer, cap, nrows, sha, nbytes, args.chunk_rows)
            records.append(record)
            grids[(tracer, cap)] = grid
            redshift_grids[(tracer, cap)] = per_z
            print(f"RANDOM_OK {tracer} {cap} rows={nrows} "
                  f"candidate={record['candidate_random_rows']} sha256={sha}", flush=True)
            # Data products are deliberately not uploaded to the Actions artifact.
            path.unlink()
        except (OSError, ValueError, KeyError, RuntimeError, MemoryError) as exc:
            errors.append(f"{filename}: {exc}")
            print("RANDOM_ERROR", errors[-1], flush=True)

    per_cap = {}
    for cap in ("NGC", "SGC"):
        k_lrg, k_elg = ("LRG", cap), ("ELG", cap)
        if k_lrg not in grids or k_elg not in grids:
            continue
        per_cap[cap] = {
            "candidate_interval": overlap(
                grids[k_lrg], grids[k_elg], MIN_RANDOMS_PER_PIXEL),
            "redshift_bins": [
                {"zlo": lo, "zhi": hi,
                 "support": overlap(redshift_grids[k_lrg][i],
                                    redshift_grids[k_elg][i],
                                    MIN_RANDOMS_PER_PIXEL)}
                for i, (lo, hi) in enumerate(zip(Z_EDGES[:-1], Z_EDGES[1:]))
            ],
        }
    status = "four_random_catalogues_checked" if len(records) == 4 and not errors else "partial"
    report = {
        "study": "eBOSS DR16 LRG–ELG pre-unblinding joint random-selection diagnostic",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "catalogue_release": BASE, "status": status,
        "candidate_interval_frozen": False,
        "diagnostic_sky_grid": {"ra_step_deg": RA_STEP,
                                "sin_dec_step": SIN_DEC_STEP,
                                "minimum_randoms_per_tracer_per_pixel": MIN_RANDOMS_PER_PIXEL},
        "random_catalogues": records, "joint_random_support_by_cap": per_cap,
        "errors": errors,
        "data_odd_vector_read": False, "pair_counts_computed": False,
        "window_matrix_computed": False, "wake_amplitude_fitted": False,
        "mock_covariance_calculated": False,
        "scope_note": (
            "These full random-file checks establish file identity, column and "
            "weight health, redshift distributions and gridded overlap. Pixel "
            "support is not an exact common angular/redshift mask. The physical "
            "joint selection, weights, mock adequacy and RR window still require "
            "independent checks before opening any odd-sector measurement."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"JOINT_RANDOM_AUDIT {out} status={status}", flush=True)
    for cap, result in per_cap.items():
        x = result["candidate_interval"]
        print(f"JOINT_SUPPORT {cap} LRG={x['lrg_supported_pixels']} "
              f"ELG={x['elg_supported_pixels']} joint={x['joint_supported_pixels']}",
              flush=True)
    return 0 if status == "four_random_catalogues_checked" else 2


if __name__ == "__main__":
    raise SystemExit(main())
