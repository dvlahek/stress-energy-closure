#!/usr/bin/env python3
"""Inspect DR16 LRG/ELG selection metadata without constructing an odd statistic.

We download only the four public clustering *data* FITS files and record file
hashes, FITS structure, redshift counts, weight-column health and coarse sky
occupancy. Random masks, pair counts and gravitational templates are excluded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.parse import urlsplit
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits

BASE = "https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/"
FILES = {
    ("LRG", "NGC"): ("eBOSS_LRG_clustering_data-NGC-vDR16.fits", 107500),
    ("LRG", "SGC"): ("eBOSS_LRG_clustering_data-SGC-vDR16.fits", 67316),
    ("ELG", "NGC"): ("eBOSS_ELG_clustering_data-NGC-vDR16.fits", 83769),
    ("ELG", "SGC"): ("eBOSS_ELG_clustering_data-SGC-vDR16.fits", 89967),
}
WEIGHTS = ("WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP")
Z_EDGES = (0.6, 0.7, 0.8, 0.9, 1.0, 1.1)


def download(name: str, root: Path, max_bytes: int, timeout: float) -> tuple[Path, str, int]:
    root.mkdir(parents=True, exist_ok=True)
    url = BASE + name
    target = root / name
    tmp = root / (name + ".part")
    digest = hashlib.sha256()
    length = 0
    request = Request(url, headers={"User-Agent": "eBOSS-DR16-selection-audit/1.0"})
    try:
        with urlopen(request, timeout=timeout) as response, tmp.open("wb") as output:
            final = urlsplit(response.geturl())
            expected = urlsplit(BASE)
            if final.scheme != "https" or final.netloc != expected.netloc or not final.path.startswith(expected.path):
                raise ValueError("The catalogue URL redirected outside the approved SDSS directory")
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                length += len(block)
                if length > max_bytes:
                    raise ValueError(f"Catalogue exceeds configured byte limit: {name}")
                output.write(block)
                digest.update(block)
            advertised = response.headers.get("Content-Length", "")
            if advertised.isdigit() and int(advertised) != length:
                raise ValueError(f"HTTP Content-Length mismatch for {name}")
        if not length:
            raise ValueError(f"Empty catalogue: {name}")
        tmp.replace(target)
        return target, digest.hexdigest(), length
    finally:
        if tmp.exists():
            tmp.unlink()


def inspect(path: Path, tracer: str, cap: str, expected_rows: int,
            checksum: str, size: int) -> tuple[dict, set[int]]:
    with fits.open(path, memmap=True) as hdus:
        hdus.verify("exception")
        table_hdus = [hdu for hdu in hdus if isinstance(hdu, fits.BinTableHDU)]
        if len(table_hdus) != 1:
            raise ValueError(f"Expected one FITS binary table, got {len(table_hdus)}: {path.name}")
        table = table_hdus[0]
        names = list(table.columns.names)
        if not all(name in names for name in ("RA", "DEC", "Z", *WEIGHTS)):
            raise ValueError(f"Missing required coordinate/redshift/weight columns in {path.name}")
        if int(table.header["NAXIS2"]) != expected_rows:
            raise ValueError(f"Catalogue row count differs from inventoried header: {path.name}")
        ra = np.asarray(table.data["RA"], dtype=np.float64)
        dec = np.asarray(table.data["DEC"], dtype=np.float64)
        z = np.asarray(table.data["Z"], dtype=np.float64)
        valid = (np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
                 & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90))
        valid_count = int(np.count_nonzero(valid))
        counts = {f"{lo:.1f}_{hi:.1f}": int(np.count_nonzero(valid & (z >= lo) & (z < hi)))
                  for lo, hi in zip(Z_EDGES[:-1], Z_EDGES[1:])}
        candidate = valid & (z >= 0.6) & (z < 1.0)
        # Data-object occupancy is a diagnostic, not a joint survey mask.
        px = np.floor(ra[candidate] / 4.0).astype(np.int32)
        py = np.floor((np.sin(np.deg2rad(dec[candidate])) + 1.0) / 0.05).astype(np.int32)
        pixels = set((px * 41 + np.clip(py, 0, 40)).tolist())
        weight_diagnostics = {}
        for column in WEIGHTS:
            values = np.asarray(table.data[column], dtype=np.float64)
            finite = np.isfinite(values)
            weight_diagnostics[column] = {
                "nonfinite_rows": int(np.count_nonzero(~finite)),
                "nonpositive_rows": int(np.count_nonzero(finite & (values <= 0))),
                "minimum_finite": float(np.min(values[finite])) if finite.any() else None,
                "maximum_finite": float(np.max(values[finite])) if finite.any() else None,
            }
        record = {
            "tracer": tracer, "cap": cap, "filename": path.name,
            "file_sha256": checksum, "downloaded_bytes": size,
            "header_rows": expected_rows, "finite_valid_coordinate_and_z_rows": valid_count,
            "valid_candidate_0p6_to_1p0_rows": int(np.count_nonzero(candidate)),
            "candidate_z_interval_is_preliminary": True,
            "raw_redshift_bin_counts": counts,
            "weight_column_diagnostics": weight_diagnostics,
            "redshift_min_valid": float(np.min(z[valid])) if valid.any() else None,
            "redshift_max_valid": float(np.max(z[valid])) if valid.any() else None,
            "coarse_candidate_data_occupied_cells": len(pixels),
        }
        return record, pixels


def self_test() -> None:
    z = np.array([0.599, 0.6, 0.71, 0.8, 0.999, 1.0, 1.099, 1.1])
    bins = [int(np.count_nonzero((z >= lo) & (z < hi)))
            for lo, hi in zip(Z_EDGES[:-1], Z_EDGES[1:])]
    assert bins == [1, 1, 1, 1, 2]
    print("Selection-only redshift bin test passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/data_selection_audit.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/public_data_catalogues")
    parser.add_argument("--max-file-mib", type=int, default=32)
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.max_file_mib < 1 or args.timeout <= 0:
        parser.error("max-file-mib and timeout must be positive")

    rows = []
    occupancy = {}
    errors = []
    for (tracer, cap), (filename, expected_rows) in FILES.items():
        try:
            local, checksum, size = download(filename, Path(args.cache_dir),
                                             args.max_file_mib * 1024 * 1024,
                                             args.timeout)
            record, pixels = inspect(local, tracer, cap, expected_rows, checksum, size)
            rows.append(record)
            occupancy[(tracer, cap)] = pixels
            print(filename, "SHA256", checksum, "candidate rows",
                  record["valid_candidate_0p6_to_1p0_rows"])
        except (OSError, ValueError, KeyError, RuntimeError) as exc:
            errors.append(f"{filename}: {exc}")
            print("Catalogue validation error:", errors[-1])

    common = {}
    for cap in ("NGC", "SGC"):
        if ("LRG", cap) in occupancy and ("ELG", cap) in occupancy:
            a, b = occupancy["LRG", cap], occupancy["ELG", cap]
            common[cap] = {
                "lrg_occupied_cells": len(a), "elg_occupied_cells": len(b),
                "both_tracers_occupied_cells": len(a & b),
                "interpretation": "Only coarse data-object occupancy; not an angular random-mask intersection.",
            }
    report = {
        "study": "eBOSS DR16 LRG–ELG pre-unblinding data selection audit",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "four_data_catalogues_checked" if len(rows) == 4 and not errors else "partial",
        "catalogue_release": BASE,
        "data_catalogues": rows,
        "preliminary_common_data_occupancy": common,
        "errors": errors,
        "random_catalogues_downloaded": False,
        "mock_catalogues_downloaded": False,
        "pair_counts_computed": False,
        "odd_multipoles_computed": False,
        "physical_template_fitted": False,
        "note": (
            "This audit records data-only metadata before fixing the common "
            "selection. It does not validate the random masks, random weights, "
            "common selection function, mock calibration or an odd-sector statistic."
        ),
    }
    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("Selection audit:", output, "status:", report["status"])
    return 0 if report["status"] == "four_data_catalogues_checked" else 2


if __name__ == "__main__":
    raise SystemExit(main())
