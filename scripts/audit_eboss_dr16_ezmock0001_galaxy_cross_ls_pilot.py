#!/usr/bin/env python3
"""Blinded matched REALISTIC EZmock0001 galaxy + own-random cross-LS pilot.

First authenticate the byte-identical uploaded four mock-galaxy SHA report,
the original source-only protocol, every full galaxy gzip SHA256 and all
four matching same-realization/cap/tracer realistic random gzip SHA256s.
Only AFTER all eight complete-source checks may mock FITS rows be read.
The fixed 600 galaxy + 1200 random per tracer/cap pilot measures pair
orientation/normalization algebra, NOT covariance, physical mask, an odd
detection or the observed eBOSS odd data vector.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from audit_eboss_dr16_ezmock0001_galaxy_raw_sha import (
    ROOT, PROTOCOL as RAW_PROTOCOL, PASS as RAW_PASS,
    ORDER, digest, load, preflight as raw_preflight,
)
from audit_eboss_dr16_9mock_fine_weighted_nz import acquire
from audit_eboss_dr16_rr_pair_closure import rr_histogram, mirrored_closure
from check_eboss_cross_ls_synthetic import (
    cross_landy_szalay, LABELS, audit as scalar_synthetic_audit,
)
from eboss_dr16_fiducial import (
    PRIMARY_GEOMETRY, WEIGHT_COLUMNS, SYSTOT_NUMERICAL_ZERO_TOL,
    comoving_mpc_over_h, validated_weight_product,
)
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path

PILOT = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_protocol_2026-09-26.json"
MANIFEST = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_raw_sha_uploaded_manifest_2026-09-26.json"
ARCHIVE = ROOT / "source_data/eboss_dr16_ezmock0001_galaxy_raw_sha_report_2026-09-26.json"
OLD_RANDOM = ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
OLD_HIGHZ = ROOT / "source_data/eboss_dr16_mock_0001_fine_rr_2026-09-24.json"
STATUS = "EZMOCK0001_MATCHED_GALAXY_CROSS_LS_PILOT_ALGEBRA_ONLY"
STOP = "EZMOCK0001_MATCHED_GALAXY_CROSS_LS_PILOT_INCOMPLETE_STOP"
S_EDGES = np.asarray([20., 40., 60., 80., 100., 120., 140.], dtype="f8")
MU_EDGES = np.linspace(-1. - 1e-7, 1. + 1e-7, 25, dtype="f8")
DATA_COUNT, RANDOM_COUNT = 600, 1200
ORDER_KEYS = tuple(cap + "/" + tracer for cap, tracer in ORDER)
RANDOM_CACHED_DIRS = (
    "eboss_workspace/local_pair_window/mock_fits",
    "eboss_workspace/local_nz/mock_fits",
    "eboss_workspace/mock_selection_catalogues",
    "eboss_workspace/quarantine_realistic_ezmock0001_randoms",
)


def sha_git_blob(raw):
    return hashlib.sha1(
        ("blob " + str(len(raw))).encode("ascii") + bytes([0]) + raw
    ).hexdigest()


def source_manifest_gate(p):
    fixed = load(MANIFEST)
    orig = load(RAW_PROTOCOL)
    raw = ARCHIVE.read_bytes()
    if (
        len(raw) != fixed["uploaded_report_bytes"]
        or hashlib.sha256(raw).hexdigest() != fixed["uploaded_report_sha256"]
        or sha_git_blob(raw) != fixed["uploaded_report_git_blob_sha1"]
        or fixed["uploaded_report_git_blob_sha1"] != fixed["archived_git_blob_sha1"]
        or hashlib.sha256(RAW_PROTOCOL.read_bytes()).hexdigest() !=
           fixed["parent_protocol_sha256"]
        or fixed["uploaded_report_sha256"] !=
           "ae9dc13ca27d7c4f4a6123cbc2526b44866aacc7db3add3f7a2b000773d86b4b"
        or p["prerequisite_source_byte_gate"] !=
           str(RAW_PROTOCOL.relative_to(ROOT))
        or p["prerequisite_source_byte_report_local"] != orig["local_report"]
        or p["mock_id"] != orig["mock_realization_id"]
        or p["mock_id"] != 1
        or p["mock_ensemble_reserved_ids"] != orig["prior_ensemble_ids"]
        or p["caps"] != orig["caps"]
        or p["tracers"] != orig["tracers"]
        or fixed["fixed_source_order"] != list(ORDER_KEYS)
        or p["exact_fixed_highz_bin"] != [0.9, 1.0]
        or p["pilot_counts_per_cap_tracer"] != {
            "mock_galaxy_D": DATA_COUNT, "mock_random_R": RANDOM_COUNT
        }
        or p["sampling"]["seed_root"] != 93127
        or p["geometry"]["s_edges_mpc_h"] != S_EDGES.tolist()
        or p["geometry"]["mu_edges"] !=
           "np.linspace(-1 - 1e-7,1 + 1e-7,25); 24 signed bins"
        or p["geometry"]["fiducial"] != PRIMARY_GEOMETRY
        or p["geometry"]["theta_min_deg"] != 0.05
        or p["geometry"]["line_of_sight"] != "midpoint"
        or p["weights"]["formula"] !=
           "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP"
        or SYSTOT_NUMERICAL_ZERO_TOL != 1e-12
        or p["observed_odd_data_vector_read"] is not False
        or p["new_selection_rule_applied"] is not False
        or p["unblinding_authorized"] is not False
    ):
        raise ValueError("Prior prospective pilot or uploaded full SHA source pin changed")
    report = json.loads(raw)
    if (
        report.get("status") != RAW_PASS
        or report.get("protocol_sha256") != fixed["parent_protocol_sha256"]
        or report.get("mock_realization_id") != 1
        or list(report.get("samples", {})) != list(ORDER_KEYS)
        or report.get("errors") != []
        or report.get("mock_galaxy_rows_read") is not False
        or report.get("mock_random_rows_read") is not False
        or report.get("observed_odd_data_vector_read") is not False
        or report.get("new_science_selection_applied") is not False
        or report.get("total_four_galaxy_compressed_bytes") !=
           fixed["expected_total_compressed_galaxy_bytes"]
    ):
        raise ValueError("Uploaded actual four-mock source SHA report failed its archived contract")
    old_random = load(OLD_RANDOM)
    highz = load(OLD_HIGHZ)
    if (
        old_random.get("status") != "mock_random_sample_selection_compatible"
        or highz.get("status") != "mock_0001_fine_rr_comparison_complete"
        or highz.get("mock_id") != 1
    ):
        raise ValueError("Previously pinned matching mock0001 random evidence incomplete")
    by_random = {
        x["cap"] + "/eBOSS_" + x["tracer"]: x
        for x in old_random["mock_random_catalogues"] if x["id"] == "0001"
    }
    if set(by_random) != set(ORDER_KEYS):
        raise ValueError("Required mock0001 random provenance incomplete")
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        earlier = orig["exact_expected_header_rows_by_cap_tracer"][key]
        frozen = fixed["per_source_frozen_from_uploaded_report"][key]
        sample = report["samples"][key]
        filename = Path(mock_path(tracer, cap, "dat", 1)).name
        expected_suffix = ("/" + orig["local_quarantine_dir"] + "/" +
                           mock_path(tracer, cap, "dat", 1))
        if (
            frozen["released_filename"] != filename
            or frozen["expected_mock_galaxy_header_rows_from_prior_header_audit"] != earlier
            or sample["previously_declared_header_rows_NOT_READ"] != earlier
            or sample["first_seen_full_compressed_sha256"] != frozen["full_compressed_sha256"]
            or sample["compressed_bytes"] != frozen["compressed_bytes"]
            or not sample["local_path"].endswith(expected_suffix)
            or sample["source_url"] != MOCK_BASE + mock_path(tracer, cap, "dat", 1)
            or sample["returned_url"] != sample["source_url"]
            or sample["HTTP_status"] != 200
            or sample["HTTP_content_length_if_available"] != frozen["compressed_bytes"]
            or sample["gzip_magic_checked"] is not True
            or sample["gzip_payload_decompressed"] is not False
            or sample["FITS_headers_or_mock_rows_read"] is not False
            or sample["matched_mock_random_sha256_FROM_PRIOR_AUDIT"] !=
               frozen["matched_same_mock_realization_random_full_compressed_sha256"]
            or by_random[key]["compressed_file_sha256"] !=
               frozen["matched_same_mock_realization_random_full_compressed_sha256"]
        ):
            raise ValueError("Exact uploaded mock galaxy or companion random SHA identity changed: " + key)
    return orig, report, fixed, by_random, highz


def verify_local_source_report(orig):
    raw = ARCHIVE.read_bytes()
    path = ROOT / orig["local_report"]
    if path.read_bytes() != raw:
        raise ValueError("Local raw-source report is not byte-identical to Git-archived user upload")
    return path


def verified_eight_sources(orig, report, fixed, by_random, *, no_download):
    verify_local_source_report(orig)
    raw_preflight(orig)  # Source/ensemble metadata only; never reads FITS.
    sources = {}
    # Full gzip digest of ALL FOUR mock galaxies before reading any FITS row.
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        sample = report["samples"][key]
        path = ROOT / orig["local_quarantine_dir"] / mock_path(tracer, cap, "dat", 1)
        sha, n = digest(path, orig["max_compressed_bytes_per_file"])
        if (
            sha != sample["first_seen_full_compressed_sha256"]
            or sha != fixed["per_source_frozen_from_uploaded_report"][key]["full_compressed_sha256"]
            or n != sample["compressed_bytes"]
        ):
            raise ValueError("Matched mock galaxy binary source SHA changed: " + key)
        sources[key + "/dat"] = path
        print("VERIFIED_MOCK_GALAXY_FULL_SHA", key, n, sha, flush=True)
    # Verify/reuse only exact pre-approved local random caches; otherwise fetch
    # the SAME already pinned SDSS URL and require its complete SHA256.
    for cap, tracer in ORDER:
        key = cap + "/" + tracer
        relative = mock_path(tracer, cap, "ran", 1)
        basename = Path(relative).name
        expected = by_random[key]["compressed_file_sha256"]
        approved_existing = [
            ROOT / folder / basename for folder in RANDOM_CACHED_DIRS
            if (ROOT / folder / basename).exists()
        ]
        if approved_existing:
            path = approved_existing[0]
            for candidate in approved_existing:
                sha, n = digest(candidate, 512 * 1024 * 1024)
                if sha != expected:
                    raise ValueError("Existing expected mock0001 random cache SHA differs: " + str(candidate))
                if candidate == path:
                    size = n
        else:
            path = ROOT / RANDOM_CACHED_DIRS[-1] / basename
            obtained, sha, size = acquire(
                path, expected, observed=False, url=MOCK_BASE + relative,
                no_download=no_download, timeout=120,
            )
            if obtained != path or sha != expected:
                raise ValueError("Pinned companion random acquisition has changed identity")
        if path.name != basename:
            raise ValueError("Random source basename changed")
        sources[key + "/ran"] = path
        print("VERIFIED_MATCHED_MOCK_RANDOM_FULL_SHA", key, size, expected, flush=True)
    if set(sources) != {key + "/" + role for key in ORDER_KEYS for role in ("dat", "ran")}:
        raise ValueError("Eight source FITS binaries not all pinned before row access")
    return sources


def sample_catalogue(path, *, expected_rows, cap, tracer, role, expected_highz, p):
    high = p["exact_fixed_highz_bin"]
    n_select = DATA_COUNT if role == "dat" else RANDOM_COUNT
    seed = (p["sampling"]["seed_root"] + 10000 * p["caps"].index(cap)
            + 100 * p["tracers"].index(tracer)
            + (0 if role == "dat" else 1))
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected exactly one matched realistic EZmock BINTABLE")
        table = tables[0]
        if int(table.header["NAXIS2"]) != expected_rows:
            raise ValueError("Expected mock source FITS header row count differs: " + str(path))
        names = set(table.columns.names)
        required = {"RA", "DEC", "Z", *WEIGHT_COLUMNS}
        if tracer == "eBOSS_ELG":
            required.add("chunk")
        if not required.issubset(names):
            raise ValueError("Required exact mock source columns missing: " + str(required - names))
        data = table.data
        z = np.asarray(data["Z"], dtype="f8")
        ra = np.asarray(data["RA"], dtype="f8")
        dec = np.asarray(data["DEC"], dtype="f8")
        if (len(z) != expected_rows or not np.isfinite(z).all()
                or not np.isfinite(ra).all() or not np.isfinite(dec).all()
                or np.any(ra < 0) or np.any(ra >= 360)
                or np.any(dec < -90) or np.any(dec > 90)):
            raise ValueError("Invalid coordinate/redshift input in complete pinned mock source")
        candidate = (z >= high[0]) & (z < high[1])
        raw_candidate = int(np.count_nonzero(candidate))
        fields = {
            name: np.asarray(data[name][candidate], dtype="f8")
            for name in WEIGHT_COLUMNS
        }
        weight, retained = validated_weight_product(fields)
        eligible = np.flatnonzero(candidate)[retained]
        eligible_weight = weight[retained]
        if expected_highz is not None and len(eligible) != expected_highz:
            raise ValueError("Earlier source-pinned mock random high-z count differs")
        if len(eligible) < n_select:
            raise ValueError("Too few eligible rows for fixed D/R count: " + str(path))
        choose = np.sort(np.random.default_rng(seed).choice(
            len(eligible), size=n_select, replace=False
        ))
        indices = eligible[choose]
        cat = tuple(np.asarray(v, dtype="f8").copy() for v in (
            ra[indices], dec[indices], z[indices], eligible_weight[choose]
        ))
        if (
            len(np.unique(indices)) != n_select
            or not all(np.isfinite(v).all() for v in cat)
            or not np.all(cat[3] > 0)
            or not np.all((cat[2] >= high[0]) & (cat[2] < high[1]))
        ):
            raise ValueError("Deterministic selected mock rows are invalid")
        digest = hashlib.sha256()
        for vec in cat:
            digest.update(np.ascontiguousarray(vec).tobytes())
        chunk_diag = None
        if tracer == "eBOSS_ELG":
            eligible_labels = np.asarray(data["chunk"])[eligible]
            sample_labels = np.asarray(data["chunk"])[indices]
            def label(value):
                if isinstance(value, (bytes, np.bytes_)):
                    value = value.decode("ascii", errors="strict")
                out = str(value).strip()
                if not out or len(out) > 128:
                    raise ValueError("Invalid ELG chunk label in official mock source")
                return out
            chunk_diag = {
                "eligible_by_exact_chunk": {
                    label(v): int(np.count_nonzero(eligible_labels == v))
                    for v in np.unique(eligible_labels)
                },
                "selected_by_exact_chunk": {
                    label(v): int(np.count_nonzero(sample_labels == v))
                    for v in np.unique(sample_labels)
                },
                "chunk_conditioned_sampling_or_weight_change": False,
            }
    return cat, {
        "cap": cap, "tracer": tracer, "role": role, "seed": seed,
        "input_FITS_rows": expected_rows,
        "candidate_highz_rows_before_weight_validation": raw_candidate,
        "numerical_zero_weight_excluded_rows": raw_candidate - len(eligible),
        "eligible_highz_rows_after_fixed_weight_gate": int(len(eligible)),
        "selected_rows": n_select,
        "selected_array_SHA256": digest.hexdigest(),
        "sum_selected_weights": float(cat[3].sum(dtype="f8")),
        "ELG_exact_chunk_diagnostic": chunk_diag,
    }


def oriented_pair_terms(d_l, d_e, r_l, r_e, *, distance):
    samples = {
        "D1D2": (d_l, d_e), "D1R2": (d_l, r_e),
        "R1D2": (r_l, d_e), "R1R2": (r_l, r_e),
    }
    hist, norms, meta = {}, {}, {}
    for term in LABELS:
        first, second = samples[term]
        h, info = rr_histogram(
            first, second, S_EDGES, MU_EDGES, 0.05, distance, block=128
        )
        independent_norm = (
            first[3].sum(dtype="f8") * second[3].sum(dtype="f8")
        )
        rel = abs(info["pair_normalization"] / independent_norm - 1.0)
        if (
            not np.isfinite(h).all() or np.any(h < 0)
            or h.shape != (6, 24) or rel >= 1e-12
            or info["accepted_pairs"] < 0
        ):
            raise ValueError("Nonfinite or misnormalized independently oriented pair term: " + term)
        hist[term], norms[term] = h, info["pair_normalization"]
        meta[term] = {
            "accepted_pairs": info["accepted_pairs"],
            "total_weighted_pairs_in_fixed_s_mu_bins": float(h.sum()),
            "independently_normalized_pair_weight": float(independent_norm),
            "normalization_relative_residual": float(rel),
            "weighted_histogram_SHA256": hashlib.sha256(
                np.ascontiguousarray(h).tobytes()
            ).hexdigest(),
        }
    return hist, norms, meta


def closure_for_cap(d_l, d_e, r_l, r_e):
    distance = lambda z: comoving_mpc_over_h(z, PRIMARY_GEOMETRY)
    fwd, norms, fmeta = oriented_pair_terms(
        d_l, d_e, r_l, r_e, distance=distance
    )
    rev, rnorms, rmeta = oriented_pair_terms(
        d_e, d_l, r_e, r_l, distance=distance
    )
    mapping = {"D1D2": "D1D2", "D1R2": "R1D2",
               "R1D2": "D1R2", "R1R2": "R1R2"}
    mirror = {}
    for label, reversed_label in mapping.items():
        # Reversal normalizations must swap together with the cross-term.
        a = dict(pair_normalization=norms[label],
                 accepted_pairs=fmeta[label]["accepted_pairs"])
        b = dict(pair_normalization=rnorms[reversed_label],
                 accepted_pairs=rmeta[reversed_label]["accepted_pairs"])
        found = mirrored_closure(fwd[label], rev[reversed_label], a, b, MU_EDGES)
        if found["closure_passed"] is not True:
            raise ValueError("True mock-galaxy cross-pair reversal check failed: " + label)
        mirror[label] = found
    xi, support = cross_landy_szalay(fwd, norms)
    xirev, reverse_support = cross_landy_szalay(rev, rnorms)
    if not np.array_equal(support, reverse_support[:, ::-1]) or not np.any(support):
        raise ValueError("Matched mock galaxy RR supported signed-mu cells differ on reverse")
    err = float(np.max(np.abs(xi[support] - xirev[:, ::-1][support])))
    if not np.isfinite(err) or err >= 1e-8:
        raise ValueError("Matched mock galaxy cross-LS field mirror closure failed")
    rr, rrrev = fwd["R1R2"], rev["R1R2"]
    parities = {}
    for ell in (1, 3):
        if ell == 1:
            integral = (MU_EDGES[1:]**2 - MU_EDGES[:-1]**2) / 2.
        else:
            integral = ((5./8.) * (MU_EDGES[1:]**4 - MU_EDGES[:-1]**4)
                        - (3./4.) * (MU_EDGES[1:]**2 - MU_EDGES[:-1]**2))
        avg = integral / np.diff(MU_EDGES)
        wa, wb = rr.sum(), rrrev.sum()
        if wa <= 0 or wb <= 0:
            raise ValueError("Empty fixed R_L R_E pair histogram")
        a = float((rr @ avg).sum() / wa)
        b = float((rrrev @ avg).sum() / wb)
        residual = float(abs(a+b))
        if not np.isfinite(residual) or residual >= 1e-10:
            raise ValueError("Mock random RR orientation odd-moment algebra failed")
        parities[str(ell)] = {"raw_random_odd_moment_mirror_sum_abs": residual}
    return {
        "status": "matched_realistic_mock_galaxy_cross_ls_orientation_checked",
        "forward_pair_terms": fmeta,
        "reverse_pair_terms": rmeta,
        "forward_reverse_pair_mirror": mirror,
        "RR_supported_s_mu_cells": int(np.count_nonzero(support)),
        "RR_total_s_mu_cells": int(support.size),
        "cross_xi_reverse_max_abs_residual": err,
        "raw_RR_odd_parity_1_and_3": parities,
        "mock_galaxy_one_realization_pilot_not_statistical_inference": True,
        "unsupported_mu_cells_not_zero_filled_for_physical_multipoles": True,
    }


def self_test(p):
    orig, report, fixed, byrandom, highz = source_manifest_gate(p)
    if not (len(report["samples"]) == 4
            and len(byrandom) == 4
            and fixed["expected_total_compressed_galaxy_bytes"] == 8307621):
        raise AssertionError("Four exact source fingerprints not available")
    test = scalar_synthetic_audit(0.05)
    if test["supported_s_mu_cells"] <= 0 or not test["forward_reverse_xi_mirror_passed"]:
        raise AssertionError("Independent scalar synthetic cross-LS check failed")
    print("EBOSS_EZMOCK0001_GALAXY_CROSS_LS_SYNTHETIC_SELF_TEST_OK", flush=True)


def atomic(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(".tmp.json")
    part.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")
    part.replace(path)


def run(p, *, no_download):
    orig, report, fixed, byrandom, highz = source_manifest_gate(p)
    paths = verified_eight_sources(
        orig, report, fixed, byrandom, no_download=no_download
    )
    result = {
        "status": STATUS,
        "pilot_mock_id": 1,
        "caps": list(p["caps"]),
        "tracers": list(p["tracers"]),
        "exact_uploaded_raw_report_sha256": fixed["uploaded_report_sha256"],
        "source_galaxy_full_compressed_SHA256": {
            key: fixed["per_source_frozen_from_uploaded_report"][key]["full_compressed_sha256"]
            for key in ORDER_KEYS
        },
        "matched_source_random_full_compressed_SHA256": {
            key: byrandom[key]["compressed_file_sha256"]
            for key in ORDER_KEYS
        },
        "fixed_geometry": p["geometry"],
        "fixed_D_R_counts_per_tracer_cap": p["pilot_counts_per_cap_tracer"],
        "cases": [],
        "mock_galaxy_rows_read": True,
        "mock_random_rows_read": True,
        "observed_galaxy_rows_read": False,
        "observed_odd_data_vector_read": False,
        "new_science_selection_applied": False,
        "physical_mask_certified": False,
        "physical_pair_window_certified": False,
        "mock_covariance_computed": False,
        "not_a_detection_or_hypothesis_test": True,
        "errors": [],
    }
    for cap in p["caps"]:
        inputs, metas = {}, {}
        for tracer in p["tracers"]:
            key = cap + "/" + tracer
            for role in ("dat", "ran"):
                sample, meta = sample_catalogue(
                    paths[key + "/" + role],
                    expected_rows=(
                        fixed["per_source_frozen_from_uploaded_report"][key][
                            "expected_mock_galaxy_header_rows_from_prior_header_audit"
                        ] if role == "dat"
                        else int(byrandom[key]["header_rows"])
                    ),
                    cap=cap, tracer=tracer, role=role,
                    expected_highz=(
                        int(highz["caps"][cap]["mock_randoms"][tracer.removeprefix("eBOSS_")])
                        if role == "ran" else None
                    ), p=p
                )
                inputs[tracer, role] = sample
                metas[tracer + "_" + role] = meta
                print("SELECTED_FIXED_MOCK_INPUT", cap, tracer, role,
                      meta["eligible_highz_rows_after_fixed_weight_gate"],
                      meta["selected_rows"], flush=True)
        closure = closure_for_cap(
            inputs["eBOSS_LRG", "dat"], inputs["eBOSS_ELG", "dat"],
            inputs["eBOSS_LRG", "ran"], inputs["eBOSS_ELG", "ran"],
        )
        result["cases"].append({"cap": cap, "input_sample_diagnostics": metas,
                                "cross_ls": closure})
        print("MATCHED_MOCK_GALAXY_CROSS_LS_CAP", cap,
              closure["RR_supported_s_mu_cells"],
              closure["cross_xi_reverse_max_abs_residual"], flush=True)
    if len(result["cases"]) != 2:
        raise ValueError("Both predeclared Galactic caps required")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--no-download", action="store_true",
                    help="Require all four already pinned mock randoms present in known caches")
    args = ap.parse_args()
    p = load(PILOT)
    if args.self_test:
        self_test(p)
        return 0
    out = ROOT / "eboss_workspace/official_mask_inventory/ezmock0001_galaxy_cross_ls_pilot.json"
    try:
        result = run(p, no_download=args.no_download)
    except Exception as exc:
        result = {
            "status": STOP,
            "errors": [str(exc)],
            "mock_galaxy_rows_may_have_been_read": True,
            "observed_galaxy_rows_read": False,
            "observed_odd_data_vector_read": False,
            "new_science_selection_applied": False,
            "physical_mask_certified": False,
            "physical_pair_window_certified": False,
        }
    atomic(out, result)
    print("EBOSS_EZMOCK0001_GALAXY_CROSS_LS_PILOT", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep="\n", flush=True)
        return 2
    print("COMPLETED_CAPS", len(result["cases"]), flush=True)
    print("OBSERVED_ODD_DATA_READ", result["observed_odd_data_vector_read"],
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
