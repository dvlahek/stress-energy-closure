#!/usr/bin/env python3
"""Blinded, fixed 0.01-z weighted n(z) control on eBOSS DR16 RANDOMS.

Process each observed/random or realistic EZmock/random FITS catalogue once
on local WSL. Input SHA256s come from existing committed observed provenance
and the fixed successful nine-mock random-only ensemble. Per-cap/tracer,
and for ELG per exact matching CHUNK label, compare the forty predeclared
weighted 0.01-redshift bins on 0.6 <= z < 1.0. No observed/mock galaxy
positions or observed odd vector are opened, and no cuts are optimized.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
from astropy.io import fits

from aggregate_eboss_dr16_9mock_window import IDS, CAPS, TRACERS
from audit_eboss_dr16_full_random_rr import ROOT
from eboss_dr16_fiducial import (
    WEIGHT_COLUMNS, SYSTOT_NUMERICAL_ZERO_TOL, validated_weight_product,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import fetch_with_retry
from inspect_eboss_dr16_selection import download

PROTOCOL = ROOT / "source_data/eboss_dr16_9mock_fine_weighted_nz_protocol_2026-09-24.json"
PROTOCOL_COMMIT = "6322adaa7b4fbac0020cd7dbb4874e07ff3b4e67"
OBS_REF = ROOT / "source_data/eboss_dr16_weighted_normalization_2026-09-24.json"
COARSE_REF = ROOT / "source_data/eboss_dr16_full_observed_rr_2026-09-24.json"
PRIOR_MOCK_REF = ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
SOURCE_ENSEMBLE_RUN = "36017670812"
SOURCE_ENSEMBLE_SHA = "054005edc6ca193176b129e1951e4bd3ce8751a7"
SOURCE_MOCK_RUN = "36015246541"
SOURCE_MOCK_SHA = "9150992ac06c03cf1b16c5fdd4e3bb3856f89ff0"
FINE_EDGES = np.linspace(0.6, 1.0, 41, dtype="f8")
BLOCK_ROWS = 100000


def sha256_bytes(path: Path):
    digest, nbytes = hashlib.sha256(), 0
    with path.open("rb") as src:
        for block in iter(lambda: src.read(1024 * 1024), b""):
            digest.update(block)
            nbytes += len(block)
    return digest.hexdigest(), nbytes


def acquire(path, expected_sha, *, observed, filename=None, url=None,
            timeout=120, no_download=False):
    if path.is_file():
        digest, size = sha256_bytes(path)
        if digest != expected_sha:
            raise ValueError(f"Invalid cached SHA256 for {path}; remove it explicitly")
        print("FINE_NZ_VERIFIED_CACHE", path, digest, flush=True)
        return path, digest, size
    if no_download:
        raise FileNotFoundError(f"No SHA-verified catalogue in cache: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    # A byte-complete HTTP transfer can still fail the pinned SHA256:
    # retry the *same pinned URL*, never replace the expected fingerprint.
    # Every mismatching payload is deleted, and every actual SHA is logged.
    # Retry only transport failures for the same exact, pinned release
    # URL. This does not relax the byte identity or accept a new source.
    attempts = 3
    last_digest, last_size = None, None
    for attempt in range(1, attempts + 1):
        try:
            if observed:
                actual, digest, size = download(
                    filename, path.parent, 1024 * 1024 * 1024, timeout)
                if actual != path:
                    raise ValueError("Observed random download landed at unexpected path")
            else:
                digest, size = fetch_with_retry(
                    url, path, 512 * 1024 * 1024, timeout)
        except (OSError, TimeoutError) as exc:
            if attempt == attempts:
                raise
            print("FINE_NZ_PINNED_INPUT_TRANSPORT_RETRY",
                  path.name, f"attempt={attempt}/{attempts}",
                  f"error={type(exc).__name__}", flush=True)
            time.sleep(2 * attempt)
            continue
        if digest == expected_sha:
            print("FINE_NZ_INPUT_SHA_OK", path.name, size, digest, flush=True)
            return path, digest, size
        last_digest, last_size = digest, size
        path.unlink(missing_ok=True)
        print("FINE_NZ_SHA_MISMATCH", path.name,
              f"attempt={attempt}/{attempts}",
              f"bytes={size}", f"expected={expected_sha}",
              f"actual={digest}", flush=True)
        if attempt < attempts:
            time.sleep(2 * attempt)
    raise ValueError(
        "Pinned random FITS SHA256 mismatch after "
        f"{attempts} attempt(s): {path}; expected={expected_sha}; "
        f"last_actual={last_digest}; last_bytes={last_size}. "
        "Do not repin or skip the fixed mock based on this failure.")


def chunk_label(value):
    if isinstance(value, (bytes, np.bytes_)):
        s = value.decode("ascii", errors="strict")
    else:
        s = str(value)
    label = s.strip()
    if not label or len(label) > 128:
        raise ValueError("Empty or excessive ELG chunk label in candidate random")
    return label


def hist_component(z, w):
    counts = np.histogram(z, bins=FINE_EDGES)[0].astype("i8")
    weighted = np.histogram(z, bins=FINE_EDGES, weights=w)[0].astype("f8")
    squared = np.histogram(z, bins=FINE_EDGES, weights=w*w)[0].astype("f8")
    return counts, weighted, squared


def add_component(destination, component):
    for key, values in zip(("count", "weighted", "weight_squared"), component):
        destination[key] += values


def empty_component():
    return {
        "count": np.zeros(40, dtype="i8"),
        "weighted": np.zeros(40, dtype="f8"),
        "weight_squared": np.zeros(40, dtype="f8"),
    }


def final_component(component):
    n = int(np.sum(component["count"], dtype="i8"))
    sw = float(np.sum(component["weighted"], dtype="f8"))
    sw2 = float(np.sum(component["weight_squared"], dtype="f8"))
    if n <= 0 or sw <= 0 or sw2 <= 0:
        raise ValueError("Empty or invalid fixed weighted n(z) component")
    if not np.isfinite(sw + sw2):
        raise ValueError("Nonfinite weighted n(z) normalization")
    return {
        "retained_rows": n,
        "sum_weights": sw,
        "sum_squared_weights": sw2,
        "kish_effective_count": sw * sw / sw2,
        "count_per_fine_bin": component["count"].tolist(),
        "sum_weights_per_fine_bin": component["weighted"].tolist(),
        "normalized_weighted_nz": (component["weighted"] / sw).tolist(),
        "normalized_unweighted_nz": (component["count"] / n).tolist(),
        "coarse_0p1_weighted_fractions": (
            component["weighted"].reshape(4, 10).sum(axis=1) / sw
        ).tolist(),
    }


def scan_fits(path, *, tracer, cap, kind, digest, size, expected_rows=None):
    if tracer not in TRACERS or cap not in CAPS or kind not in ("observed", "mock"):
        raise ValueError("Unsupported fixed catalogue type")
    whole, by_chunk = empty_component(), {}
    raw_candidate = excluded = tiny = tol = scanned = 0
    with fits.open(path, memmap=False) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected exactly one random FITS BINTABLE")
        hdu = tables[0]
        nrows = int(hdu.header["NAXIS2"])
        if expected_rows is not None and nrows != expected_rows:
            raise ValueError("Previously certified FITS header row count changed")
        required = {"RA", "DEC", "Z", *WEIGHT_COLUMNS}
        if tracer == "ELG":
            required.add("chunk")
        missing = required - set(hdu.columns.names)
        if missing:
            raise ValueError("Missing required random FITS columns: " + str(sorted(missing)))
        for begin in range(0, nrows, BLOCK_ROWS):
            batch = hdu.data[begin:min(nrows, begin + BLOCK_ROWS)]
            z = np.asarray(batch["Z"], dtype="f8")
            ra = np.asarray(batch["RA"], dtype="f8")
            dec = np.asarray(batch["DEC"], dtype="f8")
            scanned += len(z)
            if (not np.isfinite(z).all() or not np.isfinite(ra).all()
                    or not np.isfinite(dec).all()
                    or np.any(ra < 0) or np.any(ra >= 360)
                    or np.any(dec < -90) or np.any(dec > 90)):
                raise ValueError("Nonfinite or out-of-range random coordinates/redshifts")
            selected = (z >= FINE_EDGES[0]) & (z < FINE_EDGES[-1])
            if not np.any(selected):
                continue
            raw_candidate += int(np.count_nonzero(selected))
            fields = {
                name: np.asarray(batch[name][selected], dtype="f8")
                for name in WEIGHT_COLUMNS
            }
            tiny += int(np.count_nonzero(
                np.abs(fields["WEIGHT_SYSTOT"]) <= 1e-20))
            tol += int(np.count_nonzero(
                np.abs(fields["WEIGHT_SYSTOT"]) <= SYSTOT_NUMERICAL_ZERO_TOL))
            weights, retained = validated_weight_product(fields)
            excluded += int(np.count_nonzero(~retained))
            z_keep = z[selected][retained]
            w_keep = weights[retained]
            if z_keep.size == 0:
                continue
            add_component(whole, hist_component(z_keep, w_keep))
            if tracer == "ELG":
                raw_labels = np.asarray(batch["chunk"])[selected][retained]
                for value in np.unique(raw_labels):
                    label = chunk_label(value)
                    mask = raw_labels == value
                    item = by_chunk.setdefault(label, empty_component())
                    add_component(item, hist_component(z_keep[mask], w_keep[mask]))
    if scanned != nrows or tiny != tol or excluded != tol:
        raise ValueError("FITS rows or predeclared numerical-zero weight gate inconsistent")
    if int(np.sum(whole["count"])) != raw_candidate - excluded:
        raise ValueError("Fixed fine n(z) bins fail retained candidate row conservation")
    if tracer == "ELG":
        if not by_chunk:
            raise ValueError("Required ELG CHUNK labels absent")
        for field in ("count", "weighted", "weight_squared"):
            reconstructed = np.sum(
                [v[field] for v in by_chunk.values()], axis=0)
            np.testing.assert_allclose(
                reconstructed, whole[field], rtol=5e-12, atol=1e-7)
    result = {
        "status": "fine_weighted_random_catalogue_checked",
        "protocol_commit": PROTOCOL_COMMIT,
        "kind": kind, "cap": cap, "tracer": tracer,
        "source_filename": path.name, "source_sha256": digest,
        "source_file_bytes": size, "header_rows": nrows,
        "candidate_raw_rows": raw_candidate,
        "excluded_numerical_zero_systot_rows": excluded,
        "near_zero_le_1e_20_rows": tiny,
        "near_zero_le_1e_12_rows": tol,
        "fine_z_edges": FINE_EDGES.tolist(),
        "weight_convention":
            "WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP",
        "full": final_component(whole),
        "per_exact_chunk": {
            key: final_component(by_chunk[key]) for key in sorted(by_chunk)
        },
        "observed_galaxy_data_read": False,
        "mock_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
    }
    return result


def normalized_tv(p, q):
    p, q = np.asarray(p, dtype="f8"), np.asarray(q, dtype="f8")
    if p.shape != (40,) or q.shape != (40,) or not np.isfinite(p+q).all():
        raise ValueError("Invalid fixed forty-bin weighted n(z)")
    if abs(float(np.sum(p)) - 1) > 1e-10 or abs(float(np.sum(q)) - 1) > 1e-10:
        raise ValueError("Fine n(z) does not normalize to one")
    return {
        "weighted_normalized_total_variation": float(
            0.5 * np.sum(np.abs(p-q))),
        "weighted_normalized_max_abs_fine_bin_difference": float(
            np.max(np.abs(p-q))),
        "coarse_weighted_z_fraction_max_abs_difference": float(np.max(np.abs(
            p.reshape(4,10).sum(axis=1)-q.reshape(4,10).sum(axis=1)))),
    }


def compare(observed, mock, *, mid):
    if any(observed[k] != mock[k]
           for k in ("cap", "tracer", "fine_z_edges", "weight_convention")):
        raise ValueError("Observed/mock selection conventions differ")
    obs_chunk = observed["per_exact_chunk"]
    mock_chunk = mock["per_exact_chunk"]
    shared = sorted(set(obs_chunk) & set(mock_chunk))
    row = {
        "mock_id": mid, "cap": observed["cap"], "tracer": observed["tracer"],
        "observed_sha256": observed["source_sha256"],
        "mock_sha256": mock["source_sha256"],
        "full_weighted_nz_difference": normalized_tv(
            observed["full"]["normalized_weighted_nz"],
            mock["full"]["normalized_weighted_nz"]),
        "exactly_shared_elg_chunks": shared,
        "observed_only_elg_chunks": sorted(set(obs_chunk)-set(mock_chunk)),
        "mock_only_elg_chunks": sorted(set(mock_chunk)-set(obs_chunk)),
        "per_shared_chunk_difference": {
            label: {
                **normalized_tv(
                    obs_chunk[label]["normalized_weighted_nz"],
                    mock_chunk[label]["normalized_weighted_nz"]),
                "observed_retained_rows": obs_chunk[label]["retained_rows"],
                "mock_retained_rows": mock_chunk[label]["retained_rows"],
            } for label in shared
        },
    }
    return row


def json_write_atomic(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    tmp.replace(path)


def cached_record(path, sha, cap, tracer, kind):
    if not path.is_file():
        return None
    record = json.loads(path.read_text())
    if (record.get("status") != "fine_weighted_random_catalogue_checked"
            or record.get("protocol_commit") != PROTOCOL_COMMIT
            or record.get("source_sha256") != sha
            or record.get("cap") != cap
            or record.get("tracer") != tracer
            or record.get("kind") != kind
            or record.get("fine_z_edges") != FINE_EDGES.tolist()
            or record.get("observed_odd_data_vector_read") is not False):
        raise ValueError("Existing checkpoint does not match pinned run: " + str(path))
    return record


def preflight(args):
    manifest = json.loads(PROTOCOL.read_text())
    if (tuple(manifest["mock_ids"]) != IDS
            or tuple(manifest["caps"]) != CAPS
            or tuple(manifest["tracers"]) != TRACERS
            or manifest["fine_z_edges_definition"]
            != "np.linspace(0.6,1.0,41,dtype=np.float64); forty fixed bins, width 0.01"
            or manifest["observed_odd_data_vector_read"] is not False):
        raise ValueError("Fixed weighted n(z) registration changed")
    if args.ensemble_json is None:
        raise ValueError("Require downloaded pinned nine-mock ensemble JSON")
    ensemble = json.loads(Path(args.ensemble_json).read_text())
    if (ensemble.get("status") != "nine_mock_random_only_window_pilot_complete"
            or ensemble.get("revision_commit") != SOURCE_ENSEMBLE_SHA
            or tuple(ensemble.get("predeclared_mock_ids", [])) != IDS
            or ensemble.get("source_mock_shard_workflow_run") != SOURCE_MOCK_RUN
            or ensemble.get("source_mock_shard_commit") != SOURCE_MOCK_SHA
            or ensemble.get("observed_odd_data_vector_read") is not False
            or len(ensemble.get("cases", [])) != 18):
        raise ValueError("Nine-mock input ensemble source is incomplete or not blinded")
    mock_sha = {}
    for item in ensemble["cases"]:
        mid, cap = int(item["mock_id"]), item["cap"]
        for tracer in TRACERS:
            digest = item["input_file_sha256"][tracer]
            if (mid not in IDS or cap not in CAPS
                    or len(digest) != 64
                    or (mid,cap,tracer) in mock_sha):
                raise ValueError("Invalid or duplicate pinned mock catalogue SHA")
            mock_sha[mid, cap, tracer] = digest
    if len(mock_sha) != 36:
        raise ValueError("Full fixed cohort requires 36 matched mock random SHA256s")
    previous = json.loads(PRIOR_MOCK_REF.read_text())
    if previous.get("status") != "mock_random_sample_selection_compatible":
        raise ValueError("Earlier three mock inputs not certified")
    previous_files = {
        (int(item["id"]), item["cap"], item["tracer"]): item
        for item in previous["mock_random_catalogues"]
    }
    for key, item in previous_files.items():
        if mock_sha[key] != item["compressed_file_sha256"]:
            raise ValueError("Prior three-mock SHA identity differs from nine-mock source")

    obs_ref = json.loads(OBS_REF.read_text())
    obs_sha = {
        (item["cap"], item["tracer"]): item["sha256"]
        for item in obs_ref["per_file_records"]
        if item["survey"] == "observed" and item["role"] == "random"
    }
    coarse = json.loads(COARSE_REF.read_text())
    obs_counts = {
        (item["cap"], item["tracer"], item["zlo"], item["zhi"]):
            item["retained_rows"]
        for item in coarse["input_random_redshift_counts"]
    }
    if len(obs_sha) != 4 or len(obs_counts) != 16:
        raise ValueError("Observed random source count or SHA provenance incomplete")
    return mock_sha, obs_sha, obs_counts, previous_files


def self_test():
    np.testing.assert_allclose(np.diff(FINE_EDGES), np.full(40, 0.01),
                               rtol=0, atol=1e-14)
    z = np.array([0.599,0.6,0.605,0.615,0.999,1.0])
    w = np.array([0.0,2.,3.,4.,5.,0.])
    c, h, h2 = hist_component(z[1:-1], w[1:-1])
    assert int(c.sum()) == 4
    np.testing.assert_allclose(h.sum(), 14.)
    np.testing.assert_allclose(h2.sum(), 2**2+3**2+4**2+5**2)
    full = empty_component()
    add_component(full, (c,h,h2))
    x = final_component(full)
    assert abs(sum(x["normalized_weighted_nz"])-1) < 1e-14
    assert normalized_tv(x["normalized_weighted_nz"],
                         x["normalized_weighted_nz"])[
        "weighted_normalized_total_variation"] == 0
    assert chunk_label(b" test_A  ") == "test_A"
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as tmp:
        zrows = np.array([.59,.60,.61,.69,.71,.99,1.0,.75])
        wsys = np.array([1.,1.,1.,1.,1e-31,1.,1.,1.])
        cols = [
            fits.Column(name="RA", format="D",
                        array=np.arange(8,dtype="f8")+50),
            fits.Column(name="DEC", format="D",
                        array=np.arange(8,dtype="f8")),
            fits.Column(name="Z", format="D", array=zrows),
            fits.Column(name="WEIGHT_SYSTOT", format="D", array=wsys),
            fits.Column(name="WEIGHT_CP", format="D",
                        array=np.ones(8)),
            fits.Column(name="WEIGHT_NOZ", format="D",
                        array=np.ones(8)),
            fits.Column(name="WEIGHT_FKP", format="D",
                        array=2*np.ones(8)),
            fits.Column(name="chunk", format="8A",
                        array=np.array(["A","A","A","B","B","B","C","B"])),
        ]
        path = Path(tmp)/"synthetic_random.fits"
        fits.HDUList([fits.PrimaryHDU(),
                      fits.BinTableHDU.from_columns(cols)]).writeto(path)
        digest,size = sha256_bytes(path)
        record = scan_fits(path,tracer="ELG",cap="SGC",kind="mock",
                           digest=digest,size=size,expected_rows=8)
        assert record["candidate_raw_rows"] == 6
        assert record["excluded_numerical_zero_systot_rows"] == 1
        assert record["full"]["retained_rows"] == 5
        np.testing.assert_allclose(record["full"]["sum_weights"],10.)
        assert set(record["per_exact_chunk"]) == {"A","B"}
        assert record["per_exact_chunk"]["A"]["retained_rows"] == 2
        assert record["per_exact_chunk"]["B"]["retained_rows"] == 3
        np.testing.assert_allclose(
            record["full"]["normalized_weighted_nz"],
            record["full"]["normalized_unweighted_nz"])
        assert compare(record,record,mid=1)["full_weighted_nz_difference"][
            "weighted_normalized_total_variation"] == 0
    # Integrity recovery unit test: a wrong but complete first payload
    # must be rejected, and the same immutable pin is retained on retry.
    with TemporaryDirectory() as tmp:
        good, bad = b"fixed-pinned-synthetic-fits", b"wrong-response"
        expected = hashlib.sha256(good).hexdigest()
        previous_fetch = globals()["fetch_with_retry"]
        seen = []
        def mock_fetch(url, destination, byte_limit, timeout):
            seen.append(url)
            payload = bad if len(seen) == 1 else good
            destination.write_bytes(payload)
            return hashlib.sha256(payload).hexdigest(), len(payload)
        globals()["fetch_with_retry"] = mock_fetch
        try:
            local = Path(tmp) / "synthetic.ran.fits.gz"
            got = acquire(local, expected, observed=False,
                          url=MOCK_BASE + "synthetic", timeout=1)
            assert len(seen) == 2 and got[1] == expected
            np.testing.assert_array_equal(
                np.frombuffer(local.read_bytes(), dtype="u1"),
                np.frombuffer(good, dtype="u1"))
            local.unlink()
            seen.clear()
            def always_bad(url, destination, byte_limit, timeout):
                seen.append(url)
                destination.write_bytes(bad)
                return hashlib.sha256(bad).hexdigest(), len(bad)
            globals()["fetch_with_retry"] = always_bad
            try:
                acquire(local, expected, observed=False,
                        url=MOCK_BASE + "synthetic", timeout=1)
            except ValueError as exc:
                assert expected in str(exc) and "last_actual=" in str(exc)
            else:
                raise AssertionError("A persistent SHA mismatch must fail closed")
            assert len(seen) == 3 and not local.exists()
        finally:
            globals()["fetch_with_retry"] = previous_fetch
    print("EBOSS_FINE_NZ_PINNED_SHA_RETRY_SELF_TEST_OK", flush=True)
    # A remote disconnect is a transport error, not a reason to repin the
    # observed random catalogue or change its SHA256.
    with TemporaryDirectory() as tmp:
        previous_download = globals()["download"]
        expected_bytes = b"synthetic-observed-random-only"
        expected_sha = hashlib.sha256(expected_bytes).hexdigest()
        observed_tries = []
        def flaky_observed(filename, root, max_bytes, timeout):
            observed_tries.append(filename)
            if len(observed_tries) < 3:
                raise ConnectionResetError("synthetic remote disconnect")
            target = root / filename
            target.write_bytes(expected_bytes)
            return target, expected_sha, len(expected_bytes)
        globals()["download"] = flaky_observed
        try:
            observed_file = Path(tmp) / "synthetic_observed_random.fits"
            received = acquire(observed_file, expected_sha,
                               observed=True, filename=observed_file.name,
                               timeout=1)
            assert len(observed_tries) == 3
            assert received[1] == expected_sha
            assert observed_file.read_bytes() == expected_bytes
        finally:
            globals()["download"] = previous_download
    print("EBOSS_FINE_NZ_OBSERVED_TRANSPORT_RETRY_SELF_TEST_OK", flush=True)
    assert PROTOCOL.is_file()
    print("EBOSS_NINEMOCK_FINE_WEIGHTED_NZ_SELF_TEST_OK", flush=True)


def run(args):
    mock_sha, obs_sha, obs_counts, prev_files = preflight(args)
    ids = IDS if args.mock_ids == "all" else tuple(
        int(v.strip()) for v in args.mock_ids.split(","))
    caps = CAPS if args.caps == "all" else tuple(
        v.strip() for v in args.caps.split(","))
    if (not ids or len(ids) != len(set(ids)) or any(mid not in IDS for mid in ids)
            or not caps or len(caps) != len(set(caps))
            or any(c not in CAPS for c in caps)):
        raise ValueError("Run only unique preregistered IDs/caps; use 'all' for full cohort")
    out_dir = Path(args.out_dir)
    observed, cases = {}, []
    try:
        for cap in caps:
            for tracer in TRACERS:
                file_name, nrows = RANDOMS[tracer, cap]
                path = out_dir / f"observed_{cap}_{tracer}.json"
                cached = cached_record(
                    path, obs_sha[cap,tracer], cap, tracer, "observed")
                if cached is None:
                    local = Path(args.observed_cache_dir) / file_name
                    local, digest, size = acquire(
                        local, obs_sha[cap,tracer], observed=True,
                        filename=file_name, timeout=args.timeout,
                        no_download=args.no_download)
                    cached = scan_fits(
                        local, tracer=tracer, cap=cap, kind="observed",
                        digest=digest, size=size, expected_rows=nrows)
                    for iz, (lo,hi) in enumerate(zip(FINE_EDGES[::10][:-1],
                                                        FINE_EDGES[::10][1:])):
                        key = (cap,tracer,round(float(lo),1),round(float(hi),1))
                        if sum(cached["full"]["count_per_fine_bin"][iz*10:(iz+1)*10]) != obs_counts[key]:
                            raise ValueError("Observed fine n(z) differs from prior coarse input count")
                    json_write_atomic(path, cached)
                observed[cap,tracer] = cached
                print("FINE_NZ_OBS_OK",cap,tracer,cached["full"]["retained_rows"],flush=True)
        for mid in ids:
            for cap in caps:
                by_tracer = {}
                for tracer in TRACERS:
                    relative = mock_path("eBOSS_"+tracer,cap,"ran",mid)
                    name = Path(relative).name
                    path = out_dir / f"mock_{mid:04d}_{cap}_{tracer}.json"
                    cached = cached_record(
                        path,mock_sha[mid,cap,tracer],cap,tracer,"mock")
                    if cached is None:
                        local = Path(args.mock_cache_dir) / name
                        try:
                            local,digest,size = acquire(
                                local,mock_sha[mid,cap,tracer], observed=False,
                                url=MOCK_BASE+relative,timeout=args.timeout,
                                no_download=args.no_download)
                            old = prev_files.get((mid,cap,tracer))
                            cached = scan_fits(
                                local,tracer=tracer,cap=cap,kind="mock",
                                digest=digest,size=size,
                                expected_rows=old["header_rows"] if old else None)
                            if old and (
                                cached["candidate_raw_rows"] != old["candidate_random_rows"]
                                or cached["excluded_numerical_zero_systot_rows"]
                                != old["candidate_zero_weight_exclusions"]):
                                raise ValueError("Mock candidate count changed from earlier source audit")
                            json_write_atomic(path,cached)
                        finally:
                            if not args.keep_mock_cache:
                                local.unlink(missing_ok=True)
                    by_tracer[tracer] = cached
                    print("FINE_NZ_MOCK_OK",mid,cap,tracer,
                          cached["full"]["retained_rows"],flush=True)
                for tracer in TRACERS:
                    case = compare(observed[cap,tracer],by_tracer[tracer],mid=mid)
                    cases.append(case)
                    print("FINE_NZ_COMPARE",mid,cap,tracer,
                          case["full_weighted_nz_difference"][
                              "weighted_normalized_total_variation"],
                          "common_chunks",len(case["exactly_shared_elg_chunks"]),
                          flush=True)
                # Immutable-looking progress file; this is an input-only
                # status marker, never a full-ensemble completion claim.
                json_write_atomic(out_dir / "fine_weighted_nz_checkpoint.json",{
                    "status":"fine_weighted_nz_checkpoint",
                    "protocol_commit":PROTOCOL_COMMIT,
                    "completed_cases":len(cases),
                    "requested_ids":list(ids),"requested_caps":list(caps),
                    "cases":cases,
                    "observed_odd_data_vector_read":False})
                gc.collect()
    except Exception as exc:
        json_write_atomic(out_dir/"fine_weighted_nz_failure.json",{
            "status":"fine_weighted_nz_incomplete",
            "protocol_commit":PROTOCOL_COMMIT,
            "completed_cases":len(cases),
            "errors":[str(exc)],
            "observed_odd_data_vector_read":False})
        raise

    is_full = set(ids)==set(IDS) and set(caps)==set(CAPS)
    summary = {}
    if len(ids)==9:
        for cap in caps:
            for tracer in TRACERS:
                selected=[c for c in cases if c["cap"]==cap and c["tracer"]==tracer]
                if len(selected)!=9:
                    raise ValueError("Fixed full nine-mock per-tracer cohort missing")
                x=np.asarray([c["full_weighted_nz_difference"][
                    "weighted_normalized_total_variation"] for c in selected])
                summary[f"{cap}_{tracer}"]={
                    "mock_count":9,"weighted_fine_nz_TV_min":float(x.min()),
                    "weighted_fine_nz_TV_median":float(np.median(x)),
                    "weighted_fine_nz_TV_max":float(x.max()),
                    "by_fixed_id":{str(c["mock_id"]):float(v)
                                   for c,v in zip(selected,x)}}
    report={
        "status":("nine_mock_fine_weighted_nz_control_complete"
                  if is_full else "local_fine_weighted_nz_shards_complete"),
        "protocol_commit":PROTOCOL_COMMIT,
        "source_ensemble_run":SOURCE_ENSEMBLE_RUN,
        "source_mock_run":SOURCE_MOCK_RUN,
        "requested_ids":list(ids),"requested_caps":list(caps),
        "completed_comparisons":len(cases),
        "fine_z_edges":FINE_EDGES.tolist(),
        "cases":cases,"all_nine_per_cap_tracer_summary":summary,
        "errors":[],
        "observed_galaxy_data_read":False,
        "mock_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
        "exact_joint_mask_certified":False,
        "joint_mock_covariance_estimated":False,
        "physical_window_convolution_validated":False,
        "inference_protocol_frozen":False,
        "note":"Random-only normalized weighted n(z), not galaxy odd signal; ELG chunks compared only for exact shared labels. Larger full cohort, official veto mask and pair estimator remain independent requirements."
    }
    json_write_atomic(out_dir/"fine_weighted_nz_summary.json",report)
    # A successful new report supersedes any error marker from an earlier run.
    (out_dir/"fine_weighted_nz_failure.json").unlink(missing_ok=True)
    print("EBOSS_FINE_WEIGHTED_NZ",report["status"],len(cases),flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ensemble-json",type=Path)
    ap.add_argument("--mock-ids",default="1",
                    help="comma-separated fixed IDs, or all (default 1 for local pilot)")
    ap.add_argument("--caps",default="SGC",
                    help="NGC,SGC, or all (default SGC)")
    ap.add_argument("--observed-cache-dir",
                    default="eboss_workspace/local_rr/fits")
    ap.add_argument("--mock-cache-dir",
                    default="eboss_workspace/local_nz/mock_fits")
    ap.add_argument("--out-dir",default="eboss_workspace/local_nz")
    ap.add_argument("--timeout",type=float,default=120)
    ap.add_argument("--no-download",action="store_true")
    ap.add_argument("--keep-mock-cache",action="store_true")
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--preflight-only",action="store_true",
                    help="Verify pinned nine-mock ensemble and SHA manifest without FITS download")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout<=0:
        ap.error("Timeout must be positive")
    if args.preflight_only:
        mock_sha,obs_sha,obs_counts,prior = preflight(args)
        print("EBOSS_FINE_WEIGHTED_NZ_PINNED_INPUT_PREFLIGHT_OK",
              len(mock_sha),len(obs_sha),len(obs_counts),len(prior),
              flush=True)
        return 0
    run(args)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
