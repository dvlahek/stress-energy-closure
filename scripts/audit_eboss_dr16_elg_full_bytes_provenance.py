#!/usr/bin/env python3
"""SHA-pin the exact public DR16 ELG full FITS bytes without parsing table rows.

This is a BYTE PROVENANCE gate. Never call astropy, FITS table readers,
healpy or mask selection in this runner. Do not interpret observed data.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "source_data/eboss_dr16_elg_full_bytes_provenance_protocol_2026-09-25.json"
INDEX_PROTOCOL = ROOT / "source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json"
HEADER_PROTOCOL = ROOT / "source_data/eboss_dr16_elg_full_header_only_protocol_2026-09-25.json"
INDEX_REPORT = ROOT / "eboss_workspace/official_mask_inventory/bit8_reference_index.json"
HEADER_REPORT = ROOT / "eboss_workspace/official_mask_inventory/elg_full_header_only.json"
CHUNK = 1024 * 1024


def load_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require_frozen_prior(p):
    idx_p = load_json(INDEX_PROTOCOL)
    hdr_p = load_json(HEADER_PROTOCOL)
    idx = load_json(INDEX_REPORT)
    hdr = load_json(HEADER_REPORT)
    url = p["official_source_url"]
    if (
        p["method"].startswith("Only stream exact official FITS bytes") is False
        or url != p["official_source_url"]
        or url != idx_p["official_root"] + hdr_p["source_file"]
        or hdr_p["source_root"] != idx_p["official_root"]
        or idx.get("status") != "OFFICIAL_INDEX_MATCHED_REFERENCE_NOT_YET_VERIFIED"
        or idx.get("root_index_sha256") != idx_p["root_listing_sha256"]
        or hdr.get("status") != "ELG_FULL_HEADER_ONLY_PARSED_REFERENCE_NOT_VERIFIED"
        or hdr.get("source_url") != url
        or hdr.get("remote_bytes_reported") != p["expected_file_bytes"]
        or hdr.get("remote_etag") != p["expected_etag"]
        or hdr.get("remote_last_modified") != p["expected_last_modified"]
        or hdr.get("header_bytes_requested") != p["expected_header_bytes"]
        or hdr.get("header_range_sha256") != p["expected_header_sha256"]
        or hdr.get("header_declared_rows_NOT_READ")
             != p["expected_bintable_header_rows_NOT_READ"]
        or hdr.get("header_declared_row_bytes") != p["expected_bintable_row_bytes"]
        or hdr.get("column_count") != p["expected_column_count"]
        or hdr.get("has_mskbit_column") is not True
        or hdr.get("observed_galaxy_data_rows_read") is not False
        or hdr.get("observed_odd_data_vector_read") is not False
        or p["byte_cap"] != p["expected_file_bytes"]
        or p["full_source_sha256_NOT_YET_KNOWN"] is not True
        or p["observed_galaxy_rows_parsed"] is not False
        or p["observed_odd_data_vector_read"] is not False
    ):
        raise ValueError("Frozen official index / FITS header / provenance mismatch")
    candidate = [
        x for x in idx["candidate_files"]
        if x["filename"] == hdr_p["source_file"]
        and x["source_url"] == url
        and x["classification"] == "unverified_file_name_only"
        and x["is_directory"] is False
    ]
    if len(candidate) != 1:
        raise ValueError("Candidate not authenticated by pinned official index")


def digest_stream(source, target, *, expected_bytes, header_bytes, header_sha):
    """Bounded raw byte copy+SHA; no FITS decoding or table inspection."""
    if not 0 < header_bytes <= expected_bytes:
        raise ValueError("Invalid expected header/total size")
    full_hash = hashlib.sha256()
    first_hash = hashlib.sha256()
    count = 0
    while True:
        data = source.read(min(CHUNK, expected_bytes - count + 1))
        if not data:
            break
        count += len(data)
        if count > expected_bytes:
            raise ValueError("Full-file stream exceeded pinned exact byte count")
        target.write(data)
        full_hash.update(data)
        previous = count - len(data)
        if previous < header_bytes:
            n_header = min(len(data), header_bytes - previous)
            first_hash.update(data[:n_header])
    if count != expected_bytes:
        raise ValueError("Short official file download / local file")
    if first_hash.hexdigest() != header_sha:
        raise ValueError("Exact local first-header SHA256 differs from frozen HTTP-range header")
    return full_hash.hexdigest(), count


def inspect_remote_headers(response, p):
    if response.status != p["expected_http_status"]:
        raise ValueError("Server refused exact full download HTTP status 200")
    if response.geturl() != p["official_source_url"]:
        raise ValueError("Redirect or source URL differs from pinned official product")
    headers = response.headers
    if headers.get("Content-Encoding", "identity").lower() != "identity":
        raise ValueError("Server supplied content-encoded bytes, not original FITS file")
    if headers.get("Content-Length") != str(p["expected_file_bytes"]):
        raise ValueError("Content-Length does not match pinned FITS header metadata")
    if headers.get("ETag") != p["expected_etag"]:
        raise ValueError("Official ETag differs from independently pinned header request")
    if headers.get("Last-Modified") != p["expected_last_modified"]:
        raise ValueError("Official Last-Modified differs from pinned header request")


def hash_existing(path, p):
    with path.open("rb") as f:
        # verify size first to avoid reading an unexpectedly large local file
        if path.stat().st_size != p["expected_file_bytes"]:
            raise ValueError("Local file size differs from pinned exact official size")
        return digest_stream(
            f, io.BytesIO(), expected_bytes=p["expected_file_bytes"],
            header_bytes=p["expected_header_bytes"],
            header_sha=p["expected_header_sha256"],
        )


def acquire(p, path, timeout):
    if path.is_file():
        digest, nbytes = hash_existing(path, p)
        return digest, nbytes, "verified_existing_local_bytes"
    if path.exists():
        raise ValueError("Destination exists but is not a regular file")
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + ".part")
    if part.exists():
        raise ValueError("Stale .part download present: inspect manually; no silent overwrite")
    request = Request(
        p["official_source_url"],
        headers={
            "User-Agent": "eBOSS-blinded-ELG-raw-provenance-only/1.0",
            "Accept-Encoding": "identity",
        },
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            inspect_remote_headers(response, p)
            with part.open("xb") as sink:
                digest, nbytes = digest_stream(
                    response, sink, expected_bytes=p["expected_file_bytes"],
                    header_bytes=p["expected_header_bytes"],
                    header_sha=p["expected_header_sha256"],
                )
        part.replace(path)
    except BaseException:
        if part.is_file():
            part.unlink()
        raise
    return digest, nbytes, "official_http_full_bytes_SHA_pinned"


def self_test():
    header = b"synthetic header bytes"
    rest = b"not parsed table bytes"
    raw = header + rest
    h = hashlib.sha256(header).hexdigest()
    out = io.BytesIO()
    digest, total = digest_stream(
        io.BytesIO(raw), out, expected_bytes=len(raw),
        header_bytes=len(header), header_sha=h,
    )
    assert digest == hashlib.sha256(raw).hexdigest()
    assert total == len(raw) and out.getvalue() == raw
    for bad in (raw[:-1], raw + b"!", b"X" + raw[1:]):
        try:
            digest_stream(
                io.BytesIO(bad), io.BytesIO(), expected_bytes=len(raw),
                header_bytes=len(header), header_sha=h,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("Corrupted / short / overlong input accepted")
    print("EBOSS_ELG_RAW_BYTES_PROVENANCE_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--existing-file", type=Path,
                    help="Verify already downloaded official file; never read as FITS")
    ap.add_argument("--timeout", type=float, default=120)
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    p = load_json(PROTOCOL)
    out = ROOT / p["report_output"]
    try:
        if args.timeout <= 0:
            raise ValueError("timeout must be positive")
        require_frozen_prior(p)
        path = args.existing_file or (ROOT / p["local_output"])
        digest, nbytes, mode = (
            (lambda d, n: (d, n, "verified_existing_local_bytes"))(
                *hash_existing(path, p)
            ) if args.existing_file else acquire(p, path, args.timeout)
        )
        if out.is_file():
            earlier = load_json(out)
            if (earlier.get("status") != "OFFICIAL_ELG_FULL_BYTES_SHA_PINNED_ONLY"
                    or earlier.get("full_file_sha256") != digest
                    or earlier.get("file_bytes") != nbytes):
                raise ValueError("Previous local SHA pin differs: refuse silent repin")
        report = {
            "status": "OFFICIAL_ELG_FULL_BYTES_SHA_PINNED_ONLY",
            "official_source_url": p["official_source_url"],
            "local_file_path": str(path.resolve()),
            "acquisition_mode": mode,
            "file_bytes": nbytes,
            "full_file_sha256": digest,
            "expected_header_sha256": p["expected_header_sha256"],
            "header_bytes_sha256_verified": True,
            "full_source_sha256_frozen_in_repo": False,
            "official_pre_veto_bit8_reference_verified": False,
            "official_bit8_production_convention_resolved": False,
            "observed_galaxy_data_rows_parsed": False,
            "observed_galaxy_columns_parsed": False,
            "observed_odd_data_vector_read": False,
            "note": (
                "Full raw official FITS bytes were hashed, not decoded. "
                "Freeze this first-seen SHA256 separately before any mask-row "
                "audit; the hash is not a pre-existing SDSS publisher checksum."
            ),
            "errors": [],
        }
    except Exception as exc:
        report = {
            "status": "OFFICIAL_ELG_FULL_BYTES_PROVENANCE_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "official_pre_veto_bit8_reference_verified": False,
            "observed_galaxy_data_rows_parsed": False,
            "observed_odd_data_vector_read": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    temporary = out.with_suffix(".tmp.json")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(out)
    print("EBOSS_ELG_FULL_BYTES_PROVENANCE", report["status"], flush=True)
    print("REPORT", out, flush=True)
    if report.get("errors"):
        print("ERRORS", *report["errors"], sep="\n", flush=True)
        return 2
    print("BYTES", report["file_bytes"], flush=True)
    print("SHA256", report["full_file_sha256"], flush=True)
    print("OBSERVED_ODD_READ", report["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
