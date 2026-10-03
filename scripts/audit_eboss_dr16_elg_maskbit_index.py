#!/usr/bin/env python3
"""Offline source-index SHA gate for 19,381 official eBOSS DR16 ELG mask images.

This reads only the user's previously downloaded JSON directory index
and the committed published source/polygon provenance manifests. It
never downloads or parses a FITS image, galaxy, random or odd-sector
datum. A listed FITS file is not yet a file whose bytes are certified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
PIN = ROOT / "source_data/eboss_dr16_elg_maskbit_index_pin_2026-09-26.json"
SOURCE = ROOT / "source_data/eboss_dr16_official_mask_source_protocol_2026-09-25.json"
POLYGONS = ROOT / "source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json"
NAME_RE = re.compile(r"mask-(eboss21|eboss22|eboss23|eboss25)-[0-9]{4}[mp][0-9]{3}[.]fits[.]gz")


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def family_digest(names):
    return hashlib.sha256(
        "".join(name + chr(10) for name in names).encode("utf-8")
    ).hexdigest()


def verify_families(index, pinned, expected_total):
    if list(index) != list(pinned):
        raise ValueError("Official four-chunk family order or completeness changed")
    result, all_names = {}, set()
    for sample_id, (chunk, specification) in enumerate(pinned.items()):
        entry = index[chunk]
        names = entry["maskbit_filenames"]
        if (
            not isinstance(names, list)
            or entry["listed_maskbits"] != specification["count"]
            or len(names) != specification["count"]
            or names != sorted(names)
            or len(set(names)) != len(names)
            or specification["subsample_id"] != sample_id
            or family_digest(names) != specification["filename_sha256"]
            or not all(
                isinstance(name, str)
                and NAME_RE.fullmatch(name) is not None
                and name.startswith("mask-" + chunk + "-")
                for name in names
            )
            or names[0] != specification["first"]
            or names[len(names) // 2] != specification["middle"]
            or names[-1] != specification["last"]
        ):
            raise ValueError("Pinned official maskbit filenames differ: " + chunk)
        if all_names.intersection(names):
            raise ValueError("One official image name occurs in multiple chunks")
        all_names.update(names)
        result[chunk] = {
            "listed_images": len(names),
            "filename_list_sha256": specification["filename_sha256"],
            "subsample_id_from_pinned_upstream": sample_id,
            "first": names[0],
            "middle": names[len(names) // 2],
            "last": names[-1],
            "individual_image_file_bytes_sha256_verified": False,
        }
    if len(all_names) != expected_total:
        raise ValueError("Combined official maskbit inventory is incomplete")
    return result


def audit(p, path):
    source, polygon = load(SOURCE), load(POLYGONS)
    raw = path.read_bytes()
    report_sha = hashlib.sha256(raw).hexdigest()
    if (
        len(raw) != p["original_uploaded_report_bytes"]
        or report_sha != p["original_uploaded_report_sha256"]
        or p["original_uploaded_report_sha256"]
           != "84a3a163ec0b28a5392dadc95aa757e469cb9cb6e94cc731aac2fa2958fc6c18"
        or p["total_listed_maskbit_fits"] != 19381
        or p["observed_odd_data_vector_read"] is not False
        or p["physical_elg_mask_certified"] is not False
        or source["source_data_roots"]["nersc_mirror"] != p["published_official_source_root"]
        or source["source_data_roots"]["nersc_mirror"]
           + source["source_data_roots"]["elg_mask_directory"]
           != p["published_elg_directory"]
        or source["source_index_sha256"] != {
            "dr16_root": p["frozen_directory_index_sha256"]["root"],
            "lrg_directory": p["frozen_directory_index_sha256"]["lrg"],
            "elg_directory": p["frozen_directory_index_sha256"]["elg"],
        }
        or source["published_brickmask"] != {
            "repository": p["pinned_upstream_brickmask"]["repository"],
            "upstream_commit": p["pinned_upstream_brickmask"]["commit"],
            "extra_mask_source": p["pinned_upstream_brickmask"]["extra_script"],
            "source_blob_sha": p["pinned_upstream_brickmask"]["extra_script_blob_sha"],
            "program_compilation_flag": p["pinned_upstream_brickmask"]["compile_flag"],
        }
        or list(source["elg_maskbits_by_exact_chunk"]) != list(p["ordered_chunk_families"])
        or source["elg_extra_polygons"]
           != list(p["three_published_extra_polygons_ALREADY_PINNED_IN_SEPARATE_11_FILE_MANIFEST"])
    ):
        raise ValueError("Uploaded ELG source-index pin or parent protocol mismatch")
    index = json.loads(raw)
    if (
        index.get("status") != "official_elg_chunk_index_discovered"
        or index.get("approved_source_root") != p["published_official_source_root"]
        or index.get("approved_elg_directory") != p["published_elg_directory"]
        or index.get("source_filenames_confirmed") is not True
        or index.get("physical_joint_mask_certified") is not False
        or index.get("observed_galaxy_data_read") is not False
        or index.get("observed_odd_data_vector_read") is not False
        or index.get("errors") != []
        or {
            x: index[x + "_index"]["sha256"]
            for x in ("root", "elg", "lrg")
        } != p["frozen_directory_index_sha256"]
        or index.get("upstream_brickmask") != source["published_brickmask"]
    ):
        raise ValueError("Uploaded index status, source or published directory SHA mismatch")
    families = verify_families(
        index["elg_chunks"], p["ordered_chunk_families"],
        p["total_listed_maskbit_fits"],
    )
    extra = index["elg_extra_polygons"]
    products = [
        obj for obj in polygon["products"]
        if obj["role"] == "elg_extra_polygon"
    ]
    actual_sha = {obj["filename"]: obj["sha256"] for obj in products}
    if (
        polygon["status"] != "official_lrg_elg_polygon_bytes_pinned_only"
        or len(products) != 3
        or actual_sha != p[
            "three_published_extra_polygons_ALREADY_PINNED_IN_SEPARATE_11_FILE_MANIFEST"
        ]
        or set(extra) != set(actual_sha)
    ):
        raise ValueError("Previously pinned three ELG veto polygon sources changed")
    for name, v in extra.items():
        if (
            v.get("listed") is not True
            or v.get("source_url") != p["published_elg_directory"] + name
        ):
            raise ValueError("Uploaded index lacks published ELG extra polygon: " + name)
    return {
        "status": "OFFICIAL_ELG_FOUR_FAMILY_FILENAME_INDEX_SHA_PINNED_ONLY",
        "uploaded_index_sha256_reverified": report_sha,
        "uploaded_index_bytes": len(raw),
        "official_elg_directory_index_sha256":
            p["frozen_directory_index_sha256"]["elg"],
        "total_listed_maskbit_fits": p["total_listed_maskbit_fits"],
        "families": families,
        "previously_independently_pinned_extra_polygon_sha256": actual_sha,
        "maskbit_image_bytes_downloaded_or_verified": False,
        "observed_galaxy_or_random_rows_read": False,
        "observed_odd_data_vector_read": False,
        "physical_elg_mask_certified": False,
        "physical_lrg_elg_pair_window_certified": False,
        "new_selection_rule_applied": False,
        "note": (
            "Exact published filename inventory and previous extra-polygon "
            "SHA only. No brickmask FITS images or maskbit values inspected."
        ),
        "errors": [],
    }


def self_test(p):
    fake_spec, fake_index = {}, {}
    for i, (chunk, _) in enumerate(p["ordered_chunk_families"].items()):
        names = [
            "mask-" + chunk + "-0001m002.fits.gz",
            "mask-" + chunk + "-0001p000.fits.gz",
        ]
        fake_spec[chunk] = {
            "subsample_id": i, "count": 2,
            "filename_sha256": family_digest(names),
            "first": names[0], "middle": names[1], "last": names[-1],
        }
        fake_index[chunk] = {"listed_maskbits": 2, "maskbit_filenames": names}
    result = verify_families(fake_index, fake_spec, 8)
    assert list(result) == list(p["ordered_chunk_families"])
    for damaged in (
        dict(fake_index, eboss22={
            "listed_maskbits": 2,
            "maskbit_filenames": fake_index["eboss22"]["maskbit_filenames"][:1] * 2,
        }),
        {k: v for k, v in fake_index.items() if k != "eboss25"},
    ):
        try:
            verify_families(damaged, fake_spec, 8)
        except ValueError:
            pass
        else:
            raise AssertionError("Corrupted four-family synthetic index accepted")
    assert p["total_listed_maskbit_fits"] == 19381
    assert p["observed_odd_data_vector_read"] is False
    print("EBOSS_ELG_FOUR_FAMILY_INDEX_SYNTHETIC_SELF_TEST_OK", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--index", type=Path)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    p = load(PIN)
    if args.self_test:
        self_test(p)
        return 0
    index = args.index or ROOT / p["local_source_report"]
    out = args.out or ROOT / p["output"]
    try:
        result = audit(p, index)
    except Exception as exc:
        result = {
            "status": "OFFICIAL_ELG_FOUR_FAMILY_FILENAME_INDEX_INCOMPLETE_STOP",
            "errors": [str(exc)],
            "observed_odd_data_vector_read": False,
            "physical_elg_mask_certified": False,
            "new_selection_rule_applied": False,
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(result, indent=2) + chr(10), encoding="utf-8")
    tmp.replace(out)
    print("EBOSS_ELG_FOUR_FAMILY_INDEX", result["status"], flush=True)
    print("REPORT", out, flush=True)
    if result["errors"]:
        print("ERRORS", *result["errors"], sep=chr(10), flush=True)
        return 2
    print("TOTAL", result["total_listed_maskbit_fits"], flush=True)
    for chunk, obj in result["families"].items():
        print("CHUNK", chunk, obj["listed_images"],
              obj["filename_list_sha256"], flush=True)
    print("ODD_DATA_READ", result["observed_odd_data_vector_read"], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
