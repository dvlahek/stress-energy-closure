#!/usr/bin/env python3
"""Random-only weighted RR pair closure for the eBOSS LRG-to-ELG estimator.

This is a *diagnostic*, not the eBOSS production survey window. It uses a
fixed-seed subsample of published random FITS catalogues, a provisional
multiplicative column-weight rule and an explicitly non-frozen distance
mapping. Swapping the tracer order must mirror the signed midpoint-LOS mu
histogram and reverse every odd raw-RR angular moment. No observed galaxy
pair, odd data vector, wake template or covariance is read or calculated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from inspect_eboss_dr16_mock_selection import (
    NUMERICAL_ZERO_WEIGHT_TOL, fetch_with_retry,
)
from inspect_eboss_dr16_selection import BASE as REAL_BASE, WEIGHTS, download

ROOT = Path(__file__).resolve().parents[1]
REAL_REF = ROOT / "source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json"
MOCK_REF = (
    ROOT / "source_data/eboss_dr16_realistic_mock_random_selection_audit_2026-09-24.json"
)
# This is the already documented DESI legacy mapping. It is used for geometric
# code closure only; it is *not* the chosen eBOSS fiducial model.
TEST_COSMO = FlatLambdaCDM(H0=67.4, Om0=0.315, Tcmb0=2.7255)
TEST_H = 0.674
DIAGNOSTIC_Z = (0.8, 0.9)
DIAGNOSTIC_SEP_EDGES = (20., 40., 60., 80., 100., 120., 140.)


def synthetic_distance(z: np.ndarray) -> np.ndarray:
    return np.asarray(z, dtype="f8") * 2800.0


def legacy_distance(z: np.ndarray) -> np.ndarray:
    return np.asarray(TEST_COSMO.comoving_distance(z).value * TEST_H, dtype="f8")


def cartesian(cat: tuple[np.ndarray, ...], distance) -> tuple[np.ndarray, ...]:
    ra, dec, z, w = cat
    lon, lat = np.deg2rad(ra), np.deg2rad(dec)
    c = np.cos(lat)
    unit = np.column_stack((c * np.cos(lon), c * np.sin(lon), np.sin(lat)))
    radial = np.asarray(distance(z), dtype="f8")
    return radial[:, None] * unit, unit, np.asarray(w, dtype="f8")


def rr_histogram(
    first: tuple[np.ndarray, ...], second: tuple[np.ndarray, ...],
    sedges: np.ndarray, muedges: np.ndarray, theta_min_deg: float,
    distance, block: int = 128,
) -> tuple[np.ndarray, dict]:
    """Weighted orientation LRG->ELG; mu=(r2^2-r1^2)/(|s||r1+r2|)."""
    pos1, dir1, w1 = cartesian(first, distance)
    pos2, dir2, w2 = cartesian(second, distance)
    r1sq = np.einsum("ij,ij->i", pos1, pos1)
    r2sq = np.einsum("ij,ij->i", pos2, pos2)
    nsep, nmu = len(sedges) - 1, len(muedges) - 1
    hist = np.zeros((nsep, nmu), dtype="f8")
    theta_cos_limit = np.cos(np.deg2rad(theta_min_deg))
    valid_pairs = 0
    for i0 in range(0, len(pos1), block):
        i1 = min(len(pos1), i0 + block)
        dot = pos1[i0:i1] @ pos2.T
        s2 = r1sq[i0:i1, None] + r2sq[None, :] - 2.0 * dot
        m2 = r1sq[i0:i1, None] + r2sq[None, :] + 2.0 * dot
        np.maximum(s2, 0.0, out=s2)
        np.maximum(m2, 0.0, out=m2)
        denom = np.sqrt(s2 * m2)
        mu = np.zeros_like(denom)
        good = denom > 0
        radial_difference = r2sq[None, :] - r1sq[i0:i1, None]
        np.divide(radial_difference, denom, out=mu, where=good)
        cos_theta = dir1[i0:i1] @ dir2.T
        sep_index = np.searchsorted(sedges, np.sqrt(s2), side="right") - 1
        mu_index = np.searchsorted(muedges, mu, side="right") - 1
        accepted = (
            good & (sep_index >= 0) & (sep_index < nsep)
            & (mu_index >= 0) & (mu_index < nmu)
            & (cos_theta <= theta_cos_limit) & (cos_theta > -1)
        )
        if not np.any(accepted):
            continue
        valid_pairs += int(np.count_nonzero(accepted))
        flat = sep_index[accepted] * nmu + mu_index[accepted]
        pair_weights = (w1[i0:i1, None] * w2[None, :])[accepted]
        hist += np.bincount(
            flat, weights=pair_weights, minlength=nsep * nmu
        ).reshape(nsep, nmu)
    norm = float(np.sum(w1, dtype="f8") * np.sum(w2, dtype="f8"))
    if norm <= 0 or not np.isfinite(norm):
        raise ValueError("The random-pair weight normalization is nonpositive")
    return hist, {
        "sampled_first": len(w1), "sampled_second": len(w2),
        "accepted_pairs": valid_pairs,
        "sum_first_weights": float(np.sum(w1, dtype="f8")),
        "sum_second_weights": float(np.sum(w2, dtype="f8")),
        "pair_normalization": norm,
    }


def odd_rr_moments(hist: np.ndarray, muedges: np.ndarray) -> dict:
    """Bin-integrated P1 and P3 moments of RR, *not* measured galaxy multipoles."""
    if hist.shape[1] != len(muedges) - 1:
        raise ValueError("mu histogram dimensions differ")
    p1_int = (muedges[1:] ** 2 - muedges[:-1] ** 2) / 2.0
    p3_int = (
        (5.0 / 8.0) * (muedges[1:] ** 4 - muedges[:-1] ** 4)
        - (3.0 / 4.0) * (muedges[1:] ** 2 - muedges[:-1] ** 2)
    )
    widths = np.diff(muedges)
    p1_avg, p3_avg = p1_int / widths, p3_int / widths
    out = []
    for row in hist:
        total = float(np.sum(row))
        out.append({
            "weighted_rr_in_sep_bin": total,
            "weighted_rr_p1_mean": (
                float(np.dot(row, p1_avg) / total) if total else None),
            "weighted_rr_p3_mean": (
                float(np.dot(row, p3_avg) / total) if total else None),
            "rr_mu_mirror_l1_over_total": (
                float(np.sum(np.abs(row - row[::-1])) / total) if total else None),
        })
    return {"per_separation_bin": out}


def mirrored_closure(
    forward: np.ndarray, reverse: np.ndarray, forward_norm: dict,
    reverse_norm: dict, muedges: np.ndarray,
) -> dict:
    if forward.shape != reverse.shape:
        raise ValueError("Forward and reverse RR histograms have different shapes")
    if not np.allclose(muedges, -muedges[::-1], atol=1e-14, rtol=0):
        raise ValueError("mu edges must be symmetric for the reversal test")
    res = forward - reverse[:, ::-1]
    abs_sum = float(np.sum(np.abs(forward)))
    max_abs = float(np.max(np.abs(res)))
    l1 = float(np.sum(np.abs(res)) / abs_sum) if abs_sum else None
    n1 = forward_norm["pair_normalization"]
    n2 = reverse_norm["pair_normalization"]
    norm_rel = float(abs(n1 - n2) / max(n1, n2))
    passed = (
        abs_sum > 0 and l1 is not None and l1 < 1e-10
        and norm_rel < 1e-12
        and forward_norm["accepted_pairs"] == reverse_norm["accepted_pairs"]
    )
    return {
        "mirror_l1_relative": l1, "mirror_max_abs": max_abs,
        "pair_normalization_relative_residual": norm_rel,
        "accepted_pair_counts_match": (
            forward_norm["accepted_pairs"] == reverse_norm["accepted_pairs"]),
        "closure_passed": passed,
    }


def select_fixed_random(
    path: Path, sample_size: int, seed: int, zlo: float, zhi: float
) -> tuple[tuple[np.ndarray, ...], dict]:
    with fits.open(path, memmap=False if path.suffix == ".gz" else True) as h:
        h.verify("exception")
        data = h[1].data
        ra = np.asarray(data["RA"], dtype="f8")
        dec = np.asarray(data["DEC"], dtype="f8")
        z = np.asarray(data["Z"], dtype="f8")
        allowed = (
            np.isfinite(ra) & np.isfinite(dec) & np.isfinite(z)
            & (ra >= 0) & (ra < 360) & (dec >= -90) & (dec <= 90)
            & (z >= zlo) & (z < zhi)
        )
        count_before_weight = int(np.count_nonzero(allowed))
        base_candidate = allowed.copy()
        nonfinite = {}
        nonpositive = {}
        numerical_zero = {}
        for name in WEIGHTS:
            w = np.asarray(data[name], dtype="f8")
            finite = np.isfinite(w)
            floor = NUMERICAL_ZERO_WEIGHT_TOL if name == "WEIGHT_SYSTOT" else 0.0
            nonfinite[name] = int(np.count_nonzero(base_candidate & ~finite))
            nonpositive[name] = int(np.count_nonzero(
                base_candidate & finite & (w <= 0)))
            numerical_zero[name] = int(np.count_nonzero(
                base_candidate & finite &
                (np.abs(w) <= NUMERICAL_ZERO_WEIGHT_TOL)))
            if np.any(base_candidate & finite &
                      (w < -NUMERICAL_ZERO_WEIGHT_TOL)):
                raise ValueError(f"Significant negative {name} in random slice")
            allowed &= finite & (w > floor)
        eligible = np.flatnonzero(allowed)
        if len(eligible) < sample_size:
            raise ValueError(
                f"Too few eligible random rows ({len(eligible)} < {sample_size})")
        rng = np.random.default_rng(seed)
        indexes = np.sort(rng.choice(eligible, size=sample_size, replace=False))
        weight = np.ones(sample_size, dtype="f8")
        for name in WEIGHTS:
            weight *= np.asarray(data[name][indexes], dtype="f8")
        if not (np.isfinite(weight).all() and np.all(weight > 0)):
            raise ValueError("Sample contains nonfinite or nonpositive pair weights")
        cat = (
            np.asarray(ra[indexes], dtype="f8").copy(),
            np.asarray(dec[indexes], dtype="f8").copy(),
            np.asarray(z[indexes], dtype="f8").copy(),
            weight,
        )
    digest = hashlib.sha256()
    for field in cat:
        digest.update(np.ascontiguousarray(field).tobytes())
    return cat, {
        "candidate_rows_before_weight": count_before_weight,
        "candidate_eligible_rows": int(len(eligible)),
        "sample_size": sample_size,
        "seed": seed, "sample_arrays_sha256": digest.hexdigest(),
        "sample_weight_min": float(np.min(weight)),
        "sample_weight_max": float(np.max(weight)),
        "candidate_nonfinite_by_weight": nonfinite,
        "candidate_nonpositive_by_weight": nonpositive,
        "candidate_near_zero_by_weight": numerical_zero,
    }


def parse_edges(txt: str) -> np.ndarray:
    edges = np.array([float(x) for x in txt.split(",")], dtype="f8")
    if len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("Separation edges must be finite and increasing")
    return edges


def self_test() -> None:
    rng = np.random.default_rng(291)
    a = (
        rng.uniform(120, 124, 35), rng.uniform(10, 15, 35),
        rng.uniform(0.82, 0.88, 35), rng.uniform(0.8, 1.5, 35))
    b = (
        rng.uniform(121, 125, 41), rng.uniform(11, 16, 41),
        rng.uniform(0.82, 0.88, 41), rng.uniform(0.75, 1.3, 41))
    sedges = np.array([0., 10., 40., 100., 200., 400.])
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, 25)
    fw, f = rr_histogram(a, b, sedges, muedges, 0.05, synthetic_distance, block=10)
    rv, r = rr_histogram(b, a, sedges, muedges, 0.05, synthetic_distance, block=9)
    chk = mirrored_closure(fw, rv, f, r, muedges)
    assert chk["closure_passed"], chk
    odd_f, odd_r = odd_rr_moments(fw, muedges), odd_rr_moments(rv, muedges)
    for frow, rrow in zip(
        odd_f["per_separation_bin"], odd_r["per_separation_bin"]
    ):
        for label in ("weighted_rr_p1_mean", "weighted_rr_p3_mean"):
            if frow[label] is not None:
                assert abs(frow[label] + rrow[label]) < 1e-12
    assert sum(x["weighted_rr_in_sep_bin"] for x in
               odd_f["per_separation_bin"]) > 0
    print("Weighted synthetic RR mirror and odd-moment closure passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/rr_pair_closure.json")
    parser.add_argument("--cache-dir", default="eboss_workspace/rr_closure_randoms")
    parser.add_argument("--sample", type=int, default=1200)
    parser.add_argument("--mock-id", type=int, default=1)
    parser.add_argument("--seed", type=int, default=4071)
    parser.add_argument("--zlo", type=float, default=DIAGNOSTIC_Z[0])
    parser.add_argument("--zhi", type=float, default=DIAGNOSTIC_Z[1])
    parser.add_argument("--sep-edges", default=",".join(str(int(e)) for e in DIAGNOSTIC_SEP_EDGES))
    parser.add_argument("--mu-bins", type=int, default=24)
    parser.add_argument("--theta-min-deg", type=float, default=0.05)
    parser.add_argument("--block", type=int, default=128)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if (args.sample < 100 or args.mu_bins < 4 or args.mu_bins % 2
            or args.zlo < 0.6 or args.zhi > 1 or args.zlo >= args.zhi
            or args.theta_min_deg < 0 or args.block < 1 or args.timeout <= 0):
        parser.error("Invalid diagnostic sample, geometry or redshift configuration")
    sep = parse_edges(args.sep_edges)
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, args.mu_bins + 1)
    real_ref = json.loads(REAL_REF.read_text(encoding="utf-8"))
    mock_ref = json.loads(MOCK_REF.read_text(encoding="utf-8"))
    if args.mock_id not in mock_ref["sample_realization_ids"]:
        parser.error("mock-id must be one of the preselected 1, 500 or 1000")
    real_checks = {(v["tracer"], v["cap"]): v for v in real_ref["randoms"]}
    mock_checks = {
        (int(v["id"]), v["tracer"], v["cap"]): v
        for v in mock_ref["mock_random_catalogues"]
    }
    reports = []
    failures = []
    root = Path(args.cache_dir)
    for cap_number, cap in enumerate(("NGC", "SGC")):
        for survey_number, survey in enumerate(("observed", "realistic_mock")):
            cats, inputs = {}, {}
            for tracer_number, tracer in enumerate(("LRG", "ELG")):
                real = survey == "observed"
                if real:
                    filename, expected_rows = RANDOMS[(tracer, cap)]
                    path = root / survey / filename
                else:
                    relative = mock_path(
                        "eBOSS_" + tracer, cap, "ran", args.mock_id)
                    path = root / survey / Path(relative).name
                try:
                    if real:
                        path, sha, nbytes = download(
                            filename, path.parent, 1024 * 1024 * 1024, args.timeout)
                        expected = real_checks[(tracer, cap)]
                        if sha != expected["sha256"] or expected_rows != expected["rows"]:
                            raise ValueError("Observed random identity differs from audited reference")
                    else:
                        sha, nbytes = fetch_with_retry(
                            MOCK_BASE + relative, path,
                            512 * 1024 * 1024, args.timeout)
                        expected = mock_checks[(args.mock_id, tracer, cap)]
                        if sha != expected["compressed_file_sha256"]:
                            raise ValueError("Mock random identity differs from audited reference")
                    cat, selected = select_fixed_random(
                        path, args.sample,
                        args.seed + 100 * cap_number + 10 * survey_number + tracer_number,
                        args.zlo, args.zhi,
                    )
                    cats[tracer] = cat
                    inputs[tracer] = {
                        "filename": path.name, "full_file_sha256": sha,
                        "full_file_bytes": nbytes, **selected,
                    }
                    print("RR_INPUT_OK", survey, cap, tracer, selected["sample_arrays_sha256"],
                          flush=True)
                except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                    failures.append(f"{survey}/{cap}/{tracer}: {exc}")
                    print("RR_INPUT_ERROR", failures[-1], flush=True)
                finally:
                    if path.exists():
                        path.unlink()
            if set(cats) != {"LRG", "ELG"}:
                continue
            try:
                f, nf = rr_histogram(
                    cats["LRG"], cats["ELG"], sep, mu,
                    args.theta_min_deg, legacy_distance, block=args.block)
                r, nr = rr_histogram(
                    cats["ELG"], cats["LRG"], sep, mu,
                    args.theta_min_deg, legacy_distance, block=args.block)
                closure = mirrored_closure(f, r, nf, nr, mu)
                odd_f = odd_rr_moments(f, mu)
                odd_r = odd_rr_moments(r, mu)
                report = {
                    "survey": survey, "cap": cap, "inputs": inputs,
                    "forward_counts": nf, "reverse_counts": nr,
                    "pair_reversal_closure": closure,
                    "random_only_forward_odd_moments": odd_f,
                    "random_only_reverse_odd_moments": odd_r,
                    "forward_raw_weighted_rr_histogram": f.tolist(),
                    "reverse_raw_weighted_rr_histogram": r.tolist(),
                }
                reports.append(report)
                print(f"RR_CLOSURE {survey} {cap} "
                      f"accepted={nf['accepted_pairs']} "
                      f"mirror_l1={closure['mirror_l1_relative']} "
                      f"pass={closure['closure_passed']}", flush=True)
                if not closure["closure_passed"]:
                    failures.append(f"{survey}/{cap}: RR orientation mirror closure failed")
            except (OSError, ValueError, RuntimeError, MemoryError) as exc:
                failures.append(f"{survey}/{cap}: {exc}")
                print("RR_CLOSURE_ERROR", failures[-1], flush=True)
    success = len(reports) == 4 and not failures
    result = {
        "study": "eBOSS LRG-ELG *random-only* weighted RR pair-count closure",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": "diagnostic_rr_pair_closure_passed" if success else "partial",
        "fixed_source_references": [
            str(REAL_REF.relative_to(ROOT)), str(MOCK_REF.relative_to(ROOT))],
        "mock_realization_id": args.mock_id,
        "sample_size_per_tracer_per_cap": args.sample,
        "seed": args.seed,
        "diagnostic_redshift_interval": [args.zlo, args.zhi],
        "diagnostic_redshift_interval_frozen_for_inference": False,
        "separation_edges_Mpc_over_h": sep.tolist(),
        "mu_edges": mu.tolist(),
        "theta_min_deg": args.theta_min_deg,
        "weights": "provisional product WEIGHT_FKP * WEIGHT_SYSTOT * WEIGHT_CP * WEIGHT_NOZ",
        "weights_validated_as_eboss_production_convention": False,
        "distance_mapping": "legacy FlatLambdaCDM H0=67.4 Om0=0.315 Tcmb0=2.7255; Mpc/h=0.674*comoving_Mpc; geometry closure only",
        "distance_mapping_validated_as_eboss_fiducial": False,
        "line_of_sight": "midpoint",
        "orientation": "first LRG then ELG; reverse is ELG then LRG",
        "reports": reports,
        "errors": failures,
        "observed_galaxy_positions_read": False,
        "observed_odd_vector_read": False,
        "full_density_RR_window_computed": False,
        "mock_covariance_computed": False,
        "note": (
            "The two orientations must mirror the signed pair mu histogram. "
            "Reported odd angular moments are raw RR geometric diagnostics, "
            "not observed galaxy multipoles or a wake signal. Fixed-seed "
            "subsampling, non-frozen cosmology and provisional weights preclude "
            "using this output as the final survey window."
        ),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("RR_PAIR_CLOSURE", result["status"], out, flush=True)
    return 0 if success else 2


if __name__ == "__main__":
    raise SystemExit(main())
