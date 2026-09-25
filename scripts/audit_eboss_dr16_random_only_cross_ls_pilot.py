#!/usr/bin/env python3
"""Blinded high-z cross-LS algebra on disjoint splits of official random FITS.

Uses released independent DR16 LRG and ELG random selections as the
tracer-specific RR window. Fixed pseudo-D/pseudo-R subsamples are *both*
drawn without replacement from randoms, not physical galaxy mocks.
The output is solely an estimator orientation/normalization diagnostic.
Never interpret its pseudo-xi amplitudes as an observed odd signal.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from audit_eboss_dr16_9mock_fine_weighted_nz import (
    acquire, preflight, sha256_bytes, json_write_atomic,
)
from audit_eboss_dr16_rr_pair_closure import rr_histogram, mirrored_closure
from check_eboss_cross_ls_synthetic import (
    audit as synthetic_scalar_closure, cross_landy_szalay, LABELS,
)
from eboss_dr16_fiducial import (
    WEIGHT_COLUMNS, PRIMARY_GEOMETRY, comoving_mpc_over_h,
    validated_weight_product,
)
from inspect_eboss_dr16_joint_randoms import RANDOMS
from inspect_eboss_dr16_mock_headers import BASE as MOCK_BASE, mock_path
from aggregate_eboss_dr16_9mock_window import IDS, CAPS, TRACERS

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = ROOT / "source_data/eboss_dr16_random_only_cross_ls_protocol_2026-09-25.json"
EDGES_S = np.asarray([20., 40., 60., 80., 100., 120., 140.], dtype="f8")
EDGES_MU = np.linspace(-1 - 1e-7, 1 + 1e-7, 25, dtype="f8")
SPLIT_SIZE = 600
HIGH_Z = (0.9, 1.0)
SEED_BASE = 89021
ANGLE = 0.05


def distance(z):
    return comoving_mpc_over_h(z, PRIMARY_GEOMETRY)


def split_random(path, *, kind, tracer, cap, mid, expected_highz=None):
    if tracer not in TRACERS or cap not in CAPS or kind not in ("observed", "mock"):
        raise ValueError("Unsupported fixed source")
    with fits.open(path, memmap=(path.suffix != ".gz")) as hdus:
        hdus.verify("exception")
        tables = [h for h in hdus if isinstance(h, fits.BinTableHDU)]
        if len(tables) != 1:
            raise ValueError("Expected exactly one random FITS BINTABLE")
        data = tables[0].data
        if data is None or not {"RA", "DEC", "Z", *WEIGHT_COLUMNS}.issubset(
                tables[0].columns.names):
            raise ValueError("Missing validated public random fields")
        z = np.asarray(data["Z"], dtype="f8")
        candidate = (z >= HIGH_Z[0]) & (z < HIGH_Z[1])
        wfields = {
            k: np.asarray(data[k][candidate], dtype="f8")
            for k in WEIGHT_COLUMNS
        }
        w, retained = validated_weight_product(wfields)
        eligible = np.flatnonzero(candidate)[retained]
        weligible = w[retained]
        if expected_highz is not None and len(eligible) != expected_highz:
            raise ValueError("Fine n(z) checkpoint high-z count is inconsistent")
        if len(eligible) < 2 * SPLIT_SIZE:
            raise ValueError("Too few predeclared independent split random rows")
        ci, ti = CAPS.index(cap), TRACERS.index(tracer)
        seed = SEED_BASE + 10000 * ci + 100 * mid + ti
        selected = np.random.default_rng(seed).permutation(len(eligible))[
            :2 * SPLIT_SIZE]
        sets = {}
        for label, inds in (("pseudo_D", selected[:SPLIT_SIZE]),
                            ("pseudo_R", selected[SPLIT_SIZE:])):
            rows = eligible[inds]
            ra = np.asarray(data["RA"][rows], dtype="f8").copy()
            dec = np.asarray(data["DEC"][rows], dtype="f8").copy()
            zz = np.asarray(z[rows], dtype="f8").copy()
            ww = np.asarray(weligible[inds], dtype="f8").copy()
            if (not np.isfinite(ra+dec+zz+ww).all()
                    or np.any(ra < 0) or np.any(ra >= 360)
                    or np.any(dec < -90) or np.any(dec > 90)
                    or np.any(ww <= 0)):
                raise ValueError("Invalid selected high-z random values")
            sets[label] = (ra, dec, zz, ww)
        digest = hashlib.sha256()
        for key in ("pseudo_D", "pseudo_R"):
            for v in sets[key]:
                digest.update(np.ascontiguousarray(v).tobytes())
    if np.intersect1d(selected[:SPLIT_SIZE],
                      selected[SPLIT_SIZE:]).size:
        raise ValueError("Pseudo-D and pseudo-R split overlap")
    return sets, {
        "kind": kind, "cap": cap, "tracer": tracer,
        "mock_id": mid if kind == "mock" else None,
        "seed": seed, "eligible_highz_random_rows": len(eligible),
        "pseudo_D_rows": SPLIT_SIZE, "pseudo_R_rows": SPLIT_SIZE,
        "disjoint_split": True, "sample_arrays_sha256": digest.hexdigest(),
        "pseudo_D_weight_sum": float(sets["pseudo_D"][3].sum()),
        "pseudo_R_weight_sum": float(sets["pseudo_R"][3].sum()),
    }


def oriented_counts(lrg, elg):
    d_l, r_l = lrg["pseudo_D"], lrg["pseudo_R"]
    d_e, r_e = elg["pseudo_D"], elg["pseudo_R"]
    pairs = {"D1D2": (d_l,d_e), "D1R2": (d_l,r_e),
             "R1D2": (r_l,d_e), "R1R2": (r_l,r_e)}
    hist, meta = {}, {}
    for key, (first, second) in pairs.items():
        arr, info = rr_histogram(first, second, EDGES_S, EDGES_MU,
                                 ANGLE, distance, block=128)
        hist[key], meta[key] = arr, info
        if not np.isfinite(arr).all() or np.any(arr < 0):
            raise ValueError("Invalid fixed cross-pair histogram")
        independently_normalized = float(
            first[3].sum(dtype="f8") * second[3].sum(dtype="f8"))
        if abs(info["pair_normalization"] / independently_normalized - 1) > 1e-12:
            raise ValueError("Independent cross-pair normalization differs")
    return hist, meta


def algebra_case(lrg, elg, *, cap, kind, mid):
    forward, fmeta = oriented_counts(lrg, elg)
    reverse, rmeta = oriented_counts(elg, lrg)
    swap = {"D1D2":"D1D2", "D1R2":"R1D2",
            "R1D2":"D1R2", "R1R2":"R1R2"}
    residuals = {}
    for label, revname in swap.items():
        val = mirrored_closure(forward[label], reverse[revname],
            fmeta[label], rmeta[revname], EDGES_MU)
        if not val["closure_passed"]:
            raise ValueError("Independent tracer-order cross-LS pair closure failed: "
                             + label)
        residuals[label] = val
    xi, support = cross_landy_szalay(
        forward, {k:v["pair_normalization"] for k,v in fmeta.items()})
    xi_rev, support_rev = cross_landy_szalay(
        reverse, {k:v["pair_normalization"] for k,v in rmeta.items()})
    if not np.array_equal(support, support_rev[:, ::-1]):
        raise ValueError("Reversed cross-LS RR support differs")
    if not np.any(support):
        raise ValueError("No RR-supported cross-LS (s,mu) cells")
    err = float(np.max(np.abs(xi[support] - xi_rev[:, ::-1][support])))
    if not np.isfinite(err) or err >= 1e-8:
        raise ValueError("Cross-LS tracer-swap field mismatch")
    rr = forward["R1R2"]
    rr_rev = reverse["R1R2"]
    p1 = (EDGES_MU[1:]**2-EDGES_MU[:-1]**2)/(2*np.diff(EDGES_MU))
    p3 = ((5/8)*(EDGES_MU[1:]**4-EDGES_MU[:-1]**4) -
          (3/4)*(EDGES_MU[1:]**2-EDGES_MU[:-1]**2))/np.diff(EDGES_MU)
    raw_parity = {}
    for ell,p in ((1,p1),(3,p3)):
        weight = float(np.sum(rr))
        weight_rev = float(np.sum(rr_rev))
        if weight <= 0 or weight_rev <= 0:
            raise ValueError("Zero full RR histogram in fixed diagnostic bins")
        a = float(np.sum(rr @ p)/weight)
        b = float(np.sum(rr_rev @ p)/weight_rev)
        if abs(a + b) >= 1e-10:
            raise ValueError("Raw RR odd-moment swap-parity failed")
        raw_parity[str(ell)] = {
            "forward_normalized_rr_moment":a,
            "reverse_normalized_rr_moment":b,
            "reflection_sum_abs":abs(a+b),
        }
    return {
        "status":"split_random_cross_ls_algebra_checked",
        "kind":kind,"cap":cap,"mock_id":mid if kind=="mock" else None,
        "pair_terms":list(LABELS),
        "forward_reverse_pair_closure":residuals,
        "supported_s_mu_cells":int(np.count_nonzero(support)),
        "full_s_mu_cells":int(support.size),
        "reverse_cross_ls_max_abs_residual":err,
        "rr_raw_odd_parity":raw_parity,
        "pseudo_xi_not_an_observed_signal":True,
        "pseudo_D_and_R_from_same_parent_random_catalogues":True,
        "measured_galaxy_odd_data_vector_read":False,
    }


def run(args):
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    if (protocol["highz"] != [0.9, 1.0]
            or protocol["geometry"]["s_edges_mpc_h"] != EDGES_S.tolist()
            or protocol["diagnostic_sample"]["count_per_pseudo_D_or_R_per_tracer"]
            != SPLIT_SIZE or protocol["observed_odd_data_vector_read"] is not False):
        raise ValueError("Prospective random-only cross-LS protocol changed")
    if args.mock_id not in IDS or args.cap not in CAPS:
        raise ValueError("Outside preregistered cohort")
    sha_mock, sha_obs, _, _ = preflight(args)
    out = Path(args.out_dir)
    checkpoint = Path(args.fine_checkpoint_dir)
    cats = {}
    source = {}
    for kind in ("observed", "mock"):
        by_tracer, inputs = {}, {}
        for tracer in TRACERS:
            filename = (RANDOMS[tracer,args.cap][0] if kind=="observed"
                else Path(mock_path("eBOSS_"+tracer,args.cap,"ran",args.mock_id)).name)
            path = (Path(args.observed_cache_dir) if kind=="observed"
                    else Path(args.mock_cache_dir)) / filename
            sha = (sha_obs[args.cap,tracer] if kind=="observed"
                   else sha_mock[args.mock_id,args.cap,tracer])
            fine_path = checkpoint / (
                f"observed_{args.cap}_{tracer}.json" if kind=="observed"
                else f"mock_{args.mock_id:04d}_{args.cap}_{tracer}.json")
            fine = json.loads(fine_path.read_text(encoding="utf-8"))
            if (fine.get("source_sha256") != sha
                    or fine.get("status") != "fine_weighted_random_catalogue_checked"
                    or fine.get("observed_odd_data_vector_read") is not False):
                raise ValueError("Input not backed by prior fine weighted n(z) SHA checkpoint")
            expected_high = sum(fine["full"]["count_per_fine_bin"][30:])
            acquired = None
            try:
                acquired,verified_sha,size = acquire(
                    path, sha, observed=(kind=="observed"),
                    filename=filename if kind=="observed" else None,
                    url=(MOCK_BASE+mock_path("eBOSS_"+tracer,args.cap,"ran",args.mock_id)
                         if kind=="mock" else None),
                    timeout=args.timeout,no_download=args.no_download)
                pair, meta = split_random(
                    acquired,kind=kind,tracer=tracer,cap=args.cap,
                    mid=args.mock_id if kind=="mock" else 0,
                    expected_highz=expected_high)
                by_tracer[tracer] = pair
                inputs[tracer] = {
                    **meta,"source_sha256":verified_sha,
                    "full_source_file_bytes":size,
                    "fine_nz_checkpoint":str(fine_path)}
            finally:
                if kind=="mock" and acquired is not None and not args.keep_mock_cache:
                    acquired.unlink(missing_ok=True)
            print("CROSS_LS_RANDOM_INPUT_OK",kind,args.cap,tracer,
                  inputs[tracer]["eligible_highz_random_rows"],flush=True)
        cats[kind] = by_tracer
        source[kind] = inputs
    cases = []
    for kind in ("observed","mock"):
        case = algebra_case(cats[kind]["LRG"],cats[kind]["ELG"],
                            cap=args.cap,kind=kind,mid=args.mock_id)
        cases.append(case)
        print("CROSS_LS_RANDOM_PAIR_ALGEBRA_OK",kind,args.cap,
              case["supported_s_mu_cells"],
              case["reverse_cross_ls_max_abs_residual"],flush=True)
    out.mkdir(parents=True,exist_ok=True)
    report = {
        "status":"split_random_cross_ls_random_only_pilot_complete",
        "protocol_path":str(PROTOCOL_PATH.relative_to(ROOT)),
        "mock_id":args.mock_id,"cap":args.cap,
        "source_ensemble_run":36017670812,
        "cases":cases,"inputs":source,"errors":[],
        "observed_galaxy_data_read":False,
        "mock_galaxy_data_read":False,
        "observed_odd_data_vector_read":False,
        "official_lrg_elg_mask_certified":False,
        "physical_window_convolution_validated":False,
        "mock_galaxy_covariance_estimated":False,
        "not_statistical_inference":True,
        "note":"Disjoint pseudo D/R are both from the same published random catalogue. This tests cross-LS pair orientation and normalization only; it is not an independent mock-galaxy estimator closure or physical odd constraint.",
    }
    json_write_atomic(out/"random_only_cross_ls_pilot.json",report)
    print("EBOSS_RANDOM_ONLY_CROSS_LS_PILOT_OK",args.cap,args.mock_id,flush=True)


def self_test():
    for theta in (ANGLE,):
        x=synthetic_scalar_closure(theta)
        assert x["four_pair_terms_checked"] and x["forward_reverse_xi_mirror_passed"]
    # The synthetic-only test must not resolve any real FITS path.
    protocol=json.loads(PROTOCOL_PATH.read_text())
    assert protocol["observed_odd_data_vector_read"] is False
    assert tuple(protocol["mock_ids"])==IDS and tuple(protocol["caps"])==CAPS
    assert tuple(protocol["tracers"])==TRACERS
    assert np.allclose(EDGES_MU,-EDGES_MU[::-1],atol=1e-14,rtol=0)
    print("EBOSS_RANDOM_ONLY_CROSS_LS_SELF_TEST_OK",flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ensemble-json",type=Path)
    ap.add_argument("--cap",choices=CAPS,default="SGC")
    ap.add_argument("--mock-id",type=int,default=1)
    ap.add_argument("--observed-cache-dir",default="eboss_workspace/local_rr/fits")
    ap.add_argument("--mock-cache-dir",default="eboss_workspace/local_pair_window/mock_fits")
    ap.add_argument("--fine-checkpoint-dir",default="eboss_workspace/local_nz")
    ap.add_argument("--out-dir",default="eboss_workspace/local_pair_window")
    ap.add_argument("--timeout",type=float,default=120)
    ap.add_argument("--no-download",action="store_true")
    ap.add_argument("--keep-mock-cache",action="store_true")
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.timeout<=0:
        ap.error("timeout must be positive")
    if args.ensemble_json is None:
        ap.error("--ensemble-json is required")
    run(args)
    return 0


if __name__=="__main__":
    raise SystemExit(main())
