#!/usr/bin/env python3
"""Synthetic-only independent closure of oriented cross-Landy-Szalay.

No public galaxy data or eBOSS odd-sector vector is opened. Weighted
D1D2, D1R2, R1D2 and R1R2 are computed by the vectorized RR helper
and independently by scalar pair loops. We compare their normalized
cross-Landy-Szalay fields and the orientation-mirror parity of odd poles.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from audit_eboss_dr16_rr_pair_closure import (
    rr_histogram, synthetic_distance, mirrored_closure,
)
from check_eboss_rr_scalar import independent_scalar_rr

LABELS = ("D1D2", "D1R2", "R1D2", "R1R2")


def synthetic_catalogue(rng, n, ra_center, dec_center):
    return (
        rng.uniform(ra_center - 2, ra_center + 2, n),
        rng.uniform(dec_center - 2, dec_center + 2, n),
        rng.uniform(0.81, 0.89, n),
        rng.uniform(0.7, 1.7, n),
    )


def cross_counts(cats, sedges, muedges, theta, implementation):
    a, b, ra, rb = cats
    inputs = {
        "D1D2": (a, b), "D1R2": (a, rb),
        "R1D2": (ra, b), "R1R2": (ra, rb),
    }
    histograms, normalizations, accepted = {}, {}, {}
    for label in LABELS:
        first, second = inputs[label]
        if implementation == "scalar":
            h, npair, norm = independent_scalar_rr(
                first, second, sedges, muedges, theta)
        elif implementation == "vectorized":
            h, info = rr_histogram(
                first, second, sedges, muedges, theta,
                synthetic_distance, block=13)
            npair, norm = info["accepted_pairs"], info["pair_normalization"]
        else:
            raise ValueError("Unrecognized independent pair implementation")
        histograms[label] = h
        normalizations[label] = float(norm)
        accepted[label] = int(npair)
    return histograms, normalizations, accepted


def cross_landy_szalay(histograms, norms):
    normalized = {k: histograms[k] / norms[k] for k in LABELS}
    rr = normalized["R1R2"]
    support = rr > 0
    xi = np.full_like(rr, np.nan)
    xi[support] = (
        normalized["D1D2"][support] - normalized["D1R2"][support]
        - normalized["R1D2"][support] + rr[support]
    ) / rr[support]
    return xi, support


def odd_synthetic_pole(xi, muedges, ell):
    if ell == 1:
        integral = (muedges[1:] ** 2 - muedges[:-1] ** 2) / 2
    elif ell == 3:
        integral = (
            5 / 8 * (muedges[1:] ** 4 - muedges[:-1] ** 4)
            - 3 / 4 * (muedges[1:] ** 2 - muedges[:-1] ** 2)
        )
    else:
        raise ValueError("Only synthetic odd poles 1 and 3 are checked")
    averaged = integral / np.diff(muedges)
    # NaN cells have no random pairs. Zero-filling is ONLY a synthetic
    # reflection test here and is not a usable observational projection.
    return np.sum(np.nan_to_num(xi, nan=0.0) * averaged[None, :], axis=1)


def audit(theta):
    rng = np.random.default_rng(99271)
    d1 = synthetic_catalogue(rng, 37, 131.5, 13.5)
    d2 = synthetic_catalogue(rng, 41, 132.0, 14.0)
    r1 = synthetic_catalogue(rng, 160, 131.6, 13.6)
    r2 = synthetic_catalogue(rng, 180, 132.1, 14.1)
    cats = (d1, d2, r1, r2)
    sedges = np.array([0., 30., 60., 90., 120., 160., 220., 340.])
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, 25)

    hist, norm, accepted = cross_counts(
        cats, sedges, muedges, theta, "vectorized")
    independent, inorm, iaccepted = cross_counts(
        cats, sedges, muedges, theta, "scalar")
    worst = 0.0
    for label in LABELS:
        np.testing.assert_allclose(
            hist[label], independent[label],
            rtol=2e-12, atol=1e-10,
            err_msg=f"Independent {label} scalar pair-count mismatch")
        assert accepted[label] == iaccepted[label]
        assert abs(norm[label] - inorm[label]) < 1e-12 * norm[label]
        worst = max(worst, float(np.max(
            np.abs(hist[label] - independent[label]))))
    xi, support = cross_landy_szalay(hist, norm)
    ixi, isupport = cross_landy_szalay(independent, inorm)
    assert np.array_equal(support, isupport)
    assert np.any(support)
    np.testing.assert_allclose(
        xi[support], ixi[support], rtol=2e-11, atol=2e-10)
    reversed_cats = (d2, d1, r2, r1)
    revhist, revnorm, revaccepted = cross_counts(
        reversed_cats, sedges, muedges, theta, "vectorized")
    for key, mirror in (
        ("D1D2", "D1D2"), ("D1R2", "R1D2"),
        ("R1D2", "D1R2"), ("R1R2", "R1R2"),
    ):
        check = mirrored_closure(
            hist[key], revhist[mirror],
            {"pair_normalization": norm[key],
             "accepted_pairs": accepted[key]},
            {"pair_normalization": revnorm[mirror],
             "accepted_pairs": revaccepted[mirror]},
            muedges,
        )
        assert check["closure_passed"], (key, mirror, check)
    reverse, reverse_support = cross_landy_szalay(revhist, revnorm)
    assert np.array_equal(support, reverse_support[:, ::-1])
    np.testing.assert_allclose(
        xi[support], reverse[:, ::-1][support],
        rtol=2e-11, atol=2e-10)
    odd_residual = {}
    for ell in (1, 3):
        a = odd_synthetic_pole(xi, muedges, ell)
        b = odd_synthetic_pole(reverse, muedges, ell)
        residual = float(np.max(np.abs(a + b)))
        assert residual < 1e-10, (ell, residual)
        odd_residual[str(ell)] = residual
    return {
        "theta_min_deg": theta,
        "raw_pair_counts": accepted,
        "supported_s_mu_cells": int(np.count_nonzero(support)),
        "independent_scalar_max_weighted_count_absolute_residual": worst,
        "odd_pole_reflection_max_absolute_residual": odd_residual,
        "four_pair_terms_checked": True,
        "forward_reverse_xi_mirror_passed": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/synthetic_cross_ls_closure.json")
    args = parser.parse_args()
    tests = [audit(theta) for theta in (0.0, 0.05, 0.5)]
    result = {
        "study": "Synthetic-only independent eBOSS-style cross-Landy-Szalay closure",
        "status": "synthetic_scalar_and_vectorized_closure_passed",
        "synthetic_seed": 99271,
        "orientation": "LRG-to-ELG; reversed ELG-to-LRG",
        "line_of_sight": "midpoint",
        "tests": tests,
        "uses_real_observed_galaxies": False,
        "uses_observed_odd_data_vector": False,
        "physical_eboss_weights_or_distance_fiducial_validated": False,
        "real_survey_window_computed": False,
        "note": (
            "This is a synthetic algorithmic parity/normalization check only. "
            "Zero-filled unsupported mu cells are not observational multipoles."
        ),
    }
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for test in tests:
        print("SYNTHETIC_CROSS_LS_OK", test["theta_min_deg"],
              "supported", test["supported_s_mu_cells"],
              "worst_count_residual",
              test["independent_scalar_max_weighted_count_absolute_residual"],
              flush=True)
    print("SYNTHETIC_CROSS_LS_AUDIT", target, result["status"], flush=True)


if __name__ == "__main__":
    main()
