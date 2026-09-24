#!/usr/bin/env python3
"""Independent scalar benchmark for the eBOSS random-only RR histogram.

The reference deliberately does not call the vectorized Cartesian or pair
counting helpers. It performs explicit scalar pair geometry, angular cuts and
weight sums, and compares every (s, mu) bin with the separate block-vectorized
implementation. Only synthetic random positions are used.
"""
from __future__ import annotations

import math

import numpy as np

from audit_eboss_dr16_rr_pair_closure import (
    rr_histogram, mirrored_closure, synthetic_distance, odd_rr_moments,
)


def independent_scalar_rr(cat1, cat2, sedges, muedges, theta_min_deg):
    def each_position(cat):
        ra, dec, z, weight = cat
        out = []
        for lon, lat, redshift, w in zip(ra, dec, z, weight):
            lng, phi = math.radians(float(lon)), math.radians(float(lat))
            direction = np.array([
                math.cos(phi) * math.cos(lng),
                math.cos(phi) * math.sin(lng),
                math.sin(phi),
            ], dtype="f8")
            r = 2800.0 * float(redshift)
            out.append((r * direction, direction, float(w)))
        return out

    a, b = each_position(cat1), each_position(cat2)
    output = np.zeros((len(sedges) - 1, len(muedges) - 1), dtype="f8")
    count = 0
    for pos_a, dir_a, w_a in a:
        for pos_b, dir_b, w_b in b:
            theta = math.degrees(math.acos(float(np.clip(
                np.dot(dir_a, dir_b), -1, 1))))
            if theta < theta_min_deg or theta >= 180:
                continue
            difference = pos_b - pos_a
            midpoint = (pos_b + pos_a) / 2.0
            s = float(np.linalg.norm(difference))
            midpoint_norm = float(np.linalg.norm(midpoint))
            if s <= 0 or midpoint_norm <= 0:
                continue
            mu = float(np.dot(difference, midpoint) / (s * midpoint_norm))
            si = int(np.searchsorted(sedges, s, side="right") - 1)
            mi = int(np.searchsorted(muedges, mu, side="right") - 1)
            if not (0 <= si < output.shape[0] and 0 <= mi < output.shape[1]):
                continue
            output[si, mi] += w_a * w_b
            count += 1
    norm = float(sum(q[2] for q in a) * sum(q[2] for q in b))
    return output, count, norm


def main():
    rng = np.random.default_rng(72731)
    lrg = (
        rng.uniform(130.0, 136.0, 31),
        rng.uniform(10.0, 15.0, 31),
        rng.uniform(0.81, 0.89, 31),
        rng.uniform(0.8, 1.9, 31),
    )
    elg = (
        rng.uniform(131.0, 137.0, 29),
        rng.uniform(11.0, 16.0, 29),
        rng.uniform(0.80, 0.90, 29),
        rng.uniform(0.7, 1.6, 29),
    )
    sedges = np.array([0., 30., 60., 90., 120., 160., 220., 340.])
    muedges = np.linspace(-1 - 1e-7, 1 + 1e-7, 25)
    for theta_min_deg in (0.0, 0.05, 0.5):
        for a, b, direction in (
            (lrg, elg, "LRG_to_ELG"), (elg, lrg, "ELG_to_LRG")
        ):
            reference, accepted, reference_norm = independent_scalar_rr(
                a, b, sedges, muedges, theta_min_deg)
            fast, stats = rr_histogram(
                a, b, sedges, muedges, theta_min_deg, synthetic_distance,
                block=7)
            np.testing.assert_allclose(
                fast, reference, rtol=2e-12, atol=1e-10,
                err_msg=f"Independent scalar mismatch in {direction}")
            assert stats["accepted_pairs"] == accepted
            assert abs(stats["pair_normalization"] - reference_norm) <= (
                1e-12 * reference_norm)
            print(f"SCALAR_RR_OK theta={theta_min_deg} "
                  f"direction={direction} accepted={accepted}")
        forward, nf = rr_histogram(
            lrg, elg, sedges, muedges, theta_min_deg,
            synthetic_distance, block=8)
        reverse, nr = rr_histogram(
            elg, lrg, sedges, muedges, theta_min_deg,
            synthetic_distance, block=9)
        closure = mirrored_closure(forward, reverse, nf, nr, muedges)
        assert closure["closure_passed"], closure
        p1 = odd_rr_moments(forward, muedges)["per_separation_bin"]
        p1rev = odd_rr_moments(reverse, muedges)["per_separation_bin"]
        for fw, rv in zip(p1, p1rev):
            if fw["weighted_rr_p1_mean"] is not None:
                assert abs(fw["weighted_rr_p1_mean"]
                           + rv["weighted_rr_p1_mean"]) < 1e-12
            if fw["weighted_rr_p3_mean"] is not None:
                assert abs(fw["weighted_rr_p3_mean"]
                           + rv["weighted_rr_p3_mean"]) < 1e-12
    print("Independent scalar and vectorized weighted RR closure passed")


if __name__ == "__main__":
    main()
