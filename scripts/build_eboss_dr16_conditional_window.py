#!/usr/bin/env python3
"""Data-blind, finite-bin eBOSS RR response for the LS estimator.

For each output separation bin S and signed-mu bin j, the expected
cross-Landy-Szalay field for a piecewise constant fine-s theory is the
RR-weighted fine-s average of xi(s_i,mu_j), because the random selection
cancels inside the LS ratio. Only afterwards do we project uniformly in
signed mu onto the odd output poles 1 and 3.

The resulting matrix maps prescribed fine-s ell=0..4 coefficients to
coarse odd poles. It is a conditional *finite-bin* response, not a
fully validated eBOSS survey window or an observed odd-sector statistic.
All inputs are random-only, and no wake template enters this module.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import Legendre

OUT_ELLS = (1, 3)
IN_ELLS = (0, 1, 2, 3, 4)
THINNING = ("full", "half_A", "half_B")
CAPS = ("NGC", "SGC")
CONSTANT_CLOSURE_TOL = 3e-6


def legendre_product_integrals(muedges: np.ndarray,
                               ellout: int, ellin: int) -> np.ndarray:
    poly = (Legendre.basis(ellout) * Legendre.basis(ellin)).integ()
    return np.diff(poly(muedges))


def conditional_window_block(
    rr: np.ndarray, fine: np.ndarray, coarse: np.ndarray,
    muedges: np.ndarray, ellout: int, ellin: int,
) -> np.ndarray:
    """Rows output coarse-s; columns input fine-s, exact mu-bin integrals."""
    if rr.shape != (len(fine) - 1, len(muedges) - 1):
        raise ValueError("Fine RR shape disagrees with separation/mu edges")
    if (not np.isfinite(rr).all() or np.any(rr < 0)
            or not np.all(np.diff(fine) > 0)
            or not np.all(np.diff(coarse) > 0)
            or not np.all(np.diff(muedges) > 0)):
        raise ValueError("Invalid RR counts or histogram bin edges")
    if ellout not in OUT_ELLS or ellin not in IN_ELLS:
        raise ValueError("Unexpected input or output multipole")
    mu_span = float(muedges[-1] - muedges[0])
    kernel = (
        (2 * ellout + 1) / mu_span *
        legendre_product_integrals(muedges, ellout, ellin)
    )
    matrix = np.zeros((len(coarse) - 1, len(fine) - 1), dtype="f8")
    for row, (lo, hi) in enumerate(zip(coarse[:-1], coarse[1:])):
        covered = np.flatnonzero(
            (fine[:-1] >= lo - 1e-12) &
            (fine[1:] <= hi + 1e-12)
        )
        if (not covered.size
                or abs(fine[covered[0]] - lo) > 1e-12
                or abs(fine[covered[-1] + 1] - hi) > 1e-12):
            raise ValueError("Fine separation edges do not cover a coarse bin")
        window = rr[covered]
        by_mu = np.sum(window, axis=0, dtype="f8")
        if np.any(by_mu <= 0):
            raise ValueError(
                "Missing coarse-bin mu support: cannot normalize LS estimator")
        conditional_s_given_mu = window / by_mu[None, :]
        matrix[row, covered] = conditional_s_given_mu @ kernel
    return matrix


def build_blocks(rr, fine, coarse, muedges) -> dict[tuple[int, int], np.ndarray]:
    return {
        (o, i): conditional_window_block(
            rr, fine, coarse, muedges, o, i)
        for o in OUT_ELLS for i in IN_ELLS
    }


def closure_summary(blocks: dict) -> dict:
    result = {}
    for o in OUT_ELLS:
        for i in IN_ELLS:
            row_sums = np.sum(blocks[(o, i)], axis=1, dtype="f8")
            expected = 1.0 if o == i else 0.0
            residual = float(np.max(np.abs(row_sums - expected)))
            result[f"{o}<-{i}"] = {
                "constant_in_s_theory_response": row_sums.tolist(),
                "expected_orthogonality": expected,
                "max_constant_closure_abs_error": residual,
            }
    return result


def ramp_for_coarse_bin(fine, coarse, row):
    """Predeclared, unit-free within-bin radial-gradient probe."""
    center = 0.5 * (coarse[row] + coarse[row + 1])
    width = coarse[row + 1] - coarse[row]
    s = 0.5 * (fine[1:] + fine[:-1])
    keep = (fine[:-1] >= coarse[row] - 1e-12) & (
        fine[1:] <= coarse[row + 1] + 1e-12)
    return np.where(keep, (s - center) / width, 0.0)


def probe_even_to_odd(blocks, fine, coarse, ellout, ellin):
    return np.asarray([
        np.dot(
            blocks[(ellout, ellin)][row],
            ramp_for_coarse_bin(fine, coarse, row),
        )
        for row in range(len(coarse) - 1)
    ], dtype="f8")


def unit_test() -> None:
    fine = np.array([20., 25., 30., 35., 40.])
    coarse = np.array([20., 40.])
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, 25)
    centers = 0.5 * (mu[1:] + mu[:-1])
    rng = np.random.default_rng(9221)
    # Nonseparable, nonmirror-symmetric positive RR varies with fine s.
    rr = 100 + rng.uniform(0, 20, (4, 24))
    rr += 20.0 * np.arange(4)[:, None] * centers[None, :]
    assert np.min(rr) > 0
    blocks = build_blocks(rr, fine, coarse, mu)
    for val in closure_summary(blocks).values():
        assert val["max_constant_closure_abs_error"] < CONSTANT_CLOSURE_TOL
    for o in OUT_ELLS:
        for i in IN_ELLS:
            mirrored = conditional_window_block(
                rr[:, ::-1], fine, coarse, mu, o, i)
            parity = -1 if (o + i) % 2 else 1
            np.testing.assert_allclose(
                mirrored, parity * blocks[(o, i)],
                rtol=2e-11, atol=2e-11)
    even_rr = rr + rr[:, ::-1]
    for o in OUT_ELLS:
        for i in (0, 2, 4):
            neutral = conditional_window_block(
                even_rr, fine, coarse, mu, o, i)
            assert np.max(np.abs(neutral)) < 1e-12
    # A flat-in-s monopole does not leak into the dipole even with
    # strongly asymmetric random geometry.
    assert abs(np.sum(blocks[(1, 0)])) < CONSTANT_CLOSURE_TOL
    # A radial gradient can leak under finite-s, signed-mu selection.
    gradient = probe_even_to_odd(blocks, fine, coarse, 1, 0)
    assert float(np.max(np.abs(gradient))) > 1e-6
    print("EBOSS_CONDITIONAL_RR_WINDOW_SYNTHETIC_CLOSURE_OK", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fine-rr", help="Fine random-only RR NPZ artifact")
    parser.add_argument("--fine-audit", help="Matching fine random-only RR JSON")
    parser.add_argument(
        "--out", default="eboss_workspace/conditional_rr_window_audit.json")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        unit_test()
        return 0
    if not args.fine_rr or not args.fine_audit:
        parser.error("Full fine-RR NPZ and matching JSON must be supplied")
    provenance = json.loads(Path(args.fine_audit).read_text())
    if provenance["status"] != "fine_RR_high_z_geometric_closure_passed":
        raise ValueError("The prior fine-RR pre-odd geometry gate did not pass")
    if provenance["observed_odd_data_vector_read"]:
        raise ValueError("Refusing any input carrying an observed odd data vector")
    expected_source = os.environ.get("SOURCE_WORKFLOW_SHA")
    if expected_source and provenance.get("revision_commit") != expected_source:
        raise ValueError("Fine-RR provenance does not match the source workflow commit")
    counts = np.load(args.fine_rr, allow_pickle=False)
    fine = np.asarray(counts["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(counts["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(counts["muedges"], dtype="f8")
    if (len(fine) != 121 or len(coarse) != 7 or len(mu) != 241
            or not np.allclose(np.diff(fine), 1, rtol=0, atol=1e-12)):
        raise ValueError("Fine candidate histogram geometry has changed")
    arrays = {
        "fine_sedges_mpc_over_h": fine,
        "output_sedges_mpc_over_h": coarse,
        "muedges": mu,
    }
    cases, failures = [], []
    for cap in CAPS:
        try:
            matrices = {}
            raw = {}
            for name in THINNING:
                key = f"fine_rr_{cap}_{name}"
                rr = np.asarray(counts[key], dtype="f8")
                norm = float(np.asarray(
                    counts[f"fine_norm_{cap}_{name}"]).item())
                if norm <= 0 or not np.isfinite(norm):
                    raise ValueError("Invalid full or thinned random normalization")
                raw[name] = rr
                matrices[name] = build_blocks(rr, fine, coarse, mu)
                for (o, i), m in matrices[name].items():
                    arrays[f"M_{cap}_{name}_out{o}_in{i}"] = m
            closure = closure_summary(matrices["full"])
            worst_closure = max(
                row["max_constant_closure_abs_error"]
                for row in closure.values()
            )
            if worst_closure >= CONSTANT_CLOSURE_TOL:
                raise ValueError("Constant-in-s Legendre closure failed")
            mirror_max = 0.0
            for o in OUT_ELLS:
                for i in IN_ELLS:
                    reflection = conditional_window_block(
                        raw["full"][:, ::-1],
                        fine, coarse, mu, o, i)
                    parity = -1 if (o + i) % 2 else 1
                    mismatch = float(np.max(np.abs(
                        reflection - parity * matrices["full"][(o, i)])))
                    mirror_max = max(mirror_max, mismatch)
            if mirror_max > 1e-10:
                raise ValueError("Reversed signed-mu operator parity failed")
            blocks = {}
            probes = {}
            for o in OUT_ELLS:
                denom = float(np.linalg.norm(
                    matrices["full"][(o, o)]))
                if denom <= 0:
                    raise ValueError("Degenerate odd diagonal response block")
                for i in IN_ELLS:
                    k = f"{o}<-{i}"
                    reference = matrices["full"][(o, i)]
                    differences = {
                        label: {
                            "max_absolute_matrix_cell_shift": float(
                                np.max(np.abs(
                                    matrices[label][(o, i)] - reference))),
                            "matrix_frobenius_difference_over_odd_diagonal": float(
                                np.linalg.norm(
                                    matrices[label][(o, i)] - reference) / denom),
                        }
                        for label in ("half_A", "half_B")
                    }
                    blocks[k] = {
                        "full_matrix_frobenius_norm": float(
                            np.linalg.norm(reference)),
                        "half_sample_differences": differences,
                        "npz_full_key": f"M_{cap}_full_out{o}_in{i}",
                    }
                    if i in (0, 2, 4):
                        full_probe = probe_even_to_odd(
                            matrices["full"], fine, coarse, o, i)
                        probes[k] = {
                            "predeclared_radial_ramp_full_output": full_probe.tolist(),
                            "max_abs_half_sample_probe_shift": max(
                                float(np.max(np.abs(
                                    probe_even_to_odd(
                                        matrices[label], fine, coarse, o, i)
                                    - full_probe)))
                                for label in ("half_A", "half_B")
                            ),
                        }
            cases.append({
                "cap": cap,
                "high_z_candidate_interval": [0.9, 1.0],
                "full_rr_smu_cells": int(raw["full"].size),
                "positive_full_rr_smu_cells": int(
                    np.count_nonzero(raw["full"])),
                "maximum_constant_in_s_closure_error": worst_closure,
                "maximum_signed_mu_reversal_parity_error": mirror_max,
                "constant_theory_closure": closure,
                "response_blocks": blocks,
                "synthetic_even_input_radial_ramp_probes": probes,
            })
            print("EBOSS_CONDITIONAL_RR_WINDOW_OK", cap,
                  "constant_max", worst_closure,
                  "mirror_max", mirror_max, flush=True)
        except (OSError, ValueError, RuntimeError, KeyError) as exc:
            failures.append(f"{cap}: {exc}")
            print("EBOSS_CONDITIONAL_RR_WINDOW_ERROR", failures[-1], flush=True)
    success = len(cases) == len(CAPS) and not failures
    summary = {
        "study": "Data-blind eBOSS conditional finite-bin RR response",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "status": (
            "conditional_finite_bin_response_validated"
            if success else "partial"),
        "source_artifact_run": os.environ.get("SOURCE_WORKFLOW_RUN"),
        "input_fine_rr_provenance": provenance.get("revision_commit"),
        "input_z_interval": [0.9, 1.0],
        "z_interval_frozen_for_inference": False,
        "input_ells": list(IN_ELLS), "output_ells": list(OUT_ELLS),
        "orientation": "LRG->ELG", "los": "midpoint",
        "definition": (
            "For each signed mu bin, xi is averaged over fine-s subbins "
            "using RR_i,mu/sum_fine_i RR_i,mu. It is then projected with "
            "exact Legendre-product integrals across mu bin edges."),
        "synthetic_ramp_probe": (
            "Within each coarse separation bin, xi_input(s)=(s-s_center)"
            "/coarse_bin_width; input units are arbitrary and no physical "
            "template or observed odd value is used."),
        "cases": cases, "errors": failures,
        "observed_galaxy_data_read": False,
        "observed_odd_vector_read": False,
        "wake_template_fitted_or_tuned": False,
        "full_four_z_bin_window_computed": False,
        "mock_window_and_cross_covariance_computed": False,
        "complete_physical_survey_forward_model_validated": False,
        "scope_note": (
            "A signed-mu, conditional finite-separation-bin response for "
            "the declared LS estimator. A raw RR odd moment is not an "
            "observed multipole. The full physical forward model, mock "
            "window validation and protocol remain separate gates."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    if arrays:
        np.savez_compressed(out.with_suffix(".npz"), **arrays)
    print("EBOSS_CONDITIONAL_RR_WINDOW_AUDIT", summary["status"], out,
          flush=True)
    return 0 if success else 2


if __name__ == "__main__":
    raise SystemExit(main())
