#!/usr/bin/env python3
"""Build blinded conditional finite-bin eBOSS RR response in all four z bins.

Input is the successful observed-random-only 1 Mpc/h RR artifact from
audit_eboss_dr16_fine_rr_allbins.py. For each Galactic cap and candidate
redshift bin, construct matrices mapping fine-s Legendre inputs ell=0..4
to coarse observed odd outputs ell=1,3 using conditional RR weighting in
each signed-mu bin. Repeat on deterministic disjoint half-samples.

No galaxy catalogue, observed odd data vector, mock galaxy statistic,
physical wake template, covariance, or fitted nuisance parameter enters.
The output is a random-only finite-bin response diagnostic; the redshift
bins remain candidate choices until the prospective protocol is frozen.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

from audit_eboss_dr16_full_random_rr import BINS
from build_eboss_dr16_conditional_window import (
    CAPS, OUT_ELLS, IN_ELLS, THINNING, CONSTANT_CLOSURE_TOL,
    build_blocks, closure_summary, conditional_window_block,
    probe_even_to_odd,
)

EXPECTED_INPUT_STATUS = "full_four_bin_fine_RR_geometric_closure_passed"


def self_test() -> None:
    fine = np.linspace(20, 140, 121)
    coarse = np.asarray([20, 40, 60, 80, 100, 120, 140.], dtype="f8")
    mu = np.linspace(-1 - 1e-7, 1 + 1e-7, 241)
    rng = np.random.default_rng(91937)
    for iz in range(4):
        rr = rng.uniform(10, 100, (120, 240))
        rr *= 1 + (0.08 + 0.01 * iz) * (
            np.linspace(-1, 1, 240)[None, :]
            * np.linspace(-1, 1, 120)[:, None])
        blocks = build_blocks(rr, fine, coarse, mu)
        worst = max(
            v["max_constant_closure_abs_error"]
            for v in closure_summary(blocks).values())
        assert worst < CONSTANT_CLOSURE_TOL
        for o in OUT_ELLS:
            for i in IN_ELLS:
                mirrored = conditional_window_block(
                    rr[:, ::-1], fine, coarse, mu, o, i)
                parity = -1 if (o + i) % 2 else 1
                np.testing.assert_allclose(
                    mirrored, parity * blocks[(o, i)],
                    rtol=2e-11, atol=2e-11)
    print("EBOSS_ALLBIN_CONDITIONAL_WINDOW_SELF_TEST_OK", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fine-rr", help="Successful all-bin fine RR NPZ")
    ap.add_argument("--fine-audit", help="Matching all-bin fine RR JSON")
    ap.add_argument(
        "--out", default="eboss_workspace/allbin_conditional_window.json")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        self_test()
        return 0
    if not args.fine_rr or not args.fine_audit:
        ap.error("Successful all-bin fine RR NPZ and JSON are required")

    audit = json.loads(Path(args.fine_audit).read_text())
    expected_sha = os.environ.get("SOURCE_WORKFLOW_SHA")
    if (audit.get("status") != EXPECTED_INPUT_STATUS
            or audit.get("uses_observed_odd_data_vector") is not False
            or audit.get("uses_observed_galaxy_positions") is not False):
        raise ValueError("Input is not a successful blinded all-bin RR audit")
    if expected_sha and audit.get("revision_commit") != expected_sha:
        raise ValueError("All-bin RR provenance SHA does not match source run")

    data = np.load(args.fine_rr, allow_pickle=False)
    fine = np.asarray(data["fine_sedges_mpc_over_h"], dtype="f8")
    coarse = np.asarray(data["coarse_sedges_mpc_over_h"], dtype="f8")
    mu = np.asarray(data["muedges"], dtype="f8")
    zlo = np.asarray(data["zlo"], dtype="f8")
    zhi = np.asarray(data["zhi"], dtype="f8")
    if (fine.shape != (121,) or coarse.shape != (7,) or mu.shape != (241,)
            or zlo.shape != (4,) or zhi.shape != (4,)
            or not np.allclose(np.diff(fine), 1, atol=1e-12, rtol=0)):
        raise ValueError("All-bin fine RR geometry differs from declared grid")
    if not np.allclose(
        np.c_[zlo, zhi], np.asarray(BINS), atol=1e-12, rtol=0):
        raise ValueError("All-bin z-grid differs from candidate bin declaration")

    arrays = {
        "fine_sedges_mpc_over_h": fine,
        "output_sedges_mpc_over_h": coarse,
        "muedges": mu, "zlo": zlo, "zhi": zhi,
    }
    cases, errors = [], []
    for cap in CAPS:
        for iz, (lo, hi) in enumerate(BINS):
            try:
                all_matrices, raw = {}, {}
                for sample in THINNING:
                    rr = np.asarray(
                        data[f"fine_rr_{cap}_z{iz}_{sample}"], dtype="f8")
                    norm = float(np.asarray(
                        data[f"fine_norm_{cap}_z{iz}_{sample}"]).item())
                    if rr.shape != (120, 240) or norm <= 0:
                        raise ValueError("Fine RR sample shape/norm differs")
                    if np.count_nonzero(rr) != rr.size:
                        raise ValueError(
                            "Fine RR contains unsupported (s,mu) cells")
                    raw[sample] = rr
                    all_matrices[sample] = build_blocks(
                        rr, fine, coarse, mu)
                    for (o, i), matrix in all_matrices[sample].items():
                        arrays[f"M_{cap}_z{iz}_{sample}_out{o}_in{i}"] = matrix

                closure = closure_summary(all_matrices["full"])
                worst_closure = max(
                    x["max_constant_closure_abs_error"]
                    for x in closure.values())
                if worst_closure >= CONSTANT_CLOSURE_TOL:
                    raise ValueError("Constant Legendre closure failed")

                mirror_max = 0.0
                for o in OUT_ELLS:
                    for i in IN_ELLS:
                        mirror = conditional_window_block(
                            raw["full"][:, ::-1], fine, coarse, mu, o, i)
                        parity = -1 if (o + i) % 2 else 1
                        err = float(np.max(np.abs(
                            mirror - parity
                            * all_matrices["full"][(o, i)])))
                        mirror_max = max(mirror_max, err)
                if mirror_max > 1e-10:
                    raise ValueError("Signed-mu reversal parity failed")

                blocks, probes = {}, {}
                for o in OUT_ELLS:
                    diagonal = float(np.linalg.norm(
                        all_matrices["full"][(o, o)]))
                    if diagonal <= 0:
                        raise ValueError("Degenerate odd diagonal response")
                    for i in IN_ELLS:
                        key = f"{o}<-{i}"
                        full = all_matrices["full"][(o, i)]
                        blocks[key] = {
                            "full_frobenius_norm": float(
                                np.linalg.norm(full)),
                            "half_A_frobenius_difference_over_odd_diagonal":
                                float(np.linalg.norm(
                                    all_matrices["half_A"][(o, i)] - full)
                                      / diagonal),
                            "half_B_frobenius_difference_over_odd_diagonal":
                                float(np.linalg.norm(
                                    all_matrices["half_B"][(o, i)] - full)
                                      / diagonal),
                            "max_abs_half_A_cell_shift": float(np.max(np.abs(
                                all_matrices["half_A"][(o, i)] - full))),
                            "max_abs_half_B_cell_shift": float(np.max(np.abs(
                                all_matrices["half_B"][(o, i)] - full))),
                            "npz_key":
                                f"M_{cap}_z{iz}_full_out{o}_in{i}",
                        }
                        if i in (0, 2, 4):
                            reference = probe_even_to_odd(
                                all_matrices["full"], fine, coarse, o, i)
                            half_a = probe_even_to_odd(
                                all_matrices["half_A"], fine, coarse, o, i)
                            half_b = probe_even_to_odd(
                                all_matrices["half_B"], fine, coarse, o, i)
                            probes[key] = {
                                "predeclared_unitless_radial_ramp_full":
                                    reference.tolist(),
                                "max_abs_half_A_minus_full": float(
                                    np.max(np.abs(half_a - reference))),
                                "max_abs_half_B_minus_full": float(
                                    np.max(np.abs(half_b - reference))),
                            }

                cases.append({
                    "cap": cap, "bin_index": iz,
                    "zlo": lo, "zhi": hi,
                    "maximum_constant_in_s_closure_error": worst_closure,
                    "maximum_signed_mu_reversal_parity_error": mirror_max,
                    "constant_theory_closure": closure,
                    "response_blocks": blocks,
                    "synthetic_even_input_radial_ramp_probes": probes,
                })
                print(
                    "ALLBIN_CONDITIONAL_WINDOW_OK", cap,
                    f"{lo:.1f}-{hi:.1f}", "constant_max",
                    worst_closure, "mirror_max", mirror_max, flush=True)
            except (OSError, ValueError, RuntimeError, KeyError) as exc:
                errors.append(f"{cap}/z{iz}: {exc}")
                print("ALLBIN_CONDITIONAL_WINDOW_ERROR", errors[-1],
                      flush=True)

    complete = len(cases) == 8 and not errors
    report = {
        "study":
            "Blinded four-bin eBOSS conditional finite-bin RR response",
        "revision_commit": os.environ.get("GITHUB_SHA"),
        "source_fine_rr_workflow_run": os.environ.get("SOURCE_WORKFLOW_RUN"),
        "source_fine_rr_commit": audit.get("revision_commit"),
        "status": (
            "all_four_bin_conditional_finite_bin_response_validated"
            if complete else "partial"),
        "candidate_redshift_bins": [list(x) for x in BINS],
        "candidate_bins_frozen_for_inference": False,
        "caps": list(CAPS),
        "input_ells": list(IN_ELLS),
        "output_ells": list(OUT_ELLS),
        "fine_separation_step_mpc_over_h": 1.0,
        "coarse_separation_edges_mpc_over_h": coarse.tolist(),
        "mu_bins": len(mu) - 1,
        "definition": (
            "Within each signed-mu bin, average prescribed fine-s xi using "
            "RR_i,mu / sum_i RR_i,mu and then integrate exact Legendre "
            "products over the signed-mu bin."),
        "cases": cases, "errors": errors,
        "observed_galaxy_data_read": False,
        "observed_odd_data_vector_read": False,
        "mock_galaxy_data_read": False,
        "physical_wake_template_fit_or_tuning": False,
        "statistical_covariance_built": False,
        "mock_ensemble_window_validated": False,
        "scope_note": (
            "This validates the complete candidate-bin observed-random "
            "conditional finite-bin response only. It does not freeze "
            "the bins, establish a physical nuisance amplitude, validate "
            "the mock ensemble, or open the galaxy odd sector."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(out.with_suffix(".npz"), **arrays)
    print("EBOSS_ALLBIN_CONDITIONAL_WINDOW_AUDIT", report["status"], out,
          flush=True)
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
