#!/usr/bin/env python3
"""eBOSS A-03E6: SOURCE-ONLY synthetic theory-to-finite-RR bridge.

The only numerical window used here is a deterministic synthetic positive RR
array. Original eBOSS E0/E1 operator NPZ, all FITS and observed odd data are
NEVER opened. A passing CI validates algebra and STOP provenance, not a
physically normalized Einstein-Vlasov prediction or eBOSS inference.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from build_eboss_dr16_conditional_window import build_blocks, conditional_window_block

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "source_data"
PROTO = SRC / "eboss_dr16_a03_e6_physical_template_to_empirical_window_bridge_protocol_2026-09-27.json"
E01PROTO = SRC / "eboss_dr16_a03_empirical_rr_window_injection_protocol_2026-09-26.json"
WINDOW_CODE = ROOT / "scripts/build_eboss_dr16_conditional_window.py"
E4REPORT = SRC / "eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json"
E4ERRATUM = SRC / "eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json"
E5 = SRC / "eboss_dr16_a03_e5_a04_physics_first_readiness_protocol_2026-09-27.json"
DESI_SHAPE = ROOT / "code/build_lrg_elg_wake_template.py"
DESI_BASIS = ROOT / "code/build_lrg_elg_physical_odd_basis.py"
DESI_3Z = ROOT / "code/build_lrg_elg_zresolved_physical_basis.py"
PROTO_BLOB = "4df20872cff5b3d8aa09831ab50b22fa5f8b6581"
E4_REPORT_SHA = "ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f"
PARENTS = {
    E01PROTO: "d8bb2d52f9c542edc26a80a23b37d639cadf64f7",
    WINDOW_CODE: "5e183e9bd6cce74743169f91c8b235c1fabbf1e1",
    E4ERRATUM: "0baf9c5a009bf8f5e78f5725e44395bad6f91715",
    E5: "0c81483bd9c70daf7cc6341f35a857761ab964d6",
    DESI_SHAPE: "5f803071d3b3de19705893fb08e621b7f8ea06db",
    DESI_BASIS: "e08382ae3eb7c701c4219d27ef00ecbbd7a2f8af",
    DESI_3Z: "181bb1d500fb44b8c864d352cecec1f2821f862d",
}
FINE = np.arange(20., 141., 1., dtype=np.float64)
COARSE = np.arange(20., 141., 20., dtype=np.float64)
MU = np.linspace(-1. - 1e-7, 1. + 1e-7, 241, dtype=np.float64)
CAPS = ("NGC", "SGC")
INS = (0, 1, 2, 3, 4)
OUTS = (1, 3)
EVEN = (0, 2, 4)
ODD = (1, 3)


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw).hexdigest()


def source_only_gate() -> dict:
    raw = PROTO.read_bytes()
    require(git_blob(raw) == PROTO_BLOB,
            "Exact preregistered E6 mathematical/source-only protocol mutated")
    p = json.loads(raw)
    for file, expected in PARENTS.items():
        require(git_blob(file.read_bytes()) == expected,
                "Frozen original E6 parent Git blob changed: " + str(file))
    e4_raw = E4REPORT.read_bytes()
    require(hashlib.sha256(e4_raw).hexdigest() == E4_REPORT_SHA,
            "Frozen original E4 user report SHA changed")
    r = json.loads(e4_raw)
    erratum = json.loads(E4ERRATUM.read_bytes())
    e5 = json.loads(E5.read_bytes())
    e01p = json.loads(E01PROTO.read_bytes())
    require(r["completed_cases"] == 18 and r["failed_cases"] == 0 and
            r["observed_odd_data_vector_read"] is False and
            r["physical_empirical_pair_window_certified"] is False and
            r["inferential_18D_covariance_computed"] is False and
            erratum["pilot_joint_dim"] == 24 and
            erratum["after_the_original_E4_trial"] is True and
            erratum["eBOSS_inferential_covariance_computed"] is False and
            e5["current_decision"] ==
            "STOP_A03_PHYSICAL_AND_A04_EBOSS_INFERENCE_OBSERVED_ODD_SEALED",
            "Original E4/E5 post-inspection STOP/covariance scope changed")
    require(
        p["date"] == "2026-09-27" and
        p["registration_boundary"].startswith("Before new E6 synthetic maths tests") and
        p["math_contract"]["input_prewindow_ells"] == list(INS) and
        p["math_contract"]["output_odd_ells"] == list(OUTS) and
        p["math_contract"]["caps"] == list(CAPS) and
        p["math_contract"]["pilot_per_cap_dimension"] == 12 and
        p["math_contract"]["pilot_joint_dimension"] == 24 and
        p["math_contract"]["final_inference_dimension"] ==
        "UNDECIDED_NOT_AUTOMATIC_18" and
        p["physics_prediction_status"] ==
        "NO_EBOSS_ABSOLUTE_PHYSICAL_WAKE_AMPLITUDE_OR_EBOSS_FINAL_P_DIMENSIONAL_OBSERVABLE_AVAILABLE",
        "Original mathematical eBOSS 24D pilot/future physical unit guards changed")
    for field in ("new_observed_rows_read", "observed_odd_data_vector_read",
                  "new_mock_gzip_downloads", "real_E4_reexecution",
                  "new_science_selection", "author_contact", "main_mutation",
                  "physical_A03_acceptance_set", "eBOSS_A04_inference_performed",
                  "unblinding_authorized"):
        require(p[field] is False, "E6 protocol scope/observed guard changed: " + field)
    expected_keys = {"fine_sedges_mpc_over_h", "output_sedges_mpc_over_h", "muedges"}
    expected_keys |= {f"observed_M_{cap}_out{o}_in{i}"
                      for cap in CAPS for o in OUTS for i in INS}
    expected_keys |= {f"mock_M_{cap}_id{mid:04d}_out{o}_in{i}"
                      for cap in CAPS for mid in (1,125,250,375,500,625,750,875,1000)
                      for o in OUTS for i in INS}
    require(len(expected_keys) == 203 and
            e01p["geometry"]["npz_expected_array_count"] == 203 and
            e01p["source_artifact_npz_member_sha256"] ==
            p["unchanged_original_empirical_operator"]["original_A03_E0E1_npz_sha256"] and
            e01p["source_artifact_npz_member_bytes"] == 250104 and
            e01p["geometry"]["npz_expected_block_shapes"] == [6, 120],
            "E6 eBOSS existing random-only operator source/key/shape identity changed")
    shape_source = DESI_SHAPE.read_text(encoding="utf-8")
    basis_source = DESI_BASIS.read_text(encoding="utf-8")
    require("absolute_wake_prediction=False" in shape_source and
            "wake/=max(" in shape_source and
            "dop/=max(" in shape_source and
            'scope": "Pre-window physical odd-sector basis for exact DESI DR1' in basis_source and
            "volume_bin_average" in basis_source,
            "Original DESI shape-only/ideal-shell provenance changed; re-audit before eBOSS transfer")
    return p


def require_future_physical_model_metadata(info: dict) -> None:
    """Syntactic STOP only; metadata cannot itself certify a true physical model."""
    require(info.get("target_survey") == "eBOSS DR16" and
            info.get("template_type") == "absolute_physical_prewindow_xi_ell" and
            info.get("dimensions") == "dimensionless_xi_ell" and
            info.get("tracer_orientation") == "LRG_to_ELG" and
            info.get("redshift_treatment") in
            ("independently_validated_effective_z_by_cap", "certified_z_conditioned_operator") and
            info.get("LRG_ELG_response_calibrated") is True and
            info.get("per_bin_max_normalization") is False and
            info.get("upstream_model_sha256") and
            info.get("physical_amplitude_provenance") and
            info.get("model_definition_frozen_before_next_physical_test") is True and
            info.get("units_and_window_response_independently_verified") is True,
            "STOP: model has DESI shape-only units, unverified eBOSS redshift, "
            "missing physical normalization/provenance or wrong tracer orientation")


def validate_synthetic_blocks(blocks: dict, *, nfine: int = 120,
                              ncoarse: int = 6) -> None:
    require(set(blocks) == {(o, i) for o in OUTS for i in INS},
            "Missing/extra conditional RR window blocks")
    for key, value in blocks.items():
        x = np.asarray(value)
        require(x.shape == (ncoarse, nfine) and x.dtype == np.dtype("float64")
                and np.isfinite(x).all(),
                "Expected complete finite float64 6x120 original-geometry window block: " +
                str(key))


def validated_synthetic_prewindow(field: dict, *, nfine: int = 120) -> dict:
    require(set(field) == set(CAPS), "Missing/extra cap in synthetic prewindow field")
    clean = {}
    for cap in CAPS:
        v = field[cap]
        require(set(v) == set(INS),
                "Missing even or odd input ell: cannot audit physically possible leakage")
        clean[cap] = {}
        for ell in INS:
            a = np.asarray(v[ell])
            require(a.shape == (nfine,) and np.issubdtype(a.dtype, np.number)
                    and np.isfinite(a).all(),
                    "Synthetic prewindow input must contain finite 120 fine-s values per ell")
            clean[cap][ell] = a.astype(np.float64, copy=False)
    return clean


def project_pilot_24d(blocks_by_cap: dict, prewindow_by_cap: dict) -> dict:
    require(set(blocks_by_cap) == set(CAPS), "Need both NGC and SGC window blocks")
    theory = validated_synthetic_prewindow(prewindow_by_cap)
    report = {}
    joint = []
    for cap in CAPS:
        blocks = blocks_by_cap[cap]
        validate_synthetic_blocks(blocks)
        report[cap] = {}
        for o in OUTS:
            intrinsic = np.sum([blocks[(o, i)] @ theory[cap][i]
                                for i in ODD], axis=0)
            leakage = np.sum([blocks[(o, i)] @ theory[cap][i]
                              for i in EVEN], axis=0)
            full = np.sum([blocks[(o, i)] @ theory[cap][i]
                           for i in INS], axis=0)
            require(full.shape == (6,) and
                    np.allclose(full, intrinsic + leakage, rtol=0, atol=2e-14),
                    "Conditional RR window has broken even/odd linear decomposition")
            report[cap][o] = {
                "intrinsic_odd": intrinsic, "even_to_odd": leakage,
                "full": full,
            }
            joint.extend(full.tolist())
    vector = np.asarray(joint, dtype=np.float64)
    require(vector.shape == (24,) and np.isfinite(vector).all(),
            "eBOSS pilot diagnostic joint length is exactly 24, not 18")
    return {"by_cap": report, "pilot_joint_24d": vector}


def synthetic_budget_bias(template: np.ndarray, delta: np.ndarray,
                          covariance: np.ndarray,
                          nuisance: np.ndarray | None = None) -> dict:
    """Algebraic demonstration only: no eBOSS covariance or physical amplitude."""
    t = np.asarray(template, dtype=np.float64)
    d = np.asarray(delta, dtype=np.float64)
    c = np.asarray(covariance, dtype=np.float64)
    require(t.ndim == d.ndim == 1 and t.shape == d.shape and
            t.size >= 2 and c.shape == (t.size, t.size) and
            np.isfinite(t).all() and np.isfinite(d).all() and
            np.isfinite(c).all() and np.allclose(c, c.T, rtol=0, atol=1e-12),
            "Synthetic noise budget requires same-shape real finite template, delta and C")
    try:
        np.linalg.cholesky(c)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Synthetic covariance not strictly positive definite") from exc
    precision_t = np.linalg.solve(c, t)
    if nuisance is None:
        tp = t
    else:
        B = np.asarray(nuisance, dtype=np.float64)
        require(B.ndim == 2 and B.shape[0] == t.size and 0 < B.shape[1] < t.size
                and np.isfinite(B).all(),
                "Nuisance columns must be finite and have fewer columns than data vector")
        PB = np.linalg.solve(c, B)
        gram = B.T @ PB
        require(np.linalg.matrix_rank(gram) == B.shape[1],
                "Nuisance columns rank deficient; no matched filter")
        tp = t - B @ np.linalg.solve(gram, B.T @ precision_t)
    Ptp = np.linalg.solve(c, tp)
    denom = float(tp @ Ptp)
    require(math.isfinite(denom) and denom > 1e-16 * float(t @ precision_t),
            "Nuisance projected physical template is zero or degenerate")
    amplitude_bias = float(Ptp @ d) / denom
    sigma_A = denom**(-.5)
    return {
        "amplitude_bias_algebraic_only": amplitude_bias,
        "sigma_A_algebraic_only": sigma_A,
        "bias_over_sigma_algebraic_only": amplitude_bias / sigma_A,
        "not_physical_or_eBOSS_inference": True,
    }


def self_test(p: dict) -> None:
    fine_centers = .5 * (FINE[1:] + FINE[:-1])
    mu_centers = .5 * (MU[1:] + MU[:-1])
    sgradient = (fine_centers - 80.) / 60.
    rr_ngc = 1000. * (1. + .25 * sgradient[:, None] * mu_centers[None, :])
    rr_sgc = 700. * (1. - .19 * sgradient[:, None] * mu_centers[None, :])
    require(rr_ngc.shape == rr_sgc.shape == (120, 240) and
            np.min(rr_ngc) > 0 and np.min(rr_sgc) > 0,
            "Synthetic RR positivity failed")
    blocks = {
        "NGC": build_blocks(rr_ngc, FINE, COARSE, MU),
        "SGC": build_blocks(rr_sgc, FINE, COARSE, MU),
    }
    for cap, rr in (("NGC", rr_ngc), ("SGC", rr_sgc)):
        validate_synthetic_blocks(blocks[cap])
        for o in OUTS:
            for i in INS:
                rev = conditional_window_block(rr[:, ::-1], FINE, COARSE, MU, o, i)
                assert np.allclose(rev, (-1)**(o+i) * blocks[cap][(o,i)],
                                   rtol=1e-10, atol=1e-10), "Signed-mu parity failed"
    zero = {cap: {ell: np.zeros(120) for ell in INS} for cap in CAPS}
    z = project_pilot_24d(blocks, zero)
    require(np.array_equal(z["pilot_joint_24d"], np.zeros(24)),
            "Zero synthetic prewindow theory must map to zero")
    flat_even = {cap: {ell: np.zeros(120) for ell in INS} for cap in CAPS}
    for cap in CAPS:
        flat_even[cap][0] = np.full(120, .10)
        flat_even[cap][2] = np.full(120, .025)
        flat_even[cap][4] = np.full(120, -.01)
    pure_even = project_pilot_24d(blocks, flat_even)
    require(np.max(np.abs(pure_even["pilot_joint_24d"])) < 3e-6,
            "Flat-in-fine-s even input spuriously created synthetic odd")
    ramp = copy.deepcopy(flat_even)
    for cap in CAPS:
        ramp[cap][0] = .10 + .08 * sgradient
    r = project_pilot_24d(blocks, ramp)
    require(np.max(np.abs(r["pilot_joint_24d"])) > 1e-8,
            "Nonseparable RR must reveal synthetic radial even-to-odd finite-bin leakage")
    odd = copy.deepcopy(zero)
    for cap in CAPS:
        odd[cap][1] = np.full(120, .02)
        odd[cap][3] = np.full(120, -.007)
    o = project_pilot_24d(blocks, odd)
    for cap in CAPS:
        require(np.max(np.abs(o["by_cap"][cap][1]["full"] - .02)) < 3e-6 and
                np.max(np.abs(o["by_cap"][cap][3]["full"] + .007)) < 3e-6,
                "Constant synthetic odd response not recovered by operator")
    total = {cap: {ell: ramp[cap][ell] + odd[cap][ell] for ell in INS}
             for cap in CAPS}
    both = project_pilot_24d(blocks, total)
    require(np.allclose(both["pilot_joint_24d"],
                        r["pilot_joint_24d"] + o["pilot_joint_24d"],
                        rtol=0, atol=2e-14),
            "Synthetic theory/window linear superposition failed")
    toy = synthetic_budget_bias(np.array([1., 2., 0., 0.]),
                                np.array([.02, .2, 0., 0.]),
                                np.diag([1., 4., 9., 16.]),
                                np.array([[1.],[0.],[0.],[0.]]))
    require(math.isclose(toy["amplitude_bias_algebraic_only"], .1,
                         rel_tol=0, abs_tol=1e-12) and
            math.isclose(toy["sigma_A_algebraic_only"], 1.,
                         rel_tol=0, abs_tol=1e-12),
            "Toy nuisance projected matched filter bias analytical control failed")
    for bad, name in (
        ({"target_survey":"DESI DR1", "template_type":"shape_only"}, "DESI shape"),
        ({"target_survey":"eBOSS DR16", "template_type":"shape_only",
          "dimensions":"dimensionless_xi_ell"}, "unphysical max-normalized eBOSS shape"),
        ({"target_survey":"eBOSS DR16", "template_type":"absolute_physical_prewindow_xi_ell",
          "dimensions":"dimensionless_xi_ell", "tracer_orientation":"LRG_to_ELG",
          "redshift_treatment":"z_eff_assumed_0.95"}, "hardcoded z effective"),
    ):
        try:
            require_future_physical_model_metadata(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("Unauthorized absolute eBOSS physical model accepted: " + name)
    incomplete = copy.deepcopy(total)
    del incomplete["NGC"][2]
    try:
        project_pilot_24d(blocks, incomplete)
    except ValueError:
        pass
    else:
        raise AssertionError("Missing even ell 2 silently accepted")
    truncated = copy.deepcopy(blocks)
    truncated["SGC"] = {(o,i):a[:,:18] for (o,i),a in blocks["SGC"].items()}
    try:
        project_pilot_24d(truncated, total)
    except ValueError:
        pass
    else:
        raise AssertionError("Wrong 18-slice eBOSS source window accepted")
    try:
        synthetic_budget_bias(np.ones(4), np.ones(4), np.zeros((4,4)))
    except ValueError:
        pass
    else:
        raise AssertionError("Singular covariance permitted in algebraic budget")
    try:
        synthetic_budget_bias(np.ones(4), np.ones(4), np.eye(4),
                              np.ones((4,1)))
    except ValueError:
        pass
    else:
        raise AssertionError("Fully nuisance-degenerate theory vector permitted")
    print("A03E6_SOURCE_ONLY_SYNTHETIC_THEORY_RR_WINDOW_BRIDGE_OK",
          "10_BLOCKS_PER_CAP", "24D_PILOT", "EVEN_ODD_DECOMPOSITION",
          "NUISANCE_BIAS_ALGEBRA", "NEGATIVE_CONTROLS_PASS", flush=True)
    print("NO_ABSOLUTE_EBOSS_PHYSICAL_TEMPLATE NO_EBOSS_COVARIANCE "
          "NO_NEW_FITS NO_OBSERVED_ODD A03_PHYSICAL_STOP A04_INFERENCE_STOP",
          flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true",
                        help="Pure source-verified positive synthetic RR and algebra only")
    args = parser.parse_args()
    p = source_only_gate()
    require(args.self_test, "No real E6 physical run authorized; use --self-test")
    self_test(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
