#!/usr/bin/env python3
"""Attribute NGC/SGC wake-standard degeneracy to templates vs covariance metric.

This audit intentionally never reads the observed DESI data vector.  It combines
NGC/SGC signal templates, nuisance templates and covariance metrics in all
2 x 2 x 2 combinations and reports metric cosines and retained wake norm.

For one nuisance vector v and signal w:
  rho = (w^T P v) / sqrt[(w^T P w)(v^T P v)]
  retained = sqrt(1-rho^2)

The Hartlap correction is a scalar and therefore cancels from rho/retained, but
is included in provenance for consistency.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def load_template(path, signal_col, nuisance_col):
    a = np.genfromtxt(path, delimiter=",", names=True)
    keys = np.column_stack([a["zlo"], a["zhi"], a["s_Mpc_over_h"]]).astype(float)
    w = np.asarray(a[signal_col], float)
    v = np.asarray(a[nuisance_col], float)
    if not np.all(np.isfinite(np.c_[keys, w, v])):
        raise RuntimeError(f"non-finite template values in {path}")
    return keys, w, v


def load_cov(path):
    b = np.load(path, allow_pickle=False)
    keys = np.column_stack([b["zlo"], b["zhi"], b["separation"]]).astype(float)
    C = np.asarray(b["cov_xi1"], float)
    n = int(np.asarray(b["n_mocks"]).item())
    if not np.all(np.isfinite(C)):
        raise RuntimeError(f"non-finite covariance in {path}")
    p = C.shape[0]
    h = float((n - p - 2) / (n - 1)) if n > p + 2 else None
    return keys, C, n, h


def geometry(w, v, C):
    P = np.linalg.pinv(C, rcond=1e-12)
    ww = float(w @ P @ w)
    vv = float(v @ P @ v)
    wv = float(w @ P @ v)
    rho = wv / np.sqrt(ww * vv)
    rho = float(np.clip(rho, -1.0, 1.0))
    retained = float(np.sqrt(max(0.0, 1.0 - rho * rho)))
    return {
        "metric_cosine": rho,
        "abs_metric_cosine": abs(rho),
        "retained_wake_metric_norm": retained,
        "wake_metric_norm2": ww,
        "nuisance_metric_norm2": vv,
    }


def euclidean(a, b):
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ngc-cov", required=True)
    ap.add_argument("--sgc-cov", required=True)
    ap.add_argument("--ngc-template", required=True)
    ap.add_argument("--sgc-template", required=True)
    ap.add_argument("--signal-col", default="wake_shape")
    ap.add_argument("--nuisance-col", default="standard_linked_shape")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    ktn, wn, vn = load_template(args.ngc_template, args.signal_col, args.nuisance_col)
    kts, ws, vs = load_template(args.sgc_template, args.signal_col, args.nuisance_col)
    kcn, Cn, nn, hn = load_cov(args.ngc_cov)
    kcs, Cs, ns, hs = load_cov(args.sgc_cov)

    for name, k in [("SGC template", kts), ("NGC covariance", kcn), ("SGC covariance", kcs)]:
        if not np.allclose(ktn, k, rtol=0.0, atol=1e-12):
            raise RuntimeError(f"grid mismatch: NGC template vs {name}")

    signals = {"NGC": wn, "SGC": ws}
    nuisances = {"NGC": vn, "SGC": vs}
    metrics = {"NGC": Cn, "SGC": Cs}

    combos = {}
    for sw, w in signals.items():
        for sv, v in nuisances.items():
            for sm, C in metrics.items():
                key = f"signal_{sw}__nuisance_{sv}__metric_{sm}"
                combos[key] = geometry(w, v, C)

    # Baselines reported by the regional fit.
    rho_ngc = combos["signal_NGC__nuisance_NGC__metric_NGC"]["abs_metric_cosine"]
    rho_sgc = combos["signal_SGC__nuisance_SGC__metric_SGC"]["abs_metric_cosine"]

    # Controlled swaps.
    metric_effect_on_ngc_templates = (
        combos["signal_NGC__nuisance_NGC__metric_SGC"]["abs_metric_cosine"] - rho_ngc
    )
    metric_effect_on_sgc_templates = (
        rho_sgc - combos["signal_SGC__nuisance_SGC__metric_NGC"]["abs_metric_cosine"]
    )
    template_effect_under_ngc_metric = (
        combos["signal_SGC__nuisance_SGC__metric_NGC"]["abs_metric_cosine"] - rho_ngc
    )
    template_effect_under_sgc_metric = (
        rho_sgc - combos["signal_NGC__nuisance_NGC__metric_SGC"]["abs_metric_cosine"]
    )

    out = {
        "scope": "NGC/SGC template-versus-covariance attribution of wake-standard degeneracy",
        "uses_observed_data_vector": False,
        "signal_column": args.signal_col,
        "nuisance_column": args.nuisance_col,
        "inputs": {
            "ngc_covariance": {
                "path": args.ngc_cov,
                "sha256": sha256(args.ngc_cov),
                "n_mocks": nn,
                "hartlap_factor": hn,
                "condition_number": float(np.linalg.cond(Cn)),
            },
            "sgc_covariance": {
                "path": args.sgc_cov,
                "sha256": sha256(args.sgc_cov),
                "n_mocks": ns,
                "hartlap_factor": hs,
                "condition_number": float(np.linalg.cond(Cs)),
            },
            "ngc_template": {"path": args.ngc_template, "sha256": sha256(args.ngc_template)},
            "sgc_template": {"path": args.sgc_template, "sha256": sha256(args.sgc_template)},
        },
        "template_similarity": {
            "wake_NGC_vs_SGC_euclidean_cosine": euclidean(wn, ws),
            "standard_NGC_vs_SGC_euclidean_cosine": euclidean(vn, vs),
        },
        "all_combinations": combos,
        "baseline": {
            "NGC_abs_cosine": rho_ngc,
            "NGC_retained_norm": combos["signal_NGC__nuisance_NGC__metric_NGC"]["retained_wake_metric_norm"],
            "SGC_abs_cosine": rho_sgc,
            "SGC_retained_norm": combos["signal_SGC__nuisance_SGC__metric_SGC"]["retained_wake_metric_norm"],
            "absolute_cosine_gap_SGC_minus_NGC": rho_sgc - rho_ngc,
        },
        "controlled_effects": {
            "metric_swap_effect_on_NGC_template_pair_abs_cosine": metric_effect_on_ngc_templates,
            "metric_swap_effect_on_SGC_template_pair_abs_cosine": metric_effect_on_sgc_templates,
            "template_pair_swap_effect_under_NGC_metric_abs_cosine": template_effect_under_ngc_metric,
            "template_pair_swap_effect_under_SGC_metric_abs_cosine": template_effect_under_sgc_metric,
        },
        "interpretation_guide": (
            "If metric-swap effects are much larger than template-pair-swap effects, the regional "
            "sensitivity difference is dominated by the covariance/survey-noise metric. If template "
            "effects dominate under a fixed metric, cap-specific physical nuisance shapes are the main "
            "source. Comparable terms imply a mixed origin. This audit does not by itself isolate survey "
            "window geometry from all other contributors to the regional covariance."
        ),
        "guardrail": (
            "This is a data-blind geometry diagnostic. Do not use it to select a cap, alter nuisance "
            "inputs, or retune the wake template."
        ),
    }

    Path(args.out).write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
