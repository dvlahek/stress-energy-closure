#!/usr/bin/env python3
"""Select the strongest candidate in the fixed-present-day-relic-density sweep."""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from class_scalar_fixed_omega0_run import CONTROL, control_parameters

MASS_GRID = [
    ("m003", 0.03),
    ("m006", 0.06),
    ("m010", 0.10),
    ("m018", 0.18),
    ("m030", 0.30),
    ("m060", 0.60),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--z-match", type=float, default=1100.0)
    args = ap.parse_args()

    rows = []
    for tag, mass in MASS_GRID:
        d = args.root / tag
        summary_path = d / "scalar_candidate_actual_summary.json"
        if not summary_path.exists():
            raise FileNotFoundError(summary_path)
        s = json.loads(summary_path.read_text())
        winner = s["winner"]
        strength = float(s["winner_strength_ideal_SN"])
        wm = s["moderate_precision_metrics"][winner]
        cp = control_parameters(mass, args.z_match)
        rows.append({
            "tag": tag,
            "mass_eV": mass,
            "winner": winner,
            "winner_strength_ideal_SN": strength,
            "winner_TTEE_SN": float(
                wm["TTEE_TE_combined"]["ideal_fullsky_gaussian_cv_SN"]
            ),
            "winner_phiphi_SN": float(
                wm["phiphi"]["ideal_fullsky_cv_SN_auto_only"]
            ),
            "winner_Pk_max_abs_relative_difference_percent": float(
                wm["Pk_z0"]["max_abs_relative_difference_percent"]
            ),
            "omega_cdm": cp["omega_cdm"],
            "deg_ncdm": cp["deg_ncdm"],
            "controlled_omega_ncdm_z0": cp["controlled_omega_ncdm_z0"],
        })

    best = max(rows, key=lambda r: r["winner_strength_ideal_SN"])
    src = args.root / best["tag"]
    args.outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / "winner_plus.dat", args.outdir / "winner_plus.dat")
    shutil.copy2(src / "winner_minus.dat", args.outdir / "winner_minus.dat")
    shutil.copy2(
        src / "scalar_candidate_actual_summary.json",
        args.outdir / "scalar_candidate_actual_summary.json",
    )
    if (src / "scalar_optimization_summary.json").exists():
        shutil.copy2(
            src / "scalar_optimization_summary.json",
            args.outdir / "scalar_optimization_summary.json",
        )

    winner_params = control_parameters(float(best["mass_eV"]), args.z_match)
    (args.outdir / "control_parameters.json").write_text(
        json.dumps(winner_params, indent=2) + "\n"
    )

    out = {
        "control": CONTROL,
        "z_match": args.z_match,
        "mass_grid_eV": [m for _, m in MASS_GRID],
        "ranking_metric": (
            "max(ideal full-sky joint TT+TE+EE S/N, "
            "ideal full-sky phi-phi auto S/N)"
        ),
        "per_mass": rows,
        "global_winner": best,
        "global_winner_control_parameters": winner_params,
    }
    (args.outdir / "control_mass_sweep_selection.json").write_text(
        json.dumps(out, indent=2) + "\n"
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
