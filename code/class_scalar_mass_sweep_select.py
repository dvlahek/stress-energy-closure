#!/usr/bin/env python3
"""Select the strongest actual scalar/lensing hidden-state candidate across masses.

Each mass directory must contain scalar_candidate_actual_summary.json plus
winner_plus.dat and winner_minus.dat produced by class_scalar_response_optimize.py.
The script ranks masses by the same ideal full-sky strength used within each mass:
max[joint TT+TE+EE S/N, phi-phi auto S/N].
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

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
    ap.add_argument("--root", type=Path, required=True,
                    help="Directory containing m003, m006, m010, m018, m030, m060")
    ap.add_argument("--outdir", type=Path, required=True)
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
        rows.append({
            "tag": tag,
            "mass_eV": mass,
            "winner": winner,
            "winner_strength_ideal_SN": strength,
            "winner_TTEE_SN": float(wm["TTEE_TE_combined"]["ideal_fullsky_gaussian_cv_SN"]),
            "winner_phiphi_SN": float(wm["phiphi"]["ideal_fullsky_cv_SN_auto_only"]),
            "winner_Pk_max_abs_relative_difference_percent": float(
                wm["Pk_z0"]["max_abs_relative_difference_percent"]
            ),
        })

    best = max(rows, key=lambda r: r["winner_strength_ideal_SN"])
    src = args.root / best["tag"]
    args.outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src / "winner_plus.dat", args.outdir / "winner_plus.dat")
    shutil.copy2(src / "winner_minus.dat", args.outdir / "winner_minus.dat")

    # Keep the selected per-mass candidate summary under the filename expected by
    # class_scalar_response_optimize.py for the final high-precision summary.
    shutil.copy2(src / "scalar_candidate_actual_summary.json",
                 args.outdir / "scalar_candidate_actual_summary.json")
    if (src / "scalar_optimization_summary.json").exists():
        shutil.copy2(src / "scalar_optimization_summary.json",
                     args.outdir / "scalar_optimization_summary.json")

    out = {
        "mass_grid_eV": [m for _, m in MASS_GRID],
        "ranking_metric": "max(ideal full-sky joint TT+TE+EE S/N, ideal full-sky phi-phi auto S/N)",
        "per_mass": rows,
        "global_winner": best,
    }
    (args.outdir / "global_mass_sweep_selection.json").write_text(
        json.dumps(out, indent=2) + "\n"
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
