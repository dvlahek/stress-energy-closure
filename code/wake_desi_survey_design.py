#!/usr/bin/env python3
"""Survey-design map for the final nine-tracer DESI-BGS wake forecast.

The baseline hidden kinetic direction is optimized once for the validated
nine-tracer halo-mass partition using the same full N x N covariance and fixed
physical wake calibration as the production forecast.  Survey area and tracer
number density are then varied without recalibrating the wake amplitude and
without re-optimizing the hidden state.  This isolates survey-design scaling
from changes in the underlying physical signal.

Area factor changes volume at fixed number density. Density factor changes
number density at fixed volume. The reported fixed-density area threshold is
therefore an exact sqrt(area) rescaling of the validated baseline S/N.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_desi_bonvin_fisher as w
import wake_desi_multitracer_fisher as mt
import wake_desi_multitracer_fullgrid as fg
from wake_desi_split_fixedcal_fisher import make_calibration

AREA_FACTORS = [1.0, 1.25, 1.5, 2.0, 41253.0 / 14000.0]
DENSITY_FACTORS = [1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0]


def scaled_survey(survey, area, density):
    out = {}
    for z, tracers in survey.items():
        out[z] = []
        for tracer in tracers:
            row = dict(tracer)
            row["N"] = float(tracer["N"] * area * density)
            row["n"] = float(tracer["n"] * density)
            out[z].append(row)
    return out


def optimize_baseline(q, f0, weights, shapes, state0, mass, frac, survey, cal,
                      out, samples, validate, seed):
    projected = mt.linear_projected(q, f0, shapes, state0, mass, frac,
                                    survey, cal, "full")
    selected = base.candidate_pool(projected, shapes, f0, frac,
                                   samples, seed, validate)
    best = None
    validations = []
    for rank, (pred2, coeff, norm, shape) in enumerate(selected):
        fp = f0 + frac * shape
        fm = f0 - frac * shape
        pp = out / f"baseline_{rank:02d}_plus.dat"
        pm = out / f"baseline_{rank:02d}_minus.dat"
        base.write_psd(pp, q, fp)
        base.write_psd(pm, q, fm)
        statep = base.build_state(pp, mass)
        statem = base.build_state(pm, mass)
        sn_un, sn_full, perz = mt.nonlinear_projected(
            q, f0, fp, fm, state0, statep, statem, mass, survey, cal, "full"
        )
        _, sn_amp, _ = mt.nonlinear_projected(
            q, f0, fp, fm, state0, statep, statem, mass, survey, cal, "amp"
        )
        mp = cro.moments(fp, q, weights)
        mm = cro.moments(fm, q, weights)
        mismatch = np.abs(mp - mm) / np.maximum(
            0.5 * (np.abs(mp) + np.abs(mm)), 1e-300
        )
        row = {
            "rank": rank,
            "predicted_full_projected_SN": math.sqrt(max(pred2, 0.0)),
            "validated_unprojected_hidden_SN": sn_un,
            "validated_after_wake_amplitude_projection_SN": sn_amp,
            "validated_after_full_odd_nuisance_projection_SN": sn_full,
            "max_relative_moment_mismatch": float(mismatch.max()),
            "per_redshift": perz,
            "coefficients": coeff.tolist(),
        }
        validations.append(row)
        print("BASELINE", json.dumps(row))
        if best is None or sn_full > best[0]:
            best = (sn_full, row, fp, fm, shape, statep, statem)
    if best is None:
        raise RuntimeError("No baseline candidate validated")
    return best, validations


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="wake_survey_design_output")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--validate", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    refsurvey = fg.reference_geometry()
    prepared = fg.cumulative_profiles(refsurvey)
    boundaries = tuple(fg.BOUNDARIES)
    survey = mt.disjoint_tracers(boundaries, prepared)
    if len(boundaries) != 8:
        raise RuntimeError(f"Expected eight mass boundaries for nine tracers, got {boundaries}")
    if any(len(tracers) != 9 for tracers in survey.values()):
        raise RuntimeError("Nine-tracer partition construction failed")

    base.ZBINS = np.array(sorted(survey), float)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, nnull, shapes, y, moment_matrix = cro.kinetic_objects(
        q, args.mass, args.z_match
    )
    p0 = out / "fd_reference.dat"
    base.write_psd(p0, q, f0)
    state0 = base.build_state(p0, args.mass)

    w.SURVEY = list(refsurvey)
    calibration_reference = w.prepare_survey(state0)
    cal = make_calibration(q, f0, state0, args.mass, calibration_reference)

    best, validations = optimize_baseline(
        q, f0, weights, shapes, state0, args.mass, args.frac, survey, cal,
        out, args.samples, args.validate, args.seed
    )
    baseline_sn, baseline_row, fp, fm, shape, statep, statem = best
    np.savetxt(
        out / "best_pair.csv",
        np.column_stack([q, f0, fp, fm, shape]),
        delimiter=",",
        header="q,f_FD,f_plus,f_minus,normalized_delta_shape",
        comments="",
    )

    rows = []
    for area in AREA_FACTORS:
        for density in DENSITY_FACTORS:
            varied = scaled_survey(survey, area, density)
            sn_un, sn_full, _ = mt.nonlinear_projected(
                q, f0, fp, fm, state0, statep, statem,
                args.mass, varied, cal, "full"
            )
            _, sn_amp, _ = mt.nonlinear_projected(
                q, f0, fp, fm, state0, statep, statem,
                args.mass, varied, cal, "amp"
            )
            row = {
                "area_factor": float(area),
                "density_factor": float(density),
                "unprojected_SN": sn_un,
                "amplitude_projected_SN": sn_amp,
                "full_projected_SN": sn_full,
            }
            rows.append(row)
            print("DESIGN", json.dumps(row))

    baseline = next(
        row for row in rows
        if row["area_factor"] == 1.0 and row["density_factor"] == 1.0
    )
    if abs(baseline["full_projected_SN"] - baseline_sn) > 1e-10:
        raise RuntimeError("Survey-design baseline does not reproduce optimized baseline")

    lo, hi = 1.0, 20.0
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        varied = scaled_survey(survey, 1.0, mid)
        _, sn, _ = mt.nonlinear_projected(
            q, f0, fp, fm, state0, statep, statem,
            args.mass, varied, cal, "full"
        )
        if sn >= 1.0:
            hi = mid
        else:
            lo = mid
    density_1sigma = hi

    area_1sigma = (1.0 / baseline_sn) ** 2
    area_2sigma = (2.0 / baseline_sn) ** 2
    summary = {
        "mass_eV": args.mass,
        "pointwise_cap": args.frac,
        "ntracer": 9,
        "partition_boundaries": list(boundaries),
        "baseline_area_deg2": 14000.0,
        "baseline": baseline,
        "baseline_validation": baseline_row,
        "baseline_candidates": validations,
        "grid": rows,
        "density_factor_for_1sigma_at_14000deg2": density_1sigma,
        "area_factor_for_1sigma_at_fixed_density": area_1sigma,
        "area_deg2_for_1sigma_at_fixed_density": 14000.0 * area_1sigma,
        "area_factor_for_2sigma_at_fixed_density": area_2sigma,
        "full_sky_SN_at_fixed_density": baseline_sn * math.sqrt(41253.0 / 14000.0),
        "calibration_rule": (
            "Physical wake normalization is fixed once from the 13.75 two-tracer "
            "reference on the baseline k range. Survey area and density changes do "
            "not recalibrate the signal."
        ),
        "design_rule": (
            "The hidden kinetic direction is optimized once at the baseline nine-tracer "
            "configuration and then held fixed throughout the survey-design sweep."
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("SUMMARY", json.dumps(summary))


if __name__ == "__main__":
    main()
