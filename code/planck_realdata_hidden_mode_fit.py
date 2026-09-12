#!/usr/bin/env python3
"""Blind Planck-2018 Plik-lite fit of the pre-defined hidden kinetic mode.

The hidden direction, mass and fixed-present-day-abundance control are frozen
before looking at Planck data. The calculation uses the already-computed
high-precision CLASS spectra from the survey-validation campaign and the
public Planck 2018 Plik-lite TTTEEE data/covariance via planck-lite-py.

The theory response is linearized about the midpoint hidden state in the six
LambdaCDM directions already used by the survey forecast. We additionally
profile over the Planck calibration amplitude. A Gaussian tau prior with
sigma=0.007 and a calibration prior sigma=0.0025 are included. Foreground
and other high-l nuisance parameters are already marginalized in Plik-lite.

alpha = +1 and -1 denote the existing +/-30% pointwise-capped hidden-state
pair. The physically validated family is therefore |alpha| <= 1. We also
report the formal unconstrained Gaussian template amplitude and uncertainty;
values outside |alpha|<=1 are only a local template extrapolation, not a
positive kinetic distribution claim.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import scipy.linalg

TCMB_K = 2.7255
UK2_PER_DIMENSIONLESS = (TCMB_K * 1.0e6) ** 2
TAU_PRIOR_SIGMA = 0.007
CAL_PRIOR_SIGMA = 0.0025
PLANCK_LITE_COMMIT = "2c0d0f67e59ce781654cf62dd7fb10757b0e60be"

DERIVATIVE_STEPS = {
    "H0": 0.50,
    "omega_b": 1.5e-4,
    "omega_cdm": 8.0e-4,
    "lnAs": 1.0e-2,
    "n_s": 5.0e-3,
    "tau_reio": 4.0e-3,
}


def load_planck_lite(root: Path):
    module_path = root / "planck_lite_py.py"
    spec = importlib.util.spec_from_file_location("planck_lite_py_external", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {module_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.PlanckLitePy(
        data_directory=str(root / "data"), year=2018,
        spectra="TTTEEE", use_low_ell_bins=False,
    )


def read_class_lensed(path: Path) -> dict[str, np.ndarray]:
    header = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                break
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr[None, :]
    import re
    labels = {}
    for line in header:
        for num, name in re.findall(r"(\d+):([^\s]+)", line):
            labels[name.strip().lower()] = int(num) - 1

    def col(name: str) -> int:
        if name not in labels:
            raise RuntimeError(f"Missing {name} in {path}; labels={labels}")
        return labels[name]

    ell = arr[:, 0].astype(int)
    if ell[0] > 2:
        raise RuntimeError(f"Expected CLASS spectra from ell<=2, got ell_min={ell[0]}")
    return {
        "ell": ell,
        "TT": arr[:, col("tt")] * UK2_PER_DIMENSIONLESS,
        "TE": arr[:, col("te")] * UK2_PER_DIMENSIONLESS,
        "EE": arr[:, col("ee")] * UK2_PER_DIMENSIONLESS,
    }


def unique_spectrum(root: Path, stem: str) -> Path:
    hits = list(root.rglob(f"{stem}_00_cl_lensed.dat"))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one lensed spectrum for {stem}, found {hits}")
    return hits[0]


def align(spectra: list[dict[str, np.ndarray]]) -> None:
    ref = spectra[0]["ell"]
    for x in spectra[1:]:
        if not np.array_equal(ref, x["ell"]):
            raise RuntimeError("CLASS ell grids differ")


def valid_bin_indices(like, ell_max: int) -> tuple[list[int], dict[str, int]]:
    """Return Plik-lite data-vector indices fully covered by available theory."""
    keep: list[int] = []
    ntt = nte = nee = 0
    for i in range(like.nbintt):
        if int(like.blmax_TT[i] + like.plmin_TT) <= ell_max:
            keep.append(i); ntt += 1
    te0 = like.nbintt
    for i in range(like.nbinte):
        if int(like.blmax[i] + like.plmin) <= ell_max:
            keep.append(te0 + i); nte += 1
    ee0 = like.nbintt + like.nbinte
    for i in range(like.nbinee):
        if int(like.blmax[i] + like.plmin) <= ell_max:
            keep.append(ee0 + i); nee += 1
    return keep, {"TT": ntt, "TE": nte, "EE": nee}


def binned_model(like, sp: dict[str, np.ndarray], keep: list[int]) -> np.ndarray:
    ell = sp["ell"]
    ellmin, ellmax = int(ell[0]), int(ell[-1])
    ls = ell.astype(float)
    fac = ls * (ls + 1.0) / (2.0 * math.pi)
    fac = np.where(fac > 0, fac, np.inf)
    Cltt, Clte, Clee = sp["TT"] / fac, sp["TE"] / fac, sp["EE"] / fac
    full = np.full(like.nbin_tot, np.nan, dtype=float)

    for i in range(like.nbintt):
        lo = int(like.blmin_TT[i] + like.plmin_TT)
        hi = int(like.blmax_TT[i] + like.plmin_TT)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w_TT[like.blmin_TT[i]:like.blmax_TT[i] + 1]
        full[i] = float(np.sum(Cltt[sl] * w))

    off = like.nbintt
    for i in range(like.nbinte):
        lo = int(like.blmin[i] + like.plmin)
        hi = int(like.blmax[i] + like.plmin)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w[like.blmin[i]:like.blmax[i] + 1]
        full[off + i] = float(np.sum(Clte[sl] * w))

    off = like.nbintt + like.nbinte
    for i in range(like.nbinee):
        lo = int(like.blmin[i] + like.plmin)
        hi = int(like.blmax[i] + like.plmin)
        if hi > ellmax or lo < ellmin: continue
        sl = slice(lo - ellmin, hi + 1 - ellmin)
        w = like.bin_w[like.blmin[i]:like.blmax[i] + 1]
        full[off + i] = float(np.sum(Clee[sl] * w))

    vec = full[np.asarray(keep, dtype=int)]
    if not np.all(np.isfinite(vec)):
        raise RuntimeError("Non-finite binned theory vector")
    return vec


def chi2(beta, r0, A, precision, prior_precision) -> float:
    r = r0 - A @ beta
    return float(r @ precision @ r + beta @ prior_precision @ beta)


def solve_profile(r0, a_alpha, B, precision, prior_nuis, alpha):
    rr = r0 - alpha * a_alpha
    N = B.T @ precision @ B + prior_nuis
    g = B.T @ precision @ rr
    nuisance = scipy.linalg.solve(N, g, assume_a="sym")
    resid = rr - B @ nuisance
    val = float(resid @ precision @ resid + nuisance @ prior_nuis @ nuisance)
    return val, nuisance


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectra-root", type=Path, required=True)
    ap.add_argument("--planck-lite-root", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    like = load_planck_lite(args.planck_lite_root)
    stems = ["fiducial", "winner_plus", "winner_minus"]
    for p in DERIVATIVE_STEPS:
        stems.extend([f"d_{p}_plus", f"d_{p}_minus"])
    spec = {s: read_class_lensed(unique_spectrum(args.spectra_root, s)) for s in stems}
    align(list(spec.values()))
    ell_max = int(spec["fiducial"]["ell"][-1])

    keep, bin_counts = valid_bin_indices(like, ell_max)
    if bin_counts["TE"] != like.nbinte or bin_counts["EE"] != like.nbinee:
        raise RuntimeError(f"Unexpected polarization truncation: {bin_counts}")
    if bin_counts["TT"] < like.nbintt - 3:
        raise RuntimeError(f"Too many TT bins lost at ell_max={ell_max}: {bin_counts}")

    cov_full = scipy.linalg.inv(like.fisher)
    idx = np.asarray(keep, dtype=int)
    cov = cov_full[np.ix_(idx, idx)]
    precision = scipy.linalg.cho_solve(scipy.linalg.cho_factor(cov), np.eye(cov.shape[0]))
    data = np.asarray(like.X_data, dtype=float)[idx]

    xb = {s: binned_model(like, v, keep) for s, v in spec.items()}
    x0 = xb["fiducial"]
    r0 = data - x0

    cols, names, units = [], [], {}
    cols.append(0.5 * (xb["winner_plus"] - xb["winner_minus"]))
    names.append("alpha_hidden")
    units["alpha_hidden"] = "1 corresponds to the +30% hidden-state endpoint"
    for p, step in DERIVATIVE_STEPS.items():
        cols.append(0.5 * (xb[f"d_{p}_plus"] - xb[f"d_{p}_minus"]))
        names.append(p)
        units[p] = f"1 fit unit = {step:g} physical parameter shift"

    cplus = x0 / (1.0 + CAL_PRIOR_SIGMA) ** 2
    cminus = x0 / (1.0 - CAL_PRIOR_SIGMA) ** 2
    cols.append(0.5 * (cplus - cminus))
    names.append("A_planck")
    units["A_planck"] = f"1 fit unit = {CAL_PRIOR_SIGMA:g} calibration shift"

    A = np.column_stack(cols)
    p = A.shape[1]
    prior_precision = np.zeros((p, p), dtype=float)
    tau_i, cal_i = names.index("tau_reio"), names.index("A_planck")
    prior_precision[tau_i, tau_i] = (DERIVATIVE_STEPS["tau_reio"] / TAU_PRIOR_SIGMA) ** 2
    prior_precision[cal_i, cal_i] = 1.0

    N = A.T @ precision @ A + prior_precision
    g = A.T @ precision @ r0
    cov_beta = scipy.linalg.inv(N)
    beta = cov_beta @ g
    sigma = np.sqrt(np.clip(np.diag(cov_beta), 0.0, None))
    chi2_min = chi2(beta, r0, A, precision, prior_precision)

    a_alpha, B = A[:, 0], A[:, 1:]
    prior_nuis = prior_precision[1:, 1:]
    chi2_0, _ = solve_profile(r0, a_alpha, B, precision, prior_nuis, 0.0)
    delta_chi2_0 = max(chi2_0 - chi2_min, 0.0)
    alpha_hat, sigma_alpha = float(beta[0]), float(sigma[0])
    alpha_phys = float(np.clip(alpha_hat, -1.0, 1.0))
    chi2_phys, _ = solve_profile(r0, a_alpha, B, precision, prior_nuis, alpha_phys)

    profile_points = {}
    for aa in (-1.0, 0.0, 1.0):
        cc, _ = solve_profile(r0, a_alpha, B, precision, prior_nuis, aa)
        profile_points[f"{aa:+.1f}"] = {
            "chi2_profile": cc,
            "delta_chi2_from_best_physical": max(cc - chi2_phys, 0.0),
        }

    formal_68 = [alpha_hat - sigma_alpha, alpha_hat + sigma_alpha]
    formal_95 = [alpha_hat - 1.95996398454 * sigma_alpha,
                 alpha_hat + 1.95996398454 * sigma_alpha]

    nuisance = {}
    for i, name in enumerate(names[1:], start=1):
        entry = {"fit_units": float(beta[i]), "sigma_fit_units": float(sigma[i])}
        if name in DERIVATIVE_STEPS:
            entry["physical_shift"] = float(beta[i] * DERIVATIVE_STEPS[name])
            entry["physical_sigma"] = float(sigma[i] * DERIVATIVE_STEPS[name])
        elif name == "A_planck":
            entry["physical_shift"] = float(beta[i] * CAL_PRIOR_SIGMA)
            entry["physical_sigma"] = float(sigma[i] * CAL_PRIOR_SIGMA)
        nuisance[name] = entry

    chi2_fid = float(r0 @ precision @ r0)
    out = {
        "experiment": "Planck 2018 real-data blind hidden-mode fit",
        "likelihood": "public Planck 2018 Plik-lite high-l TTTEEE",
        "planck_lite_py_commit": PLANCK_LITE_COMMIT,
        "source_class_spectra_run_id": 34677322082,
        "theory_units_conversion": {
            "CLASS_input": "dimensionless D_ell",
            "Planck_likelihood_input": "microK^2 D_ell before Plik binning",
            "T_CMB_K": TCMB_K,
            "multiplicative_factor": UK2_PER_DIMENSIONLESS,
        },
        "data_selection": {
            "available_theory_ell_max": ell_max,
            "retained_bins": bin_counts,
            "official_full_bins": {"TT": int(like.nbintt), "TE": int(like.nbinte), "EE": int(like.nbinee)},
            "note": "TT bins requiring ell>2500 are omitted; TE/EE retain the full Plik-lite range to ell=1996.",
        },
        "pre_registered_hidden_template": {
            "mass_eV": 0.60,
            "control": "fixed present-day relic omega_ncdm",
            "alpha_definition": "alpha=+1/-1 are the existing +/-30% pointwise-capped winning hidden-state pair; alpha=0 is the midpoint",
            "physical_validated_range": [-1.0, 1.0],
            "no_data_driven_mass_or_direction_optimization": True,
        },
        "nuisance_model": {
            "linearized_parameters": list(DERIVATIVE_STEPS),
            "derivative_steps": DERIVATIVE_STEPS,
            "tau_gaussian_prior_sigma": TAU_PRIOR_SIGMA,
            "A_planck_gaussian_prior_sigma": CAL_PRIOR_SIGMA,
            "foreground_nuisance_status": "already marginalized in Plik-lite",
            "fit_units": units,
        },
        "formal_unconstrained_linear_template_fit": {
            "alpha_hat": alpha_hat,
            "sigma_alpha": sigma_alpha,
            "alpha_over_sigma": alpha_hat / sigma_alpha if sigma_alpha > 0 else None,
            "delta_chi2_alpha0": delta_chi2_0,
            "formal_68_percent_interval": formal_68,
            "formal_95_percent_interval": formal_95,
            "chi2_min": chi2_min,
            "chi2_alpha0_profiled": chi2_0,
        },
        "physical_family_profile": {
            "alpha_best_in_minus1_plus1": alpha_phys,
            "chi2_best_physical": chi2_phys,
            "profile_at_endpoints_and_zero": profile_points,
            "max_delta_chi2_across_endpoints_vs_best_physical": max(v["delta_chi2_from_best_physical"] for v in profile_points.values()),
            "interpretation_guardrail": "If the formal interval extends outside |alpha|<=1, it only quantifies local template sensitivity; it is not a positive-distribution constraint beyond the validated family.",
        },
        "nuisance_best_fit": nuisance,
        "diagnostics": {
            "chi2_fiducial_no_refit": chi2_fid,
            "number_of_retained_data_bins": int(len(keep)),
            "number_of_fit_parameters_including_alpha": int(p),
            "normal_matrix_condition_number": float(np.linalg.cond(N)),
        },
        "scope": "One pre-defined hidden mode confronted with actual Planck 2018 Plik-lite data. This is a local linearized profile-likelihood test, not a full nonlinear MCMC and not a detection claim.",
    }

    (args.outdir / "planck2018_hidden_mode_fit.json").write_text(json.dumps(out, indent=2) + "\n")
    grid = np.linspace(-1.0, 1.0, 81)
    rows = []
    for aa in grid:
        cc, _ = solve_profile(r0, a_alpha, B, precision, prior_nuis, float(aa))
        rows.append((aa, cc, cc - chi2_phys))
    np.savetxt(
        args.outdir / "planck2018_hidden_mode_profile.csv", np.asarray(rows), delimiter=",",
        header="alpha,chi2_profile,delta_chi2_from_best_physical", comments="",
    )
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
