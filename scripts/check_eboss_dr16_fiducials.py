#!/usr/bin/env python3
"""Reference eBOSS redshift-distance mapping and published-convention comparison.

This *input-only* audit compares two published eBOSS fiducial backgrounds
without looking at data or random catalogues. Flat matter-plus-Lambda
distances are evaluated independently with scipy quadrature and astropy.
Radiation/neutrino details are NOT included in this simplified background,
so a future final RR implementation must explicitly freeze that convention.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.integrate import quad

C_KM_S = 299792.458
GRID = (0.6, 0.7, 0.8, 0.9, 1.0)
MODELS = {
    "published_configuration_space_LRG_ELG_multitracer": {
        "h": 0.6777, "Omega_m": 0.307,
        "source": "https://academic.oup.com/mnras/article/498/3/3470/5897383",
        "scope": "Reference mapping for a future eBOSS cross-tracer RR calculation",
    },
    "published_DR16_ELG_BAO": {
        "h": 0.676, "Omega_m": 0.31,
        "source": "https://repositorio.uam.es/bitstream/handle/10486/702297/8787449.pdf",
        "scope": "Predeclared reference-cosmology sensitivity mapping only",
    },
    "legacy_DESI_debug_mapping": {
        "h": 0.674, "Omega_m": 0.315,
        "source": "code/desi_dr1_lrg_elg_exact.py",
        "scope": "Legacy diagnostic mapping; not the eBOSS production fiducial",
    },
}


def distance_mpc_over_h_quad(z: float, h: float, omega_m: float) -> float:
    if z < 0 or not (0 < h < 1.5 and 0 < omega_m < 1):
        raise ValueError("Invalid redshift or flat background")
    integral, err = quad(
        lambda t: 1.0 / math.sqrt(omega_m * (1.0 + t) ** 3
                                  + (1.0 - omega_m)),
        0.0, z, epsrel=2e-12, epsabs=1e-11)
    if not math.isfinite(err) or err > 1e-7:
        raise ValueError("Background distance quadrature did not converge")
    return (C_KM_S / (100.0 * h)) * h * integral


def distance_mpc_over_h_astropy(z: np.ndarray, h: float,
                                omega_m: float) -> np.ndarray:
    model = FlatLambdaCDM(H0=100 * h, Om0=omega_m, Tcmb0=0)
    return np.asarray(model.comoving_distance(z).value, dtype="f8") * h


def audit() -> dict:
    z = np.asarray(GRID, dtype="f8")
    traces = {}
    for name, params in MODELS.items():
        h, om = params["h"], params["Omega_m"]
        quadrature = np.asarray(
            [distance_mpc_over_h_quad(t, h, om) for t in z], dtype="f8")
        astropy = distance_mpc_over_h_astropy(z, h, om)
        rel = float(np.max(np.abs(quadrature - astropy) / astropy))
        if rel > 3e-11:
            raise RuntimeError(f"Independent distance mismatch for {name}: {rel}")
        traces[name] = {
            **params, "z": z.tolist(),
            "comoving_Mpc_over_h": quadrature.tolist(),
            "independent_astropy_max_relative_residual": rel,
        }
        print("FIDUCIAL_DISTANCE_OK", name,
              "z_0p8_Mpc_over_h", quadrature[2],
              "independent_rel", rel, flush=True)
    base = traces["published_configuration_space_LRG_ELG_multitracer"]
    comparison = {}
    for name, record in traces.items():
        if name == "published_configuration_space_LRG_ELG_multitracer":
            continue
        x = np.asarray(record["comoving_Mpc_over_h"])
        y = np.asarray(base["comoving_Mpc_over_h"])
        comparison[name] = {
            "distance_fractional_difference_at_z": ((x - y) / y).tolist(),
            "max_abs_fractional_difference_on_grid": float(
                np.max(np.abs((x - y) / y))),
        }
    # In the simplified radiation-free flat background, h cancels in Mpc/h.
    test = distance_mpc_over_h_quad(0.9, 0.5, 0.307)
    reference = distance_mpc_over_h_quad(0.9, 0.9, 0.307)
    if not math.isclose(test, reference, rel_tol=2e-14):
        raise RuntimeError("Unexpected H0 scaling in Mpc/h")
    return {
        "study": "Published eBOSS fiducial distance preflight",
        "status": "two_published_fiducial_backgrounds_independently_cross_checked",
        "model": "flat matter plus Lambda, negligible radiation, no massive neutrino correction",
        "model_is_exact_published_eboss_distance_implementation": False,
        "z_grid": list(GRID),
        "cosmologies": traces,
        "relative_to_multitracer_reference": comparison,
        "observed_galaxy_or_random_files_read": False,
        "observed_odd_vector_read": False,
        "RR_pair_window_computed": False,
        "fiducial_distance_convention_frozen_for_inference": False,
        "note": (
            "Both literature backgrounds are independently verified in a "
            "simplified flat-LCDM distance implementation. A final production "
            "distance implementation must explicitly settle any required "
            "radiation and neutrino contributions before the RR calculation. "
            "No selection, bins or physics template is tuned here."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="eboss_workspace/fiducial_distance_audit.json")
    args = parser.parse_args()
    result = audit()
    target = Path(args.out)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("FIDUCIAL_DISTANCE_AUDIT", result["status"], target, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
