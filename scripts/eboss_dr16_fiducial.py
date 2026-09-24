#!/usr/bin/env python3
"""Declared eBOSS LRG-ELG pre-odd distance and column-weight conventions.

The geometry compares the published multitracer fiducial (h=0.6777,
Om0=0.307) against the DR16 ELG fiducial (h=0.676, Om0=0.310).
Only the homogeneous flat matter+Lambda comoving-distance mapping is used,
with Tcmb0=0 so radiation/neutrino corrections are explicitly omitted.
The two literature conventions are not a mock-versus-data calibration:
the *same* chosen mapping must apply to both observed and mock inputs.
Neither a distance nor a weight option may be chosen from the odd signal.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from astropy.cosmology import FlatLambdaCDM

WEIGHT_COLUMNS = (
    "WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"
)
SYSTOT_NUMERICAL_ZERO_TOL = 1e-12
DISTANCE_CONVENTIONS = {
    "published_multitracer": {
        "h": 0.6777, "om0": 0.307,
        "citation_url": "https://academic.oup.com/mnras/article/498/3/3470/5897383",
        "context": "Wang et al. DR16 LRG/ELG multitracer fiducial",
    },
    "dr16_elg_bao": {
        "h": 0.676, "om0": 0.310,
        "citation_url": "https://academic.oup.com/mnras/article/500/3/3254/5942664",
        "context": "DR16 ELG BAO fiducial, explicitly secondary here",
    },
}
PRIMARY_GEOMETRY = "published_multitracer"
SECONDARY_GEOMETRY = "dr16_elg_bao"


def comoving_mpc_over_h(redshift: Any, convention: str = PRIMARY_GEOMETRY
                       ) -> np.ndarray:
    if convention not in DISTANCE_CONVENTIONS:
        raise ValueError("Unknown, undeclared eBOSS fiducial convention")
    z = np.asarray(redshift, dtype="f8")
    if not np.all(np.isfinite(z)) or np.any(z < 0):
        raise ValueError("Redshifts must be finite and nonnegative")
    cfg = DISTANCE_CONVENTIONS[convention]
    geometry = FlatLambdaCDM(
        H0=100.0 * cfg["h"], Om0=cfg["om0"], Tcmb0=0.0
    )
    return np.asarray(
        geometry.comoving_distance(z).value * cfg["h"], dtype="f8")


def validated_weight_product(columns: dict[str, Any]
                             ) -> tuple[np.ndarray, np.ndarray]:
    """Return selected total weights and the identical data/random row gate.

    Reject nonfinite, significant negative and nonpositive non-systematic
    weights. Exact or negligible WEIGHT_SYSTOT rows are zero-weight excluded.
    No CMASS additive-correction formula is used for pure eBOSS LRG/ELG.
    """
    missing = set(WEIGHT_COLUMNS) - set(columns)
    if missing:
        raise ValueError(f"Missing required catalogue weights: {sorted(missing)}")
    arrays = [np.asarray(columns[name], dtype="f8") for name in WEIGHT_COLUMNS]
    if not arrays or any(a.shape != arrays[0].shape for a in arrays):
        raise ValueError("Weight columns do not have identical shapes")
    if any(a.ndim != 1 for a in arrays):
        raise ValueError("Weight columns must be one-dimensional")
    wsys, *other = arrays
    if not all(np.isfinite(a).all() for a in arrays):
        raise ValueError("A required catalogue weight is nonfinite")
    if np.any(wsys < -SYSTOT_NUMERICAL_ZERO_TOL):
        raise ValueError("Significantly negative WEIGHT_SYSTOT")
    if any(np.any(w <= 0) for w in other):
        raise ValueError("A non-systematic catalogue weight is nonpositive")
    retained = wsys > SYSTOT_NUMERICAL_ZERO_TOL
    total = np.ones(wsys.shape, dtype="f8")
    for array in arrays:
        total *= np.where(retained, array, 1.0)
    total[~retained] = 0.0
    if not np.isfinite(total).all() or np.any(total[retained] <= 0):
        raise ValueError("A combined catalogue weight is invalid")
    return total, retained


def self_test() -> None:
    from scipy.integrate import quad

    a = np.array([0.0, 0.6, 0.8, 1.0])
    for name, cfg in DISTANCE_CONVENTIONS.items():
        distances = comoving_mpc_over_h(a, name)
        assert distances[0] == 0 and np.all(np.diff(distances) > 0)
        # Independent dimensionless flat-matter+Lambda distance integral.
        for z, got in zip(a[1:], distances[1:]):
            expected = 2997.92458 * quad(
                lambda t: 1 / np.sqrt(
                    cfg["om0"] * (1 + t) ** 3 + (1 - cfg["om0"])),
                0, z, epsabs=1e-12, epsrel=1e-12,
            )[0]
            assert abs(got / expected - 1) < 1e-10
    x = {
        "WEIGHT_SYSTOT": [1.0, 1e-30, 2.0],
        "WEIGHT_CP": [1.2, 1.2, 1.0],
        "WEIGHT_NOZ": [1.1, 1.1, 1.0],
        "WEIGHT_FKP": [0.6, 0.6, 0.5],
    }
    product, gate = validated_weight_product(x)
    np.testing.assert_allclose(product, [1.0 * 1.2 * 1.1 * 0.6, 0.0, 1.0])
    assert gate.tolist() == [True, False, True]
    x["WEIGHT_SYSTOT"][1] = -1e-9
    try:
        validated_weight_product(x)
    except ValueError:
        pass
    else:
        raise AssertionError("Significant negative weights must be rejected")
    print("EBOSS_FIDUCIAL_AND_WEIGHTS_SELF_TEST_OK")


if __name__ == "__main__":
    self_test()
