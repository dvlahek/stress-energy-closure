#!/usr/bin/env python3
"""Build a pre-window physical odd-sector basis for DESI LRG->ELG.

The standard dipole follows Bonvin et al. (MNRAS 525, 4611, 2023), Eq. (8):
a relativistic/Doppler-like nu_1(d,z) contribution plus the leading
wide-angle d/r * mu_2(d,z) contribution.  The hidden-state wake uses the
frozen publication direction already used by build_lrg_elg_wake_template.py.

This builder is deliberately pre-window.  It improves on the old single H/k
toy nuisance by exposing the physically distinct radial basis functions and
by volume-averaging all templates over the same separation bins as the data.
A publication-level prediction still requires DESI-window convolution and
externally calibrated tracer bias / magnification-bias / evolution-bias
parameters.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from classy import Class
from scipy.special import spherical_jn

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_phase7_template as wpt


def trapz_hankel(k, pk, kernel, ell, s):
    k = np.asarray(k, float)
    pk = np.asarray(pk, float)
    kernel = np.asarray(kernel, float)
    s = np.asarray(s, float)
    kk = k[:, None]
    ss = s[None, :]
    integ = (kernel * pk)[:, None] * spherical_jn(ell, kk * ss)
    return np.trapezoid(integ, k, axis=0) / (2.0 * np.pi**2)


def volume_bin_average(sfine, values, edges):
    sfine = np.asarray(sfine, float)
    values = np.asarray(values, float)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mask = (sfine >= lo) & (sfine <= hi)
        if np.count_nonzero(mask) < 2:
            raise RuntimeError(f"insufficient fine-grid support in bin {lo}-{hi}")
        x = sfine[mask]
        y = values[mask]
        # Ideal-shell average.  The exact DESI window upgrade should replace
        # s^2 by the measured random-pair radial weight in each bin.
        out.append(float(np.trapezoid(x*x*y, x) / np.trapezoid(x*x, x)))
    return np.asarray(out, float)


def normalized(x):
    x = np.asarray(x, float)
    scale = float(np.max(np.abs(x)))
    if not np.isfinite(scale) or scale <= 0.0:
        raise RuntimeError("degenerate template")
    return x / scale, scale


def class_geometry(psd, mass, z):
    c = Class()
    c.set(base.class_params(psd, mass))
    c.compute()
    try:
        h = float(c.h())
        r_h = (1.0 + z) * float(c.angular_distance(z)) * h
        Hphys_h = float(c.Hubble(z)) / h
        H0_h = float(c.Hubble(0.0)) / h
        Hconf_h = Hphys_h / (1.0 + z)

        dz = 2.0e-3 * (1.0 + z)
        zm = max(0.0, z - dz)
        zp = z + dz
        Hcm = (float(c.Hubble(zm)) / h) / (1.0 + zm)
        Hcp = (float(c.Hubble(zp)) / h) / (1.0 + zp)
        dlnHc_dz = (np.log(Hcp) - np.log(Hcm)) / (zp - zm)
        Hdot_over_H2 = -(1.0 + z) * dlnHc_dz
    finally:
        c.struct_cleanup()
        c.empty()
    return {
        "h": h,
        "r_Mpc_over_h": float(r_h),
        "Hconf_h_over_Mpc": float(Hconf_h),
        "H0_h_over_Mpc": float(H0_h),
        "Hdot_over_H2": float(Hdot_over_H2),
    }


def physical_coefficients(b1, b2, s1, s2, fe1, fe2, f, geom):
    rH = geom["r_Mpc_over_h"] * geom["Hconf_h_over_Mpc"]
    inv_rH = 1.0 / rH
    Arel = (
        (b1 - b2) * f * (2.0 * inv_rH + geom["Hdot_over_H2"])
        + 3.0 * (s2 - s1) * f*f * (1.0 - inv_rH)
        + 5.0 * (b1*s2 - b2*s1) * f * (1.0 - inv_rH)
        + (3.0/5.0) * (fe1 - fe2) * f*f
        + (b2*fe1 - b1*fe2) * f
    )
    Awa = -(2.0/5.0) * (b1 - b2) * f
    return float(Arel), float(Awa)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--z", type=float, required=True)
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--sep-edges", default="20,40,60,80,100,120,140")
    ap.add_argument("--fine-step", type=float, default=0.25)

    # Optional externally calibrated tracer parameters.  If all are supplied,
    # also write the linked standard GR odd-dipole shape.
    ap.add_argument("--b-lrg", type=float)
    ap.add_argument("--b-elg", type=float)
    ap.add_argument("--s-lrg", type=float)
    ap.add_argument("--s-elg", type=float)
    ap.add_argument("--fevo-lrg", type=float)
    ap.add_argument("--fevo-elg", type=float)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    edges = np.asarray([float(v) for v in args.sep_edges.split(",")], float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    sfine = np.arange(edges[0], edges[-1] + 0.5*args.fine_step, args.fine_step)

    z = float(args.z)
    base.ZBINS = np.asarray([z], float)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, Nnull, shapes, y, M = cro.kinetic_objects(q, args.mass, args.z_match)
    if shapes.shape[0] != len(wpt.COEFF):
        raise RuntimeError("production kinetic basis changed")
    raw = wpt.COEFF @ shapes
    shape = raw / float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
    fp = f0 + args.frac * shape
    fm = f0 - args.frac * shape

    p0 = out / "fd_reference.dat"
    pp = out / "hidden_plus.dat"
    pm = out / "hidden_minus.dat"
    base.write_psd(p0, q, f0)
    base.write_psd(pp, q, fp)
    base.write_psd(pm, q, fm)
    st0 = base.build_state(p0, args.mass)
    stp = base.build_state(pp, args.mass)
    stm = base.build_state(pm, args.mass)

    info = st0["z"][z]
    k = np.asarray(base.KOBS, float)
    pk = np.asarray(info["pk"], float)
    fg = np.asarray(info["f"], float)
    fscalar = float(np.median(fg[k <= 0.05])) if np.any(k <= 0.05) else float(np.median(fg))
    geom = class_geometry(p0, args.mass, z)

    hid = wpt.signed_hidden_response(q, f0, fp, fm, st0, stp, stm, args.mass, z)
    wake_fine = wpt.hankel_dipole(k, pk, hid, sfine)

    # Bonvin Eq. (9): nu1 = integral dk k H0 P j1 / (2 pi^2).
    nu1_fine = trapz_hankel(
        k, pk, k * geom["H0_h_over_Mpc"], 1, sfine
    )
    # Bonvin Eq. (10): mu2 = integral dk k^2 P j2 / (2 pi^2).
    mu2_fine = trapz_hankel(k, pk, k*k, 2, sfine)
    wa_fine = (sfine / geom["r_Mpc_over_h"]) * mu2_fine

    wake = volume_bin_average(sfine, wake_fine, edges)
    nu1 = volume_bin_average(sfine, nu1_fine, edges)
    wa = volume_bin_average(sfine, wa_fine, edges)
    wake_n, wake_scale = normalized(wake)
    nu1_n, nu1_scale = normalized(nu1)
    wa_n, wa_scale = normalized(wa)

    cols = [centers, wake_n, nu1_n, wa_n]
    names = [
        "s_Mpc_over_h",
        "wake_shape",
        "relativistic_nu1_shape",
        "wide_angle_shape",
    ]

    pars = [
        args.b_lrg, args.b_elg, args.s_lrg, args.s_elg,
        args.fevo_lrg, args.fevo_elg,
    ]
    physical = None
    if all(v is not None for v in pars):
        Arel, Awa = physical_coefficients(
            args.b_lrg, args.b_elg, args.s_lrg, args.s_elg,
            args.fevo_lrg, args.fevo_elg, fscalar, geom
        )
        # Bonvin Eq. (8), with tracer-1=LRG and tracer-2=ELG.
        standard = (
            (geom["Hconf_h_over_Mpc"] / geom["H0_h_over_Mpc"]) * Arel * nu1
            + Awa * wa
        )
        standard_n, standard_scale = normalized(standard)
        cols.append(standard_n)
        names.append("standard_odd_shape")
        physical = {
            "b_lrg": args.b_lrg,
            "b_elg": args.b_elg,
            "s_lrg": args.s_lrg,
            "s_elg": args.s_elg,
            "fevo_lrg": args.fevo_lrg,
            "fevo_elg": args.fevo_elg,
            "growth_f": fscalar,
            "A_relativistic": Arel,
            "A_wide_angle": Awa,
            "standard_shape_scale": standard_scale,
        }

    np.savetxt(
        out / "lrg_elg_physical_odd_basis.csv",
        np.column_stack(cols),
        delimiter=",",
        header=",".join(names),
        comments="",
    )

    mp, mm = cro.moments(fp, q, weights), cro.moments(fm, q, weights)
    mis = np.abs(mp-mm) / np.maximum(0.5*(np.abs(mp)+np.abs(mm)), 1e-300)
    summary = {
        "scope": "Pre-window physical odd-sector basis for exact DESI DR1 LRG->ELG analysis",
        "z_effective": z,
        "mass_eV": args.mass,
        "z_match": args.z_match,
        "pointwise_cap": args.frac,
        "separation_edges_Mpc_over_h": edges.tolist(),
        "bin_average": "ideal-shell s^2 weighted; replace by DESI random-pair/window weighting for production",
        "orientation": "tracer1=LRG -> tracer2=ELG; reversing tracer order flips all odd templates",
        "geometry": geom,
        "growth_f_reference": fscalar,
        "growth_f_k_minmax": [float(np.min(fg)), float(np.max(fg))],
        "basis_definitions": {
            "relativistic_nu1_shape": "Bonvin et al. 2023 Eq. (9), normalized after bin averaging",
            "wide_angle_shape": "d/r times Bonvin et al. 2023 Eq. (10), normalized after bin averaging",
            "wake_shape": "frozen hidden-state Hankel dipole, volume-bin-averaged and normalized",
        },
        "normalization_scales": {
            "wake": wake_scale,
            "relativistic_nu1": nu1_scale,
            "wide_angle": wa_scale,
        },
        "physical_standard_shape": physical,
        "max_relative_moment_mismatch": float(np.max(mis)),
        "guardrail": (
            "The two basis columns may be used for geometry diagnostics. Do not select "
            "their amplitudes from the observed odd data. A physical standard_odd_shape "
            "requires externally calibrated b, s and f_evo values and eventual DESI-window convolution."
        ),
    }
    (out / "basis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
