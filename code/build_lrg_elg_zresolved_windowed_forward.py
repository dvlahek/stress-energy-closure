#!/usr/bin/env python3
"""Forward-model frozen DESI LRGxELG templates through the measured RR window.

The model is data-blind: it uses only fine random-pair counts plus the frozen
theory definition.  It follows the correlation-function RR-window construction
used in the DESI/desilike stack, generalized here to odd output multipoles.

For every redshift slice we build response blocks M_{ell_out,ell_in} from
fine R_LRG R_ELG(s,mu) counts.  We then propagate:
  * the frozen hidden-state wake dipole (ell=1);
  * the linked relativistic + leading wide-angle standard dipole (ell=1);
  * linear Kaiser LRGxELG even multipoles ell=0,2,4.

The even multipoles are included because a cross-tracer survey window and
finite separation binning can leak even clustering into an observed odd
multipole.  No observed odd-sector values are read by this script.
"""
from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path

import numpy as np
from scipy.special import legendre

import class_response_optimize as cro
import wake_two_tracer_fisher as base
import wake_phase7_template as wpt
import build_lrg_elg_physical_odd_basis as phys


OUT_ELLS = (1, 3)
IN_ELLS = (0, 1, 2, 3, 4)


def normalized(x):
    x = np.asarray(x, float)
    scale = float(np.max(np.abs(x)))
    if not np.isfinite(scale) or scale <= 0:
        raise RuntimeError("Cannot normalize degenerate vector")
    return x / scale, scale


def cosine(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / den) if den > 0 else None


def volume_effective_centers(edges):
    lo = np.asarray(edges[:-1], float)
    hi = np.asarray(edges[1:], float)
    return 0.75 * (hi**4 - lo**4) / (hi**3 - lo**3)


def volume_average_from_fine(values, fine_edges, output_edges):
    values = np.asarray(values, float)
    vol = fine_edges[1:]**3 - fine_edges[:-1]**3
    out = []
    for lo, hi in zip(output_edges[:-1], output_edges[1:]):
        m = (fine_edges[:-1] >= lo - 1e-12) & (fine_edges[1:] <= hi + 1e-12)
        if not np.any(m):
            raise RuntimeError(f"No fine bins in output bin {lo}-{hi}")
        out.append(float(np.sum(values[m] * vol[m]) / np.sum(vol[m])))
    return np.asarray(out)


def window_block(fine_edges, muedges, wcounts, output_edges, ellout, ellin):
    """Return matrix mapping fine xi_ellin(s) to coarse observed xi_ellout."""
    ns = len(fine_edges) - 1
    nout = len(output_edges) - 1
    if wcounts.shape != (ns, len(muedges) - 1):
        raise ValueError("RR array shape does not match fine s/mu edges")

    polyint = (legendre(int(ellout)) * legendre(int(ellin))).integ()
    lint = polyint(muedges[1:]) - polyint(muedges[:-1])
    dmu = np.diff(muedges)
    M = np.zeros((nout, ns), dtype="f8")

    for ib, (lo, hi) in enumerate(zip(output_edges[:-1], output_edges[1:])):
        idx = np.flatnonzero(
            (fine_edges[:-1] >= lo - 1e-12) & (fine_edges[1:] <= hi + 1e-12)
        )
        if not idx.size:
            raise RuntimeError(f"No RR fine bins inside output bin {lo}-{hi}")
        wc = np.asarray(wcounts[idx], float)
        wcmu = np.sum(wc, axis=0)
        valid = wcmu > 0
        if not np.any(valid):
            raise RuntimeError(f"No nonzero RR mu cells in output bin {lo}-{hi}")

        cond = wc[:, valid] / wcmu[valid][None, :]
        mu_norm = float(np.sum(dmu[valid]))
        mu_kernel = (2 * int(ellout) + 1) * lint[valid] / mu_norm
        M[ib, idx] = cond @ mu_kernel
    return M


def lookup_linked_row(summary, lo, hi):
    for row in summary["per_bin"]:
        if np.isclose(float(row["zlo"]), lo, atol=1e-12, rtol=0) and np.isclose(
            float(row["zhi"]), hi, atol=1e-12, rtol=0
        ):
            return row
    raise KeyError(f"No linked-standard metadata for z={lo}-{hi}")


def build_kinetic_states(zgrid, mass, z_match, frac, workdir):
    base.ZBINS = np.asarray(zgrid, float)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, Nnull, shapes, y, M = cro.kinetic_objects(q, mass, z_match)
    if shapes.shape[0] != len(wpt.COEFF):
        raise RuntimeError("production kinetic basis changed")
    raw = wpt.COEFF @ shapes
    shape = raw / float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
    fp = f0 + frac * shape
    fm = f0 - frac * shape

    workdir = Path(workdir)
    p0 = workdir / "fd_reference.dat"
    pp = workdir / "hidden_plus.dat"
    pm = workdir / "hidden_minus.dat"
    base.write_psd(p0, q, f0)
    base.write_psd(pp, q, fp)
    base.write_psd(pm, q, fm)
    st0 = base.build_state(p0, mass)
    stp = base.build_state(pp, mass)
    stm = base.build_state(pm, mass)
    return q, f0, fp, fm, p0, st0, stp, stm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rr-counts", required=True)
    ap.add_argument("--prewindow-template", required=True)
    ap.add_argument("--linked-summary", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    args = ap.parse_args()

    R = np.load(args.rr_counts, allow_pickle=False)
    T = np.genfromtxt(args.prewindow_template, delimiter=",", names=True)
    linked_summary = json.loads(Path(args.linked_summary).read_text())

    zlo_rr = np.asarray(R["zlo"], float)
    zhi_rr = np.asarray(R["zhi"], float)
    fine_edges = np.asarray(R["sedges"], float)
    muedges = np.asarray(R["muedges"], float)
    W = np.asarray(R["wcounts"], float)
    outedges = np.asarray(R["output_sedges"], float)

    # Frozen coarse template grid.
    zlo_t = np.asarray(T["zlo"], float)
    zhi_t = np.asarray(T["zhi"], float)
    ze_t = np.asarray(T["z_effective"], float)
    s_t = np.asarray(T["s_Mpc_over_h"], float)
    unique_bins = list(dict.fromkeys((float(a), float(b)) for a, b in zip(zlo_t, zhi_t)))
    if len(unique_bins) != len(zlo_rr):
        raise RuntimeError("RR and template redshift-bin counts differ")

    zgrid = []
    for iz, (lo, hi) in enumerate(unique_bins):
        if not (
            np.isclose(zlo_rr[iz], lo, atol=1e-12, rtol=0)
            and np.isclose(zhi_rr[iz], hi, atol=1e-12, rtol=0)
        ):
            raise RuntimeError("RR and template z-bin ordering differs")
        m = np.isclose(zlo_t, lo, atol=1e-12) & np.isclose(zhi_t, hi, atol=1e-12)
        zgrid.append(float(np.mean(ze_t[m])))

    sfine = volume_effective_centers(fine_edges)
    nout = len(outedges) - 1
    if not all(np.count_nonzero(np.isclose(zlo_t, lo) & np.isclose(zhi_t, hi)) == nout for lo, hi in unique_bins):
        raise RuntimeError("Template does not have one row per frozen output separation bin")

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    matrices = np.zeros(
        (len(unique_bins), len(OUT_ELLS), len(IN_ELLS), nout, len(sfine)), dtype="f8"
    )
    for iz in range(len(unique_bins)):
        for io, ellout in enumerate(OUT_ELLS):
            for ii, ellin in enumerate(IN_ELLS):
                matrices[iz, io, ii] = window_block(
                    fine_edges, muedges, W[iz], outedges, ellout, ellin
                )

    # Constant-in-s closure: with full symmetric mu support, each block should
    # integrate to delta_{ellout,ellin}. This is a window-only diagnostic.
    constant_response = {}
    for io, ellout in enumerate(OUT_ELLS):
        for ii, ellin in enumerate(IN_ELLS):
            vals = matrices[:, io, ii].sum(axis=-1)
            constant_response[f"{ellout}<-{ellin}"] = vals.tolist()

    wake_ideal, std_ideal = [], []
    wake_win1, stdodd_win1, even_win1, totalstd_win1 = [], [], [], []
    wake_win3, stdodd_win3, even_win3, totalstd_win3 = [], [], [], []
    raw_rows = []
    per_bin_meta = []

    with tempfile.TemporaryDirectory(prefix="lrg_elg_window_theory_") as td:
        q, f0, fp, fm, p0, st0, stp, stm = build_kinetic_states(
            zgrid, args.mass, args.z_match, args.frac, td
        )

        for iz, ((lo, hi), z) in enumerate(zip(unique_bins, zgrid)):
            info = st0["z"][z]
            k = np.asarray(base.KOBS, float)
            pk = np.asarray(info["pk"], float)
            linked = lookup_linked_row(linked_summary, lo, hi)
            f = float(linked["growth_f"])
            geom = phys.class_geometry(p0, args.mass, z)

            hid = wpt.signed_hidden_response(q, f0, fp, fm, st0, stp, stm, args.mass, z)
            xi1_wake = wpt.hankel_dipole(k, pk, hid, sfine)

            nu1 = phys.trapz_hankel(k, pk, k * geom["H0_h_over_Mpc"], 1, sfine)
            mu2_odd = phys.trapz_hankel(k, pk, k * k, 2, sfine)
            wa = (sfine / geom["r_Mpc_over_h"]) * mu2_odd
            xi1_std = (
                (geom["Hconf_h_over_Mpc"] / geom["H0_h_over_Mpc"])
                * float(linked["A_relativistic"])
                * nu1
                + float(linked["A_wide_angle"]) * wa
            )

            b1 = float(linked["b_lrg"])
            b2 = float(linked["b_elg"])
            c0 = b1 * b2 + (b1 + b2) * f / 3.0 + f * f / 5.0
            c2 = 2.0 * (b1 + b2) * f / 3.0 + 4.0 * f * f / 7.0
            c4 = 8.0 * f * f / 35.0
            h0 = phys.trapz_hankel(k, pk, k * k, 0, sfine)
            h2 = phys.trapz_hankel(k, pk, k * k, 2, sfine)
            h4 = phys.trapz_hankel(k, pk, k * k, 4, sfine)
            xi0 = c0 * h0
            xi2 = -c2 * h2  # i^2 = -1 in the Fourier-to-configuration transform
            xi4 = c4 * h4   # i^4 = +1

            Min = {
                (ellout, ellin): matrices[
                    iz, OUT_ELLS.index(ellout), IN_ELLS.index(ellin)
                ]
                for ellout in OUT_ELLS
                for ellin in IN_ELLS
            }
            w1 = Min[(1, 1)] @ xi1_wake
            s1 = Min[(1, 1)] @ xi1_std
            e1 = Min[(1, 0)] @ xi0 + Min[(1, 2)] @ xi2 + Min[(1, 4)] @ xi4
            w3 = Min[(3, 1)] @ xi1_wake
            s3 = Min[(3, 1)] @ xi1_std
            e3 = Min[(3, 0)] @ xi0 + Min[(3, 2)] @ xi2 + Min[(3, 4)] @ xi4

            wi = volume_average_from_fine(xi1_wake, fine_edges, outedges)
            si = volume_average_from_fine(xi1_std, fine_edges, outedges)

            wake_ideal.extend(wi.tolist())
            std_ideal.extend(si.tolist())
            wake_win1.extend(w1.tolist())
            stdodd_win1.extend(s1.tolist())
            even_win1.extend(e1.tolist())
            totalstd_win1.extend((s1 + e1).tolist())
            wake_win3.extend(w3.tolist())
            stdodd_win3.extend(s3.tolist())
            even_win3.extend(e3.tolist())
            totalstd_win3.extend((s3 + e3).tolist())

            for ib, s in enumerate(0.5 * (outedges[:-1] + outedges[1:])):
                raw_rows.append(
                    (
                        lo, hi, z, s,
                        wi[ib], si[ib],
                        w1[ib], s1[ib], e1[ib], s1[ib] + e1[ib],
                        w3[ib], s3[ib], e3[ib], s3[ib] + e3[ib],
                    )
                )

            per_bin_meta.append(
                {
                    "zlo": lo,
                    "zhi": hi,
                    "z_effective": z,
                    "b_lrg": b1,
                    "b_elg": b2,
                    "growth_f": f,
                    "kaiser_coefficients": {"ell0": c0, "ell2": c2, "ell4": c4},
                    "geometry": geom,
                }
            )

    wake_ideal = np.asarray(wake_ideal)
    std_ideal = np.asarray(std_ideal)
    wake_win1 = np.asarray(wake_win1)
    stdodd_win1 = np.asarray(stdodd_win1)
    even_win1 = np.asarray(even_win1)
    totalstd_win1 = np.asarray(totalstd_win1)

    wake_ideal_n, _ = normalized(wake_ideal)
    std_ideal_n, _ = normalized(std_ideal)
    wake_win1_n, wake_win_scale = normalized(wake_win1)
    stdodd_win1_n, stdodd_win_scale = normalized(stdodd_win1)
    totalstd_win1_n, totalstd_win_scale = normalized(totalstd_win1)

    old_wake = np.asarray(T["wake_shape"], float)
    old_std = np.asarray(T["standard_linked_shape"], float)
    if len(old_wake) != len(wake_win1_n):
        raise RuntimeError("Pre-window template length mismatch")

    # Output row ordering is the same frozen z-major, s-minor ordering.
    centers = 0.5 * (outedges[:-1] + outedges[1:])
    rows = []
    i = 0
    for (lo, hi), z in zip(unique_bins, zgrid):
        for s in centers:
            rows.append(
                (
                    lo, hi, z, s,
                    wake_win1_n[i],
                    stdodd_win1_n[i],
                    totalstd_win1_n[i],
                    even_win1[i],
                    wake_win1[i],
                    stdodd_win1[i],
                    totalstd_win1[i],
                    np.asarray(wake_win3)[i],
                    np.asarray(stdodd_win3)[i],
                    np.asarray(even_win3)[i],
                    np.asarray(totalstd_win3)[i],
                )
            )
            i += 1

    np.savetxt(
        out / "lrg_elg_zresolved_windowed_templates.csv",
        np.asarray(rows, float),
        delimiter=",",
        header=(
            "zlo,zhi,z_effective,s_Mpc_over_h,"
            "wake_windowed_shape,standard_odd_windowed_shape,standard_total_windowed_shape,"
            "even_to_odd_window_leakage_raw,wake_windowed_raw,standard_odd_windowed_raw,"
            "standard_total_windowed_raw,wake_to_xi3_raw,standard_odd_to_xi3_raw,"
            "even_to_xi3_window_leakage_raw,standard_total_xi3_raw"
        ),
        comments="",
    )
    np.savetxt(
        out / "lrg_elg_zresolved_windowed_raw_components.csv",
        np.asarray(raw_rows, float),
        delimiter=",",
        header=(
            "zlo,zhi,z_effective,s_Mpc_over_h,wake_ideal_raw,standard_odd_ideal_raw,"
            "wake_xi1_windowed_raw,standard_odd_xi1_windowed_raw,even_xi1_leakage_raw,"
            "standard_total_xi1_windowed_raw,wake_xi3_windowed_raw,"
            "standard_odd_xi3_windowed_raw,even_xi3_leakage_raw,standard_total_xi3_windowed_raw"
        ),
        comments="",
    )

    np.savez_compressed(
        out / "lrg_elg_zresolved_window_matrix.npz",
        output_ells=np.asarray(OUT_ELLS, dtype="i8"),
        input_ells=np.asarray(IN_ELLS, dtype="i8"),
        fine_s=sfine,
        fine_sedges=fine_edges,
        muedges=muedges,
        output_sedges=outedges,
        matrix=matrices,
    )

    rms = lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))
    summary = {
        "scope": "Data-blind RR-window-convolved forward model for frozen DESI DR1 LRGxELG templates",
        "uses_observed_odd_data_vector": False,
        "window_method": (
            "fine R_LRG R_ELG(s,mu) conditional radial weighting within each frozen output s-bin; "
            "Legendre-product projection generalized to odd output multipoles"
        ),
        "output_multipoles": list(OUT_ELLS),
        "input_multipoles": list(IN_ELLS),
        "theory": {
            "wake": "frozen hidden-state ell=1 template",
            "standard_odd": "linked Bonvin relativistic + leading wide-angle ell=1 benchmark",
            "standard_even": "linear Kaiser LRGxELG ell=0,2,4 using the same fixed b_LRG, b_ELG and per-bin growth f",
            "mass_eV": float(args.mass),
            "z_match": float(args.z_match),
            "pointwise_cap": float(args.frac),
        },
        "per_zbin": per_bin_meta,
        "constant_radial_response_by_zbin": constant_response,
        "prewindow_recomputation_closure": {
            "wake_cosine_vs_frozen": cosine(wake_ideal_n, old_wake),
            "wake_max_abs_normalized_difference": float(np.max(np.abs(wake_ideal_n - old_wake))),
            "standard_cosine_vs_frozen": cosine(std_ideal_n, old_std),
            "standard_max_abs_normalized_difference": float(np.max(np.abs(std_ideal_n - old_std))),
        },
        "window_effect_on_frozen_shapes": {
            "wake_prewindow_vs_windowed_cosine": cosine(old_wake, wake_win1_n),
            "wake_max_abs_normalized_change": float(np.max(np.abs(old_wake - wake_win1_n))),
            "standard_prewindow_vs_windowed_oddonly_cosine": cosine(old_std, stdodd_win1_n),
            "standard_prewindow_vs_windowed_total_cosine": cosine(old_std, totalstd_win1_n),
            "standard_oddonly_vs_total_windowed_cosine": cosine(stdodd_win1_n, totalstd_win1_n),
        },
        "even_to_odd_leakage": {
            "rms_even_leakage_raw": rms(even_win1),
            "rms_windowed_physical_odd_standard_raw": rms(stdodd_win1),
            "rms_ratio_even_leakage_to_physical_odd_standard": (
                rms(even_win1) / rms(stdodd_win1) if rms(stdodd_win1) > 0 else None
            ),
            "max_abs_even_leakage_raw": float(np.max(np.abs(even_win1))),
            "max_abs_physical_odd_standard_raw": float(np.max(np.abs(stdodd_win1))),
        },
        "normalization_scales": {
            "wake_windowed": wake_win_scale,
            "standard_odd_windowed": stdodd_win_scale,
            "standard_total_windowed": totalstd_win_scale,
        },
        "analysis_scope": (
            "We derive the response from random-pair counts and prescribed theoretical inputs, independently "
            "of the measured odd-sector amplitudes. The even sector uses linear Kaiser multipoles; "
            "nonlinear corrections are outside this benchmark."
        ),
    }
    (out / "windowed_forward_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
