#!/usr/bin/env python3
"""Calibration-hierarchy test for the parity-odd hidden-state wake forecast.

Re-optimizes the matched n-rho-P null space under four nuisance models:
  global:       one common wake normalization across all redshift bins
  linear:       global normalization plus linear redshift calibration drift
  quadratic:    global + linear + quadratic redshift calibration drift
  perbin:       independent wake normalization in every redshift bin

All models additionally project conservative k^-1, k, and k^2 odd broadband
shapes independently in each redshift bin. Candidate directions from all four
linear Fisher optimizations are pooled, then each unique candidate is validated
once with full CLASS F+/F- transfer functions and scored under every hierarchy.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path
import numpy as np

import class_response_optimize as cro
import wake_two_tracer_fisher as base

LEVELS = ("global", "linear", "quadratic", "perbin")
ZMEAN = float(np.mean(base.ZBINS))
ZSCALE = float(np.max(np.abs(base.ZBINS - ZMEAN)))


def amp_terms(z: float, level: str):
    x = (z - ZMEAN) / max(ZSCALE, 1e-30)
    if level == "global":
        return [1.0]
    if level == "linear":
        return [1.0, x]
    if level == "quadratic":
        return [1.0, x, x*x]
    if level == "perbin":
        return None
    raise ValueError(level)


def nuisance_size(level: str):
    nz = len(base.ZBINS)
    namp = nz if level == "perbin" else len(amp_terms(float(base.ZBINS[0]), level))
    return namp + 3*nz


def nuisance_shapes(iz: int, z: float, k: float, level: str):
    nz = len(base.ZBINS)
    vals = []
    if level == "perbin":
        vals.extend([1.0 if j == iz else 0.0 for j in range(nz)])
    else:
        vals.extend(amp_terms(z, level))
    # Conservative independent odd-shape nuisances per redshift bin.  No
    # constant term here, because the wake amplitude is handled above.
    for j in range(nz):
        a = 1.0 if j == iz else 0.0
        vals.extend([a*(0.01/k), a*(k/0.01), a*(k/0.01)**2])
    return np.asarray(vals, float)


def calibrated_linear_fisher(qgrid, f0, shapes, state0, mass, frac, volume_total, level):
    nd = shapes.shape[0]
    nn = nuisance_size(level)
    F = np.zeros((nd, nd)); G = np.zeros((nd, nn)); N = np.zeros((nn, nn))
    volume_bin = volume_total / len(base.ZBINS)
    perz = []
    for iz, zv in enumerate(base.ZBINS):
        z = float(zv); info = state0["z"][z]
        rows = base.mode_noise_weight(info, volume_bin)
        cache = {}; raw_base = 0.0
        for ik, im, w in rows:
            if ik not in cache:
                cache[ik] = base.linear_basis_moments(qgrid, f0, shapes, state0, mass, z, ik, frac)
            m00, m0, mm = cache[ik]
            raw_base += w*m00
        target = base.desired_fd_sn2(mass, z, volume_bin)
        cal = target / max(raw_base, 1e-300)
        perz.append({"z": z, "FD_SN": math.sqrt(target), "calibration": cal})
        for ik, im, w0 in rows:
            w = w0*cal; m00, m0, mm = cache[ik]
            ns = nuisance_shapes(iz, z, base.KOBS[ik], level)
            F += w*mm
            G += w*np.outer(m0, ns)
            N += w*m00*np.outer(ns, ns)
    P = F - G @ np.linalg.pinv(N, rcond=1e-11) @ G.T
    return 0.5*(P+P.T), 0.5*(F+F.T), perz


def nonlinear_fisher(qgrid, f0, fp, fm, state0, statep, statem, mass, volume_total, level):
    nn = nuisance_size(level)
    Fdd = 0.0; g = np.zeros(nn); N = np.zeros((nn, nn))
    volume_bin = volume_total / len(base.ZBINS)
    perz = []
    for iz, zv in enumerate(base.ZBINS):
        z = float(zv); info = state0["z"][z]
        rows = base.mode_noise_weight(info, volume_bin)
        cache = {}; raw_base = 0.0
        for ik, im, w in rows:
            if ik not in cache:
                cache[ik] = base.stochastic_moments(qgrid, f0, fp, fm, state0, statep, statem, mass, z, ik)
            m00, mdd, m0d = cache[ik]
            raw_base += w*m00
        target = base.desired_fd_sn2(mass, z, volume_bin)
        cal = target/max(raw_base, 1e-300)
        zun = 0.0
        for ik, im, w0 in rows:
            w = w0*cal; m00, mdd, m0d = cache[ik]
            ns = nuisance_shapes(iz, z, base.KOBS[ik], level)
            Fdd += w*mdd; zun += w*mdd
            g += w*m0d*ns
            N += w*m00*np.outer(ns, ns)
        perz.append({"z": z, "unprojected_hidden_SN": math.sqrt(max(zun, 0.0))})
    proj = max(Fdd - float(g @ np.linalg.pinv(N, rcond=1e-11) @ g), 0.0)
    return math.sqrt(max(Fdd, 0.0)), math.sqrt(proj), perz


def candidate_pool(Bproj, shapes, f0, samples, seed, keep):
    nd = shapes.shape[0]
    evals, evecs = np.linalg.eigh(0.5*(Bproj+Bproj.T))
    cand = [evecs[:, j] for j in np.argsort(evals)[::-1]]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0; cand.extend([e, -e])
    rng = np.random.default_rng(seed)
    rr = rng.normal(size=(samples, nd)); rr /= np.linalg.norm(rr, axis=1, keepdims=True)
    cand.extend(rr)
    scored = []
    for c in cand:
        raw = c @ shapes
        mr = float(np.max(np.abs(raw)/np.maximum(f0, 1e-300)))
        if not np.isfinite(mr) or mr <= 0: continue
        norm = 1.0/mr
        s2 = float(norm*norm*(c @ Bproj @ c))
        scored.append((s2, c.copy(), norm, raw*norm))
    scored.sort(key=lambda x: x[0], reverse=True)
    out = []
    for item in scored:
        u = item[1]/np.linalg.norm(item[1])
        if all(abs(float(u @ (v[1]/np.linalg.norm(v[1])))) < 0.9995 for v in out):
            out.append(item)
        if len(out) >= keep: break
    return out


def dedupe(items, max_keep):
    out = []
    for source_level, item in items:
        u = item[1]/np.linalg.norm(item[1])
        if all(abs(float(u @ (x[1][1]/np.linalg.norm(x[1][1])))) < 0.9995 for x in out):
            out.append((source_level, item))
        if len(out) >= max_keep: break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="wake_calibration_hierarchy_output")
    ap.add_argument("--mass", type=float, default=0.06)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--frac", type=float, default=0.30)
    ap.add_argument("--volume", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=40000)
    ap.add_argument("--per-level", type=int, default=4)
    ap.add_argument("--validate", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, Nnull, shapes, y, M = cro.kinetic_objects(q, args.mass, args.z_match)
    p0 = out/"fd_reference.dat"; base.write_psd(p0, q, f0)
    state0 = base.build_state(p0, args.mass)

    linear = {}; pool = []
    for i, level in enumerate(LEVELS):
        Bp, Br, perz = calibrated_linear_fisher(q, f0, shapes, state0, args.mass, args.frac, args.volume, level)
        sel = candidate_pool(Bp, shapes, f0, args.samples, args.seed+i, args.per_level)
        linear[level] = {
            "best_predicted_projected_SN": math.sqrt(max(sel[0][0], 0.0)),
            "nuisance_dimension": nuisance_size(level),
        }
        pool.extend((level, x) for x in sel)

    pool.sort(key=lambda x: x[1][0], reverse=True)
    selected = dedupe(pool, args.validate)
    validations = []
    best = {lev: None for lev in LEVELS}

    for rank, (source_level, (pred2, coeff, norm, shape)) in enumerate(selected):
        fp = f0 + args.frac*shape; fm = f0 - args.frac*shape
        pp = out/f"cand_{rank:02d}_plus.dat"; pm = out/f"cand_{rank:02d}_minus.dat"
        base.write_psd(pp, q, fp); base.write_psd(pm, q, fm)
        statep = base.build_state(pp, args.mass); statem = base.build_state(pm, args.mass)
        mp, mm = cro.moments(fp, q, weights), cro.moments(fm, q, weights)
        mismatch = np.abs(mp-mm)/np.maximum(0.5*(np.abs(mp)+np.abs(mm)), 1e-300)
        scores = {}
        sn_un = None
        for level in LEVELS:
            un, proj, perz = nonlinear_fisher(q, f0, fp, fm, state0, statep, statem, args.mass, args.volume, level)
            sn_un = un
            scores[level] = proj
            if best[level] is None or proj > best[level][0]:
                best[level] = (proj, rank, fp.copy(), fm.copy(), shape.copy())
        row = {
            "rank": rank,
            "source_linear_level": source_level,
            "source_predicted_SN": math.sqrt(max(pred2,0.0)),
            "validated_unprojected_hidden_SN": sn_un,
            "validated_projected_SN": scores,
            "max_relative_moment_mismatch": float(mismatch.max()),
            "coefficients": coeff.tolist(),
        }
        validations.append(row)
        print(json.dumps(row))

    best_summary = {}
    for level in LEVELS:
        proj, rank, fp, fm, shape = best[level]
        best_summary[level] = {
            "validated_projected_SN": float(proj),
            "candidate_rank": int(rank),
            "volume_1sigma_hminus3_Gpc3": float(args.volume/max(proj*proj, 1e-30)),
            "volume_2sigma_hminus3_Gpc3": float(4.0*args.volume/max(proj*proj, 1e-30)),
            "volume_3sigma_hminus3_Gpc3": float(9.0*args.volume/max(proj*proj, 1e-30)),
        }
        np.savetxt(out/f"best_{level}_pair.csv", np.column_stack([q,f0,fp,fm,shape]), delimiter=",",
                   header="q,f_FD,f_plus,f_minus,normalized_delta_shape", comments="")

    summary = {
        "class_commit": cro.CLASS_COMMIT,
        "mass_eV": args.mass,
        "fractional_distortion_cap": args.frac,
        "survey_volume_hminus3_Gpc3": args.volume,
        "redshift_bins": base.ZBINS.tolist(),
        "calibration_models": {
            "global": "one common wake normalization across all redshift bins",
            "linear": "common normalization plus linear drift in centered/scaled redshift",
            "quadratic": "common normalization plus linear and quadratic redshift drift",
            "perbin": "independent wake normalization in each redshift bin",
        },
        "common_odd_nuisances": "independent k^-1, k, and k^2 odd broadband templates in each redshift bin",
        "linear_optimization": linear,
        "nonlinear_candidates_validated": len(validations),
        "best_by_model": best_summary,
        "all_validations": validations,
        "interpretation": "Literature-calibrated Okoli-type two-tracer wake Fisher. Each calibration hierarchy is optimized in the exact matched-moment null space, then pooled candidates are validated with full CLASS F+/F- transfer functions.",
    }
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print("SUMMARY\n"+json.dumps(summary["best_by_model"], indent=2))

if __name__ == "__main__":
    main()
