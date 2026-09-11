#!/usr/bin/env python3
"""Optimize stress-energy-matched relic distributions for scalar CMB/lensing response.

Experimental diagnostic only. Uses the same smooth 10-function deformation basis and
exact (n,rho,P) null space as class_response_optimize.py, but targets scalar TT/TE/EE
or phi-phi response instead of tensor BB.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
import numpy as np

from class_response_optimize import (
    H0, OMEGA_B, OMEGA_CDM, A_S, N_UR, T_NCDM, CLASS_COMMIT,
    kinetic_objects, pair_stats, write_psd,
)

N_S = 0.9649
TAU_REIO = 0.054


def write_scalar_ini(path: Path, psd: Path, root: Path, mass_eV: float,
                     lmax: int, pkmax: float) -> None:
    text = f"""# Scalar/lensing null-space response optimization.
output = tCl,pCl,lCl,mPk
modes = s
lensing = yes
H0 = {H0}
omega_b = {OMEGA_B}
omega_cdm = {OMEGA_CDM}
A_s = {A_S:.12e}
n_s = {N_S}
tau_reio = {TAU_REIO}
N_ur = {N_UR}
N_ncdm = 1
use_ncdm_psd_files = 1
ncdm_psd_filenames = {psd.resolve()}
m_ncdm = {mass_eV}
T_ncdm = {T_NCDM}
deg_ncdm = 1.0
l_max_scalars = {lmax}
P_k_max_h/Mpc = {pkmax}
z_pk = 0
root = {root}
headers = yes
write warnings = yes
"""
    path.write_text(text, encoding="utf-8")


def read_table(path: Path):
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
    labels = {}
    import re
    for line in header:
        for num, name in re.findall(r"(\d+):([^\s]+)", line):
            labels[name.strip().lower()] = int(num) - 1
    return arr, labels, header


def col(labels, *names):
    for name in names:
        if name.lower() in labels:
            return labels[name.lower()]
    return None


def load_pair(d: Path, stem: str):
    lp = d / f"{stem}_plus_00_cl_lensed.dat"
    lm = d / f"{stem}_minus_00_cl_lensed.dat"
    up = d / f"{stem}_plus_00_cl.dat"
    um = d / f"{stem}_minus_00_cl.dat"
    for p in (lp, lm, up, um):
        if not p.exists():
            raise FileNotFoundError(p)
    AL, LL, _ = read_table(lp)
    BL, LM, _ = read_table(lm)
    AU, UL, HU = read_table(up)
    BU, UM, _ = read_table(um)
    if not np.array_equal(AL[:, 0], BL[:, 0]) or not np.array_equal(AU[:, 0], BU[:, 0]):
        raise RuntimeError(f"ell grids differ for {stem}")
    iTT, iEE, iTE = col(LL, "tt"), col(LL, "ee"), col(LL, "te")
    if None in (iTT, iEE, iTE):
        raise RuntimeError(f"Missing T/E columns for {stem}: {LL}")
    jPP = col(UL, "pp", "phiphi", "phi-phi", "phi_phi")
    if jPP is None:
        if AU.shape[1] >= 6:
            jPP = 5
        else:
            raise RuntimeError("Could not identify phi-phi column: " + " | ".join(HU))
    return {
        "ell": AL[:, 0].astype(int),
        "TTp": AL[:, iTT], "TTm": BL[:, iTT],
        "EEp": AL[:, iEE], "EEm": BL[:, iEE],
        "TEp": AL[:, iTE], "TEm": BL[:, iTE],
        "ellu": AU[:, 0].astype(int),
        "PPp": AU[:, jPP], "PPm": BU[:, jPP],
    }


def load_pk_pair(d: Path, stem: str):
    ap, _, _ = read_table(d / f"{stem}_plus_00_pk.dat")
    am, _, _ = read_table(d / f"{stem}_minus_00_pk.dat")
    kp, pp = ap[:, 0], ap[:, 1]
    km, pm = am[:, 0], am[:, 1]
    lo, hi = max(kp.min(), km.min()), min(kp.max(), km.max())
    kp_sel = kp[(kp >= lo) & (kp <= hi)]
    km_sel = km[(km >= lo) & (km <= hi)]
    k = kp_sel if len(kp_sel) >= len(km_sel) else km_sel
    k = k[k > 0]
    if np.any(pp <= 0) or np.any(pm <= 0):
        raise RuntimeError("Non-positive P(k), cannot log-interpolate")
    lpp = np.interp(np.log(k), np.log(kp), np.log(pp))
    lpm = np.interp(np.log(k), np.log(km), np.log(pm))
    return k, np.exp(lpp), np.exp(lpm)


def prepare_probes(d: Path, mass_eV: float, z_match: float, probe_frac: float,
                   lmax: int, pkmax: float):
    d.mkdir(parents=True, exist_ok=True)
    q = np.linspace(0.0, 20.0, 4000)
    f0, weights, basis, N, shapes, y, M = kinetic_objects(q, mass_eV, z_match)
    np.save(d / "q.npy", q)
    np.save(d / "f0.npy", f0)
    np.save(d / "weights.npy", weights)
    np.save(d / "null_shapes.npy", shapes)

    for j, shape in enumerate(shapes):
        fp, fm = f0 + probe_frac * shape, f0 - probe_frac * shape
        if fp.min() <= 0 or fm.min() <= 0:
            raise RuntimeError(f"Probe positivity failed for direction {j}")
        pp, pm = d / f"probe_{j:02d}_plus.dat", d / f"probe_{j:02d}_minus.dat"
        write_psd(pp, q, fp)
        write_psd(pm, q, fm)
        write_scalar_ini(d / f"probe_{j:02d}_plus.ini", pp, d / f"probe_{j:02d}_plus_",
                         mass_eV, lmax, pkmax)
        write_scalar_ini(d / f"probe_{j:02d}_minus.ini", pm, d / f"probe_{j:02d}_minus_",
                         mass_eV, lmax, pkmax)

    meta = {
        "mass_eV": mass_eV, "z_match": z_match, "mass_over_Tnu_at_match": y,
        "probe_fractional_distortion": probe_frac,
        "null_dimension": int(shapes.shape[0]), "basis_dimension": int(basis.shape[0]),
        "moment_matrix_rank": int(np.linalg.matrix_rank(M)),
        "lmax": lmax, "pkmax_h_Mpc": pkmax, "class_commit": CLASS_COMMIT,
        "targets": ["ideal full-sky Gaussian TT+TE+EE S/N", "ideal full-sky phi-phi auto S/N"],
    }
    (d / "scalar_probe_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))


def build_gram(d: Path, lmin: int, lmax: int):
    shapes = np.load(d / "null_shapes.npy")
    nd = shapes.shape[0]
    pairs = [load_pair(d, f"probe_{j:02d}") for j in range(nd)]
    ell = pairs[0]["ell"]
    ellu = pairs[0]["ellu"]
    for p in pairs[1:]:
        if not np.array_equal(ell, p["ell"]) or not np.array_equal(ellu, p["ellu"]):
            raise RuntimeError("Probe ell grids differ")

    T0 = np.mean([0.5 * (p["TTp"] + p["TTm"]) for p in pairs], axis=0)
    E0 = np.mean([0.5 * (p["EEp"] + p["EEm"]) for p in pairs], axis=0)
    X0 = np.mean([0.5 * (p["TEp"] + p["TEm"]) for p in pairs], axis=0)
    P0 = np.mean([0.5 * (p["PPp"] + p["PPm"]) for p in pairs], axis=0)

    dT = np.column_stack([p["TTp"] - p["TTm"] for p in pairs])
    dE = np.column_stack([p["EEp"] - p["EEm"] for p in pairs])
    dX = np.column_stack([p["TEp"] - p["TEm"] for p in pairs])
    dP = np.column_stack([p["PPp"] - p["PPm"] for p in pairs])

    Gs = np.zeros((nd, nd))
    for i, l in enumerate(ell):
        if l < lmin or l > lmax or T0[i] <= 0 or E0[i] <= 0:
            continue
        n = 2.0 * l + 1.0
        cov = np.array([
            [2*T0[i]*T0[i], 2*X0[i]*X0[i], 2*T0[i]*X0[i]],
            [2*X0[i]*X0[i], 2*E0[i]*E0[i], 2*E0[i]*X0[i]],
            [2*T0[i]*X0[i], 2*E0[i]*X0[i], X0[i]*X0[i] + T0[i]*E0[i]],
        ]) / n
        inv = np.linalg.pinv(cov, rcond=1e-12)
        B = np.vstack([dT[i], dE[i], dX[i]])
        Gs += B.T @ inv @ B

    Gp = np.zeros((nd, nd))
    for i, l in enumerate(ellu):
        if l < 8 or l > lmax or P0[i] == 0:
            continue
        v = dP[i] / P0[i]
        Gp += 0.5 * (2.0 * l + 1.0) * np.outer(v, v)

    return Gs, Gp


def optimize(d: Path, final_frac: float, lmin: int, lmax: int,
             random_samples: int, seed: int):
    meta = json.loads((d / "scalar_probe_meta.json").read_text())
    probe_frac = float(meta["probe_fractional_distortion"])
    mass_eV = float(meta["mass_eV"])
    pkmax = float(meta["pkmax_h_Mpc"])
    shapes = np.load(d / "null_shapes.npy")
    q = np.load(d / "q.npy")
    f0 = np.load(d / "f0.npy")
    weights = np.load(d / "weights.npy")
    nd = shapes.shape[0]
    Gs, Gp = build_gram(d, lmin, lmax)

    def top_eig(G):
        w, v = np.linalg.eigh(0.5 * (G + G.T))
        return v[:, np.argmax(w)]

    seeds = [top_eig(Gs), top_eig(Gp)]
    for j in range(nd):
        e = np.zeros(nd); e[j] = 1.0
        seeds.extend([e, -e])

    rng = np.random.default_rng(seed)
    batch = rng.normal(size=(random_samples, nd))
    batch /= np.linalg.norm(batch, axis=1, keepdims=True)

    best = {"scalar": None, "lensing": None}

    def score(c):
        raw = c @ shapes
        maxrel = float(np.max(np.abs(raw) / np.maximum(f0, 1e-300)))
        if not np.isfinite(maxrel) or maxrel <= 0:
            return None
        norm = 1.0 / maxrel
        scale = (final_frac / probe_frac) * norm
        ss = float(scale * np.sqrt(max(c @ Gs @ c, 0.0)))
        sp = float(scale * np.sqrt(max(c @ Gp @ c, 0.0)))
        return norm, ss, sp, raw * norm

    for c in list(seeds) + list(batch):
        s = score(c)
        if s is None:
            continue
        norm, ss, sp, shape = s
        if best["scalar"] is None or ss > best["scalar"]["score"]:
            best["scalar"] = {"score": ss, "coeff": c.copy(), "norm": norm,
                              "cross": sp, "shape": shape.copy()}
        if best["lensing"] is None or sp > best["lensing"]["score"]:
            best["lensing"] = {"score": sp, "coeff": c.copy(), "norm": norm,
                               "cross": ss, "shape": shape.copy()}

    out = dict(meta)
    out.update({
        "final_fractional_distortion_cap": final_frac,
        "optimization_ell_min": lmin, "optimization_ell_max": lmax,
        "random_search_samples": random_samples, "random_seed": seed,
        "predicted_scalar_opt_TTEE_SN": best["scalar"]["score"],
        "predicted_scalar_opt_phiphi_SN": best["scalar"]["cross"],
        "predicted_lensing_opt_phiphi_SN": best["lensing"]["score"],
        "predicted_lensing_opt_TTEE_SN": best["lensing"]["cross"],
        "scalar_opt_coefficients": best["scalar"]["coeff"].tolist(),
        "lensing_opt_coefficients": best["lensing"]["coeff"].tolist(),
    })

    for name, item in best.items():
        shape = item["shape"]
        fp, fm = f0 + final_frac * shape, f0 - final_frac * shape
        if fp.min() <= 0 or fm.min() <= 0:
            raise RuntimeError(f"{name} optimized pair failed positivity")
        stats = pair_stats(f0, fp, fm, q, weights)
        out[f"{name}_pair_stats"] = stats
        pp, pm = d / f"{name}_final_plus.dat", d / f"{name}_final_minus.dat"
        write_psd(pp, q, fp); write_psd(pm, q, fm)
        write_scalar_ini(d / f"{name}_final_plus.ini", pp, d / f"{name}_final_plus_",
                         mass_eV, int(meta["lmax"]), pkmax)
        write_scalar_ini(d / f"{name}_final_minus.ini", pm, d / f"{name}_final_minus_",
                         mass_eV, int(meta["lmax"]), pkmax)

    (d / "scalar_optimization_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def symrel(a, b, floor_frac=1e-10):
    den = 0.5 * (np.abs(a) + np.abs(b))
    floor = max(float(np.nanmax(den)) * floor_frac, 1e-300)
    ok = np.isfinite(a) & np.isfinite(b) & (den > floor)
    r = np.full_like(den, np.nan, dtype=float)
    r[ok] = (a[ok] - b[ok]) / den[ok]
    return r, ok


def cv_auto(ell, a, b, lmin, lmax):
    r, ok = symrel(a, b)
    s = (ell >= lmin) & (ell <= lmax) & ok
    return float(np.sqrt(np.sum(0.5 * (2*ell[s] + 1) * r[s]**2)))


def cv_joint(p, lmin, lmax):
    sn2, used = 0.0, 0
    for l,T1,E1,X1,T2,E2,X2 in zip(p["ell"],p["TTp"],p["EEp"],p["TEp"],p["TTm"],p["EEm"],p["TEm"]):
        if l < lmin or l > lmax:
            continue
        T,E,X = 0.5*(T1+T2),0.5*(E1+E2),0.5*(X1+X2)
        if T <= 0 or E <= 0:
            continue
        cov = np.array([[2*T*T,2*X*X,2*T*X],
                        [2*X*X,2*E*E,2*E*X],
                        [2*T*X,2*E*X,X*X+T*E]])/(2*l+1)
        dv = np.array([T1-T2,E1-E2,X1-X2])
        val = float(dv @ np.linalg.pinv(cov, rcond=1e-12) @ dv)
        if np.isfinite(val) and val >= 0:
            sn2 += val; used += 1
    return float(np.sqrt(sn2)), used


def actual_metrics(d: Path, stem: str, lmin: int, lmax: int, pkmax: float):
    p = load_pair(d, stem)
    out = {}
    for key,a,b,ell,lo in [
        ("TT",p["TTp"],p["TTm"],p["ell"],lmin),
        ("EE",p["EEp"],p["EEm"],p["ell"],lmin),
        ("phiphi",p["PPp"],p["PPm"],p["ellu"],8),
    ]:
        r, ok = symrel(a,b)
        s = (ell >= lo) & (ell <= lmax) & ok
        ii = np.where(s)[0]
        im = ii[np.argmax(np.abs(r[s]))]
        out[key] = {
            "max_abs_relative_difference_percent": float(100*abs(r[im])),
            "signed_relative_difference_percent_at_max": float(100*r[im]),
            "location": int(ell[im]),
            "ideal_fullsky_cv_SN_auto_only": cv_auto(ell,a,b,lo,lmax),
        }
    sj, used = cv_joint(p, lmin, lmax)
    out["TTEE_TE_combined"] = {"ideal_fullsky_gaussian_cv_SN": sj, "multipoles_used": used}
    k, pp, pm = load_pk_pair(d, stem)
    r, ok = symrel(pp, pm)
    s = (k >= 1e-3) & (k <= pkmax) & ok
    ii = np.where(s)[0]
    im = ii[np.argmax(np.abs(r[s]))]
    out["Pk_z0"] = {
        "max_abs_relative_difference_percent": float(100*abs(r[im])),
        "signed_relative_difference_percent_at_max": float(100*r[im]),
        "location_k_h_Mpc": float(k[im]),
    }
    return out


def summarize_candidates(d: Path, lmin: int, lmax: int, pkmax: float):
    results = {}
    for name in ("scalar", "lensing"):
        results[name] = actual_metrics(d, f"{name}_final", lmin, lmax, pkmax)
    def strength(m):
        return max(m["TTEE_TE_combined"]["ideal_fullsky_gaussian_cv_SN"],
                   m["phiphi"]["ideal_fullsky_cv_SN_auto_only"])
    winner = max(results, key=lambda n: strength(results[n]))
    shutil.copy2(d / f"{winner}_final_plus.dat", d / "winner_plus.dat")
    shutil.copy2(d / f"{winner}_final_minus.dat", d / "winner_minus.dat")
    out = {"moderate_precision_metrics": results, "winner": winner,
           "winner_strength_ideal_SN": strength(results[winner])}
    (d / "scalar_candidate_actual_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def prepare_winner_highprec(d: Path, mass_eV: float, lmax: int, pkmax: float):
    hp = d / "highprec"
    hp.mkdir(exist_ok=True)
    write_scalar_ini(hp / "winner_plus.ini", d / "winner_plus.dat", hp / "winner_plus_",
                     mass_eV, lmax, pkmax)
    write_scalar_ini(hp / "winner_minus.ini", d / "winner_minus.dat", hp / "winner_minus_",
                     mass_eV, lmax, pkmax)


def summarize_winner_highprec(d: Path, lmin: int, lmax: int, pkmax: float):
    hp = d / "highprec"
    m = actual_metrics(hp, "winner", lmin, lmax, pkmax)
    cand = json.loads((d / "scalar_candidate_actual_summary.json").read_text())
    out = {"winner": cand["winner"], "high_precision_metrics": m}
    (d / "scalar_winner_highprec_summary.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=[
        "prepare-probes", "optimize", "summarize-candidates",
        "prepare-winner-highprec", "summarize-winner-highprec"])
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mass-eV", type=float, default=0.10)
    ap.add_argument("--z-match", type=float, default=1100.0)
    ap.add_argument("--probe-frac", type=float, default=0.08)
    ap.add_argument("--final-frac", type=float, default=0.30)
    ap.add_argument("--lmin", type=int, default=30)
    ap.add_argument("--lmax", type=int, default=800)
    ap.add_argument("--pkmax", type=float, default=0.30)
    ap.add_argument("--random-samples", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=314159)
    args = ap.parse_args()
    if args.mode == "prepare-probes":
        prepare_probes(args.outdir,args.mass_eV,args.z_match,args.probe_frac,args.lmax,args.pkmax)
    elif args.mode == "optimize":
        optimize(args.outdir,args.final_frac,args.lmin,args.lmax,args.random_samples,args.seed)
    elif args.mode == "summarize-candidates":
        summarize_candidates(args.outdir,args.lmin,args.lmax,args.pkmax)
    elif args.mode == "prepare-winner-highprec":
        prepare_winner_highprec(args.outdir,args.mass_eV,args.lmax,args.pkmax)
    else:
        summarize_winner_highprec(args.outdir,args.lmin,args.lmax,args.pkmax)


if __name__ == "__main__":
    main()
