#!/usr/bin/env python3
"""DESI DR1 luminosity-rank odd-octupole consistency test.

We measure the dipole and octupole on the same pair sample and evaluate the
18-component octupole against the within-stratum permutation distribution.
We do not fit a physical wake template to the octupole. The dipole comparison
uses the recorded 32-permutation reference from this control analysis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

import desi_dr1_phase7_lss as p
import desi_dr1_phase7_multitracer_fullsample as fs


def legendre_weight(mu: np.ndarray, ell: int) -> np.ndarray:
    mu = np.asarray(mu, float)
    if ell == 1:
        return 3.0 * mu
    if ell == 3:
        return 3.5 * (5.0 * mu**3 - 3.0 * mu)
    raise ValueError("only ell=1 and ell=3 are supported")


def hist_pair_multipole(tab, A, B, markA, markB, nb, ell, drop=None):
    keep = np.ones(len(tab["a"]), bool)
    if drop is not None:
        keep &= (A["jk"][tab["a"]] != drop) & (B["jk"][tab["b"]] != drop)
    a = tab["a"][keep]
    b = tab["b"][keep]
    ib = tab["bin"][keep]
    pw = A["w"][a] * B["w"][b] * tab["imp"][keep]
    h0 = np.bincount(ib, weights=pw, minlength=nb).astype(float)
    odd = pw * (markB[b] - markA[a]) * legendre_weight(tab["mu"][keep], ell)
    hell = np.bincount(ib, weights=odd, minlength=nb).astype(float)
    return h0, hell


def estimate_block_multipole(block, ell, markD=None, drop=None):
    D, R = block["D"], block["R"]
    markD = D["mark"] if markD is None else markD
    markR = R["mark"]
    nb = block["nb"]

    dd0, ddell = hist_pair_multipole(block["DD"], D, D, markD, markD, nb, ell, drop)
    dr0, drell = hist_pair_multipole(block["DR"], D, R, markD, markR, nb, ell, drop)
    rr0, rrell = hist_pair_multipole(block["RR"], R, R, markR, markR, nb, ell, drop)

    WD, NDD = fs.object_norm(D, drop)
    WR, NRR = fs.object_norm(R, drop)
    NDR = max(WD * WR, 1e-300)

    dd0 /= NDD
    ddell /= NDD
    dr0 /= NDR
    drell /= NDR
    rr0 /= NRR
    rrell /= NRR

    den = np.maximum(rr0, 1e-300)
    xi0 = (dd0 - 2.0 * dr0 + rr0) / den
    xiell = (ddell - 2.0 * drell + rrell) / den
    return xi0, xiell


def vector_multipole(blocks, ell, marks=None, drop=None):
    x0, xell = [], []
    for iz, block in enumerate(blocks):
        md = None if marks is None else marks[iz]
        a, c = estimate_block_multipole(block, ell, md, drop)
        x0.extend(a.tolist())
        xell.extend(c.tolist())
    return np.asarray(x0), np.asarray(xell)


def regularized_jackknife_cov(jk):
    jk = np.asarray(jk, float)
    n, pdim = jk.shape
    jm = jk.mean(axis=0)
    raw = (n - 1.0) / n * ((jk - jm).T @ (jk - jm))
    ev = np.linalg.eigvalsh(raw)
    maxev = float(max(ev[-1], 1e-300))
    ridge = max(maxev * 1e-7, float(np.median(np.diag(raw))) * 1e-6, 1e-14)
    cov = raw + np.eye(pdim) * ridge
    cinv = np.linalg.pinv(cov, rcond=1e-10)
    rank = int(np.linalg.matrix_rank(raw, tol=maxev * 1e-10))
    cond = float(np.linalg.cond(cov))
    return raw, cov, cinv, rank, cond, ridge


def null_test(data, null, cov, cinv):
    null = np.asarray(null, float)
    null_mean = null.mean(axis=0)
    y = np.asarray(data, float) - null_mean
    tdata = float(y @ cinv @ y)
    tperm = np.asarray([(v - null_mean) @ cinv @ (v - null_mean) for v in null], float)
    pemp = float((1 + np.sum(tperm >= tdata)) / (len(tperm) + 1))
    pasym = float(chi2.sf(tdata, len(y)))
    sig = np.sqrt(np.maximum(np.diag(cov), 0.0))
    zbin = np.divide(y, sig, out=np.zeros_like(y), where=sig > 0)
    imax = int(np.argmax(np.abs(zbin)))
    return {
        "null_mean": null_mean,
        "null_corrected": y,
        "mahalanobis_data": tdata,
        "permutation_mahalanobis": tperm,
        "empirical_pvalue": pemp,
        "asymptotic_chi2_pvalue_diagnostic": pasym,
        "max_abs_single_bin_z_diagnostic": float(abs(zbin[imax])),
        "max_abs_single_bin_index": imax,
        "bin_z": zbin,
    }


def reference_dipole_check(path, y):
    pth = Path(path)
    if not pth.is_file():
        return {"available": False, "path": str(pth)}
    ref = np.genfromtxt(pth, delimiter=",", names=True)
    if "xi1_null_corrected" not in ref.dtype.names:
        raise RuntimeError(f"{pth} lacks xi1_null_corrected")
    r = np.asarray(ref["xi1_null_corrected"], float)
    if len(r) != len(y):
        raise RuntimeError(f"reference dipole length mismatch: {len(r)} != {len(y)}")
    diff = np.asarray(y, float) - r
    corr = float(np.corrcoef(y, r)[0, 1]) if np.std(y) > 0 and np.std(r) > 0 else float("nan")
    return {
        "available": True,
        "path": str(pth),
        "max_abs_difference": float(np.max(np.abs(diff))),
        "rms_difference": float(np.sqrt(np.mean(diff * diff))),
        "correlation": corr,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", nargs="+", required=True)
    ap.add_argument("--random", nargs="+", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--reference-dipole", default="source_data/phase7_octupole_control/reference_dipole_32perm.csv")
    ap.add_argument("--zmin", type=float, default=0.10)
    ap.add_argument("--zmax", type=float, default=0.40)
    ap.add_argument("--dz-proxy", type=float, default=0.02)
    ap.add_argument("--ntracer", type=int, default=5)
    ap.add_argument("--analysis-z-edges", default="0.10,0.20,0.30,0.40")
    ap.add_argument("--sep-edges", default="20,40,60,80,100,120,140")
    ap.add_argument("--random-factor", type=float, default=2.0)
    ap.add_argument("--neighbors-per-anchor", type=int, default=48)
    ap.add_argument("--theta-min-deg", type=float, default=0.05)
    ap.add_argument("--jackknife", type=int, default=30)
    ap.add_argument("--permutations", type=int, default=32)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    edges = np.asarray([float(x) for x in args.sep_edges.split(",")])
    s = 0.5 * (edges[:-1] + edges[1:])
    nb = len(s)
    zedges = np.asarray([float(x) for x in args.analysis_z_edges.split(",")])
    nz = len(zedges) - 1
    if args.ntracer < 3:
        raise ValueError("ntracer must be >=3")

    catD = p.read_catalog(args.data, "data")
    catR = p.read_catalog(args.random, "random")

    iD, lD, mD, sD, proxy_info = fs.make_data_proxy(
        catD, args.zmin, args.zmax, args.dz_proxy, args.ntracer
    )
    iR, lR, mR, sR = fs.make_random_proxy(
        catR, args.zmin, args.zmax, args.dz_proxy, args.ntracer,
        proxy_info, args.random_factor, rng
    )

    regD, njd = p.sky_jackknife_regions(catD["ra"][iD], catD["region"][iD], args.jackknife)
    regR, njr = p.sky_jackknife_regions(catR["ra"][iR], catR["region"][iR], args.jackknife)
    nj = min(njd, njr)
    if nj <= nz * nb + 2:
        raise RuntimeError("insufficient jackknife regions")

    blocks, zmeta = [], []
    for iz, (lo, hi) in enumerate(zip(zedges[:-1], zedges[1:])):
        qd = np.where((catD["z"][iD] >= lo) & (catD["z"][iD] < hi))[0]
        qr = np.where((catR["z"][iR] >= lo) & (catR["z"][iR] < hi))[0]
        D = fs.prepare(catD, iD[qd], mD[qd], sD[qd], regD[qd])
        R = fs.prepare(catR, iR[qr], mR[qr], sR[qr], regR[qr])
        if min(len(D["w"]), len(R["w"])) < 100:
            raise RuntimeError(f"too few objects in z bin {lo}-{hi}")

        DD = fs.sample_pairs(
            D, D, edges, args.neighbors_per_anchor, rng,
            args.theta_min_deg, auto=True, label=f"DD_z{iz}"
        )
        DR = fs.sample_pairs(
            D, R, edges, args.neighbors_per_anchor, rng,
            args.theta_min_deg, auto=False, label=f"DR_z{iz}"
        )
        RR = fs.sample_pairs(
            R, R, edges, args.neighbors_per_anchor, rng,
            args.theta_min_deg, auto=True, label=f"RR_z{iz}"
        )
        blocks.append({"D": D, "R": R, "DD": DD, "DR": DR, "RR": RR, "nb": nb})
        ze = float(np.average(D["z"], weights=D["w"]))
        zmeta.append({
            "zlo": float(lo), "zhi": float(hi), "z_effective": ze,
            "N_data": int(len(D["w"])), "N_random": int(len(R["w"]))
        })

    xi0, xi1 = vector_multipole(blocks, 1)
    xi0_b, xi3 = vector_multipole(blocks, 3)
    if not np.allclose(xi0, xi0_b, rtol=0.0, atol=1e-14):
        raise RuntimeError("internal xi0 mismatch between ell=1 and ell=3 paths")
    pdim = len(xi1)

    jk1, jk3 = [], []
    for j in range(nj):
        _, v1 = vector_multipole(blocks, 1, drop=j)
        _, v3 = vector_multipole(blocks, 3, drop=j)
        jk1.append(v1)
        jk3.append(v3)
        print("JK_OCT", j, json.dumps({"xi1": v1.tolist(), "xi3": v3.tolist()}))

    jk1 = np.asarray(jk1)
    jk3 = np.asarray(jk3)
    raw1, cov1, cinv1, rank1, cond1, ridge1 = regularized_jackknife_cov(jk1)
    raw3, cov3, cinv3, rank3, cond3, ridge3 = regularized_jackknife_cov(jk3)

    null1, null3 = [], []
    for ip in range(args.permutations):
        pm = fs.permuted_marks(blocks, rng)
        _, v1 = vector_multipole(blocks, 1, pm)
        _, v3 = vector_multipole(blocks, 3, pm)
        null1.append(v1)
        null3.append(v3)
        print("PERM_OCT", ip, json.dumps({"xi1": v1.tolist(), "xi3": v3.tolist()}))

    null1 = np.asarray(null1)
    null3 = np.asarray(null3)
    test1 = null_test(xi1, null1, cov1, cinv1)
    test3 = null_test(xi3, null3, cov3, cinv3)

    dipole_ref = reference_dipole_check(args.reference_dipole, test1["null_corrected"])

    rows = []
    k = 0
    sig1 = np.sqrt(np.maximum(np.diag(cov1), 0.0))
    sig3 = np.sqrt(np.maximum(np.diag(cov3), 0.0))
    for meta in zmeta:
        for si in s:
            rows.append((
                meta["zlo"], meta["zhi"], meta["z_effective"], si,
                xi0[k],
                xi1[k], test1["null_mean"][k], test1["null_corrected"][k], sig1[k],
                xi3[k], test3["null_mean"][k], test3["null_corrected"][k], sig3[k],
            ))
            k += 1

    np.savetxt(
        out / "data_vector_octupole_control.csv",
        np.asarray(rows), delimiter=",",
        header=(
            "zlo,zhi,z_effective,s_Mpc_over_h,xi0_proxy,"
            "xi1_proxy_odd,xi1_permutation_null_mean,xi1_null_corrected,xi1_jackknife_sigma,"
            "xi3_proxy_odd,xi3_permutation_null_mean,xi3_null_corrected,xi3_jackknife_sigma"
        ),
        comments=""
    )
    np.savetxt(out / "jackknife_covariance_dipole.csv", cov1, delimiter=",")
    np.savetxt(out / "jackknife_covariance_octupole.csv", cov3, delimiter=",")
    np.savetxt(out / "permutation_vectors_dipole.csv", null1, delimiter=",")
    np.savetxt(out / "permutation_vectors_octupole.csv", null3, delimiter=",")

    p3 = float(test3["empirical_pvalue"])
    status = "CONSISTENT_WITH_NULL" if p3 >= 0.05 else "OCTUPOLE_EXCESS_INVESTIGATE_SYSTEMATICS"

    summary = {
        "scope": (
            "Odd-octupole consistency test on the DESI DR1 BGS "
            "five-tracer luminosity-rank sample. The ell=3 octupole is tested "
            "against its permutation distribution; no wake template is fitted to ell=3."
        ),
        "predeclared_before_data_inspection": {
            "primary_observable": "odd octupole ell=3",
            "basis": "7 P_3(mu) = (7/2)(5 mu^3 - 3 mu)",
            "primary_statistic": "global permutation-calibrated Mahalanobis distance",
            "decision_rule": (
                "empirical p >= 0.05: consistent with null; empirical p < 0.05: "
                "investigate survey/estimator systematics before any physical interpretation"
            ),
            "dipole_role": "reference-vector reproduction test",
            "no_octupole_wake_template_fit": True
        },
        "configuration": {
            "seed": int(args.seed),
            "ntracer": int(args.ntracer),
            "zmin": float(args.zmin),
            "zmax": float(args.zmax),
            "dz_proxy": float(args.dz_proxy),
            "analysis_z_edges": zedges.tolist(),
            "separation_edges_Mpc_over_h": edges.tolist(),
            "random_factor": float(args.random_factor),
            "neighbors_per_anchor": int(args.neighbors_per_anchor),
            "theta_min_deg": float(args.theta_min_deg),
            "jackknife_regions": int(nj),
            "permutations": int(args.permutations),
            "data_count": int(len(iD)),
            "random_count": int(len(iR))
        },
        "dipole_reproduction": {
            "raw_covariance_rank": rank1,
            "regularized_condition_number": cond1,
            "ridge_added": ridge1,
            "global_empirical_pvalue": float(test1["empirical_pvalue"]),
            "global_mahalanobis": float(test1["mahalanobis_data"]),
            "reference_check": dipole_ref
        },
        "octupole_null_test": {
            "status": status,
            "raw_covariance_rank": rank3,
            "regularized_condition_number": cond3,
            "ridge_added": ridge3,
            "global_mahalanobis": float(test3["mahalanobis_data"]),
            "global_empirical_pvalue": p3,
            "permutation_pvalue_resolution": float(1.0 / (len(null3) + 1)),
            "asymptotic_chi2_pvalue_diagnostic": float(test3["asymptotic_chi2_pvalue_diagnostic"]),
            "max_abs_single_bin_z_diagnostic": float(test3["max_abs_single_bin_z_diagnostic"]),
            "max_abs_single_bin_index": int(test3["max_abs_single_bin_index"])
        },
        "proxy_strata": proxy_info,
        "analysis_scope": (
            "The octupole is a separate odd-multipole consistency observable. "
            "The dipole reference uses the corresponding 32-permutation realization; "
            "the primary DESI dipole coefficient uses 256 permutations."
        )
    }
    (out / "summary_octupole_control.json").write_text(json.dumps(summary, indent=2) + "\n")

    print("PHASE7_OCTUPOLE_CONTROL", json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
