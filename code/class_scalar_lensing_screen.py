#!/usr/bin/env python3
"""Screen an existing stress-energy-matched relic pair in scalar CMB/lensing/P(k).

This is a channel-screening calculation, not a data fit. It keeps the optimized
m_ncdm=0.10 eV F_+/F_- pair fixed and asks whether scalar observables provide a
larger separation than the tensor B-mode channel used in the manuscript.
"""
from __future__ import annotations
import argparse, json, re
from pathlib import Path
import numpy as np

H0 = 67.36
OMEGA_B = 0.02237
OMEGA_CDM = 0.1200
A_S = float(np.exp(3.044) * 1e-10)
N_S = 0.9649
TAU_REIO = 0.054
N_UR = 2.0328
T_NCDM = 0.71611
CLASS_COMMIT = "e85808324f51fc694d12e3ed7439552a3c3f9540"


def write_ini(path: Path, psd: Path, root: Path, mass_eV: float, lmax: int, pkmax: float) -> None:
    text = f"""# Scalar/lensing channel screen for a fixed stress-energy-matched relic pair.
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
    for line in header:
        for num, name in re.findall(r"(\d+):([^\s]+)", line):
            labels[name.strip().lower()] = int(num) - 1
    return arr, labels, header


def col(labels, *names):
    for name in names:
        if name.lower() in labels:
            return labels[name.lower()]
    return None


def symmetric_relative(a, b, floor_frac=1e-10):
    den = 0.5 * (np.abs(a) + np.abs(b))
    floor = max(float(np.nanmax(den)) * floor_frac, 1e-300)
    mask = np.isfinite(a) & np.isfinite(b) & (den > floor)
    rel = np.full_like(den, np.nan, dtype=float)
    rel[mask] = (a[mask] - b[mask]) / den[mask]
    return rel, mask


def max_metric(x, rel, sel):
    ok = sel & np.isfinite(rel)
    if not np.any(ok):
        return {"max_abs_relative_difference_percent": None, "location": None}
    ii = np.where(ok)[0]
    j = ii[np.argmax(np.abs(rel[ok]))]
    return {
        "max_abs_relative_difference_percent": float(100.0 * abs(rel[j])),
        "signed_relative_difference_percent_at_max": float(100.0 * rel[j]),
        "location": float(x[j]),
    }


def cv_sn_auto(ell, a, b, lmin, lmax):
    rel, mask = symmetric_relative(a, b)
    sel = (ell >= lmin) & (ell <= lmax) & mask
    return float(np.sqrt(np.sum(0.5 * (2.0 * ell[sel] + 1.0) * rel[sel] ** 2)))


def cv_sn_te_block(ell, tp, ep, xp, tm, em, xm, lmin, lmax):
    """Ideal full-sky Gaussian T/E covariance S/N, no noise or nuisance parameters."""
    sn2 = 0.0
    used = 0
    for l, T1, E1, X1, T2, E2, X2 in zip(ell, tp, ep, xp, tm, em, xm):
        if l < lmin or l > lmax:
            continue
        T = 0.5 * (T1 + T2); E = 0.5 * (E1 + E2); X = 0.5 * (X1 + X2)
        if T <= 0 or E <= 0 or not np.all(np.isfinite([T,E,X])):
            continue
        n = 2.0 * l + 1.0
        cov = np.array([
            [2*T*T, 2*X*X, 2*T*X],
            [2*X*X, 2*E*E, 2*E*X],
            [2*T*X, 2*E*X, X*X + T*E],
        ]) / n
        d = np.array([T1-T2, E1-E2, X1-X2])
        try:
            val = float(d @ np.linalg.pinv(cov, rcond=1e-12) @ d)
        except np.linalg.LinAlgError:
            continue
        if np.isfinite(val) and val >= 0:
            sn2 += val; used += 1
    return float(np.sqrt(max(sn2, 0.0))), used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["prepare", "summarize"], required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--plus-psd", type=Path)
    ap.add_argument("--minus-psd", type=Path)
    ap.add_argument("--mass-eV", type=float, default=0.10)
    ap.add_argument("--lmax", type=int, default=2500)
    ap.add_argument("--pkmax", type=float, default=1.0)
    args = ap.parse_args()
    d = args.outdir
    d.mkdir(parents=True, exist_ok=True)

    if args.mode == "prepare":
        if args.plus_psd is None or args.minus_psd is None:
            raise SystemExit("--plus-psd and --minus-psd are required")
        pplus = d / "final_plus.dat"; pminus = d / "final_minus.dat"
        pplus.write_bytes(args.plus_psd.read_bytes()); pminus.write_bytes(args.minus_psd.read_bytes())
        write_ini(d/"scalar_plus.ini", pplus, d/"scalar_plus_", args.mass_eV, args.lmax, args.pkmax)
        write_ini(d/"scalar_minus.ini", pminus, d/"scalar_minus_", args.mass_eV, args.lmax, args.pkmax)
        meta = {"mass_eV":args.mass_eV,"lmax":args.lmax,"pkmax_h_Mpc":args.pkmax,"class_commit":CLASS_COMMIT,
                "pair_provenance":"GitHub Actions run 33857363622 artifact class-response-mass-010",
                "purpose":"scalar/lensing/P(k) screen of the fixed BB-optimized pair; not a data fit"}
        (d/"screen_meta.json").write_text(json.dumps(meta,indent=2)+"\n")
        print(json.dumps(meta,indent=2)); return

    lp = d/"scalar_plus_cl_lensed.dat"; lm = d/"scalar_minus_cl_lensed.dat"
    up = d/"scalar_plus_cl.dat"; um = d/"scalar_minus_cl.dat"
    pp = d/"scalar_plus_pk.dat"; pm = d/"scalar_minus_pk.dat"
    for p in [lp,lm,up,um,pp,pm]:
        if not p.exists(): raise FileNotFoundError(p)

    AL, LL, _ = read_table(lp); BL, LM, _ = read_table(lm)
    AU, UL, HU = read_table(up); BU, UM, HM = read_table(um)
    APK, _, _ = read_table(pp); BPK, _, _ = read_table(pm)
    if not np.array_equal(AL[:,0], BL[:,0]) or not np.array_equal(AU[:,0], BU[:,0]):
        raise RuntimeError("ell grids differ")
    ell = AL[:,0].astype(int)
    iTT=col(LL,"tt"); iEE=col(LL,"ee"); iTE=col(LL,"te")
    if None in (iTT,iEE,iTE): raise RuntimeError(f"Missing T/E columns in lensed file: {LL}")
    jPP=col(UL,"pp","phiphi","phi-phi","phi_phi")
    if jPP is None:
        if AU.shape[1] >= 6:
            jPP=5
        else:
            raise RuntimeError("Could not identify phi-phi column. Header: " + " | ".join(HU))

    TTp,EEp,TEp=AL[:,iTT],AL[:,iEE],AL[:,iTE]
    TTm,EEm,TEm=BL[:,iTT],BL[:,iEE],BL[:,iTE]
    ellu=AU[:,0].astype(int); PPp,PPm=AU[:,jPP],BU[:,jPP]

    metrics={}
    for name,a,b,lo,hi in [("TT",TTp,TTm,30,2500),("EE",EEp,EEm,30,2500),("phiphi",PPp,PPm,8,2000)]:
        x=ell if name!="phiphi" else ellu
        rel,mask=symmetric_relative(a,b)
        sel=(x>=lo)&(x<=hi)&mask
        metrics[name]=max_metric(x,rel,sel)
        metrics[name]["ideal_fullsky_cv_SN_auto_only"] = cv_sn_auto(x,a,b,lo,hi)

    sn_te,nused=cv_sn_te_block(ell,TTp,EEp,TEp,TTm,EEm,TEm,30,2500)
    metrics["TTEE_TE_combined"]={"ideal_fullsky_gaussian_cv_SN":sn_te,"multipoles_used":nused,
        "note":"optimistic screening statistic; no instrument noise, foregrounds or parameter/nuisance degeneracies"}

    denom=np.sqrt(np.maximum(0.5*(TTp+TTm),1e-300)*np.maximum(0.5*(EEp+EEm),1e-300))
    te_norm=(TEp-TEm)/denom
    sel=(ell>=30)&(ell<=2500)&np.isfinite(te_norm)
    metrics["TE"]={"max_abs_delta_TE_over_sqrt_TT_EE_percent":float(100*np.max(np.abs(te_norm[sel]))),
                   "ell_at_max":int(ell[np.where(sel)[0][np.argmax(np.abs(te_norm[sel]))]])}

    if APK.shape != BPK.shape or not np.allclose(APK[:,0],BPK[:,0],rtol=0,atol=1e-14):
        raise RuntimeError("P(k) grids differ")
    k=APK[:,0]; Pp=APK[:,1]; Pm=BPK[:,1]
    rpk,mpk=symmetric_relative(Pp,Pm)
    sel=(k>=1e-3)&(k<=args.pkmax)&mpk
    metrics["Pk_z0"]=max_metric(k,rpk,sel)
    metrics["Pk_z0"]["k_units"]="h/Mpc (CLASS output for P_k_max_h/Mpc input)"

    out={"meta":json.loads((d/"screen_meta.json").read_text()),"metrics":metrics}
    (d/"scalar_lensing_screen_summary.json").write_text(json.dumps(out,indent=2)+"\n")
    np.savetxt(d/"scalar_lensed_pair.csv",np.column_stack([ell,TTp,TTm,EEp,EEm,TEp,TEm]),delimiter=",",
               header="ell,TT_plus,TT_minus,EE_plus,EE_minus,TE_plus,TE_minus",comments="")
    np.savetxt(d/"scalar_phiphi_pair.csv",np.column_stack([ellu,PPp,PPm]),delimiter=",",header="ell,phiphi_plus,phiphi_minus",comments="")
    np.savetxt(d/"scalar_pk_z0_pair.csv",np.column_stack([k,Pp,Pm]),delimiter=",",header="k_h_Mpc,P_plus,P_minus",comments="")
    print(json.dumps(out,indent=2))

if __name__ == "__main__": main()
