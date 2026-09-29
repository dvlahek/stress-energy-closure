#!/usr/bin/env python3
"""E28R1 posthoc registered finer grid: original E28 conditional 3D force ONLY.

First E28 had unresolved signed Fpm contrast. Never quietly accept its sign.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil as old

PRO=ROOT/"source_data/eboss_dr16_a03_e28r1_posthoc_signed_recoil_contrast_and_band_convergence_protocol_2026-09-29.json"
PRO_BLOB="5d0a0ac96364916e95328c0f9ef6d30fb528c229"
RUNNER_BLOB="c02dc66f4c235427d6f3e7dac85b5d617ed7717e"
GRID=((24,257),(48,513))
STATES=("FD","plus","minus")

def need(ok,msg):
    if not ok:
        raise ValueError("E28R1_POSTHOC_GRID_STOP "+msg)

def gitblob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def sources():
    need(gitblob(PRO.read_bytes())==PRO_BLOB,"R1 protocol changed")
    p=json.loads(PRO.read_bytes())
    need(p["original_E28_run_id"]==36601783674 and
         p["original_E28_code_blob"]==RUNNER_BLOB and
         p["new_grid"]["coarse_Simpson_intervals_PER_segment"]==24 and
         p["new_grid"]["fine_Simpson_intervals_PER_segment"]==48 and
         p["new_grid"]["coarse_z_nodes"]==257 and
         p["new_grid"]["fine_z_nodes"]==513 and
         all(p["guards"].values()),"R1 posthoc frozen grid or hard STOP drift")
    need(gitblob(Path(old.__file__).read_bytes())==RUNNER_BLOB,
         "immutable original E28 runner changed")
    _,base,profiles=old.sources()  # all original E8/E9/B1/E27 plus DM14 pins
    return p,base,profiles

def grid_source(prep,background,nseg,nz):
    k=old.k_grid(nseg)
    z,H,hf,ang,dx=old.kernels_for_state(prep,background,k,nz)
    return {"k":k,"z":z,"H":H,"hf":hf,"ang":ang,"dx":dx}

def one_halo(packed,nseg,halo,profile,rho):
    M0,alpha=halo
    k,z,H,hf,ang=[packed[q] for q in ("k","z","H","hf","ang")]
    c=profile["c200c_population_median_DM14"]
    rs=profile["r_s_comoving_Mpc_DM14"]
    nodes,weights=np.polynomial.legendre.leggauss(old.NFW_N)
    u=old.nfw_u(k,rs,c,nodes,weights)
    growth=np.exp(-alpha*(z-old.old.ZOBS)/(1+old.old.ZOBS))
    integ=np.trapz(hf*ang*growth[None,:]/H[None,:],z,axis=1)
    pref=-4*math.pi*old.old.G/(old.old.C**2)*(old.old.MNU/old.old.TNU0)*(M0/k)
    inlog=k*k*pref*integ*u*u
    bands=old.simpson_bands(k,inlog,nseg)
    return (old.old.G/(1+old.old.ZOBS)*rho/math.pi)*bands

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    p,base,profiles=sources()
    _,preps,backgrounds,halos=base
    k0=old.k_grid(24)
    k1=old.k_grid(48)
    need(len(k1)==2*len(k0)-1 and
         np.allclose(k1[::2],k0,rtol=1e-14,atol=0),
         "new R1 k grid must nest the registered original E28 fine grid")
    data={}
    for state in STATES:
        rho=old.q_mass_density(preps[state],old.old.ZOBS)
        packed={(nseg,nz):grid_source(preps[state],backgrounds[state],nseg,nz)
                for nseg,nz in GRID}
        cases={}
        for i,halo in enumerate(halos):
            profile=profiles[i//2]
            co=one_halo(packed[GRID[0]],GRID[0][0],halo,profile,rho)
            fi=one_halo(packed[GRID[1]],GRID[1][0],halo,profile,rho)
            need(np.all(np.isfinite(co)) and np.all(np.isfinite(fi)) and
                 np.all(co<0) and np.all(fi<0),
                 "R1 nonfinite or wrong signed conditional recoil")
            gap=np.abs(fi-co)/np.maximum(np.abs(fi),1e-80)
            need(float(gap[-1])<p["acceptance"]["physical_FD_each_band"].count(
                 "0.01-1")*1e2, "impossible internal finite band QA")
            contributions=np.diff(np.r_[0.,fi])
            cases[str(i)]={"alpha":halo[1],"M0_Msun":halo[0],
                 "fine_prefix_accel_km2_s2_per_Mpc":list(map(float,fi)),
                 "coarse_prefix_accel_km2_s2_per_Mpc":list(map(float,co)),
                 "each_fine_band_accel_km2_s2_per_Mpc":list(map(float,contributions)),
                 "last_band_fraction_of_fine_accel":float(contributions[-1]/fi[-1]),
                 "coarse_fine_prefix_relative_gap":list(map(float,gap))}
        data[state]={"rho_nu_rest_Msun_Mpc3":rho,"cases":cases}
        print("E28R1_STATE",state,"FIRST_HALO_TOTAL_BAND",
              format(cases["0"]["fine_prefix_accel_km2_s2_per_Mpc"][-1],".15g"),
              "LAST_BAND_FRACTION",
              format(cases["0"]["last_band_fraction_of_fine_accel"],".7g"),
              flush=True)
    comparisons={}
    for i in range(4):
        def val(state,grid):
            return data[state]["cases"][str(i)][grid][-1]
        fine=val("minus","fine_prefix_accel_km2_s2_per_Mpc")-val("plus","fine_prefix_accel_km2_s2_per_Mpc")
        coarse=val("minus","coarse_prefix_accel_km2_s2_per_Mpc")-val("plus","coarse_prefix_accel_km2_s2_per_Mpc")
        gap=abs(fine-coarse)/max(abs(fine),1e-80)
        resolved=(fine!=0 and coarse!=0 and (fine>0)==(coarse>0) and gap<.01)
        comparisons[str(i)]={"fine_Fminus_minus_Fplus_accel":fine,
            "coarse_Fminus_minus_Fplus_accel":coarse,
            "signed_difference_relative_to_FD":fine/val("FD","fine_prefix_accel_km2_s2_per_Mpc"),
            "coarse_vs_fine_gap_over_signed_difference":gap,
            "signed_Fpm_difference_numerically_resolved":resolved}
        print("E28R1_SIGNED_FPM",i,"REL_FD",
              format(comparisons[str(i)]["signed_difference_relative_to_FD"],".12g"),
              "SIGNED_REL_GRID_GAP",format(gap,".9g"),
              "RESOLVED",resolved,flush=True)
    # No claim of all-k convergence: explicit uncomputed high-k tail is unknown.
    print("E28R1_ALL_K_NOT_CALCULATED_ZERO_INCOMING_WAKE_AT_Z1",flush=True)
    report={"stage":p["stage"],"registered_protocol_blob":PRO_BLOB,
        "original_E28_CI":p["original_E28_run_id"],
        "original_E28_source_runner_git_blob":RUNNER_BLOB,
        "conditioned_original_3F_4halo_bands":data,"signed_Fpm":comparisons,
        "only_label_contrast_resolved_if_registered_one_percent_pass":True,
        "k_band_comoving_Mpc_inverse":list(old.BANDS),
        "no_all_k_drag_no_real_halo_history_no_beta_chi_no_observed_odd":True}
    if args.out:
        need(args.out.suffix==".json" and not args.out.exists(),"new JSON only")
        with args.out.open("x",encoding="utf8") as fh:
            json.dump(report,fh,indent=2,allow_nan=False);fh.write("\n")
        print("E28R1_REGISTERED_POSTHOC_FINER_NUMERIC_JSON",args.out,flush=True)
    print("E28R1_POSTHOC_FINER_GRID_COMPLETED_WITH_HONEST_CONTRAST_STATUS",flush=True)

if __name__=="__main__":main()
