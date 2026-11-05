#!/usr/bin/env python3
"""A03E8b: SI-normalized *conditional static halo* neutrino-wake Green kernel.

Published Okoli et al., MNRAS 468 (2017) Eq.(19)-(20), physical k, fixed
potential and parallel relative speed. Ratio of isotropic matched F+/- to the
published Fermi-Dirac resonance is determined by the frozen CLASS-normalized
kinetic distribution. NO evolving CLASS solution or absolute eBOSS galaxy xi.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import build_eboss_dr16_a03_e8_frozen_resonant_occupation as e8

PROTOCOL=ROOT/"source_data/eboss_dr16_a03_e8b_static_halo_green_kernel_protocol_2026-09-27.json"
P_BLOB="7865f059d312a0a97753f7c13b6faa5bb07a0e16"
ORIGINAL_E8_SOURCE_BLOB="3db381a2942832b6500a60291ac864eff929b200"
OUT_DEFAULT=ROOT/"eboss_workspace/a03_physics_source/e8_frozen_kinetic_resonance/conditional_static_halo_green_kernel.json"

def require(condition,message):
    if not condition:raise ValueError(message)

def protocol():
    require(e8.git_blob(PROTOCOL.read_bytes())==P_BLOB,
            "E8b source-only SI physical protocol Git blob changed")
    require(e8.git_blob(Path(e8.__file__).read_bytes())==ORIGINAL_E8_SOURCE_BLOB,
            "E8a frozen source code changed after SI physical protocol")
    p=json.loads(PROTOCOL.read_bytes())
    require(p["parent_frozen_E8_source_code"]["git_blob_sha1"]==ORIGINAL_E8_SOURCE_BLOB and
            p["original_frozen_source_csv_sha256"]==
            "bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0" and
            p["guards"]["reinterpret_as_eBOSS_galaxy_xi"] is False and
            p["guards"]["observed_odd_read"] is False and
            p["physical_reference"]["physical_not_comoving_k"] is True and
            p["physical_reference"]["one_ncdm_species"]==1,
            "E8b physical SI source scope changed")
    e8.source_gate()
    ref=p["physical_reference"]
    require(math.isclose(float(e8.cro.H0)/100.,float(ref["h"]),rel_tol=0,abs_tol=1e-14) and
            float(ref["mass_eV"])==e8.MASS_EV and
            float(ref["reference_z"])==.95 and
            float(ref["reference_abs_v_parallel_kms"])==200,
            "Historical original CLASS H0 or source physical constants changed")
    return p

def conditional_Green(src,p,*,z,v_parallel_kms,k_physical_h_per_Mpc=.05):
    """Per-Phi Fourier gravitational-field coefficient along +k, in SI m^-1.

    Eq.(20) k and Phi are from static PHYSICAL-coordinates halo response.
    Do not feed k_com from E6 without conversion k_phys=(1+z) k_com.
    The physical neutrino occupancy at resonance enters exactly as
    F_FD(q_res)=1/(exp(q_res)+1), and F+/- rescale it by the same-source ratio.
    """
    ref=p["physical_reference"]
    require(np.isfinite(k_physical_h_per_Mpc) and k_physical_h_per_Mpc>0 and
            np.isfinite(v_parallel_kms) and np.isfinite(z) and z>=0,
            "Nonphysical SI Green-function inputs")
    k=float(k_physical_h_per_Mpc)*float(ref["h"])/float(ref["Mpc_m"])
    v=float(v_parallel_kms)*1000.
    m=float(ref["mass_eV"])*float(ref["rest_mass_conversion_eVc2_to_kg"])
    G=float(ref["Newton_G_m3_kg_s2"])
    hb=float(ref["hbar_Js"])
    speed=float(ref["c_m_per_s"])
    require(k>0 and m>0 and hb>0 and G>0 and speed>0,
            "Invalid dimensional constants or physical k")
    need=e8.conditional_resonant_ratio(src["q"],src["f0"],src["fp"],src["fm"],
              e8.MASS_EV,z,abs(v_parallel_kms))
    qr=float(need["q_res"])
    fd=1./(1.+math.exp(qr))
    K=2.*G*m**4*v/(hb**3*k)*fd # m^-1 signed along k at fixed Phi_k
    kp=K*float(need["f_plus_over_FD"])
    km=K*float(need["f_minus_over_FD"])
    require(all(np.isfinite([K,kp,km])) and
            abs((kp+km)/2.-K)<=1e-12*max(abs(K),1e-300) and
            abs((kp-km)-K*float(need["f_plus_minus_f_minus_over_FD"]))<=
            1e-12*max(abs(K),1e-300),
            "Conditional absolute Green source loses original FD-relative closure")
    return {
        "z_example_not_measured_eBOSS_zeff":float(z),
        "v_parallel_kms_example_not_measured_eBOSS":float(v_parallel_kms),
        "k_PHYSICAL_h_per_Mpc_not_comoving":float(k_physical_h_per_Mpc),
        "k_PHYSICAL_SI_per_m":k,
        "m_nu_kg":m,
        "q_res":qr,
        "FD_occupation_without_CLASS_phase_space_prefactor":fd,
        "signed_gnu_over_Phi_FD_per_m":K,
        "signed_gnu_over_Phi_plus_per_m":kp,
        "signed_gnu_over_Phi_minus_per_m":km,
        "signed_gnu_over_Phi_plus_minus_minus_per_m":kp-km,
        "exact_pair_average_equals_FD":True,
        "no_galaxy_xi_or_survey_detection":True,
    }

def self_test():
    p=protocol()
    src=e8.kinetic_source()
    x=conditional_Green(src,p,z=.95,v_parallel_kms=200.,
                        k_physical_h_per_Mpc=.05)
    minus=conditional_Green(src,p,z=.95,v_parallel_kms=-200.,
                            k_physical_h_per_Mpc=.05)
    zero=conditional_Green(src,p,z=.95,v_parallel_kms=0.,
                           k_physical_h_per_Mpc=.05)
    require(x["k_PHYSICAL_SI_per_m"]>0 and x["signed_gnu_over_Phi_FD_per_m"]>0 and
            abs(x["signed_gnu_over_Phi_FD_per_m"]+
                minus["signed_gnu_over_Phi_FD_per_m"])<
                abs(x["signed_gnu_over_Phi_FD_per_m"])*1e-13 and
            abs(x["signed_gnu_over_Phi_plus_per_m"]+
                minus["signed_gnu_over_Phi_plus_per_m"])<
                abs(x["signed_gnu_over_Phi_plus_per_m"])*1e-13 and
            zero["signed_gnu_over_Phi_FD_per_m"]==0. and
            zero["signed_gnu_over_Phi_plus_per_m"]==0.,
            "Physical signed velocity/zero-force controls failed")
    for k in (0.,-1.):
        try:conditional_Green(src,p,z=.95,v_parallel_kms=200.,
                              k_physical_h_per_Mpc=k)
        except ValueError:pass
        else:raise AssertionError("Invalid SI physical k accepted")
    require(x["signed_gnu_over_Phi_plus_per_m"]<
            x["signed_gnu_over_Phi_FD_per_m"]<
            x["signed_gnu_over_Phi_minus_per_m"] and
            abs(x["signed_gnu_over_Phi_plus_minus_minus_per_m"]/
                x["signed_gnu_over_Phi_FD_per_m"]-
                (-.3704048092478419))<2e-12,
            "Frozen E8a source and exact SI Eq20 Green coefficient disagree")
    print("E8B_SOURCE_ONLY_EXACT_STATIC_HALO_SI_GREEN_COEFFICIENT_OK",
          "KFD_PER_M",x["signed_gnu_over_Phi_FD_per_m"],
          "KPLUS_PER_M",x["signed_gnu_over_Phi_plus_per_m"],
          "KMINUS_PER_M",x["signed_gnu_over_Phi_minus_per_m"],
          "DIFF_PER_M",x["signed_gnu_over_Phi_plus_minus_minus_per_m"],
          "NO_CLASS NO_FITS NO_OBSERVED",flush=True)
    return p,src

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--self-test",action="store_true")
    a.add_argument("--output",type=Path,default=OUT_DEFAULT)
    args=a.parse_args()
    p,src=self_test()
    if args.self_test:return 0
    response=[conditional_Green(src,p,z=z,v_parallel_kms=v,
                                k_physical_h_per_Mpc=.05)
              for z in (.90,.95,1.0) for v in (100.,200.,400.)]
    summary={
        "status":"CONDITIONAL_STATIC_PHYSICAL_K_NEUTRINO_WAKE_GREEN_SI_SOURCE_READY_GALAXY_XI_UNAVAILABLE",
        "scope":"Published Eq20 per-unit-potential STATIC halo linear collisionless response for one 0.06eV massive neutrino species, same potential/v/k and original CLASS kinetic-source normalization",
        "model_reference":"Okoli et al. MNRAS 468 (2017) Eq.14--20; exact resonant FD occupation, no fitted mu(v)",
        "physical_wave_number_convention":"PHYSICAL k, not E6 comoving k. k_phys=(1+z)*k_com for future transfer to comoving window.",
        "original_E8_source_CSV_sha256":p["original_frozen_source_csv_sha256"],
        "original_E8_source_protocol_git_blob_sha1":p["parent_frozen_E8_protocol"]["git_blob_sha1"],
        "source_E8B_protocol_git_blob_sha1":P_BLOB,
        "SI_unit":"m^-1, signed longitudinal gravitational Fourier field per unit Phi_k; physical static halo Green coefficient",
        "illustrative_cases_not_fit_or_predeclared_science_acceptance":response,
        "eBOSS_LRG_ELG_absolute_xi_ell_computed":False,
        "CLASS_time_dependent_phi_or_halo_velocity_response_computed":False,
        "eBOSS_tracer_bias_and_HOD_calibrated":False,
        "eBOSS_pair_redshift_window_certified":False,
        "eBOSS_independent_inference_covariance_computed":False,
        "observed_odd_data_vector_read":False,
        "new_FITS_or_mock_download":False,
    }
    raw=(json.dumps(summary,indent=2,allow_nan=False)+"\n").encode()
    e8.atomic_reproducible(args.output,raw)
    print("E8B_CONDITIONAL_SI_STATIC_HALO_GREEN_SOURCE_ONLY_SAVED",
          args.output,"SHA256",e8.sha(raw),flush=True)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
