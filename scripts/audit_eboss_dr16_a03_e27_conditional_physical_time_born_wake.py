#!/usr/bin/env python3
"""E27 physical-time, EXTERNALLY FORCED, conditional point-halo Born wake.
Not self-consistent Einstein-Vlasov; not LRG/ELG halo drag or an eBOSS odd.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
import audit_eboss_dr16_a03_e17d0_external_potential_vlasov_source as e17
PRO=ROOT/"source_data/eboss_dr16_a03_e27_conditional_physical_time_external_halo_wake_protocol_2026-09-29.json"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E9=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E9"
B1=ROOT/"source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28"
PINS={"FD":"3d21ccdce4d3547be927a1fa76626d69e861109d",
      "plus":"e4765dbab6d9814e999ff2178e020fb7df547c9c",
      "minus":"106acd41328482dde887c1bda0976f34ad44902e"}
B1_SHA=("10fefba505916ecb173352f17da21929b118978a",
        "72707304484a24c361bcff09ad5061bb82dca434",
        "f0a478caf2c8385ff845d0aa3d8cc9353a2a41d2",
        "12e7ae9401dd3e9ea2e95a96d5daa680a00cb3d4")
S=("FD","plus","minus")
COL={"FD":"F0_CLASS_normalized","plus":"Fplus_CLASS_normalized",
     "minus":"Fminus_CLASS_normalized"}
C=299792.458               # km/s
G=4.30091e-9              # Mpc (km/s)^2 / physical solar mass
TCMB=2.7255                # original CLASS code K
KB=8.617333262e-5         # original source eV/K
TNU0=.71611*TCMB*KB        # original CLASS T_ncdm
MNU=.06                   # original neutrino mass eV
KCOM=.05*.6736            # original frozen short k: 1/comoving Mpc
ZOBS=.95; ZINIT=1.0; VHALO=200.0
NCOARSE=129; NFINE=257

def need(ok,why):
    if not ok: raise ValueError("E27_CONDITIONAL_PHYSICAL_TIME_STOP: "+why)

def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def read_pinned(path,sha):
    need(path.is_file() and not path.is_symlink(),"missing "+str(path))
    raw=path.read_bytes()
    need(git_blob(raw)==sha,"changed immutable source "+str(path))
    return json.loads(raw)

def read_sources():
    p=read_pinned(PRO,"d9c9a71841d00cd04db1016ceb2ce1adb5cfc7d5")
    need(p["stage"]=="E27_FIRST_PHYSICAL_TIME_CONDITIONAL_EXTERNAL_HALO_BORN_WAKE"
         and all(p["physical_nonclaims"].values()) and
         p["immutable_numerical_inputs"]["z_initial_zero_response"]==ZINIT
         and p["immutable_numerical_inputs"]["z_observation"]==ZOBS
         and p["immutable_numerical_inputs"]["k_comoving_h_Mpc"]==.05,
         "prospective protocol or physical nonclaim altered")
    need(git_blob(Path(e17.__file__).read_bytes())==
         p["source_pins"]["E17D0_impulse_code"],"original source kernel changed")
    raw=E8.read_bytes()
    need(git_blob(raw)==p["source_pins"]["E8_frozen_4000q_git_blob"] and
         hashlib.sha256(raw).hexdigest()==p["source_pins"]["E8_frozen_4000q_SHA256"],
         "original E8 4000q state source changed")
    with E8.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    need(len(rows)==4000,"original 4000q source cohort")
    q=np.asarray([float(z["q_dimensionless"]) for z in rows])
    need(q[0]==0 and q[-1]==20 and np.all(np.diff(q)>0),
         "original q interval changed")
    prep={s:e17.prepare(q,np.asarray([float(z[COL[s]]) for z in rows]))
          for s in S}
    h={}
    for s in S:
        j=read_pinned(E9/("class_R16_wind_"+s+".json"),PINS[s])
        cp=j["CLASS_params_from_original_wake_two_tracer_fisher"]
        need(j["state"]==s and cp["m_ncdm"]==MNU and
             cp["T_ncdm"]==.71611 and cp["H0"]==67.36 and
             j["observed_odd_data_vector_read"] is False and
             j["CLASS_state_galaxy_xi_not_computed"] is True,
             "original CLASS E9 parameters or observed seal")
        zz=np.asarray([v["z"] for v in j["z_records"]],dtype=float)
        HH=np.asarray([v["H_CLASS_inverse_Mpc"] for v in j["z_records"]])
        need(np.all(np.diff(zz)>0) and zz[0]<=ZOBS and zz[-1]>=ZINIT
             and np.all(HH>0),"original E9 H(z) support")
        h[s]=(zz,HH)
    halos=[]
    for i,sha in enumerate(B1_SHA):
        j=read_pinned(B1/("e17d2b1_original_conditional_wechsler_exterior_0%d.json"%i),sha)
        need(j["redshift_anchor_z0"]==ZOBS and
             j["Wechsler_form_alpha_equals_two_ac_math_QA_NOT_fitted_to_halo"]
                ==(.4 if i%2==0 else .8) and
             j["halo_M200c_use_of_Wechsler_virial_M_history_UNCALIBRATED_ansatz"],
             "original CONDITIONAL halo reference altered")
        halos.append((float(j["original_DM14_anchor_mass_physical_Msun"]),
                      float(j["Wechsler_form_alpha_equals_two_ac_math_QA_NOT_fitted_to_halo"])))
    need(halos[2][0]/halos[0][0]>9.999999999 and
         halos[2][0]/halos[0][0]<10.000000001,
         "original two physical halo mass references changed")
    return p,prep,h,halos

def physical_grid(zh,zH,n):
    z=np.linspace(ZOBS,ZINIT,n)
    H=np.interp(z,zh,zH)
    need(np.all(H>0),"bad H interpolation")
    inv=1/H
    tau=np.zeros(n); fs=np.zeros(n)
    dz=np.diff(z)
    tau[1:]=np.cumsum(dz*(inv[1:]+inv[:-1])/2)
    fs[1:]=np.cumsum(dz*((1+z[1:])*inv[1:]+
                             (1+z[:-1])*inv[:-1])/2)
    phase=KCOM*(TNU0/MNU)*fs
    lookback=VHALO/C*tau  # comoving Mpc, constant physical wind
    return z,H,phase,lookback,tau[-1]

def response(prep,background,halo,n,mu=1.,wind=VHALO):
    M0,alpha=halo
    z,H,s,dx,tau=physical_grid(*background,n)
    if wind!=VHALO: dx*=wind/VHALO
    H_F=np.asarray([e17.source(prep,float(x))[1] for x in s])
    weight=np.exp(-alpha*(z-ZOBS)/(1+ZOBS))*H_F/H
    phase=np.exp(1j*KCOM*mu*dx)
    integral=np.trapz(weight*phase,z)
    pref=-4*math.pi*G/(C*C)*(MNU/TNU0)*(M0/KCOM)
    return complex(pref*integral),tau,phase[-1]

def evaluate(x):
    _,prep,h,halos=x
    result={}
    for s in S:
        runs={}
        for i,halo in enumerate(halos):
            a,tau,ph=response(prep[s],h[s],halo,NFINE)
            c,_,_=response(prep[s],h[s],halo,NCOARSE)
            need(a.real>0 and a.imag>0 and
                 abs(a-c)/abs(a)<5e-3 and
                 abs(a.imag-c.imag)/abs(a.imag)<5e-3,
                 "physical-time integral or independent grid convergence "+s+str(i))
            neg,_,_=response(prep[s],h[s],halo,NFINE,mu=-1)
            zero,_,_=response(prep[s],h[s],halo,NFINE,wind=0)
            transverse,_,_=response(prep[s],h[s],halo,NFINE,mu=0)
            need(abs(a-neg.conjugate())<1e-12*abs(a) and
                 zero.imag==0 and transverse.imag==0 and
                 abs(a.real-transverse.real)<1e-4*abs(a.real),
                 "wind/LOS symmetry or transverse phase "+s+str(i))
            runs[str(i)]={"anchor_M0_physical_Msun":halo[0],"alpha":halo[1],
                          "Re_delta_nu_over_nbar_k_Mpc3":a.real,
                          "Im_delta_nu_over_nbar_k_Mpc3":a.imag,
                          "redshift_quadrature_complex_relative_gap":abs(a-c)/abs(a),
                          "redshift_quadrature_imag_relative_gap":
                              abs(a.imag-c.imag)/abs(a.imag),
                          "lookback_conformal_Mpc":tau,
                          "end_halo_displacement_Fourier_phase_rad":float(np.angle(ph))}
        for i in (0,1):
            lo=runs[str(i)];hi=runs[str(i+2)]
            need(math.isclose(hi["Re_delta_nu_over_nbar_k_Mpc3"]/
                              lo["Re_delta_nu_over_nbar_k_Mpc3"],10,rel_tol=1e-12)
                 and math.isclose(hi["Im_delta_nu_over_nbar_k_Mpc3"]/
                                  lo["Im_delta_nu_over_nbar_k_Mpc3"],10,rel_tol=1e-12),
                 "point-mass Born halo mass scale closure")
        need(runs["0"]["Re_delta_nu_over_nbar_k_Mpc3"]>
             runs["1"]["Re_delta_nu_over_nbar_k_Mpc3"],
             "earlier-alpha conditional mass source should yield less response")
        result[s]=runs
    return result

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path)
    args=a.parse_args()
    x=read_sources()
    res=evaluate(x)
    print("E27_PHYSICAL_TIME_CAUSAL_EXTERNAL_HALO_BORN_3F_4H_PASS",flush=True)
    for s in S:
        for key in ("0","1","2","3"):
            w=res[s][key]
            print("E27_PHYSICAL_BORN_SOURCE",s,"HALO_CASE",key,
                  "RE_MPC3",format(w["Re_delta_nu_over_nbar_k_Mpc3"],".15g"),
                  "IM_MPC3",format(w["Im_delta_nu_over_nbar_k_Mpc3"],".15g"),
                  "D_Z_REL",format(w["redshift_quadrature_imag_relative_gap"],".5g"),
                  flush=True)
    out={"stage":"E27_FIRST_PHYSICAL_TIME_CONDITIONAL_EXTERNAL_HALO_BORN_WAKE",
         "status":"CONDITIONAL_PHYSICAL_TIME_POINT_HALO_BORN_WAKE_MODE_PASS_NOT_GALAXY",
         "protocol_git_blob":"d9c9a71841d00cd04db1016ceb2ce1adb5cfc7d5",
         "z_obs":ZOBS,"z_initial_zero_wake":ZINIT,"k_comoving_h_Mpc":.05,
         "reference_halo_v_relative_neutrinos_km_s":VHALO,
         "external_halo_conditioning":"Original E17D2B1 alpha=0.4/0.8 and physical M0, point-mass one-k",
         "original_F_results":res,"source_only_previous_E8_E9_E17D0_B1_untouched":True,
         "not_full_primordial_halo_wake":True,"not_halo_drag_or_tracer_beta":True,
         "no_galaxy_xi_or_observed_odd_or_new_CLASS_ABACUS_WSL":True}
    if args.out:
        need(args.out.suffix==".json" and not args.out.exists(),"output collision")
        with args.out.open("x") as f:
            json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
        print("E27_NEW_CONDITIONAL_NUMERIC_RESULT",args.out,flush=True)
if __name__=="__main__":main()
