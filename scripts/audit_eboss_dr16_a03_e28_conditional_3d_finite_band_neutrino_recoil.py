#!/usr/bin/env python3
"""E28: finite-comoving-k-band gravitational recoil from E27 conditional wake.

This is a FIRST-ORDER externally forced short-history Born MODEL, not total
halo drag, a measured M200c history, a galaxy response or an eBOSS odd.
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
import audit_eboss_dr16_a03_e27_conditional_physical_time_born_wake as old

PRO=ROOT/"source_data/eboss_dr16_a03_e28_conditional_finite_band_3d_wake_recoil_protocol_2026-09-29.json"
PINS={
 "protocol":"e8feedc4d09965a24d8eb60c16454d6e03b2fbff",
 "old_runner":"851df1faf5897464665ba79f4534c4b6de234eca",
 "old_report":"87380ff6bdcb7cda3e244b0b4a716c84df0c4034",
 "low_halo":"776322703b3f93071e7c53f2cb803a667048c870",
 "high_halo":"f6fc9e34b0382e6f8d596a2a5c3721adc9dcd5f8"}
REPORT=ROOT/"source_data/eboss_dr16_a03_e27_archived_conditional_physical_time_born_report_2026-09-29.json"
HALOS=ROOT/"source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28"
HFILES=("e17d2b0_original_DM14_conditional_snapshot_04.json",
        "e17d2b0_original_DM14_conditional_snapshot_13.json")
MPC_M=3.0856775814913673e22
HBARC_EVM=1.973269804e-7
EV_J=1.602176634e-19
C_MS=299792458.
MSUN_KG=1.98847e30
BANDS=(.01,.1,1.,8.)  # 1/comoving Mpc, analysis bands ONLY, not survey cuts
STATES=("FD","plus","minus")
NFW_N=64

def need(ok,msg):
    if not ok:raise ValueError("E28_PHYSICAL_MODEL_FAIL_CLOSED: "+msg)

def blob(path,sha):
    need(path.is_file() and not path.is_symlink(),"missing original "+str(path))
    raw=path.read_bytes()
    b=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    need(b==sha,"original source Git blob changed "+str(path))
    return raw

def sources():
    p=json.loads(blob(PRO,PINS["protocol"]))
    need(p["stage"]=="E28_CONDITIONAL_FINITE_BAND_3D_NEUTRINO_WAKE_HALO_RECOIL"
         and p["frozen_physical_benchmark"]["z_obs"]==old.ZOBS
         and p["frozen_physical_benchmark"]["z_init"]==old.ZINIT
         and p["frozen_physical_benchmark"]["k_comoving_Mpc_inverse_band"]==list(BANDS)
         and p["protocol_clarification_before_any_E28_execution"].startswith(
             "First protocol text mistakenly said physical r_s fixed")
         and all(p["limits"].values()),"E28 preregistration/science boundary")
    blob(Path(old.__file__),PINS["old_runner"])
    old_report=json.loads(blob(REPORT,PINS["old_report"]))
    need(old_report["original_Fpm_halos_independently_calibrated"] is False
         and old_report["full_neutrino_halo_drag_or_3D_k_integral_calculated"] is False
         and old_report["observed_odd_SEALED"],"E27 original scientific STOP")
    base=old.read_sources()  # Pins original 4000q, E9 CLASS H, E17D0, E17D2B1
    halodata=[]
    for file,key in zip(HFILES,("low_halo","high_halo")):
        h=json.loads(blob(HALOS/file,PINS[key]))["published_population_halo_snapshot"]
        need(h["z_original_CLASS_math_node_only"]==old.ZOBS
             and h["illustrative_log10_concentration_offset_dex_NOT_z095_posterior"]==0
             and h["c200c_population_median_DM14"]>0
             and h["r_s_comoving_Mpc_DM14"]>0,"original conditional NFW mass/profile")
        halodata.append(h)
    for i,h in enumerate(halodata):
        need(math.isclose(h["M200c_physical_Msun_DM14"],base[3][2*i][0],
                          rel_tol=1e-13),"E27 halo mass anchor vs DM14 snapshot")
    return p,base,halodata

def q_mass_density(prep,z):
    # Original F(q) already includes the two nu/anti spin states, in
    # 2/(2pi)^3 convention. This is nonrelativistic REST mass density only.
    n_per_m3=(4*math.pi*(old.TNU0*(1+z)/HBARC_EVM)**3)*prep[4]
    rho=n_per_m3*(old.MNU*EV_J/C_MS**2)*MPC_M**3/MSUN_KG
    need(0<rho<1e11,"F-state rest-mass density physical unit")
    return rho

def k_grid(nseg):
    need(nseg in (12,24),"unregistered k integration grid")
    pieces=[np.geomspace(BANDS[i],BANDS[i+1],nseg+1)
            for i in range(3)]
    k=np.concatenate([pieces[0],pieces[1][1:],pieces[2][1:]])
    need(len(k)==3*nseg+1 and np.all(np.diff(k)>0),"fixed k band")
    return k

def simpson_bands(k,values,nseg):
    result=[]
    for ib in range(3):
        start=ib*nseg
        q=k[start:start+nseg+1]
        f=values[start:start+nseg+1]
        dlog=np.log(q[-1]/q[0])/nseg
        integral=(dlog/3)*(f[0]+f[-1]+
                 4*np.sum(f[1:-1:2])+2*np.sum(f[2:-1:2]))
        result.append(float(integral))
    return np.cumsum(result)

def nfw_u(k,rscom,c,nodes,weights):
    # Exact normalized spherical truncated NFW density transform,
    # integrated independently of original E17D2a Si/Ci implementation.
    y=(nodes+1)*c/2
    w=(weights*c/2)*y/(1+y)**2
    A=math.log1p(c)-c/(1+c)
    u=np.sum(w[None,:]*np.sinc(k[:,None]*rscom*y[None,:]/math.pi),
             axis=1)/A
    need(np.all(np.isfinite(u)) and np.all(u>0) and
         max(abs(float(u[0])-1),abs(float(u[-1])-1))<.9,
         "truncated NFW radial transform unphysical")
    return u

def standalone_E27_mode_check(prep,background,halo):
    # Same unfiltered U=1 one-k, mu=+1 E27 parent, independent of k-band.
    z,H,s,dx,_=old.physical_grid(*background,257)
    hf=np.asarray([old.e17.source(prep,float(v))[1] for v in s])
    pref=-4*math.pi*old.G/(old.C**2)*(old.MNU/old.TNU0)*(halo[0]/old.KCOM)
    here=pref*np.trapz(np.exp(-halo[1]*(z-old.ZOBS)/(1+old.ZOBS))*
                       hf/H*np.exp(1j*old.KCOM*dx),z)
    ref=old.response(prep,background,halo,257,mu=1)[0]
    need(abs(here-ref)<1e-12*max(1,abs(ref)),"E27 original single mode drift")
    return abs(here-ref)

def kernels_for_state(prep,background,k,nz):
    z,H,s,dx,_=old.physical_grid(*background,nz)
    # E27 s was computed at original short K; here only change k itself.
    phases=k[:,None]*s[None,:]/old.KCOM
    hf=np.empty_like(phases)
    for ik in range(len(k)):
        for iz in range(nz):
            hf[ik,iz]=old.e17.source(prep,float(phases[ik,iz]))[1]
    ang=2*old.e17.jp(k[:,None]*dx[None,:])
    need(np.all(np.isfinite(hf)) and np.all(np.isfinite(ang))
         and np.max(abs(ang+2*old.e17.jp(-k[:,None]*dx[None,:])))<2e-15,
         "exact angular conjugation/reversal")
    return z,H,hf,ang,dx

def evaluate(prep,background,halo,k,nz,nseg,rho,nodes,weights):
    M0,alpha,c,rscom=halo
    profile=nfw_u(k,rscom,c,nodes,weights)
    z,H,hf,ang,dx=kernels_for_state(prep,background,k,nz)
    growth=np.exp(-alpha*(z-old.ZOBS)/(1+old.ZOBS))
    integral=np.trapz(hf*growth[None,:]*ang/H[None,:],z,axis=1)
    pref=-4*math.pi*old.G/(old.C**2)*(old.MNU/old.TNU0)*(M0/k)
    # Two profile factors: forcing density U_h and mass-averaged recoil U_h.
    angular_integrated=pref*integral*profile**2
    # ∫dk k * angular integrated scalar = ∫dlnk k² * scalar.
    per_log_k=k*k*angular_integrated
    three=simpson_bands(k,per_log_k,nseg)
    a=1/(1+old.ZOBS)
    acceleration=old.G*a*rho/math.pi*three
    need(np.all(np.isfinite(acceleration)) and
         np.all(acceleration<0),"unexpected sign / divergence of conditional drag")
    return {"accel_prefix":acceleration,"u_end":float(profile[-1]),
            "u_original_short":float(nfw_u(np.array([old.KCOM]),rscom,c,
                                                  nodes,weights)[0])}

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--out",type=Path)
    args=arg.parse_args()
    p,base,profiles=sources()
    _,preps,backgrounds,original_halos=base
    nodes,weights=np.polynomial.legendre.leggauss(NFW_N)
    k_coarse=k_grid(12);k_fine=k_grid(24)
    need(all(math.isclose(k_fine[2*i],k_coarse[i],rel_tol=1e-14)
             for i in range(len(k_coarse))),"non-nested independent k QA")
    report={}
    for state in STATES:
        rho=q_mass_density(preps[state],old.ZOBS)
        refgap=standalone_E27_mode_check(preps[state],backgrounds[state],
                                        original_halos[0])
        cases={}
        # Cache response kernels: do NOT re-compute 4000q F for every halo mass.
        cached={}
        for nseg,nz,k in ((12,129,k_coarse),(24,257,k_fine)):
            cached[nseg]=kernels_for_state(preps[state],backgrounds[state],k,nz)
        for i,(M0,alpha) in enumerate(original_halos):
            prof=profiles[i//2]
            c=prof["c200c_population_median_DM14"]
            rs=prof["r_s_comoving_Mpc_DM14"]
            grid_results={}
            for nseg,nz,k in ((12,129,k_coarse),(24,257,k_fine)):
                z,H,hf,ang,dx=cached[nseg]
                u=nfw_u(k,rs,c,nodes,weights)
                growth=np.exp(-alpha*(z-old.ZOBS)/(1+old.ZOBS))
                time=np.trapz(hf*growth[None,:]*ang/H[None,:],z,axis=1)
                pref=-4*math.pi*old.G/(old.C**2)*(old.MNU/old.TNU0)*(M0/k)
                log_integrand=k*k*pref*time*u*u
                prefix=simpson_bands(k,log_integrand,nseg)
                accel=(old.G/(1+old.ZOBS)*rho/math.pi)*prefix
                need(np.all(np.isfinite(accel)) and np.all(accel<0),
                     "finite-band conditional drag lost odd sign")
                grid_results[nseg]={"band_prefix_accel_km2_s2_per_Mpc":list(map(float,accel)),
                                    "u_at_kmax":float(u[-1])}
            fine=np.asarray(grid_results[24]["band_prefix_accel_km2_s2_per_Mpc"])
            coarse=np.asarray(grid_results[12]["band_prefix_accel_km2_s2_per_Mpc"])
            gap=np.abs(fine-coarse)/np.maximum(np.abs(fine),1e-80)
            need(gap[-1]<p["acceptance"]["FD_coarse_fine_full_band_relative_gap_at_most"],
                 "combined time and log k convergence not sufficient")
            cases[str(i)]={"M0_Msun":M0,"alpha":alpha,"r_s_comoving_Mpc":rs,
                           "c200c_conditional_DM14":c,
                           "rho_nu_rest_Msun_per_Mpc3":rho,
                           "u_kmax":grid_results[24]["u_at_kmax"],
                           "prefix_end_k_comoving_Mpc_inv":list(BANDS[1:]),
                           "coarse_prefix_accel_km2_s2_per_Mpc":list(map(float,coarse)),
                           "fine_prefix_accel_km2_s2_per_Mpc":list(map(float,fine)),
                           "coarse_vs_fine_prefix_relative_gap":list(map(float,gap))}
        report[state]={"rho_nu_rest_Msun_per_Mpc3":rho,
                       "E27_original_one_k_257_node_replay_absolute_gap":refgap,
                       "cases":cases}
        print("E28_NUMERIC_STATE",state,"RHO_REST",format(rho,".12g"),
              "MLOW_ALPHA04_ACCEL_BAND_0p01_TO_8",
              format(cases["0"]["fine_prefix_accel_km2_s2_per_Mpc"][-1],".13g"),
              flush=True)
    contrast={}
    for i in range(4):
        fd=report["FD"]["cases"][str(i)]
        plus=report["plus"]["cases"][str(i)]
        minus=report["minus"]["cases"][str(i)]
        fine=minus["fine_prefix_accel_km2_s2_per_Mpc"][-1]-plus["fine_prefix_accel_km2_s2_per_Mpc"][-1]
        coarse=minus["coarse_prefix_accel_km2_s2_per_Mpc"][-1]-plus["coarse_prefix_accel_km2_s2_per_Mpc"][-1]
        gap=abs(fine-coarse)/max(abs(fine),1e-80)
        resolved=(fine!=0 and (fine>0)==(coarse>0) and gap<.01)
        contrast[str(i)]={"Fminus_minus_Fplus_km2_s2_Mpc_inv":fine,
                          "relative_to_FD_band_recoil":fine/fd["fine_prefix_accel_km2_s2_per_Mpc"][-1],
                          "combined_coarse_fine_signed_difference_relative_gap":gap,
                          "resolved_to_registered_one_percent_of_difference":resolved}
        print("E28_SIGNED_CONTRAST",i,"REL_FD",
              format(contrast[str(i)]["relative_to_FD_band_recoil"],".12g"),
              "DIFFERENCE_RESOLVED",resolved,"GAP",format(gap,".6g"),flush=True)
    # Exact analytic angular nulls on the frozen k domain: no odd force at
    # zero relative velocity, opposite sign for reversed halo wind.
    test=np.array([0.,.003,.5,1.3])
    need(np.all(old.e17.jp(np.zeros_like(test))==0) and
         np.max(abs(old.e17.jp(test)+old.e17.jp(-test)))<1e-15,
         "wind reversal or no-wind recoil null")
    print("E28_3D_ANGULAR_SINC_DERIVATIVE_RECOIL_SIGN_AND_WIND_NULL_PASS",flush=True)
    out={"stage":p["stage"],"status":"CONDITIONAL_FINITE_K_BAND_3D_RECOIL_NOT_TOTAL_HALO_DRAG",
         "protocol_git_blob":PINS["protocol"],"frozen_families":report,
         "signed_Fpm_contrast":contrast,
         "profile":"DM14 median truncated NFW COMOVING SHAPE HELD FIXED z1 to z.95; normalization original E17D2B1 M_alpha",
         "rho_nu":"one massive nu/anti species nonrel rest mass from original F CLASS normalization",
         "force_definition":"halo-averaged direct gravity from neutrino wake Psi_nu included in Psi_total; not independent E25 epsilon",
         "units":"a_parallel=(km/s)^2/comoving Mpc converted from gravitational potential derivative to physical acceleration",
         "numerical_k_bands_only_not_all_k_total_force":list(BANDS),
         "no_pre_z1_neutrino_wake":True,
         "E27_external_source_Born_and_DM14_CLASS_cosmology_mismatch":True,
         "no_physical_highz_LRG_ELG_beta_chi_or_observed_odd":True,
         "observed_odd_SEALED_main_untouched_PR_draft":True}
    if args.out:
        need(args.out.suffix==".json" and not args.out.exists(),"new JSON only")
        with args.out.open("x",encoding="utf8") as fh:
            json.dump(out,fh,indent=2,allow_nan=False);fh.write("\n")
        print("E28_ORIGINAL_CONDITIONAL_3D_RECOIL_JSON",args.out,flush=True)
    print("E28_CONDITIONAL_FINITE_BAND_NUMERICAL_RECOIL_PASS_NOT_TOTAL_FORCE",flush=True)

if __name__=="__main__":main()
