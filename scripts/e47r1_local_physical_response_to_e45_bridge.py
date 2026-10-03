#!/usr/bin/env python3
"""
E47 local source-to-eBOSS response bridge.

Purpose:
  Derive an absolute dimensionless eBOSS odd-vector *basis* for the frozen
  Fminus-Fplus source difference under the restricted E25 linear local response
  model, parameterized by

      g0 = b_L d_E - b_E d_L
      g2 = f (d_E - d_L)

  with d_a dimensionless in the E25 convention. No values of g0/g2 are fitted
  or assumed. The script returns A_E45 = alpha*g0 + beta*g2 per mock/cap.

Data scope:
  - exact archived E8 4000q F+/F- CSV
  - CLASS source-only transfer calculation
  - exact small E0/E1 random-window NPZ (250104 bytes)
  - archived E3 JSON for the synthetic E45 matched-projection q vectors
  - NO FITS, NO ASDF, NO observed galaxy rows, NO observed odd vector
"""
from pathlib import Path
import io, os, json, math, hashlib, zipfile, statistics
import numpy as np
from scipy.special import spherical_jn
from classy import Class

ROOT=Path(__file__).resolve().parents[1] if Path(__file__).resolve().parent.name=="scripts" else Path.cwd()
SRC=ROOT/"source_data"

E8CSV=SRC/"eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E3=SRC/"eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json"
E45=SRC/"e45_archive_only_synthetic_odd_amplitude_recovery.json"
OUT=SRC/"e47_local_physical_response_basis_to_e45.json"

NPZ_NAME="ninemock_window_ensemble.npz"
NPZ_SHA="4f63e04ce4c5cfdb9404e5a56dee44907cb9a28cc015c6ed27ac718195feae6b"
NPZ_BYTES=250104

E8_CSV_SHA="bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0"
C_KMS=299792.458
Z=0.95
MASS=.06
KCHECK=np.array([.001,.002,.003,.005],float)
E21_PLUS=np.array([4328422.243894628,11606037.777232643,19862015.365340475,33329438.921102364],float)
E21_MINUS=np.array([4297706.35960601,11550362.014643177,19824645.207705196,33268689.17565131],float)

def sha(b): return hashlib.sha256(b).hexdigest()

def require(c,msg):
    if not c: raise RuntimeError(msg)

def find_npz_bytes():
    # Exact extracted member first.
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in {".git",".venv","venv","__pycache__"}]
        if NPZ_NAME in filenames:
            p=Path(dirpath)/NPZ_NAME
            b=p.read_bytes()
            if len(b)==NPZ_BYTES and sha(b)==NPZ_SHA:
                return b, str(p), "direct_npz"
    # Then small/medium local zip archives only.
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in {".git",".venv","venv","__pycache__"}]
        for fn in filenames:
            if not fn.lower().endswith(".zip"): continue
            p=Path(dirpath)/fn
            try:
                if p.stat().st_size > 200_000_000: continue
                with zipfile.ZipFile(p) as z:
                    names=z.namelist()
                    candidates=[n for n in names if n.endswith("/"+NPZ_NAME) or n==NPZ_NAME]
                    for n in candidates:
                        b=z.read(n)
                        if len(b)==NPZ_BYTES and sha(b)==NPZ_SHA:
                            return b, f"{p}::{n}", "zip_member"
            except (OSError, zipfile.BadZipFile):
                pass
    raise RuntimeError(
        "MISSING exact ninemock_window_ensemble.npz (250104 B, SHA "+NPZ_SHA+"). "
        "Do not substitute another operator."
    )

def load_e8():
    require(E8CSV.exists(),"MISSING "+str(E8CSV))
    raw=E8CSV.read_bytes()
    require(sha(raw)==E8_CSV_SHA,"Frozen E8 CSV SHA changed")
    a=np.genfromtxt(E8CSV,delimiter=",",names=True)
    q=np.asarray(a["q_dimensionless"],float)
    fp=np.asarray(a["Fplus_CLASS_normalized"],float)
    fm=np.asarray(a["Fminus_CLASS_normalized"],float)
    require(len(q)==4000 and np.min(fp)>0 and np.min(fm)>0,"Invalid frozen E8 arrays")
    return q,fp,fm

def write_psd(path,q,f):
    txt="\n".join(f"{qi:.14e} {fi:.14e}" for qi,fi in zip(q,f))+"\n"
    path.write_text(txt,encoding="utf-8")

R16=16.0

def tophat(x):
    x=np.asarray(x,float)
    out=np.ones_like(x)
    m=np.abs(x)>1e-5
    xm=x[m]
    out[m]=3.0*(np.sin(xm)-xm*np.cos(xm))/xm**3
    return out

def dense_transfer_source(q,f,label):
    # Reuse the frozen project CLASS parameterization, changing only output to include vTk.
    # IMPORTANT: preserve the original E12/E13 construction order:
    # interpolate *transfer functions first*, then form the nonlinear cross product.
    import sys
    sys.path.insert(0,str(ROOT/"code"))
    import wake_two_tracer_fisher as base
    import class_response_optimize as cro

    work=ROOT/"eboss_workspace/a03_physics_source/e47_dense_linear_response"
    work.mkdir(parents=True,exist_ok=True)
    psd=work/f"frozen_{label}.dat"
    write_psd(psd,q,f)

    base.ZBINS=np.array([.9,1.],dtype=float)
    pars=base.class_params(psd,MASS)
    pars["output"]="mPk,dTk,vTk"
    require(pars["gauge"]=="newtonian","Unexpected CLASS gauge")
    require(float(pars["P_k_max_h/Mpc"])==.25,"Frozen P_k_max changed")

    c=Class(); c.set(pars); c.compute()
    try:
        h=float(c.h())
        t=c.get_transfer(z=Z,output_format="class")
        names=set(t.keys())
        need={"k (h/Mpc)","t_ncdm[0]","t_cdm","d_b","d_cdm"}
        require(need.issubset(names),"Missing direct CLASS transfer columns")
        kin=np.asarray(t["k (h/Mpc)"],float)
        require(np.all(np.diff(kin)>0),"CLASS k grid not increasing")
        theta=np.asarray(t["t_ncdm[0]"],float)-np.asarray(t["t_cdm"],float)
        db=np.asarray(t["d_b"],float)
        dc=np.asarray(t["d_cdm"],float)
        cb=(float(cro.OMEGA_B)*db+float(cro.OMEGA_CDM)*dc)/(float(cro.OMEGA_B)+float(cro.OMEGA_CDM))
        require(np.isfinite(theta).all() and np.isfinite(cb).all(),"Nonfinite CLASS transfer source")
        return {
            "h":h,
            "kin":kin,
            "theta":theta,
            "cb":cb,
            "As":float(cro.A_S),
            "ns":float(base.NS),
        }
    finally:
        c.struct_cleanup(); c.empty()

def filtered_Cv_at(kh,source):
    # E21 physical-unit source is the E13 R16-filtered cross:
    # C_v = W_R16(k) * P_R(k) * delta_cb(k) * c * theta_rel(k) / k.
    # Do NOT interpolate the already-multiplied C_v product; E12/E13 interpolate
    # d_cb and theta_rel individually to the target k before multiplying.
    kh=np.asarray(kh,float)
    kin=source["kin"]
    require(np.min(kh)>=kin[0] and np.max(kh)<=kin[-1],"Requested k outside CLASS transfer grid")
    cb=np.interp(kh,kin,source["cb"])
    theta=np.interp(kh,kin,source["theta"])
    k_mpc=kh*source["h"]
    primordial_delta2=source["As"]*(k_mpc/0.05)**(source["ns"]-1.0)
    P=2*np.pi**2/k_mpc**3*primordial_delta2
    raw=P*cb*C_KMS*theta/k_mpc
    return tophat(kh*R16)*raw

def replay_e21(source,target,label):
    got=filtered_Cv_at(KCHECK,source)
    rel=np.abs(got-target)/np.maximum(np.abs(target),1e-300)
    m=float(rel.max())
    require(m<2e-5,f"{label} exact-transfer E21 replay fails: {m}")
    return got,m

def hankel(k_mpc, coeff_mpc3, ell, s_h, h):
    # s_h is Mpc/h; transform in physical Mpc because coeff is Mpc^3.
    s_mpc=np.asarray(s_h,float)/h
    kk=k_mpc[:,None]
    ss=s_mpc[None,:]
    integrand=(k_mpc*k_mpc*coeff_mpc3)[:,None]*spherical_jn(ell,kk*ss)
    return np.trapezoid(integrand,k_mpc,axis=0)/(2*np.pi**2)

def e3_q_vectors():
    x=json.loads(E3.read_text(encoding="utf-8"))
    require(x.get("observed_odd_data_vector_read") is False,"Observed odd guard changed")
    out={}
    for case_name,case in x["cases"].items():
        inj=case["levels"]["full_1200"]["DD_pair_level_synthetic_odd_injections"]["0.02"]
        m=inj["full_ell0to3_injected_increment_only_if_144_supported"]
        require(m is not None,"Missing full E3 injected multipoles "+case_name)
        qv=[]
        for ell in ("1","3"):
            ev=np.asarray(m[ell]["expected_finite_bin_pair_analytic_increment_by_s_bin"],float)/.02
            require(ev.shape==(6,),"Bad E3 q shape")
            qv.extend(ev.tolist())
        out[case_name]=np.asarray(qv,float)
    require(len(out)==18,"Expected 18 E3 q vectors")
    return out

def get_blocks(npz,cap,mid):
    d={}
    for o in (1,3):
        for i in (1,3):
            key=f"mock_M_{cap}_id{mid:04d}_out{o}_in{i}"
            require(key in npz.files,"Missing NPZ block "+key)
            a=np.asarray(npz[key],float)
            require(a.shape==(6,120),"Bad block shape "+key)
            d[(o,i)]=a
    return d

def project(blocks,xi1,xi3):
    y=[]
    for o in (1,3):
        v=blocks[(o,1)]@xi1 + blocks[(o,3)]@xi3
        y.extend(v.tolist())
    return np.asarray(y,float)

def recovery_fraction(rows,cap,coef_key,B):
    rr=[r for r in rows if r["cap"]==cap]
    n=len(rr)
    okp=okm=0
    for r in rr:
        a=float(r[coef_key])*B
        base=float(r["baseline_A0"])
        # symmetric +/- physical coefficient sign
        if a!=0 and (base+a)*a>0: okp+=1
        am=-a
        if am!=0 and (base+am)*am>0: okm+=1
    return 0.5*(okp+okm)/n

def threshold_scan(rows,cap,coef_key,target):
    # logarithmic search; descriptive nine-mock threshold only.
    grid=np.logspace(-6,6,24001)
    for B in grid:
        if recovery_fraction(rows,cap,coef_key,float(B)) >= target:
            return float(B), recovery_fraction(rows,cap,coef_key,float(B))
    return None,None

def main():
    require(E3.exists(),"MISSING "+str(E3))
    require(E45.exists(),"MISSING "+str(E45))
    e45=json.loads(E45.read_text(encoding="utf-8"))
    require(e45["status"]=="PASS_DESCRIPTIVE_NOT_PHYSICAL_WAKE_NOT_INFERENCE","E45 status changed")
    A0={r["case"]:float(r["baseline_apparent_amplitude"]) for r in e45["rows"]}
    require(len(A0)==18,"E45 needs 18 cases")

    npzb, npz_source, npz_mode=find_npz_bytes()
    npz=np.load(io.BytesIO(npzb),allow_pickle=False)
    require(len(npz.files)==203,"Wrong E0/E1 NPZ array count")

    q,fp,fm=load_e8()
    sp=dense_transfer_source(q,fp,"Fplus")
    sm=dense_transfer_source(q,fm,"Fminus")
    require(abs(sp["h"]-sm["h"])<1e-13,"F states changed h")
    h=sp["h"]
    gotp,gapp=replay_e21(sp,E21_PLUS,"Fplus")
    gotm,gapm=replay_e21(sm,E21_MINUS,"Fminus")

    # Common deterministic integration grid within frozen CLASS k coverage.
    # Form C_v only after separately interpolating the transfer fields.
    klo=max(float(sp["kin"].min()),float(sm["kin"].min()),1e-4)
    khi=min(float(sp["kin"].max()),float(sm["kin"].max()),.25)
    require(klo<.001 and khi>=.249,"Dense CLASS source does not cover required range")
    kh=np.geomspace(klo,khi,2048)
    kp=kh*h
    cvp=filtered_Cv_at(kh,sp)
    cvm=filtered_Cv_at(kh,sm)
    dCU=(cvm-cvp)/C_KMS  # Mpc^3; Fminus-Fplus R16-smoothed dimensionless-velocity cross source

    sedges=np.asarray(npz["fine_sedges_mpc_over_h"],float)
    require(sedges.shape==(121,),"Unexpected fine separation edges")
    sf=.5*(sedges[:-1]+sedges[1:])

    H1=hankel(kp,dCU,1,sf,h)
    H3=hankel(kp,dCU,3,sf,h)

    # Fourier odd contrast:
    # i*C_U [g0*mu + g2*mu^3]
    # = i*C_U [(g0+3g2/5)P1 + (2g2/5)P3].
    # xi_l = i^l integral P_l j_l => xi1=-..., xi3=+...
    xi1_g0=-H1
    xi3_g0=np.zeros_like(H1)
    xi1_g2=-(3/5)*H1
    xi3_g2=(2/5)*H3

    qvec=e3_q_vectors()
    rows=[]
    for case_name in sorted(qvec):
        mid_s,cap=case_name.split("/")
        mid=int(mid_s)
        blocks=get_blocks(npz,cap,mid)
        y0=project(blocks,xi1_g0,xi3_g0)
        y2=project(blocks,xi1_g2,xi3_g2)
        qv=qvec[case_name]
        qq=float(qv@qv)
        require(qq>0,"Zero E45 q norm")
        alpha=float(qv@y0/qq)
        beta=float(qv@y2/qq)
        rows.append({
            "case":case_name,"cap":cap,"mock_id":mid,
            "baseline_A0":A0[case_name],
            "alpha_A_per_unit_g0":alpha,
            "beta_A_per_unit_g2":beta,
            "unit_g0_postwindow_12d":y0.tolist(),
            "unit_g2_postwindow_12d":y2.tolist(),
        })

    # Descriptive required coefficient magnitudes on two one-parameter axes.
    summaries={}
    for cap in ("NGC","SGC"):
        rc=[r for r in rows if r["cap"]==cap]
        al=[r["alpha_A_per_unit_g0"] for r in rc]
        be=[r["beta_A_per_unit_g2"] for r in rc]
        s={
          "n":len(rc),
          "alpha_g0_median":statistics.median(al),
          "alpha_g0_minmax":[min(al),max(al)],
          "beta_g2_median":statistics.median(be),
          "beta_g2_minmax":[min(be),max(be)],
        }
        for key in ("alpha_A_per_unit_g0","beta_A_per_unit_g2"):
            tag="g0" if key.startswith("alpha") else "g2"
            for targ,label in ((.8333333333333333,"83pct"),(.9444444444444444,"94pct")):
                B,rec=threshold_scan(rows,cap,key,targ)
                s[f"required_abs_{tag}_for_{label}_symmetric_sign_recovery"]=B
                s[f"actual_discrete_{tag}_{label}_recovery"]=rec
        summaries[cap]=s

    out={
      "stage":"E47_LOCAL_FROZEN_FPM_LINEAR_RESPONSE_TO_EBOSS_E45_BRIDGE",
      "status":"PASS_PARAMETERIZED_PHYSICAL_RESPONSE_BASIS_NOT_CALIBRATED_GALAXY_MODEL",
      "response_parameters":{
        "g0":"b_L*d_E - b_E*d_L",
        "g2":"f*(d_E-d_L)",
        "d_a":"E25 dimensionless local neutrino-relative-velocity number-count response coefficient",
        "prediction":"A_E45(case)=alpha_case*g0 + beta_case*g2 for the Fminus-Fplus source contrast, under state-independent tracer response"
      },
      "source_replay":{
        "Fplus_E21_at_four_K":gotp.tolist(),
        "Fminus_E21_at_four_K":gotm.tolist(),
        "max_relative_gap_Fplus":gapp,
        "max_relative_gap_Fminus":gapm,
      },
      "integration":{
        "source_construction":"transfer-first interpolation, then nonlinear cross product, then W_R16 filter; exact E12/E13 ordering",
        "R16_Mpc_over_h":R16,
        "z":Z,"h":h,"k_h_min":float(kh[0]),"k_h_max":float(kh[-1]),
        "k_nodes":len(kh),"fine_s_centers_Mpc_over_h":sf.tolist()
      },
      "operator":{
        "source":npz_source,"mode":npz_mode,"sha256":NPZ_SHA,"bytes":NPZ_BYTES
      },
      "cap_summary":summaries,
      "rows":rows,
      "guardrails":{
        "observed_odd_used":False,
        "observed_galaxy_rows_used":False,
        "FITS_used":False,
        "ASDF_used":False,
        "physical_g0_g2_calibrated":False,
        "full_number_count_model":False,
        "even_F_state_changes_included":False,
        "state_dependent_response_included":False,
        "inference_or_sigma":False
      }
    }
    OUT.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")

    print("E47_WSL_PHYSICAL_RESPONSE_BASIS_PASS")
    print("E21_REPLAY_MAX_REL_GAP_PLUS",gapp)
    print("E21_REPLAY_MAX_REL_GAP_MINUS",gapm)
    print("K_H_RANGE",kh[0],kh[-1],"N",len(kh))
    for cap in ("NGC","SGC"):
        s=summaries[cap]
        print("CAP",cap)
        print("ALPHA_G0_MEDIAN",s["alpha_g0_median"],"RANGE",s["alpha_g0_minmax"])
        print("BETA_G2_MEDIAN",s["beta_g2_median"],"RANGE",s["beta_g2_minmax"])
        print("REQ_G0_83",s["required_abs_g0_for_83pct_symmetric_sign_recovery"])
        print("REQ_G0_94",s["required_abs_g0_for_94pct_symmetric_sign_recovery"])
        print("REQ_G2_83",s["required_abs_g2_for_83pct_symmetric_sign_recovery"])
        print("REQ_G2_94",s["required_abs_g2_for_94pct_symmetric_sign_recovery"])
    print("NO_FITS_ASDF_OBSERVED_ODD",True)
    print("REPORT",OUT)

if __name__=="__main__":
    main()
