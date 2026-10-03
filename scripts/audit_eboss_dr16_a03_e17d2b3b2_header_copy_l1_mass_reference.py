#!/usr/bin/env python3
"""E17D2b3b2: source-only precise copied-header CompaSO L1 threshold audit.

Reads only immutable 9.5MB PUBLIC GIT-BUNDLED original-header COPY from
abacusorg/abacusutils. NO concrete halo_info ASDF, particle/PID, merger tree,
original E8/F sources, new CLASS, observed galaxy/odd or astrophysical M200c.
The returned threshold ratio does NOT determine any individual halo mass.
"""
from __future__ import annotations
import argparse
from decimal import Decimal,localcontext
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_official_header_copy_mass_reference_prereg_2026-09-28.json"
P_BLOB="324e722200997b21a0fa97d965b8ef096e1587e6"
B1DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b3b1_archived_CI_2026_09_28"
B1=B1DIR/"e17d2b3b1_pinned_official_embedded_header_copy_state_z095.json"
B1MAN=B1DIR/"archive_manifest.json"
B1PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b1_official_embedded_header_metadata_prereg_2026-09-28.json"
B1RUN=ROOT/"scripts/audit_eboss_dr16_a03_e17d2b3b1_official_header_copy.py"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3b2_header_copy_so"
OUTNAME="e17d2b3b2_pinned_header_copy_l1_threshold_vs_200crit.json"

def check(ok,why):
    if not ok: raise ValueError("E17D2B3B2_FAIL_CLOSED: "+why)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def gitblob(raw):return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def save(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():check(path.read_bytes()==raw,"refuse to overwrite changed original report")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3b2_",delete=False) as f:
            tmp=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(tmp,path)
        finally:tmp.unlink(missing_ok=True)
    print("E17D2B3B2_SOURCE_COPY_REPORT",path,"SHA256",sha(raw),flush=True)

def preflight():
    check(gitblob(P.read_bytes())==P_BLOB,"precommitted protocol Git blob")
    p=json.loads(P.read_bytes());lock=p["prior_immutable"]
    for f,key in ((E8,"E8_csv_sha256"),(E16,"E16_geometry_sha256"),
                  (E17A,"E17A_joint_sha256"),
                  (B1,"E17D2b3b1_source_copy_original_sha256")):
        check(sha(f.read_bytes())==lock[key],"original parent SHA "+key)
    for f,key in ((B1MAN,"E17D2b3b1_manifest_git_blob"),
                  (B1PRE,"E17D2b3b1_amended_prereg_git_blob"),
                  (B1RUN,"E17D2b3b1_code_git_blob")):
        check(gitblob(f.read_bytes())==lock[key],"original parent Git blob "+key)
    previous=json.loads(B1.read_bytes())
    check(previous["observed_odd_SEALED"] is True and
          previous["real_halo_info_ASDF_header_and_GNU_cksum_SHA256_verified"] is False
          and previous["actual_halo_or_merger_tree_or_particle_data_downloaded"] is False
          and previous["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
          and previous["upstream_source_file_sha256"]==
              p["official_upstream"]["metadata_sha256"]
          and previous["exact_source_copy_state_label"]=="z0.950"
          and all(p["absolute_STOP"].values()),
          "original b1 source-copy and hard real-data scope")
    check(p["frozen_source_copy_inputs"]["require_scalar_fields"]==
          ["Redshift","ScaleFactor","OmegaNow_m","SODensityL1",
           "ParticleMassHMsun","ParticleMassMsun"],
          "six exact fields were changed")
    return p,previous

def scalar(v):
    if type(v) in (int,float):
        q=float(v)
    elif hasattr(v,"item"):
        q=v.item()
        check(type(q) in (int,float),"copy field not scalar numeric")
        q=float(q)
    else:
        raise ValueError("E17D2B3B2_FAIL_CLOSED: nonnumeric source-copy state field")
    check(math.isfinite(q),"copy field not finite")
    return q

def validate_fields(v,p,previous):
    expected=p["frozen_source_copy_inputs"]
    for k in expected["require_scalar_fields"]:
        check(k in v,"missing source-copy header field "+k)
    data={k:scalar(v[k]) for k in expected["require_scalar_fields"]}
    check(data["Redshift"]>=0 and data["ScaleFactor"]>0
          and data["ScaleFactor"]<=1 and 0<data["OmegaNow_m"]<1
          and data["SODensityL1"]>0
          and data["ParticleMassHMsun"]>0
          and data["ParticleMassMsun"]>0,
          "unphysical copy epoch, L1 threshold or particle mass")
    for k,previous_key,limit in (
        ("Redshift","expected_previous_copy_redshift",
         p["prereg_QA"]["copied_redshift_prior_absolute_diff_tolerance"]),
        ("ScaleFactor","expected_previous_copy_scale_factor",
         p["prereg_QA"]["copied_scale_factor_prior_absolute_diff_tolerance"]),
        ("OmegaNow_m","expected_previous_copy_omega_now_m",
         p["prereg_QA"]["copied_Omega_now_m_prior_absolute_diff_tolerance"])
    ):
        check(abs(data[k]-expected[previous_key])<limit,
              "copied header deviates from prereg prior: "+k)
    check(abs(1/data["ScaleFactor"]-1-data["Redshift"])<
          p["prereg_QA"]["redshift_scale_consistency_abs_tolerance"],
          "copied source z and a disagree")
    check(data["Redshift"]==previous["copied_header_z_value_if_unambiguous"],
          "source-copy Redshift not identical to previously SHA-pinned result")
    return data

def inspect(p,previous,upstream):
    import asdf
    import msgpack
    orig=p["official_upstream"]
    m=upstream/orig["metadata_path"]
    extpath=upstream/orig["pinned_blosc_source_path"]
    decpath=upstream/orig["pinned_upstream_decoder_path"]
    check(m.is_file() and m.stat().st_size==orig["metadata_exact_bytes"]
          and gitblob(m.read_bytes())==orig["metadata_git_blob"]
          and sha(m.read_bytes())==orig["metadata_sha256"]
          and gitblob(extpath.read_bytes())==orig["blosc_source_git_blob"]
          and gitblob(decpath.read_bytes())==orig["decoder_git_blob"],
          "official upstream source-copy byte size/Git blob/SHA or decoder drift")
    spec=importlib.util.spec_from_file_location("e17d2b3b2_upstream_blosc",extpath)
    check(spec and spec.loader,"official pinned Blosc extension not importable")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sim=orig["simname"];zkey=orig["state_key"]
    with asdf.open(m,extensions=[module.AbacusExtension()],memmap=False) as af:
        check(sim in af.tree,"source-copy sim missing")
        source=af.tree[sim]
        params=msgpack.loads(bytes(source["param"].data),strict_map_key=False)
        states=msgpack.loads(bytes(source["state"].data),strict_map_key=False)
    check(params.get("SimName")==sim and zkey in states,
          "source-copy sim or state identity mismatch")
    data=validate_fields(states[zkey],p,previous)
    # Independent decoding and Decimal arithmetic on the SAME copied bytes,
    # not a second astrophysical source or an independent N-body simulation.
    with asdf.open(m,extensions=[module.AbacusExtension()],memmap=False) as af:
        rerun=msgpack.loads(bytes(af.tree[sim]["state"].data),
                             strict_map_key=False)
    data_second=validate_fields(rerun[zkey],p,previous)
    check(data_second==data,"two immutable upstream source-copy decodes differ")
    threshold_c=200/data["OmegaNow_m"]
    ratio=data["SODensityL1"]/threshold_c
    with localcontext() as ctx:
        ctx.prec=40
        D=Decimal
        source_ratio=D(str(data["SODensityL1"]))*D(str(data["OmegaNow_m"]))/D("200")
        source_target=D("200")/D(str(data["OmegaNow_m"]))
        source_z=D("1")/D(str(data["ScaleFactor"]))-D("1")
        relative_arithmetic_gap=abs(float(source_ratio)-ratio)/max(1.,abs(ratio))
        target_arithmetic_gap=abs(float(source_target)-threshold_c)/max(1.,abs(threshold_c))
    check(math.isfinite(ratio) and ratio>0
          and relative_arithmetic_gap<1e-13
          and target_arithmetic_gap<1e-13
          and abs(float(source_z)-data["Redshift"])<1e-10,
          "copied threshold arithmetic or redshift independent Decimal mismatch")
    negative=[]
    for edit,label in (({"SODensityL1":0},"ZERO_L1_THRESHOLD"),
                       ({"Redshift":.95},"FABRICATED_DIRECTORY_REDSHIFT")):
        wrong=dict(data);wrong.update(edit)
        try:validate_fields(wrong,p,previous)
        except ValueError:negative.append(label)
        else:raise AssertionError("deliberately invalid "+label+" was accepted")
    check(len(negative)==2,"two preregistered negative controls not rejected")
    return {
      "date":"2026-09-28",
      "status":"E17D2B3B2_PINNED_OFFICIAL_EMBEDDED_HEADER_COPY_L1_THRESHOLD_AND_PARTICLE_MASS_MATCH_ORIGINAL_SOURCE_NO_ACTUAL_HALO_ASDF",
      "prospective_protocol_git_blob":P_BLOB,
      "original_E17D2b3b1_source_copy_sha256":p["prior_immutable"][
          "E17D2b3b1_source_copy_original_sha256"],
      "official_commit":orig["commit"],
      "official_metadata_git_blob":orig["metadata_git_blob"],
      "official_metadata_sha256":orig["metadata_sha256"],
      "official_metadata_bytes":m.stat().st_size,
      "simname":sim,"state_key":zkey,
      "official_source_copy_six_exact_numeric_fields":data,
      "copy_delta_200critical_in_mean_cosmic_density_units":threshold_c,
      "copy_L1_SO_threshold_over_200critical_SO_threshold":ratio,
      "ratio_interpretation":"ONLY threshold-density comparison using official mean cosmic density and source-copy OmegaNow_m; NOT an individual halo mass, radius conversion, neutrino wake or measured LRG/ELG quantity.",
      "Decimal_40digit_independent_copy_threshold_ratio_gap":relative_arithmetic_gap,
      "Decimal_40digit_independent_copy_target_gap":target_arithmetic_gap,
      "two_same_pinned_official_copy_decodes_agree":True,
      "synthetic_tampered_copy_field_negative_tests_rejected":negative,
      "actual_halo_info_ASDF_byte_checksum_header_verified":False,
      "actual_provider_halo_file_GNU_cksum_or_merger_tree_read":False,
      "same_object_true_M200c_200critical_remeasured":False,
      "original_frozen_E8_E16_and_observed_odd_unchanged":True,
      "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
      "observed_odd_SEALED":True,
      "main_untouched_PR_draft":True
    }

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preflight",action="store_true")
    ap.add_argument("--upstream-root",type=Path,default=Path("official-abacusutils"))
    ap.add_argument("--output-dir",type=Path,default=OUT)
    a=ap.parse_args()
    p,prev=preflight()
    if a.preflight:
        print("E17D2B3B2_SOURCE_COPY_SHA_PREREG_AND_PARENT_SEAL_PASS",
              "OFFICIAL_DATA_NOT_YET_READ",flush=True)
        return
    o=inspect(p,prev,a.upstream_root)
    save(a.output_dir/OUTNAME,o)
    print("E17D2B3B2_C000_COPIED_HEADER_SO_THRESHOLD_PASS",
          "z",o["official_source_copy_six_exact_numeric_fields"]["Redshift"],
          "SODensityL1",o["official_source_copy_six_exact_numeric_fields"]["SODensityL1"],
          "Mparticle_Msun_over_h",o["official_source_copy_six_exact_numeric_fields"]["ParticleMassHMsun"],
          "L1_TO_200CRIT_DENSITY_THRESHOLD_RATIO",
          o["copy_L1_SO_threshold_over_200critical_SO_threshold"],
          "REAL_HALO_ASDF_AND_M200C_BLOCKED",flush=True)
if __name__=="__main__":main()
