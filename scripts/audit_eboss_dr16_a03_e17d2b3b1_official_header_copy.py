#!/usr/bin/env python3
"""E17D2b3b1: OFFICIAL EMBEDDED HEADER-DERIVED metadata only.

Decode the Git-blob-verified 9.5MB abacusutils header COPY for exact candidate
AbacusSummit_base_c000_ph000 state key z0.950. Never access concrete halo ASDF,
particles/PIDs/merger trees, mock, observed odd, or CLASS. Even successful
extraction from a bundled header COPY is NOT actual per-file halo ASDF QA.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PROTO=ROOT/"source_data/eboss_dr16_a03_e17d2b3b1_official_embedded_header_metadata_prereg_2026-09-28.json"
PROTO_BLOB="a282b157a80f88a1e0bf95c01942453d9a59efe4"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
E17A=ROOT/"source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json"
B2=ROOT/"source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28/e17d2b2_original_abacus_c000_documented_metadata_bridge.json"
B3A=ROOT/"source_data/eboss_dr16_a03_e17d2b3a_archived_CI_2026_09_28/e17d2b3a_original_l1_200critical_nonidentifiability.json"
B0D=ROOT/"source_data/eboss_dr16_a03_e17d2b3b0_archived_CI_2026_09_28"
B0=B0D/"e17d2b3b0_offline_real_abacus_access_contract_BLOCKED.json"
B0MAN=B0D/"archive_manifest.json"
B0PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b0_abacus_real_data_access_dryrun_prereg_2026-09-28.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3b1_official_header_copy"
REPORT="e17d2b3b1_pinned_official_embedded_header_copy_state_z095.json"
OFFICIAL_COMMIT="24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6"

def require(ok,msg):
    if not ok:raise ValueError("E17D2B3B1_FAIL_CLOSED: "+msg)
def sha(b):return hashlib.sha256(b).hexdigest()
def blob(b):return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def write_once(path,obj):
    raw=(json.dumps(obj,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"immutable copied-header report would change")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3b1_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B3B1_HEADER_COPY_REPORT",path,"SHA256",sha(raw),flush=True)

def gate(upstream_root,skip_upstream=False):
    require(blob(PROTO.read_bytes())==PROTO_BLOB,"pre-numerical official metadata protocol Git blob")
    p=json.loads(PROTO.read_bytes());locked=p["immutable_project_parents"]
    for f,key in ((E8,"E8_original_4000q_sha256"),
                  (E16,"E16_original_576_triangle_sha256"),
                  (E17A,"E17A_original_CLASS_joint_sha256"),
                  (B2,"E17D2b2_official_metadata_original_sha256"),
                  (B3A,"E17D2B3A_original_mass_nonidentifiability_sha256"),
                  (B0,"E17D2B3B0_offline_report_sha256")):
        require(sha(f.read_bytes())==locked[key],"immutable original SHA: "+key)
    for f,key in ((B0MAN,"E17D2B3B0_archive_manifest_git_blob"),
                  (B0PRE,"E17D2B3B0_protocol_git_blob")):
        require(blob(f.read_bytes())==locked[key],"immutable original protocol/manifest Git blob "+key)
    parent=json.loads(B0.read_bytes())
    a=json.loads(E17A.read_bytes())
    require(parent["offline_decision"]["real_ABACUS_data_ready"] is False
            and parent["no_actual_ASDF_halo_or_tree_downloaded"] is True
            and parent["observed_odd_SEALED"] is True
            and parent["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and a["eBOSS_observed_odd_read"] is False
            and a["retarded_halo_tracer_physical_bispectrum"]=="BLOCKED"
            and all(p["absolute_STOP"].values()),
            "parent original no-read/observed seal")
    src=p["upstream_public_source_pins"]
    require(src["public_commit"]==OFFICIAL_COMMIT
            and src["metadata_exact_size_bytes"]==9486279
            and src["metadata_git_blob"]=="0065029cdd515ffa76c61a27b81e9f16b2656f38",
            "upstream revision and bounded public metadata source")
    if skip_upstream:return p,None
    m=upstream_root/src["metadata_filename"]
    extension=upstream_root/src["blosc_extension_python"]
    decoder=upstream_root/src["official_metadata_decoder_python"]
    require(m.exists() and m.is_file() and m.stat().st_size==9486279,
            "not exact official 9.5MB SOURCE-BUNDLED metadata copy")
    require(blob(m.read_bytes())==src["metadata_git_blob"]
            and blob(extension.read_bytes())==src["blosc_extension_git_blob"]
            and blob(decoder.read_bytes())==src["official_metadata_decoder_git_blob"],
            "official compressed Git metadata or decompressor version mismatch")
    # No provider directory probes, no actual files, no halo catalog.
    return p,m

def scalar(x):
    try:
        import numpy as np
        if isinstance(x,np.generic):x=x.item()
    except ImportError:pass
    return type(x) in (int,float) and math.isfinite(float(x))

def read_official_embedded_metadata(p,m,upstream):
    import asdf
    import msgpack
    src=p["upstream_public_source_pins"]
    extfile=upstream/src["blosc_extension_python"]
    spec=importlib.util.spec_from_file_location("e17d2b3b1_pinned_upstream_blosc",extfile)
    require(spec is not None and spec.loader is not None,"official pinned blosc extension missing")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    ext=module.AbacusExtension()
    sim=p["prospective_locked_access"]["exact_simname"]
    state_label=p["prospective_locked_access"][
        "exact_state_label_from_official_get_meta_z_normalization"]
    with asdf.open(m,extensions=[ext],memmap=False) as af:
        require(sim in af.tree,"exact c000 ph000 absent from upstream embedded metadata")
        simtree=af.tree[sim]
        params=msgpack.loads(bytes(simtree["param"].data),strict_map_key=False)
        states=msgpack.loads(bytes(simtree["state"].data),strict_map_key=False)
    require(sim in [params.get("SimName"),sim] and params.get("SimName")==sim,
            "embedded source metadata simulation identity absent")
    require(state_label in states,
            "exact z0.950 source HEADER-COPY state is absent; do not substitute other epoch")
    state=states[state_label]
    require(isinstance(state,dict) and isinstance(params,dict),
            "official source state/param decoder returned wrong types")
    scalar_candidates={}
    for key,value in state.items():
        k=str(key)
        if any(token in k.lower() for token in ("redshift","scalefactor","scale_factor","anow","a_now","timenow")):
            if scalar(value):
                scalar_candidates[k]=float(value)
            else:
                scalar_candidates[k]="NON_SCALAR_OR_NONFINITE_NOT_USED"
    exact_candidate_keys=("Redshift","redshift","RedshiftNow","CurrentRedshift","z","z_now","zNow")
    detected={k:float(state[k]) for k in exact_candidate_keys if k in state and scalar(state[k])}
    def derive_from_a():
        for k in ("aNow","a_now","ScaleFactor","scale_factor"):
            if k in state and scalar(state[k]) and 0<float(state[k])<=1.0:
                return {k:1./float(state[k])-1.}
        return {}
    from_a=derive_from_a()
    copied_z_status="BLOCKED_NO_UNAMBIGUOUS_NUMERIC_RED_SHIFT_KEY_IN_SOURCE_COPY"
    copied_z=None
    if detected:
        values=list(detected.values())
        if all(abs(z-values[0])<1e-10 for z in values):
            copied_z=values[0];copied_z_status="SOURCE_BUNDLED_COPIED_HEADER_NUMERIC_Z_UNAMBIGUOUS_NOT_REAL_HALO_ASDF"
        else:
            copied_z_status="BLOCKED_MULTIPLE_INCONSISTENT_COPIED_HEADER_Z_VALUES"
    elif len(from_a)==1:
        copied_z=list(from_a.values())[0]
        copied_z_status="SOURCE_BUNDLED_COPIED_HEADER_SCALE_FACTOR_DERIVED_Z_NOT_REAL_HALO_ASDF"
    if copied_z is not None:
        require(0.<=copied_z<10. and math.isfinite(copied_z),
                "unphysical or nonfinite embedded copied-header z")
        if from_a:
            require(all(abs(val-copied_z)<1e-8 for val in from_a.values()),
                    "copied-header Redshift vs scale factor mismatch")
    original_meta_bundle_sha=sha(m.read_bytes())
    # independent decompression: bytewise reopen under separately constructed
    # ASDF container and verify same selected exact scalar fields, not new data.
    with asdf.open(m,extensions=[module.AbacusExtension()],memmap=False) as af2:
        verify=msgpack.loads(bytes(af2.tree[sim]["state"].data),strict_map_key=False)
    require(set(verify.keys())==set(states.keys())
            and verify[state_label]==state,
            "second full official source-copy metadata decode drift")
    out={
      "date":"2026-09-28",
      "status":"E17D2B3B1_OFFICIAL_GIT_PINNED_ABACUS_EMBEDDED_HEADER_COPY_DECODED_NO_ACTUAL_HALO_ASDF_VERIFIED",
      "prospective_protocol_git_blob":PROTO_BLOB,
      "upstream_public_official_commit":OFFICIAL_COMMIT,
      "upstream_source_git_blob":src["metadata_git_blob"],
      "upstream_source_size_bytes":m.stat().st_size,
      "upstream_source_file_sha256":original_meta_bundle_sha,
      "official_pinned_decoder_git_blob":src["official_metadata_decoder_git_blob"],
      "official_pinned_blosc_extension_git_blob":src["blosc_extension_git_blob"],
      "upstream_simulation":sim,
      "exact_source_copy_state_label":state_label,
      "upstream_parameter_SimName":params.get("SimName"),
      "upstream_state_key_inventory_at_z095":sorted(str(k) for k in state),
      "upstream_copy_redshift_related_scalar_candidates":scalar_candidates,
      "explicit_redshift_scalar_candidates":detected,
      "independently_derived_redshift_from_copied_scalefactor_if_present":from_a,
      "copied_header_z_extraction_status":copied_z_status,
      "copied_header_z_value_if_unambiguous":copied_z,
      "state_in_official_bundle_exists":True,
      "full_upstream_state_second_decode_identical":True,
      "real_halo_info_ASDF_header_and_GNU_cksum_SHA256_verified":False,
      "actual_halo_or_merger_tree_or_particle_data_downloaded":False,
      "actual_M200c_200critical_same_halo_history_verified":False,
      "original_CLASS_As_tau_Fplus_Fminus_nonlinear_wake_matched":False,
      "real_LRG_ELG_HOD_or_3pt_window_covariance_verified":False,
      "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
      "observed_odd_SEALED":True,
      "main_untouched_PR_draft":True}
    return out

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--upstream-root",type=Path,default=Path("official-abacusutils"))
    ap.add_argument("--output-dir",type=Path,default=OUT)
    ap.add_argument("--preflight",action="store_true")
    a=ap.parse_args()
    p,m=gate(a.upstream_root,skip_upstream=a.preflight)
    if a.preflight:
        print("E17D2B3B1_ORIGINAL_FROZEN_SHA_AND_OFFICIAL_PIN_PREFLIGHT_PASS",
              "UPSTREAM_METADATA_COPY_BYTES_MAX",9486279,
              "NO_HALO_ASDF_NO_CLASS_NO_ODD",flush=True)
        return
    out=read_official_embedded_metadata(p,m,a.upstream_root)
    write_once(a.output_dir/REPORT,out)
    print("E17D2B3B1_OFFICIAL_BUNDLED_HEADER_COPY_SOURCE_PASS",
          "STATE",out["exact_source_copy_state_label"],
          "COPIED_HEADER_Z_STATUS",out["copied_header_z_extraction_status"],
          "Z_IF_PROVEN",out["copied_header_z_value_if_unambiguous"],
          "ACTUAL_HALO_ASDF_AND_M200C_UNVERIFIED",flush=True)
if __name__=="__main__":main()
