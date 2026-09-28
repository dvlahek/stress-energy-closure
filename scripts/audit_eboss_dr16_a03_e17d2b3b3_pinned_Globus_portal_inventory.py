#!/usr/bin/env python3
"""E17D2b3b3: OFFLINE Git-pinned public Abacus portal inventory, NO HALO ASDF.

Reads ONLY the upstream source-controlled static 5.1 MB aggregate portal
manifest, its public config/build/view Python source, and previously frozen
local science metadata. DOES NOT contact Globus, load ASDF/catalog/particles
or authenticate a provider file. Static portal manifest != live availability.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from urllib.parse import urlencode

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"source_data/eboss_dr16_a03_e17d2b3b3_official_Globus_portal_product_inventory_prereg_2026-09-28.json"
P_BLOB="42dab76d5e4fb7c8325ece9ea4465775b384da21"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
B1=ROOT/"source_data/eboss_dr16_a03_e17d2b3b1_archived_CI_2026_09_28/e17d2b3b1_pinned_official_embedded_header_copy_state_z095.json"
B2=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28/e17d2b3b2_pinned_header_copy_l1_threshold_vs_200crit.json"
B2MAN=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28/archive_manifest.json"
B2PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_official_header_copy_mass_reference_prereg_2026-09-28.json"
R1PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2r1_real_ASDF_byte_header_verifier_prereg_2026-09-28.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3b3_official_portal"
OUTNAME="e17d2b3b3_original_pinned_portal_c000_ph000_z095_aggregate_inventory.json"

def need(ok,what):
    if not ok:raise ValueError("E17D2B3B3_FAIL_CLOSED: "+what)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def preflight():
    need(blob(P.read_bytes())==P_BLOB,"prospectively registered portal source protocol Git blob")
    p=json.loads(P.read_bytes())
    pin=p["immutable_local_parents"]
    for path,key in ((E8,"E8_csv_sha256"),(E16,"E16_triangle_sha256"),
                     (B1,"E17D2b3b1_header_copy_sha256"),
                     (B2,"E17D2b3b2_source_copy_mass_reference_sha256")):
        need(sha(path.read_bytes())==pin[key],"immutable science parent SHA "+key)
    for path,key in ((B2MAN,"E17D2b3b2_archive_manifest_git_blob"),
                     (B2PRE,"E17D2b3b2_protocol_git_blob"),
                     (R1PRE,"E17D2b3b2r1_real_ASDF_checker_prereg_blob")):
        need(blob(path.read_bytes())==pin[key],"immutable science source Git blob "+key)
    a=json.loads(B1.read_bytes());b=json.loads(B2.read_bytes())
    need(a["copied_header_z_value_if_unambiguous"]==0.952838237036305
         and a["real_halo_info_ASDF_header_and_GNU_cksum_SHA256_verified"] is False
         and b["actual_halo_info_ASDF_byte_checksum_header_verified"] is False
         and b["actual_provider_halo_file_GNU_cksum_or_merger_tree_read"] is False
         and b["same_object_true_M200c_200critical_remeasured"] is False
         and a["observed_odd_SEALED"] is True
         and b["observed_odd_SEALED"] is True
         and a["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
         and b["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
         and all(p["absolute_STOP"].values()),
         "observed odd, original frozen science or no-real-halo gate violated")
    return p

def upstream(p,source):
    up=p["upstream_public_portal_pins"]
    expected=[
        (up["manifest_path"],up["manifest_git_blob"],up["manifest_reported_bytes"]),
        (up["manifest_builder_path"],up["manifest_builder_git_blob"],None),
        (up["portal_config_path"],up["portal_config_git_blob"],None),
        (up["portal_view_path"],up["portal_view_git_blob"],None)
    ]
    reports={}
    for name,expectedblob,nbytes in expected:
        f=source/name
        need(f.is_file() and not f.is_symlink(),"upstream pinned public portal file missing "+name)
        raw=f.read_bytes()
        need(blob(raw)==expectedblob,"portal upstream source Git blob mismatch "+name)
        if nbytes is not None:need(len(raw)==nbytes,"portal manifest bytes mismatch")
        reports[name]={"git_blob":expectedblob,"size_bytes":len(raw),"sha256":sha(raw)}
    m=json.loads((source/up["manifest_path"]).read_bytes())
    config=(source/up["portal_config_path"]).read_text(encoding="utf-8")
    view=(source/up["portal_view_path"]).read_text(encoding="utf-8")
    builder=(source/up["manifest_builder_path"]).read_text(encoding="utf-8")
    need("urlencode(dict(origin_id=endpoint_id,origin_path=endpoint_path))" in view
         and "DATASET_ENDPOINT_ID" in view
         and "source_endpoint_base / pathfmt.format(sim['root']) / zstr / ftype" in view,
         "pinned portal source does not document claimed origin-id/path or product tree")
    need("DEFAULT_PRODUCTS = dict(halos=dict(path='{}/halos'" in builder
         and "cleaning=dict(path='cleaning/{}'" in builder,
         "pinned builder product tree semantics drift")
    ids=re.findall(r"^\s*DATASET_ENDPOINT_ID\s*=\s*'([^']+)'",config,re.M)
    bases=re.findall(r"^\s*DATASET_ENDPOINT_BASE\s*=\s*'([^']+)'",config,re.M)
    need(ids==[up["endpoint_uuid_from_pinned_portal_source"]]
         and bases==[up["endpoint_base_from_pinned_source"]],
         "upstream official source Globus collection UUID/base mismatch")
    return m,reports,ids[0],bases[0]

def inventory(j,p,endpoint,base):
    up=p["upstream_public_portal_pins"]
    need(isinstance(j,dict) and isinstance(j.get("data"),list)
         and isinstance(j.get("products"),dict),"not a published table manifest")
    name=up["exact_simulation"]
    rows=[r for r in j["data"] if r.get("name")==name]
    need(len(rows)==1,"exact frozen simulation not unique in public portal snapshot")
    row=rows[0]
    need(row.get("root")==up["expected_source_tree_root"]
         and row.get("id")==0,
         "exact simulation path or archived row identity drift")
    need(up["exact_directory_label_z"]=="z0.950"
         and up["expected_source_manifest_z_key"]=="0.95"
         and "0.95" in row.get("halos",{})
         and "0.95" in row.get("cleaning",{}),
         "no exact archived candidate epoch product")
    want=p["preobserved_public_manifest_inventory_NOT_new_science_measurement"]
    pairs=[
        ("halos","halo_info","halo_info_files","halo_info_aggregate_bytes"),
        ("halos","halo_pid_A","halo_pid_A_files","halo_pid_A_aggregate_bytes"),
        ("cleaning","cleaned_halo_info","cleaned_halo_info_files","cleaned_halo_info_aggregate_bytes"),
        ("cleaning","cleaned_rvpid","cleaned_rvpid_files","cleaned_rvpid_aggregate_bytes")
    ]
    counts={}
    for family,product,fn,bn in pairs:
        y=row[family]["0.95"][product]
        need(type(y)==list and len(y)==2
             and type(y[0])==int and type(y[1])==int
             and y[0]==want[fn] and y[1]==want[bn],
             "archived upstream exact file count/aggregate bytes drift "+product)
        counts[product]={"file_count_in_archived_snapshot":y[0],
                         "total_bytes_for_all_files_in_archived_snapshot":y[1],
                         "NOTE":"NOT per-file byte size; NOT current live Globus item listing"}
    products=j["products"]
    need(products["halos"]["path"]=="{}/halos"
         and products["cleaning"]["path"]=="cleaning/{}"
         and "halo_info" in products["halos"]["ftypes"]
         and "cleaned_halo_info" in products["cleaning"]["ftypes"],
         "raw and cleaned portal product templates drift")
    raw="/"+products["halos"]["path"].format(row["root"])+"/z0.950/halo_info/"
    clean="/"+products["cleaning"]["path"].format(row["root"])+"/z0.950/cleaned_halo_info/"
    need(endpoint==up["endpoint_uuid_from_pinned_portal_source"]
         and base=="/"
         and raw=="/AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/"
         and clean=="/cleaning/AbacusSummit_base_c000_ph000/z0.950/cleaned_halo_info/",
         "official portal source product/collection path mismatch")
    return counts,{
        "official_portal_pinned_Globus_collection_UUID":endpoint,
        "official_source_raw_halo_info_directory":raw,
        "official_source_cleaned_halo_info_directory":clean,
        "browse_raw_halo_info_Globus_URL":"https://app.globus.org/file-manager?"+
                        urlencode({"origin_id":endpoint,"origin_path":raw}),
        "browse_cleaned_halo_info_Globus_URL":"https://app.globus.org/file-manager?"+
                        urlencode({"origin_id":endpoint,"origin_path":clean}),
        "source_URL_is_browse_ONLY_and_live_collection_availability_UNVERIFIED":True
    }

def qa_negative(j,p,endpoint,base):
    for mutation,label in [
        (lambda x:x["data"][0]["halos"]["0.95"]["halo_info"].__setitem__(0,34),"COUNT"),
        (lambda x:x["data"][0].__setitem__("name","ABACUS_FAKE"),"SIM_NAME"),
        (lambda x:x["products"]["halos"].__setitem__("path","wrong/path"),"PRODUCT_ROOT")
    ]:
        x=copy.deepcopy(j);mutation(x)
        try:inventory(x,p,endpoint,base)
        except ValueError:pass
        else:raise AssertionError("fake portal metadata accepted: "+label)
    try:inventory(j,p,"00000000-0000-0000-0000-000000000000",base)
    except ValueError:pass
    else:raise AssertionError("fake Globus collection UUID accepted")
    print("E17D2B3B3_FOUR_METADATA_TAMPER_NEGATIVE_CONTROLS_PASS",flush=True)

def write(path,j):
    raw=(json.dumps(j,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():need(path.read_bytes()==raw,"source-only report already exists with different bytes")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3b3_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B3B3_ORIGINAL_PUBLIC_PORTAL_INVENTORY_SHA256",sha(raw),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--upstream-root",type=Path,default=Path("official-abacus-portal"))
    ap.add_argument("--output-dir",type=Path,default=OUT)
    a=ap.parse_args();p=preflight()
    if a.self_test:
        print("E17D2B3B3_FROZEN_E8_E16_B1_B2_PORTAL_PREREG_SHA_PASS",
              "NO_GLOBUS_OR_HALO_FILE_ACCESS",flush=True)
        return
    m,raw,endpoint,base=upstream(p,a.upstream_root)
    q,links=inventory(m,p,endpoint,base)
    qa_negative(m,p,endpoint,base)
    result={
       "date":"2026-09-28",
       "status":"E17D2B3B3_PINNED_PUBLIC_PORTAL_SNAPSHOT_C000_PH000_Z095_35_HALO_INFO_ARCHIVED_NO_LIVE_FILE_ATTESTATION",
       "prospective_protocol_git_blob":P_BLOB,
       "upstream_source_commit":p["upstream_public_portal_pins"]["immutable_commit"],
       "upstream_commit_date_UTC":p["upstream_public_portal_pins"]["commit_date_UTC"],
       "exact_pinned_portal_source_file_git_blobs_and_SHA256":raw,
       "selection":{"simname":p["upstream_public_portal_pins"]["exact_simulation"],
                    "manifest_row_id":0,"manifest_epoch_key":"0.95",
                    "directory_label_ONLY":"z0.950",
                    "actual_source_COPY_z_from_prior_not_real_halo_ASDF":0.952838237036305},
       "archived_provider_portal_source_inventory_not_live":q,
       "official_source_code_browse_links_NOT_tested_with_live_Globus":links,
       "four_tampered_manifest_and_endpoint_negative_controls_REJECTED":True,
       "single_real_halo_info_file_or_supplied_cksum_selected":False,
       "actual_Globus_source_file_per_item_size_sha_POSIX_cksum_or_live_status_read":False,
       "actual_halo_ASDF_provider_file_or_tree_particle_data_downloaded":False,
       "real_M200c_same_halo_mass_history_verified":False,
       "original_Fplus_Fminus_nonlinear_neutrino_wake_from_Abacus":False,
       "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
       "observed_odd_SEALED":True,
       "main_untouched_PR_draft":True
    }
    write(a.output_dir/OUTNAME,result)
    print("E17D2B3B3_PINNED_PUBLIC_PORTAL_SOURCE_C000_PH000_Z095_PASS",
          "RAW_HALO_INFO_FILES",q["halo_info"]["file_count_in_archived_snapshot"],
          "CLEANED_HALO_INFO_FILES",q["cleaned_halo_info"]["file_count_in_archived_snapshot"],
          "DIRECT_GLOBUS_BROWSE_ONLY",links["browse_raw_halo_info_Globus_URL"],
          "LIVE_AVAILABILITY_UNKNOWN_REAL_M200C_B_BLOCKED",flush=True)
if __name__=="__main__":main()
