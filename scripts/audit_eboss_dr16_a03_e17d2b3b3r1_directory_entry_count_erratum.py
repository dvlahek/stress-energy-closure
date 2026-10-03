#!/usr/bin/env python3
"""E17D2b3b3R1: erratum to old portal '35 files' wording.

Official pinned upstream build_manifest.py counts ALL entries in product
directories via fpath.iterdir(), not only ASDF superslabs. Preserve old
archive byte-for-byte, write a separately SHA-pinned correction and PROHIBIT
assertions about exact real .asdf counts without provider directory listing.
No provider/Globus/halo/particle/observed data read or downloaded.
"""
from __future__ import annotations
import argparse
import ast
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"source_data/eboss_dr16_a03_e17d2b3b3r1_portal_directory_entries_not_ASDF_counts_erratum_prereg_2026-09-28.json"
P_BLOB="534142340197c34160cb970a53cc86d8f7bcf57c"
OLD_DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b3b3_archived_CI_2026_09_28"
OLD=OLD_DIR/"e17d2b3b3_original_pinned_portal_c000_ph000_z095_aggregate_inventory.json"
OLD_MAN=OLD_DIR/"archive_manifest.json"
OUT=ROOT/"eboss_workspace/a03_physics_source/e17d2b3b3r1_count_semantics"
OUTNAME="e17d2b3b3r1_corrected_directory_entry_count_not_ASDF_count.json"

def require(flag,message):
    if not flag:raise ValueError("E17D2B3B3R1_FAIL_CLOSED: "+message)
def sha(raw):return hashlib.sha256(raw).hexdigest()
def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def check_original():
    require(blob(P.read_bytes())==P_BLOB,"corrective prospective prereg source blob")
    p=json.loads(P.read_bytes());pin=p["prior_original_AUDIT_preserve"]
    require(sha(OLD.read_bytes())==pin["original_result_SHA256"],
            "old original portal source report must remain immutable")
    require(blob(OLD_MAN.read_bytes())==pin["original_archive_manifest_git_blob"],
            "old original portal archive manifest must remain immutable")
    old=json.loads(OLD.read_bytes())
    man=json.loads(OLD_MAN.read_bytes())
    require(man["CI_run_id"]==str(pin["original_CI_run_id"])
            and old["prospective_protocol_git_blob"]==
            pin["original_E17D2b3b3_prospective_prereg_git_blob"]
            and old["observed_odd_SEALED"] is True
            and old["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and old["actual_halo_ASDF_provider_file_or_tree_particle_data_downloaded"] is False
            and all(p["absolute_STOP"].values()),
            "old first-class source report and observed odd preservation")
    return p,old

def builder_count_semantics(source):
    tree=ast.parse(source)
    comps=[]
    for n in ast.walk(tree):
        if not isinstance(n,ast.ListComp):continue
        for gen in n.generators:
            f=gen.iter
            if (isinstance(f,ast.Call)
                and isinstance(f.func,ast.Attribute)
                and f.func.attr=="iterdir"
                and isinstance(f.func.value,ast.Name)
                and f.func.value.id=="fpath"):
                comps.append((n,gen))
    require(len(comps)==1,"unexpected change in official directory walker")
    n,g=comps[0]
    require(not g.ifs,"upstream published count now has extension or other filtering")
    require(isinstance(n.elt,ast.Attribute) and n.elt.attr=="st_size",
            "manifest product aggregate is not directly stat().st_size")
    foundlen=False;foundsum=False
    for a in ast.walk(tree):
        if (isinstance(a,ast.Assign) and len(a.targets)==1
            and isinstance(a.targets[0],ast.Name)):
            t=a.targets[0].id;val=a.value
            if (isinstance(val,ast.Call) and isinstance(val.func,ast.Name)
                and len(val.args)==1 and isinstance(val.args[0],ast.Name)
                and val.args[0].id=="du"):
                if t=="nfile" and val.func.id=="len":foundlen=True
                if t=="du" and val.func.id=="sum":foundsum=True
    require(foundlen and foundsum,
            "source builder no longer sets nfile=len(unfiltered entries), du=sum(unfiltered sizes)")
    return {"source_builder_AST_listcomp_uses_direct_fpath_iterdir_without_extension_filter":True,
            "source_builder_AST_nfile_is_len_all_directory_entries":True,
            "source_builder_AST_directory_total_bytes_includes_all_entry_sizes":True}

def run(upstream,p,old):
    src=p["upstream_public_pins"]
    builder=upstream/src["builder_path"]
    table=upstream/src["portal_snapshot_path"]
    require(builder.is_file() and table.is_file()
            and blob(builder.read_bytes())==src["builder_blob"]
            and blob(table.read_bytes())==src["portal_snapshot_blob"],
            "official upstream source Git pin mismatch or unavailable")
    code=builder.read_text(encoding="utf-8")
    semantics=builder_count_semantics(code)
    probe=code.replace(
        "du = [fn.stat().st_size for fn in fpath.iterdir()]",
        "du = [fn.stat().st_size for fn in fpath.iterdir() if fn.suffix == '.asdf']")
    require(probe!=code,"upstream expected unfiltered expression not found")
    try:builder_count_semantics(probe)
    except ValueError:pass
    else:raise AssertionError("synthetic fake ASDF-only source incorrectly accepted")
    j=json.loads(table.read_bytes())
    selected=[x for x in j["data"] if x.get("name")==src["expected_simname"]]
    require(len(selected)==1,"official source table sim not unique")
    row=selected[0]
    require(row.get("root")==src["expected_simname"]
            and src["expected_epoch_key"] in row["halos"]
            and src["expected_epoch_key"] in row["cleaning"],
            "source table simulation, z0.950 snapshot or product metadata drift")
    lookup={
       "halo_info":row["halos"]["0.95"]["halo_info"],
       "halo_pid_A":row["halos"]["0.95"]["halo_pid_A"],
       "cleaned_halo_info":row["cleaning"]["0.95"]["cleaned_halo_info"],
       "cleaned_rvpid":row["cleaning"]["0.95"]["cleaned_rvpid"]
    }
    expected=src["expected_counts_and_aggregate_dir_bytes"]
    oldentries=old["archived_provider_portal_source_inventory_not_live"]
    corrected={}
    for prod,nums in lookup.items():
        require(nums==expected[prod],"official unchanged source directory-count/total bytes drift "+prod)
        require(nums[0]==oldentries[prod]["file_count_in_archived_snapshot"]
                and nums[1]==oldentries[prod]["total_bytes_for_all_files_in_archived_snapshot"],
                "original 'file' count cannot be reinterpreted without matching exact byte values")
        corrected[prod]={
          "directory_entry_count_in_archived_source_snapshot_NOT_verified_ASDF_file_count":nums[0],
          "directory_total_bytes_including_any_auxiliary_entries":nums[1],
          "actual_XX_ASDF_superslab_count":None,
          "individual_superslab_names_sizes_and_GNU_POSIX_cksum":None,
          "live_Globus_provider_listing_verified":False
        }
    tampered=copy.deepcopy(lookup)
    tampered["halo_info"][0]-=1
    require(tampered["halo_info"]!=expected["halo_info"],
            "negative control must reject fabricated portal entry count")
    return {
     "date":"2026-09-28",
     "status":"E17D2B3B3R1_ERRATUM_ORIGINAL_35_COUNTS_ARE_UNFILTERED_DIRECTORY_ENTRIES_REAL_ASDF_COUNT_UNVERIFIED",
     "prospective_corrective_protocol_git_blob":P_BLOB,
     "preserved_old_report_sha256":p["prior_original_AUDIT_preserve"]["original_result_SHA256"],
     "preserved_old_manifest_git_blob":p["prior_original_AUDIT_preserve"]["original_archive_manifest_git_blob"],
     "original_run_id":p["prior_original_AUDIT_preserve"]["original_CI_run_id"],
     "upstream_portal_source_commit":src["commit"],
     "upstream_manifest_git_blob":src["portal_snapshot_blob"],
     "upstream_builder_git_blob":src["builder_blob"],
     "source_code_proven_directory_entry_semantics":semantics,
     "products_CORRECTED_as_directory_entries_NOT_ASDF_file_counts":corrected,
     "specific_real_ASDF_halo_info_file_count":None,
     "official_checksums_crc32_present_in_each_product_directory_per_docs_not_individually_verified":True,
     "no_claim_35_or_34_real_ASDFs_without_live_listing":True,
     "synthetic_extension_filter_tamper_rejected":True,
     "synthetic_35_to_34_source_count_tamper_rejected":True,
     "original_immutable_source_report_and_manifest_preserved":True,
     "external_real_Globus_directory_listing_or_ASDF_header_opened":False,
     "actual_real_halo_file_or_merger_tree_downloaded":False,
     "real_M200c_same_halo_mass_history_verified":False,
     "full_physical_finiteK_galaxy_bispectrum":"BLOCKED",
     "observed_odd_SEALED":True,
     "main_untouched_PR_draft":True
    }

def write(path,j):
    raw=(json.dumps(j,indent=2,allow_nan=False)+"\n").encode()
    path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():require(path.read_bytes()==raw,"old corrective output differs")
    else:
        with tempfile.NamedTemporaryFile(dir=path.parent,prefix=".e17d2b3b3r1_",delete=False) as f:
            t=Path(f.name);f.write(raw);f.flush();os.fsync(f.fileno())
        try:os.link(t,path)
        finally:t.unlink(missing_ok=True)
    print("E17D2B3B3R1_DIRECTORY_COUNT_ERRATUM_SHA256",sha(raw),flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--upstream-root",type=Path,default=Path("official-abacus-portal"))
    ap.add_argument("--output-dir",type=Path,default=OUT)
    a=ap.parse_args();p,old=check_original()
    if a.self_test:
        fake="du = [fn.stat().st_size for fn in fpath.iterdir()]\nnfile = len(du)\ndu = sum(du)\n"
        fake="def example(fpath):\n    "+fake.replace("\n","\n    ")
        require(builder_count_semantics(fake)["source_builder_AST_nfile_is_len_all_directory_entries"],
                "standalone AST source semantic self-test")
        print("E17D2B3B3R1_ORIGINAL_ARCHIVE_SHA_AND_AST_UNFILTERED_ITERDIR_PREFLIGHT_PASS",
              "NO_REAL_GLOBUS_OR_ASDF",flush=True)
        return
    out=run(a.upstream_root,p,old)
    write(a.output_dir/OUTNAME,out)
    print("E17D2B3B3R1_OFFICIAL_SOURCE_SEMANTICS_CORRECTED",
          "HALO_INFO_DIRECTORY_ENTRIES",out[
            "products_CORRECTED_as_directory_entries_NOT_ASDF_file_counts"][
                "halo_info"]["directory_entry_count_in_archived_source_snapshot_NOT_verified_ASDF_file_count"],
          "ACTUAL_ASDF_COUNT_UNKNOWN_REAL_HALO_UNREAD",flush=True)
if __name__=="__main__":main()
