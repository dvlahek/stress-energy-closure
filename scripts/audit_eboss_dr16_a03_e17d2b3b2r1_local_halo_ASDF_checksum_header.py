#!/usr/bin/env python3
"""E17D2b3b2r1: fail-closed local Abacus halo ASDF byte+header checker.

CI ONLY runs --self-test on its OWN temporary synthetic bytes. This module
contains NO downloader and does not access real Abacus data unless a user
separately supplies --file, --checksum-manifest and --official-extension.
Even success on a real file verifies LOCAL BYTES against SUPPLIED checksums
and copied header, NOT the external provider's attestation, halo M200c,
merger trees, original F+/- wake, or eBOSS galaxy bispectrum.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
from collections.abc import Mapping

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2r1_real_ASDF_byte_header_verifier_prereg_2026-09-28.json"
PRE_BLOB="b5e4bed15c9841bdfa8933a4d6960c9eda701e36"
B2DIR=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28"
B2=B2DIR/"e17d2b3b2_pinned_header_copy_l1_threshold_vs_200crit.json"
B2MAN=B2DIR/"archive_manifest.json"
B2PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_official_header_copy_mass_reference_prereg_2026-09-28.json"
B1=ROOT/"source_data/eboss_dr16_a03_e17d2b3b1_archived_CI_2026_09_28/e17d2b3b1_pinned_official_embedded_header_copy_state_z095.json"
E8=ROOT/"source_data/eboss_dr16_a03_e8_e11_archived_CI_2026_09_27/E8/frozen_matched_distributions_4000q.csv"
E16=ROOT/"source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json"
HALO_PATTERN=re.compile(r"halo_info_[0-9]{3}[.]asdf")
LINE=re.compile(r"^([0-9]+)\s+([0-9]+)\s+(.+)$")
SIM="AbacusSummit_base_c000_ph000"
REF_FIELDS=("Redshift","ScaleFactor","OmegaNow_m","SODensityL1",
            "ParticleMassHMsun","ParticleMassMsun")

def require(ok,reason):
    if not ok:raise ValueError("E17D2B3B2R1_FAIL_CLOSED: "+reason)
def sha(data):return hashlib.sha256(data).hexdigest()
def blob(data):
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()

def pinned():
    require(blob(PRE.read_bytes())==PRE_BLOB,"prospective verifier prereg Git blob drift")
    p=json.loads(PRE.read_bytes())
    parents=p["parents"]
    for f,key in ((B2,"E17D2B3B2_original_SHA256"),
                  (B1,"E17D2B3B1_original_SHA256"),
                  (E8,"E8_original_SHA256"),
                  (E16,"E16_original_SHA256")):
        require(sha(f.read_bytes())==parents[key],"frozen original SHA "+key)
    for f,key in ((B2MAN,"E17D2B3B2_manifest_git_blob"),
                  (B2PRE,"E17D2B3B2_protocol_git_blob")):
        require(blob(f.read_bytes())==parents[key],"frozen original Git blob "+key)
    reference=json.loads(B2.read_bytes())
    require(reference["observed_odd_SEALED"] is True
            and reference["actual_halo_info_ASDF_byte_checksum_header_verified"] is False
            and reference["same_object_true_M200c_200critical_remeasured"] is False
            and reference["full_physical_finiteK_galaxy_bispectrum"]=="BLOCKED"
            and all(p["absolute_STOP"].values()),
            "original source-only no-real-ASDF and odd seal")
    fields=reference["official_source_copy_six_exact_numeric_fields"]
    require(set(fields)==set(REF_FIELDS),"expected copied source header reference incomplete")
    require(reference["copy_L1_SO_threshold_over_200critical_SO_threshold"]>0,
            "frozen CompaSO L1 200crit threshold comparison missing")
    return p,fields

def path_gate(target,manifest):
    f=Path(target).absolute()
    m=Path(manifest).absolute()
    require(f.is_file() and not f.is_symlink(),"local candidate file missing or symlink")
    require(m.is_file() and not m.is_symlink(),"local checksum manifest missing or symlink")
    require(HALO_PATTERN.fullmatch(f.name) is not None,
            "candidate basename must be exact halo_info_NNN.asdf")
    require(f.parent.name=="halo_info" and f.parent.parent.name=="z0.950"
            and f.parent.parent.parent.name=="halos"
            and f.parent.parent.parent.parent.name==SIM,
            "candidate not nested in exact chosen simulation/halo epoch")
    require(m.parent==f.parent and m.name=="checksums.crc32",
            "checksum manifest must be supplied from same exact halo_info dir")
    require(0<f.stat().st_size<8*1024**3,
            "missing or unexpectedly huge individual superslab; no 6.6TB transfer")
    return f,m

def parse_manifest(raw,filename):
    try:lines=raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as e:raise ValueError("E17D2B3B2R1_FAIL_CLOSED: checksum list not ASCII/UTF8") from e
    selected=[]
    for line in lines:
        if not line.strip():continue
        m=LINE.fullmatch(line.strip())
        require(m is not None,"unknown provider checksum file line syntax: fail closed")
        label=m.group(3).strip()
        if label in (filename,"./"+filename):
            selected.append((int(m.group(1)),int(m.group(2))))
    require(len(selected)==1,"exactly one checksum entry for chosen actual filename required")
    return selected[0]

def checksums(f,m):
    expected=parse_manifest(m.read_bytes(),f.name)
    try:
        result=subprocess.run(["cksum",str(f)],capture_output=True,text=True,
                              check=True,timeout=240)
    except (OSError,subprocess.SubprocessError) as e:
        raise ValueError("E17D2B3B2R1_FAIL_CLOSED: GNU/POSIX cksum not available/successful") from e
    match=LINE.fullmatch(result.stdout.strip())
    require(match is not None,"GNU cksum returned unexpected syntax")
    actual=(int(match.group(1)),int(match.group(2)))
    require(actual==expected and actual[1]==f.stat().st_size,
            "GNU cksum CRC or byte count differs from supplied file manifest")
    with f.open("rb") as source:
        h=hashlib.sha256()
        for block in iter(lambda:source.read(1024*1024),b""):
            h.update(block)
    return {"file_size_bytes":actual[1],
            "GNU_POSIX_cksum_crc":actual[0],
            "file_local_computed_sha256":h.hexdigest(),
            "checksum_matches_SUPPLIED_manifest_only":True,
            "supplied_manifest_SHA256":sha(m.read_bytes()),
            "provider_attestation_of_supplied_manifest_independently_verified":False}

def number(v,name):
    if type(v) in (float,int):x=float(v)
    elif hasattr(v,"item"):
        y=v.item()
        require(type(y) in (int,float),"header nonnumeric numpy scalar "+name)
        x=float(y)
    else:
        raise ValueError("E17D2B3B2R1_FAIL_CLOSED: actual header field missing or nonnumeric "+name)
    require(math.isfinite(x),"actual ASDF header field nonfinite "+name)
    return x

def compare_header(header,reference):
    require(isinstance(header,Mapping),"actual ASDF tree header must be a mapping")
    require(header.get("SimName")==SIM,"actual ASDF SimName not selected c000 ph000")
    out={}
    for k in REF_FIELDS:
        require(k in header,"missing actual ASDF header key "+k)
        q=number(header[k],k);v=reference[k]
        tol=1e-10 if k in ("Redshift","ScaleFactor") else max(1e-10,2e-7*abs(v))
        require(abs(q-v)<=tol,"per-file actual header mismatch from official pinned copied state "+k)
        out[k]=q
    require(abs(1/out["ScaleFactor"]-1-out["Redshift"])<1e-10,
            "actual ASDF Redshift and ScaleFactor inconsistent")
    return out

def actual_header(f,extfile,reference,parents):
    import asdf
    extension=Path(extfile).absolute()
    require(extension.is_file() and not extension.is_symlink()
            and blob(extension.read_bytes())==
                parents["upstream_official_Blosc_extension_git_blob"],
            "official ASDF Blosc decoder not original pinned upstream source")
    spec=importlib.util.spec_from_file_location("e17d2b3b2r1_original_blosc",extension)
    require(spec and spec.loader,"official ASDF decoder cannot be loaded")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    # Lazy load only the header. Do not access halo columns or particle arrays.
    with asdf.open(f,extensions=[module.AbacusExtension()],
                   lazy_load=True,memmap=True) as af:
        require("header" in af.tree,"not an actual ASDF halo file with header")
        return compare_header(af.tree["header"],reference)

def synthetic_tests(p,reference):
    # Tiny CREATED-BY-THIS-TEST bytes, NEVER an actual provider halo_info ASDF.
    with tempfile.TemporaryDirectory(prefix="e17d2b3b2r1_SYNTHETIC_ONLY_") as t:
        root=Path(t)/SIM/"halos"/"z0.950"/"halo_info"
        root.mkdir(parents=True)
        f=root/"halo_info_002.asdf";f.write_bytes(b"FAKE ASDF SYNTHETIC CHECKSUM FIXTURE NOT A HALO")
        c=subprocess.run(["cksum",str(f)],capture_output=True,text=True,
                         check=True).stdout.split()
        m=root/"checksums.crc32"
        m.write_text(c[0]+" "+c[1]+" "+f.name+"\n")
        file_,manifest=path_gate(f,m)
        j=checksums(file_,manifest)
        require(j["file_local_computed_sha256"]==sha(f.read_bytes())
                and j["provider_attestation_of_supplied_manifest_independently_verified"] is False,
                "synthetic local SHA or provider-origin STOP failed")
        bad=[(c[0]+" "+str(int(c[1])+1)+" "+f.name+"\n","BYTECOUNT_CORRUPTED"),
             ("0 "+c[1]+" "+f.name+"\n","CRC_CORRUPTED"),
             (c[0]+" "+c[1]+" "+f.name+"\n"*2,"DUPLICATE")]
        negatives=[]
        for text,label in bad:
            if label=="DUPLICATE":text=(c[0]+" "+c[1]+" "+f.name+"\n")*2
            m.write_text(text)
            try:checksums(f,m)
            except ValueError:negatives.append(label)
            else:raise AssertionError("accepted synthetic altered "+label)
        m.write_text(c[0]+" "+c[1]+" "+f.name+"\n")
        wrong=dict(reference);wrong["Redshift"]=.95
        try:compare_header({"SimName":SIM,**wrong},reference)
        except ValueError:negatives.append("DIRECTORY_Z_AS_REAL_HEADER")
        else:raise AssertionError("accepted guessed directory redshift")
        try:compare_header({"SimName":"OTHER",**reference},reference)
        except ValueError:negatives.append("WRONG_SIMNAME")
        else:raise AssertionError("accepted wrong simulation identity")
        f_link=root/"halo_info_003.asdf";f_link.symlink_to(f)
        try:path_gate(f_link,m)
        except ValueError:negatives.append("SYMLINK_REFUSED")
        else:raise AssertionError("accepted synthetic file symlink")
        require(len(negatives)==6,"missing prereg negative controls")
        print("E17D2B3B2R1_SYNTHETIC_GNU_CKSUM_SHA_HEADER_AND_6_NEGATIVE_QA_PASS",
              "NO_REAL_ABACUS_ASDF_READ",flush=True)
        return negatives

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    mode=ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--self-test",action="store_true")
    mode.add_argument("--file",type=Path)
    ap.add_argument("--checksum-manifest",type=Path)
    ap.add_argument("--official-extension",type=Path)
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    p,reference=pinned()
    if a.self_test:
        negatives=synthetic_tests(p,reference)
        print("E17D2B3B2R1_SELFTEST_CERTIFIES_CODE_ONLY_REAL_DATA_BLOCKED",
              len(negatives),flush=True)
        return
    require(a.file and a.checksum_manifest and a.official_extension and a.output,
            "real-file mode requires exact local file, supplied provider manifest, pinned official decoder and a LOCAL output")
    f,m=path_gate(a.file,a.checksum_manifest)
    file_q=checksums(f,m)
    header=actual_header(f,a.official_extension,reference,p["parents"])
    # Provenance is never authenticated purely by a supplied local checksum
    # file. This report stays local. It MUST NOT be auto-pushed to GitHub.
    result={
      "stage":"E17D2B3B2R1",
      "status":"LOCAL_BYTES_MATCH_SUPPLIED_CKSUM_AND_OFFICIAL_COPIED_HEADER_PROVIDER_MANIFEST_ORIGIN_UNATTESTED",
      "file_basename":f.name,
      "simname":SIM,"directory_label_only":"z0.950",
      "source_copy_reference_original_sha256":p["parents"]["E17D2B3B2_original_SHA256"],
      "actual_local_file_byte_integrity":file_q,
      "local_ASDF_header_fields_matching_pinned_source_copy":header,
      "independent_provider_manifest_origin_authenticated":False,
      "provider_item_uri_verified":False,
      "actual_M200c_200critical_remeasured":False,
      "merger_tree_progenitor_branches_verified":False,
      "original_E8_Fplus_Fminus_neutrino_halo_wake_verified":False,
      "observed_odd_SEALED":True,
      "full_physical_finiteK_galaxy_bispectrum":"BLOCKED"}
    raw=(json.dumps(result,indent=2,allow_nan=False)+"\n").encode()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    require(not a.output.exists(),"local report already exists; refuse overwrite")
    a.output.write_bytes(raw)
    print("E17D2B3B2R1_LOCAL_BYTES_SUPPLIED_CHECKSUM_AND_COPIED_HEADER_MATCH",
          "REPORT_SHA256",sha(raw),
          "PROVIDER_MANIFEST_ORIGIN_UNATTESTED_M200C_AND_GALAXY_B_BLOCKED",flush=True)
if __name__=="__main__":main()
