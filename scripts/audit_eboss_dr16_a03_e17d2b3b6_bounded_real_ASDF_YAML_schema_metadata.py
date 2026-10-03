#!/usr/bin/env python3
"""E17D2b3b6: ONLY ASDF YAML descriptors; never read/decompress real halo arrays."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import yaml
from yaml.nodes import MappingNode, ScalarNode, SequenceNode

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b6_official_schema_real_ASDF_YAML_descriptors_only_prereg_2026-09-28.json"
PARENT=ROOT/"source_data/eboss_dr16_a03_e17d2b3b5_user_reported_real_halo_000_byte_and_dual_header_local_acceptance_2026-09-28.json"
PRE_BLOB="e1f43da1a73aad3ef85c76a114fbffadc76d4761"
PARENT_BLOB="c0ac072263f3d8b19d01caae091cb743c22a4037"
LIMIT=262144
SIM="AbacusSummit_base_c000_ph000"
MANDATORY=("id","N","SO_radius","SO_central_particle")
OPTIONAL=("SO_central_density","x_L2com","r100_L2com","r50_L2com_i16")
FIELDS=MANDATORY+OPTIONAL

def need(ok,reason):
    if not ok: raise ValueError("E17D2B3B6_FAIL_CLOSED: "+reason)

def sha(b): return hashlib.sha256(b).hexdigest()
def blob(b): return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def pinned():
    raw=PRE.read_bytes()
    need(blob(raw)==PRE_BLOB,"prereg Git blob drift")
    p=json.loads(raw)
    need(all(p["absolute_STOP"].values()),"absolute STOP gate")
    previous=PARENT.read_bytes()
    need(blob(previous)==PARENT_BLOB,"prior evidence Git blob drift")
    return p,json.loads(previous)

def mapping(n,label):
    need(isinstance(n,MappingNode),label+" is not mapping")
    out={}
    for k,v in n.value:
        need(isinstance(k,ScalarNode),label+" has non-scalar key")
        need(k.value not in out,label+" duplicate "+k.value)
        out[k.value]=v
    return out

def scalar(n,label):
    need(isinstance(n,ScalarNode),label+" is not scalar")
    return n.value

def prefix(f):
    buf=bytearray()
    with f.open("rb") as h:
        while len(buf)<LIMIT:
            line=h.readline(LIMIT-len(buf)+1)
            need(bool(line),"YAML terminator absent")
            need(len(buf)+len(line)<=LIMIT,"YAML > 256 KiB")
            buf.extend(line)
            if line.strip()==b"...":return bytes(buf)
    raise ValueError("E17D2B3B6_FAIL_CLOSED: YAML > 256 KiB")

def parse(raw):
    need(raw.startswith(b"#ASDF ") and len(raw)<=LIMIT,"ASDF prefix/byte cap")
    need(raw.splitlines()[-1].strip()==b"...","no terminal YAML line")
    doc=yaml.compose(raw.decode("utf-8"),Loader=yaml.SafeLoader)
    top=mapping(doc,"top")
    need("header" in top and "data" in top,"header/data missing")
    header=mapping(top["header"],"header")
    need(scalar(header.get("SimName"),"SimName")==SIM,"wrong sim")
    need(abs(float(scalar(header.get("Redshift"),"Redshift"))-0.952838237036305)<=1e-10,
         "wrong redshift")
    data=mapping(top["data"],"data")
    need(all(n in data for n in MANDATORY),"mandatory fields absent")
    selected={}
    for name in FIELDS:
        if name not in data:continue
        desc=mapping(data[name],name)
        need("shape" in desc and "datatype" in desc,"descriptor missing shape/datatype "+name)
        sh=desc["shape"]
        need(isinstance(sh,SequenceNode),"shape not sequence "+name)
        dims=[int(scalar(v,"shape dim")) for v in sh.value]
        need(1<=len(dims)<=4 and all(0<v<10**10 for v in dims),"invalid shape "+name)
        dtype=desc["datatype"]
        if isinstance(dtype,ScalarNode):d=dtype.value
        elif isinstance(dtype,SequenceNode):d="structured_dtype_metadata"
        else:raise ValueError("E17D2B3B6_FAIL_CLOSED: invalid datatype "+name)
        need(bool(d),"empty datatype "+name)
        item={"shape_METADATA_ONLY":dims,"datatype_METADATA_ONLY":d}
        for k in ("byteorder","source"):
            if k in desc:
                v=scalar(desc[k],k)
                need(k!="source" or v.isdecimal(),"invalid source block index")
                item[k+"_METADATA_ONLY"]=int(v) if k=="source" else v
        selected[name]=item
    return {
        "YAML_prefix_bytes":len(raw),
        "YAML_prefix_SHA256":sha(raw),
        "raw_column_count_NAMES_ONLY":len(data),
        "raw_column_names_NAMES_ONLY":sorted(data),
        "preselected_descriptors_METADATA_ONLY":selected,
        "halo_row_values_accessed":False,
        "binary_block_bytes_read":0,
        "M200c_measured":False
    }

def selftest():
    header=("#ASDF 1.0.0\n#ASDF_STANDARD 1.5.0\n%YAML 1.1\n"
            "%TAG ! tag:stsci.edu:asdf/\n--- !core/asdf-1.1.0\n"
            "header:\n  SimName: "+SIM+"\n  Redshift: 0.952838237036305\n"
            "data:\n")
    def col(n,i):
        return "  "+n+": !core/ndarray-1.0.0\n    shape: [12]\n    datatype: uint32\n    source: "+str(i)+"\n"
    good=(header+"".join(col(n,i) for i,n in enumerate(MANDATORY))+"...\n").encode()
    with tempfile.TemporaryDirectory(prefix="e17d2b3b6_synthetic_") as t:
        f=Path(t)/"fake.asdf"
        f.write_bytes(good+b"SYNTHETIC_BINARY_NOT_A_REAL_HALO")
        r=parse(prefix(f))
        need(r["raw_column_count_NAMES_ONLY"]==4 and r["binary_block_bytes_read"]==0,
             "positive synthetic")
    cases={
      "MISSING_N":(header+"".join(col(n,i) for i,n in enumerate(MANDATORY) if n!="N")+"...\n").encode(),
      "DUPLICATE_ID":good.replace(b"...\n",col("id",9).encode()+b"...\n"),
      "WRONG_SIM":good.replace(SIM.encode(),b"OTHER_SIM",1),
      "MISSING_DTYPE":good.replace(b"    datatype: uint32\n",b"",1),
      "NO_END":good.replace(b"...\n",b"",1),
      "OVER_LIMIT":good+b"#"*(LIMIT+1)
    }
    for k,v in cases.items():
        try:parse(v)
        except (ValueError,TypeError,UnicodeError,yaml.YAMLError):print("NEGATIVE_REJECT",k,flush=True)
        else:raise AssertionError("synthetic negative accepted "+k)
    print("E17D2B3B6_SYNTHETIC_SCHEMA_AND_SIX_NEGATIVES_PASS",flush=True)

def actual(f,m,out,p,previous):
    need(f.is_file() and not f.is_symlink(),"file absent or symlink")
    need(m.is_file() and not m.is_symlink(),"manifest absent or symlink")
    need(f.name=="halo_info_000.asdf" and m.name=="checksums.crc32" and f.parent==m.parent,
         "exact file/manifest path")
    need(f.parent.name=="halo_info" and f.parent.parent.name=="z0.950"
         and f.parent.parent.parent.name=="halos"
         and f.parent.parent.parent.parent.name==SIM,"wrong simulation/epoch path")
    need(f.stat().st_size==2223483833,"file size drift")
    need(sha(m.read_bytes())==p["only_allowed_actual_file"]["same_directory_manifest_sha256"],
         "manifest SHA drift")
    need(not out.exists(),"refuse report overwrite")
    report=parse(prefix(f))
    report.update({
      "stage":p["stage"],
      "source":"USER_LOCAL_ACTUAL_ASDF_YAML_TEXT_DESCRIPTORS_ONLY",
      "pre_registered_protocol_git_blob":PRE_BLOB,
      "previous_user_stdout_evidence_git_blob":PARENT_BLOB,
      "previous_full_file_SHA256_USER_REPORTED_NOT_REHASHED_HERE":
          previous["local_user_reported_byte_preflight"]["reported_computed_full_ASDF_SHA256"],
      "source_provider_independently_attested":False,
      "full_ASDF_binary_file_not_loaded":True,
      "no_halo_arrays_decompressed":True,
      "observed_odd_SEALED":True,
      "M200c_and_physical_neutrino_wake_B":"BLOCKED"
    })
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_bytes((json.dumps(report,indent=2,allow_nan=False)+"\n").encode())
    print("E17D2B3B6_REAL_ASDF_YAML_METADATA_ONLY_PASS",flush=True)
    print("YAML_BYTES",report["YAML_prefix_bytes"],flush=True)
    print("RAW_COLUMN_COUNT",report["raw_column_count_NAMES_ONLY"],flush=True)
    for k,v in report["preselected_descriptors_METADATA_ONLY"].items():
        print("SCHEMA",k,v,flush=True)
    print("REPORT",out,"SHA256",sha(out.read_bytes()),flush=True)
    print("NO_BINARY_BLOCK_NO_HALO_ROWS_OBSERVED_ODD_SEALED",flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    g=a.add_mutually_exclusive_group(required=True)
    g.add_argument("--self-test",action="store_true")
    g.add_argument("--file",type=Path)
    a.add_argument("--checksum-manifest",type=Path)
    a.add_argument("--output",type=Path)
    args=a.parse_args()
    p,previous=pinned()
    if args.self_test:selftest()
    else:
        need(args.checksum_manifest is not None and args.output is not None,
             "file mode requires checksum manifest and local output")
        actual(args.file,args.checksum_manifest,args.output,p,previous)

if __name__=="__main__":main()
