#!/usr/bin/env python3
"""E17D2b3b7: prospectively locked first real L1 id/N columns, never observed eBOSS."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import sys
from collections.abc import Mapping

ROOT = Path(__file__).resolve().parents[1]
SIM = "AbacusSummit_base_c000_ph000"
PRE = ROOT / "source_data/eboss_dr16_a03_e17d2b3b7_first_real_L1_id_N_columns_prereg_2026-09-29.json"
EVIDENCE = ROOT / "source_data/eboss_dr16_a03_e17d2b3b6_user_reported_real_ASDF_YAML_schema_local_acceptance_2026-09-29.json"
B6_SCRIPT = ROOT / "scripts/audit_eboss_dr16_a03_e17d2b3b6_bounded_real_ASDF_YAML_schema_metadata.py"
PRE_BLOB = "296bc22e3d38a414017cc018e1f8ba8658899278"
EVIDENCE_BLOB = "570366acf3c9db90a3d7542b2ce80bc21195a6ba"
B6_SCRIPT_BLOB = "718d2dde75f29dfbf2e32ab54cb624fc2a6688f8"
DECODER_BLOB = "8c5e5135736409bb0fab77bff06ad1259951189d"
ROWS = 11676687
SELECTED = (0, 5838343, 11676686)
MEMORY_LIMIT = 6442450944
MAX_RAW_ARRAY_BYTES = 268435456
MASS_REF = 2109081520.453063

def need(ok, why):
    if not ok:
        raise ValueError("E17D2B3B7_FAIL_CLOSED: " + why)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def gitblob(raw):
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()

def pins():
    need(gitblob(PRE.read_bytes()) == PRE_BLOB, "prospective B7 Git blob drift")
    p = json.loads(PRE.read_bytes())
    need(gitblob(EVIDENCE.read_bytes()) == EVIDENCE_BLOB, "B6 user-reported evidence archive drift")
    need(gitblob(B6_SCRIPT.read_bytes()) == B6_SCRIPT_BLOB, "original B6 YAML reader drift")
    need(all(p["absolute_STOP"].values()), "absolute STOP prereg")
    need(p["fields_prospectively_locked"]["source_binary_blocks_exact"] == [0, 7], "binary blocks modified")
    need(p["fields_prospectively_locked"]["allowed_raw_fields_exact"] == ["id", "N"], "field selector modified")
    need(p["selector"]["rows_for_individual_reporting"] == list(SELECTED), "preselected row selector drift")
    need(p["identity"]["expected_rows"] == ROWS, "row-count parent drift")
    return p

def check_schema(report):
    need(report.get("YAML_prefix_bytes") == 14416, "B6 YAML length changed")
    need(report.get("raw_column_count_NAMES_ONLY") == 70, "B6 column count changed")
    need(report.get("binary_block_bytes_read") == 0, "B6 unexpectedly read binary blocks")
    need(report.get("halo_row_values_accessed") is False, "B6 real rows were already accessed")
    desc = report.get("preselected_descriptors_METADATA_ONLY")
    need(isinstance(desc, dict), "missing B6 descriptors")
    for field, dtype, source, shape in (
        ("id", "uint64", 0, [ROWS]), ("N", "uint32", 7, [ROWS])
    ):
        d = desc.get(field)
        need(isinstance(d, dict), "missing B6 descriptor " + field)
        need(d.get("shape_METADATA_ONLY") == shape, "B6 shape mismatch " + field)
        need(d.get("datatype_METADATA_ONLY") == dtype, "B6 type mismatch " + field)
        need(d.get("byteorder_METADATA_ONLY") == "little", "B6 byteorder mismatch " + field)
        need(d.get("source_METADATA_ONLY") == source, "B6 binary source mismatch " + field)

def check_report_bytes(raw, expected):
    need(sha(raw) == expected, "user-local B6 report SHA256 drift")

def report_preflight(path):
    need(path.is_file() and not path.is_symlink(), "B6 local report missing or symlink")
    raw = path.read_bytes()
    need(len(raw) < 65536, "B6 report unexpectedly large")
    check_report_bytes(raw, "70820fb64965849ec692a035d4c2d6a1dedf38a8de492f19c88a4082dd42cedd")
    j = json.loads(raw)
    check_schema(j)
    need(j.get("M200c_measured") is False and j.get("observed_odd_SEALED") is True,
         "B6 physical/observational STOP invalid")
    return j

def input_preflight(f, m, b6, output, p):
    need(f.is_file() and not f.is_symlink(), "actual ASDF missing or symlink")
    need(m.is_file() and not m.is_symlink(), "checksum manifest missing or symlink")
    need(f.name == "halo_info_000.asdf" and m.name == "checksums.crc32"
         and f.parent == m.parent, "wrong exact local file/manifest")
    need(f.parent.name == "halo_info"
         and f.parent.parent.name == "z0.950"
         and f.parent.parent.parent.name == "halos"
         and f.parent.parent.parent.parent.name == SIM, "wrong simulation/epoch path")
    need(f.stat().st_size == p["identity"]["expected_file_bytes"], "ASDF byte-size drift")
    need(sha(m.read_bytes()) == p["identity"]["expected_manifest_SHA256"], "manifest SHA drift")
    expected_b6_dir = ROOT / "eboss_workspace/e17d2b3b6_reports"
    need(b6 == expected_b6_dir / "halo_info_000_schema_only.json", "B6 report path drift")
    need(output == ROOT / "eboss_workspace/e17d2b3b7_reports/halo_info_000_id_N_locked_rows.json",
         "output path drift")
    need(not output.exists(), "do not overwrite any earlier local output")

def limit_memory():
    need(sys.platform.startswith("linux"), "this protocol requires Linux WSL resource cap")
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    need(hard == resource.RLIM_INFINITY or hard >= MEMORY_LIMIT,
         "current hard RLIMIT_AS lower than fixed protocol")
    resource.setrlimit(resource.RLIMIT_AS, (MEMORY_LIMIT, MEMORY_LIMIT))
    need(resource.getrlimit(resource.RLIMIT_AS)[0] == MEMORY_LIMIT, "memory cap failed")

def read_b6_prefix_again(f, b6):
    spec = importlib.util.spec_from_file_location("b6_pinned_metadata_only", B6_SCRIPT)
    need(spec is not None and spec.loader is not None, "original B6 script unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    now = module.parse(module.prefix(f))
    for name in ("YAML_prefix_SHA256", "YAML_prefix_bytes", "raw_column_count_NAMES_ONLY",
                 "preselected_descriptors_METADATA_ONLY"):
        need(now.get(name) == b6.get(name), "B6/current ASDF prefix mismatch " + name)
    need(now.get("binary_block_bytes_read") == 0, "prefix probe accessed binary data")

def pinned_extension(path):
    need(path.is_file() and not path.is_symlink(), "official Abacus ASDF extension absent/symlink")
    need(gitblob(path.read_bytes()) == DECODER_BLOB, "official upstream blsc source Git blob drift")
    import asdf
    spec = importlib.util.spec_from_file_location("e17d2b3b7_official_abacus_blosc", path)
    need(spec is not None and spec.loader is not None, "cannot load pinned blsc source")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    asdf.get_config().add_extension(mod.AbacusExtension())
    return asdf

def evaluate_arrays(ids, n, mass=MASS_REF, indexes=SELECTED, expected_rows=ROWS, synthetic_fixture=False):
    import numpy as np
    need(sys.byteorder == "little", "unexpected machine byteorder")
    need(expected_rows == ROWS or (synthetic_fixture and expected_rows == 3),
         "unexpected row contract")
    need(getattr(ids, "shape", None) == (expected_rows,), "unexpected id array shape")
    need(getattr(n, "shape", None) == (expected_rows,), "unexpected N array shape")
    for name, arr, kind, size in (("id", ids, "u", 8), ("N", n, "u", 4)):
        dt = np.dtype(arr.dtype)
        need(dt.kind == kind and dt.itemsize == size and dt.byteorder in ("<", "="),
             "unexpected raw dtype " + name)
        need(arr.nbytes == expected_rows * size, "unexpected materialized byte count " + name)
    need(ids.nbytes + n.nbytes <= MAX_RAW_ARRAY_BYTES,
         "two columns exceed preregistered raw-array memory bound")
    need(len(indexes) == 3 and
         tuple(indexes) == ((0, 1, 2) if synthetic_fixture else SELECTED),
         "row index selector was altered")
    need(math.isfinite(float(mass)) and abs(float(mass) - MASS_REF) <= 1e-7,
         "unapproved mass reference")
    records = []
    for i in indexes:
        count = int(n[i])
        records.append({"row_index":i, "raw_id":int(ids[i]), "raw_N":count,
                        "L1_assigned_mass_Msun_per_h":count * MASS_REF})
    need(len({r["raw_id"] for r in records}) == len(records),
         "selected raw halo ids are not unique")
    total = int(np.sum(n, dtype=np.int64))
    return {
        "rows_in_exact_superslab_000":expected_rows,
        "minimum_N_full_superslab":int(np.min(n)),
        "maximum_N_full_superslab":int(np.max(n)),
        "sum_N_full_superslab":total,
        "zero_N_count_full_superslab_NO_CUT":int(np.count_nonzero(n == 0)),
        "selected_fixed_rows":records,
        "raw_columns_loaded_exact":["id", "N"],
        "aggregate_is_only_selected_superslab_000":True,
        "M200c_measured":False,
        "merger_tree_measured":False,
        "physical_neutrino_wake_or_galaxy_B_measured":False,
        "observed_odd_SEALED":True
    }

def selftest(p):
    import copy
    import numpy as np
    fake = {
        "YAML_prefix_bytes":14416, "raw_column_count_NAMES_ONLY":70,
        "binary_block_bytes_read":0,"halo_row_values_accessed":False,
        "preselected_descriptors_METADATA_ONLY":{
            "id":{"shape_METADATA_ONLY":[ROWS],"datatype_METADATA_ONLY":"uint64",
                  "byteorder_METADATA_ONLY":"little","source_METADATA_ONLY":0},
            "N":{"shape_METADATA_ONLY":[ROWS],"datatype_METADATA_ONLY":"uint32",
                 "byteorder_METADATA_ONLY":"little","source_METADATA_ONLY":7}
        }
    }
    check_schema(fake)
    ids = np.array([3, 5, 9], dtype="<u8")
    n = np.array([36, 40, 80], dtype="<u4")
    synthetic_result = evaluate_arrays(ids, n, indexes=(0,1,2),
                                       expected_rows=3, synthetic_fixture=True)
    need(synthetic_result["sum_N_full_superslab"] == 156
         and [q["L1_assigned_mass_Msun_per_h"] for q in synthetic_result["selected_fixed_rows"]]
             == [36*MASS_REF,40*MASS_REF,80*MASS_REF]
         and synthetic_result["minimum_N_full_superslab"] == 36
         and synthetic_result["maximum_N_full_superslab"] == 80,
         "same actual evaluation function synthetic numeric happy-path")
    def reject(label, action):
        try:action()
        except (ValueError, TypeError, KeyError):
            print("NEGATIVE_REJECT",label,flush=True)
        else:raise AssertionError("negative control accepted: " + label)
    def changed(field,key,value):
        j=copy.deepcopy(fake);j["preselected_descriptors_METADATA_ONLY"][field][key]=value
        return j
    reject("WRONG_ROW_COUNT",lambda:check_schema(changed("N","shape_METADATA_ONLY",[ROWS+1])))
    reject("WRONG_DTYPE",lambda:check_schema(changed("id","datatype_METADATA_ONLY","int64")))
    reject("WRONG_SOURCE_INDEX",lambda:check_schema(changed("N","source_METADATA_ONLY",8)))
    reject("B6_REPORT_HASH_DRIFT",lambda:check_report_bytes(b"tampered",p["identity"]["local_b6_report_required_exact_SHA256"]))
    reject("MASS_REFERENCE_DRIFT",lambda:need(abs(MASS_REF + 10 - MASS_REF)<=1e-7,"mass mismatch"))
    class OnlyAllowed:
        def __init__(self):self.data={"id":ids,"N":n}
        def __getitem__(self,k):
            need(k in ("id","N"),"unregistered extra binary column access")
            return self.data[k]
    only=OnlyAllowed()
    need(only["id"] is ids and only["N"] is n,"synthetic whitelist")
    reject("EXTRA_BINARY_COLUMN",lambda:only["SO_radius"])
    reject("OUTPUT_ALREADY_EXISTS",lambda:need(not True,"no overwrite"))
    reject("UNAPPROVED_ROW_SELECTOR",lambda:need([0,1,2]==list(SELECTED),"selector changed"))
    print("E17D2B3B7_SYNTHETIC_ID_N_AND_EIGHT_NEGATIVES_PASS",flush=True)
    print("E17D2B3B7_NO_REAL_ASDF_NO_OBSERVED_ODD_CI_ONLY",flush=True)

def actual(args, p):
    import numpy as np
    f = args.file.absolute()
    manifest = args.checksum_manifest.absolute()
    local_b6 = args.b6_report.absolute()
    output = args.output.absolute()
    input_preflight(f, manifest, local_b6, output, p)
    b6 = report_preflight(local_b6)
    limit_memory()
    read_b6_prefix_again(f, b6)  # bounded YAML ONLY, before any binary block
    asdf = pinned_extension(args.official_extension.absolute())
    with asdf.open(f, lazy_load=True, memmap=True) as af:
        header=af.tree["header"]
        need(header.get("SimName")==SIM, "ASDF header simulation drift")
        need(abs(float(header["Redshift"])-0.952838237036305)<1e-10, "ASDF epoch drift")
        need(abs(float(header["SODensityL1"])-228.52306365966797)<1e-7,
             "ASDF L1 SO threshold drift")
        need(abs(float(header["ParticleMassHMsun"])-MASS_REF)<=1e-7,
             "ASDF particle mass reference drift")
        data=af.tree["data"]
        # ONLY these two named columns. The ASDF format stores them as
        # separate binary blocks; no SO/radius/particle array is touched.
        ids=np.asarray(data["id"])
        n=np.asarray(data["N"])
        result=evaluate_arrays(ids,n,float(header["ParticleMassHMsun"]))
    result.update({
        "stage":p["stage"],"status":"USER_LOCAL_FIRST_REAL_TWO_L1_COLUMNS_ONLY",
        "b7_prereg_git_blob":PRE_BLOB,
        "b6_user_local_report_sha256_verified":p["identity"]["local_b6_report_required_exact_SHA256"],
        "supplied_manifest_origin_provider_signed":False,
        "original_real_ASDF_full_sha256_recomputed_here":False,
        "source_file":str(f),
        "no_other_ASDF_columns_or_particle_files_accessed":True,
        "no_new_science_cuts_or_eBOSS_unblinding":True
    })
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("x",encoding="utf-8") as h:
        json.dump(result,h,indent=2,allow_nan=False)
        h.write("\n")
        h.flush()
        os.fsync(h.fileno())
    print("E17D2B3B7_LOCAL_REAL_ID_N_ONLY_PASS",flush=True)
    print("SUM_N",result["sum_N_full_superslab"],"MIN_N",result["minimum_N_full_superslab"],
          "MAX_N",result["maximum_N_full_superslab"],flush=True)
    for r in result["selected_fixed_rows"]:
        print("LOCKED_ROW",json.dumps(r,sort_keys=True),flush=True)
    print("REPORT",output,"SHA256",sha(output.read_bytes()),flush=True)
    print("ONLY_L1_CATALOG_MASS_NOT_M200C_NO_WAKE_NO_ODD",flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--self-test",action="store_true")
    group.add_argument("--file",type=Path)
    parser.add_argument("--checksum-manifest",type=Path)
    parser.add_argument("--b6-report",type=Path)
    parser.add_argument("--official-extension",type=Path)
    parser.add_argument("--output",type=Path)
    a=parser.parse_args()
    p=pins()
    if a.self_test:selftest(p)
    else:
        need(all(x is not None for x in
                 (a.checksum_manifest,a.b6_report,a.official_extension,a.output)),
             "actual mode requires original manifest, B6 report, pinned decoder, output")
        actual(a,p)
if __name__=="__main__":
    main()
