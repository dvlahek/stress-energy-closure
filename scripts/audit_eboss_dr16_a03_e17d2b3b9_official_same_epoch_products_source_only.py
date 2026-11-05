#!/usr/bin/env python3
"""E17D2b3b9: exact Git-pinned official TEXT-only same-epoch M200c availability.

CI reads five sparse-checkout official text sources and two existing tiny
tracked audit JSONs. It MUST NOT read any real ASDF, particle, merger tree,
observed galaxy data, or user-local WSL file. "Not identified" is scoped
to the exact pinned official documented products, not all possible providers.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b9_pinned_same_epoch_M200c_public_product_feasibility_prereg_2026-09-29.json"
B8=ROOT/"source_data/eboss_dr16_a03_e17d2b3b8_official_source_only_N35_M200c_feasibility_QA_result_2026-09-29.json"
B7=ROOT/"source_data/eboss_dr16_a03_e17d2b3b7r1_user_uploaded_local_id_N_result_and_existing_output_failure_diagnostic_2026-09-29.json"
PINS=(("prereg",PRE,"8ed61aa3927a4ab89cf423c424aa9c24563637a9"),
      ("B8",B8,"1c1de601c657af7c20590a86945d2cb39abe72e9"),
      ("B7",B7,"f27ed284c95e48548e50336da04c1a1d4b137186"))
OFFICIAL=(
    ("simulation_document",ROOT/"official-abacussummit/docs/simulations.rst","b5291009d2bb84dfae0c95da2a9734300841f84f"),
    ("data_access_document",ROOT/"official-abacussummit/docs/data-access.rst","6aad7833ab4482ebdd7cd2a23a3011f25e16e329"),
    ("data_products_document",ROOT/"official-abacussummit/docs/data-products.rst","f7c847b174365744b23cb6d1ebaeb010a5bd6ca7"),
    ("simulation_table",ROOT/"official-abacussummit/Simulations/simulations.csv","53d32e3fbec68ab64da1d52a3bf904fc3eb06c97"),
    ("loader",ROOT/"official-abacusutils/abacusnbody/data/compaso_halo_catalog.py","adb16aee1cbac863db5301c6e380938ee7f76547"),
)
SIM="AbacusSummit_base_c000_ph000"
SECONDARY_Z=0.952838237036305
EXPECTED_FULL=(3.0,2.5,2.0,1.7,1.4,1.1,0.8,0.5,0.4,0.3,0.2,0.1)

def require(ok,message):
    if not ok:raise ValueError("E17D2B3B9_SOURCE_ONLY_FAIL_CLOSED: "+message)

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def pinned_parent():
    out={}
    for label,path,sha in PINS:
        require(path.is_file() and not path.is_symlink() and blob(path.read_bytes())==sha,
                "original pinned "+label+" Git blob changed")
        out[label]=json.loads(path.read_bytes())
    pre,b8,b7=out["prereg"],out["B8"],out["B7"]
    require(pre["status"]=="SOURCE_ONLY_PREREG_BEFORE_B9_CI_NO_REAL_HALO_IO"
            and all(pre["absolute_STOP"].values()),"original B9 prereg STOP drift")
    require(b8["status"]=="PINNED_SOURCE_AND_EIGHT_SYNTHETIC_NEGATIVES_CI_PASS_NO_NEW_REAL_DATA"
            and b8["observed_odd_SEALED"],"B8 original source-only outcome drift")
    require(b7["data_scope"]["real_redshift_from_prev_metadata"]==SECONDARY_Z
            and b7["limitations"]["halo_M200c_measured"] is False
            and b7["limitations"]["whole_34_file_catalogue_sampled"] is False,
            "user reported single-epoch evidence promoted or drifted")
    return pre,b8,b7

def pinned_sources(pre,args):
    supplied=(args.simulations,args.access,args.products,args.table,args.loader)
    expect=tuple(item[1] for item in OFFICIAL)
    require(tuple(p.resolve() for p in supplied)==tuple(p.resolve() for p in expect),
            "only five exact sparse-checkout TEXT source paths allowed")
    out={}
    for (key,path,expected),inp in zip(OFFICIAL,supplied):
        require(inp.is_file() and not inp.is_symlink() and inp.suffix in (".rst",".csv",".py"),
                "official source not regular text or missing "+key)
        raw=inp.read_bytes()
        require(blob(raw)==expected,"original upstream Git blob changed "+key)
        out[key]=raw.decode("utf-8")
    p=pre["source_pins"]
    require(p["repository"]=="abacusorg/AbacusSummit"
            and p["commit"]=="4b1959c710cb0c49aa305c6213a228aa2a4587ff"
            and p["loader_repository"]=="abacusorg/abacusutils"
            and p["loader_commit"]=="24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6",
            "upstream repository/commit Git pin drift")
    return out

def documented_product_checks(s):
    sim,access,prod,table,loader=(s[k] for k in
        ("simulation_document","data_access_document","data_products_document","simulation_table","loader"))
    match=re.search(r"(?m)^\* \*\*Full\*\*: \*z\* = ([\d., ]+)$",sim)
    require(match is not None,"official Full redshift source field absent")
    full=tuple(float(x.strip()) for x in match.group(1).split(","))
    rows=list(csv.reader(table.splitlines(),skipinitialspace=True))
    require(len(rows)>1 and "Full Outputs" in rows[0][6],"official CSV header drift")
    selected=[r for r in rows[1:] if r and r[0].strip()==SIM]
    require(len(selected)==1,"exact c000 ph000 simulation table identity missing/duplicate")
    table_mode=selected[0][6].strip()
    schema_segments={}
    for a,b in (("clean_dt = np.dtype(","clean_dt_progen = np.dtype("),
                ("clean_dt_progen = np.dtype(","halo_lc_dt = np.dtype("),
                ("user_dt = np.dtype(","align=True,")):
        require(a in loader and b in loader.split(a,1)[1],"official loader dtype section absent "+a)
        region=loader.split(a,1)[1].split(b,1)[0]
        schema_segments[a]=re.findall(r"\('([^']+)',\s*np\.\w+",region)
    schema_fields=set(q for fields in schema_segments.values() for q in fields)
    checks={
        "EXACT_C000_PH000_OFFICIAL_TABLE_IS_FULL":table_mode=="Full",
        "OFFICIAL_FULL_PRIMARY_EPOCHS_EXACT":full==EXPECTED_FULL,
        "SECONDARY_EXACT_Z_ABSENT_FROM_FULL_SNAPSHOTS":
            all(abs(SECONDARY_Z-z)>1e-3 for z in full) and 0.95 not in full,
        "NEIGHBORING_FULL_Z_1p1_AND_0p8_ONLY_AS_DIFFERENT_EPOCHS":1.1 in full and 0.8 in full,
        "SECONDARY_STANDARD_PRODUCTS_HALO_AND_PID_NOT_RV_FIELD":
            bool(re.search(r"For the 21 secondary redshifts, we output the halo catalogs and the halo\s+subsample particle IDs.*?only, so not the\s+positions/velocities nor the field particles",prod,re.S)),
        "STANDARD_NERSC_DISK_EXCLUDES_100PCT_FULL_OUTPUTS":
            "- the 100% time slice outputs." in access,
        "TAPE_ACCESS_IS_COARSE_HIGH_LATENCY":
            "magnetic-tape-backed" in access and "6.6 TB" in access,
        "CLEANED_HALO_OFFICIAL_DOI_PRESENT":
            "10.13139/OLCF/1828535" in access and "*cleaned halo catalogs*" in access,
        "MERGER_TREES_NOT_ALWAYS_IN_PORTAL_WEB_INTERFACE":
            "merger trees" in access and "not yet exposed via the web interface" in access,
        "CLEANED_DTYPE_HAS_MEMBERSHIP_N_TOTAL_AND_PROGEN_FIELDS":
            {"N_total","N_mainprog","haloindex_mainprog"}.issubset(schema_fields),
        "NO_EXPLICIT_M200C_IN_PINNED_STANDARD_LOADER_SCHEMA":
            not any(re.search(r"(?i)m_?200c|mass_?200c",f) for f in schema_fields),
        "LOADER_N_TOTAL_DIFFERS_FROM_RAW_N_ON_CLEANING":
            "If we load cleaned, 'N' no longer has meaning" in loader,
        "SOURCE_ONLY_NOT_A_LIVE_PORTAL_INVENTORY":True
    }
    require(all(checks.values()),"official source assertions failed: "+
            ",".join(k for k,v in checks.items() if not v))
    return {"checks":checks,"full_primary_redshifts_documented":list(full),
            "exact_simulation_csv_Full_Outputs":table_mode,
            "documented_loader_schema_fields_checked":sorted(schema_fields),
            "documented_M200c_field_found_in_pinned_loader":False}

def negatives(s):
    def reject(label,callback):
        try:callback()
        except (ValueError,TypeError,KeyError):
            print("B9_SYNTHETIC_NEGATIVE_REJECT",label,flush=True)
        else:raise AssertionError("accepted tampered source fixture "+label)
    sim=s["simulation_document"]
    reject("FAKE_FULL_Z095",lambda:documented_product_checks(
        {**s,"simulation_document":sim.replace("1.1, 0.8","1.1, 0.95, 0.8")}))
    reject("WRONG_SIMULATION_ID",lambda:documented_product_checks(
        {**s,"simulation_table":s["simulation_table"].replace(SIM,"AbacusSummit_base_c000_ph999",1)}))
    reject("WRONG_FULL_DESIGNATION",lambda:documented_product_checks(
        {**s,"simulation_table":s["simulation_table"].replace("0.1     , Full","0.1     , none",1)}))
    reject("FAKE_SECONDARY_FULL_RV",lambda:documented_product_checks(
        {**s,"data_products_document":s["data_products_document"].replace(
            "positions/velocities nor the field particles","complete same-epoch field and RV",1)}))
    reject("CLEANED_N_TOTAL_RENAMED_TO_M200C",lambda:documented_product_checks(
        {**s,"loader":s["loader"].replace("('N_total', np.uint32)","('M200c', np.float32)",1)}))
    reject("PRIMARY_NEIGHBOR_RECAST_AS_SECONDARY",lambda:documented_product_checks(
        {**s,"simulation_document":sim.replace("1.1, 0.8","0.95, 0.8",1)}))
    reject("FULL_OUTPUTS_FALSELY_ON_DISK",lambda:documented_product_checks(
        {**s,"data_access_document":s["data_access_document"].replace(
            "- the 100% time slice outputs.","- the 100% time slice outputs are on disk.",1)}))
    reject("FAKE_M200C_SCHEMA",lambda:documented_product_checks(
        {**s,"loader":s["loader"].replace(
            "('N_total', np.uint32)","('M200c', np.float32),\n        ('N_total', np.uint32)",1)}))
    print("E17D2B3B9_EIGHT_SYNTHETIC_NEGATIVES_PASS_NO_REAL_ASDF",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ("simulations","access","products","table","loader"):
        ap.add_argument("--"+name,type=Path,required=True)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    pre,b8,b7=pinned_parent()
    sources=pinned_sources(pre,args)
    checked=documented_product_checks(sources)
    if args.self_test:negatives(sources)
    report={
      "stage":pre["stage"],"status":"PINNED_STANDARD_OFFICIAL_DOCUMENTED_PRODUCT_GATE_PASS_NO_SAME_EPOCH_M200C_IDENTIFIED",
      "prospective_prereg_git_blob":PINS[0][2],"source_pin_commit":pre["source_pins"]["commit"],
      "official_documented_product_checks":checked,
      "scope":"Only these FIVE pinned official text docs/table/loader. No live Globus listing, no claim excluding third-party or special provider-created data.",
      "c000_ph000_z095_full_particle_snapshot_in_pinned_official_standard_products":False,
      "nearest_primary_full_outputs_are_NOT_same_z":["z1.100","z0.800"],
      "official_standard_cleaned_fields_are_NOT_measured_M200c":True,
      "independent_same_epoch_M200c_in_pinned_standard_catalog_schema_documented":False,
      "actual_verified_same_epoch_M200c_independent_product_identifier":None,
      "future_real_M200c_a_gate":"BLOCKED: independent official exact-cosmology exact-phase exact-epoch mass definition+provenance required before any further local transfers or computation",
      "original_local_B7_JSON_preserved":True,"user_WSL_run":False,
      "real_ASDF_or_PID_RV_cleaned_tree_snapshot_data_read":False,
      "observed_odd_SEALED":True,"original_E8_E16_48k_frozen":True
    }
    if args.output:
        require(args.output.suffix==".json" and not args.output.exists(),
                "new source-only JSON path needed, refuse overwrite")
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print("E17D2B3B9_SOURCE_ONLY_REPORT",args.output,flush=True)
    print("E17D2B3B9_PINNED_SAME_EPOCH_AVAILABILITY_SOURCE_ONLY_PASS",flush=True)
    print("E17D2B3B9_Z095_FULL_SNAPSHOT_NOT_IN_DOCUMENTED_C000_PH000_TIMESLICES",flush=True)
    print("E17D2B3B9_CLEANED_L1_FIELDS_NOT_INDEPENDENT_M200C",flush=True)
    print("E17D2B3B9_NO_NEW_REAL_ASDF_NO_WSL_OBSERVED_ODD_SEALED",flush=True)
if __name__=="__main__":
    main()
