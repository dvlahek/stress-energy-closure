#!/usr/bin/env python3
from __future__ import annotations
import csv, io, json, math, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"source_data/e63_cosmodc2_outerrim_lightcone_source_preflight_result.json"
TAP="https://irsa.ipac.caltech.edu/TAP/sync"
TABLE="cosmodc2mockv1"
REQ=["galaxy_id","ra_true","dec_true","redshift_true","position_x","position_y","position_z","velocity_x","velocity_y","velocity_z","halo_mass","is_central"]
RA0,DEC0,RAD=55.0,-41.0,2.0
Z0,Z1=0.9,1.0
TOP=4097

def need(c,m):
    if not c: raise RuntimeError(m)

def query(adql):
    data=urllib.parse.urlencode({
      "REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":adql
    }).encode()
    req=urllib.request.Request(TAP,data=data,headers={"User-Agent":"EinsteinVlasovNP-E63/1.0"})
    with urllib.request.urlopen(req,timeout=180) as r:
        raw=r.read().decode("utf-8")
    if raw.lstrip().startswith("<?xml") or "QUERY_STATUS" in raw[:1000]:
        raise RuntimeError("IRSA TAP returned an error document instead of CSV")
    return raw

def parse_csv(raw):
    return list(csv.DictReader(io.StringIO(raw)))

def main():
    schema_q=("SELECT column_name,datatype,unit,description FROM TAP_SCHEMA.columns "
              f"WHERE table_name='{TABLE}'")
    schema=parse_csv(query(schema_q))
    names={r["column_name"] for r in schema}
    missing=[x for x in REQ if x not in names]
    need(not missing,f"Missing required CosmoDC2 columns: {missing}")

    cols=",".join(REQ)
    sample_q=(f"SELECT TOP {TOP} {cols} FROM {TABLE} WHERE "
              f"redshift_true>={Z0} AND redshift_true<{Z1} AND "
              f"1=CONTAINS(POINT('ICRS',ra_true,dec_true),"
              f"CIRCLE('ICRS',{RA0},{DEC0},{RAD}))")
    rows=parse_csv(query(sample_q))
    need(len(rows)>=256,f"Bounded CosmoDC2 probe returned only {len(rows)} rows")

    numeric=[x for x in REQ if x!="galaxy_id"]
    arr={k:np.asarray([float(r[k]) for r in rows],dtype=float) for k in numeric}
    for k,a in arr.items():
        need(np.isfinite(a).all(),f"Nonfinite values in {k}")
    need(np.all((arr["redshift_true"]>=Z0)&(arr["redshift_true"]<Z1)),"redshift_true range drift")

    pos=np.column_stack([arr["position_x"],arr["position_y"],arr["position_z"]])
    vel=np.column_stack([arr["velocity_x"],arr["velocity_y"],arr["velocity_z"]])
    r=np.linalg.norm(pos,axis=1)
    v=np.linalg.norm(vel,axis=1)

    result={
      "stage":"E63_COSMODC2_OUTERRIM_LIGHTCONE_SOURCE_PREFLIGHT",
      "date":"2026-10-01",
      "status":"PASS_PUBLIC_OUTERRIM_LIGHTCONE_SOURCE_AVAILABLE",
      "tap_endpoint":"https://irsa.ipac.caltech.edu/TAP",
      "table":TABLE,
      "required_columns":REQ,
      "schema_column_count":len(schema),
      "bounded_probe":{
        "center_ra_deg":RA0,"center_dec_deg":DEC0,"radius_deg":RAD,
        "redshift_true_range":[Z0,Z1],"rows_returned":len(rows),"top_limit":TOP
      },
      "sample_diagnostics":{
        "redshift_true_min":float(arr["redshift_true"].min()),
        "redshift_true_max":float(arr["redshift_true"].max()),
        "position_radius_median":float(np.median(r)),
        "position_radius_min":float(r.min()),
        "position_radius_max":float(r.max()),
        "velocity_norm_median_km_s":float(np.median(v)),
        "velocity_norm_min_km_s":float(v.min()),
        "velocity_norm_max_km_s":float(v.max()),
        "central_fraction":float(np.mean(arr["is_central"]>0.5)),
        "halo_mass_median":float(np.median(arr["halo_mass"]))
      },
      "decision":{
        "source_preflight_pass":True,
        "interpretation":"PUBLIC_OUTERRIM_LIGHTCONE_WITH_VELOCITY_TRUTH_AVAILABLE; NEXT_PREREGISTER_E63_SCIENCE_TRANSFER_TEST",
        "important_limit":"This bounded query is a source/schema preflight only, not the final lightcone tracer or geometry."
      },
      "observed_eBOSS_rows_used":False,
      "observed_odd_used":False
    }
    OUT.parent.mkdir(parents=True,exist_ok=True)
    OUT.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("E63_COSMODC2_SOURCE_PREFLIGHT_PASS")
    print("ROWS",len(rows))
    print("CENTRAL_FRACTION",result["sample_diagnostics"]["central_fraction"])
    print("Z_RANGE",result["sample_diagnostics"]["redshift_true_min"],result["sample_diagnostics"]["redshift_true_max"])
    print("DECISION",result["decision"]["interpretation"])
    print("REPORT",OUT)

if __name__=="__main__":
    main()
