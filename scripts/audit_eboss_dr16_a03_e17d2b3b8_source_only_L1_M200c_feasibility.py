#!/usr/bin/env python3
"""E17D2b3b8 pinned SOURCE ONLY. Never open real ASDF, PID, RV, odd or WSL data."""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
PRE=ROOT/"source_data/eboss_dr16_a03_e17d2b3b8_official_source_only_L1_threshold_and_M200c_feasibility_prereg_2026-09-29.json"
B7=ROOT/"source_data/eboss_dr16_a03_e17d2b3b7r1_user_uploaded_local_id_N_result_and_existing_output_failure_diagnostic_2026-09-29.json"
B2=ROOT/"source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28/e17d2b3b2_pinned_header_copy_l1_threshold_vs_200crit.json"
BLOBS=("40abc1616ea85d96ec9cb2ba2ca26dc6a8d6fe22","f27ed284c95e48548e50336da04c1a1d4b137186","cb911e0f7864738695f20fe3228826374d343378")
DOCS=ROOT/"official-abacussummit/docs/data-products.rst"
COMP=ROOT/"official-abacussummit/docs/compaso.rst"
LOADER=ROOT/"official-abacusutils/abacusnbody/data/compaso_halo_catalog.py"
SOURCE_BLOBS=("f7c847b174365744b23cb6d1ebaeb010a5bd6ca7","a5170574467501e0b9adb4a71e8e432a32843dbf","adb16aee1cbac863db5301c6e380938ee7f76547")

def need(ok,msg):
    if not ok: raise ValueError("E17D2B3B8_SOURCE_ONLY_FAIL_CLOSED: "+msg)

def blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()

def parent_inputs():
    out=[]
    for p,sha in zip((PRE,B7,B2),BLOBS):
        need(p.is_file() and not p.is_symlink() and blob(p.read_bytes())==sha,"parent Git blob drift "+p.name)
        out.append(json.loads(p.read_bytes()))
    pre,b7,b2=out
    need(pre["status"]=="PREREG_BEFORE_B8_SOURCE_ONLY_CI_NO_MORE_REAL_HALO_IO"
         and pre["no_new_real_data_transfer"] and all(pre["absolute_STOP"].values()),"protocol drift")
    need(pre["parents"]["b7_prospective_data_field_protocol_git_blob"]
         =="296bc22e3d38a414017cc018e1f8ba8658899278"
         and pre["parents"]["b7_original_numeric_reader_git_blob"]
         =="f3df40e44a87aa0784fe1c017ce431da089920d8","original B7 pins drift")
    return pre,b7,b2

def source_inputs(paths):
    need(tuple(p.resolve() for p in paths)==tuple(p.resolve() for p in (DOCS,COMP,LOADER)),
         "only exact sparse checkout official three text files allowed")
    out=[]
    for p,h in zip(paths,SOURCE_BLOBS):
        need(p.is_file() and not p.is_symlink() and p.suffix in (".rst",".py")
             and blob(p.read_bytes())==h,"official source Git blob drift "+p.name)
        out.append(p.read_text(encoding="utf-8"))
    return out

def source_checks(d,c,l):
    v={
      "N_DOCUMENTED_GT35":"We output properties for all L1 halos with more than 35 particles." in c,
      "COMPETITIVE_L1_ASSIGNMENT":bool(re.search(r"Now a particle is assigned to the new group\s+if is previously unassigned",c)),
      "L1_BASED_ON_L0_PARTICLES":bool(re.search(r"all L1/L2 finding and all halo statistics are based\s+solely on the particles in the L0 halo",c)),
      "N_IS_ASSIGNED_PRIMARY_FIELD":"The number of particles in this halo.  This is the primary halo mass field." in d,
      "SODENSITY_EPOCH_RELATIVE_TO_MEAN":bool(re.search(r"spherical overdensity definition\s+is a function of epoch.*?SODensityL1.*?mean cosmic density",d,re.S)),
      "SO_RADIUS_HAS_NONCROSSING_SENTINEL":"or a constant if the SO crossing is not reached" in d,
      "SO_RADIUS_SENTINEL_EXPLAINED":bool(re.search(r"Some fraction of low-mass halos have.*?SO_radius.*?same\s+value.*?does not drop below the SO threshold",d,re.S)),
      "L2COM_PERCENTILES_ARE_CONDENSED":bool(re.search(r"Radii of this percentage of mass, relative to L2\s+center.*?Expressed as ratios of r100 and condensed",d,re.S)),
      "SECONDARY_NO_SAME_EPOCH_RV_FIELD":bool(re.search(r"For the 21 secondary redshifts, we output the halo catalogs and the halo\s+subsample particle IDs \(w/densities and sticky L2 tag\) only, so not the\s+positions/velocities nor the field particles",d)),
      "LOADER_CLEANED_TRUE_DEFAULT":bool(re.search(r"def __init__\(\s*self,\s*path,\s*cleaned=True",l)),
      "LOADER_CLEANED_N_NOT_RAW_N":"If we load cleaned, 'N' no longer has meaning" in l,
      "LOADER_SO_RADIUS_BOX_SCALE":"r'SO(?:_L2max)?(?:_central_particle|_radius)'" in l,
      "LOADER_PERCENTILE_RATIO_I16":"raw[m[0] + '_i16'] * raw['r100' + m['com']] / INT16SCALE * box" in l,
      "HEADER_NOT_DIRECTORY_REDSHIFT":"Always use the redshift in the header" in d}
    need(all(v.values()),"upstream source checks: "+", ".join(k for k,passed in v.items() if not passed))
    return v

def boundary(b7,c):
    v=b7["data_scope"]
    need("more than 35 particles" in c and v["min_N"]==35
         and v["total_rows"]==11676687 and v["raw_file"]=="halo_info_000.asdf",
         "raw N=35 / documentation boundary parent mismatch")
    need(b7["limitations"]["halo_M200c_measured"] is False
         and b7["limitations"]["whole_34_file_catalogue_sampled"] is False,
         "source-only raw result promoted to physics")
    return "DOC_GT35_VS_USER_RAW_MIN35_UNRESOLVED_NO_CUT"

def thresholds(b2):
    src=b2["official_source_copy_six_exact_numeric_fields"]
    t=200/src["OmegaNow_m"]
    need(math.isclose(t,b2["copy_delta_200critical_in_mean_cosmic_density_units"],rel_tol=1e-12)
         and not math.isclose(src["SODensityL1"],t,rel_tol=1e-6)
         and b2["same_object_true_M200c_200critical_remeasured"] is False,
         "wrong L1/200critical threshold identity or evidence")
    return {"z_source_copy":src["Redshift"],"SODensityL1_over_mean":src["SODensityL1"],
            "200critical_over_mean":t,"L1_over_200critical":src["SODensityL1"]/t,
            "NOT_an_individual_M200c_conversion":True}

def synthetic_rank_toy():
    low=[i/50 for i in range(1,51)]
    tail=[2+i/100 for i in range(1,34)]
    a=low+[1.1]*16+[2.0]+tail
    b=low+[1.9]*16+[2.0]+tail
    ranks=(10,25,33,50,67,75,90,95,98,100)
    need(len(a)==len(b)==100 and all(a[k-1]==b[k-1] for k in ranks),"toy quantiles changed")
    ac=sum(x<=1.5 for x in a)
    bc=sum(x<=1.5 for x in b)
    need(ac==66 and bc==50,"toy enclosed counts should disagree")
    return {"kind":"SYNTHETIC_MATHEMATICAL_RANK_TOY_NOT_ABACUS",
            "same_percentile_ranks":list(ranks),"case_A_count_at_radius_1p5":ac,
            "case_B_count_at_radius_1p5":bc,"same_quantiles_do_not_fix_full_enclosed_profile":True}

def neg(label,fn):
    try: fn()
    except (ValueError,KeyError,TypeError): print("B8_NEGATIVE_REJECT",label,flush=True)
    else: raise AssertionError("Accepted negative "+label)

def synthetic_negatives(pre,b7,b2,d,c,l):
    broken=copy.deepcopy(b7);broken["data_scope"]["min_N"]=36
    neg("N35_PARENT_DRIFT",lambda:boundary(broken,c))
    neg("GT35_DOCUMENTATION_WORDING_DRIFT",lambda:source_checks(d,c.replace("more than 35 particles","35 or more particles"),l))
    broken=copy.deepcopy(b2)
    broken["official_source_copy_six_exact_numeric_fields"]["SODensityL1"]=broken["copy_delta_200critical_in_mean_cosmic_density_units"]
    neg("FALSE_EQUAL_SO_AND_200CRITICAL",lambda:thresholds(broken))
    neg("FALSE_SECONDARY_RV_SOURCE",lambda:source_checks(d.replace("positions/velocities nor the field particles","unlimited full RV and field"),c,l))
    broken=copy.deepcopy(b7);broken["limitations"]["halo_M200c_measured"]=True
    neg("L1_N_PROMOTED_TO_MEASURED_M200C",lambda:boundary(broken,c))
    neg("SO_RADIUS_SENTINEL_ERASED",lambda:source_checks(d.replace("or a constant if the SO crossing is not reached","always measured R200c"),c,l))
    neg("REAL_ASDF_NOT_ALLOWED_SOURCE",lambda:need(Path("halo_info_000.asdf").suffix in (".py",".rst"),"binary input"))
    neg("LOADER_CLEANED_DEFAULT_DRIFT",lambda:source_checks(d,c,l.replace("cleaned=True,","cleaned=False,")))
    print("E17D2B3B8_EIGHT_SYNTHETIC_NEGATIVES_PASS",flush=True)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--docs",type=Path,required=True)
    ap.add_argument("--compaso",type=Path,required=True)
    ap.add_argument("--loader",type=Path,required=True)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--output",type=Path)
    a=ap.parse_args()
    pre,b7,b2=parent_inputs()
    d,c,l=source_inputs((a.docs,a.compaso,a.loader))
    checks=source_checks(d,c,l)
    classification=boundary(b7,c)
    comparison=thresholds(b2)
    toy=synthetic_rank_toy()
    if a.self_test: synthetic_negatives(pre,b7,b2,d,c,l)
    report={"stage":pre["stage"],"status":"PINNED_SOURCE_ONLY_AND_SYNTHETIC_PASS_NO_NEW_REAL_DATA",
       "prospective_prereg_git_blob":BLOBS[0],"official_source_checks":checks,
       "N35_boundary":classification,"fixed_source_threshold_comparison":comparison,
       "synthetic_rank_profile_example":toy,
       "secondary_available":"halo_info plus selected 10-percent PID/density/sticky tag; no same-time RV/field",
       "future_measured_M200c_requires":"same-epoch enclosed TOTAL mass profile about defined center and 200critical crossing; absent in current raw id/N and finite assigned-L1 percentiles",
       "user_local_numeric_result_not_replayed_by_CI":True,"no_new_real_ASDF_or_PID_RV_read":True,
       "M200c_measured":False,"merger_tree_measured":False,
       "physical_original_Fplus_Fminus_wake_or_galaxy_B_measured":False,"observed_odd_SEALED":True}
    if a.output:
        need(a.output.suffix==".json" and not a.output.exists(),"source report must be a new JSON")
        with a.output.open("x",encoding="utf-8") as f:
            json.dump(report,f,indent=2,allow_nan=False);f.write("\n")
    print("E17D2B3B8_PINNED_OFFICIAL_SOURCE_ONLY_FEASIBILITY_PASS",flush=True)
    print("E17D2B3B8_N35_DOC_DISCREPANCY_UNRESOLVED_NO_CUT",flush=True)
    print("E17D2B3B8_NO_CERTIFIED_M200C_FROM_CURRENT_L1_CATALOG",flush=True)
    print("E17D2B3B8_NO_NEW_REAL_HALO_READ_OBSERVED_ODD_SEALED",flush=True)
if __name__=="__main__":
    main()
