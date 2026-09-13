#!/usr/bin/env python3
"""Full hidden-state wake Fisher for alternative DESI-BGS halo-mass splits.

The parent BGS-BRIGHT selection is the independently validated r<19.5,
representative g-r=1.0 HOD proxy.  Run 34754436382 reconstructed the same
parent sample for several halo-mass thresholds.  This script inserts those
reconstructed counts and ST99 biases into the *same* full survey-specific
hidden-state Fisher used for the 10^13.75 benchmark.

The purpose is deliberately narrower than the simple threshold-ranking
diagnostic: determine if the nonlinear CLASS-validated, nuisance-projected
hidden-state S/N itself prefers a threshold above 10^13.75 h^-1 Msun.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import wake_desi_bonvin_fisher as w

SURVEYS = {
    "14.00": [
        {"z":0.075,"zlo":0.05,"zhi":0.10,"N_faint":1389965,"N_bright":86088,"b_faint":0.9424776852,"b_bright":2.9957900109},
        {"z":0.125,"zlo":0.10,"zhi":0.15,"N_faint":2288003,"N_bright":128009,"b_faint":0.9764313329,"b_bright":3.0803858158},
        {"z":0.175,"zlo":0.15,"zhi":0.20,"N_faint":2459939,"N_bright":126251,"b_faint":1.0239334323,"b_bright":3.1678236524},
        {"z":0.225,"zlo":0.20,"zhi":0.25,"N_faint":2167828,"N_bright":105095,"b_faint":1.0873671917,"b_bright":3.2586221239},
        {"z":0.275,"zlo":0.25,"zhi":0.30,"N_faint":1613560,"N_bright":78067,"b_faint":1.1755548032,"b_bright":3.3540005075},
        {"z":0.325,"zlo":0.30,"zhi":0.35,"N_faint":993704,"N_bright":52844,"b_faint":1.3007762506,"b_bright":3.4561447674},
        {"z":0.375,"zlo":0.35,"zhi":0.40,"N_faint":480073,"N_bright":32441,"b_faint":1.4747796176,"b_bright":3.5637956567},
    ],
    "14.25": [
        {"z":0.075,"zlo":0.05,"zhi":0.10,"N_faint":1423407,"N_bright":52946,"b_faint":0.9728482192,"b_bright":3.4612934505},
        {"z":0.125,"zlo":0.10,"zhi":0.15,"N_faint":2339617,"N_bright":76907,"b_faint":1.0061752371,"b_bright":3.5706389604},
        {"z":0.175,"zlo":0.15,"zhi":0.20,"N_faint":2512876,"N_bright":73899,"b_faint":1.0533231208,"b_bright":3.6849292234},
        {"z":0.225,"zlo":0.20,"zhi":0.25,"N_faint":2213716,"N_bright":59776,"b_faint":1.1170256109,"b_bright":3.8049062420},
        {"z":0.275,"zlo":0.25,"zhi":0.30,"N_faint":1649078,"N_bright":43046,"b_faint":1.2067085611,"b_bright":3.9324347572},
        {"z":0.325,"zlo":0.30,"zhi":0.35,"N_faint":1018757,"N_bright":28191,"b_faint":1.3359332344,"b_bright":4.0714253872},
        {"z":0.375,"zlo":0.35,"zhi":0.40,"N_faint":496188,"N_bright":16634,"b_faint":1.5191947509,"b_bright":4.2271139619},
    ],
    "14.50": [
        {"z":0.075,"zlo":0.05,"zhi":0.10,"N_faint":1448667,"N_bright":28037,"b_faint":1.0034824319,"b_bright":4.0905906869},
        {"z":0.125,"zlo":0.10,"zhi":0.15,"N_faint":2377584,"N_bright":39511,"b_faint":1.0355565997,"b_bright":4.2332051090},
        {"z":0.175,"zlo":0.15,"zhi":0.20,"N_faint":2550648,"N_bright":36744,"b_faint":1.0816883183,"b_bright":4.3830417181},
        {"z":0.225,"zlo":0.20,"zhi":0.25,"N_faint":2245347,"N_bright":28707,"b_faint":1.1449275670,"b_bright":4.5409146816},
        {"z":0.275,"zlo":0.25,"zhi":0.30,"N_faint":1672625,"N_bright":19953,"b_faint":1.2352335814,"b_bright":4.7088889520},
        {"z":0.325,"zlo":0.30,"zhi":0.35,"N_faint":1034638,"N_bright":12641,"b_faint":1.3672250648,"b_bright":4.8914928300},
        {"z":0.375,"zlo":0.35,"zhi":0.40,"N_faint":505785,"N_bright":7254,"b_faint":1.5572041188,"b_bright":5.0973093238},
    ],
}

def main():
    p=argparse.ArgumentParser(add_help=False)
    p.add_argument("--split",required=True,choices=sorted(SURVEYS))
    known,rest=p.parse_known_args()
    split=known.split
    w.SURVEY=SURVEYS[split]
    # Remove our wrapper-only argument before handing control to the common parser.
    sys.argv=[sys.argv[0]]+rest
    w.main()
    # Annotate provenance after the common forecast finishes.
    out=Path("wake_desi_bonvin_output")
    if "--outdir" in rest:
        out=Path(rest[rest.index("--outdir")+1])
    sf=out/"summary.json"
    if sf.exists():
        s=json.loads(sf.read_text())
        s["survey_source"]=(
            "Validated DESI BGS-BRIGHT HOD proxy (r<19.5, representative g-r=1.0, "
            "SMT HMF, ST99 bias, 14,000 deg^2), halo-mass split log10(M/h^-1 Msun)="+split
        )
        s["halo_mass_split_log10"]=float(split)
        s["parent_selection_validation"]={
            "target_density_deg2":852.6064,
            "DESI_BGS_BRIGHT_reference_deg2":864.0,
            "dNdz_peak_z":0.15,
            "profile_source_run":34754436382,
        }
        sf.write_text(json.dumps(s,indent=2)+"\n")

if __name__=="__main__":
    main()
