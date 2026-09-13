#!/usr/bin/env python3
"""Hidden-state wake Fisher for a validated DESI-BGS HOD halo-mass split.

The underlying BGS selection is the public MXXL/hodpy HOD evaluated at r<19.5
with representative observer-frame g-r=1.0.  Our independent selection check
finds 852.6 galaxies/deg^2 over 0.05<z<0.4 and a redshift-distribution peak at
z=0.15, closely matching the final DESI BGS-BRIGHT target density (864/deg^2)
and the Ge, Pasquini & Tan forecast peak.  Galaxies are split at
M_h=10^13.75 h^-1 Msun, the split optimized in Ge et al. JCAP 05 (2024) 108.

Seven Delta-z=0.05 bins use the directly reconstructed HOD dN/dz and ST99
biases.  Counts below are midpoint integrations of the saved colour_1.0 profile
for 14,000 deg^2.  The Fisher machinery is identical to the published-split
Bonvin benchmark and retains the conservative independent odd nuisances per
redshift bin.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import wake_desi_bonvin_fisher as w

# low = M_h < 10^13.75 h^-1 Msun, high = M_h >= threshold
# Columns originate from run 34753916185, bgs_profiles_colour_1.0.csv.
w.SURVEY = [
    {"z":0.075,"zlo":0.05,"zhi":0.10,"N_faint":1349967,"N_bright":125859,"b_faint":0.9149368069,"b_bright":2.6446367510},
    {"z":0.125,"zlo":0.10,"zhi":0.15,"N_faint":2224823,"N_bright":190774,"b_faint":0.9490235822,"b_bright":2.7100639629},
    {"z":0.175,"zlo":0.15,"zhi":0.20,"N_faint":2393243,"N_bright":192429,"b_faint":0.9963186571,"b_bright":2.7760945720},
    {"z":0.225,"zlo":0.20,"zhi":0.25,"N_faint":2107932,"N_bright":164437,"b_faint":1.0588273359,"b_bright":2.8428318856},
    {"z":0.275,"zlo":0.25,"zhi":0.30,"N_faint":1565080,"N_bright":126003,"b_faint":1.1446827618,"b_bright":2.9103092237},
    {"z":0.325,"zlo":0.30,"zhi":0.35,"N_faint":957383,"N_bright":88656,"b_faint":1.2645443573,"b_bright":2.9778921133},
    {"z":0.375,"zlo":0.35,"zhi":0.40,"N_faint":454336,"N_bright":57705,"b_faint":1.4256685142,"b_bright":3.0363199304},
]

if __name__ == "__main__":
    # Run the common survey-specific machinery.
    w.main()
    # Correct provenance text inherited from the common implementation.
    try:
        oi=sys.argv.index("--outdir")
        out=Path(sys.argv[oi+1])
    except Exception:
        out=Path("wake_desi_bonvin_output")
    p=out/"summary.json"
    if p.exists():
        s=json.loads(p.read_text())
        s["survey_source"]=(
            "Validated DESI BGS-BRIGHT HOD selection: public amjsmith/hodpy r<19.5, "
            "representative g-r=1.0, SMT HMF, ST99 bias, 14,000 deg^2, with halo-mass "
            "split M_h=10^13.75 h^-1 Msun following Ge, Pasquini & Tan JCAP 05 (2024) 108."
        )
        s["selection_validation"]={
            "integrated_density_deg2":852.6064,
            "DESI_BGS_BRIGHT_target_density_deg2":864.0,
            "redshift_distribution_peak":0.15,
            "midpoint_total_count_used":sum(x["N_faint"]+x["N_bright"] for x in w.SURVEY),
            "note":"This is a forecast for a halo-mass-selected multi-tracer design, not an already delivered DESI group catalogue."
        }
        s["interpretation"]=(
            "The BGS parent selection is validated against the observed target density and redshift-distribution scale. "
            "The high/low halo-mass split is a survey-design forecast. Absolute wake normalization remains literature-calibrated to Okoli et al. Eq.35."
        )
        p.write_text(json.dumps(s,indent=2)+"\n")
