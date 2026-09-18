# DESI DR1 genuine BGS × LRG boundary diagnostic — frozen protocol

Date frozen: 2026-09-18  
Branch: `dr1-genuine-multitracer`

This is the secondary diagnostic declared in `DESI_DR1_GENUINE_MULTITRACER_PREREG.md`.

- Tracers: LRG -> BGS_BRIGHT-21.5.
- Orientation: LRG is tracer A and BGS is tracer B, fixed before inspecting the odd statistic.
- DESI DR1 LSS v1.5 clustering catalogs, NGC + SGC.
- LRG selection: 0.40 <= z < 0.50.
- BGS selection: 0.30 <= z < 0.40.
- Pair separations: 20–140 h^-1 Mpc in 20 h^-1 Mpc bins.
- Angular pair cut: theta >= 0.05 deg.
- Pair compression: maximum 48 eligible neighbours per anchor with inverse-probability weights.
- Random density target: 2 times each selected data density in cap/redshift strata.
- Jackknife: 30 common sky regions.
- Seed: 20260918.
- Wake matched-filter shape: frozen hidden-state direction, shape-only, free amplitude.
- Standard Doppler plus s^-1 and s^-2 odd broadband components are nuisances.
- Odd octupole is an orthogonal control.

This is deliberately a boundary diagnostic rather than a headline measurement: the BGS and LRG clustering selections meet at z=0.4 and do not share a broad common redshift volume. No settings above are to be changed after inspecting the result.