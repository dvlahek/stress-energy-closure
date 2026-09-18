# DESI DR1 genuine cross-population odd-dipole test — frozen protocol

Date frozen: 2026-09-18  
Branch: `dr1-genuine-multitracer`

## Purpose

The publication forecast uses a genuine multi-tracer structure, while the current DESI DR1 deployment uses luminosity-ranked subdivisions of BGS. This test closes that structural mismatch with distinct DESI tracer populations before inspecting the resulting odd statistic.

## Primary test

- Tracers: LRG × ELG_LOPnotqso.
- Orientation: LRG -> ELG. The orientation is fixed by tracer identity before the odd statistic is evaluated.
- Survey release: public DESI DR1 LSS clustering catalogs, v1.5.
- Caps: NGC + SGC.
- Shared redshift range: 0.80 <= z < 1.10.
- Analysis redshift bins: [0.80, 0.95), [0.95, 1.10).
- Separation bins: 20, 40, 60, 80, 100, 120, 140 h^-1 Mpc.
- Angular pair cut: theta >= 0.05 deg.
- Pair compression: at most 48 eligible neighbours per anchor, sampled uniformly with inverse-probability pair weights.
- Random density target: 2 times the selected data density for each tracer/cap/redshift stratum.
- Jackknife: 30 common sky regions, defined from pooled tracer geometry and then frozen for data and randoms.
- Seed: 20260918.

The estimator is the cross-population Landy-Szalay odd dipole
`(D1D2 - D1R2 - R1D2 + R1R2)/R1R2`, with the odd pair factor `3 mu`.
No luminosity or halo-mass split is applied inside either tracer.

The matched-filter wake shape is generated from the already frozen hidden-state direction used in the manuscript. Its overall amplitude is free. The high-z shape is computed directly from the hidden response and is not given the BGS-specific absolute forecast calibration. Standard Doppler and smooth odd broadband terms are nuisance components.

## Secondary diagnostic

BGS_BRIGHT-21.5 × LRG is evaluated only as a low-redshift boundary diagnostic, using BGS 0.30 <= z < 0.40 and LRG 0.40 <= z < 0.50. This is not the primary test because the two clustering selections meet at the z=0.4 boundary rather than occupying a broad common volume.

BGS × ELG is not tested because their DR1 clustering redshift supports do not produce a useful 20–140 h^-1 Mpc cross-pair sample.

## Fixed interpretation

The primary result is the LRG × ELG coefficient and its jackknife uncertainty. It is an exploratory public-data cross-population test, not a DESI collaboration likelihood and not a detection claim. No binning, orientation, nuisance basis, pair-sampling rule, seed, or tracer choice will be changed after inspecting the result.

A stronger publication claim would additionally require cross-population mock/window validation.