# eBOSS DR16 LRG×ELG common-footprint and radial-selection gate

Status: **not yet certified for unblinding**. This is an input-only gate,
separate from the completed random-only RR counting and nine-realization
window pilot. Do not open observed LRG×ELG odd multipoles to decide cuts.

## Two tracer-specific selections define the cross window

The released DR16 *clustering* LRG and ELG random catalogues are
tracer-specific and are already constructed to sample their respective
survey selection functions. The Landy–Szalay cross estimator uses
`D_L D_E - D_L R_E - R_L D_E + R_L R_E` with the corresponding
independently normalized cross-pair counts; the expected pair window
is sampled by the LRG×ELG `R_L R_E`. A forced identical angular mask
is **not** a prerequisite for the cross estimator. It would be a
new analysis cut requiring prior registration and consistent
application to both tracer data and randoms.

Official MANGLE/BRICKMASK products are independently checked as
selection provenance, not automatically reapplied to already masked
DR16 clustering samples. A literal geometrical intersection may be
useful for a separate systematic-control test only if its precise
definition, sample loss and estimator treatment are predeclared.

## Why the angular joint mask is not a Jaccard threshold

The DR16 clustering randoms sample the tracer-specific selection functions
after footprint, veto and completeness treatment. However, the sampled
angular support of a finite random catalogue, even with exact-count
resampling, is not a geometrically exact continuous intersection.
The observed LRG and ELG randoms have different number densities and
selection prescriptions. A pixel threshold on their coarse sky grids
therefore does not reproduce the published tracer mask.

The published multi-tracer EZmock construction documents a necessary
difference between the input masks: LRG/QSO vetoes are MANGLE polygons
(applied with MAKE_SURVEY/MPLY_TRIM), but ELG vetoes are associated with
DECaLS brick pixels and applied with BRICKMASK. The released BRICKMASK
documentation also identifies extra eBOSS ELG mask operations beyond its
main maskbits assignment. The LRG×ELG *pair selection* is the product
of the two tracer-specific selection functions, sampled by the
published cross-random pair counts R_L R_E. Replacing it with a single
coarse thresholded joint sky map or automatically applying a new
intersection cut is not equivalent.

References:
- Zhao et al., MNRAS 503 (2021) 1149, DR16 matched multi-tracer EZmocks:
  https://academic.oup.com/mnras/article/503/1/1149/6149160
- Ross et al., MNRAS 498 (2020) 2354, DR16 LRG LSS catalogue and vetoes:
  https://academic.oup.com/mnras/article/498/2/2354/5900562
- SDSS DR16 LSS catalogue documentation:
  https://www.sdss4.org/dr17/spectro/lss/
- BRICKMASK code and eBOSS ELG-specific notes:
  https://github.com/cheng-zhao/brickmask

## Current separately auditable input-only checks

The complete observed random-only fine RR is available in both caps and
all four candidate bins with positive support in all fine cells and an
independent coarse rebin check. Nine preregistered matched mock IDs have
completed a high-z random-only input/window audit. The next stage,
`scripts/audit_eboss_dr16_9mock_density_radial.py`, compares all nine
mock random pixel maps with the observed maps after *exact-count* matching
separately by cap, tracer and candidate redshift interval. It retains
all fixed thresholds 1, 5, 10, 15 and 25; normalized sky-map total
variation does not require a chosen support threshold. Two fixed-seed
multinomial replicates of each pixel map provide a finite-input
resampling reference. These are correlated through the same parent
catalogue and must not be interpreted as independent mocks.

The preregistered local fine weighted n(z) audit has now completed all
36 comparisons for nine fixed realistic mocks, both caps and both
tracers on forty 0.01-wide z bins on [0.6,1.0). Exact shared ELG
chunk labels were audited separately. This closes the fixed fine
random n(z) *input check*, not the physical pair-window validation.
The ELG radial random assignment depends on imaging depth in the
DR16 LSS prescription. NGC/SGC and ELG chunk normalizations cannot
be replaced by one pooled n(z).

The three SHA-pinned published ELG extra MANGLE polygons have also been
evaluated on prospectively stratified samples of the released ELG
clustering RANDOMS. For observed NGC and SGC and mock 0001 NGC and
SGC, each case has eight exact chunk-by-z strata of 2500 random rows:
0 of 80000 sampled random rows falls inside any of the three extra
polygons. This does not prove the full catalogues contain no membership.
The additional upstream HEALPix bit 8 and the complete BRICKMASK
image files are separate open checks. Published clustering random
catalogues must not be newly vetoed solely because this diagnostic
was run.

The complete official eleven-polygon source-byte cohort is now SHA256-
pinned: seven published LRG/QSO veto polygons, the QSO+LRG noveto
footprint polygon, and three additional ELG veto polygons. The exact
source filenames, byte lengths and SHA256s are retained in
`source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json`.
The three input-index SHA256s and the local eleven-file manifest agree.
This closes the polygon *byte-provenance* gate only. The published
tracer-specific LRG veto composition, ELG brickmask image families,
ELG extra pixel-bit coordinate convention, and continuous pair-window
membership are separate open gates.

## Required production checks before opening observed odd data

1. Pin the official LRG MANGLE footprint/veto/sector products and
   the ELG DECaLS brick-mask and extra-mask products. Record exact
   source path, release version, SHA256, veto flags and completeness
   semantics. Obtain these files from the official release; do not
   synthesize their geometry from coarse occupancy.
2. Audit the two *tracer-specific* selection functions and their
   pair-level cross window, testing NGC and SGC separately, including
   exact ELG chunk definitions. Validate official mask membership
   and weighted accepted-random fractions against the already pinned
   released DR16 clustering randoms. Do not reapply a veto blindly:
   first determine if the published clustering sample already
   excludes that region. Any additional intersection cut must be
   prospectively frozen and applied consistently to both tracers'
   data and random catalogues, without using odd measurements.
3. Independently compare **weighted**, finer n(z), including ELG depth
   or chunk-conditioned distributions, for each observed/matched-mock
   random catalogue and both caps. Declare the fine z-grid, weights,
   tolerances and treatment of near-zero weights *before* this test.
4. Validate pair-level cross-LS normalization, analytic even-to-odd
   leakage, mock-galaxy estimator closure and joint 18D covariance
   under the two published tracer-specific masks and their cross
   window. Nine mock randoms are a window/selection pilot only. The inferential mock membership, finite-mock likelihood
   and stopping rule must be separately preregistered.
5. Only then freeze and commit the full prospective eBOSS analysis
   protocol and unblind the actual observed odd-sector vector.

This checklist records completed fine weighted random n(z) and
published-polygon input checks. It does **not** assert an exact
physical angular mask, validated pair-level convolution,
mock-galaxy covariance or any observed galaxy odd statistic.

## Additional ELG pixel-bit coordinate-contract check (open)

The byte-identical published `scripts/eBOSS_ELG_extra.py` currently
forms `theta = radians(90-dec)` and `phi = radians(360-ra)` and then
calls `healpy.ang2pix(..., theta, phi, lonlat=True)` when constructing
its extra bit `2**8`. The documented `healpy.ang2pix` convention is:
`lonlat=True` requires longitude and latitude **in degrees**;
`lonlat=False` takes colatitude and longitude **in radians**.
Those two statements are not syntactically consistent as written.
Therefore the extra pixel bit is intentionally **not applied** in
the current source and polygon membership pilots. Preserve the
upstream file byte-for-byte, but resolve and independently test the
actual reference-mask pixel convention before implementing bit 8.
Do not silently replace the published formula with an assumed fix;
compare original versus alternative pixel indexing with official
reference bitmaps or a verified output catalogue.

Primary sources:
- https://github.com/cheng-zhao/brickmask/blob/b9eb684a579b56ec3dbdb46549224be7e3fa2830/scripts/eBOSS_ELG_extra.py
- https://healpy.readthedocs.io/en/stable/generated/healpy.pixelfunc.ang2pix.html
