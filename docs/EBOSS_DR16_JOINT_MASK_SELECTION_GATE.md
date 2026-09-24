# eBOSS DR16 LRG×ELG common-footprint and radial-selection gate

Status: **not yet certified for unblinding**. This is an input-only gate,
separate from the completed random-only RR counting and nine-realization
window pilot. Do not open observed LRG×ELG odd multipoles to decide cuts.

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
main maskbits assignment. The actual LRG×ELG common footprint must
combine the two *published* selection rules, including their vetoes,
rather than intersecting two thresholded random maps.

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

The existing mock random shard reports retain only four broad unweighted
redshift-bin counts for each tracer/cap. Thus the new radial check can
compare coarse normalized redshift fractions but cannot yet certify
weighted sub-bin n(z). The ELG radial random assignment is dependent
on imaging depth in the DR16 LSS prescription. NGC/SGC and ELG chunk
normalizations cannot be replaced by one pooled n(z).

## Required production checks before opening observed odd data

1. Pin the official LRG MANGLE footprint/veto/sector products and
   the ELG DECaLS brick-mask and extra-mask products. Record exact
   source path, release version, SHA256, veto flags and completeness
   semantics. Obtain these files from the official release; do not
   synthesize their geometry from coarse occupancy.
2. Construct the common *tracer-pair* selection under the published
   LRG and ELG angular rules, testing NGC and SGC separately, including
   the ELG chunk definitions. Validate membership and weighted
   accepted-random fractions against the already pinned DR16 catalogues
   without changing a cut based on any odd measurement.
3. Independently compare **weighted**, finer n(z), including ELG depth
   or chunk-conditioned distributions, for each observed/matched-mock
   random catalogue and both caps. Declare the fine z-grid, weights,
   tolerances and treatment of near-zero weights *before* this test.
4. Validate pair-level LS normalization, analytic even-to-odd leakage,
   mock-galaxy estimator closure and joint 18D covariance under the
   exact common footprint. Nine mock randoms are a window/selection
   pilot only. The inferential mock membership, finite-mock likelihood
   and stopping rule must be separately preregistered.
5. Only then freeze and commit the full prospective eBOSS analysis
   protocol and unblind the actual observed odd-sector vector.

This checklist does not assert that official mask files have been
retrieved or that any final physical mask, weighted n(z), mock-galaxy
covariance or observed galaxy odd statistic has been certified.
