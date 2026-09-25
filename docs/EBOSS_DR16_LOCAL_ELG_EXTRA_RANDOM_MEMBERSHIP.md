# Blinded DR16 extra ELG veto membership on published randoms

This is a **source-verified additional-polygon diagnostic**, not a
reconstructed BRICKMASK, a new cut on released clustering catalogues,
an exact LRG×ELG common footprint, or a galaxy odd measurement.

## Inputs already certified

The predeclared three ELG extra polygons have SHA256s in
`source_data/eboss_dr16_elg_three_extra_polygon_sha_2026-09-25.json`.
Their full-byte source identities were obtained from official DR16
`ELGmasks/` in workflow run `36100783598`. The upstream
`scripts/eBOSS_ELG_extra.py` is a byte-identical vendored copy of
`cheng-zhao/brickmask` commit
`b9eb684a579b56ec3dbdb46549224be7e3fa2830`.

This diagnostic implements only the upstream MANGLE polygon bits
`9,10,11` for `ELG_centerpost.ply`,
`ELG_TDSSFES_62arcsec.pix.snap.balk.ply`, and
`ebosselg_badphot.26Aug2019.ply`. The script deliberately does
**not** implement the upstream HEALPix discrepancy bit 8,
DECaLS brick-image masks, LRG veto union, or any new veto exclusion.

Each exact ELG chunk is analyzed separately in NGC (`eboss23/25`)
and SGC (`eboss21/22`) in the four fixed candidate bins
`[0.6,0.7),[0.7,0.8),[0.8,0.9),[0.9,1.0)`.
An independently seeded sample of exactly 2500 eligible published
clustering RANDOMS per chunk and z-bin is drawn without replacement.
Only SHA-verified observed and predeclared realistic EZmock RANDOMS
may enter. Membership is descriptive; no threshold is fitted.

## Run the observed SGC pilot in WSL

```bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate
python -m pip install pymangle

mkdir -p eboss_workspace/official_mask_inventory
gh run download 36100783598 \
  --repo dvlahek/stress-energy-closure \
  --name eboss-dr16-elg-three-extra-polygon-sha \
  --dir eboss_workspace/official_mask_inventory

python scripts/audit_eboss_dr16_elg_extra_random_membership.py --self-test

python scripts/audit_eboss_dr16_elg_extra_random_membership.py \
  --ensemble-json eboss_workspace/local_nz/source/ninemock_window_ensemble.json \
  --polygon-manifest eboss_workspace/official_mask_inventory/official_polygon_sha_manifest.json \
  --observed-cache-dir eboss_workspace/local_rr/fits \
  --polygon-dir eboss_workspace/official_mask_inventory/official_polygons \
  --mock-ids none --caps SGC \
  --out-dir eboss_workspace/official_mask_inventory/membership \
  2>&1 | tee eboss_workspace/official_mask_inventory/observed_SGC_elg_extra.log
```

The three source polygons are downloaded only from exact official
URLs and accepted only if their bytes match pinned SHA256s. Published
ELG clustering random FITS are also SHA-verified. Successfully
audited catalogue cases are checkpointed without reopening any
galaxy data. The source sample is fixed before any polygon check.

## One fixed mock, then full previously registered cohort

Run the same command with `--mock-ids 1 --caps SGC` to include
realistic EZmock random 0001. Its sample is independently drawn
from its own fixed parent random catalogue; it is not a galaxy mock
covariance or a pseudo-detection. After checking the pilot,
replace those options with `--mock-ids all --caps all` to evaluate
all nine fixed mocks and both caps. By default any downloaded
compressed mock input is retained in the mock cache. Use
`--purge-new-mock-fits` only when disk space requires deletion;
previously existing cache files are never deleted by that option.

The output `elg_extra_membership_summary.json` contains each
per-cap/case checkpoint and sampled per-polygon maskword fractions.
No sample is removed by the polygon test; it measures membership
inside veto regions on the *already published randoms*.

## Required additional gates

A nonzero membership count is a diagnostic to investigate against
official boundary rules, not proof of a physical odd signal.
Even zero sampled memberships does not prove a mask valid over
the full catalogue. The ELG bit-8 correction and all four brickmask
image families, exact published LRG mask rules, cross-pair normalizations
and mock-galaxy estimator must still be addressed separately.
Do not use a thresholded gridded random mask or an
unregistered common-footprint intersection in their place.
