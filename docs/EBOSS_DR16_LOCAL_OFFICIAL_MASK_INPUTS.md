# Local eBOSS DR16 LRG×ELG official-mask input gate

The completed nine-mock fine weighted random-only `n(z)` control is
**not** an exact joint-footprint certification. This stage discovers
and byte-pins official DR16 angular-mask products before applying any
selection to the *observed galaxy* catalogue or constructing an odd
statistic. The work uses the user's local WSL environment.

## Published products and why the tracers differ

The official DR16 catalogues supply independently constructed LRG and
ELG clustering randoms. The matched EZmock methodology uses MANGLE
polygons for LRG/QSO veto geometry but DECaLS brick-pixel masks for the
ELG veto geometry. The eBOSS implementation of Cheng Zhao's
`brickmask` requires compilation with `-DEBOSS` for this tracer and
a *separate* additional mask pass.

Pin the upstream `cheng-zhao/brickmask` commit
`b9eb684a579b56ec3dbdb46549224be7e3fa2830`.
Its `scripts/eBOSS_ELG_extra.py` identifies three required
additional veto polygons:

- `ELG_centerpost.ply`
- `ELG_TDSSFES_62arcsec.pix.snap.balk.ply`
- `ebosselg_badphot.26Aug2019.ply`

That script also assigns its additional pixel-discrepancy bit
`2**8`; copying just the three polygon filenames is **not**
equivalent to reproducing its ELG mask. The original script and its
bit semantics must be independently tested before imposing any
exclusion.

The four per-chunk ELG brickmask file families are exactly
`mask-eboss21-*.fits.gz`, `mask-eboss22-*.fits.gz`,
`mask-eboss23-*.fits.gz`, `mask-eboss25-*.fits.gz`.
Do not merge the chunk families before determining the final
selection rule. Upstream `CONFIG.md` uses the Legacy Survey DR7
`survey-bricks.fits.gz` brick table with these files and assigns
subsample IDs `0,1,2,3` in that exact chunk order.

Source links:
- SDSS DR16 LSS release documentation:
  https://www.sdss4.org/dr17/spectro/lss/
- Official mask directory mirrored by NERSC:
  https://portal.nersc.gov/project/cosmo/data/sdss/dr17/eboss/lss/catalogs/DR16/ELGmasks/
- Brickmask pinned upstream extra-mask code:
  https://github.com/cheng-zhao/brickmask/blob/b9eb684a579b56ec3dbdb46549224be7e3fa2830/scripts/eBOSS_ELG_extra.py
- Brickmask pinned configuration details:
  https://github.com/cheng-zhao/brickmask/blob/b9eb684a579b56ec3dbdb46549224be7e3fa2830/CONFIG.md

## First local run: only discover public official filenames

```bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate

python scripts/audit_eboss_dr16_official_mask_sources.py --self-test
python scripts/audit_eboss_dr16_official_mask_sources.py \
  --out-dir eboss_workspace/official_mask_inventory
```

The source-inventory script accepts only the published eBOSS DR16
catalogue tree on the official SDSS/NERSC hosts; it does not read
any observed or mock galaxy data. It parses the root index and
`ELGmasks/` index, reports the exact filenames of each fixed ELG
chunk family, checks for all three extra polygon files and records
candidate LRG-related root names plus root subdirectories.

```bash
python - <<'PY'
import json
from pathlib import Path
p=Path('eboss_workspace/official_mask_inventory/official_mask_index.json')
r=json.loads(p.read_text())
print('status:',r['status'])
print('root directories:',r.get('root_subdirectories',[]))
for chunk,item in r.get('elg_chunks',{}).items():
    print(chunk, 'listed brick maskbits:', item['listed_maskbits'])
for name,item in r.get('elg_extra_polygons',{}).items():
    print('extra veto:',name,'listed:',item['listed'])
print('LRG candidates:',r.get('lrg_root_candidates',[]))
print('errors:',r['errors'])
PY
```

On complete inventory the status is
`official_elg_chunk_index_discovered`. If some official file cannot
be located, it writes a machine-readable incomplete report and
exits nonzero rather than guessing filenames. No large maskbit FITS
family is downloaded automatically.

## Optional: byte-pin the three listed extra veto polygons

Only after the first index run has verified all three are listed:

```bash
python scripts/audit_eboss_dr16_official_mask_sources.py \
  --out-dir eboss_workspace/official_mask_inventory \
  --download-extra
```

This records the exact SHA256 and byte length of each extra polygon
in the JSON report. Re-running will fail if previously pinned bytes
change. These SHA values are *not* predetermined and must not be
fabricated. The status
`official_elg_extra_bytes_pinned_and_chunk_indexed`
means only that inputs are pinned, **not** that a joint mask has
been reconstructed.

## Still needed to certify the physical joint footprint

Identify the exact official LRG MANGLE sector/completeness and veto
files, including north/south conventions and all relevant veto
components; pin all their exact release paths and SHA256. Pin and
process all four ELG brickmask families and the upstream extra-mask
script with its bit semantics. Reproduce each tracer's independent
accept/reject rule, then evaluate their pair-level intersection
and weighted completeness separately in NGC and SGC against pinned
random-only inputs. Predeclare numeric tolerances and mask membership
checks before running them.

The `official_mask_index.json` file must never be treated as an
exact common selection mask, estimator closure, mock covariance or
observational significance calculation. The actual eBOSS odd vector
remains sealed until those separate checks pass.

## Verified LRG index and eleven-file polygon provenance

The official DR16 LRGandQuasarmasks/ directory contains seven published
polygons: `allsky_bright_star_mask_pix.ply`,
`badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply`,
`bright_object_mask_rykoff_pix.ply`, `brightstarmask_tiling_final.ply`,
`centerpost_mask_eboss_DR16_new.ply`,
`collision_priority_mask_QSO_eboss_DR16_v9_singletiles.ply`, and
`collision_priority_mask_lrg_eboss_DR16_new.ply`. The official root
also contains `eBOSS_QSOandLRG_fullfootprintgeometry_noveto.ply`.
These eight files and the three ELG extra polygons form the declared
eleven-file *byte-provenance* set. Do not infer that all seven veto
polygons should be indiscriminately unioned; the exact LRG selection
rule and sector completeness still require independent certification.

After obtaining `official_mask_index.json` as described above, run:

```bash
python scripts/audit_eboss_dr16_official_polygon_sha.py --self-test

python scripts/audit_eboss_dr16_official_polygon_sha.py \
  --inventory-json eboss_workspace/official_mask_inventory/official_mask_index.json \
  --out-dir eboss_workspace/official_mask_inventory \
  2>&1 | tee eboss_workspace/official_mask_inventory/polygon_sha.log
```

The runner verifies the three already pinned official directory SHA256
indexes, streams the eleven exact published polygons with SHA256 and
byte counts, rejects previously pinned byte changes, and records
`official_polygon_sha_manifest.json`. A complete status
`official_lrg_elg_polygon_bytes_pinned_only` does **not** certify
any continuous angular mask, ELG brickmask FITS selection, exact
LRG×ELG pair footprint or odd statistic. The four ELG brickmask
image families must still be processed separately.

## Resumable scientific provenance when the official LRG server times out

The original eleven-file GitHub job did **not** pass: the official
host timed out before completion. Do not claim SHA certification for
any LRG file whose JSON manifest was not written. The separate
official three-ELG-polygon job passed with pinned SHA256 records.

To isolate large LRG transfers, each of the eight exact published
LRG/QSO polygon files now has an independent GitHub job in
`eboss_dr16_lrg_polygon_shards.yml`, run `36102933048`. Failure of a
large file does not invalidate independently completed SHA manifest
artifacts for other filenames. Two transfers run concurrently to
limit source load. This is source-byte provenance only and not mask
application.

To reproduce a single file locally using the already downloaded
index, choose its exact published basename, for example:

```bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate

python scripts/audit_eboss_dr16_official_polygon_sha.py \
  --inventory-json eboss_workspace/official_mask_inventory/official_mask_index.json \
  --only-filename collision_priority_mask_lrg_eboss_DR16_new.ply \
  --out-dir eboss_workspace/official_mask_lrg/collision_priority_mask_lrg_eboss_DR16_new.ply \
  --timeout 125
```

Successful per-file status is
`official_single_polygon_sha_pinned_only`. Every downloaded polygon
must agree with the same pinned official source-directory indexes.
Transfer bytes are SHA256-hashed and checked against their advertised
HTTP Content-Length when available. An incomplete file is never
accepted as a certified mask input. Exact mask veto composition is a
separate scientific gate.
