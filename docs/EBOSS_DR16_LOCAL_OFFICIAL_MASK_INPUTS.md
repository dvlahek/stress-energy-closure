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
