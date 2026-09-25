# eBOSS DR16 LRG sector/veto provenance: next local header-only gate

Status (2026-09-25): **published selection categories identified; exact LRG veto composition, sector weights and physical mask not yet certified.** This follows the successful blinded official ELG bit-8 label audit. The exact 269,178-row official ELG full-catalogue report was archived byte-identically in `source_data/eboss_dr16_elg_full_bit8_label_geometry_audit_2026-09-25.json`. The native-RA source-pixel mapping reproduces all 15 official bit-8 positive labels with zero FP/FN; no science selection changed.

## What the published LRG studies establish

Ross et al. (2020), MNRAS 498, 2354, Sections 5.2 and 5.4, Table 2
(https://academic.oup.com/mnras/article/498/2/2354/5900562)
describe an LRG footprint formed from sector-level survey geometry, tracer-specific vetoes and completeness. The paper distinguishes LRG collision-priority veto from the QSO collision-priority veto. It also lists LRG bad-field, bright-star, infrared-bright-star, bright-object and centrepost veto categories. Their veto areas overlap; summing nominal areas or blindly unioning all seven published LRG/QSO polygons is not a validated selection algorithm.

For the released LRG clustering construction, the published cuts are
`C_eBOSS > 0.5` and `C_z > 0.5`, with strict inequality. Ross et al.
define `C_eBOSS` per sector and describe completeness as a random
*weight*, not a new deterministic subsampling of released random
positions. Zhao et al. (2021), MNRAS 503, 1149, Section 2.3.2
(https://academic.oup.com/mnras/article/503/1/1149/6149160)
describe the corresponding LRG/QSO sector trimming in matched
EZmocks. Their ELG prescription is different. These literature
statements **do not** by themselves identify the correct FITS
`COMP_BOSS`, `sector_TSR` or `sector_SSR` mappings for this exact
release, and must not be re-applied as new cuts to published clustering
data or randoms.

The existing official SHA-pinned eleven-polygon manifest is
`source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json`.
Six published LRG/QSO polygon filenames are *candidates* for the six
LRG veto classes; their exact optical-versus-infrared bright-star
mapping and composition have not been authenticated against the
production pipeline. The seventh filename,
`collision_priority_mask_QSO_eboss_DR16_v9_singletiles.ply`,
is QSO-specific in the published categorization, not an automatic
LRG veto. The eighth LRG/QSO file is the no-veto footprint polygon.
All these published bytes were previously SHA256-pinned; this next
gate does **not** download polygons again.

## One local test: only the official LRG full-catalogue FITS header

The source and expected published categories are frozen in
`source_data/eboss_dr16_lrg_sector_header_protocol_2026-09-25.json`.
The runner `scripts/audit_eboss_dr16_lrg_full_sector_header.py`
requires the user's earlier SHA-matched official release index and
the repository's complete 11-file polygon provenance. It requests only
exact 2880-byte HTTP 206 FITS header blocks from the single official
`eBOSS_LRG_full_ALLdata-vDR16.fits`, and stops at the first BINTABLE
END card. If the server responds with HTTP 200, or the source URL,
Content-Range or header contract differs, it fails closed before
requesting table data.

The runner records the official column names relevant to sectors and
completeness, including available `sector*` and `COMP*` columns,
their FITS types, file ETag, source size and the bounded header SHA256.
It **does not read a single LRG data row**, inspect a redshift or
weight value, construct a mask, or evaluate a pair count.

Run in the *same audit branch* in the user's WSL checkout:

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate

echo "=== LRG FULL SECTOR HEADER ==="
python -u scripts/audit_eboss_dr16_lrg_full_sector_header.py --self-test
python -u scripts/audit_eboss_dr16_lrg_full_sector_header.py
echo "=== KRAJ ==="
```

The synthetic CI test checks both the immutable official polygon
manifest and the header-only, no-row-access behavior. A passed real
header run closes **only** the official FITS-column source inventory.
Next, pin the complete exact LRG file SHA and independently freeze a
sector-only metadata inspection protocol before reading any rows, then
authenticate published veto/sector composition against production
code or other official release products.

The published tracer-specific randoms still define the cross-LS
random-pair window `R_L R_E`. No common hard angular cut is inferred
from a coarse occupancy grid, no veto is newly applied to released
data/randoms, and the observed eBOSS odd-sector vector stays unopened.
