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

## User's official full-LRG header result and frozen byte-only follow-up

The official header-only run completed locally and its JSON was uploaded.
The uploaded 4,565-byte JSON has SHA256
`901c5486e817da57a162d9bec74171a0fbc09f4e0901d72b87d7f4495ff51308`.
This exact report fingerprint is separately frozen in
`source_data/eboss_dr16_lrg_full_sector_header_uploaded_manifest_2026-09-25.json`
and required by the next local raw-byte-only runner. The original report
remains at `eboss_workspace/official_mask_inventory/lrg_full_sector_header_only.json`
in the user's WSL checkout.

The official LRG full FITS reports 196,485,120 bytes, 20,160 header
bytes (header SHA256
`a48b0d1055b5f7638dc77a86eb5ff5c861b629d19ff75cd8ae8122caa7e9014b`),
311,848 declared rows at 630 bytes each, and 72 columns. It contains
`SECTOR` (FITS J), `sector_SSR` (D), `sector_TSR` (D),
and `COMP_BOSS` (D). These column identities are now verified as
released **header structure** but their source-to-paper completeness
equivalences are not inferred from names alone. No galaxy table row,
sector value or odd statistic was read at this stage.

The frozen follow-up protocol is
`source_data/eboss_dr16_lrg_full_bytes_provenance_protocol_2026-09-25.json`.
Its runner `scripts/audit_eboss_dr16_lrg_full_bytes_provenance.py`
checks the precise previous uploaded report SHA, the earlier
SHA-pinned release index, and all eight official LRG/QSO MANGLE polygon
fingerprints, then streams the exact 196-MB official
`eBOSS_LRG_full_ALLdata-vDR16.fits` to a quarantined local folder.
It checks HTTP status 200, exact ETag, Last-Modified, byte count
and the previously frozen first-header SHA256. It writes the **first-seen
full-file SHA256 without interpreting a single FITS data row**. An
already-downloaded exact local copy can alternatively be passed by
`--existing-file /absolute/path/to/file.fits`.

Run in the same audit-branch checkout:

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python -u scripts/audit_eboss_dr16_lrg_full_bytes_provenance.py --self-test
python -u scripts/audit_eboss_dr16_lrg_full_bytes_provenance.py
```

[The synthetic byte-only CI passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36161956132).
After the local full-file checksum is returned, freeze it in a separate
sector-only protocol **before** reading any values of `SECTOR`,
`sector_TSR`, `sector_SSR` or `COMP_BOSS`. Those fields may then
be audited with aggregate sector-level output and official semantic
source references, without opening individual LRG redshifts, weights,
identifiers or the eBOSS odd vector.

Do not redownload the 11 SHA-pinned polygons. Do not automatically
apply all seven LRG/QSO veto names, insert new completeness cuts into
published clustering catalogues, or infer a single joint LRG/ELG mask.
The resulting science-pair window must remain the published
tracer-specific LRG×ELG random-pair selection.

## Completed LRG raw-byte SHA gate and separately frozen four-field audit

The user's actual local LRG raw-byte-only runner completed successfully.
The official `eBOSS_LRG_full_ALLdata-vDR16.fits` has
196,485,120 bytes and the first-seen complete-file SHA256
`39b831801adec04fe6dc5d6ab76a4b303aa58bfd303548cb7fddfae7d7f9331d`.
The local JSON status is
`OFFICIAL_LRG_FULL_BYTES_SHA_PINNED_ONLY`, with no FITS row or sector
field parsed and `OBSERVED_ODD_READ=False`. The user's exact uploaded
1,382-byte JSON was independently byte-hashed (SHA256
`a8d515ff27c9368f65f93ee1803ae30b94eca29d945820b91807d1063c7b1d0d`)
and archived unchanged as
`source_data/eboss_dr16_lrg_full_bytes_provenance_local_2026-09-25.json`;
its Git blob SHA1
`b6bd5db807db2f9bae892d1eb0e5356ab35eec42` was verified.

The full-file SHA and exact user-uploaded report bytes were registered
**before any LRG table-row sector data were read** in the separate
`source_data/eboss_dr16_lrg_sector_metadata_protocol_2026-09-25.json`.
The next audit code,
`scripts/audit_eboss_dr16_lrg_sector_metadata.py`, checks the entire
local FITS file's SHA256 and exact first-header SHA256 BEFORE
memory-mapping its data. It decodes ONLY `SECTOR` (FITS J),
`sector_TSR`, `sector_SSR`, `COMP_BOSS` (all FITS D) with
precise FITS strides. It does not decode galaxy RA/DEC, redshifts,
individual IDs, WEIGHT columns, randoms, mocks, any pair count or odd
statistic. It never prints individual sector labels or galaxy values.

All output is aggregate: total rows, number of distinct integer sector
labels, row-count distribution across sectors, per-column health and
distribution, within-sector numerical constancy under predeclared
absolute tolerance `1e-10`, and three FIXED candidate numerical
relations `COMP_BOSS-sector_TSR`,
`COMP_BOSS-sector_TSR*sector_SSR`, and
`sector_TSR-sector_SSR`. These are mathematical diagnostics ONLY,
not a data-driven choice of the published `C_eBOSS` or `C_z`
definition. In particular, no sector is excluded on the basis of
`0.5`, no MANGLE veto is newly applied, and no science-selection
rule changes.

The runner uses the same quarantined official LRG FITS downloaded in
the previous step; **no repeat download**:

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python -u scripts/audit_eboss_dr16_lrg_sector_metadata.py --self-test
python -u scripts/audit_eboss_dr16_lrg_sector_metadata.py
```

The synthetic CI is
https://github.com/dvlahek/stress-energy-closure/actions/runs/36169472480.
Passage of the local four-field audit will close an official catalogue
sector-*metadata* consistency check, not the independent source
interpretation of sector completeness or the historical production
MANGLE polygon-composition contract. The released LRG/ELG
tracer-specific random pair window, complete ELG brickmask and
mock-galaxy estimator/covariance gates stay open. Observed eBOSS odd
data remain unopened.

## Actual local four-column sector aggregate result — semantics still open

The user ran the exact preregistered local metadata runner successfully on
2026-09-25. The original full-file SHA256 was registered before this
row-level inspection. The reported output status is
`OFFICIAL_LRG_SECTOR_METADATA_AGGREGATES_ONLY_NOT_SEMANTICS_CERTIFIED`.
This is currently **user-reported terminal stdout**, captured separately
in `source_data/eboss_dr16_lrg_sector_metadata_reported_outcome_2026-09-25.json`.
The full JSON report and its exact SHA256 are **not yet uploaded**.
Do not claim to have independently regenerated the row-level result.

The published full-LRG catalogue contains 311,848 inspected rows and
5,813 distinct integer `SECTOR` labels. At the preregistered
absolute tolerance `1e-10`, all three numerical fields
(`sector_TSR`, `sector_SSR`, `COMP_BOSS`) are **exactly
constant** within each sector (zero within-sector span everywhere).
All their reported numerical values are finite and in [0,1].
The following numbers count **catalogue rows**, not independent
sectors or sky area:

| Diagnostic field | Exact-zero rows | Rows with field <= 0.5 | Median |
|---|---:|---:|---:|
| `COMP_BOSS` | 9,413 | 58,575 | 0.9855072463768116 |
| `sector_TSR` | 9,413 | 59,250 | 0.9421487603305785 |
| `sector_SSR` | 9,424 | 9,528 | 0.9855072463768116 |

The three **predeclared** formula comparisons also show that the
published fields cannot simply be substituted for one another:
`COMP_BOSS = sector_TSR` agrees to `1e-10` on only 107,090 rows,
`COMP_BOSS = sector_TSR * sector_SSR` on 74,143 rows, and
`sector_TSR = sector_SSR` on 25,676 rows. A failed algebraic
identity is evidence against that identity on this catalogue, **not**
evidence for an alternative formula chosen after seeing the data.
Do not declare any of the logged fields to be paper `C_eBOSS` or
`C_z` on numeric correlation alone.

Ross et al. (2020), MNRAS 498, 2354, Sec. 5.4, Eqs. (10)–(11),
give independent **source-defined** target-observation completeness
`C_eBOSS` and redshift-success completeness `C_z` in terms of
per-sector categories `N_z,eboss`, `N_cp`, `N_badclass`,
`N_star`, `N_zfail`, and `N_missed`. Their documented
clustering threshold is `C_eBOSS > 0.5` and `C_z > 0.5`;
random positions are completeness-weighted instead of being
subsampled by `C_eBOSS`. See
https://academic.oup.com/mnras/article/498/2/2354/5900562 .
The full-catalogue four-field audit does not contain the underlying
per-sector category counts or an authenticated production-code mapping
from these equations to the release's `COMP_BOSS`,
`sector_TSR` or `sector_SSR` fields. An official data-model
or source-code provenance check is still required, and must be
distinguished from observing numerical coincidences.

**Next provenance action:** the user should upload the exact local
`eboss_workspace/official_mask_inventory/lrg_full_sector_metadata_audit.json`
or supply its raw SHA256. Archive the JSON unchanged and freeze its
content digest. Only then commit and run any further separate
sector-only comparison. No new veto union, sky mask, catalogue cut,
LRG×ELG random selection, or observed odd measurement was carried
out in this metadata diagnostic.
