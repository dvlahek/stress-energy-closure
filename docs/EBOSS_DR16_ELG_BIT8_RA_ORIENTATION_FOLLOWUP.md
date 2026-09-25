# eBOSS DR16 ELG bit 8: independent RA-orientation follow-up

Status (2026-09-25): **source-only geometry checked; official production bit-8 convention NOT identified.** This is separate from the existing source/API unit discrepancy and does not certify a physical mask, reapply any veto, open observed odd data or modify the DESI results.

## Problem and limitations of the previous roundtrip

The byte-pinned upstream `cheng-zhao/brickmask` script at commit `b9eb684a579b56ec3dbdb46549224be7e3fa2830` forms
`theta=radians(90-dec)` and `phi=radians(360-ra)`, then calls `healpy.ang2pix(1024,theta,phi,nest=False,lonlat=True)`.
That call passes radians to an option that expects degrees. The previous frozen source-only test recovered all 37 listed pixel indices using `lonlat=False` after constructing test RA from `360-phi`, but that demonstrates the internal roundtrip of the *reflected-RA* transformation, **not** its compatibility with celestial RA in the public survey footprint.

A separate source-only follow-up is frozen in
`source_data/eboss_dr16_elg_bit8_ra_orientation_protocol_2026-09-25.json`.
It records the exploratory arithmetic that motivated the follow-up **before** the CI replay. The pinned upstream script is read and AST-parsed to extract the entire list of 37 source pixels, and its existing Git blob and integer-array SHA256 are verified. The follow-up never opens a random or galaxy catalogue.

## Independent native-RA versus reflected-RA check

For each published HEALPix pixel, derive its centre with
`hp.pix2ang(1024, pixel, nest=False, lonlat=True)`. The returned
longitude is used as *native HEALPix longitude*; treating that as celestial RA is a separate, explicitly testable hypothesis. At these same 37 native-longitude centres, evaluate four alternatives: the upstream call as written, the units-only repair with `360-RA`, native-RA spherical radians, and native-RA longitude/latitude degrees. Native-RA radians and native-RA degrees agree and recover all source indices. The reflected-RA units-only mapping returns no pixel in the 37-pixel source list at these native-longitude centres.

As a **diagnostic only**, the published approximate ELG NGC rectangle
`126<RA<166`, `13.8<DEC<32.5` and SGC rectangle
`-43<RA<45`, `-5<DEC<5` provide an external sky-location check.
With native HEALPix longitude as RA, 10 listed pixel centres are in the approximate NGC rectangle, 26 in SGC and one is just outside the rounded NGC declination limit (DEC approximately 32.531 degrees). With reflected RA, zero centres are in the approximate NGC rectangle, 26 in SGC and 11 outside. These rectangles are **not** the exact DR16 footprint; their counts must never define a physical cut, accepted area or pair window.

The 37 equal-area nside-1024 pixels have a combined geometric area of approximately `0.121304 deg^2`, compatible with the published one-decimal `0.1 deg^2`. Raichoor et al. (MNRAS 500, 3254–3274, 2021, doi:10.1093/mnras/staa3336), Sec. 3.2 and Table 5, describe the rejection of 37 pixels for bit 8 and list 15 removed targets. We have not independently checked those 15 catalogue targets.

Code: `scripts/audit_eboss_dr16_elg_bit8_ra_orientation.py`.
The frozen synthetic self-test and all-pixel source-only CI replay passed at commit
`5f3d88e7a9b9ba46ad0b78e627d92435477fbf4e`:
https://github.com/dvlahek/stress-energy-closure/actions/runs/36151465886

## What this establishes and what remains open

The existing script has a documented radian/`lonlat=True` API discrepancy. Its `360-RA` transformation also requires an independent physical-coordinate check: the former 37/37 reflected roundtrip is insufficient to validate the actual sky orientation of bit 8. These are **source-code and source-list** checks. They do not show how the official released DR16 clustering catalogues or published production mask were generated, so do not claim that those catalogues were incorrectly masked.

The independent reference gate requires an exact-version official production bitmap, or an independently verified **pre-veto bit-coded** output with coordinates, bit-8 flags, source provenance and SHA256. Test both positive and negative positions that discriminate the literal, reflected-RA and native-RA conventions; include released-mask and chunk handling only after the production reference is authenticated. An aggregate area, the number 15, or post-veto clustering catalogue survivors alone cannot distinguish candidate algorithms. If the official production reference is unavailable, record that unresolved limitation; do **not** manufacture it or silently replace the upstream call.

The LRG MANGLE veto/sector rules and four ELG brickmask FITS families can be audited in parallel with this reference search. The published LRG and ELG clustering randoms retain their own selection functions, and the cross-LS pair window remains the corresponding LRG×ELG `R_L R_E`; do not invent a shared mask from coarse random occupancy.

The observed eBOSS odd vector remains unopened. Random-only LS closure, 11/11 polygon SHA provenance and the 36/36 fine weighted n(z) check are already complete and do not need repetition.

## Next local WSL gate: official *reference filename* discovery only

The separately frozen protocol is
`source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json`.
The runner `scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py` reads
only the official **DR16 root HTML index** using a 16-MiB bound. It checks
the previously recorded exact SHA256 of that index before accepting any
filenames. It does not repeat any polygon downloads or FITS parsing.

Run in the user's existing WSL checkout of this audit branch:

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py --self-test
python scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py \
  --out eboss_workspace/official_mask_inventory/bit8_reference_index.json
```

A successful index match means **filenames only**, not that any independent
bit-8 reference was found. A changed official HTML listing fails closed and
must be reviewed separately; do not silently substitute its SHA. Inspect the
reported candidates for a genuinely bit-coded *pre-veto* official production
output. Post-veto `eBOSS_ELG_clustering_*` files and the ordinary SDSS
`EBOSS_TARGET1` ELG selection flags cannot certify the published ELG
BRICKMASK extra bit 8.

The networking-free synthetic self-test passed in
[workflow 36152677996](https://github.com/dvlahek/stress-energy-closure/actions/runs/36152677996).
The official root listing has **not** been retrieved by that CI test.
No reference bitmap was authenticated, and no physical mask/odd-data gate
was opened.

## Next local gate: single ELG full-catalogue FITS header only

The SHA-matched official release index lists the candidate
`eBOSS_ELG_full_ALLdata-vDR16.fits`. An earlier interpretation labelled `full_ALLdata` as **post-mask** without
verifying its actual row-level retention. That inference is withdrawn.
Raichoor et al. (2021, MNRAS 500, 3254, Table 5) explicitly state that
all angular veto masks except the two bad eboss22 plates are bit-coded
in the catalogues' `mskbit` column. This does NOT establish that this
particular `full_ALLdata` file retains bit-8-positive rows. Its header may
reveal `MSKBIT` or cross-reference columns but cannot identify the
production bit-8 convention until provenance and separately frozen
mask-only row inspection are completed.

The frozen protocol is
`source_data/eboss_dr16_elg_full_header_only_protocol_2026-09-25.json`.
Its runner
`scripts/audit_eboss_dr16_elg_full_header_only.py`
requires the prior local index report, verifies the pinned source URL,
and requests **only exact 2880-byte HTTP 206 FITS header blocks** from
the official NERSC mirror. It stops after the primary and first
BINTABLE `END` cards, before the first data row. An HTTP 200
full-body response, an unexpected redirect or a missing range header
fails closed. No whole-file FITS SHA is claimed from header bytes.

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python -u scripts/audit_eboss_dr16_elg_full_header_only.py --self-test
python -u scripts/audit_eboss_dr16_elg_full_header_only.py \
  --index-report eboss_workspace/official_mask_inventory/bit8_reference_index.json \
  --out eboss_workspace/official_mask_inventory/elg_full_header_only.json
```

The synthetic CI test checks that no header parser request crosses the
start of the first data row. The real official header metadata remains
to be checked locally. In every case the independent pre-veto
production reference gate stays **unresolved** until an authenticated
reference with discriminating positive and negative bit-8 locations
is available. The observed eBOSS odd sector remains blinded.

## Actual WSL header result and next byte-provenance-only gate

The user's local official HTML index matched the frozen root SHA and
the official FITS header-only runner completed. Its report
`eboss_workspace/official_mask_inventory/elg_full_header_only.json`
records HTTP header metadata: 188,985,600 remote bytes, 20,160 requested
header bytes, BINTABLE declaration of 269,178 rows, 702 bytes/row,
85 columns and an `mskbit` column (FITS TFORM `I`). No table rows were
read. The `mskbit` field makes the full catalogue a **candidate** for
an independent production mask-reference check, not yet a verified
pre-veto reference.

Before looking at any row, use the separately frozen byte-only protocol
`source_data/eboss_dr16_elg_full_bytes_provenance_protocol_2026-09-25.json`.
The runner
`scripts/audit_eboss_dr16_elg_full_bytes_provenance.py`
accepts only the same exact official source, HTTP ETag, Last-Modified,
file length and first-20,160-byte header SHA. It streams the raw file
into a quarantined local folder, computes whole-file SHA256 and records
that first-seen checksum. It does **not** interpret a FITS column or row.
An already-downloaded local file can instead be verified with
`--existing-file /absolute/path/to/file.fits`, without copying it.

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python -u scripts/audit_eboss_dr16_elg_full_bytes_provenance.py --self-test
python -u scripts/audit_eboss_dr16_elg_full_bytes_provenance.py
```

The public source file is approximately 189 MB. Preserve the local
provenance JSON with its whole-file SHA256. It is a hash observed and
recorded **before** any row audit, not an independently published SDSS
checksum. Freeze the hash in a separate follow-up protocol, then define
a strict isolated audit of only `RA`, `DEC`, `mskbit`; do not read
redshifts, galaxy weights or compute any odd/cross-correlation.
Production bit 8 remains unresolved until positives, negatives and
coordinate-contract tests agree. Main/PR retain the previous frozen
selection and the observed odd vector stays unopened.

## Independently frozen catalogue-byte SHA and isolated bit8 label audit

The user's local byte-only run finished successfully on 2026-09-25:
`188985600` official ELG full-catalogue bytes, whole-file SHA256
`8806699e14422904ef91efb6f1632184171fb740c2e074456c0534b7d2dc28b3`,
previously pinned first-header SHA256
`138e9ef0ba23b0cc934c087e44564229822944579903223d41d3823cb737657f`.
Its `OBSERVED_ODD_READ` status is `False`. This is a local byte
fingerprint established **before** any FITS table column was decoded.
The corresponding report remains in the user's local
`eboss_workspace/official_mask_inventory/elg_full_bytes_provenance.json`.
The byte-only runner has not inferred a catalogue label.

Those values are now immutably specified for the next test in
`source_data/eboss_dr16_elg_full_bit8_label_protocol_2026-09-25.json`.
The test runner `scripts/audit_eboss_dr16_elg_full_bit8_labels.py`
rechecks the entire local file's SHA256 and the exact prior header SHA
*before reading any data row*. It decodes exactly three memory-mapped,
FITS big-endian, strided observed columns: `RA`, `DEC`, `mskbit`.
No other table field is decoded. The 37 upstream pixels, their script
Git blob and the four source-only candidate mappings are verified against
the earlier frozen protocol. Output comprises only aggregate labelled
bit8-positive count and fixed all-row TP/FP/FN/TN counts for the four
predeclared coordinate interpretations.

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python -u scripts/audit_eboss_dr16_elg_full_bit8_labels.py --self-test
python -u scripts/audit_eboss_dr16_elg_full_bit8_labels.py
```

The synthetic three-column, byte-stride and source-geometry
[CI test passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36159868366).
The full observed catalogue label audit has **not** run in CI; it is
for the user's existing WSL checkout and exact quarantined official
file. The full catalogue might have zero or a different number of
bit8-positive rows; do **not** assume the paper's 15 removed targets
are present. That reported number is only a published comparison.

Even exact per-row correspondence is limited to bit8 labels on the
sampled released catalogue. It does not independently authenticate
which production executable generated the labels, certify masks
between observed positions, provide an ELG×LRG pair window, or authorize
the odd-sector analysis. No alternative mask is applied; the PR
remains draft and `main` unchanged.
