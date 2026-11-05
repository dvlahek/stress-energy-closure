# eBOSS DR16 ELG 19,381-maskbit four-family source-index gate

**Status (2026-09-26): source filenames SHA-pinned in a separate protocol; real FITS image bytes and continuous ELG selection NOT certified.** This follows the LRG two-stage Table 2 source-informed consistency replay. No LRG/ELG science catalogue, random, pair-count or observed odd-sector data are inspected here.

## Exact uploaded official index

The user uploaded an 803,716-byte `official_mask_index.json` (the chat attachment is renamed `official_mask_index(1).json`). Independently verified SHA256:
`84a3a163ec0b28a5392dadc95aa757e469cb9cb6e94cc731aac2fa2958fc6c18`.
This SHA pins the uploaded local index as a whole; its underlying original NERSC directory-listing SHA256 fingerprints are separately pinned in `source_data/eboss_dr16_official_mask_source_protocol_2026-09-25.json`.

The uploaded index includes four **ordered** ELG `mask-<chunk>-<brick>.fits.gz` filename families: `eboss21`: **2,904**; `eboss22`: **7,315**; `eboss23`: **6,159**; `eboss25`: **3,003**. Total **19,381**. All filenames match their exact family and brick naming convention, are lexicographically sorted and unique, with no cross-family duplicates. All three official extra veto polygon filenames are listed. Their SHA256 **file-byte** fingerprints were already recorded in `source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json`; the uploaded filename index itself does not include those file-byte SHA256 values.

The exact full-index SHA256, each chunk's independently computed newline-delimited filename-list SHA256, counts, ordered upstream `SUBSAMPLE_ID` values `0,1,2,3`, pinned official source-index SHA256s and already archived extra-polygon file SHA256s are frozen in `source_data/eboss_dr16_elg_maskbit_index_pin_2026-09-26.json`. Pinning a filename list **does not** pin the 19,381 FITS gzip image bytes.

## Offline next gate

`scripts/audit_eboss_dr16_elg_maskbit_index.py` verifies the local index's **exact uploaded raw SHA256**, parent official directory index fingerprints, all four full filename lists and digests, pinned upstream brickmask commit and `-DEBOSS` configuration, and the three previously pinned extra polygons. It reads no FITS images, target or random rows, pixel bits, or observed odd data. On a mismatch it writes an explicit incomplete report and fails closed. Its result is a compact JSON:
`eboss_workspace/official_mask_inventory/elg_maskbit_index_pin.json`.

Run in the **draft audit branch** from WSL (the index is already present locally; no file download):

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_elg_maskbit_index.py --self-test && python -u scripts/audit_eboss_dr16_elg_maskbit_index.py`

No repeat polygon download or bulk download of 19,381 FITS files is authorized by this gate. The synthetic CI is
https://github.com/dvlahek/stress-energy-closure/actions/runs/36219621893 .

## Subsequent predeclared physics-input checks

The pinned upstream `cheng-zhao/brickmask` commit `b9eb684a579b56ec3dbdb46549224be7e3fa2830`, `CONFIG.md`, and `io/read_fits.c` explicitly use the four eBOSS chunks in that order; the reader requires a two-dimensional first image HDU and `RA---TAN`/`DEC--TAN` WCS. This is a **source contract**, not evidence that each particular published gzip file has been inspected.

After the full filename-list gate passes locally and its exact output JSON is archived, freeze a separate test of one preselected FITS image per chunk, with bounded transfers, exact SHA256, FITS-header/WCS validation, and no pixel-level selection or observed science data. Extend toward all files only with an explicit resource budget and published bitmask semantics. Existing bit8 label geometry and three extra-polygon bytes are distinct checks; neither certifies the complete ELG BRICKMASK production path. The eventual LRG×ELG window must be validated using tracer-specific published clustering randoms, not a fabricated common hard sky mask.

## Exact offline report archive and preregistered four-sample raw gzip SHA gate

The user's *actual local* offline four-family index check passed with
`OFFICIAL_ELG_FOUR_FAMILY_FILENAME_INDEX_SHA_PINNED_ONLY`,
19,381 files across the frozen four chunk lists and
`ODD_DATA_READ=False`. The exact uploaded report is 2,789 bytes,
SHA256 `32104c056ab86637df1d236bdde04263ab18393b66302b690d70ce331102f2b7`;
it was archived unchanged as
`source_data/eboss_dr16_elg_maskbit_index_local_report_2026-09-26.json`
and verified by the matching Git blob
`ca9555978b4198515ebd3f907cfbc72f859f5aa1`.
The accompanying manifest is
`source_data/eboss_dr16_elg_maskbit_index_uploaded_manifest_2026-09-26.json`.

Before opening **any individual mask FITS.gz image bytes**, a separate
bounded source-byte-only protocol was frozen in
`source_data/eboss_dr16_elg_maskbit_four_sample_raw_sha_protocol_2026-09-26.json`.
Its rule is the exact middle element of each SHA-pinned sorted filename
list, not a replacement picked after observing outcomes. The four
preselected official products are
`mask-eboss21-3386m005.fits.gz`,
`mask-eboss22-0226m005.fits.gz`,
`mask-eboss23-1415p232.fits.gz`, and
`mask-eboss25-1538p312.fits.gz`.

The runner `scripts/audit_eboss_dr16_elg_four_sample_raw_sha.py`
requires byte-identical archived and local previous index reports, the
full original 803,716-byte index SHA256, the four exact filename lists
and their published directory provenance before downloading. It streams
only these four HTTP 200 official gzip file responses, requiring exact
bounded `Content-Length`, unchanged source URL and gzip magic.
Compressed bytes are hashed (SHA256) with a maximum of 64 MiB per
selected file and 256 MiB combined. Completed files are quarantined,
with a resumable exact-checksum checkpoint; interrupted `.part` files
are never treated as certified source input. The four first-seen
individual full gzip SHA256s are **not** SDSS-published checksums and
must be frozen in a subsequent protocol before even FITS-header/WCS
inspection. This first step does **not decompress gzip or read a single
FITS header or pixel**, and cannot prove a physical ELG mask.

The synthetic-only CI workflow is
https://github.com/dvlahek/stress-energy-closure/actions/runs/36222451283.
After CI passes, execute one WSL line:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_elg_four_sample_raw_sha.py --self-test && python -u scripts/audit_eboss_dr16_elg_four_sample_raw_sha.py`

The local report will be
`eboss_workspace/official_mask_inventory/elg_maskbit_four_sample_raw_sha.json`.
Upload the JSON unchanged after the run. Keep the PR draft and
observed eBOSS odd-vector sealed. Do not treat these four
sampling checks as independent validation of all 19,381 images.

## Completed actual four-gzip byte gate and separately frozen FITS/WCS-only gate

The local source-byte runner completed for all **four preregistered
middle-of-list images**. The user's exact 4,338-byte report has SHA256
`f324da9232a121e72e9d6d35b7d4e7c2e0f7133007a01db8d065d4ba84356665`;
it is archived byte-identically as
`source_data/eboss_dr16_elg_four_sample_raw_sha_report_2026-09-26.json`
(verified Git blob
`9c543517be72e1622390a7164b6ab003a1ae7045`).
The uploaded-report identity is independently pinned in
`source_data/eboss_dr16_elg_four_sample_raw_sha_uploaded_manifest_2026-09-26.json`.
The archived report documents raw full-gzip SHA256, exact official
URL, HTTP 200, content length and source validators for each sample:

| Chunk | Exact selected gzip source | Compressed bytes | First-seen compressed SHA256 |
|---|---|---:|---|
| eboss21 | mask-eboss21-3386m005.fits.gz | 172308 | 064e535aba4caa7351f5ec64802421968724c1ea664b7eafc142a018a13049a9 |
| eboss22 | mask-eboss22-0226m005.fits.gz | 188582 | bec27e6ce17fba5336a48f87682ce491abcd603d1599d914d333dfacc8c2d0d8 |
| eboss23 | mask-eboss23-1415p232.fits.gz | 139821 | faf180362b182aa9fb56221af4848462ff3ba801a437a4e1ceb2a43cdb2f569e |
| eboss25 | mask-eboss25-1538p312.fits.gz | 175454 | ce2b3dc302a3998b1ac4e061989af7916ded6221a80a3643724063e6ed19fc27 |

Combined raw compressed size: **676165 bytes**. This source-byte
check did **not** decompress any sample, read a FITS header/pixel, inspect
galaxies/randoms or open the observed odd signal. The four file SHA256s
are first-seen independent local checksums, **not published SDSS
checksums** and do not certify the remaining 19377 files.

The exact four checksums were committed **before any FITS
decompression/header inspection** in a separate
`source_data/eboss_dr16_elg_four_sample_fits_header_wcs_protocol_2026-09-26.json`.
The new strictly header-level
`scripts/audit_eboss_dr16_elg_four_sample_fits_header_wcs.py`
first independently hashes all four existing compressed files and
checks the exact unchanged earlier local and archived JSON reports.
It opens each local gzip stream only after all four raw SHA checks
pass, requests at most 57600 FITS header bytes per image as
2880-byte blocks and validates first image-HDU metadata:
two-dimensional integer image with `RA---TAN`/`DEC--TAN` WCS,
CRVAL/CRPIX/CD matrix entries present and finite, and nonsingular
CD. An empty primary followed immediately by an IMAGE extension is
also allowed. It never requests image pixel values, applies mask
bit codes, reads any catalogue rows or uses odd data.
The gzip implementation may buffer/decompress some subsequent
bytes internally, so this is not a claim that absolutely no pixels
are decompressed inside a library or that the entire gzip CRC
was verified. Only the FITS header bytes are returned to the
analysis code.

The synthetic CI is
https://github.com/dvlahek/stress-energy-closure/actions/runs/36222762663 .
After it passes, the local command (no download) is:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_elg_four_sample_fits_header_wcs.py --self-test && python -u scripts/audit_eboss_dr16_elg_four_sample_fits_header_wcs.py`

The resulting JSON is
`eboss_workspace/official_mask_inventory/elg_four_sample_fits_header_wcs.json`.
Upload it unchanged to freeze its exact SHA before any pixel-level
checks. This small per-family sample is a source/WCS structural check,
not an all-image mask validation or a LRG×ELG physical-window
certificate.

## Complete local WCS-header result and frozen four-image pixel-byte QA

The user completed the strictly bounded WCS audit on all four
preselected SHA-pinned FITS.gz samples. Each is a primary 8-bit
`3600 x 3600` image with `RA---TAN`/`DEC--TAN`
WCS, `CRPIX=(1800.5,1800.5)`, a nonsingular CD matrix,
and **5760** FITS header bytes. The compressed SHA256s
were reverified before any gzip opening. The user's exact
5,705-byte `elg_four_sample_fits_header_wcs.json` has
SHA256 `233339946dda20dc7b81f83c70645761dfc728f562396fecd261750e3eec6e04`,
archived byte-identically as
`source_data/eboss_dr16_elg_four_sample_fits_header_wcs_report_2026-09-26.json`
and independently verified against Git blob
`b0f515d9bdb5bb0c7ccde7b64d982d33c9e10361`.
The upload fingerprint and four header-block SHA256s are
independently recorded in
`source_data/eboss_dr16_elg_four_sample_fits_header_wcs_uploaded_manifest_2026-09-26.json`.
No image pixel values, galaxy/random rows, new selection or
observed odd statistic were read in this completed header check.

The pinned upstream `cheng-zhao/brickmask` source at commit
`b9eb684a579b56ec3dbdb46549224be7e3fa2830`
defines `EBOSS_MASK_VALID(bit)` as `bit & 1`,
`EBOSS_XYBUG_BIT=4` and `EBOSS_XYBUG_VALID(bit)` as
`bit & 4` in `src/define.h` (blob
`a6b85bfbc7745388a076fba822a0dc27d52d9809`).
The eBOSS path of `src/bit_code.c` (blob
`0ab07730c491324fac3222ebac1ae1c763efe980`)
checks the rounded pixel for bit0 and handles bit2 by comparing
rounded and truncated pixel locations. This is a documented
production *source-code rule*; its historical invocation with the
specific released catalogue and any coordinate boundary
conventions are not yet independently authenticated.

A distinct four-image **pixel aggregate** protocol was frozen
*after source checks but before any image pixel value was decoded*:
`source_data/eboss_dr16_elg_four_sample_pixel_bitcode_protocol_2026-09-26.json`.
The runner
`scripts/audit_eboss_dr16_elg_four_sample_pixel_bitcode.py`
checks exact archived+local WCS report bytes and all four
full compressed SHA256s before opening gzip. It checks each
already SHA-pinned 5760-byte FITS header again, then streams
the full 12,960,000 uint8 primary-image pixel values per sample,
reporting ONLY aggregate 256-bin byte histograms, bit0/bit2
counts, decompressed-image SHA256, and verified gzip EOF/CRC.
No pixel coordinates, reconstructed sky maps, galaxy/random
locations, science cuts, or odd statistics are produced.
Reading the full sample image is a NEW, narrowly permitted
source-input QA stage and must never be confused with the
earlier header-only stage. These FOUR samples do not establish
full 19,381-image source coverage or a physical ELG mask.

The synthetic-only CI is
https://github.com/dvlahek/stress-energy-closure/actions/runs/36223326539 .
After it passes, run in the same WSL audit branch without
any new FITS download:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_elg_four_sample_pixel_bitcode.py --self-test && python -u scripts/audit_eboss_dr16_elg_four_sample_pixel_bitcode.py`

Expected report:
`eboss_workspace/official_mask_inventory/elg_four_sample_pixel_bitcode_aggregates.json`.
The exact output must be SHA-pinned separately after upload.
The prerequisite for physical ELG×LRG pair analysis remains
source-authenticated tracer-specific catalogue/random mask
production; neither these four histograms nor a hard
common sky-mask intersection is sufficient.

## Completed four-image pixel-bitcode source QA and interpretation boundary

The user's **actual local** source-input pixel aggregate audit completed
with `OFFICIAL_ELG_FOUR_SAMPLE_PIXEL_AGGREGATES_SOURCE_QA_ONLY`,
after it reverified all four exact previously SHA-pinned compressed
source files and all four FITS header SHAs. The exact uploaded
16,601-byte `elg_four_sample_pixel_bitcode_aggregates.json` was
independently SHA256-verified
(`ed83fcaad2b299a987e1bc759f0156dbbcfb974242ecdedbc0fdb894534f0e7c`)
and archived byte-identically as
`source_data/eboss_dr16_elg_four_sample_pixel_bitcode_aggregates_2026-09-26.json`.
Its Git blob SHA1
`26a0851516543da1449694ecde022185c112296e`
was independently verified against the uploaded bytes; the separate
file-fingerprint and first-seen results manifest is
`source_data/eboss_dr16_elg_four_sample_pixel_bitcode_uploaded_manifest_2026-09-26.json`.

All 256 histogram bins in each of the four images add up to exactly
12,960,000 pixels. The sample's bit0/bit2 totals can be independently
reconstructed from the uploaded histograms:

| Exact fixed sample chunk | bit0 set | bit2 set | both bits set | Nonzero byte codes |
|---|---:|---:|---:|---:|
| eboss21 | 12,960,000 | 1,345 | 1,345 | 12 |
| eboss22 | 12,960,000 | 3,542 | 3,542 | 28 |
| eboss23 | 12,960,000 | 2,080,939 | 2,080,939 | 22 |
| eboss25 | 12,960,000 | 0 | 0 | 9 |

The complete gzip CRC/ISIZE and EOF check passed for each sample.
Every image pixel in the four chosen images has bit0 set. This is a
**source-image histogram diagnostic only**. In upstream
`src/bit_code.c` under `-DEBOSS`, a galaxy/random coordinate first
maps to a **rounded** WCS image pixel, the bit0 test controls whether
that pixel contributes, and bit2 (mask value 4) is adjusted using
the **truncated** pixel index. Because this audit deliberately does
not output pixel positions, read target/random coordinates, or test
coordinate-to-pixel mappings, even the comparatively large
2,080,939 bit2-pixel count in the selected eboss23 image is
**not itself evidence of an anomalous released mask or galaxy
selection**. Do not invent an adaptive new brick sample or infer the
historical selection rule from this count.

These four files are four **predeclared samples**, not source-byte or
semantic certification of the remaining 19,377 images. A full
published ELG mask requires the actual historical
code/configuration of target and random acceptance, the extra
bit8 RA mapping (still unresolved by the existing independent
reference gate), all 3 extra polygon rules (bits 9–11), and the
tracer-specific released random window. No new cuts, full mask,
pair estimator or observed odd-sector vector are permitted from
these four histograms. The next substantive work should close
historical production/reference provenance, not continue
post-hoc image sampling.

## Published, post-result xybug qualitative cross-check

Raichoor et al., *MNRAS* **500** (2021) 3254–3274, Sec. 3.2,
Table 5, independently describe **bit 2 (xybug)** as affecting
`eboss23` most strongly, `eboss21` and `eboss22` to a lesser
extent, with `eboss25` unaffected. Table 5 gives 49.7 deg²
removed area and **0** removed targets for bit2; the text explains
that bits 1–5 were used at the target-selection stage.
Source: https://repositorio.uam.es/bitstreams/9c5e261f-a786-4897-95a5-e84b3fc3d5ed/download
(pdf page 6, Sec. 3.2 and Table 5). The pattern is consistent with
our four fixed sample-image histograms: 2,080,939 bit2 pixels in
`eboss23`, 1,345 in `eboss21`, 3,542 in `eboss22`, and zero
in `eboss25`. This literature comparison was made **after**
the sample histograms were observed and is explicitly a
*post-result qualitative corroboration*, not a preregistered
independent test. Sample pixel counts are neither survey-area
totals nor paper target-veto counts. Source-code bit2 handling
depends on WCS rounded and truncated sky-coordinate locations;
histograms alone do not implement that handling. No change
to a catalog's acceptance cuts follows from this comparison.

The paper also says in Sec. 3.2 that published veto bits 1–7
are reproducible with `brickmask`, while bits 8–11 and the
two eboss22 low-quality plates require custom Python handling.
Our independently documented unresolved *bit8 RA
interpretation* still needs a release-authenticated production
reference containing both positive and negative bit8
positions before any observed pair-window work.
