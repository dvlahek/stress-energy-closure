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
