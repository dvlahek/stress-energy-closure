# Nine fixed realistic EZmock GALAXY source hashes: next blinded code-transport gate

**Status (2026-09-26):** one realistic EZmock realization `0001` has passed matched LRG×ELG cross-Landy–Szalay orientation closure and a separate descriptive mock-only `ell=0,1,2,3` projection. The exact successfully uploaded output is frozen. No observed odd-sector data or real-galaxy rows have been accessed in this stage.

## Completed, exactly archived `0001` mock multipoles

The user's completed local `ezmock0001_galaxy_odd_projection.json` is **8,736 bytes**, SHA256 `7900385533cd2d79768333ab4ca9f4c8f38ee65f26f063544541457d7c35cbe6`. Its exact bytes are archived in `source_data/eboss_dr16_ezmock0001_galaxy_odd_projection_report_2026-09-26.json`, Git blob `3a629edcf0194ee17838b23476306f83d44ac9cd`. The file's exact identity and scope are independently recorded in `source_data/eboss_dr16_ezmock0001_galaxy_odd_projection_uploaded_manifest_2026-09-26.json`. The previous `KeyError('pair_normalization')` failure remains archived separately, not overwritten.

All eight original 0001 gzip source SHA256s, original sample-array SHA256s and prior forward/reversed weighted pair histogram SHA256s were reproduced. Both caps have 144/144 RR-supported `(s,mu)` cells. All even and odd forward/reverse `ell=0,1,2,3` projected parity residuals are below `1e-15` in the exact report. The first `s=20–40 h^-1 Mpc` bin has mock NGC `P1=0.42460724478109463`, `P3=0.5888369397917503`, and mock SGC `P1=0.3431104851431689`, `P3=0.4109635954637549`. These are finite-sample, **one-mock** descriptive multipoles with only `600D/1200R` per tracer/cap. Reversal is an algebraic identity, not a physical odd null; the nonzero values are not evidence of a cosmological signal. They must not be used to change redshift/separation/`mu` bins, weights, angular vetoes, sampling seeds or hypothesis templates.

## Frozen nine-ID source-only extension

The nine mock realization IDs had already been selected **before the first 0001 mock galaxy odd projection**, in `source_data/eboss_dr16_predeclared_9mock_window_pilot_2026-09-24.json`: `0001,0125,0250,0375,0500,0625,0750,0875,1000`. A new input-only protocol registered **before downloading any of the other eight IDs' galaxy data** is `source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_protocol_2026-09-26.json`. It retains the four already frozen `0001` galaxy sources and fixes exactly **32 additional** official realistic EZmock `dat.fits.gz` files: the remaining eight IDs × NGC/SGC × eBOSS_LRG/eBOSS_ELG. Each newly acquired file has a **16 MiB** full compressed-byte cap and the additional 32 together have a **512 MiB** cap.

`scripts/audit_eboss_dr16_nine_ezmock_galaxy_raw_sha.py` first verifies the byte-identical archived successful 0001 odd report, four original 0001 source SHA pins, and the prior nine-ID/random-only ensemble registration. In actual WSL mode, it rehashes the four old local 0001 source files, then acquires only the 32 exact other official HTTPS gzip DATA URLs and records their **first-seen full compressed-byte SHA256s**. Each successful file is checkpointed. A previously checkpointed file is reused only after complete SHA and exact filename/size verification. Existing unpinned files, changed sources, changed realization IDs and alternate mirrors **fail closed**. Nothing is decompressed. Neither mock galaxy nor random FITS rows are read, no pair count is made, and observed odd data remain sealed.

[Source-only/synthetic CI passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36233490532). It checks the frozen 9-ID protocol, uploaded 0001 report/manifest identities and rejection of synthetic corrupted gzip transfers; it does **not** acquire the 32 real files in CI.

One WSL command (on the existing draft audit branch):

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_nine_ezmock_galaxy_raw_sha.py --self-test && python -u scripts/audit_eboss_dr16_nine_ezmock_galaxy_raw_sha.py`

The result is `eboss_workspace/official_mask_inventory/ezmock_nine_galaxy_source_raw_sha.json`. Upload the **unchanged JSON**, including any fail-closed checkpoint result. Its exact bytes and 32 new first-seen SHA256s must be archived in a separate immutable manifest *before* any of their mock galaxy rows are read.

## What a completed source inventory does not establish

The existing nine-mock **random-only** audit previously validated its own data sources and estimator algebra. It does not automatically authenticate the companion 32 new mock-random compressed files on the user's local runtime. Before pairing any of the new mock galaxies, rehash and freeze **all same-realization, tracer- and cap-matched mock randoms** under the earlier random source protocol. Only then can a separately preregistered nine-ID galaxy *code-transport* run use the same `[0.9,1.0)`, `600D/1200R`, 6×24 `(s,mu)`, angular cut and midpoint LOS without after-the-fact tuning. Missing RR support must be reported, not fabricated by zero-fill or discarding an unfavorable mock.

Nine independent mock realizations cannot yield an invertible, scientifically qualified **18-dimensional** sample covariance: the centered sample covariance has rank at most eight. This cohort is a preliminary source/implementation and fixed-response diagnostic only. A meaningful physical mock null and covariance require a separately designed, much larger source-pinned ensemble and a validated published LRG×ELG empirical random pair window. Neither this audit nor mock odd multipoles authorize reading the observed odd signal.

## Completed actual 36 mock-galaxy source inventory

The user completed the frozen local source-only nine-ID run.
The exact uploaded `ezmock_nine_galaxy_source_raw_sha.json`
is **37,105 bytes**, SHA256
`2c6b55bf5611dc71449ecbbfc2dedf545f2d22250147ed08f6a072536b50f06c`.
It was archived byte-identically as
`source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_report_2026-09-26.json`
(Git blob SHA1 `43be7aae5b06580867bee4fff08292ebeee468bf`).
Its separately frozen all-36-data-source input manifest is
`source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_uploaded_manifest_2026-09-26.json`.
Four original `0001` full-gzip source SHAs were reverified;
the 32 newly registered exact source filenames/HTTPS URLs,
gzip magic, Content-Length and full SHA256s were verified.
The extra 32 sources total **64,937,949 compressed bytes**
and the 36-galaxy-source inventory totals **73,245,570 bytes**.
No new mock-galaxy rows or random rows were read during
this source stage. The 32 first-seen source checksums are
local byte provenance, **not** independently released
official SDSS checksums.

## Exact original nine-ID random SHA reference recovered

The original successful fixed nine-mock GitHub Actions
[window ensemble, run 36017670812](https://github.com/dvlahek/stress-energy-closure/actions/runs/36017670812)
contains artifact `eboss-dr16-nine-mock-window-ensemble`,
artifact ID `10815013322`. Its archived ZIP SHA256 is
`a3b23502073e2a2e46211582d13ee4c1e45cd4984fdedf1ef36df4e0c5a8cda0`.
The exact `ninemock_window_ensemble.json` ZIP member
is 337,280 bytes, SHA256
`8702adaf1284068c26e79555feb76c7a709df9a265f9e2194409cc4161602e01`,
with parent revision commit
`054005edc6ca193176b129e1951e4bd3ce8751a7`.
Its 18 cap-by-realization cases contain **36 pre-existing**
full-compressed-file LRG/ELG mock-random SHA256s.
All 12 values from the earlier independently committed
`0001/0500/1000` mock-random inventory match.
The extracted, per-file original SHA references are pinned
in `source_data/eboss_dr16_nine_ezmock_random_reference_sha_from_20260924_artifact.json`.
These are **earlier source identities**, not first-seen SHAs
selected using the 2026-09-26 mock odd values.

## New separate random-only full-source verification

The new
`source_data/eboss_dr16_nine_ezmock_matched_random_sha_protocol_2026-09-26.json`
freezes the original 36 random-source keys, original
artifact identity, exact successful local 36-galaxy-source
uploaded report/manifest, approved local cache paths, official
same-realization/cap/tracer `ran.fits.gz` URLs,
128 MiB compressed cap **per random** and 4 GiB overall.
`scripts/audit_eboss_dr16_nine_ezmock_matched_random_sha.py`
checks all original source/manifests and the byte-identical
user's local 36-galaxy source report **before** touching
any new random binary. It then rehashes complete random
gzip bytes against the immutable **2026-09-24 reference**,
reusing known cache locations and checkpointing each success.
It does **not** decompress FITS, open mock/observed rows
or compute pairs/odd statistics. A prior-hash mismatch
stops the audit; no repinning or substitute ID is allowed.

By default, this runner **does not download missing
random files**. It produces an incomplete checkpoint
listing the first missing fixed source. Missing official
inputs can only be acquired using the separate opt-in
`--download-missing` flag; such downloads must match
the already frozen 2026-09-24 SHA values. These files can
be collectively several gigabytes, so the cached-only
pass avoids unexpected bandwidth or storage.

[Archive/manifest and synthetic CI passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36234497932);
that CI opened no real random FITS. A single WSL cached-only
check is:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_nine_ezmock_matched_random_sha.py --self-test && python -u scripts/audit_eboss_dr16_nine_ezmock_matched_random_sha.py`

Upload the unchanged result
`eboss_workspace/official_mask_inventory/ezmock_nine_matched_random_source_sha.json`
even on an incomplete run. A future completed nine-ID
galaxy estimator pilot must first freeze the exact successful
random-only upload. Do **not** infer a physical mask, observed
odd detection or 18D covariance from this source provenance;
nine mock realizations yield sample covariance rank at most 8.
