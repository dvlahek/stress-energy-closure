# eBOSS DR16: matched EZmock GALAXY cross-Landy–Szalay pilot (blinded)

**Status, 2026-09-26:** the next stage is frozen at the source-byte gate, not yet run on four local mock data files. All previously completed synthetic and random-only algebra tests remain valid. Do not compute observed galaxy odd multipoles or change the published tracer-specific angular selection.

## Why this is the next independent gate

The existing nine-realization random-only pseudo-D/R pilot established that the oriented four-term cross-estimator, independently weighted normalization and LRG→ELG versus ELG→LRG signed-μ reflection close on the fixed random-only samples. It did **not** involve galaxies drawn from the realistic EZmock density fields, so it cannot establish mock-galaxy estimator behaviour or mock-galaxy covariance. This next pilot pairs **realistic EZmock 0001 mock GALAXIES** with the corresponding **same-realization, same-tracer, same-cap released realistic EZmock RANDOMS**. It never pairs realistic mock galaxies with observed real-data galaxies or a different realization's randoms. The fixed nine-realization set remains `0001,0125,0250,0375,0500,0625,0750,0875,1000`; choosing 0001 first reprises an earlier fixed pilot and is not an odd-data or mock-outcome-driven choice.

The published LRG and ELG mask and redshift selections differ, so there is **no forced common sky mask**. Keep NGC and SGC independent. ELG exact `chunk` labels are retained and reported; source eligibility is not newly altered using the four example BRICKMASK pixel histograms or the still-unproven full historical production chain.

## Stage A: pin mock GALAXY source bytes before any row access

`source_data/eboss_dr16_ezmock0001_galaxy_raw_bytes_protocol_2026-09-26.json` fixes four exact realistic EZmock `v1_0_0` released gzip DATA files before this pilot:

| Cap | Tracer | Exact released gzip filename | Independently previously checked FITS-header rows |
|---|---|---|---:|
| NGC | eBOSS_LRG | EZmock_realistic_eBOSS_LRG_NGC_v7_0001.dat.fits.gz | 129262 |
| NGC | eBOSS_ELG | EZmock_realistic_eBOSS_ELG_NGC_v7_0001.dat.fits.gz | 83747 |
| SGC | eBOSS_LRG | EZmock_realistic_eBOSS_LRG_SGC_v7_0001.dat.fits.gz | 82921 |
| SGC | eBOSS_ELG | EZmock_realistic_eBOSS_ELG_SGC_v7_0001.dat.fits.gz | 90650 |

Approved base: `https://data.sdss.org/sas/dr17/eboss/lss/EZmocks/v1_0_0/realistic/` with the expected `eBOSS_LRG/dat/` or `eBOSS_ELG/dat/` subpath. The source-only runner `scripts/audit_eboss_dr16_ezmock0001_galaxy_raw_sha.py` checks the earlier nine-mock pilot, matching 0001 random SHAs and prior three-ID mock data/header inventories before a single input transfer. It downloads at most 16 MiB **per exact gzip** and 64 MiB combined, records first-seen whole-compressed-file SHA256 and source validators, and checkpoints the four files. A completed, locally byte-pinned input is only reused if the checkpoint and exact full SHA match. No gzip, FITS header, galaxy row, random row, pairing or observed odd vector is opened in Stage A. The new first-seen galaxy file SHA256 values are **not** previously published official checksums or invented from the earlier compact selection report.

[Stage A synthetic and earlier-manifest CI](https://github.com/dvlahek/stress-energy-closure/actions/runs/36225171171) passed without downloading any real file. One WSL command:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_raw_sha.py --self-test && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_raw_sha.py`

Upload the **unchanged** resulting `eboss_workspace/official_mask_inventory/ezmock0001_galaxy_data_raw_sha.json`; archive exact bytes and freeze four new source SHA256s in a separate stage-B manifest before mock FITS row decoding. Do not substitute another mock ID or source mirror after a checksum mismatch.

## Stage B already specified, not yet executable

The separate prior-to-galaxy-rows protocol is `source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_protocol_2026-09-26.json`. It predefines one high-z pilot bin `[0.9,1.0)`, 600 eligible mock DATA and 1200 **matched** mock RANDOM rows per tracer and cap (uniform no-replacement source subsamples, fixed seed root 93127), midpoint line of sight, published multitracer fiducial, `theta_min=0.05 deg`, `s=20,40,...,140 h^-1 Mpc` coarse edges and 24 signed μ cells. Use the already declared `WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP` rule and numerical-zero threshold `1e-12`. If eligibility is insufficient or source identity fails, stop without sample-size, mock-ID or redshift-bin substitution.

Compute all four independently normalized cross pair histograms, independently reverse tracer order, and verify full weighted-pair reflection and cross-LS reflection only on positive-RR-supported cells. Zero RR cells are missing, **not zeros in an observational odd multipole**. Preserve exact ELG chunk labels without forcing equal composition. Raw pair odd-parity reversal is an algorithmic check. A finite `600D/1200R` sample of one mock is **not** a physical mock covariance, an odd-null significance test, an observational detection, an all-19,381-image mask certification or an authorization to unblind. The eventual inferential 18D covariance must use separately fixed many-realization mock galaxies under matched selections and explicitly check mock adequacy.

The open historical ELG production-source and exact `eboss22` bad-plate veto geometry remain separate provenance gates. The already released tracer-specific random pairs sample the empirical cross window for this estimator test, without pretending they certify a continuous physical intersection.

## Stage A completed locally; exact upload archived and Stage B now available

The user completed the fixed four-mock-galaxy **source-byte-only**
gate for `0001`. The raw output
`eboss_workspace/official_mask_inventory/ezmock0001_galaxy_data_raw_sha.json`
has exact size **5,609 bytes** and independently verified SHA256
`ae9dc13ca27d7c4f4a6123cbc2526b44866aacc7db3add3f7a2b000773d86b4b`.
It is archived **byte-identically** in
`source_data/eboss_dr16_ezmock0001_galaxy_raw_sha_report_2026-09-26.json`,
confirmed by Git blob
`406617081d3f50ad872de41b53a2fe4e91acd4b2`.
Its separate post-source/pre-row frozen input manifest is
`source_data/eboss_dr16_ezmock0001_galaxy_raw_sha_uploaded_manifest_2026-09-26.json`.

| Same-ID fixed source | Full released compressed data bytes | Full gzip SHA256 first seen in local source audit |
|---|---:|---|
| NGC eBOSS_LRG `0001` | 2576333 | 67862ecb7860eef939f03c459ba424709cd37e8bbf00aa916f2bb2d351896536 |
| NGC eBOSS_ELG `0001` | 1933801 | aa399797762153b0ccff4960928c55336b35dc1131a44fbc735a7e668fbeb31e |
| SGC eBOSS_LRG `0001` | 1675016 | ba271958cf01a10c293b7c2876a5cad688c4618228c82378273cdb5451ae8a46 |
| SGC eBOSS_ELG `0001` | 2122471 | a7df8700eed4e4dbcc97cd184007efff2fcbe7e354001ed1334b1bf572d9595a |

Combined compressed bytes **8,307,621**. The report also
retains the four independently earlier SHA-pinned same-ID
realistic mock random file checksums; all source HTTP transfers
completed and no mock galaxy or random FITS rows, observed galaxy
rows or observed odd vector were read during Stage A. These
mock DATA first-seen SHA256s are local full-source
fingerprints, not independently published SDSS checksums.

A new `scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py`
first checks that the user's local source-only JSON is **byte-for-byte
identical** to the archived uploaded report, all four frozen
mock galaxy gzip bytes retain their full SHA256s, and all
four exact corresponding mock random files are similarly
verified against the previous official released-mock-random
audit. If a released mock random is not present in the
known local caches, it can fetch only the **same previously
SHA-pinned official URL**. It never acquires observed galaxy
files, compares to the observed odd vector, applies a new
sky mask or substitutes another mock ID.

Only after **all eight** source SHA checks pass does the runner
decode mock galaxy and mock random rows. It samples exactly
600 eligible mock DATA and 1,200 matched mock RANDOM rows
per tracer and cap in the frozen high-z interval,
retaining ELG exact chunk labels as diagnostics.
For each cap it calculates all four oriented
independently normalized cross terms and independently
recomputes their reverse-orientation reflection, RR
support and odd raw-RR parity. Unsupported RR cells
remain missing, not zero-imputed physical multipoles.
A one-realization finite sampled cross-LS mirror test
is NOT a galaxy covariance or physical null detection.

Synthetic-only CI, which audits the archived SHA manifest
and independent scalar pair-count closure but opens NO
real FITS sources:
https://github.com/dvlahek/stress-energy-closure/actions/runs/36226848080 .
Run locally on the existing draft audit branch:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py --self-test && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py`

Local result:
`eboss_workspace/official_mask_inventory/ezmock0001_galaxy_cross_ls_pilot.json`.
Upload the **unchanged** JSON, including any fail-closed
status, for a separate exact-byte archive before
further mock expansion. Preserve all fixed source
IDs/caps/bin edges/sampling seeds and the sealed observed
odd data vector.
