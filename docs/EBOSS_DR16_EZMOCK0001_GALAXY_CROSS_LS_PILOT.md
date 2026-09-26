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
https://github.com/dvlahek/stress-energy-closure/actions/runs/36226935218 .
Run locally on the existing draft audit branch:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py --self-test && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_cross_ls_pilot.py`

Local result:
`eboss_workspace/official_mask_inventory/ezmock0001_galaxy_cross_ls_pilot.json`.
Upload the **unchanged** JSON, including any fail-closed
status, for a separate exact-byte archive before
further mock expansion. Preserve all fixed source
IDs/caps/bin edges/sampling seeds and the sealed observed
odd data vector.

The first CI attempt failed because its source-only manifest check incorrectly compared the **user's absolute WSL file path** to GitHub runner's different absolute checkout location. The archive, raw source SHA256 pins, and actual local binary gate were unaffected. The fixed CI instead validates the SHA-pinned report's **relative source suffix**, then the actual user-run stage must rehash the entire local binary file at its own checkout path. The updated synthetic/source-only CI [passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36226935218); no real FITS galaxy/random rows are accessed in CI.

## Completed fixed realistic mock0001 galaxy-pair algebra and exact report archive

The actual source-verified local four-`dat` plus four matched-`ran`
mock0001 run finished for both Galactic caps. The user's exact
`ezmock0001_galaxy_cross_ls_pilot.json` is **17,995 bytes**,
SHA256 `b8ea4545d433900c1364d1c12e2638a2d2288e35f0a327c22b8ad1ec9ae8510c`.
It is archived byte-identically in
`source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_report_2026-09-26.json`,
with independently verified Git blob
`4f35e17a44a226235d9f45ff3afb5adc5f0c3579`.
The independently registered report fingerprint, exact eight input
source SHA256s, accepted pair totals and no-unblinding assertions are
`source_data/eboss_dr16_ezmock0001_galaxy_cross_ls_pilot_uploaded_manifest_2026-09-26.json`.

All eight full source gzip SHAs matched their existing prior pins
*before* any FITS row reads. The previously locked high-z bin
`0.9 <= z < 1.0` contained eligible galaxy rows:
NGC LRG `7537`, NGC ELG `17149`; SGC LRG `4159`, SGC ELG
`15860`. The fixed sample uses `600D` and `1200R`
per tracer per cap, without numerical-zero exclusions or
chunk-conditioned sampling. In NGC, the four weighted cross
pair terms have accepted raw pair counts
`DD=1486, DR=2825, RD=3946, RR=7151`.
In SGC, those counts are
`DD=5637, DR=11159, RD=10414, RR=20102`.
All four independently normalized pair terms reproduce their
reverse-orientation signed-`mu` mirrors in **both** caps,
with forward/reverse `xi(s,mu)` residuals
`8.881784197001252e-16` (NGC) and
`1.3322676295501878e-15` (SGC), at 144/144
RR-supported cells in each cap. This is a true
**realistic-mock-galaxy estimator algebra closure**, unlike the
earlier nine-mock random-only pseudo-D/R test.
It is **not** a measurement that the mock odd
`xi_1`/`xi_3` vanish, not an observed odd-sector constraint
and not a mock covariance.

## Follow-up: fixed descriptive mock `xi_ell`, no observed odd analysis

Before any mock galaxy odd multipole was derived, a separate,
**post-first-pilot descriptive** protocol was recorded:
`source_data/eboss_dr16_ezmock0001_galaxy_odd_projection_protocol_2026-09-26.json`.
The new runner
`scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py`
uses the existing `0001` mock files, same original
`600D/1200R` sample indices/weights, both caps, same
fiducial, LOS, angular cut and 6x24 signed grid.
It first SHA-compares the local pilot report with its exact
Git archive, then rehashes **all eight complete gzip
files**, independently reproduces eight original selected
sample-array hashes and every original forward/reverse
weighted-pair histogram SHA256. It fails closed if any
source/sample/pair differs. Only then does it project the
reconstructed **mock-only** `xi(s,mu)` onto fixed
`ell=0,1,2,3` using exact bin-integrated Legendre
weights clipped to physical `-1<=mu<=1`.
All 144 random-supported cells must remain present in
both caps, with no NaN-to-zero substitution, and the
independently reversed mock projection must obey
even/odd parity.

The projection protocol was registered **after**
the successful 0001 estimator pilot and its support
counts were observed. It is a post-pilot implementation
diagnostic, **not** an independent preregistered
scientific null or hypothesis test. A mock's `ell=1`
and `ell=3` can be nonzero because of finite
sampling and tracer/window effects. Their measured
values must not be used to tune a science cut or
select mock IDs. For actual inference the next
stage requires an adequate many-realization
source-pinned mock galaxy + matched-random ensemble
and a separately qualified joint covariance;
the previous nine IDs are preliminary source/window
coverage, not a credible full 18D covariance.

The exact-archive + synthetic projection CI
[passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36232109867).
The first CI attempt failed because GitHub's
runner lacks the user's local JSON; CI now checks the
committed exact archive, while **the real WSL
runner still requires the local JSON to match the
archive byte-for-byte**, with all eight binary SHAs
checked before any rows. Run locally (no additional
files if prior mock random quarantine remains present):

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py --self-test && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py`

Upload the unchanged local
`eboss_workspace/official_mask_inventory/ezmock0001_galaxy_odd_projection.json`
even on a fail-closed run. At no point open real
observed galaxies, odd `xi`, wake matched filters or
derive detection significance from this single mock.

## First local 0001 mock-only multipole run stopped on a metadata-key regression

The first local follow-up run reported
`EZMOCK0001_MOCK_GALAXY_ODD_PROJECTION_INCOMPLETE_STOP`
with `KeyError('pair_normalization')`. All four complete
mock-galaxy and four complete same-realization mock-random compressed
source SHA256 checks had already passed. The failure occurred in
`verify_pair_reproduction`: the `oriented_pair_terms` API returns
`(histograms, norms, metadata)`. Metadata includes
`independently_normalized_pair_weight`, while the actual
`pair_normalization` is in the separate `norms` dict.
The mock-only projection never completed; **no P1/P3 numerical
results from the failed run exist**. It did not read observed galaxies
or the observed odd-sector vector.

The exact uploaded fail-closed 334-byte JSON report,
SHA256
`4423ab407f5f5da67db4736623e18d29c5d77177e07c1c21e6d5023743e962e7`,
was archived byte-identically to
`source_data/eboss_dr16_ezmock0001_galaxy_odd_projection_failure_2026-09-26.json`;
its Git blob is
`fedd41a9fac42de262b264375851917270699af9`.
It is retained as a failed-attempt provenance artifact, **not**
a substitute for a complete run.

The implementation was corrected to pass the actual independently
returned `fnorm[label]` / `rnorm[label]` into the
SHA-gated archived pair-reproduction check. It also validates
the separate metadata's independent weight normalization,
archived accepted-pair count, weighted histogram SHA256,
and recorded normalization residual. A synthetic regression
constructs the same public return shape with **no**
`pair_normalization` key in metadata; its positive control
must succeed, while tampered norms and histogram bytes
must fail. [Corrected source-only/synthetic CI passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36232522998).
No source hash, catalogue row count, sample seed,
fiducial, redshift/angle/separation/μ bin or science rule
was changed.

The real local rerun retains exact source and prior-pilot
report gates. Re-run:

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py --self-test && python -u scripts/audit_eboss_dr16_ezmock0001_galaxy_odd_projection.py`

Upload the complete or fail-closed unchanged
`eboss_workspace/official_mask_inventory/ezmock0001_galaxy_odd_projection.json`.
Only a successfully SHA-reproduced **mock-only**
report may be used for the already declared descriptive
single-mock P1/P3 projection; no statistical significance
or observed odd access follows.
