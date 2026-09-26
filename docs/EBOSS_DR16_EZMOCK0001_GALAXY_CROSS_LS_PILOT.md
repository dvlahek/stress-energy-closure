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
