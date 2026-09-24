# eBOSS DR16 LRG×ELG: catalogue inventory and prospective replication

Status: pre-analysis catalogue and mock validation, 2026-09-24. We have parsed all eight public real-data and random FITS headers, inspected the four full data catalogues and all four full random catalogues for selection metadata, audited 1000 matched EZmock filename sets, and checked 24 bounded realistic-mock FITS headers. The exact common mask, full mock selection and cross-covariance remain unvalidated. No eBOSS odd-sector vector has been calculated or inspected, and the analysis protocol is not yet frozen.

## Scientific question

We use a second survey to test the prescribed odd-sector wake response without changing the DESI DR1 results. The eBOSS estimator, selection function, covariance and survey window must be validated independently before a real-data measurement. A similar measured amplitude alone would not establish an independent detection, particularly if the surveys share cosmological volume.

## Published catalogue links

The read-only inventory workflow enumerated the public [SDSS DR16 LSS catalogue directory](https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/) on 2026-09-24. It identified one clustering data file and one clustering random file for each tracer and Galactic cap:

| Tracer | Cap | Data FITS | Random FITS |
| --- | --- | --- | --- |
| LRG | NGC | `eBOSS_LRG_clustering_data-NGC-vDR16.fits` | `eBOSS_LRG_clustering_random-NGC-vDR16.fits` |
| LRG | SGC | `eBOSS_LRG_clustering_data-SGC-vDR16.fits` | `eBOSS_LRG_clustering_random-SGC-vDR16.fits` |
| ELG | NGC | `eBOSS_ELG_clustering_data-NGC-vDR16.fits` | `eBOSS_ELG_clustering_random-NGC-vDR16.fits` |
| ELG | SGC | `eBOSS_ELG_clustering_data-SGC-vDR16.fits` | `eBOSS_ELG_clustering_random-SGC-vDR16.fits` |

Each file has the URL formed by appending its published filename to the SDSS directory above. The [inventory workflow run](https://github.com/dvlahek/stress-energy-closure/actions/runs/35965943010) retains the complete read-only listing as artifact `eboss-dr16-release-index`. This inventory establishes file locations only. A successful filename inventory does not establish FITS validity, redshift or footprint overlap, compatible random selection, weight conventions or mock calibration.

The SDSS [DR16 LSS documentation](https://www.sdss4.org/dr17/spectro/lss/) describes the released clustering data and random catalogues. The catalogue-header and mock-filename checks below establish the published file structure. They do not establish the physical overlap, survey selection or a calibrated cross-tracer covariance.

## Header and mock-filename validation

We inspect FITS metadata through bounded HTTP range requests without reading catalogue rows. All eight public clustering files have a binary-table extension with `RA`, `DEC`, `Z`, `WEIGHT_SYSTOT`, `WEIGHT_CP`, `WEIGHT_NOZ` and `WEIGHT_FKP`. The declared row counts are:

| Tracer | Cap | Data header rows | Random header rows |
| --- | --- | ---: | ---: |
| LRG | NGC | 107,500 | 5,460,719 |
| LRG | SGC | 67,316 | 3,453,453 |
| ELG | NGC | 83,769 | 3,728,363 |
| ELG | SGC | 89,967 | 3,609,460 |

These are FITS `NAXIS2` values, not galaxy counts after quality, redshift or common-footprint cuts. We have not checked the contents or determined the formula for combining eBOSS weights. The [header inventory workflow](https://github.com/dvlahek/stress-energy-closure/actions/runs/35966541598) and `source_data/eboss_dr16_fits_header_audit_2026-09-24.json` record the metadata.

### Data-only selection checks

We downloaded the four public clustering data FITS files into a temporary CI workspace and calculated a full SHA256 digest for each. The selection audit verifies FITS structure and the finiteness and positivity of the four published weight columns. For the *candidate* interval `0.6 <= z < 1.0`, before common-footprint selection, the files contain 107,500 and 67,316 LRG objects in NGC and SGC, and 76,395 and 82,596 ELG objects. The 0.6–1.0 candidate contains all objects in these two eBOSS LRG clustering data files, but excludes the ELG objects above z=1.0. These counts do not define final analysis bins or a statistical test.

A coarse data-object sky-occupancy diagnostic finds 37 jointly occupied cells out of 66 ELG-occupied cells in NGC, and 70 out of 70 in SGC, on the diagnostic grid. This is **not** a measurement of the common survey mask: sparse samples, cell boundaries and the absence of a joint random-mask selection can affect it. The full-random diagnostic below now checks the gridded intersection by cap. An exact common selection and survey window still need to be established before choosing the final sample.

The [data-only selection workflow](https://github.com/dvlahek/stress-energy-closure/actions/runs/35967074244) and `source_data/eboss_dr16_data_selection_audit_2026-09-24.json` retain the full-data SHA256 digests, raw redshift counts and coarse diagnostic. The workflow did not calculate pair counts, multipoles, matched-filter scores or any wake amplitude.

We also enumerate the public [eBOSS DR16 EZmock v1_0_0 release](https://data.sdss.org/sas/dr17/eboss/lss/EZmocks/v1_0_0/). Both `complete` and `realistic` branches provide eBOSS LRG and ELG data files for IDs `0001`–`1000` in NGC and SGC. The complete mocks have per-realization shuffled randoms and separate shared random catalogues; the realistic mocks have per-realization randoms. The filename audit finds all 1000 joint LRG–ELG IDs in both caps, with no missing data or matched-random filenames. The [pairing-audit workflow](https://github.com/dvlahek/stress-energy-closure/actions/runs/35966739277) and `source_data/eboss_dr16_mock_filename_audit_2026-09-24.json` retain the counts and the listing fingerprint.

The published mock-construction study ([Zhao et al., 2021](https://academic.oup.com/mnras/article/503/1/1149/6149160)) states that the tracers in a given realization share initial conditions and the underlying density field. A published [LRG–ELG cross-correlation analysis](https://academic.oup.com/mnras/article/511/4/5492/6527584) uses the joint EZmock ensemble. A later [multi-tracer EFT analysis](https://academic.oup.com/mnras/article/532/1/783/7693747) notes mock–data cross-power discrepancies at approximately `k > 0.15 h/Mpc`, because the tracer populations are sampled independently from their shared underlying density field. We have verified matching released filenames and a bounded sample of FITS headers, not the full mock random selection or cross-covariance. We therefore retain a mock-adequacy check on the final configuration-space scales rather than treating the existence of 1000 files as sufficient calibration.

The release also contains `eBOSS_LRGpCMASS` products, which are distinct from the `eBOSS_LRG` sample inventoried above. We have not selected between those tracer definitions; the choice and the corresponding data/mock provenance must be recorded before unblinding.

## Full random-catalogue selection diagnostic

The [full-random selection workflow](https://github.com/dvlahek/stress-energy-closure/actions/runs/35967630362) downloaded and verified all four public random FITS files. The retained numerical summary is `source_data/eboss_dr16_joint_random_selection_audit_2026-09-24.json`; its workflow artifact contains the full selection-only output. The declared row counts and SHA256 digests were recorded for all four catalogues. All random rows have valid coordinates and redshifts, and the four inspected weight columns have finite positive values.

For the candidate `0.6 <= z < 1.0` interval, we formed a *diagnostic* random support on a grid with `0.5 deg` in RA and `0.01` in `sin(dec)`. Requiring at least 25 LRG randoms and 25 ELG randoms in the same cell gives 965 joint cells in NGC (of 2,013 supported ELG cells) and 2,228 in SGC (of 2,238 supported ELG cells). The corresponding joint-cell counts in the four `0.6-0.7`, `0.7-0.8`, `0.8-0.9`, `0.9-1.0` intervals are `938, 944, 922, 818` in NGC and `2045, 2206, 2204, 1812` in SGC. These are thresholded, finite-grid random supports, **not** an exact common survey mask, effective area or RR-window matrix. We retain the distinction between NGC and SGC; pixel threshold, angular boundary and selection-function robustness still require validation.

A separate [bounded realistic-mock header audit](https://github.com/dvlahek/stress-energy-closure/actions/runs/35967950354) inspected the LRG and ELG data and random FITS headers of realizations `0001`, `0500` and `1000`, in both caps. All 24 headers parsed, and each tracer/role schema is consistent across the inspected IDs and caps. The realistic mock files provide `RA`, `DEC`, `Z` and the four weight columns, but the header audit does not validate the full mock selection. The compact metadata summary is `source_data/eboss_dr16_realistic_mock_header_audit_2026-09-24.json`. For example, LRG NGC realization `0001` declares 129,262 data rows, compared with 107,500 in the public real LRG NGC catalogue. The discrepancy needs a redshift/selection comparison before mock covariance is accepted.

## Validation gates before measuring the odd signal

1. **Catalogue integrity and provenance.** Download and record exact catalogue versions, byte sizes and SHA256 digests. Validate FITS structure, row counts, angular coordinates, redshift and quality fields, and the catalogue-specific data and random weights. Confirm the published LSS data model and the selection and failure-correction conventions separately for LRG and ELG. Do not assume that DESI DR1 weight or random-catalogue conventions apply to eBOSS.

2. **Common selection.** Determine the actual joint sky mask and usable redshift overlap from catalogue metadata, sky geometry, data/random distributions and quality cuts. Check NGC and SGC separately. The previously discussed approximate interval `0.6 < z < 1.0` is a candidate only; it is not a frozen redshift selection or binning. Record all inclusion criteria and exclusions before examining the odd-sector result.

3. **Mock eligibility.** Verify that any proposed LRG and ELG mock realizations form matched joint realizations with compatible selection, redshift distributions, angular geometry and random catalogues. Establish a fixed, sufficient mock ensemble and its cross-tracer covariance. Do not substitute separate single-tracer mock covariances for a validated cross-covariance. If appropriate joint mocks are unavailable, document the limitation and redesign the calibration before unblinding.

4. **Estimator and forward-model closure.** Implement exact oriented LRG→ELG cross-Landy–Szalay pair counting, using midpoint line of sight if the catalogues and selection support the same definition. Validate pair counting against an independent implementation; test angular resolution, random density, cap combination and the adopted small-angle exclusion. Construct the actual eBOSS random-pair survey window and its even-to-odd leakage. Use the same prescribed physical response as in DESI where its assumptions remain valid, with independently justified eBOSS tracer inputs. Do not use the real eBOSS odd signal to tune the template.

5. **Prospective statistical protocol.** Before opening the real odd-sector vector, commit the final catalogue hashes, quality cuts, common footprint, redshift and separation bins, multipoles, nuisance basis, wake template, mock membership and count, covariance estimator, finite-mock likelihood and empirical tail statistic. Specify the two-sided test, treatment of multiple analyses and stopping rule. Complete null, injection, mock and survey-window checks under this protocol.

6. **Real-data measurement and reporting.** Open the real odd-sector vector only after the preceding checks and the protocol commit. Report nominal and empirical finite-mock results, null controls and sensitivity to the declared modeling assumptions. Keep the DESI DR1 primary null result and secondary exact-pair consistency result unchanged. Do not combine surveys as statistically independent without an appropriate cross-survey covariance assessment.

This document records the design and outstanding checks. It does not certify that these gates have been completed.
