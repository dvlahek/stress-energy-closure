# eBOSS DR16 LRG×ELG: catalogue inventory and prospective replication

Status: catalogue-link inventory only, 2026-09-24. No eBOSS FITS files have been downloaded or validated in this repository, and no eBOSS odd-sector data vector has been measured or inspected. The analysis definition below is a prospective plan, not a frozen protocol.

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

The SDSS [DR16 LSS documentation](https://www.sdss4.org/dr17/spectro/lss/) describes the released clustering data and random catalogues. SDSS also publishes a [DR16 multi-tracer EZmock collection](https://sdss.org/dr20/data_access/value-added-catalogs/?vac_id=68). Its listing does not, by itself, establish a jointly usable LRG×ELG mock ensemble or its cross-covariance.

## Validation gates before measuring the odd signal

1. **Catalogue integrity and provenance.** Download and record exact catalogue versions, byte sizes and SHA256 digests. Validate FITS structure, row counts, angular coordinates, redshift and quality fields, and the catalogue-specific data and random weights. Confirm the published LSS data model and the selection and failure-correction conventions separately for LRG and ELG. Do not assume that DESI DR1 weight or random-catalogue conventions apply to eBOSS.

2. **Common selection.** Determine the actual joint sky mask and usable redshift overlap from catalogue metadata, sky geometry, data/random distributions and quality cuts. Check NGC and SGC separately. The previously discussed approximate interval `0.6 < z < 1.0` is a candidate only; it is not a frozen redshift selection or binning. Record all inclusion criteria and exclusions before examining the odd-sector result.

3. **Mock eligibility.** Verify that any proposed LRG and ELG mock realizations form matched joint realizations with compatible selection, redshift distributions, angular geometry and random catalogues. Establish a fixed, sufficient mock ensemble and its cross-tracer covariance. Do not substitute separate single-tracer mock covariances for a validated cross-covariance. If appropriate joint mocks are unavailable, document the limitation and redesign the calibration before unblinding.

4. **Estimator and forward-model closure.** Implement exact oriented LRG→ELG cross-Landy–Szalay pair counting, using midpoint line of sight if the catalogues and selection support the same definition. Validate pair counting against an independent implementation; test angular resolution, random density, cap combination and the adopted small-angle exclusion. Construct the actual eBOSS random-pair survey window and its even-to-odd leakage. Use the same prescribed physical response as in DESI where its assumptions remain valid, with independently justified eBOSS tracer inputs. Do not use the real eBOSS odd signal to tune the template.

5. **Prospective statistical protocol.** Before opening the real odd-sector vector, commit the final catalogue hashes, quality cuts, common footprint, redshift and separation bins, multipoles, nuisance basis, wake template, mock membership and count, covariance estimator, finite-mock likelihood and empirical tail statistic. Specify the two-sided test, treatment of multiple analyses and stopping rule. Complete null, injection, mock and survey-window checks under this protocol.

6. **Real-data measurement and reporting.** Open the real odd-sector vector only after the preceding checks and the protocol commit. Report nominal and empirical finite-mock results, null controls and sensitivity to the declared modeling assumptions. Keep the DESI DR1 primary null result and secondary exact-pair consistency result unchanged. Do not combine surveys as statistically independent without an appropriate cross-survey covariance assessment.

This document records the design and outstanding checks. It does not certify that these gates have been completed.
