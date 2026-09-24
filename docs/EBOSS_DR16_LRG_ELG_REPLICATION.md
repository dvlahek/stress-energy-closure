# eBOSS DR16 LRG–ELG odd cross-correlation

## Objective

We test the parity-odd wake response specified by the Einstein–Vlasov model in the eBOSS DR16 LRG–ELG cross-correlation. This is an external-survey replication of the DESI DR1 exact-pair test, not a re-optimization of its measured signal. The DESI data vector, templates, covariance, analysis cuts and reported significance remain unchanged.

eBOSS and DESI may cover common comoving volumes. We therefore report their results separately unless the cross-survey covariance can be estimated. The use of different spectroscopic selections does not by itself make their cosmic-variance fluctuations independent.

## Samples and external inputs

We adopt the published eBOSS DR16 *pre-reconstruction* clustering catalogues. The intended red-tracer sample is the combined LRGpCMASS catalogue used in the published eBOSS LRG–ELG cross-correlation analysis; the blue-tracer sample is the DR16 ELG clustering catalogue. We verify the exact released filenames, mask definitions, weight columns and matching multi-tracer mock selection from catalogue metadata before constructing an odd data vector. We do not substitute pure eBOSS LRGs, another catalogue version or different mocks on the basis of the measured odd signal.

We analyse NGC and SGC with their respective data and random catalogues. The eBOSS randoms, completeness corrections and FKP/systematic weights are *not* assumed to have the same column names or normalization as DESI DR1. The final data–random weighting convention is determined from the released eBOSS catalogues and their accompanying survey documentation before odd-sector measurements.

The common redshift interval is `0.60 <= z < 1.00`. The four analysis bins are `[0.60,0.70)`, `[0.70,0.80)`, `[0.80,0.90)`, and `[0.90,1.00)`. We use six separation bins with edges `20,40,60,80,100,120,140` in `Mpc/h`, 240 angular bins, the midpoint line of sight and LRGpCMASS-to-ELG pair orientation. The primary statistic is the 24-component `ell=1` vector. The `ell=3` vector provides a control. We adopt the 0.05-degree angular exclusion used by the DESI exact-pair method, subject to verifying that the eBOSS catalogue cuts and randoms implement it consistently. Any infeasible bin is identified from catalogue and mock support before an odd-sector data vector is evaluated; a failed feasibility check is reported, not repaired by choosing bins from the observed signal.

The fiducial redshift-to-distance relation is taken from the published eBOSS clustering convention and applied identically to data, randoms and mocks. Its numerical implementation and source are recorded with the catalogue manifest.

## Physical model

We keep the Einstein–Vlasov hidden-state construction fixed: `m_nu=0.06 eV`, `z_match=1100` and a 30 per cent pointwise deformation cap. The theoretical redshift dependence is evaluated at the eBOSS effective redshifts. Numerical DESI template arrays are not reused as eBOSS templates.

The standard odd component includes the linked relativistic and leading wide-angle benchmark with eBOSS-specific tracer inputs obtained independently of the measured odd vector. The even sector is propagated through the measured eBOSS random-pair window. We report the validity domain of the linear Kaiser approximation, the treatment of the selection function and any unmodelled wide-angle or nonlinear corrections. We do not interpret a benchmark-template excess as a model-independent detection of the hidden state.

## Covariance and calibration

The intended covariance sample is the public eBOSS DR16 multi-tracer EZmock ensemble, with LRG and ELG catalogues paired by realization and sharing the same underlying mock density field. The released ensemble contains 1000 realizations. We verify that the paired files reproduce the selected data sample, Galactic caps, random selection and systematic treatment before choosing the usable ensemble. A pilot subset may be used to validate file ingestion and pair-counting; it is not used for significance-based stopping.

Data and mocks use the same estimator and window convention. We examine mock mean, covariance eigenvalues, condition number, finite-sample precision correction, cap-specific behavior and estimator closure on non-data catalogues. We calibrate the prescribed one-dimensional nuisance-projected matched filter against held-out mocks and report empirical two-sided `p=(1+n_ge)/(1+N)`. We also evaluate the Sellentin–Heavens likelihood-ratio statistic. The complete usable ensemble and any exclusions are documented before calculating the observed odd-sector statistic. The target is all 1000 released paired realizations if their catalogue selection and availability can be confirmed.

The published multi-tracer mock comparisons identify limitations in the LRG–ELG small-scale cross power. We assess their relevance for our separation interval using mock-only and published clustering diagnostics before treating the ensemble as a reliable null calibration.

## Analysis order

1. Record released catalogue versions, headers, mask and weight fields, caps, effective redshift support and matched mock realization IDs. This step does not measure data odd multipoles.
2. Validate redshift selection, random normalization, orientation reversal and independent small-subset pair counting. Inspect random-pair geometry and calculate the window operator.
3. Compute and document the eBOSS-specific wake and standard templates from prescribed theory and external tracer inputs.
4. Validate the covariance and matched-filter distribution with paired mocks, including pseudo-signal injection and leave-one-out calibration. Record the final model inputs and covariance sample.
5. Measure the observed eBOSS odd vector once with the preceding configuration and report the result, including a null outcome, without changing bins or templates in response.

The primary DESI DR1 result is not refitted. A joint DESI–eBOSS significance is not quoted without an estimated cross-survey covariance.

## Data and methods references

- SDSS, [DR16 LSS catalogues and clustering-ready data/random definitions](https://www.sdss4.org/dr17/spectro/lss/).
- SDSS, [public DR16 LSS catalogue directory](https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/).
- Zhao et al., [DR16 multi-tracer EZmock catalogues](https://academic.oup.com/mnras/article/503/1/1149/6149160).
- Wang et al., [joint eBOSS DR16 LRG and ELG clustering in configuration space](https://academic.oup.com/mnras/article/498/3/3470/5897383).
- eBOSS multitracer full-shape study, [1000 EZmock comparison and its cross-spectrum limitations](https://academic.oup.com/mnras/article/532/1/783/7693747).
