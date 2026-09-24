# DESI DR1 luminosity-rank octupole control

We measure the odd octupole on the same DESI DR1 BGS tracer sample and sampled pair geometry used for the luminosity-rank dipole. This provides a higher-multipole consistency test without fitting a second wake amplitude.

## Observable and null distribution

The dipole and octupole angular weights are $3\mu$ and $7P_3(\mu)=7(5\mu^3-3\mu)/2$, respectively. We construct an 18-component octupole vector with three redshift intervals and six separation bins. Its global test statistic is the Mahalanobis distance from the permutation-null mean, calibrated with 32 within-stratum luminosity-mark permutations.

We use 30 delete-one jackknife regions. For inversion, the jackknife covariance is regularized as

$$
C_{\mathrm{reg}}=C_{\mathrm{JK}}+rI,\qquad
r=\max\left[10^{-7}\lambda_{\max}(C_{\mathrm{JK}}),
10^{-6}\operatorname{median}\operatorname{diag}(C_{\mathrm{JK}}),10^{-14}\right].
$$

The reported condition number, $440.8145$, refers to $C_{\mathrm{reg}}$. The unregularized covariance has rank 18; its condition number was not reported in this run.

## Results

The global octupole Mahalanobis statistic is 15.6936. Its empirical permutation $p$-value is 0.81818; the asymptotic chi-square diagnostic gives $p=0.61392$. The largest single-bin fluctuation has $|z|=1.3221$. These results are consistent with the permutation null.

The same calculation reproduces the **corresponding 32-permutation dipole reference vector** exactly: the maximum absolute difference and RMS difference are zero, and the correlation is one. This is an implementation check, not an additional independent measurement. The reference vector is stored in `source_data/phase7_octupole_control/reference_dipole_32perm.csv`.

The primary reported DESI dipole coefficient uses a later 256-permutation null calibration, $A_{\rm wake}=-0.0768409\pm0.0851660$ with empirical two-sided $p=0.37354$. The octupole reproduction check is not a direct equality test against that later null-corrected vector.

## Source data and reproduction

The measured vector, covariance, compact numerical summary, reference dipole and provenance files are in `source_data/phase7_octupole_control/`. The local workflow is:

```bash
bash scripts/run_phase7_octupole_control_local.sh
```

The workflow uses the same source catalogues and recorded pair-sampling seed as the original octupole calculation. Large external DESI catalogues are not stored in this repository.
