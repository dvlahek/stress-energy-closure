# Development history

## 2026-09-15 — Observable-retention deep dive

The observable hierarchy was extended beyond the CMB/RSD/wake comparison to test how source-matched hidden kinetic directions survive in other response kernels.

### Screening branches

- Weak-lensing wake screening retained a large fractional hidden-state contrast (about 0.46 in the tested reference configurations), but no absolute weak-lensing survey detectability claim was made because line-of-sight projection is expected to strongly suppress the usable signal.
- kSZ-tagged wake screening retained about 0.50 of the hidden-state contrast and improved the wake weighting only modestly, by roughly 8–13% in the tested configurations.
- CLASS relative-velocity alignment showed strong directional alignment between the baryon/kSZ velocity and the neutrino-CDM relative velocity, with `|r|` about 0.87–0.93 across the tested redshift/smoothing grid. This is an alignment result, not a kSZ detection forecast.
- A parity-odd linear relative-velocity screening calculation found a much smaller hidden-state response, with mode-count RMS fractions of order `1e-4`–`1e-3`. This motivated a direct same-pair comparison against the resonant wake.

### Deep operator diagnostic

A seven-dimensional matched-moment response operator was constructed from small symmetric CLASS probes. The linearized optimization predicted a large `~3.17e-2` density-weighted proxy response for one direction, but direct nonlinear CLASS validation gave only `~6e-4`. The prediction therefore overfit the tiny small-probe response and is not used quantitatively.

The direct endpoint response initially appeared stable under changes in redshift finite-difference step, k-grid density and denominator mask. A later precision-convergence test showed that the response-optimized direction is not robust to tightened CLASS integration tolerances. It is therefore excluded from quantitative results.

### Direct CLASS velocity-transfer cross-check

The OOM-safe worker was extended to request `mPk,dTk,vTk`, exposing `t_cdm` and `t_ncdm[0]` directly. This allowed a direct relative velocity-divergence control, `theta_ncdm-theta_cdm`, without reconstructing velocity from redshift derivatives of density transfers.

For the reference/CREF deformation at `m_nu = 0.06 eV`, `z = 0.3`, a 30% pointwise cap and `0.003 <= k <= 0.2 h/Mpc`, the `nk=48` direct-theta run gave

- density-derived `v_(nu-cdm) P_cb` half-pair RMS: `5.369172869387564e-4`,
- direct `theta_(nu-cdm)` half-pair RMS: `3.9505423084350484e-4`,
- direct `theta_(nu-cdm) P_cb` half-pair RMS: `5.369438898120348e-4`,
- resonant wake half-pair response: `0.23084420034493347`.

The density-derived and direct-theta-times-`P_cb` proxies agree at the `~5e-5` relative level in this run, showing that the earlier density-derived construction was not a finite-difference artifact.

### k-resolved decomposition

The direct-theta diagnostic was then resolved over 96 logarithmic `k` modes. The exact identity

`delta(theta P) = bar(theta) delta(P) + bar(P) delta(theta)`

closed to absolute residuals of order `1e-16`. Pure `theta` and `theta P_cb` are therefore distinct transfer-level quantities, not two normalizations of one quantity.

For CREF at standard precision (`nk=96`):

- direct `theta P_cb` RMS: `5.281598155627618e-4`,
- density-derived `v P_cb` RMS: `5.280865049182505e-4`,
- wake / direct-`theta P_cb`: `437.0726313189232`,
- three logarithmic `k`-third RMS values: `4.409099614642714e-3`, `1.925236732732645e-3`, `4.7520094922303546e-4`.

The response is scale dependent, so no claim of a uniform pointwise `~400x` contrast is made across the full `k` interval.

### Moderate-precision convergence

The original aggressive high-precision preset repeatedly exhausted the available WSL memory. A lower-memory precision test was therefore defined by tightening only integration tolerances while leaving the ncdm hierarchy size and momentum grid unchanged.

For CREF, the moderate-precision `nk=96` run gave

- direct `theta P_cb` RMS: `5.322814846689793e-4`,
- density-derived `v P_cb` RMS: `5.322828757582342e-4`,
- wake / direct-`theta P_cb`: `433.68820256540243`.

The main `theta P_cb` RMS changed by only about `+0.78%` from standard to moderate precision, and the wake/linear contrast changed by about `-0.77%`. The pure `theta` fractional RMS is more precision sensitive and is not used for the main quantitative comparison.

For the response-selected direction, the moderate-precision result changed qualitatively:

- standard direct `theta P_cb` RMS: `6.464617200898742e-4`,
- moderate direct `theta P_cb` RMS: `1.0669648447319053e-5`.

The middle and high-`k` responses collapse under tighter tolerances and the number of zero crossings increases. The response-selected direction is therefore excluded from quantitative results. Its previously quoted wake/linear factors (`~285`–`288`) are retained only as development history.

### Response-independent orthogonal control

A second null direction was then defined without using any CLASS response amplitude: the first canonical coefficient-space basis axis was Gram-Schmidt orthogonalized against normalized CREF, normalized under the same pointwise cap, and propagated through the same pipeline. Its coefficient-space dot product with CREF is `-1.91e-17`, and its matched-moment mismatch is `2.39e-16`.

At standard precision (`nk=96`):

- direct `theta P_cb` RMS: `1.1713198944007974e-3`,
- density-derived `v P_cb` RMS: `1.1713282196477968e-3`,
- wake response: `0.1101394327438084`,
- wake / direct `theta P_cb`: `94.03019044609631`,
- zero crossings: `0`.

At moderate precision:

- direct `theta P_cb` RMS: `1.3876806621063746e-3`,
- density-derived `v P_cb` RMS: `1.3876847572915057e-3`,
- wake / direct `theta P_cb`: `79.36943689668965`,
- zero crossings: `0`.

The standard-to-moderate change in the integrated `theta P_cb` RMS is `+18.47%`, so this control is not precision-converged at the percent level. The change is driven mainly by the high-`k` third, whose RMS increases by about `26.1%`; the low and middle thirds change by about `0.44%` and `2.34%`. The qualitative separation remains strong: the integrated wake/linear contrast stays between about `79` and `94`, and the direct-theta-times-`P_cb` and density-derived proxies agree at the `1e-5` relative level or better in both precision settings.

This control is therefore retained as a response-independent qualitative robustness check, not as a second precision-grade coefficient.

### Current retained conclusion

The precision-grade quantitative control remains CREF: its integrated density-weighted direct velocity-divergence response is `5.28e-4` at standard precision and `5.32e-4` at moderate precision, compared with a resonant wake half-pair response of `0.230844`. The corresponding integrated wake-to-linear contrast remains about `4.3e2`.

The orthogonal response-independent control gives a separate qualitative confirmation that the resonant kernel can retain much more of the hidden kinetic perturbation than the integrated linear velocity kernel, with an integrated contrast of about `79`–`94` across the two precision settings. Because that second control shifts by about `18%`, it is not used as an exact second coefficient.

This is a transfer-level response diagnostic, not a kSZ or RSD survey forecast and not a theorem over all source-matched distributions.

### Resource note

The aggressive high-precision CLASS profile is not used for this diagnostic because it repeatedly exceeded the available development-machine memory. The moderate convergence test is intentionally narrower: it changes integration tolerances without simultaneously changing hierarchy size or momentum-grid resolution.

See `docs/HIDDEN_CHANNEL_RETENTION.md` and `source_data/hidden_channel_retention_direct.json` for the current reproducible status.

## 2026-09-21 — DESI LRGxELG regional robustness and estimator convergence

The z-resolved LRGxELG odd-sector validation was extended before the frozen 200-mock full-sky production run.

### Cap-specific NGC/SGC robustness

The regional split used separate NGC and SGC covariance matrices and cap-specific linked-standard nuisance templates. The external magnification/Doppler inputs were kept distinct between caps, and the DR1 evolution-bias benchmark was estimated separately from each cap's weighted selection function.

For the 18-component dipole with 40 r4 mocks per cap:

- NGC: wake amplitude `0.00670552 +/- 0.00184800`, nominal `Z = 3.62853`, `Delta chi2 = 13.1662`.
- SGC: wake amplitude `-0.000361950 +/- 0.00346529`, nominal `Z = -0.10445`, `Delta chi2 = 0.01091`.
- The approximate independent-cap amplitude difference is `1.80 sigma`, so the split is treated as inconclusive rather than a clean agreement or disagreement test.
- After nuisance projection, the retained wake metric norm is `0.76175` in NGC and `0.42275` in SGC. The SGC split is therefore substantially less informative for this fixed matched filter.

The 40-mock leave-one-out calibration places the NGC data near the mock tail but cannot resolve the nominal Gaussian significance: one mock is more extreme than the NGC data, giving the +1 empirical value `p = 2/41 = 0.04878` (`1.9705 sigma` two-sided). SGC is fully consistent with the mock ensemble.

### Extreme NGC mock audit

The only NGC mock more extreme than the data, mock `0008` with LOO matched-filter `Z = -4.55206`, was audited rather than removed.

An independent r4 rerun reproduced the stored dipole to

- maximum absolute difference `1.04e-16`,
- RMS difference `4.81e-17`.

Catalog and random counts were non-pathological, and the reverse-orientation closure remained numerically small. No objective catalog, random, estimator, or reproducibility failure was found. Mock `0008` is therefore retained as a valid tail realization.

### Angular-discretization convergence

The frozen production estimator uses 240 mu bins. Full-sky r4 reruns at 120 and 480 mu bins were compared without changing redshift bins, separation bins, random density, angular cut, cosmology, weights, or tracer definitions.

For the 18D dipole:

- `120 - 240`: median component shift `0.00445` mock sigma, maximum `0.01058`, covariance-metric shift `sqrt(Delta chi2) = 0.03545`.
- `480 - 240`: median component shift `0.00222` mock sigma, maximum `0.02971`, covariance-metric shift `sqrt(Delta chi2) = 0.02873`.
- `480 - 120`: covariance-metric shift `sqrt(Delta chi2) = 0.03725`.

The octupole is slightly more sensitive to angular discretization but remains negligible relative to mock scatter; the largest 120/240/480 covariance-metric shift is `0.11868`.

The mu-bin robustness test therefore passes. The production choice remains frozen at 240 bins; no significance-based retuning is performed.

### Production status

A new 200-realization full-sky r4 EZmock campaign has been started. The inspection order remains frozen:

`covariance diagnostics -> zero null -> nuisance-only -> wake`.

Intermediate wake significances are not inspected during production.

See `source_data/lrg_elg_ngc_sgc_regional_checkpoint_2026-09-21.json` and `source_data/lrg_elg_mubin_convergence_r4.json`.
