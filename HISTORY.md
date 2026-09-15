# Development history

## 2026-09-15 — Observable-retention deep dive

The observable hierarchy was extended beyond the published CMB/RSD/wake comparison to test how the same source-matched hidden kinetic directions survive in other response kernels.

### Screening branches

- Weak-lensing wake screening retained a large fractional hidden-state contrast (about 0.46 in the tested reference configurations), but no absolute weak-lensing survey detectability claim was made because line-of-sight projection is expected to strongly suppress the usable signal.
- kSZ-tagged wake screening retained about 0.50 of the hidden-state contrast and improved the wake weighting only modestly, by roughly 8–13% in the tested configurations.
- CLASS relative-velocity alignment showed strong directional alignment between the baryon/kSZ velocity and the neutrino-CDM relative velocity, with `|r|` about 0.87–0.93 across the tested redshift/smoothing grid. This is an alignment result, not a kSZ detection forecast.
- A parity-odd linear relative-velocity screening calculation found a much smaller hidden-state response, with mode-count RMS fractions of order `1e-4`–`1e-3`. This motivated a direct same-pair comparison against the resonant wake.

### Deep operator diagnostic

A seven-dimensional matched-moment response operator was constructed from small symmetric CLASS probes. The linearized optimization predicted a large `~3.17e-2` density-weighted proxy response for one direction, but direct nonlinear CLASS validation gave only `~6e-4`. The prediction therefore overfit the tiny small-probe response and is not used quantitatively.

The direct endpoint response initially appeared stable under changes in redshift finite-difference step, k-grid density and denominator mask. A later precision-convergence test showed that the response-optimized direction is not robust to tightened CLASS integration tolerances. It is therefore excluded from manuscript-level quantitative claims.

### Direct CLASS velocity-transfer cross-check

The OOM-safe worker was extended to request `mPk,dTk,vTk`, exposing `t_cdm` and `t_ncdm[0]` directly. This allowed a direct relative velocity-divergence control, `theta_ncdm-theta_cdm`, without reconstructing velocity from redshift derivatives of density transfers.

For the reference/CREF deformation at `m_nu = 0.06 eV`, `z = 0.3`, a 30% pointwise cap and `0.003 <= k <= 0.2 h/Mpc`, the `nk=48` direct-theta run gave

- density-derived `v_(nu-cdm) P_cb` half-pair RMS: `5.369172869387564e-4`,
- direct `theta_(nu-cdm)` half-pair RMS: `3.9505423084350484e-4`,
- direct `theta_(nu-cdm) P_cb` half-pair RMS: `5.369438898120348e-4`,
- resonant wake half-pair response: `0.23084420034493347`.

The density-derived and direct-theta-times-`P_cb` proxies agree at the `~5e-5` relative level in this run, showing that the earlier density-derived construction was not a finite-difference artifact.

A response-selected direction gave the same internal agreement at standard precision (`6.547802239846086e-4` versus `6.547375721265696e-4`), but this direction later failed precision convergence and is not retained as a quantitative control.

### k-resolved decomposition

The direct-theta diagnostic was then resolved over 96 logarithmic `k` modes. The exact identity

`delta(theta P) = bar(theta) delta(P) + bar(P) delta(theta)`

closed to absolute residuals of order `1e-16` for both tested directions. The apparently larger relative residual for the selected direction came only from division by a smaller `theta P` profile amplitude; its absolute closure residual was slightly smaller than the CREF residual.

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

The main `theta P_cb` RMS therefore changed by only about `+0.78%` from standard to moderate precision, and the wake/linear contrast changed by about `-0.77%`. The direct-theta-times-`P_cb` and density-derived proxies remain essentially identical. The three moderate-precision `k`-third RMS values are `4.504055217659639e-3`, `2.19849459583823e-3`, and `4.62497953857375e-4`.

The pure `theta` fractional RMS is more precision sensitive and is not used as the headline quantitative comparison.

For the response-selected direction, the moderate-precision result changed qualitatively:

- standard direct `theta P_cb` RMS: `6.464617200898742e-4`,
- moderate direct `theta P_cb` RMS: `1.0669648447319053e-5`,
- standard density-derived `v P_cb` RMS: `6.464779244086267e-4`,
- moderate density-derived `v P_cb` RMS: `1.0670810511705817e-5`.

The middle and high-`k` responses collapse under tighter tolerances and the number of zero crossings increases. This shows that the response-optimized direction exploited numerically delicate cancellations. Its previously quoted wake/linear factors (`~285`–`288`) are therefore retired from manuscript-level claims.

### Current retained conclusion

The robust quantitative control is the reference/CREF deformation. At `nk=96`, its density-weighted direct velocity-divergence response is `5.28e-4` at standard precision and `5.32e-4` at moderate precision, compared with a resonant wake half-pair response of `0.230844`. The corresponding integrated wake-to-linear contrast remains about `4.3e2`.

This is a transfer-level response diagnostic, not a kSZ or RSD survey forecast and not a theorem over all source-matched distributions. A final response-independent coefficient-space-orthogonal null direction is added as a separate robustness control; it is not selected using CLASS response amplitudes.

### Resource note

The aggressive publication high-precision CLASS profile is not used for this diagnostic because it repeatedly exceeded the available development-machine memory. The moderate convergence test is intentionally narrower: it changes integration tolerances without simultaneously changing hierarchy size or momentum-grid resolution.

See `docs/HIDDEN_CHANNEL_RETENTION.md` and `source_data/hidden_channel_retention_direct.json` for the current reproducible status.
