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

The direct response itself was numerically stable under changes in redshift finite-difference step, k-grid density and denominator mask. The failure was in the extrapolated operator optimization, not in the final direct endpoint response.

### Final OOM-safe direct pair checks

The decisive comparison was repeated with `code/hidden_channel_pair_check.py`, evaluating one direction and one redshift at a time with short-lived CLASS subprocesses.

At `m_nu = 0.06 eV`, `z = 0.3`, a 30% pointwise deformation cap and `0.003 <= k <= 0.2 h/Mpc`:

| direction | linear proxy half-pair RMS | wake half-pair response | wake / linear |
|---|---:|---:|---:|
| reference/CREF | `5.369172869387564e-4` | `0.23084420034493347` | `429.94369144099613` |
| independently selected direction | `6.547802239846086e-4` | `0.18644923053132711` | `284.75085792406253` |

Both pairs preserve the matched source moments to about `3.58e-16` relative accuracy. An independent 96-mode reference screening result (`~5.28e-4`) agrees with the 48-mode direct checkpoint at the percent level.

The resulting manuscript-level statement is limited and direct: for two distinct tested source-matched kinetic deformations, resonant wake sampling retains roughly 285–430 times more fractional response than the integrated linear relative-velocity proxy under the same half-pair convention. This is a response-level diagnostic, not a survey signal-to-noise forecast and not a theorem over all distributions.

### Resource note

A more aggressive high-precision CLASS profile repeatedly exceeded the available memory in the development WSL environment and was terminated by the OS. It is not part of the quoted result. The retained standard-precision values are supported by the independent stability checks above.

See `docs/HIDDEN_CHANNEL_RETENTION.md` and `source_data/hidden_channel_retention_direct.json` for the final reproducible numbers.
