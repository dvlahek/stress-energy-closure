# Hidden-channel observable-retention control

Status: direct response diagnostic, not a survey forecast.

This note records the final direct-pair checks used to compare an integrated linear neutrino-CDM relative-velocity proxy with the resonant wake response for the same source-matched kinetic states.

## Setup

- relic mass: `m_nu = 0.06 eV`
- source matching redshift: `z_match = 1100`
- evaluation redshift: `z = 0.3`
- pointwise deformation cap: `30%`
- linear-transfer range: `0.003 <= k <= 0.2 h/Mpc`
- direct checkpoint runs: `nk = 48`
- linear proxy: `X(k) = v_(nu-cdm)(k) P_cb(k)`
- half-pair fractional RMS weighting: `w_k = k^3 Delta ln k`
- wake comparison: same half-pair convention, `sigma_v = 200 km/s`

The CLASS transfer output used in these runs does not expose direct non-cold/CDM velocity-divergence fields, so the relative velocity is reconstructed from the redshift derivative of `delta_ncdm - delta_cdm`. No theta-based cross-check is claimed.

## Direct results

| direction | linear proxy half-pair RMS | wake half-pair response | wake / linear | max relative matched-moment mismatch |
|---|---:|---:|---:|---:|
| reference/CREF | `5.369172869387564e-4` | `0.23084420034493347` | `429.94369144099613` | `3.579686991212226e-16` |
| independently selected direction | `6.547802239846086e-4` | `0.18644923053132711` | `284.75085792406253` | `3.5796869912122257e-16` |

An independent 96-mode screening calculation for the reference direction gave a linear proxy response of about `5.28e-4` at the same redshift, within roughly 1.7% of the 48-mode checkpoint value.

The second direction came from an independent linear-response search. Its small-probe linearized optimizer substantially overpredicted the absolute final response, so the optimizer prediction is not used quantitatively. Only the direct full-endpoint CLASS values above are retained.

## Interpretation

For two distinct source-matched deformations, the resonant wake retains about `2.8e2` to `4.3e2` times more fractional hidden-state response than the integrated linear relative-velocity proxy under the same half-pair normalization. This supports the manuscript claim that parity sensitivity alone does not determine information retention; the observable kernel matters.

The result is restricted to the tested smooth matched-moment class and is not a theorem over all distributions. It is also not an absolute bispectrum or survey signal-to-noise forecast.

## Reproduction

The OOM-safe direct checks are run one pair at a time:

```bash
python code/hidden_channel_pair_check.py --direction cref --z 0.3 --nk 48
python code/hidden_channel_pair_check.py --direction selected --z 0.3 --nk 48 \
  --pair-csv hidden_channel_operator_diagnostics/best_proxy_pair.csv
```

The heavy high-precision profile is not part of the quoted result. On the development WSL environment it exceeded the available memory, while the standard-precision response had already passed independent `dz`, k-grid and denominator-mask stability checks.
