# Hidden-channel observable-retention control

Status: transfer-level response diagnostic, not a survey forecast.

This note records the direct same-pair checks used to compare an integrated linear neutrino-CDM relative-velocity response with the resonant wake response for source-matched kinetic states.

## Setup

- relic mass: `m_nu = 0.06 eV`
- source matching redshift: `z_match = 1100`
- evaluation redshift: `z = 0.3`
- pointwise deformation cap: `30%`
- linear-transfer range: `0.003 <= k <= 0.2 h/Mpc`
- half-pair fractional RMS weighting: `w_k = k^3 Delta ln k`
- wake comparison: same half-pair convention, `sigma_v = 200 km/s`
- CLASS transfer output for the direct velocity check: `mPk,dTk,vTk`

The direct worker exposes `t_cdm` and `t_ncdm[0]`. The density-weighted direct control is

`X_theta(k) = [theta_ncdm(k)-theta_cdm(k)] P_cb(k)`.

An independent density-derived control reconstructs relative velocity from the redshift derivative of `delta_ncdm-delta_cdm` and multiplies it by the same `P_cb`. Their agreement is used as an internal validation.

## Reference/CREF result

### Standard precision, `nk = 96`

- direct `theta_(nu-cdm)` half-pair RMS: `4.078367440126135e-4`
- direct `theta_(nu-cdm) P_cb` half-pair RMS: `5.281598155627618e-4`
- density-derived `v_(nu-cdm) P_cb` half-pair RMS: `5.280865049182505e-4`
- wake half-pair response: `0.23084420034493347`
- wake / direct `theta P_cb`: `437.0726313189232`
- wake / density-derived proxy: `437.13330712866616`
- max relative matched-moment mismatch: `3.579686991212226e-16`

The direct-theta-times-`P_cb` and density-derived proxies agree to much better than one percent.

The three logarithmic `k`-third RMS values for `theta P_cb` are:

1. `0.003 <= k < 0.0121644 h/Mpc`: `4.409099614642714e-3`
2. `0.0121644 <= k < 0.0493242 h/Mpc`: `1.925236732732645e-3`
3. `0.0493242 <= k <= 0.2 h/Mpc`: `4.7520094922303546e-4`

The suppression therefore persists over the tested interval but is scale dependent. The wake fraction used in this diagnostic has no `k` dependence after common potential factors cancel, so a bin-specific wake/linear number would reuse the same wake numerator and is not treated as an independent wake observable.

### Moderate precision, `nk = 96`

The convergence check tightens only integration tolerances while keeping the ncdm hierarchy size and momentum-grid resolution unchanged.

- direct `theta_(nu-cdm)` half-pair RMS: `5.434333009379483e-4`
- direct `theta_(nu-cdm) P_cb` half-pair RMS: `5.322814846689793e-4`
- density-derived `v_(nu-cdm) P_cb` half-pair RMS: `5.322828757582342e-4`
- wake / direct `theta P_cb`: `433.68820256540243`
- wake / density-derived proxy: `433.68706914739107`

Relative to standard precision, the headline `theta P_cb` RMS shifts by about `+0.78%`, and the wake/linear contrast by about `-0.77%`.

The moderate-precision `theta P_cb` RMS values in the same three logarithmic `k` thirds are:

1. `4.504055217659639e-3`
2. `2.19849459583823e-3`
3. `4.62497953857375e-4`

The middle subrange is more precision sensitive than the integrated result, so no claim of a constant pointwise contrast is made.

## Exact decomposition check

For each `k`, the code verifies

`delta(theta P) = bar(theta) delta(P) + bar(P) delta(theta)`.

For CREF at standard precision, the maximum absolute decomposition residual is `9.259086553026208e-17`, corresponding to `1.845680512271886e-14` relative to the peak `theta P` profile. The selected-direction residual is of the same absolute size. The difference in the quoted relative residual is caused by the smaller selected-direction profile amplitude, not by a failure of the identity.

This confirms that pure `theta` and `theta P_cb` are genuinely different transfer-level quantities, not two normalizations of one quantity.

## Response-selected direction: diagnostic only

At standard precision the response-selected direction gave

- direct `theta P_cb` RMS: `6.464617200898742e-4`
- density-derived `v P_cb` RMS: `6.464779244086267e-4`
- wake response: `0.18644923053132711`
- wake / direct `theta P_cb`: `288.4149590565178`.

However, under the moderate-precision tolerance test the same direction gives

- direct `theta P_cb` RMS: `1.0669648447319053e-5`
- density-derived `v P_cb` RMS: `1.0670810511705817e-5`.

The middle and high-`k` responses collapse and the zero-crossing count increases. The response-selected direction is therefore precision sensitive and is excluded from manuscript-level quantitative claims. Its earlier `~285`–`288` wake/linear ratio is retained only as development history.

## Retained interpretation

The robust quantitative result is the reference/CREF deformation. Its integrated density-weighted direct velocity-divergence response remains about `5.3e-4` under the standard-to-moderate precision change, while the same pair has a resonant wake response of `0.230844`. The wake-to-linear contrast therefore remains about `4.3e2`.

This supports the limited statement that, for the tested reference source-matched deformation, the resonant kernel retains substantially more hidden kinetic information than the integrated linear relative-velocity proxy.

The result is not a theorem over all source-matched distributions, not an absolute bispectrum amplitude, and not a kSZ or RSD survey signal-to-noise forecast.

A separate response-independent coefficient-space-orthogonal null direction is provided as the final robustness control. Its definition does not use CLASS response amplitudes.

## Reproduction

Direct `theta` pair check:

```bash
python code/hidden_channel_theta_pair_check.py --direction cref --z 0.3 --nk 48
```

`k`-resolved standard and moderate checks:

```bash
python code/hidden_channel_theta_kprofile.py --direction cref --z 0.3 --nk 96 --precision standard
python code/hidden_channel_theta_kprofile.py --direction cref --z 0.3 --nk 96 --precision moderate
```

The original aggressive high-precision profile is not part of this result because it repeatedly exceeded the available memory in the development WSL environment. The moderate profile changes only integration tolerances and is the retained precision-convergence control.
