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

The three logarithmic `k`-third RMS values for `theta P_cb` are `4.409099614642714e-3`, `1.925236732732645e-3`, and `4.7520094922303546e-4`. The suppression therefore persists over the tested interval but is scale dependent.

### Moderate precision, `nk = 96`

The convergence check tightens only integration tolerances while keeping the ncdm hierarchy size and momentum-grid resolution unchanged.

- direct `theta_(nu-cdm)` half-pair RMS: `5.434333009379483e-4`
- direct `theta_(nu-cdm) P_cb` half-pair RMS: `5.322814846689793e-4`
- density-derived `v_(nu-cdm) P_cb` half-pair RMS: `5.322828757582342e-4`
- wake / direct `theta P_cb`: `433.68820256540243`
- wake / density-derived proxy: `433.68706914739107`

Relative to standard precision, the headline `theta P_cb` RMS shifts by about `+0.78%`, and the wake/linear contrast by about `-0.77%`. The pure `theta` fractional RMS is more precision sensitive and is not used as the headline quantitative comparison.

## Exact decomposition check

For each `k`, the code verifies

`delta(theta P) = bar(theta) delta(P) + bar(P) delta(theta)`.

The maximum absolute residuals are of order `1e-16`. This confirms that pure `theta` and `theta P_cb` are genuinely different transfer-level quantities, not two normalizations of one quantity.

## Response-selected direction: diagnostic only

At standard precision the response-selected direction gave direct `theta P_cb` RMS `6.464617200898742e-4` and wake/direct-`theta P_cb = 288.4149590565178`.

Under the moderate-precision tolerance test the same direction gives direct `theta P_cb` RMS `1.0669648447319053e-5`. The middle and high-`k` responses collapse and the zero-crossing count increases. This direction is therefore precision sensitive and is excluded from manuscript-level quantitative claims. Its earlier `~285`–`288` wake/linear ratio is retained only as development history.

## Response-independent orthogonal control

A separate null direction is constructed without using CLASS response amplitudes. The first canonical null-basis coefficient vector is Gram-Schmidt orthogonalized against normalized CREF, then normalized under the same pointwise cap and 30% deformation amplitude.

The resulting coefficient vector is

`[0.9946511694474491, 0.006310711676011596, -0.017566041166059896, -0.050820407796465655, -0.07125215639341684, -0.04975278129587607, -0.013628548959082238]`.

Its dot product with normalized CREF is `-1.9082402331634368e-17`, and the matched-moment mismatch is `2.386457994141484e-16`.

### Standard precision, `nk = 96`

- direct `theta P_cb` RMS: `1.1713198944007974e-3`
- density-derived `v P_cb` RMS: `1.1713282196477968e-3`
- wake response: `0.1101394327438084`
- wake / direct `theta P_cb`: `94.03019044609631`
- zero crossings in `theta P_cb`: `0`
- `k`-third `theta P_cb` RMS: `1.1711868750611771e-2`, `5.588669319796887e-3`, `9.56803625365203e-4`
- same-wake-numerator contrasts by `k` third: `9.4041`, `19.7076`, `115.1118`

### Moderate precision, `nk = 96`

- direct `theta P_cb` RMS: `1.3876806621063746e-3`
- density-derived `v P_cb` RMS: `1.3876847572915057e-3`
- wake response: `0.1101394327438084`
- wake / direct `theta P_cb`: `79.36943689668965`
- zero crossings in `theta P_cb`: `0`
- `k`-third `theta P_cb` RMS: `1.1763629823371717e-2`, `5.7193927604100294e-3`, `1.2065386778474646e-3`
- same-wake-numerator contrasts by `k` third: `9.3627`, `19.2572`, `91.2855`

The integrated `theta P_cb` RMS changes by about `+18.47%` from standard to moderate precision, so this second direction is not precision-converged at the percent level. The change is concentrated at high `k`: the low, middle and high thirds shift by about `+0.44%`, `+2.34%`, and `+26.10%`, respectively.

The qualitative conclusion is nevertheless stable. The response-independent control keeps zero crossings at zero, the direct and density-derived proxies agree at the `1e-5` relative level or better, and the integrated wake/linear contrast remains between about `79` and `94`.

This direction is therefore retained as a qualitative response-independent robustness control, not as a second precision-grade headline coefficient.

## Retained interpretation

The precision-grade quantitative result is the reference/CREF deformation. Its integrated density-weighted direct velocity-divergence response remains about `5.3e-4` under the standard-to-moderate precision change, while the same pair has a resonant wake response of `0.230844`. The wake-to-linear contrast therefore remains about `4.3e2`.

The orthogonal response-independent control provides a separate qualitative confirmation with an integrated contrast of about `79`–`94`, but its `~18%` precision shift prevents using one exact number as a second headline coefficient.

The result is not a theorem over all source-matched distributions, not an absolute bispectrum amplitude, and not a kSZ or RSD survey signal-to-noise forecast.

The wake fraction used here has no `k` dependence after common potential factors cancel. Bin-specific wake/linear numbers therefore reuse the same wake numerator and are not independent wake observables.

## Reproduction

```bash
python code/hidden_channel_theta_kprofile.py --direction cref --z 0.3 --nk 96 --precision standard
python code/hidden_channel_theta_kprofile.py --direction cref --z 0.3 --nk 96 --precision moderate
python code/hidden_channel_theta_kprofile.py --direction orthogonal --z 0.3 --nk 96 --precision standard
python code/hidden_channel_theta_kprofile.py --direction orthogonal --z 0.3 --nk 96 --precision moderate
```

The original aggressive high-precision profile is not part of this result because it repeatedly exceeded the available memory in the development WSL environment. The moderate profile changes only integration tolerances and is the retained precision-convergence control.
