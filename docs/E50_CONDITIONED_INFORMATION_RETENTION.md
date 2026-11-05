# E50 — what E49 means for the conditioned F+/F− test

E49 passed the strong predeclared alignment gate: the minimum direct-CLASS baryon/relative-velocity correlation across FD/F+/F− is about `0.9217`. The ideal sign tag therefore preserves about 75% of a sign-weighted conditioned signal.

The frozen E19 vectors are

`F+ = (0.000222788767, 0.00014150297)`

`F- = (0.00031476164, 0.000210614561)`.

They differ in two useful ways:

- absolute vector separation is `43.59%` of the F+ norm;
- after arbitrary amplitude normalization, the residual shape rotation is only `1.366 deg`, equivalent to the known ~5.35% difference in `T3/T1`.

A noisy velocity sign tag multiplies both odd multipoles by the same sign-correlation attenuation. Therefore **the ratio `T3/T1` and normalized F+/F− angle survive tag noise**; the price is only lower total SNR.

This produces two sharply different scientific regimes.

## 1. Amplitude calibrated independently

If the conditional amplitude is fixed by independently determined tracer bias, `P_cb`, and velocity-tag normalization, then the full F+/F− vector difference can be used. The required 3-sigma underlying conditional-signal SNR is approximately:

- ideal baryon sign tag: `9.2`
- reconstruction correlation 0.8: `13.0`
- reconstruction correlation 0.7: `15.4`
- reconstruction correlation 0.5: `22.6`.

This is not obviously impossible.

## 2. Overall amplitude free

If a free amplitude is marginalized, only the normalized shape difference remains. The F+/F− angle is `0.023841 rad`, so a 3-sigma shape-only discrimination needs tagged-signal SNR about `125.8`. After tag dilution this becomes:

- ideal tag: `169`
- r_rec=0.8: `238`
- r_rec=0.7: `282`
- r_rec=0.5: `413`.

That is much harder.

Therefore E50 does not yet forecast eBOSS significance. It identifies the decisive design choice **before** building the mock estimator: an amplitude-calibrated conditioned test may be viable, while a pure amplitude-free shape fingerprint requires very high SNR.

The next valid survey step is a pre-registered conditioned Landy–Szalay estimator with an external velocity-tag field evaluated at galaxies and random positions, first on mocks only. The two inference modes must be frozen separately before looking at observations.
