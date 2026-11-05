# E52 — archive-only compression of the E51 conditioned estimator

E51 has two separate questions:

1. does the eBOSS pair/window operator destroy the F+/F− fingerprint?
2. is the current 600D/1200R mock implementation quiet enough to use that fingerprint?

E52 answers both without another FITS read.

## Window survival

The source-only E19 F+/F− angle is `0.02384118` rad.

After the actual E51 pair operator:

- NGC median angle = `0.02178866` rad;
- SGC median angle = `0.02344515` rad.

The post-window F−−F+ vector separation remains about
`43.56%`
of the F+ norm in NGC and
`43.49%`
in SGC.

So the survey operator does **not** erase the F-state fingerprint.

## Two one-dimensional, theory-defined compressions

For amplitude-calibrated inference define

`d = q_minus - q_plus`

and

`m_d = (d dot b)/(d dot d)`.

A technical injection of `lambda*d` shifts `m_d` by exactly `lambda`.

For amplitude-free shape-only inference first fit out the F+ amplitude,

`rho = (q_plus dot q_minus)/(q_plus dot q_plus)`

`r = q_minus - rho q_plus`

and then use

`m_r = (r dot b)/(r dot r)`.

This avoids a 12D inverse covariance. Nine mocks are enough to *describe* a
1D background, though not to make a robust tail-probability claim.

## What the current 600D/1200R mocks say

Amplitude-calibrated 83% symmetric sign-recovery thresholds in technical
`lambda` are:

- NGC: A `0.499`,
  B `0.414`,
  C `0.381`,
  D `0.335`;
- SGC: A `0.094`,
  B `0.093`,
  C `0.276`,
  D `0.211`.

The shape-only 83% thresholds are multiple units in every field/cap, confirming
that shape-only discrimination is far harder.

These `lambda` values are technical E51 pair-modulation amplitudes, not an
absolute physical eBOSS wake prediction.

## Decision

The conditioned route survives the real pair/window geometry. The immediate
limitation is the noise of the deliberately small 600D/1200R A02 samples.

The next scientifically useful test is therefore **sample-size scaling**, not
another source calculation: replay the exact same marked estimator on the
full eligible mock galaxies with denser already-defined random samples. Only
then can we tell if the remaining background is mostly finite-sample shot
noise or an intrinsic limitation of the conditioned estimator.
