# E53 — full-eligible sample-size scaling of the conditioned estimator

E51/E52 established that the F+/F− fingerprint survives the real eBOSS pair
operator, but the deliberately small 600D/1200R implementation has substantial
background. E53 tests the simplest explanation: finite-sample noise.

The science is deliberately unchanged. E53 uses exactly the same four fixed
technical sign fields and the same E19 F+/F− injection basis. It changes only
sample size.

Stage 1 uses all source-eligible mock0001 galaxies and the exact original E4
4800 randoms. This already increases the galaxy samples from 600 to the full
eligible populations while quadrupling the random sample.

Stage 2 is opt-in and uses the previously defined nested 48000 randoms, with
the same full galaxies. This separates galaxy-sampling improvement from
finite-random improvement.

Because only mock0001 is used, E53 cannot estimate mock-to-mock variance.
Its diagnostic is within-realization scaling: do the conditioned baseline
coordinates and 12D norm shrink materially from 600D/1200R -> fullD/4800R
and then stabilize from 4800R -> 48000R?

The run is fail-closed and checkpoints after every pair term. It first
rehashes the already pinned 72 mock sources, reconstructs the exact frozen
full-eligible catalogues, and reproduces the old 600D/1200R sparse/dense pair
closure. The original E4 4800 RR count and weighted sum must also replay.

No observed catalogue or observed odd vector is read.
