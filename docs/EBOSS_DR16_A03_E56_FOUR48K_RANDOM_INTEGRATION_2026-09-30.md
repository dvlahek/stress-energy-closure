# E56 — four-replica random integration over the nine full-eligible mocks

E55 completed the full nine-mock conditioned background, but its frozen
single-replica random-integration gate failed globally because NGC retained too
much 48k-random Monte-Carlo variation. E56 does not change the scientific
estimator. It only adds R3/R4 and uses a four-replica random mean.

For each cap and tag field, E56 estimates the pooled within-mock random variance
from 9 mocks x 4 independent subset draws, converts it to the Monte-Carlo
standard error of the four-replica mean, and compares that error with the
between-mock SD of the four-replica mock means.

The gate is not chosen from the E55 numerical outcome. It is inherited
algebraically from the E55 preregistration: the E55 single-rep thresholds
0.5 (median across fields) and 1.0 (maximum field) are divided by sqrt(4).
Therefore E56 requires <=0.25 median and <=0.50 maximum in EACH cap.

A methodological clarification is recorded before E56: separately seeded
without-replacement subsets are independent random draws conditional on the
fixed eligible pool even if they happen to reuse some of the same rows.
Realized overlap therefore does not itself invalidate replica averaging.
However, all replicas share the same finite parent pool, so they cannot measure
a common finite-pool approximation bias relative to an ideal/infinite random
catalogue.

Passing E56 closes only the numerical random-integration gate. It does not
authorize unblinding. The next gate would be a real survey velocity-tag model
and absolute conditioned-amplitude calibration, with the observed odd vector
still sealed.
