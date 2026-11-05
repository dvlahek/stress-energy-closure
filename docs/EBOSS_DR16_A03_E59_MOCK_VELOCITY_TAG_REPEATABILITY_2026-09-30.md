# E59 — mock-only survey-derived velocity-tag repeatability (v2)

E58 closes the frozen pair-estimator random-integration gate in both caps. E59 moves from synthetic external sign fields to a velocity-direction proxy reconstructed from mock galaxy density.

The reconstruction uses the exact real-space Green-function kernel for a linear velocity/gravity field after spherical top-hat smoothing at R=16 Mpc/h. LRG and ELG are reconstructed separately.

For computationally bounded E59, the selection field is the exact frozen E58 R1 48k random subset for each mock/cap/tracer. The SAME fixed random field is subtracted from both galaxy halves. Thus split-half repeatability isolates galaxy sampling noise conditional on that numerical mask. It is not a claim that one 48k random realization is a converged physical reconstruction mask.

The field is evaluated at 512 deterministic full-mock LRGxELG DD midpoints in 20--140 Mpc/h. Primary integration extends to 256 Mpc/h; 192 versus 256 Mpc/h is the fixed cutoff-stability diagnostic.

Passing E59 authorizes only a truth-labelled N-body/lightcone calibration, including reconstruction-mask robustness. It does not authorize observed galaxy rows or the sealed odd vector.
