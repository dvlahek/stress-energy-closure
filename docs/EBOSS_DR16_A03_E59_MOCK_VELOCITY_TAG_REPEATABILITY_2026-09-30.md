# E59 — mock-only survey-derived velocity-tag repeatability

E58 closes the frozen random-integration gate in both caps. E59 therefore moves from synthetic external sign fields to a velocity-direction proxy reconstructed from the mock galaxy density itself.

The reconstruction uses the exact real-space Green-function kernel for a linear velocity/gravity field after spherical top-hat smoothing at R=16 Mpc/h: K(r)=r/R^3 for r<R and K(r)=r/r^3 outside. Each tracer's galaxy field is selection-subtracted with its own full mock random catalogue. LRG and ELG are reconstructed separately.

The field is evaluated at 2048 deterministic full-mock LRGxELG DD midpoints per mock/cap in the frozen 20--140 Mpc/h separation range. Primary integration extends to 256 Mpc/h; 192 versus 256 Mpc/h is a fixed cutoff-stability diagnostic. Data and random catalogues are independently split in half to measure internal reconstruction repeatability.

The preregistered screen requires in both caps: median half-half Pearson >=0.324503, median sign agreement >=0.60, median 192/256 cutoff correlation >=0.90, and median cutoff sign agreement >=0.90. If both tracers pass, choose the tracer with the larger minimum-cap median half-half Pearson.

This is not a true-velocity calibration. Shared survey geometry and cosmic modes can make split-half agreement optimistic. An E59 pass authorizes only a truth-labelled N-body/lightcone calibration before any observed velocity tag or unblinding.
