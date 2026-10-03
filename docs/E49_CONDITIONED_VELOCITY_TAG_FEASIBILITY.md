# E49 — velocity-tag feasibility before building the conditioned estimator

E48 closes the original unconditioned linear channel for ordinary order-unity tracer geometry: the response required to rise above the E45 mock background is far outside the linear local-response regime.

The next route is the E10/E19 conditional statistic. Its key obstacle is operational, not algebraic: the original eBOSS 24D observable has no neutrino–CDM wind-sign tag.

E49 therefore asks the cheapest prerequisite question first:

**Could an independent baryon/electron bulk-velocity reconstruction provide a useful sign proxy for the neutrino–CDM relative velocity?**

The local E49 script uses direct CLASS `vTk`, not the older density-derivative proxy, and computes the R16-filtered field correlation

`r = <v_rel v_b> / sqrt(<v_rel^2><v_b^2>)`

at z=.95 separately for frozen F0/F+/F−. It must first reproduce the archived E21 direct-vTk R16 LOS velocity dispersion in every F state.

For a zero-mean jointly Gaussian pair of velocity fields, the ideal sign-tag correlation is

`E[sign(v_rel) sign(v_tag)] = (2/pi) asin(r)`,

so E49 also reports the ideal correct-sign probability and illustrative degradation if a survey velocity reconstruction has correlation 0.8, 0.7 or 0.5 with the true baryon velocity.

The predeclared decision is:

- minimum `|r_rel,b| >= .7`: continue to a conditioned mock estimator;
- `.4 <= minimum |r_rel,b| < .7`: possible, but reconstruction noise is central;
- minimum `< .4`: do not build the conditioned survey analysis around a baryon/kSZ sign tag.

No kSZ observation, eBOSS observed galaxy rows, observed odd vector, FITS or ASDF is used.
