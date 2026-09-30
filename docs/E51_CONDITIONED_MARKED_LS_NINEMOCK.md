# E51 — conditioned survey-estimator contract and nine-mock technical test

E49 showed that the direct baryon velocity is strongly aligned with the
neutrino–CDM relative velocity. E50 showed that the conditioned F+/F− test is
potentially useful if amplitude information is retained. E51 therefore moves
from source physics to a concrete survey estimator.

The estimator is a **marked Landy–Szalay cross-correlation**. Every LRG–ELG
pair is assigned an external sign tag at its pair midpoint. The same external
field is evaluated in DD, DR, RD and RR:

`xi_tau = (DD_tau/NDD - DR_tau/NDR - RD_tau/NRD + RR_tau/NRR)/(RR/NRR)`.

This form is preferable to splitting the random catalogue into positive and
negative branches because the denominator remains the ordinary positive RR
support. It also subtracts tag-dependent survey geometry through the same
four-term structure.

For E51 only, four smooth fixed Cartesian sign fields are used. They are
**technical null fields, not a reconstruction of the neutrino wind**. Their
purpose is to measure how much background the new marked estimator produces
on the exact nine realistic eBOSS EZmock samples.

The script also accumulates, in the same pair pass, the linear response to a
DD-level conditioned injection

`1 + lambda tau h_F(mu)`

for the frozen E19 F+ and F− angular fingerprints. This provides an exact
technical injection basis without selecting lambda after seeing the mocks.

Important limits remain:

- the nine mocks still do not support an unrestricted 12D covariance inverse;
- the synthetic tag fields are not physically correlated with the mock density;
- no actual velocity reconstruction noise or kSZ optical-depth scatter is in E51;
- the injection is not an absolute eBOSS Einstein–Vlasov prediction;
- observed galaxies and observed odd multipoles remain sealed.

A successful E51 run answers a narrower but useful question: does the
**conditioned marked estimator itself** have a manageable mock background and
does the E19 fingerprint survive the actual eBOSS pair geometry/operator?
