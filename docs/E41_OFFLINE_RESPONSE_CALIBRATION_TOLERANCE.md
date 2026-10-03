# E41 — how accurately must the unknown tracer response be controlled before the frozen F-state source differences cannot be absorbed?

This is a **source-level robustness calculation**, not an eBOSS forecast. It uses only frozen E21 and E19 values.

For one positive E21 component, suppose the physical tracer response multiplying the source is uncertain by a symmetric fractional amount `±eta` under each F-state hypothesis. The F+ and F− prediction intervals remain disjoint only if

`eta < (C_plus-C_minus)/(C_plus+C_minus)`.

For the four frozen long modes `K=.001,.002,.003,.005 h/Mpc`, the thresholds are:

- 0.3561%
- 0.2404%
- 0.0942%
- 0.0912%.

Thus keeping **all four** source components individually non-overlapping would require response control below about **0.0912%** in this deliberately simple interval model. At least one component remains non-overlapping only below the loosest threshold, 0.3561%.

The conditioned E19 source uses the amplitude-free ratio `T3/T1`: F+ = 0.635144094 and F− = 0.669123978. Under the same symmetric multiplicative interval logic for this ratio, non-overlap persists for

`eta < 0.026053`

or about **2.605%**.

So the E19 conditional angular ratio is about **28.6 times less ill-conditioned** by this narrow source-level criterion than demanding componentwise separation in all four E21 long-K points.

This does **not** establish observational superiority. E19 requires an independently reconstructed/conditioned neutrino-CDM wind and a new estimator; E21 still lacks the physical high-z LRG/ELG response `chi_F`. The point is that the original unconditioned linear source difference is small enough that sub-percent response uncertainty can absorb it, consistent with E40's near-collinearity result.
