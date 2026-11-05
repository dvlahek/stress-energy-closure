# EinsteinVlasovNP — E30–E53 catch-up and conditioned-route checkpoint

**Date:** 30. 9. 2026  
**Branch:** `audit/eboss-elg-bit8-ra-orientation-20260925`  
**Scope:** work performed after the GitHub E29 checkpoint and before E53 Stage 2.  

This checkpoint synchronizes the offline/local work that had not yet been committed to GitHub. The complete reproducibility archive contains **111 files** (scripts, protocols, reports, result JSONs, runbooks and E29 local provenance), SHA256 `2cfbab2ec2bbc62939b67245e078f2e794743cea4005364711a63afd05c265c6`. The first E47 implementation is preserved only in an explicit `archive_buggy_do_not_use/` path; the corrected E47R1 is the runnable version.

## 1. E29 local two-epoch result and the physical closure sequence E30–E44

The bounded local E29 exercise produced 20 PID-supported two-epoch candidate pairs, but **zero globally certified progenitor links**. Seventeen pairs came from the first-100k early-row window and three from a bounded spatial search. Median late PID overlap is `0.9091`, median relative L1 particle-count change `3.091%`, approximate snapshot spacing `0.307957 Gyr`, median comoving centre displacement `0.2311 h^-1 Mpc`, and median finite-difference centre-speed proxy `547.5 km/s`. These are descriptive local proxies, not M200c histories, instantaneous winds or certified merger-tree trajectories. A later targeted profile read crashed Ubuntu and was not repeated.

E30 sharpened the sign-conditioned estimand: marginal `P(S=+)=P(S=-)` does not force a general selected mean to zero; the sufficient local condition is `E[S|X]=0` for the response-relevant state `X`. E32 proved that the Abacus L2 radial percentiles do not uniquely determine the Fourier profile `u(k)` and supplied a conservative enclosure theorem. E33 quantified pure M200c pseudo-evolution from the two snapshot redshifts. E34–E36 then formalized two independent non-identifiabilities: two temporal endpoints do not determine the retarded causal history, low-k profile control needs physical radial moments, and matching low moments of the initial kinetic state does not identify later free-streamed density.

E37–E41 converted those observations into quantitative closure/conditioning bounds. E39 shows that the frozen E28R1 low-mass `alpha=.4` F-/F+ contrast is only **0.6242% of |a_FD|**, while 94.22589% of the fiducial finite-band acceleration is from `k>1/Mpc`; an adversarial high-k differential model error below about **0.662% of that high-k FD contribution** would be sufficient to preserve the frozen sign. E40 proves exact non-identifiability of the E21 channel under free state-specific tracer response `chi_F`; E41 finds that the optimistic common-response E21 source geometry needs sub-per-mille response control (about 0.0912%) to keep all four long-K components separated, while the conditioned E19 ratio is much better conditioned.

E42 combined incoming-state, history, profile, Born and tracer-response errors into one operator-level sign-preservation inequality. E43 concluded that none of the required physical blocks is currently bounded tightly enough to certify the E28 F-state sign. E44 verified from official Abacus source/docs that cleaned merger-tree auxiliary products are the right next multi-epoch object, but the local inventory did not contain the required cleaned/tree product. We therefore stopped blind Abacus mining rather than repeating broad ASDF reads.

## 2. E45–E48: the original unconditioned two-point route

E45 used only the archived E3 synthetic odd injection. The injected `|A|=0.02` is small compared with the nine-mock apparent-amplitude scatter: SD `0.09495` NGC and `0.06424` SGC. Exact +/-0.02 recovery is an algebraic estimator check, not a physical sensitivity forecast.

E46 therefore defined the only valid comparison to that axis, `A_phys=(q^T y_phys)/(q^T q)`, and explicitly rejected any direct identification of the E28 halo acceleration with E45 `A`.

E47R1 built the restricted E25 linear physical-response basis after correcting the first E47 implementation (the initial source replay failed because it interpolated the already-multiplied cross product and omitted the required R16 filter). The corrected replay matches the frozen E21 source points to machine precision. The eBOSS response is nevertheless tiny: median `alpha_g0` `7.513e-08` NGC / `7.055e-08` SGC. About 83% descriptive recovery requires `|g0|` `5.55e+05` NGC and `6.65e+05` SGC.

E48 tested self-consistency of that linear response. For diagnostic `b_L=2`, `b_E=1`, `f=1`, even the most favorable anti-aligned response requires `D=max(|d_L|,|d_E|)` `1.33e+05` NGC and `1.58e+05` SGC for ~83% recovery, corresponding to one-sigma local perturbations `60.5` and `71.4`. This is far outside the linear regime. **Decision: retain the unconditioned E21/E47 channel as a negative/limitation result and do not optimize it further.**

## 3. E49–E53: velocity-conditioned route

E49 asked the cheapest prerequisite question for the E19 conditioned observable. Direct CLASS `vTk` gives a strong R16 correlation between baryon bulk velocity and the neutrino-CDM relative velocity: minimum `|r|=0.9217` across FD/F+/F-, with minimum Gaussian sign-correlation factor `0.7463`. This passed the predeclared strong gate.

E50 then separated two inference regimes. The F-/F+ conditional vector difference is about `43.59%` of the F+ norm, but the amplitude-free shape angle is only `1.366 deg`. With an illustrative velocity-reconstruction correlation 0.8, an independently amplitude-calibrated 3-sigma comparison needs underlying conditioned S/N of order 13, while shape-only needs order 238. Therefore the conditioned route is only attractive if amplitude information is retained or externally constrained.

E51 implemented the first actual survey-level **pair-midpoint externally marked cross-Landy-Szalay estimator** on the exact nine eBOSS EZmocks, using four fixed smooth synthetic sign fields only as technical null tags. All 18 cap/mock cases passed parent sample/pair replay, and the F+/F- fingerprint survives the real pair/window operator: median post-window angle `0.02179` rad NGC and `0.02345` rad SGC. No physical velocity reconstruction, observed galaxies or observed odd vector was used.

E52 compressed the fixed theory directions to one dimension, avoiding an unjustified 12D inverse covariance from only nine mocks. For the amplitude-calibrated F--F+ direction, descriptive 83% technical-lambda thresholds are NGC `0.499, 0.414, 0.381, 0.335` for fields A-D and SGC `0.094, 0.093, 0.276, 0.211`. Shape-only thresholds are multi-unit and remain much harder.

E53 Stage 1 changed **only sample size** for mock0001: full eligible galaxies plus the exact frozen 4800 randoms. It passed all parent replay gates. The fingerprint remains stable: NGC angle `0.02499` rad and separation `43.43%`; SGC angle `0.02418` rad and separation `43.49%`. Full eligible data-row counts are 7,537 LRG + 17,149 ELG in NGC and 4,159 LRG + 15,860 ELG in SGC. The baseline 12D norm drops strongly relative to E51 600D/1200R: NGC tag-field ratios A-D `0.082, 0.070, 0.116, 0.147`, SGC `0.335, 0.112, 0.140, 0.090`. This is strong evidence that much of the E51 background was finite-sample noise.

## 4. Frozen interpretation and next action

- Original observed eBOSS galaxy rows and 24D odd vector remain **SEALED**.
- No 12D/24D inverse-covariance significance claim is authorized from the nine mocks.
- E53 is still only one mock realization, so it cannot measure mock-to-mock scatter.
- E49's baryon-velocity alignment is a source-level feasibility result, not an actual survey reconstruction.
- A03 remains **PHYSICAL_UNCERTIFIED** as an absolute observed eBOSS Einstein-Vlasov prediction; A04 remains **BLOCKED** for robust inference.
- `main` must remain untouched and PR #1 remains draft/unmerged.

**Next:** run E53 Stage 2 on mock0001 with the already-defined nested 48k randoms. If full-D/4800R -> full-D/48000R changes are small, expand the **same frozen conditioned estimator** to all nine mocks. Only after that should a real external velocity-reconstruction/tag model, amplitude calibration, coverage/covariance plan and final preregistration be built. Observed odd stays sealed until those gates are passed.
