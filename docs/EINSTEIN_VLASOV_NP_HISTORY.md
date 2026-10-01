# EinsteinVlasovNP — povijest i handoff za novi chat

## 1. 10. 2026. — E62 pre-run reader fix: vendored pure-Python readfof

The first local E62 environment setup attempted `pip install Pylians` under Python 3.12 and failed while compiling unrelated Cython extensions (notably `HI_clusters_library`) because the build environment pulled an incompatible Cython/NumPy combination. This failure occurred before any Quijote data read or E62 truth calculation.

Inspection of upstream Pylians3 confirms that `library/readfof.py`, the only component E62 needs, is a standalone pure-Python reader depending on NumPy plus standard-library modules. The audit branch therefore vendors the upstream MIT-licensed reader as `scripts/readfof.py` and E62 now imports that file directly. No E62 scientific parameter, dataset choice, halo count, smoothing scale, cutoff, probe seed, gate threshold or observed-data guardrail changed. Full Pylians installation is no longer required for E62.

## 30. 9. 2026. — E61 full-pool lock PASS; E62 truth-velocity method gate frozen

E61 completed all 18 mock/cap full-pool reconstructions using the exact E59 probes and kernel. The complete source-eligible ELG random pool is now the final deterministic survey-selection operator. The seven-replica E60 mean is diagnostic only, but it is an excellent surrogate for the full-pool reconstruction across all nine mocks in both caps.

NGC aggregate seven-replica mean vs full-pool metrics: Pearson median 0.9919635, minimum 0.9877627; sign-agreement median 0.9628906, minimum 0.9492188; RMS-difference/full-pool-field-RMS median 0.1297709, maximum 0.1606560. SGC: Pearson median 0.9920690, minimum 0.9903052; sign-agreement median 0.9589844, minimum 0.9394531; RMS ratio median 0.1264548, maximum 0.1402768. Both prospectively frozen surrogate gates pass. E61 final decision is `FULLPOOL_ELG_RECONSTRUCTION_LOCKED; SEVENREP_SURROGATE_VALID; NEXT_TRUTH_VELOCITY_CALIBRATION`. Observed galaxy rows and observed odd were not used; true velocity and absolute EV amplitude remain uncalibrated.

For the next stage, do not reinterpret internal mock repeatability as truth. The current eBOSS EZmock galaxy files do not provide the required true peculiar-velocity labels. Earlier local Abacus metadata confirmed `x_L2com` and `v_L2com`, but a single superslab cannot support the unchanged 256 Mpc/h E59 kernel without missing-neighbour/boundary bias. To avoid a large multi-superslab download and a geometrically incomplete calibration, E62 is frozen as a lightweight independent N-body method calibration on one complete periodic Quijote box: fiducial realization 0, FoF snapshot 2 (z=1).

E62 uses exactly the 165107 most massive FoF halos, unit number weights, R=16 Mpc/h smoothing, 256/192 Mpc/h cutoffs, periodic minimum-image geometry, and 4096 deterministic probe halos (seed 202609620000). The primary truth metrics are reconstructed-vs-true halo velocity Pearson and sign agreement separately for x/y/z, with frozen per-axis thresholds 0.70/0.70, plus 192-vs-256 cutoff stability thresholds 0.90/0.90. No mass-threshold or kernel tuning after truth is allowed. Passing E62 authorizes only a later end-to-end survey/lightcone truth-transfer test. It does not validate an ELG HOD, the custom F+/F− neutrino states, or an absolute Einstein–Vlasov amplitude. Observed eBOSS rows/odd remain SEALED.

Artifacts synchronized to the audit branch:
- `source_data/e61_elg_fullpool_reconstruction_lock_compact_summary_2026-09-30.json`
- `source_data/e62_quijote_z1_truth_velocity_calibration_prereg_2026-09-30.json`
- `scripts/e62_quijote_z1_truth_velocity_calibration.py`
- `docs/E62_QUIJOTE_Z1_TRUTH_VELOCITY_CALIBRATION_2026-09-30.md`

## 30. 9. 2026. — E60 mask robustness fails; full-pool operator chosen

E60 completed the ELG-only random-mask robustness test across all nine fixed EZmocks and both caps. The frozen individual-48k-replica gate fails in both caps. NGC cap medians: pairwise Pearson 0.889595 < 0.95, pairwise sign agreement 0.849609 < 0.90, and random-replica RMS / seven-replica-mean field RMS 0.346451 > 0.25 (maximum 0.370690). SGC: 0.899388, 0.843750, and 0.332598 (maximum 0.384874), also failing. Therefore one 48k random mask, and the individual 48k-replica stability criterion, are not adequate for the final velocity reconstruction.

At the same time the direct complete-random-pool anchor on mock0001 is excellent in both caps: NGC seven-replica mean vs full eligible ELG random pool gives Pearson 0.991964, sign agreement 0.966797 and RMS difference / full-pool field RMS 0.129771; SGC gives 0.992069, 0.970703 and 0.126455. E60 remains a preregistered FAIL and is not retroactively relabelled.

Downstream consequence: use the complete source-eligible ELG random pool as the deterministic survey-selection operator. E61 is frozen as an all-18-case full-pool lock at the same E59 probes/kernel. Seven-replica fields are diagnostic only. Observed galaxies and observed odd remain SEALED. Truth-velocity calibration remains the next physical gate.


## 30. 9. 2026. — E59 selects ELG; E60 mask-robustness gate frozen

E59 completed all 18 mock/cap cases. Under the prospectively frozen cap-median screen, **eBOSS_ELG passes and eBOSS_LRG fails**. ELG minimum-cap medians are split-half Pearson 0.455661, sign agreement 0.667969, 192/256 cutoff Pearson 0.999202 and cutoff sign agreement 0.96875. NGC/SGC ELG split-half Pearson medians are 0.487084/0.455661; sign-agreement medians are 0.675781/0.667969. The equal-independent-noise diagnostic gives implied full-vs-latent r medians 0.809374/0.791236, but this remains diagnostic only and is not correlation with true baryon/halo velocity. LRG minimum-cap medians are only 0.176041 in split-half Pearson and 0.564453 in sign agreement. The frozen selection rule therefore carries ELG forward. One ELG individual case (mock0125/NGC) has sign agreement 0.585938 below 0.60, but the preregistered decision rule is cap-median based; 17/18 ELG cases satisfy all four individual thresholds, while 0/18 LRG cases do.

E59 still shares one fixed E58 R1 48k random field between split halves. E60 is therefore frozen before any E60 output as an ELG-only random-mask robustness test using exact E58 R1-R7 48k subsets, the same E59 probes, and the same R=16 Mpc/h kernel. The cap gates are median pairwise Pearson >=0.95, median pairwise sign agreement >=0.90, median random-RMS/mean-field-RMS <=0.25, and maximum case ratio <=0.50. Mock0001 in both caps additionally uses the complete eligible ELG random pool as a common-bias anchor, requiring Pearson >=0.98, sign agreement >=0.95, and relative RMS difference <=0.25. Passing E60 authorizes only the next truth-labelled velocity-calibration stage. No observed galaxy rows or observed odd vector are opened.



## 30. 9. 2026. — E59 pre-run implementation refinement v2

Before any E59 output was produced, the local implementation was reduced to a computationally bounded version while preserving the same physical question and the already-frozen E59 thresholds. The exact R=16 Mpc/h top-hat-smoothed real-space velocity kernel and 192/256 Mpc/h cutoff check are unchanged. The probe count is fixed at 512 deterministic full-mock LRGxELG DD midpoints per mock/cap. For survey-selection subtraction, E59 v2 uses the exact E58 R1 48k random subset for each mock/cap/tracer and subtracts that same fixed random field from both galaxy halves. Therefore split-half agreement measures galaxy-sampling repeatability conditional on one fixed numerical mask; it does not establish random-mask convergence or truth-velocity correlation. The preregistered thresholds remain split-half Pearson >=0.324503, sign agreement >=0.60, cutoff Pearson >=0.90 and cutoff sign agreement >=0.90 in both caps. If E59 passes, the next gate remains truth-labelled N-body/lightcone calibration, explicitly including mask/reconstruction robustness. No observed rows or sealed odd vector are opened.


## 30. 9. 2026. — E58 random-integration closure and E59 velocity-tag reconstruction preregistration

E58 completed the seventh 48k random replica for every fixed mock/cap, reusing E57 R1-R6 and DD. The unchanged final-mean gate now passes in both caps. NGC field ratios A-D are 0.144669, 0.171511, 0.290179, 0.296087, with median/max 0.230845/0.296087. SGC ratios are 0.118636, 0.202948, 0.262174, 0.252714, with median/max 0.227831/0.262174. The frozen limits remain median <=0.25 and max <=0.50. The F-state fingerprint remains stable: mean ||q_- - q_+||/||q_+|| = 0.434823 NGC and 0.434951 SGC; mean post-window angle = 0.024213 and 0.024058 rad. E58 final decision: `SEVENREP_RANDOM_INTEGRATION_GATE_CLOSED; ADVANCE_TO_PHYSICAL_VELOCITY_TAG_AND_ABSOLUTE_AMPLITUDE_CALIBRATION`. This closes only conditional random-subset Monte-Carlo precision, not finite-parent-pool bias, physical tagging, absolute EV amplitude, covariance inference, or significance.

E59 is preregistered before its result. It moves from the four synthetic E51 sign fields to a mock-survey-derived linear velocity-direction proxy, reconstructed separately from LRG and ELG density fields with their own full random selections. The kernel is the exact real-space top-hat-smoothed linear velocity/gravity kernel at R=16 Mpc/h, evaluated to 256 Mpc/h with a fixed 192 Mpc/h cutoff-stability check. Probes are 2048 deterministic full-mock LRGxELG DD midpoints per cap/mock in 20--140 Mpc/h. Split-half repeatability is the primary internal diagnostic. The fixed screen requires, in both caps, median half-half Pearson >=0.324503, sign agreement >=0.60, cutoff correlation >=0.90 and cutoff sign agreement >=0.90. If more than one tracer passes, select the candidate maximizing the minimum-cap median half-half Pearson. If neither passes, do not use an observed tag. A pass authorizes only a subsequent truth-labelled N-body/lightcone calibration; split-half agreement is not true-velocity correlation and can be optimistic because the halves share mask and cosmic modes. Observed galaxy rows and the observed odd vector remain SEALED.

## 30. 9. 2026. — E30–E53 catch-up, unconditioned STOP and conditioned-route pivot

This checkpoint records the local/offline work that was performed after the GitHub E29 checkpoint and had not yet been synchronized to the audit branch.

- **E29 local two-epoch audit:** 20 PID-supported candidate pairs were examined, but zero were promoted to globally certified progenitor links. Median PID overlap was 0.9091, median relative L1 particle-count change 3.091%, median comoving centre displacement 0.2311 h^-1 Mpc, and median finite-difference centre-speed proxy 547.5 km/s. These are descriptive proxies, not M200c histories, instantaneous winds or certified merger-tree trajectories.
- **E30–E44 physical closure:** E30 sharpened the conditional-sign estimand to the sufficient condition E[S|X]=0. E32 showed radial percentiles do not uniquely determine u(k). E33 quantified M200c pseudo-evolution. E34–E38 formalized temporal/profile/initial-kinetic non-identifiability and closure bounds. E39 found the frozen E28R1 low-mass alpha=.4 F-/F+ contrast is only 0.6242% of |a_FD| while 94.22589% of the finite-band reference comes from k>1/Mpc. E40–E43 showed response/source/operator closure is not currently strong enough to certify the E28 sign. E44 stopped further blind Abacus ASDF mining because the required cleaned/tree multi-epoch object was not locally available.
- **E45–E48 unconditioned route:** E45 synthetic injection A=0.02 is small compared with the nine-mock apparent-amplitude scatter (SD 0.09495 NGC, 0.06424 SGC). E46 defined the exact physical-to-E45 projection and rejected direct E28 acceleration -> E45 amplitude identification. E47R1 corrected the original E47 interpolation/R16 implementation, replayed E21 to machine precision, and found median physical transfer coefficients of order 7e-8 (g0) and 4e-8 (g2); about 83% descriptive recovery requires |g0| ~5.55e5 NGC and 6.65e5 SGC. E48 showed the required order-unity tracer response violates the linear regime by tens to hundreds at one-sigma relative velocity. **Decision: retain E21/E47 as a negative/limitation result and do not optimize the original unconditioned linear 2pt channel further.**
- **E49–E53 conditioned route:** E49 direct CLASS vTk gives min |r(v_nu-v_cdm,v_b)|=0.92165 across F0/F+/F-, passing the preregistered strong velocity-tag feasibility gate. E50 shows amplitude-calibrated conditioned discrimination is much easier than amplitude-free shape-only discrimination. E51 implements the first pair-midpoint externally marked cross-Landy-Szalay estimator on the exact nine eBOSS EZmocks; all 18 cap/mock cases replay their parent samples/pairs and the F+/F- fingerprint survives the real pair/window operator (median angle 0.02179 rad NGC, 0.02345 rad SGC). E52 defines theory-fixed 1D compressions and confirms shape-only is much noisier. E53 Stage 1 changes only sample size for mock0001 to full eligible galaxies + frozen 4800R: NGC angle 0.02499 rad, separation 43.43%; SGC angle 0.02418 rad, separation 43.49%. Baseline 12D norms fall strongly relative to E51 600D/1200R (NGC A-D ratios 0.082,0.070,0.116,0.147; SGC 0.335,0.112,0.140,0.090), indicating much of the E51 background was finite-sample noise.
- **Frozen interpretation:** observed eBOSS galaxy rows and the original 24D odd vector remain SEALED. Nine mocks do not authorize unrestricted 12D/24D inverse-covariance inference. E49 is a source-level tag-feasibility result, not an actual survey velocity reconstruction. A03 remains PHYSICAL_UNCERTIFIED for an absolute observed eBOSS EV prediction; A04 remains BLOCKED for robust inference. main remains untouched and PR #1 remains draft/unmerged.
- **Next:** E53 Stage 2 on mock0001 with the already-defined nested 48000 randoms. If full-D/4800R -> full-D/48000R is stable, expand the exact same frozen conditioned estimator to all nine mocks before any observed-data access.


### E55 preregistration — nine full-eligible mocks

After E54 passed the 48k random-floor engineering gate, E55 is frozen before any E55 result. It applies the same E51/E19 conditioned estimator to all nine fixed EZmock IDs with all source-eligible galaxies and two fresh 48k random replicas per cap/tracer. R1/R2 are averaged only after separate marked-LS evaluation. The primary output is the 1D calibrated F--F+ coordinate. A diagnostic gate compares the paired random-integration scale with the between-mock SD: median ratio <=0.5 and every field <=1.0. Nine mocks remain descriptive only; no p-values, covariance inversion, observed rows/odd access or physical velocity reconstruction are authorized.


### E57 — six-replica nine-mock random integration

E57 reused E56 R1-R4 and added R5/R6 for all nine frozen mock IDs in both caps, retaining the E51 tag fields, E19 F+/F- basis, all source/selection cuts and the E56 final-mean engineering gate unchanged. The run completed with status `PASS_NINEMOCK_FULL_ELIGIBLE_SIX48K_BACKGROUND_QUANTIFIED`; observed galaxy rows, the observed odd vector, physical velocity reconstruction, covariance inversion and p-value/detection calculations were all absent.

NGC effective-random/between-mock ratios A-D = 0.163631, 0.190350, 0.312476, 0.314130, giving median 0.251413 and max 0.314130. The max criterion <=0.50 passes, but the median misses <=0.25 by 0.001413, so NGC formally fails. SGC ratios = 0.129260, 0.221248, 0.299603, 0.269269, giving median 0.245259 and max 0.299603; SGC passes. Global E57 interpretation therefore remains `SIXREP_RANDOM_INTEGRATION_STILL_MATERIAL; DO_NOT_ADVANCE_PHYSICAL_TAG_GATE`.

The F-state geometry remains highly stable across the nine mocks. NGC mean ||q_- - q_+||/||q_+|| = 0.434819 (SD 2.09e-4), mean angle = 0.024217 rad (SD 3.65e-4). SGC = 0.434950 (SD 1.02e-4), angle = 0.024058 rad (SD 1.54e-4). Thus the remaining issue is numerical random integration, not collapse of the conditioned F+/F- fingerprint.

E58 is preregistered before its result. It adds exactly R7 to every mock/cap and reuses E57 R1-R6 and DD. R=7 is the minimum integer obtained from the unchanged NGC gate under standard 1/sqrt(R) scaling: ceil[6*(0.2514132726578298/0.25)^2]=7. The gate stays median <=0.25 and field max <=0.50 in each cap. If E58 fails, stop replica-by-replica chasing and use deterministic or higher-density random integration. Observed rows/odd remain SEALED.


### E56 — four-replica nine-mock random integration

E56 reused E54/E55 parents and completed the frozen conditioned estimator with four 48k random replicas for all nine mock IDs in both caps. Provenance/observation guardrails remained closed. The F+/F- geometry stayed stable: NGC mean relative separation 0.434833 (SD 2.10e-4) and mean post-window angle 0.024202 rad (SD 3.75e-4); SGC 0.434952 (SD 1.07e-4) and 0.024050 rad (SD 1.52e-4).

The already-preregistered E56 final-mean random-integration gate was not met, but only narrowly on its median criterion. NGC effective-random/between-mock ratios A-D = 0.232, 0.186, 0.443, 0.346, giving median 0.289 (>0.25) and max 0.443 (<0.50). SGC = 0.149, 0.224, 0.424, 0.353, giving median 0.289 (>0.25) and max 0.424 (<0.50). Final E56 interpretation: `FOURREP_RANDOM_INTEGRATION_STILL_MATERIAL; USE_HIGHER_REPLICA_OR_DETERMINISTIC_RANDOM_INTEGRATION`.

Standard 1/sqrt(R) scaling of the measured E56 random component gives projected five-replica medians ~0.259 in both caps, still above target, and six-replica medians ~0.236 with projected maxima ~0.362 (NGC) and ~0.346 (SGC). Therefore E57 is preregistered as the minimum-integer six-replica engineering closure: reuse R1-R4 exactly, add R5/R6 only, and keep the E56 final-mean gate unchanged at median <=0.25 and per-field max <=0.50 in each cap. This uses no observed galaxy rows or observed odd vector and does not alter science cuts, E19 or E51.

### E55 — nine full-eligible mocks with paired 48k randoms

E55 completed all nine frozen EZmock IDs in both NGC and SGC using all source-eligible galaxies and two fresh 48k random replicas per tracer/cap, with the E51 tag fields and E19 F+/F- basis unchanged. The run passed its computational/provenance gates and kept observed galaxy rows and the observed odd vector sealed. The F-state geometry remains extremely stable across mocks: mean ||q_- - q_+||/||q_+|| = 0.434831 (NGC, SD 2.03e-4) and 0.434953 (SGC, SD 1.23e-4); mean post-window angle = 0.024214 rad (NGC, SD 3.58e-4) and 0.024045 rad (SGC, SD 1.37e-4).

The preregistered random-integration-subdominance gate **does not pass globally**. NGC median paired-random/between-mock ratio = 0.6141 (>0.5), max = 1.0643 (>1.0), driven especially by field C; field ratios A-D = 0.522, 0.291, 1.064, 0.706. SGC passes its cap-level gate with median = 0.4645 and max = 0.9364; field ratios = 0.244, 0.420, 0.936, 0.509. Therefore the E55 decision is `RANDOM_INTEGRATION_NOT_SUBDOMINANT; INCREASE_REPLICA_AVERAGING_BEFORE_PHYSICAL_TAG_GATE`. Nine mocks are still descriptive 1D background only, with no robust tail probability, covariance inverse, detection or exclusion.

Post-E55 methodological correction: independently generated without-replacement 48k subsets are independent Monte-Carlo draws **conditional on the fixed eligible random pool**, even if the realized subsets overlap. Overlap quantifies reuse/finite-pool coverage but does not by itself create statistical dependence between independently seeded subset draws. Replica averaging and SD/sqrt(R) are therefore valid for conditional random-subset Monte-Carlo noise; they do not measure any common bias from the finite eligible pool relative to an ideal/infinite random catalogue. This correction does not retroactively change the frozen E55 gate or its fail decision.

Next gate: E56 adds R3/R4 to every non-0001 mock/cap, reuses E54 R1-R4 for mock0001 and E55 R1/R2 elsewhere, and evaluates the four-replica mean. The effective random-integration standard error of the four-replica mean will be compared to between-mock scatter with the threshold derived prospectively from the E55 rule: median ratio <=0.25 and per-field max <=0.50. Observed rows/odd remain SEALED.

### E54 — four fresh 48k random replicas

E54 held mock0001 full eligible galaxies, the E51 four technical sign fields, the E19 F+/F- injection basis and the pair/window operator fixed, and changed only the 48k random Monte-Carlo realization. Four fresh deterministic 48k replicas per cap/tracer were drawn from the same frozen source-eligible random pools. The preregistered engineering gate passed in both caps. NGC: calibrated F--F+ coordinate SDs A-D = 0.00222, 0.01291, 0.00352, 0.00436; median 0.00394, max 0.01291; angle spread 1.92e-4 rad; relative-separation spread 9.49e-5. SGC: SDs = 0.00762, 0.00694, 0.01746, 0.00682; median 0.00728, max 0.01746; angle spread 1.72e-4 rad; relative-separation spread 1.13e-4. Therefore the 48k random floor is small on the unit technical-lambda scale and the conditioned estimator is cleared for a nine-mock **engineering expansion**. Important limits remain: this is not a physical eBOSS detection threshold, the four replicas share a finite source pool, and SGC LRG replicas overlap by about 51.4% pairwise because only 93,220 eligible random rows exist. The original E53 nested 48k realization is atypical relative to the fresh NGC replicas for several tag directions, so downstream nine-mock runs should use the preregistered fresh-replica averaging/summary policy rather than treating one nested 48k realization as canonical noise-free integration.

### E53 Stage 2 — full-D / 48000R

E53 Stage 2 completed on mock0001 for NGC and SGC with the already-defined nested 48000 randoms. The run passed with observed galaxies/odd still sealed and no physical velocity reconstruction or p-value inference. Signal geometry is highly stable: NGC post-window F+/F- angle 0.0245001 rad and ||q_- - q_+||/||q_+||=0.434696; SGC 0.0240902 rad and 0.435085. However, the finite-random background is **not yet converged at 4800R**. Increasing 4800R -> 48000R reduces the baseline 12D norm by factors A-D = 0.352,0.443,0.328,0.351 in NGC (median 0.351) and 0.294,0.536,0.448,0.798 in SGC (median 0.492). The calibrated F--F+ coordinate also drops substantially: NGC ratios 0.462,0.589,0.513,0.569 (median 0.541); SGC 0.376,0.416,0.178,0.185 (median 0.280). Therefore the predeclared criterion for immediate nine-mock expansion is not satisfied. Interpretation: the conditioned fingerprint survives and the background continues to fall with random density, but finite-random noise is still material. Next gate should quantify the random-noise floor with preregistered independent high-density random replicates (or an equivalent deterministic/infinite-random integration), before scaling the expensive conditioned estimator to all nine mocks.


### E54 preregistration — 48k random-replica floor

After E53 Stage 2 showed stable F+/F- geometry but material 4800R -> 48000R background reduction, E54 was frozen before execution. It keeps mock0001 full galaxies, E51 tag fields and E19 injection fixed and runs four fresh deterministic 48k random replicas per cap/tracer from the same eligible random pools. Engineering gate: median calibrated-difference replica SD <=0.01, max field SD <=0.02, angle range <=0.002 rad and F-state difference-ratio range <=0.01. Local script SHA256 658e91ad7d3434385844f149762d88f336566fd0928dfe248236b7d9b33b7f87; package SHA256 fa091a20df15989fc6ac7efb10f2029f58faaeeb1843802a539358a88429df73. No result exists yet; observed rows/odd remain sealed.

Canonical catch-up document: `docs/EBOSS_DR16_A03_E30_E53_CATCHUP_AND_CONDITIONED_ROUTE_2026-09-30.md`.


## 30. 9. 2026. — E30–E53 catch-up, unconditioned no-go and conditioned-estimator pivot

**Repository checkpoint before catch-up:** `c140bef3e16508246fb0fec0adb4e56fb93b477c`. Work performed locally/offline after E29 is now synchronized to the audit branch. Detailed checkpoint: [EBOSS_DR16_A03_E30_E53_CATCHUP_AND_CONDITIONED_ROUTE_2026-09-30.md](EBOSS_DR16_A03_E30_E53_CATCHUP_AND_CONDITIONED_ROUTE_2026-09-30.md). Machine-readable summary: `source_data/eboss_dr16_a03_e30_e53_conditioned_route_catchup_summary_2026-09-30.json`. Reproducibility artifacts, including scripts, protocols, local reports/results and the explicitly archived buggy E47 v0, are under `artifacts/e30_e53_catchup/` and the synchronized individual paths.

### E29–E44 physical-closure line

- Local E29 found 20 bounded PID-supported two-epoch candidates, but **zero globally certified progenitor links**. Median late PID overlap is 0.9091, median relative L1 particle-count change 3.091%, median comoving centre displacement 0.2311 h^-1 Mpc and median finite-difference centre-speed proxy 547.5 km/s. These are descriptive proxies, not M200c histories, instantaneous winds or certified merger-tree trajectories. A later targeted profile read crashed Ubuntu and was not repeated.
- E30 sharpened the sign-conditioned estimand: marginal P(S=+)=P(S=-) is not sufficient by itself; the relevant sufficient local condition is E[S|X]=0 for the response-relevant state X.
- E32–E36 established radial-profile, mass-definition, causal-history and hidden kinetic-state non-identifiabilities. E37–E43 converted these into explicit closure/error budgets. Frozen E28R1 low-mass alpha=.4 has |F−−F+|/|FD|=0.6242%, while 94.22589% of the finite-band FD acceleration comes from k>1/Mpc. None of the incoming-state/history/profile/Born/tracer-response terms is bounded tightly enough for a physical sign certificate.
- E44 verified the cleaned/tree Abacus direction as the right next multi-epoch object but found no required local cleaned/tree product. Broad ASDF mining was stopped. **A03 remains PHYSICAL_UNCERTIFIED.**

### E45–E48 original unconditioned 2-point route

- E45 archive-only synthetic odd injection: |A|=.02 is small relative to nine-mock apparent-amplitude scatter, SD 0.09495 NGC and 0.06424 SGC. Exact +/- recovery is an estimator-algebra check, not a physical forecast.
- E46 defined the valid E45 comparison A_phys=(q^T y_phys)/(q^T q) and explicitly rejected direct E28 acceleration -> E45 amplitude identification.
- Corrected E47R1 reproduces the frozen E21 source to machine precision. The physical response is nevertheless tiny: median alpha_g0 7.51e-8 NGC / 7.05e-8 SGC; ~83% descriptive recovery needs |g0| about 5.55e5 / 6.65e5.
- E48 shows the required response is outside the linear regime by orders of magnitude. For diagnostic b_L=2, b_E=1, f=1, the most favorable anti-aligned response needs D=max(|d_L|,|d_E|)=1.33e5 NGC / 1.58e5 SGC, yielding one-sigma local perturbations 60.5 / 71.4.
- **Decision:** retain E21/E47 as a negative/limitation result and do not optimize the unconditioned channel further.

### E49–E53 conditioned route

- E49 direct CLASS vTk passes the strong preregistered tag gate: min |r(v_nu-v_cdm,v_b)|=0.9217 across FD/F+/F−, with minimum Gaussian sign-correlation factor 0.7463.
- E50 separates amplitude-calibrated and shape-only regimes. The conditional F−/F+ vector difference is ~43.6% of the F+ norm, while the amplitude-free angle is only 1.366 deg. With illustrative reconstruction correlation .8, amplitude-calibrated 3-sigma discrimination needs underlying conditioned S/N of order 13, shape-only order 238.
- E51 implements the first pair-midpoint externally marked cross-Landy-Szalay estimator on all nine exact eBOSS EZmocks using four fixed synthetic technical sign fields. All 18 cap/mock cases pass parent sample/pair replay. Median post-window F+/F− angle remains 0.02179 rad NGC and 0.02345 rad SGC. No observed galaxies, physical velocity reconstruction or observed odd vector were used.
- E52 freezes one-dimensional theory directions without a 12D inverse covariance. Descriptive 83% technical-lambda thresholds for the amplitude-calibrated F−−F+ direction are NGC 0.499/0.414/0.381/0.335 and SGC 0.094/0.093/0.276/0.211 for fields A–D. Shape-only remains much harder.
- E53 Stage 1 changes **only sample size** for mock0001: full eligible galaxies + frozen 4800 randoms. Parent replay passes; fingerprint remains stable, angle 0.02499 NGC / 0.02418 SGC and separation ~43.4%. Baseline 12D norms fall strongly relative to E51 600D/1200R, indicating that much of the E51 background was finite-sample noise. E53 remains one realization only and is not an inference result.

### Frozen status and next action

Observed eBOSS galaxy rows and the original 24D odd vector remain **SEALED**. A03 remains **PHYSICAL_UNCERTIFIED** as an absolute observed Einstein–Vlasov prediction; A04 remains **BLOCKED** for robust inference. Nine mocks are not sufficient for an unrestricted 12D/24D covariance inverse. No detection significance or p-value is authorized.

**Next:** E53 Stage 2 on mock0001 with the already-defined nested 48k randoms. If full-D/4800R -> full-D/48000R changes are small, expand the exact same frozen conditioned estimator to all nine mocks. Only after that build a real external velocity reconstruction/tag model, amplitude calibration, coverage/covariance plan and final preregistration. Observed odd remains sealed.

## 29. 9. 2026. — E29 two-time halo profile / Born validity gate

- [E29 fizikalna bilješka](EBOSS_DR16_A03_E29_HALO_HISTORY_TWO_TIME_PROFILE_AND_BORN_VALIDITY_GATE_2026-09-29.md), commit `9091a6c3080bddc982f3d4f9cfb085acd3eda2b6`: u uvjetnom E28 izvornom modelu fiksni komovirajući NFW daje `u(k)^2` izvan vremenskog integrala. Za stvarni evoluirajući halo potrebno je `u_test(k,z_obs) × ∫dz' [M(z')/M_obs] u_src(k,z') K_F(k,z')`; izvorni izraz vraća se samo pod E28 pretpostavkom identičnog oblika na svim vremenima. To je fizikalna strukturna korekcija, **nije** izračun novog halo drag-a ni validacija UV-a.
- Drugi red Vlasovljeva odziva uključuje konvoluciju Fourierovih modova i gradijente po impulsu. Ne može se opravdano aproksimirati množenjem E27 jednog-k kernela proizvoljnim faktorom. Nužni ulazi: konzistentna same-object masa/profil/centar kroz vrijeme, environment/wind, raniji wake ili njegov kontrolirani bound, faznoprostorna kontrola pogreške i neovisni fizikalni UV kriterij. L1 `N` nije `M200c(a)`; ne interpolirati medijane kao povijest objekta.
- Trenutačni audit HEAD nakon E29 bilješke bio je `9091a6c...`; PR #1 ostao draft i unmerged, main netaknut. Na tom commitu CI još nije bio potvrđen pri prvoj provjeri. Ovaj history update ne mijenja fizikalne izvore, zamrznute F, opažene podatke ni ranije numeričke rezultate. **A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED, observed odd SEALED.** Nema nove lokalne WSL naredbe ni downloada bez konkretnog physical data gatea.

## 29. 9. 2026. — E28/E28R1, master plan v1.62

**Namjena:** trajni projektni checkpoint. Kanonski plan je [EINSTEIN_VLASOV_NP_MASTER_PLAN.md](../EINSTEIN_VLASOV_NP_MASTER_PLAN.md), a ovaj zapis olakšava nastavak bez ponavljanja ranijih koraka. Polazni verificirani audit HEAD prije ovoga history zapisa: `de2a15b2360f9738d56d91a8ca00ec7b9a49407b`. Aktualni HEAD i plan uvijek ponovno pročitati prije novih izmjena.

### Repozitorij i nepromjenjiva pravila

- Repo: `dvlahek/stress-energy-closure`; branch `audit/eboss-elg-bit8-ra-orientation-20260925`; draft [PR #1](https://github.com/dvlahek/stress-energy-closure/pull/1). `main` je `e50549a5c974772eb9576d301df7a5ea6e6d87d2` pri checkpointu i ne smije se mijenjati niti PR mergeati.
- Originalni F0/F+/F− iz E8, 4000q i fizikalna izvorna normalizacija, E16 576 trokuta, E4 originalnih devet mockova, izvorni E7 48k, postojeći cuts, seedovi i LRG→ELG midpoint-LOS orijentacija ostaju zamrznuti. Nema prilagodbe teorije opaženome odd rezultatu.
- Originalni opaženi eBOSS LRG×ELG 24D odd vektor i observed galaxy rows ostaju **SEALED**. Izvorni E0/E1 high-z estimand je `0.9≤z<1.0` (NGC/SGC zasebno), pa je E8/E27/E28 z=.95 matematički čvor, **ne** izmjereni pair-weighted `z_eff`. Objavljeni broad-sample `z_eff≈.77` nije naš estimand.
- A-03 fizički kalibrirani predwindow `ξ_{1,3}` i stvarni signed pair-z/LOS prozor **PHYSICAL_UNCERTIFIED**; A-04 nezavisna inferencijska kovarijanca **BLOCKED**. Originalnih 9 mockova daje za 24D najviše rang 8. Nema S/N, detekcijske tvrdnje ili otvaranja opaženih podataka.
- Nema novog CLASS/Abacus/ASDF/FITS/mock/WSL downloada samo da bi se ponovio izvorni audit. Za stvarni lokalni zadatak dati jednu potpunu kopirljivu Ubuntu/WSL naredbu te točne izlazne putanje; ne tražiti korisnika već poznate odgovore. Ne kontaktirati autore.

### Što je već gotovo prije stvarne numerike

E20–E26 odvojili su restriktivni linearni neutrinski LOS odziv od kvadratnog wake-mixed-bispectrum kanala. U E21 vrijedi `Im P_LE=μ(A_L c_E−A_E c_L)C_v`, ali fizički `c_L,c_E` nisu određeni. Standardni Kaiserov gradijent brzine sam daje parni μ², ne automatski odd. E24 je u ograničenom dvopoljnom modelu dokazao egzaktno neodređivanje globalnoga predznaka `c_L,c_E` iz samih parnih auto/cross spektara. E25 zahtijeva ukupni metric `Ψ_total=Ψ_bg+Ψ_nuwake`: istu neutrinsku gravitacijsku silu ne smije se brojiti ponovno kao neovisni Eulerov ostatak. E26 je izvornim E12/E13/E21 potvrdio rang jedan linearnoga jednog adijabatskog density–wind moda (3 stanja × 4 K), ali to ne daje galaktički halo-odziv niti numerički standard-Doppler nuisance rank. Ne vraćati se na iste source-only provjere.

### NOVO: E27/E27R1 — prvi uvjetni fizičko-vremenski Bornov neutrinski wake

[Znanstvena bilješka E27](EBOSS_DR16_A03_E27_CONDITIONAL_PHYSICAL_TIME_BORN_WAKE_AND_FPM_PRECISION_2026-09-29.md), [originalni CI 36599144647](https://github.com/dvlahek/stress-energy-closure/actions/runs/36599144647) i [posthoc-prospektivni numerički R1 CI 36599391672](https://github.com/dvlahek/stress-energy-closure/actions/runs/36599391672) su SUCCESS. Koriste zamrznute 4000q F0/F+/F−, originalni E9 CLASS `H_F(z)`, E17D0 linearni Vlasovljev impuls, dvije **nekalibrirane** Wechsler-like E17D2B1 `α=.4/.8` povijesti pri dva uvjetna masena sidra, `z=1→.95`, početni raniji wake nula, `v_h=+200 km/s`, jedan `k=.05 h_CLASS/Mpc`. R1 je razlučio samo unutar toga benchmarka `(Im_F−−Im_F+)/Im_FD≈+1.024×10^−5` za nižu masu. To nije ukupna sila, halo velocity bias, galaktički observable niti fizička F± granica.

### NAJNOVIJE: E28/E28R1 — stvarna numerika 3D gravitacijskoga recoil-a i UV STOP

[Znanstvena bilješka E28](EBOSS_DR16_A03_E28_CONDITIONAL_3D_HALO_RECOIL_AND_UV_DOMINANCE_2026-09-29.md); [trajni izvorno/CI-grounded E28/R1 sažetak](../source_data/eboss_dr16_a03_e28_conditional_3d_recoil_and_posthoc_r1_summary_2026-09-29.json). [Prvi originalni E28 CI 36601783674](https://github.com/dvlahek/stress-energy-closure/actions/runs/36601783674) SUCCESS za finite band, ali njegov signed F± kontrast **nije** prošao 1 % zahtjev preciznosti. Nakon rezultata zasebno je registriran E28R1 isti model, finija log-k/z mreža: [E28R1 CI 36602291554](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602291554) SUCCESS, sva četiri kontrasta imaju relativni numerički jaz <1 % **same potpisane razlike**.

Model: E27 finite-time Born neutrinski wake; izvorni E17D2B0 median DM14 truncirani NFW `M,c,r_s` s **nekalibriranim fiksnim komovirajućim oblikom** kroz z=1→.95; E17D2B1 uvjetne `α=.4/.8` mase; jedna izvorna masivna neutrinska vrsta mν=.06 eV i originalne F raspodjele; k pojas `[.01,8] Mpc^−1` podijeljen na .01/.1/1/8, originalna ista halo trajektorija +200 km/s, bez ranijega wakea, samogravitacije i shared environmenta. To nije isti kozmo-kalibrirani stvarni LRG/ELG halo.

R1 za **niže maseno sidro, α=.4**: FD `−0.00201850894575251`, F+ `−0.00201220917438763`, F− `−0.00202480872683428` u `(km/s)^2/Mpc`. Potpisani `(a_F−−a_F+)/a_FD` preko četiri uvjetna modela iznosi `+0.00624201,+0.00615897,−0.00344891,−0.00347544`: **predznak ovisi o uvjetnom masenom/profile sidru**.

**Najvažnija fizička granica:** `94.22589 %` FD reference za nižu masu i α=.4 dolazi iz posljednjega k pojasa `1<k≤8 Mpc^−1`. Ovo je numerički razlučen FINITE-BAND integral, ali **nije dokaz UV konvergencije ukupnoga drag-a** niti validnosti linearnog Borna ili trunciranoga unutarnjeg halo profila na tim skalama. Ne mijenjati slobodno `k_max` radi dobivanja pojačanoga signala; niti jedan F± postotak nije physical eBOSS detection/upper bound.

Pri početnom auditu HEAD `de2a15b...`, na istom commitu [repo integrity 36602934708](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602934708), [E27 36602928891](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602928891), [E28 36602928802](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602928802) i [E28R1 36602928919](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602928919) svi su SUCCESS, kao i relevantni E18–E26 i B8/B9 source-only testovi.

### Stvarni sljedeći znanstveni zadatak

1. **Ne raditi E29 kao još jedan sintetički source-only ili proizvoljni `k_max` sweep.** Najprije tražiti neovisno opravdan fizički high-k validity/regularization model, uz kvantifikaciju nelinearnog Born breakdown-a i konačnoga/prostorno usrednjenog halo responsea. Ako se to iz postojećih izvora ne može dokazati, iskreno ostaviti ukupnu silu kao uncertified i razlikovati band-limited benchmark od stvarnoga rezultata.
2. Potrebna je stvarna ili nezavisno opravdana *same-object/same-epoch* host masa/profil/formation history i zajednički halo–environment wake u jednoj konzistentnoj kozmologiji. Postojeći korisnikov Abacus raw superslab-000 lokalni `N` predstavlja assigned L1 particle count, nije `M200c` ni merger tree; službeni standardni z≈.953 sekundarni proizvodi ne daju potvrđeni same-epoch potpuni `M200c(a)`. Ne otvarati nove ASDF stupce ili preuzimati velike podatke bez zasebnoga fizikalno opravdanog data gatea.
3. Tek nakon gornjega izračunati fizički halo/galaxy velocity transfer `β_a` i eventualni stvarno zasebni total-metric Euler residual `ε_a` u skladu s E25, high-z LRG/ELG HOD, Doppler/magnification/evolution i wind selection. Iz njih dobiti `χ_F=A_L c_E−A_E c_L`, puni predwindow `ξ_1,ξ_3(s,z)` i originalni signed pair-z/LOS survey prozor. E20 quadratic mixed B/LOS je zasebna, još nezatvorena ruta.
4. Samo nakon nezavisno certificiranoga fizičkog predloška i A04 kovarijance razmotriti novi predeclared empirijski gate. Observed 24D odd ostaje SEALED.

**Za nastavak:** pročitati zadnji GitHub master plan (v1.62 na ovom checkpointu), E27/E28 bilješke, E28 original/R1 JSON i aktualni PR HEAD; nastaviti s fizičkim UV/halo-history problemom. Ne ponavljati E8–E26 i ne tražiti ponovno već donesene odluke. Bez korisničkog WSL zadatka u ovoj history fazi.
