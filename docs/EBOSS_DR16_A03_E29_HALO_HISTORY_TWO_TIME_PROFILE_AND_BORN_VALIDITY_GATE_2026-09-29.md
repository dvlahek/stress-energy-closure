# E29 — two-time halo-profile consistency and Born-validity gate (29 September 2026)

**Status:** analytic physical-model gate, not a new numerical detection, UV convergence proof, or calibrated halo force. Parent audit HEAD before this note: `0b926bd2607cc1ab2210b5b57fe4702cda156e41`. Main and original E8/E16/E27/E28 source files remain unchanged.

## 1. Concrete limitation found in the E28 integrand

The original E28 runner `scripts/audit_eboss_dr16_a03_e28_conditional_3d_finite_band_neutrino_recoil.py` computes the frozen finite-band Born recoil with `u(k)**2` outside its redshift integral. Its registered protocol explicitly fixes the **comoving** NFW shape throughout z=1 to z=.95, while multiplying the source by the E17D2B1 conditional mass history. This is internally consistent **only for that declared shape-history ansatz**. It does not establish a physical same-object evolving halo. In particular, the population DM14 median c(M,z=.95) and two illustrative alpha histories do not provide a measured c(a), r_s(a), truncation history or environmental potential for any LRG/ELG host.

For a time-dependent spherical halo, write its normalized source profile as `u_src(k,z')` and the observed recoil mass-weighting profile as `u_test(k,z_obs)`. In the same externally forced linear Born model, the E28 factorization

`u(k)^2 * integral_{z_obs}^{z_init} dz' [M(z')/M_obs] K_F(k,z')`

must be replaced by

`u_test(k,z_obs) * integral_{z_obs}^{z_init} dz' [M(z')/M_obs] u_src(k,z') K_F(k,z')`.

Here `K_F` denotes exactly the original E28 free-streaming, H(z), angular wind and time kernel with all original E8/E9 F-state inputs and units preserved. The expression is a **structural identity within the external-Born approximation**, not a calibrated prediction. The original E28 expression is recovered if `u_src(k,z')=u_test(k,z_obs)=u(k)` at all source times.

An evolving profile must be derived from a dynamically compatible source. Simply inserting independently fitted NFW snapshots or interpolating population medians does not enforce the continuity equation, accretion flux, halo-centre trajectory, mass definition or environment. The changing total M200c can also include pseudo-evolution due to changing reference density. No real halo formation history is inferred from the existing Abacus L1 particle count N.

## 2. What a Born-validity calculation must actually test

The linearized Vlasov operator `D_0` and gravitational forcing `L_Psi` give
`D_0 f1 = L_Psi F`, `D_0 f2 = L_Psi f1`.
The second-order source couples spatial Fourier modes and phase-space momentum derivatives. It is not obtained by multiplying the existing one-k E27 kernel by a scalar correction. For a meaningful comparison, the same specified halo-plus-environment potential, history, initial distribution and initial-time wake must enter f1 and f2. The appropriate target is the **halo-averaged signed acceleration** and its scale-resolved second-order correction, not only a density contrast or a single Fourier mode.

The z=1 zero-wake initial condition in E27/E28 is a benchmark boundary condition, not evidence that the physical earlier wake vanishes. The background total metric and neutrino wake must remain consistently accounted for under E25; do not add the same wake acceleration again as an independent Euler residual.

## 3. Prospective physical/data gate before any E29 numerical execution

Required: (i) a documented same-object, same-cosmology halo mass/profile/centre trajectory or an independently justified dynamical halo model; (ii) a specified treatment of halo-environment shared wake and relative wind; (iii) a justified initial neutrino response before z=1 or a controlled bound on the omitted earlier response; (iv) a phase-space numerical method resolving the nonlinear Vlasov correction and its force with independently stated error controls; (v) a physical criterion for short-scale validity or regularization fixed independently of F+ versus F− contrast and observed galaxy odd data.

No arbitrary k_max sweep, NFW core insertion, hand-chosen softening, population-median interpolation claimed as same-object history, or synthetic source-only CI can satisfy these requirements. Numerical convergence at fixed finite k band cannot certify UV convergence of the full force. If the required inputs are unavailable, leave the physical total halo acceleration **UNCERTIFIED** and do not calculate beta, epsilon, chi_F, galaxy xi1/xi3, significance or an eBOSS upper bound from E28.

## 4. Explicit boundaries

Original E8 4000q F0/F+/F−, E16 576 triangles, E7 48k, E4 nine mocks, source identities, cuts, seeds and signed LRG→ELG midpoint LOS remain frozen. Observed eBOSS 24D odd and observed galaxy rows remain SEALED. A03 PHYSICAL_UNCERTIFIED; A04 BLOCKED. No new CLASS/ASDF/Abacus/FITS/mock download or WSL execution is authorized by this note. No changes to main or merge of draft PR #1.
