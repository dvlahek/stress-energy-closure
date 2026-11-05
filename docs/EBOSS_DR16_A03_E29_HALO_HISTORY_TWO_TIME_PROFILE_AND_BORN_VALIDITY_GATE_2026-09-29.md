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

## 5. Source-verified Abacus radial-percentile feasibility and mass-definition correction (29 September 2026)

The pinned official [AbacusSummit data-products specification](https://github.com/abacusorg/AbacusSummit/blob/4b1959c710cb0c49aa305c6213a228aa2a4587ff/docs/data-products.rst) explicitly defines raw `r100_L2com` as the radius enclosing 100% of **L2 halo mass relative to its L2 center**, and `r10_L2com`, `r25_L2com`, `r33_L2com`, `r50_L2com`, `r67_L2com`, `r75_L2com`, `r90_L2com`, `r95_L2com`, `r98_L2com` as *int16 ratios to r100*, condensed to [0,30000]. These are L2 catalog mass percentiles, **not** spherical-overdensity `M200c` or an independently reconstructed total-matter `M200c(r)` profile. Raw radii have unit-box conventions and must be decoded according to the official loader/specification before comparing physical/comoving lengths. Do not interpret a raw int16 as a length or equate L2 mass percentiles with the L1 host's 200critical enclosed profile.

The same official source documents secondary approximate `z=0.95,1.025` halo catalogs and subsample PIDs (no corresponding secondary halo/field RV), with actual redshifts to be read from the headers. PID association makes a *prospective* same-object catalog-shape test possible, but no tree linkage or second-epoch radial values have yet been measured for the user's fixed raw halo IDs. The pinned [data-access specification](https://github.com/abacusorg/AbacusSummit/blob/4b1959c710cb0c49aa305c6213a228aa2a4587ff/docs/data-access.rst) notes that merger trees may require browsing the Globus directory rather than the portal web interface; this is **not** confirmation of a specific c000/ph000 file path, download availability or matching halo IDs.

**Prospective narrowly defined test:** only after a verifiable same-simulation tree/PID linkage and documented radial-column decoder, compare for the same halo the dimensionless `r10/r100`, `r50/r100`, `r90/r100` at two actual header redshifts. Record linkage quality and the raw N and tree's documented matching criteria; the user's previously fixed B7 N=44/53/61 QA rows are not presumed eligible. A detected change would reject the fixed-*catalog-L2-percentile* shape assumption for those matched objects only. No change would not validate NFW, M200c(a), the environment, the short-scale Born approximation or the full halo drag. No new halo files were downloaded or read for this source-only correction; observed odd remains SEALED, A03/A04 unchanged, main untouched and PR draft.
