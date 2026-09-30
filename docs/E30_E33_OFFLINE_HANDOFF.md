# EinsteinVlasovNP — offline E30/E32/E33 handoff (29. 9. 2026.)

## Zaključeni mali rezultati

**E30 / estimand:** prethodni E10 E19–E26 ostaju izvorni. Novi eksplicitni uvjet `m_F=E_q[t_F(X) E_q(S|X)]` razlikuje *conditional* sign balance od samo marginalnoga `P(+)=P(-)`. U središnje simetričnom kontrolnom modelu `E[delta v]=1`, ali `E[v delta delta]=0`: linearni E21 i kvadratni E20 nisu ista tvrdnja. Ne postoji fizički izračunan originalni high-z `chi_F=A_L c_E-A_E c_L`, `eta_F(X)` ili 24D EV mean. QA 21/21.

**E32 / radial u(k):** dva pozitivna sferna toy profila imaju svih 10 istih Abacus L2 percentile radijusa/r100, a različite Fourierove transformacije. Konstruirana je opća nenametnuta ograda `|u(k)-Σ_i w_i sinc(k m_i)|≤|k| Σ_i w_i Δr_i /4`, ograničena na sferni L2 profile. Bez stvarnih dekodiranih percentila nije primjenjena na 20 haloa, a ni s njima ne daje M200c, puni anisotropni ili UV-validated force. QA 41/41.

**E33 / mass definition:** iz E29 stvarnih redshiftova uz flat matter+Lambda približnu pozadinu `rho_c(late)/rho_c(early)=0.915023317906`. Čisto hipotetski nepromjenjivi fizikalni power-law haloi daju prividni M200c rast 4.5403% (gamma2) ili 1.7920% (gamma2.5) samo zbog pomaka overdensity reference. Ovo nije mjerenje 20 haloa ni procjena stvarnog N rasta. QA 12/12.

## Cjelina stvarnog projekta

- Abacus E29: 20/20 jedinstvenih **u testiranim lokalnim prozorima** PID-podržanih dvovremenskih parova; 0 globalno certificiranih merger-tree veza, 0 validiranih M200c(a) ili stvarnih `u_src/u_test(k,z)`. Uzorak je uska x/y traka superslaba 000, ne nasumična populacija. Ne otvarati novi ASDF profil nakon WSL rušenja samo radi percentila.
- E27/E28R1: pravi uvjetni finite-time Bornov rezultat i numerički razlučen finite-band signed halo recoil, ali 94.22589% FD reference niže mase alpha.4 iz k>1/Mpc, **UV/physical total drag uncertified**. Odvojiti E25 total-metric force i ne duplicirati wake u Eulerovu ostatku.
- A03: fizički EV 24D predwindow i wind-conditioned pair-window **PHYSICAL_UNCERTIFIED**, originalni observed eBOSS galaxy rows/odd SEALED; A04: eBOSS 24D joint covariance BLOCKED, devet originalnih mockova daje max rank osam; DESI 18D/120-mock rezultat ne prenositi.
- GitHub `dvlahek/stress-energy-closure`, samo `audit/eboss-elg-bit8-ra-orientation-20260925`, draft PR #1. Ovaj paket nije ništa pisao na GitHub. `main` netaknut. F0/F+/F- i originalni 4000q, E16, E7 48k, cuts, seeding, originalna orijentacija ostaju zamrznuti; nema kontakta autora.

## Sljedeći fizički prioriteti bez lovljenja rezultata

1. **Fizički linearni transfer**: za originalni `[.9,1)` LRG i ELG iz zasebno dokumentiranih halo–environment povijesti, neutrinske gravitacije u ukupnom potencijalu i opažajne number-count projekcije izvesti `beta_a,F`, mogući zasebni `epsilon_a,F`, standardni `B_a` i stvarni wind-selection `s^w_a,F` → time `chi_F(k,z)`. Ako model opravdano daje `chi=0`, zabilježiti no-go uvjet. Bez toga E21 `C_v` nije galaktički EV template.
2. **Alternativni conditional/3pt transfer**: samo uz joint sign-conditioned `p_F(S,X)` i objavljeni tracer-specifični empirijski prozor razmotriti nenulti unweighted intrinsic mean; ili zasebno izvesti fizikalni mixed `B_vδδ`, `Gamma` i LOS-odd projekciju. Velocity-conditioned observable je NOVI estimator, ne retroaktivna oznaka staroga 24D.
3. **UV/source dynamics**: prije ponavljanja E28 k-integrala definirati nezavisan high-k Born-validity/physical regularization kriterij i same-object consistent mass/profile/environment history. Percentilne ograde E32 mogu samo kvantificirati dio profilne nesigurnosti ako stvarni podaci budu dostupni bez rušenja.
4. **A03/A04 tek nakon 1–3**: fizički `xi_1,xi_3` sa standardnim GR/nuisance kroz high-z signed LRG×ELG pair window, neovisna mock kovarijanca/coverage u 24D i pisano odobrenje prije ijednoga observed odd pristupa.

## Reprodukcija

Iz ovoga paketa u običnom Pythonu 3 (bez vanjskih paketa, WSL-a ili mreže) pokrenuti `python e30_offline_estimand_qa.py`, `python e32_offline_radial_quantile_qa.py`, `python e33_offline_mass_definition_qa.py`. E33 očekuje mali parent `e29_20_pair_posthoc_offline_kinematic_audit.json`, uključen u paket. Postojeći JSON izlazi se **ne prepisuju**: replay zahtijeva identične sadržaje. Kontrolni ispis je `E30_E33_OFFLINE_REPLAY.log`; potpuna SHA256 lista `E30_E33_SHA256_MANIFEST.json`.
