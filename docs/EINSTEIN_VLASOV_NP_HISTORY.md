# EinsteinVlasovNP — povijest i handoff za novi chat

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
