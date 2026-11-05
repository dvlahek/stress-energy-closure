# A-03E9 — zamrznuti Einstein–Vlasov izvor, tri CLASS evolucije i uvjetna wake faza

**27. 9. 2026. | Status: E9 numerički dovršen, eBOSS galaktička A-03/A-04 fizika i dalje STOP.** Izvorne pozitivne distribucije `F_0,F_+,F_-` iz E8 nisu ponovno optimirane. U tri zasebna CLASS procesa i istom [izvornom CLASS commitu](https://github.com/lesgourg/class_public/commit/e85808324f51fc694d12e3ed7439552a3c3f9540) računat je svaki izvorni `nu-CDM` linearni transfer te vremenski promjenjiva, `R=16 h^{-1} Mpc` filtrirana relativna brzina. Na njoj je izračunana **uvjetna kvazistatička** Eq20/33 faza za isti komovirajući `k=0.05 h/Mpc` i dvije unaprijed zaključane numeričke derivacije.

To je novi numerički korak nakon E8c, a nije puni vremenski-retardirani halo Einstein–Vlasov izračun, izvor izmjerenog eBOSS galaktičkog vjetra, kalibrirani `xi_1/xi_3` ni inferencijska značajnost. Opaženi odd vektor ostao je zatvoren.

## Problem i zamrznuti ulazi

[Okoli et al. (2017), Eq. (32)–(34)](https://doi.org/10.1093/mnras/stx560) povezuju imaginarni dvotracerski RSD signal s vremenskom promjenom wake faze. E8c je odredio samo statičku fazu pri proizvoljno zadanom `v_parallel=200 km/s`. Za njezinu *vremensku derivaciju* trebaju nam zasebni originalni `nu-CDM` transferi i njihove promjene za sva tri već zamrznuta stanja.

**Prije novog CLASS računanja** zaključani su [E9 protokol](../source_data/eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase_prereg_2026-09-27.json), `m_nu=0.06 eV`, `z_match=1100`, 30% deformacija, isti `q=linspace(0,20,4000)`, izvorni [E8 frozen CSV SHA256 `bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0`](../source_data/eboss_dr16_a03_e8_frozen_kinetic_resonance_and_static_SI_green_result_2026-09-27.json) i CLASS `e85808324f51fc694d12e3ed7439552a3c3f9540`. Nisu odabrani novi koeficijenti iz fullmock odstupanja. [E9 implementacija](../scripts/build_eboss_dr16_a03_e9_frozen_class_wind_quasistatic_phase.py) koristi originalni Newtonian-gauge `d_ncdm[0]−d_cdm` te istu source-only fizičku Eq20 occupanciju kao E8c, **ne** neprovjereni objavljeni `mu=1` approximation. Izvorni [DOI i `mu` erratum](EBOSS_DR16_A03_E8C_EXACT_WAKE_MODE_PHASE_AND_EQ20_CONVENTION_ERRATUM_2026-09-27.md) nije poništen ovim novim rezultatom.

## Što CLASS računa, a što je naš uvjetni model vjetra

U originalnoj `code/wake_two_tracer_fisher.py` transfer konvenciji,

\[
v_{\nu c}(k,z)=-\,\frac{c H_{\rm CLASS}(z)}{k_{\rm Mpc}}\,
\frac{\partial (\delta_\nu-\delta_c)}{\partial z}\,
\sqrt{A_s(k_{\rm Mpc}/0.05\,{\rm Mpc}^{-1})^{n_s-1}} .
\]

Ovo je **signed transfer-derived mode amplitude**, ne deterministička brzina konkretnoga haloa. Originalni `KVEL=geomspace(10^{-4},0.15,180)\ h/Mpc`, originalni sferni top-hat `W(k R)` s `R=16 h^{-1}Mpc`, i predeclared kumulativni završetak na `k=0.10 h/Mpc` daju

\[
\sigma_{\nu c,R16}^2(z)=
\int_{10^{-4}}^{0.1}d\ln k\,|v_{\nu c}(k,z)|^2\,W^2(k R).
\]

Za predeclared ilustrativni pozitivan LOS `+1 sigma` koherentni rang koristimo `v_\parallel(z)=+\sigma_{\nu c,R16}(z)/\sqrt 3`. **Isti standardizirani rang** drži se kroz sedam redshift čvorova i sva tri CLASS stanja. Pretpostavka koherentnosti, izotropnosti i cutoff `0.1 h/Mpc` su dio modela, ne izmjereni eBOSS winds. Za svako od tri `F` stanja radi se *vlastita* CLASS evolucija i vlastita `sigma(z)`. Ne odbacivati promjenu pozadinskoga neutrinskog transfera umjetnim držanjem FD brzine, osim kao zasebni eksplicitni kontrolni slučaj.

Kvazistatički fazni koeficijent za **isti komovirajući** `k`, `a=1/(1+z)`, jest

\[
\phi_s(k_{\rm com},z)=
\frac{2Gm_\nu^4}{\hbar^3 k_{\rm com}^2}\,
a^2 v_{\parallel,s}(z)\,
f_s\!\left(\frac{m_\nu |v_{\parallel,s}(z)|a}{cT_{\nu0}}\right).
\]

`f_{\rm FD}(q)=1/(e^q+1)`. Za `F_\pm` vrijedi `f_\pm=f_{\rm FD}\times [F_\pm(q)/F_0(q)]`, uz istu **originalnu E8 interpolaciju na 4000-q gridu**, bez max-normalizacije predloška. Ovo je prvi točni FD-occupancy oblik objavljene Eq20, uz kvazistatičku mapu Eq33. **Nema** zajedničkoga, fizički dokazanoga halo-potential transfera ni vremenski retardiranog Vlasov kernela u ovoj formuli.

Predefinirani z-čvorovi su `[.9,.94,.945,.95,.955,.96,1.0]`. Centrirana derivacija na matematičkom midpointu `z=.95` (NE izmjereni `z_eff`) računa

\[
\frac{d\phi}{d\ln a} =
-(1+z)\frac{\phi(z+h)-\phi(z-h)}{2h},
\qquad h=0.01,\ 0.005 .
\]

`d\phi/dt=H_s(z)d\phi/d\ln a`, uz zasebni CLASS `H_s(z)` za svako stanje. Velocity-transfer `d\delta/dz` dodatno je ponovljen s unaprijed zaključanim relativnim koracima `h_v/(1+z)=0.004` i `0.002`. Za tehničku QA unaprijed je upisano 5% upozorenje, **ne** post-hoc fizička tolerancija estimatoru ili signalu.

## Dovršeni rezultati, originalni `h_v/(1+z)=0.004`, finiji fazni `h_z=0.005`

| `z=.95, k_com=.05 h/Mpc` | FD | `F_+` | `F_-` |
|---|---:|---:|---:|
| `v_parallel=+sigma_R16/sqrt 3` (km/s) | 135.68201584 | 135.86610262 | 135.48939940 |
| Trenutačna uvjetna `phi_s` | 0.00021385370943 | 0.00018071080762 | 0.00024685926277 |
| `dphi_s/dln a` | 0.00044316472341 | 0.00036005025246 | 0.00052587375696 |
| `dphi_s/dt` (s⁻¹) | 1.68203520874e−21 | 1.36657359083e−21 | 1.99595797784e−21 |

Ovdje `phi_plus−phi_minus=-6.614845514749782e−5` i `d(phi_plus−phi_minus)/dln a=-1.6582350450162594e−4` (`d/dt` razlika oko `−6.293843870120874e−22 s^{-1}` uz vlastiti `H_s`). Sljedeće potpuno **algebarsko** dijeljenje čuva isti originalni `F_\pm`:

\[
\left.\frac{d\Delta\phi}{d\ln a}\right|_{\rm total} =
\underbrace{-1.6718070025205302\times10^{-4}}_
{\rm F_\pm\ at\ common\ FD\ wind}
+\underbrace{1.3571957504271057\times10^{-6}}_
{\rm change\ due\ to\ CLASS\ winds}.
\]

Za ovaj specifični `+1 sigma` model drugi član čini `0.00818458` apsolutne ukupne razlike: dominira izravna occupancy razlika, a vlastiti F± linearni relativni vjetar djelomično ju kompenzira. To **ne** znači da je udio 0,82% univerzalna činjenica za stvarne haloe, redshift raspone, bias/HOD ili opažene galaksije.

[CI run `36325438404` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36325438404), arhivirani GitHub Actions [artefakt `10933798674`](https://github.com/dvlahek/stress-energy-closure/actions/runs/36325438404/artifacts/10933798674), puni 26 062-bajtni JSON SHA256 `450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890`. Povezani [git-auditirani sažetak](../source_data/eboss_dr16_a03_e9_CLASS_quasistatic_phase_independently_audited_summary_2026-09-27.json) čuva tri individualna CLASS-state JSON SHA, izvornu verziju, točnu metodu i 14 osnovnih slučajeva. Neovisno ponovno računanje iz user-upload E8 4000-q CSV i sva tri CI JSON-a rekonstruira svih 14 faza i sva četiri derivative stencil izračuna s maksimalnim apsolutnim odstupanjem faze `4.88e−19` i **devet od devet** upozoravajućih engineering QA kontrole bez upozorenja (maksimalni relativni numerički korak `3.73212e−5`). Ovaj test **nije** dodatna nezavisna kozmološka realizacija ili konvergencija dodatnih CLASS ncdm hierarchy parametara.

## Fizički stop i korak koji može proizvesti eBOSS `xi_1`

E9 je stvarna CLASS evolucija originalnog linearnog relativnog transfera pod tri poznata izvorna `F` stanja i uvjetna kvazistatička transformacija Eq20. Međutim, ne rješava punu vremenski-retardiranu Vlasovovu perturbaciju oko realnoga halo potencijala. Ne računamo selekciju halo masa, LRG/ELG response i njihovu evoluciju, Doppler/magnifikaciju/evolution bias, niti realan par-redshift survey-window (NGC/SGC, originalni ELG chunk/depth ovisnosti). E9 ne pokazuje da `k_R16≤0.1` i pozitivan jedan-σ LOS rang predstavljaju realnu brzinu svakoga haloa.

Zato `dphi/dln a` još **nije** eBOSS `xi_1` ili `xi_3`. A-03 empirijski parni prozor je `PHYSICAL_UNCERTIFIED` i A-04 eBOSS inferencijska `C` ostaje `BLOCKED`; devet originalnih mockova ne može invertirati zajednički 24D covariance. Izvorni original v1 SGC independent reverse kvar ostaje otvoren i 48k random apsolutna konvergencija nije fizikalno certificirana. Originalni `main`, originalni parent JSON-i, rezovi i izvori nisu dirani. Opažene galaksije, opaženi randomi i observed odd vektor **nisu čitani**; nema novih mock preuzimanja, seedova ili prilagodbe maski.

**Sljedeći fizikalni posao:** jasno definirani i izvan opaženog odd-a kalibrirani wake-to-halo/tracer response i njegov redshift/pair-window prijenos, uz zasebni null/physical-injection protokol i eBOSS nezavisnu covariance. Bez toga se iz ove faze ne smije tvrditi `>3sigma`, isključenje ni numerički fizički acceptance.
