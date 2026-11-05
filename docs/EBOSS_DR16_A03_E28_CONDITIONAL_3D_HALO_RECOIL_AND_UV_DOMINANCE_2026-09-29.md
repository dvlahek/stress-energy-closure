# E28 — 3D halo recoil od neutrinskoga wakea i UV-dominacija

**29. 9. 2026.** Nova fizički dimenzionirana numerika integrira originalni E27 neutrinski wake po kutovima i po konačnom pojasu komovirajućih valnih brojeva te računa gravitacijsko ubrzanje po jedinici mase haloa. Nije puni ukupni halo drag, kalibrirani LRG/ELG odziv, puni Einstein–Vlasov ni eBOSS detekcija.

## Model i jednadžba

Izvorni [E27](EBOSS_DR16_A03_E27_CONDITIONAL_PHYSICAL_TIME_BORN_WAKE_AND_FPM_PRECISION_2026-09-29.md) već daje Bornov odziv originalnih 4000q F0/F+/F− po stvarnom konformalnom vremenu, ali uz uvjetne Wechsler-like halo mase, z=1→.95, nulti raniji wake i konstantni v_h=+200 km/s. Ovdje preuzimamo ista dva originalna sidra mase (1e12 i 1e13 h_DM14^-1 Msun), α=.4/.8 te originalne E17D2B0 median DM14 NFW koncentracije 5.43987438 i 4.56340350. NFW oblik u komovirajućim koordinatama i koncentracija ostaju fiksni na ovom kratkom intervalu, normalizacija prati M_α(z). To je dodatna, nekalibrirana pretpostavka profila. DM14 halo cosmology i originalni CLASS H(z) nisu identični.

Srednja fizička rest-mass gustoća jednog masivnog neutrino+antineutrino stanja izravno se računa iz originalnog CLASS-normaliziranog F(q) sa mν=.06 eV. Originalni FD pri z=.95 daje 1.32559455834e9 Msun/Mpc³. Kinetička korekcija pozadinskoj gustoći i puni relativistički metric closure nisu uključeni.

Uz Fourierov neutrinski wake δν(k) po halou, fizičko halo-usrednjeno ubrzanje jest

\[
a_\parallel^{w}=4\pi G a\bar\rho_\nu\int\frac{d^3k}{(2\pi)^3}
\frac{ik_\parallel}{k^2}\delta_\nu^w(\boldsymbol k)\,u_h(k).
\]

Konačni NFW u_h(k) pojavljuje se jednom u izvoru E27 i jednom u halo-usrednjenju sile. Kutna integracija je egzaktna: integral od −1 do 1 od iμ exp(i k μ Δx) iznosi 2 sinc'(k Δx). Novi numerički pojas iznosi k=(.01,.1,1,8) Mpc^-1; to je samo numerički k-band, NIJE novi galaxy science cut. Konačan rezultat nosi jedinice (km/s)²/Mpc. Kod E25 ovu silu treba nositi kroz isti ukupni neutrinski potencijal; ne dodavati je opet kao neovisan Eulerov ostatak. [Okoli i sur., MNRAS 468 (2017) 2164](https://doi.org/10.1093/mnras/stx560) u punijem halo-modelu također uključuju zajednički two-halo wake koji E28 ne računa.

## Originalni E28 i registrirana posthoc preciznost

[E28 protokol](../source_data/eboss_dr16_a03_e28_conditional_finite_band_3d_wake_recoil_protocol_2026-09-29.json), Git blob e8feedc4d09965a24d8eb60c16454d6e03b2fbff. Njegova prva verzija govorila je o fiksnom fizičkom r_s uz jednadžbu s fiksnim komovirajućim u(k)²; unutarnja nedosljednost ispravljena je na fiksni komovirajući oblik PRIJE prvoga E28 izvođenja. [E28 prvi CI 36601783674](https://github.com/dvlahek/stress-energy-closure/actions/runs/36601783674) SUCCESS za osnovni pojasni force i 3 F × 4 modelna halo slučaja. Originalni niži maseni FD α=.4 dao je −0.002018508621961 (km/s)²/Mpc na k≤8. Izvorni E28 signed Fminus−Fplus kontrasti NISU prošli raniji zahtjev relativnog jaza manjeg od 1 % SAME razlike: gapovi 13.518 %, 13.649 %, 1.7368 %, 1.7101 %.

[E28R1 posthoc protokol](../source_data/eboss_dr16_a03_e28r1_posthoc_signed_recoil_contrast_and_band_convergence_protocol_2026-09-29.json), Git blob 5d0a0ac96364916e95328c0f9ef6d30fb528c229, zaključao je prije dodatnih izračuna isti model i pojas uz mreže 24/48 Simpson intervala po pojasu i 257/513 fizičkih z čvorova. [Prvi R1 run 36602177580](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602177580) fail-closed bio je tehnički: stari originalni E28 helper namjerno odbija nov R1 grid 48. Samo je R1 dobio vlastiti preregistrirani k-grid; originalni E28 protokol, kod i fizika ostali su nepromijenjeni. [Završni R1 CI 36602291554](https://github.com/dvlahek/stress-energy-closure/actions/runs/36602291554) SUCCESS. Puni izvorni JSON pohranjuje Actions, artefakt 11050535270; bajtovi ZIP-a nisu ponovno verificirani pri pisanju ove bilješke.

R1 niži M0≈1.490313e12 Msun α=.4 daje FD −0.00201850894575251, plus −0.00201220917438763, minus −0.00202480872683428 u (km/s)²/Mpc. R1 signed (a_Fminus−a_Fplus)/a_FD za četiri zaključana modela redom: +0.006242009714 (Mlow α=.4), +0.006158969241 (Mlow α=.8), −0.003448908126 (Mhigh α=.4), −0.003475436854 (Mhigh α=.8). Jaz SAME razlike između dviju finih mreža: 0.4195 %, 0.4238 %, 0.08143 %, 0.08014 %. Sva četiri potpisana pojasna kontrasta su NUMERIČKI razlučena unutar tog uvjetnog modela; predznak se čak mijenja s nekalibriranim masenim/profilnim sidrom. To nije robustan predznak stvarne LRG/ELG fizike.

## Znanstveni rezultat koji ograničava sljedeće korake

**94.22589 % FD pojasnoga ubrzanja za niže sidro α=.4 nastaje samo na 1<k≤8 Mpc^-1.** Precizan integral na izabranom pojasu NIJE dokaz konvergencije ukupne sile kad k_max raste. Najveći dio računa dolazi iz skala na kojima su neprovjereni unutarnji profil, povijest, halo environment i primjena linearnoga Borna najvažniji. Ne prenositi ni pojasnu amplitudu ni F± kontrast na fizički kalibrirani total drag ili galaktički odd.

Za fizički nastavak treba same-object M(a), konačni, kozmo-konzistentni halo profil i okoliš, duža i fizički specificirana povijest neutrinskoga wakea, opravdana kratkovalna regularizacija/nelinearna kontrola, shared halo force i high-z LRG/ELG HOD, Doppler i selekcija. Lokalni Abacus L1 N nije M200c ni merger tree. A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED; opaženi 24D odd SEALED. Nema novih CLASS/FITS/ASDF/Abacus/mocks/WSL, cutova/seedova, retuninga E8/F±, kontakta autora ili main mergea; PR #1 ostaje draft.

[Trajni E28/E28R1 CI-grounded sažetak](../source_data/eboss_dr16_a03_e28_conditional_3d_recoil_and_posthoc_r1_summary_2026-09-29.json) bilježi izvorni i naknadni rezultat, uključujući neuspjele međufaze i numeričke granice.
