# A-03E17D2b0 — nezavisno kalibrirana koncentracija populacije haloa, fizički radijus i granica povijesti nastanka

**28. 9. 2026. · Originalni i neovisni uvjetni prostorni halo model PASS. Stvarna individualna halo povijest, fizički kombiniran CLASS–halo izvor i puni galaktički bispektar BLOCKED. Opaženi odd SEALED.**

## Problem

E17D2a je izveo prostorno konzistentan, pozitivni NFW profil odrezan na r200, s konačnom masom i Newtonovim potencijalom koji zadovoljava Poissonovu jednadžbu. Tada su koncentracije c=(2,4,8) i numeričke skale α=(.1,1) služile **isključivo numeričkoj provjeri**. Za interpretaciju izvornoga E16 valnog područja kao stvarnog haloa treba izvana kalibrirati koncentraciju, masu i radijus uz preciznu definiciju masene konvencije.

E17D2b0 dodaje objavljenu, numeričkim simulacijama kalibriranu *populacijsku* c200c(M200c,z) relaciju, bez čitanja opaženoga odd vektora ili naknadnoga podešavanja na eBOSS. Odvojeno provjerava pripadnu 200ρcrit masu/radijus i statičku prostornu NFW težinu na originalna oba E16 kratka kraka i dugi K. Ovo nije povijest konkretne galaksije niti istovremena neutrino–halo simulacija.

## Vanjska znanstvena kalibracija i ograničenja njezine primjene

[Dutton i Macciò, MNRAS 441 (2014) 3359, DOI 10.1093/mnras/stu742](https://academic.oup.com/mnras/article/441/4/3359/1209689), objavili su koncentracijsko-masenu relaciju za NFW opise relaksiranih *dark-matter-only* haloa u Planck-2013 kozmologiji. E17D2b0 koristi točno njihov **c200c/M200c**, a ne različitu cvir/Mvir ili Einasto relaciju:

\[
\log_{10} c_{200c}(M,z)
=a(z)+b(z)\log_{10}
\left(\frac{M_{200c}}{10^{12}\,h_{\rm DM14}^{-1}M_\odot}\right),
\]

\[
a(z)=0.520+0.385 e^{-0.617z^{1.21}},
\qquad
b(z)=-0.101+0.026z.
\]

Iz izvorne tablice njihove simulacijske kozmologije u ovom *uvjetnom halo modelu* uzeti su Ωm=0.3175, hDM14=0.671, σ8=0.8344, ns=0.9624 i ravna pozadina ΩΛ=0.6825. Objavljena relacija vrijedi kao statistički opis populacije u kalibriranom masenom i redshift rasponu; ne daje poznatu koncentraciju, masu ili individualni merger tree eBOSS LRG/ELG galaksija. U članku je opisan intrinzični rasap koncentracije približno 0.11 dex oko rezultata na z=0; E17D2b0 koristi ±0.11 dex na z≈.95 **samo kao unaprijed zaključanu ilustrativnu osjetljivost**, ne kao dokazanu z≈.95 1σ distribuciju ili HOD prior. Autori također navode da Einasto bolje opisuje neke radijalne profile, pa je sama NFW forma zaseban modelni sustavak.

[Preregistrirani E17D2b0 protokol](../source_data/eboss_dr16_a03_e17d2b0_dm14_population_halo_mass_concentration_prereg_2026-09-28.json), Git blob 8faaf167403e09f321a6c62f452120f8b96898a2, commit 8313645, nastao je **prije novoga numeričkog E17D2b0 testa**. SHA-zaključava originalni E8/E16, E17A, E17D1, E17D2a originalni i nezavisni certifikat te manifeste/protokole. F0/F+/F−, E14/E15 originalnih 72 source i 24 contrast, prvotnih 576 E16 trokuta i LRG→ELG redoslijed ostaju neizmijenjeni. Nema novih opaženih kataloga, mockova, CLASS pokretanja, cutova, seedova ili unblindanja.

## Fizička masena konvencija i važna kozmološka razlika

Uvjetna DM14 pozadina, samo za izračun pojedinačnih halo **snapshotova**, definira

\[
H_{\rm DM14}(z)=100h_{\rm DM14}
\sqrt{\Omega_{m,\rm DM14}(1+z)^3+\Omega_{\Lambda,\rm DM14}}
\quad{\rm km\,s^{-1}\,Mpc^{-1}},
\]

\[
\rho_{\rm crit,DM14}(z)=\frac{3H_{\rm DM14}(z)^2}{8\pi G},
\qquad
r_{200c,\rm phys}=
\left(\frac{3M_{200c,\rm phys}}{4\pi\,200\rho_{\rm crit,DM14}(z)}\right)^{1/3},
\]

\[
r_{s,\rm phys}=r_{200c,\rm phys}/c_{200c},
\qquad
r_{s,\rm com}=(1+z)r_{s,\rm phys}.
\]

G=4.30091×10⁻⁹ Mpc (km/s)²/M☉ je izrijekom zaključana fizička konstanta u računu. Masa koja se prikazuje kao hDM14⁻¹ M☉ pravilno se pretvara u fizičke M☉. Ova uvjetna pozadina koristi plosnati matter+Λ H(z) izraz i ne uključuje puni CLASS relativistički prijenos.

**Izvorni CLASS u našem projektu koristi hCLASS=0.6736, dok DM14 koristi hDM14=0.671.** Razlika omjera hCLASS/hDM14−1 = 0.003874813710879277. E17D2b0 tu razliku ne skriva: za svaki originalni E16 valni broj k u hCLASS/Mpc uzima fizički komoving k=k_hCLASS hCLASS/Mpc, a radijalnu skalu r_s dobiva u zasebnoj DM14 pozadini i računa x=k_phys r_s,com. To je transparentno označena **cross-cosmology uvjetna prostorna osjetljivost**, a ne samokonzistentno spajanje izvornog CLASS neutrinskog polja i DM14 haloa u fizički B. Za to će trebati jedinstvena kozmologija i relevantna halo kalibracija s masivnim neutrinima.

## Unaprijed zaključani uzorci: bez galaktičkih science cuts

E17D2b0 računa točno dvije masene **numeričke reference** 10¹² i 10¹³ hDM14⁻¹ M☉, odnosno objavljeni pivot i usporedbu za jednu dekadu veće mase. Niti jedna nije proglašena stvarnom ili preferiranom masom LRG ili ELG haloa. Koristi točno tri izvorna E17A matematička z čvora (.945,.95,.955), koji nisu ponovno određeni eBOSS z_eff, i tri unaprijed zaključana log10 c pomaka (−.11,0,+.11) dex. Ukupno 2×3×3 = **18** uvjetnih snapshotova, svaki s originalnih 576 zatvorenih E16 trokuta i točno tri *zasebna* k moda: |k1|, |k2| i K, ukupno **31.104** prostornih Fourierovih vrijednosti.

Za svaku snapshot konfiguraciju izračunava se originalna E17D2a pozitivna NFW težina u(x,c), uvjetni statički oblik −u/x² za nenulti Fourierov mod i uvjetni potencijalni radijalni faktor GM200c/r_s. Ništa od toga nije umnoženo u galaktički B, dodijeljeno LRG/ELG HOD-u, pretvoreno u neutrinsku silu ili interpolirano kao individualna vremenska evolucija. Tri različita z snapshota iste masene reference **nisu tri vremena jednoga praćenog haloa**.

## Originalni CI, neovisni audit i rezultati

[Završni GitHub Actions CI 36386292471](https://github.com/dvlahek/stress-energy-closure/actions/runs/36386292471) **SUCCESS**. [Izvorni program](../scripts/audit_eboss_dr16_a03_e17d2b0_dm14_halo_snapshots.py) izračunao je svih 18 pop. snapshotova i svih 31.104 originalnih prostornih Fourierovih modova. Svi originalni 200ρcrit radius–mass closure testovi, pozitivnost c, originalni trokutni closure, Fourierov znak/povratni smjer i deklarirane SHA roditeljske kontrole prošli su bez promjene masenoga/z/c skupa.

[Neovisni pure-stdlib audit](../scripts/audit_eboss_dr16_a03_e17d2b0_independent_dm14_radial_replay.py), bez izvornog E17D2b0 runnera, NumPyja, SciPyja ili CLASS-a, ponovno izvodi Dutton–Macciò jednadžbe, fizičku kritičnu gustoću i 200ρcrit radijus drugim redoslijedom aritmetike te svih 576 E16 kratko/kratko/dugo geometrija analitičkim skalarnim produktom. Za svaki stvarni originalni k·r_s ponovno direktno integrira pozitivnu NFW radialnu gustoću 4096-panelnom Simpsonovom kvadraturom i zasebno provjerava 2048/4096 konvergenciju. Na svih 31.104 spremljenih vrijednosti najveći skalirani jaz naspram izvorne analitičke Si/Ci formule iznosi **7.172623019743967×10⁻¹³**. SHA-tamper negativna kontrola namjerno izmijenjenoga izvornog izvještaja također je PASS.

Za izvorni matematički z=.95 i **medijan** DM14 koncentracije dvije referentne mase daju:

| Ulazna masena referencija — nije HOD | c200c medijan | r_s,com u Mpc | Najveće apsolutno odstupanje u od 1 na E16 |
|---|---:|---:|---:|
| 10¹² hDM14⁻¹ M☉ | 5.439874383959189 | 0.05999406475939125 | 2.43292582162713×10⁻⁵ |
| 10¹³ hDM14⁻¹ M☉ | 4.563403500551233 | 0.1540783506639813 | 1.1817690358140176×10⁻⁴ |

Na svim unaprijed zaključanim 18 mass/z/illustrative-c slučajevima najveći pojedinačni prostorni |1−u| jest **1.2605532853371404×10⁻⁴**, približno 0.0126 %. To je *uvjetni* rezultat samo za dvije deklarirane masene reference, izvorne E16 modove i DM14 halo porodicu uz eksplicitnu razliku od izvornog CLASS-a. Ne znači da je stvarni LRG/ELG bispektar 0.0126 % ili da je halo-neutrinski wake opažen. U ovom dugovalnom području jednoga navedenog uvjetnog haloa radijalna unutarnja struktura slabo mijenja njegov normalizirani prostorni Fourierov faktor.

[Trajni originalni joint](../source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28/e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json) SHA256 848d3377830724593d75a02972f209085cb9913d03e8b73789f9b6478e6c11e6; [nezavisni certifikat](../source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28/e17d2b0_independent_stdlib_dm14_full_radial_replay.json), Git blob 1c8f1a803f52bf6128ec685afee4f667d11326e1. [Trajni manifest sa svih 20 izvornih/neovisnih JSON izvještaja](../source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28/archive_manifest.json), Git blob 81ed26202502d1539c3c128b17f5b61ba2530929, čuva SHA256 svih 18 punih izvornih snapshotova i dvaju zajedničkih certifikata u source_data/ na audit grani.

## Što je i dalje nedostajući fizički ulaz

Populacijska median c(M,z) kalibracija daje ensemble snapshot statistiku, ne informaciju o individualnom rastu M(a) ili potencijalu Ψ(k,η) istoga haloa; ne dopušta procjenu dΨ/dη iz tri gotovo susjedne maseno-fiksne z točke. Za pravi E17D2b potreban je fizički i kozmološki usklađen vanjski halo mass/assembly prior, nezavisna LRG/ELG halo mass/HOD/selection kalibracija ili realizacije haloa iz usklađene simulacije, zatim samokonzistentno dinamičko gravitacijsko Vlasovljevo rješenje i oba finite-K short-leg tracer couplinga. E17D1 dva jednako-intezivna simbolička izvora ne mogu zamijeniti tu povijest.

Također nedostaju stvarni eBOSS triple-selection/window i nezavisna 3pt kovarijanca estimatora. A03 pair-z RR, A04, originalni SGC reverse i 48k finite-random konvergencija ostaju otvoreni. Originalni 24D odd dvotočkasti vektor **nije** bispektar. **Nema novoga fizičkog B, eBOSS ξ, S/N, detekcijske značajnosti, unblindanja ni novoga opaženog čitanja. Opaženi odd SEALED, originalni E8–E17D2a zamrznuti, main netaknut i PR #1 draft.**
