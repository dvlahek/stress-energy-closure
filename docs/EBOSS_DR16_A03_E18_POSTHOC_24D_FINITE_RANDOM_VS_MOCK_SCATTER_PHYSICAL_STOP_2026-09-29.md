# A-03E18 — promjena 24D odd pilota s gustoćom randoma, originalnih devet E4 mockova

**29. 9. 2026. · ORIGINAL E4 MOCK-ONLY POSTHOC DESCRIPTIVE RESULT.** Ovo je dodatna retrospektivna numerička dijagnostika već završenoga E4 nakon što su E4 i njegova 0001 shema već bili poznati. **Nije** prospektivni test znanstvene hipoteze, mjera eBOSS fizičkoga prozora, covariance/inference QA, nova selekcijska odluka, detekcija ili ekskluzija. Bez WSL-a, novih FITS redaka, novih randoma, maske, seedova, science cuts, opaženih galaksija i opaženog odd vektora.

## Problem: E4 relativna ξ-mreža nije bila u jedinicama odd opservablea

Prethodni E4 mock-galaxy 1200→2400→4800 random test zadržao je isti originalni 600D uzorak u svakom od devet originalno fiksiranih matched mock ID-jeva, po dvije kape i oba tracera. Svih 18 ID/cap slučajeva i sve tri razine imaju 144/144 RR podršku. Veliki relativni pomaci ukupnog 6×24 \`ξ\` estimatora sami po sebi nisu fizički budžet odd komponente. Za usporedbu u istim, već postojećim dimenzijama koristimo *isključivo* postojeći 24D pilot: \`NGC ℓ1 (6 s), NGC ℓ3 (6 s), SGC ℓ1 (6 s), SGC ℓ3 (6 s)\`; originalni \`s\` centri \`30,50,70,90,110,130 Mpc/h\`.

[Retrospektivni scope/frozen source ugovor](../source_data/eboss_dr16_a03_e18_posthoc_nine_mock_paired_random_24d_stability_protocol_2026-09-29.json) Git blob \`9b13de770d4c9b051f8b6af7c5a6e2a342df77e3\` eksplicitno kaže da je 0001 shema već pregledana prije njega: bez lažne preregistracije. [Izvorni E4 444242-bajtni JSON](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json) Git blob \`32246a495506deac5b9dfe71bea51c3db33202b1\` i [raniji neovisni signed 24D erratum](../source_data/eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json) Git blob \`0baf9c5a009bf8f5e78f5725e44395bad6f91715\` ostaju nepromijenjeni. Svih 27 parnih 24D vektora \`9 IDs × 3 level contrasts\` rekonstruirano je iz E4 originalnih \`ℓ=1,3\` binova i verificirano prema erratumu na maksimalnu razliku \`≤1e-12\` isključivo kao aritmetički audit.

## Metrika, prije bilo kakve statističke interpretacije

Za i=1,...,9 i fiksnu komponentu j definiraj \`d_ij = x_ij(high R) - x_ij(low R)\`. U grupi od q komponenti \`RMS_paired = sqrt(sum_ij d_ij²/(9q))\`. Za \`x(high R)\` opisni centrirani mock rasap je \`RMS_between = sqrt(sum_ij [x_ij(high R)-mean_i x_ij(high R)]²/(8q))\`. Omjer \`RMS_paired/RMS_between\` je **isključivo opisni omjer dviju skala**. Ne predstavlja standardiziranu test statistiku, \`S/N\`, mjerenje fizičkog leakagea niti nezavisne komponente varijance: randomi su ugniježđeni, a high-R rasap također uključuje finite-random varijaciju. Nema naknadno zadanog prihvatnog praga.

| Ista originalna 24D pilot-grupa, svi 9 ID-jeva | RMS parnoga pomaka | RMS među mockovima na višoj R razini | Opisni omjer, nije sigma |
|---|---:|---:|---:|
| 1200 → 2400 | 0.183151 | 0.209158 | 0.875660 |
| 2400 → 4800 | 0.094742 | 0.180712 | 0.524269 |
| 1200 → 4800 | 0.201284 | 0.180712 | 1.113834 |

Najbliži E4 korak \`2400→4800\` zasebno: NGC \`ℓ1\` 0.085309/0.176219 (omjer 0.484106), NGC \`ℓ3\` 0.148239/0.258862 (0.572657); SGC \`ℓ1\` 0.054397/0.102062 (0.532979), SGC \`ℓ3\` 0.060768/0.148824 (0.408321). Sva 24 pojedinačna \`j\` i svih 9 izvornih ID-jeva bez filtriranja nalaze se u [malom arhiviranom rezultatu](../source_data/eboss_dr16_a03_e18_posthoc_nine_mock_24d_paired_random_vs_mock_scatter_descriptive_result_2026-09-29.json). Opisni RMS pomaka joint-24D pri \`2400→4800\` je manji od \`1200→2400\`, ali to **nije** dokaz asimptotske konvergencije ili fizički opravdane tolerancije.

## Dokaz izvedbe i granica

Izvorna Galaksijska \`D1D2\` weighted histogram SHA mora biti ista na sve tri razine za isti ID/cap; sva 144 RR moraju biti pozitivno podržana, bez zero-fill. Novo računanje radi samo nad već arhiviranim malim E4 JSON-om, a svi signed vektori odgovaraju ranijem erratumu. [Neovisni standard-library CI 36541413140 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36541413140) zasebno izračunava 3 kontrasta, 7 grupa/kontrast, svih 24 komponenti i devet ID-jeva; pet negativnih kontrola odbija RR rupu, promijenjen DD, pokvaren signed erratum, krivi broj i nedostajući cap. Prvi CI neuspjeh 36541359859 posljedica je isključivo pogrešnoga naziva polja u neovisnom čitaču; originalni B7/E4/erratum i E18 izračunati rezultat nisu prepisani.

**Ne kombinirati ovaj 600D E4 s 48k E7:** raniji E7 \`mock0001 4800→48000\` koristio je puni galaktički uzorak. Stoga njegovo \`DD\` i njegov 48k odd vektor nisu upareni random-only referent za E4 600D.

Za fizičku A-03 treba apsolutno normiran EV/standardni prewindow even+odd teorijski vektor, provjerljiv tracer-specific released \`R_L R_E\` operator i odvojen finite-random/selection budget kroz 24D ili drugi *prije* finalne analize fiksiran observable. Za A-04 je potreban zaseban opravdano nezavisan same-release matched eBOSS mock ansambl i validna covariance/tail kalibracija. Devet mockova imaju centrirani 24D sample rank \`≤8\`; njihova međusobna disperzija nije inferencijska kovarijanca. Nastavak ne smije biti ponovno računanje istih 600D parova s post-hoc većim seedovima niti otvaranje opaženog odd vektora. Zadržati originalnu korisnikovu empirijsku released-random rutu bez kontaktiranja autora, glavnu granu netaknutom i PR #1 u draftu.
