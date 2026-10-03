# A-03E4 — ispravak izvedenog odd sažetka i dijagnostika 24D pilot vektora

**Datum:** 27. 9. 2026. · **Doseg:** isključivo retrospektivni source-only audit *mock-only* podataka nakon završenog E4. Ne mijenja preregistraciju, izvorni E4 JSON, izvorni WSL log, postojeći manifest, mock ID-jeve, random seedove, uzorke, znanstvene binove ni `main`.

## Otkrivena greška izvedenog manifesta

Ranije stvoreni [E4 upload manifest](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json) u polju `odd_ell_median_case_max_abs_difference_by_density_comparison` sadrži prazne objekte za ključeve `"1"` i `"3"`, a usporedbe razina napisane su na istoj razini i zadnjom iteracijom su **prebrisale sažetke dipola sažecima oktupola**. Pogreška je u neovisno izvedenom *manifestu*, ne u originalnom [E4 korisničkom JSON-u od 444 242 bajta](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json). Originalni 18/18 checkpointi i `6×24 ξ` vrijednosti ostaju nepromijenjeni. Umjesto prepisivanja povijesnog izvornika, dodan je [zaseban izvorno vezani E4 ispravak sa svih 18 slučajeva](../source_data/eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json).

## Što je ponovno izračunano

Za svaki originalni ID, NGC i SGC, `ℓ=1` i `ℓ=3`, svih šest fiksnih separacijskih binova i tri originalne usporedbe (2400−1200, 4800−2400, 4800−1200) zasebno je izveden **predznačeni** `Δξ_ℓ(s)` iz dva originalna arhivirana multipola. Ispravak sadrži svih **108 šest-komponentnih** vektora, originalni per-case checkpoint SHA, maksimum apsolutnog pomaka po slučaju, medijan po šest binova i medijane/maksimume/minimume kroz **svih devet** fiksnih mock ID-jeva svake kape. Nijedan slučaj nije uklonjen ili ponovno izvršen.

| Kapa | Medijan `max_s |Δξ₁|`, 2400→4800 | Medijan `max_s |Δξ₃|`, 2400→4800 |
|---|---:|---:|
| NGC | 0,13558150 | 0,28188369 |
| SGC | 0,09229188 | 0,11484821 |

To su **apsolutne razlike izvornog mock-only estimatora**, ne relativne pogreške fizičkog predloška ili izmjerena selekcijska sistematika. Relativni medijani promjene izvorne `6×24 ξ` mreže iz E4 i dalje iznose 0,438741 NGC i 0,434716 SGC pri 2400→4800. Sve 54 originalne RR mreže ostaju 144/144; njihov potpuni support nije dokaz konvergencije.

## 24D retrospektivni pilot, ne novi 18D observable

Ispravak čuva i zajednički opisni vektor predznačenih pomaka za **NGC ℓ1 × 6 s; NGC ℓ3 × 6 s; SGC ℓ1 × 6 s; SGC ℓ3 × 6 s**, točno 24 komponente za svaki od devet izvornih ID-jeva i sve tri usporedbe gustoće. Pri 2400→4800 medijan devet maksimuma apsolutne komponente 24D pomaka jest **0,28188369**, a medijan euklidske norme **0,42216742**. Ovo su samo deskriptivni tehnički pomaci i nemaju jedinice/signifikanciju fizičkog wake predloška; nisu kovarijanca iz devet mockova i ne određuju finalni eBOSS 18D (koji je još nedefiniran).

## Neovisna provjera i sljedeća granica

[Source-only neovisni auditor](../scripts/audit_eboss_dr16_a03_e4_odd_erratum_24d_source_only.py) ponovno izvodi svih 108 originalnih šest-binskih vektora, dvostruko uvjetovane sažetke po kapi i `ℓ`, sve zajedničke 24D vektore, preslikavanje ranijih oktupolnih sažetaka te originalne SHA otiske izvornog izvještaja, manifesta i zasebnog ispravka. [CI `36299976302` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36299976302), uključujući namjerno izmijenjeni dipolni sažetak, predznačeni oktupolni pomak, krivu 18D dimenziju i pokušaj promjene observed-odd oznake. CI **ne otvara** lokalni gzip/FITS, opažene galaksije/randome ili odd vektor.

**Zaključak:** ispravljena je sljedivost brojčanih odd dijagnostika, ali E4 finite-random konvergencija do 4800 nije pokazana. A-03 fizički validirani released-random prozor i A-04 eBOSS kovarijanca ostaju otvoreni. Autore ne kontaktirati, bez novih bulk preuzimanja; observed odd ostaje SEALED.
