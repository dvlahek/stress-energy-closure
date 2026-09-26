# EinsteinVlasovNP — glavni istraživački plan

**Verzija:** 1.0 · **Datum odluke:** 26. 9. 2026. · **Status:** aktivan plan, ne znanstveni rezultat  
**Kanonski dokument:** \`EINSTEIN_VLASOV_NP_MASTER_PLAN.md\` u repozitoriju \`dvlahek/stress-energy-closure\`.  
**Radna grana:** \`audit/eboss-elg-bit8-ra-orientation-20260925\`, isključivo **draft PR #1**. Ne mijenjati \`main\` dok za to ne postoji zasebna odluka.

## 0. Središnja odluka i pravilo protiv “lovljenja sigme”

Einstein–Vlasov identifikabilnost razvijamo kroz **dvije neovisne istraživačke linije**:

1. **A — Opservacijski test:** dovršiti eBOSS DR16 LRG×ELG odd-paritetni/wake test, uz transparentnu mogućnost detekcije, granice isključenja ili statistički neodlučnog nalaza. DESI je zaseban, budući podatkovni skup sa zasebnim selekcijskim i kovarijancijskim auditom. Nijedno unaprijed postavljeno ograničenje nije “garantirani veliki broj sigma”.
2. **B — Fizikalno pojačan kauzalni odziv:** tražiti fizikalno dopuštene parove skrivenih Einstein–Vlasov stanja s istim zadanim početnim izvorom i grav. početnim podacima, ali različitim kasnijim kauzalnim odgovorom. Ispitati može li stvarni opservabilni imati veliki, robustan frakcijski efekt, **bez optimizacije prema opaženim podacima**.

**Najvažnije razdvajanje:** relativna razlika teorijskog kernela, numerička algebra estimatora, statističko isključenje fiksnog modela i detekcija signala četiri su različite tvrdnje. Ne prevoditi jednu u drugu bez zasebne kalibracije.

## 1. Nezaobilazna ograničenja projekta

- Opaženi eBOSS/DESI **odd-\(k\) / odd multipolni vektor ostaje SEALED** dok ne budu ispunjeni svi unaprijed dokumentirani source, selection/window, estimator, mock/covariance i inference uvjeti. Izvještaj “odd data read: false” dio je svakog audita.
- Ne uvoditi naknadne redshift, angular, \(\mu\), separacijske, uzoračke ili tracer cuts prema opaženoj odd amplitudi ili prema povoljnim mock rezultatima. Ne prilagođavati teoriju opažanjima.
- LRG i ELG imaju **zasebne selekcijske funkcije i randome**. Objavljeni odgovarajući \(R_{\rm L}R_{\rm E}\) parovi daju empirijski par-prozor u svojem dokazivom opsegu; oni nisu dokaz jedne univerzalne fizičke tvrde maske.
- Ne zaključivati geometriju nepročitanog mask-piksela iz naziva datoteke, niti ELG semantiku iz LRG MANGLE sektora. Četiri auditirana ELG FITS uzorka nisu svih 19 381 slika.
- SHA256 prvog lokalnog preuzimanja znači **identitet lokalnih izvorišnih bajtova**, ne službeno objavljen SDSS kontrolni otisak. Usporediti s ranijim neovisnim referencama gdje postoje.
- Svaki novi ulaz i svaki novi rezultat: prvo protokol, zatim fiksni identitet izvora, potom izvršavanje, izvorni JSON, neovisna provjera, verzionirana arhiva i tek onda interpretacija. Ne prepisivati neuspjele izvještaje.
- CI/sintetički test, random-only test, test na mock galaksijama, kovarijancijski audit i opservacijski test imaju zasebne oznake statusa. Ne nazivamo algebraički closure “\(10\sigma\)” ili detekcijom.
- Devet unaprijed odabranih mock realizacija služi *code-transport* auditu. Centrirana uzoračka kovarijanca iz njih ima rang najviše \(8\), pa **ne može dati punorangovnu 18D inferencijsku kovarijancu**.
- Cijeli rad odvija se u draft PR #1 dok se ne odluči drukčije; **nikad automatski ne mijenjati \`main\`**.

## 2. Stanje zabilježeno pri donošenju plana (26. 9. 2026.)

**Završeno i arhivirano:**

- Sintetička i random-only provjera četverotermnog, neovisno normaliziranog cross-Landy–Szalay estimatora.
- Stvarna uparena realistic EZmock \`0001\` LRG/ELG mock-galaxy provjera u NGC i SGC s fiksnih \(600D/1200R\) po traceru i kapi u \([0.9,1.0)\). Reverzija tracera zatvara \(\xi(s,\mu)\) numerički na oko \(10^{-15}\) uz 144/144 RR-podržanih ćelija po kapi.
- Opisna \(\ell=0,1,2,3\) projekcija **jedne** realizacije \`0001\`, s rekonstruiranim izvorima, uzorcima i histogramima. Nenulti \(\ell=1,3\) nisu fizički null test niti opservacijska detekcija.
- Puni SHA-only inventar svih **36 mock-galaxy gzip izvora** za prethodno fiksirane ID-jeve \`0001,0125,0250,0375,0500,0625,0750,0875,1000\`. Priloženi JSON i 32 nova prva SHA otiska su arhivirani.
- Iz izvornog uspješnog GitHub Actions audita od 24. 9. izdvojene su prethodne reference punog SHA256 za svih **36 odgovarajućih mock-random gzip izvora**.
- Povijesni \`brickmask\` \`v1.0\` izdanje iz rada i kasniji javni kod razlikovani su na razini izvora; početni ELG bit-8 katalog-level test reproducirao je 15 oznaka, ali povijesno izvršena puna proizvodna maska nije autentificirana.

**Aktivni, još otvoreni zadatak:** lokalni audit svih 36 odgovarajućih random datoteka. Posljednji priloženi međukorak potvrdio je **17/36**; korisnik je nastavio \`--download-missing\` posao. To je *snimka*, ne tvrdnja o njegovom konačnom statusu. Najnoviji nepromijenjeni checkpoint i završni ispis imaju prednost pred ovom brojkom.

**Arhive:** detalji, originalni JSON-ovi, protokoli i prethodne neuspjele provjere nalaze se u \`source_data/\`; operativni slijed je opisan u \`docs/EBOSS_DR16_NINE_EZMOCK_GALAXY_SOURCE_FOLLOWUP.md\`. Ne duplicirati u ovom dokumentu promjenjive hash-popise. Pri povratku na projekt najprije učitati navedene izvorne manifeste, ne pamtiti njihove otiske napamet.

## 3. Linija A — opservacijski test, s detekcijom ili legitimnom granicom isključenja

### A1. Dovršiti ulaznu i algoritamsku provjeru

**Sljedeći neposredni korak:** dovršiti upravo pokrenuti 36-random source-only job. Sačuvati završni ispis i neizmijenjeni \`eboss_workspace/official_mask_inventory/ezmock_nine_matched_random_source_sha.json\`; neovisno potvrditi svih 36 starih referentnih SHA256s, njegovu punu duljinu i status, potom arhivirati točne bajtove i zaključati zaseban manifest. Ako job stane, koristiti postojeći checkpoint, bez promjene ID-ja i bez re-pinninga. Nema FITS redaka u ovom koraku.

**Zatim:** zasebno preregistrirati transport istoga cross-LS koda na ranije izabranih devet realistic mock-galaxy realizacija sa *svojim* odgovarajućim randomima. Zadržati \([0.9,1.0)\), točne prethodne geometrijske binove, LOS, kutni uvjet, produkt težina, \(600D/1200R\), sjemenke i orijentaciju. Izvještavati sva odstupanja i nepodržane RR ćelije; ne izbacivati nepovoljne realizacije. Prije prvog dekodiranja redaka potvrditi **svih 72 odgovarajuća puna mock-galaxy/random gzip identiteta** u zaključenim manifestima.

**Uvjet prolaza:** deterministički identitet ulaza/uzorka, korektne četiri normalizacije \(DD,D R,R D,RR\), ogledalni paritet pri zamjeni tracera, transparentna RR potpora, potpuni statusi po svih devet ID-jeva. To je reproducibilnost implementacije, **ne** 18D inferencija.

### A2. Zatvoriti selekcijski prozor u dokazivom opsegu

- Odvojeno održati eBOSS DR16 LRG i ELG angular/radial/random contract.
- ELG: dokumentirati stvarnu DR16 proizvodnu verziju supplemental bitova 8–11, bit-8 koordinatnu konvenciju i točnu geometriju i originalnu primjenu zasebno isključenih \`eboss22\` ploča \`9430-58112\` i \`9395-58113\`. Njihovi identifikatori su poznati; puni prostorni veto **nije** izveden samo iz oznaka.
- Ako izvorni proizvodni kod nije dostupan, jasno odvojiti ograničeni **empirijski prozor objavljenih tracer-specifičnih randoma** od tvrdnje da je rekonstruirana potpuna kontinuirana fizička maska. Ne stvarati novu zajedničku binarnu masku.
- Napraviti unaprijed definirane blind synthetic-injection i mock-galaxy tests za moguću selekcijsku lažnu odd komponentu i osjetljivost na dostupne documented random/n(z) varijante. Ništa ne podešavati prema opaženom odd rezultatu.

**Uvjet prolaza:** dokazivo ograničen estimator/window model i dokumentiran raspon nesigurnosti. Ako puni produkcijski dokaz ostane nedostupan, ograničiti opservacijsku tvrdnju ili zadržati opaženi vektor zatvorenim.

### A3. Dizajnirati stvarnu inferencijsku kovarijancu

- Unaprijed definirati **opservacijski podatkovni vektor** (trenutačni cilj 18 dimenzija), redoslijed njegovih elemenata, prijelaz mock→observable, fiducijal i odgovarajuće nuisance parametre.
- Izraditi *zaseban* dovoljno velik, nezavisan, source-SHA-pinned realistic mock-galaxy ansambl s istom selekcijom i uparenim randomima. Broj mockova definirati na temelju konvergencije eigen-spektra, stabilnosti inverzije i kalibracije distribucije test-statistike, ne minimalnim zahtjevom da matrica bude formalno invertibilna.
- Dokumentirati eventualnu regularizaciju, shrinkage, korekciju nesigurnosti procijenjene kovarijance, tail calibration i mock-to-data systematics. Svaka odluka ide u protokol **prije** otvaranja opaženog odd vektora.
- Izvesti end-to-end mock-only false-positive, injection/recovery, cap/selection konzistenciju, perturbacije random težina i audit estimator-window degeneracija.

**Uvjet prolaza:** validirana 18D kovarijanca i unaprijed kalibriran test s dokumentiranom osjetljivošću/pokrivenošću. Devet mockova nije dovoljno.

### A4. Zaključati dva moguća znanstvena ishoda *prije* opaženih podataka

**Detekcijski test:** fiksirati predložak, skalu/amplitudu, nuisance, test-statistiku, null distribuciju, sve intervale i višestruka testiranja. Ako postoje provjere više konfiguracija, uključiti look-elsewhere ili unaprijed objasniti da su samo kontrolne. Detekciju tvrditi samo ako stvarni test i sustavne provjere zadovolje unaprijed definirani prag.

**Test isključenja:** definirati koju točno Einstein–Vlasov/wake predikciju \(A_\star\) testiramo i kako se mapira na opaženi vektor. Izračunati granice i pokrivenost kalibriranim postupkom. U jednoparametarskoj Gaussovoj aproksimaciji udaljenost fiksnog modela može biti
\[
Z_\star = \frac{|A_\star-\widehat A|}{\sigma_A},
\]
ali to je smislena sigma **samo** uz valjanu kovarijancu, test, template uncertainty i nuisance/tail kalibraciju. Granica na **95 % CL nije 5σ**. Null rezultat nije dokaz da su svi mogući wake odgovori isključeni; isključujemo samo fiksirani parametarski model ili područje parametara.

**Treća dopuštena mogućnost:** statistički neodlučan nalaz uz pošten sensitivity limit. Ne reinterpretirati ga naknadno kao detekciju ili potvrdu teorije.

### A5. eBOSS i DESI ostaju zasebne podatkovne etape

Prvo kvalificirati eBOSS DR16 u navedenom opsegu. Puni DESI katalog zahtijeva vlastiti public-release/version SHA, tracer-specifične angular/radial/random contracts, veto/prozor, mock ansambl, kovarijancu i nezavisni blind gate. Nikakva eBOSS kalibracija automatski ne vrijedi za DESI.

**Isporuka linije A:** auditabilan opservacijski rukopis s reproducibilnim podacima, estimatorom i jasno imenovanim jednim od tri ishoda: detekcija, ograničenje ili neodlučan rezultat.

## 4. Linija B — pojačan, ali fizički dopušten kauzalni odziv

### B1. Definirati identifikabilnost bez nejasnog “source matchinga”

Konstruirati dva različita početna kinetička stanja \(f_A\neq f_B\) koja zadovoljavaju točno navedene početne matching uvjete (npr. isti \(T_{\mu\nu}(t_0)\) i isti dopušteni gravitacijski Cauchyjevi podaci), ali imaju različite više momente/anisotropiju distribucije. Uvjeti moraju zadovoljavati Einsteinove constraint jednadžbe, odgovarajući gauge i Einstein–Vlasov evoluciju. **Ne tvrditi** da različita buduća metrika slijedi iz dva izvora koji su identični u cijeloj povijesti; matching je definiran na navedenoj početnoj površini ili konačnom dostupnom skupu opažljivih informacija.

Kriteriji: nenegativna fizikalna \(f\), ispravna normalizacija i podrška, maseni/energijski budžet, dopuštena anizotropija, evolucijska regularnost/stabilnost i eksplicitna područja valjanosti perturbativnog ili numeričkog opisa.

### B2. Sistematski skenirati unaprijed definirane fizičke režime

Kandidati za **testiranje**, ne unaprijed potvrđeni izvori velikog signala: anizotropni samogravitirajući bezsudarni sustavi; kolektivni modovi i phase mixing; ograničeno pojačanje blizu rezonancije; vrijeme-konačan odziv pri promjeni distribucije uz isti početni niži moment. Za svaki režim prvo specificirati astrofizički prior, relevantni vremenski/prostorni raspon i stabilnost. Ne pojačavati signal nefizičkim negativnim distribucijama, proizvoljnim forcingom ili izlaskom iz režima valjanosti.

Početni rezultati oko **3,81 %** relativne razlike kauzalnog Greenova kernela i **0,80 %** metričkog odziva odnose se na ilustrativni konstruirani par iz postojećeg teorijskog rada. Nisu univerzalna gornja granica. Brojeve od **10–50 %** tretirati kao *istraživačku hipotezu o mogućem režimu*, ne cilj koji mora biti ostvaren ili uvjet za odbacivanje slabijih valjanih rezultata.

### B3. Tri odvojene razine rezultata

1. **Teorijska identifikabilnost:** dokazati ili numerički validirati da isti zadani početni niži moment ne određuje cijeli budući response.
2. **Modelni odziv:** na fiksiranoj normi, vremenskom prozoru i gaugeu procijeniti \(\|\Delta G\|/\|G\|\), \(\|\Delta g\|/\|g\|\) i apsolutne amplitude. Odvojiti relativno povećanje zbog malog nazivnika od fizički velikog odgovora.
3. **Opservabilna razlučivost:** mapirati oba stanja na unaprijed odabran stvarni observable \(\mathcal O\), modelirati instrument/survey prozor i nuisance parametre, te procijeniti \((\Delta\mathcal O)^{\mathsf T}C_\mathcal O^{-1}(\Delta\mathcal O)\) **samo** kada je stvarna kovarijanca definirana i kalibrirana. Velika razlika kernela sama po sebi ne znači veliki \(S/N\).

### B4. Pravila za dokaz robusnosti i odluku

Zaključati fizikalnu klasu, sken i mjere prije usporedbe s opaženim podacima. Provjeriti grid/time-step konvergenciju, conservation/constraints, mode/gauge artefakte, perturbativno protiv punog režima i neovisnu rekonstrukciju responsea gdje je izvedivo. Posebno testirati matched izvor i numerički integrirane momente tijekom evolucije.

Ako robustan veliki observable postoji: izvesti jasnu predikciju i **novi, nezavisni blind opservacijski protokol**, ne podešavati postojeći DESI/eBOSS wake test. Ako ne postoji: objaviti/iskoristiti identifikabilnost i rigoroznu granicu pojačanja u ispitanoj fizičkoj klasi; ne umjetno stvarati 10–50 %.

**Isporuka linije B:** dokazani fizički dopušten source-matched par ili klasa, kontrolirana dinamika, kvantificiran kauzalni i metrički odziv, te eksplicitna predikcija i izvedivost opservabilnog kanala ili jasna ograničenja njegove razlučivosti.

## 5. Zajednička integracija i rukopis

Radna priča slijedi: **problem identifikabilnosti → granica rekonstrukcije iz nižih momenata → kauzalni mehanizam → dokaz/modelni odziv → realistični observable → blind empirijski test ili poštena granica**.

Liniju B ne vezati uz rezultat linije A post hoc. Ako pronađeni jači opservabil nije isti kao eBOSS/DESI wake vektor, voditi ga kao drugi, neovisno preregistrirani test. Ako opservacijska linija daje samo ograničenje, prikazati to kao fizički informativno ograničenje precizno definirane klase. Ne predstavljati cilj časopisa ili “veliku sigmu” kao rezultat.

## 6. Evidencija odluka i uvjeti prelaska

| Oznaka | Zadatak / artefakt | Dokaz potreban za zatvaranje |
|---|---|---|
| A-01 | 36 uparenih random SHA za fiksnih devet ID-jeva | Neizmijenjeni završni lokalni JSON, svih 36 ranije fiksiranih SHA podudarnosti, arhiva i manifest |
| A-02 | Devet-mock real-galaxy code transport | Preregistrirani protokol, svi 72 input SHA, deterministički estimator i nepobrisani fail slučajevi |
| A-03 | Empirijski/proizvodni LRG×ELG prozor | Tracer-specifičan provenijencijski contract, testovi uvjetne random selekcije i poštena ograničenja fizičke maske |
| A-04 | 18D mock inferencija | Dovoljno nezavisnih mockova, konvergencija/regularizacija/tail kalibracija, blind injection i false-positive kontrole |
| A-05 | Opaženi odd unblinding | Pisani zaključani protokol, preduvjeti A-01–A-04 i zapis dopuštenja prije prvog pristupa |
| A-06 | Detekcija/isključenje/neodlučno | Jedinstveni unaprijed zadani test, predložak i interval s odvojenim sustavnim zaključkom |
| B-01 | Matched Einstein–Vlasov klasa | Pozitivna \(f\), matching i gravitacijske constraints, početna/astrophysical plausibility |
| B-02 | Kauzalni response sken | Predefinirani prostor, numerical convergence, apsolutne i relativne razlike bez optimizacije prema podacima |
| B-03 | Realni opservabil i noise | Fizikalna projekcija, instrument/survey, nuisance, validirana \(C_\mathcal O\), neovisni test |
| B-04 | Integracija u znanstveni rad | Rezultat samo u dokazivom dosegu, obje linije pravilno odvojene |

**Trenutačni prioritet:** A-01 do njegovog potpunog izvještaja. Linija B može se razvijati paralelno teorijski i sintetički, ali ne smije koristiti još zapečaćeni A-05 odd vektor za podešavanje modela.

## 7. Kako ovaj plan ponovno koristiti i održavati

- U sljedećem razgovoru unutar projekta dovoljno je napisati **“Otvori glavni EinsteinVlasovNP plan i nastavi od trenutačnog gatea”**. Referentni dokument je ova datoteka u \`dvlahek/stress-energy-closure\`; prije konkretnog sljedećeg poteza provjeriti najnoviji draft PR, povezane izvorne JSON manifeste i zadnji korisnički checkpoint.
- Novo saznanje upisati u odgovarajući protokol/arhivu; **ovdje ažurirati status, odluku i datum**, bez retroaktivnog mijenjanja značenja prvobitno preregistriranih testova. U slučaju konceptualne izmjene napraviti verziju 1.1 s bilješkom što je promijenjeno i kada.
- GitHub dokument je trajni projektni *source of truth*, ne tvrdnja da je interna ChatGPT memorija ručno izmijenjena. Za pristup i izvan repozitorija korisnik može istu Markdown datoteku priložiti u ChatGPT Project files.
- Obavezna završna provjera prije izjave o rezultatu: **što je zapravo izmjereno, na kojoj populaciji/realizaciji, kojim prozorom, s kakvom kovarijancom, koji su nezavisni testovi prošli i je li opaženi odd vektor bio otvoren?**

### Dnevnik verzija

- **v1.0 — 26. 9. 2026.** Zabilježena odluka o dvjema linijama, razdvojene detekcija/isključenje/closure/modelni response, zamrznut trenutačni eBOSS/eZmock gate i definirani prolazni uvjeti za budući observable i inferenciju.
