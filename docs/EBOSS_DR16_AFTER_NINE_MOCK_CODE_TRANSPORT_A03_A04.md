# eBOSS DR16 — stanje nakon zatvorenog devet-mock A-02 i preduvjeti za inferenciju

**Datum:** 26. 9. 2026.  
**Kanonski projektni plan:** `EINSTEIN_VLASOV_NP_MASTER_PLAN.md` v1.2.  
**Radna grana:** `audit/eboss-elg-bit8-ra-orientation-20260925`, isključivo draft PR #1.  
**Opaženi odd-k/odd-multipolni vektor:** SEALED; nema dopuštenja za otvaranje.

## Što je sada reproducirano

**A-01:** svih 36 mock-galaxy i 36 same-realization, same-cap, tracer-specific random gzip izvora identificirano je punim SHA256; randomi odgovaraju *ranijim*, neovisnim referencama iz originalnog audita 24. 9. Arhive/manifeste čitati iz `source_data/`.

**A-02:** realni devet-ID cross-Landy–Szalay code-transport na `0001,0125,0250,0375,0500,0625,0750,0875,1000` za NGC i SGC. Sva 72 gzip izvora ponovno su heširana prije prvog mock FITS retka. Oba `0001` izvora reproduciraju originalne uzorčne i histogram SHA otiske. Svih 18 ID/cap slučajeva ima puni RR support `144/144` u fiksnom `6×24` s–μ gridu, zasebne normalizacije sva četiri cross člana i forward↔reverse ξ closure s najviše `1.7763568394002505e-15` apsolutne razlike. Nema neuspjelih slučajeva.

**Neizmijenjeni izvještaj:** `source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_report_2026-09-26.json`; točnih `151892` bajta, SHA256 `15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27`, Git blob `ffd5d801840f90dacded291a7cfaeb1ae70593f6`. Odvojeni manifest: `source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_uploaded_manifest_2026-09-26.json`. Sva su četiri predregistrirana LS člana i sample/hist fingerprints fiksirana po slučaju.

## Sljedeći znanstveni gate A-03 — točan opseg prozora

1. Dokumentirati provenijenciju LRG MANGLE i ELG BRICKMASK izdvojenih selekcija. Zaključane objavljene random datoteke omogućuju *empirijski* tracer-specifični parni prozor u dokumentiranom opsegu, ali ne impliciraju potpunu rekonstruiranu kontinuiranu tvrdu LRG×ELG masku.
2. Za ELG odvojeno dokumentirati izvor i verziju stvarno izvršene produkcijske supplemental bit 8–11 implementacije, konvenciju RA te stvarnu geometriju/primjenu izbačenih `eboss22` ploča `9430-58112` i `9395-58113`. Imena ploča ili javni checkout drugog datuma nisu sami po sebi dovoljni.
3. Ako izvorni proizvodni artefakti ostanu nedostupni, unaprijed ograničiti claim na empirijski prozor objavljenih LRG/ELG randoma i kvantificirati relevantne nepoznanice uz *blind* synthetic/mock-only kontrole. Ne stvarati zajedničku novu masku, novi science cut ili masku iz niza imena datoteka.
4. Prije bilo kojeg novog čitanja *opaženih* galaksija ili odd vektora zahtijevati posebno zaključani ulazni/opservacijski protokol i dopuštenje. A-02 uspjeh nije automatsko dopuštenje za unblinding.

## Paralelni gate A-04 — kalibrirana statistika

Prije velikih preuzimanja ili novih pair runs fiksirati 18D observable, binove, redoslijed, nuisance, teorijski wake predložak i valjanu projekciju through tracer-specific pair window. Zatim odabrati odvojeni adekvatan ansambl *nezavisnih matched realistic mock galaxies i randoms*, s točnim izvorima/SHA, istim produkcijskim selekcijskim modelom, stabilnošću eigen-spektra i kalibracijom pokrivenosti i repa null statistike. Broj realizacija odabrati prema stabilnosti i pokrivenosti, ne samo formalnoj inverzibilnosti. Procijeniti storage/transfer/CPU prije odobravanja većeg preuzimanja.

**Zašto devet nije kovarijanca za hipotezu:** centrirana 18D uzoračka kovarijanca iz devet realizacija ima rang najviše osam. U A-02 su samo kontrolirane 600D/1200R poduzorke i algebarski closure; nema 18D mock odd vektora ni fizičkog modelnog null testa. Nema dopuštenog p-value, `>5σ` exclusion ili wake detection claim-a.

## Opcionalni mock-only opisni QA

Ako želimo vidjeti kako se `ℓ=1,3` vrijednosti mijenjaju kroz devet unaprijed fiksiranih realizacija, prije toga treba **zasebno registrirati** projekciju u istim z/s/μ/LOS/weights/seed pravilima i rehashati svih 72 izvora. Reproducirati svih 18 archived pair-histogram i sample fingerprints prije projekcije. Objaviti sve fiksne ID/cap rezultate, bez re-selekcije, p-vrijednosti, težinski improviziranog coadda ili fizičke interpretacije uzorka od devet. To je dodatni QA, ne preduvjet za tvrdnju o fizičkom prozoru.

## Linija B ostaje odvojena

Paralelni Einstein–Vlasov matched-state causal-response program iz v1.2 plana ostaje aktivan. Eventualni povećani teorijski odziv ne koristi opaženi odd vektor za optimizaciju i ne pretvara numeričku closure preciznost u statističku značajnost.
