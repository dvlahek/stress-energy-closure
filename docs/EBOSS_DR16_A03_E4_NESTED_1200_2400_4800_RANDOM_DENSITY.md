# EinsteinVlasovNP A-03E4 — ugniježđena provjera gustoće mock randoma

## Erratum izvornog izvedenog manifesta — v1.13 (27. 9. 2026.)

Dodatnim auditanjem ustanovljeno je da prethodni *izvedeni* E4 upload manifest ima prazne objekte za `ℓ=1` i `ℓ=3`, a oznake usporedbe gustoće sadrže samo vrijednosti posljednje iteracije (`ℓ=3`): prvotno izračunati dipolni sažeci bili su prebrisani. **Originalni E4 korisnički JSON, log, protokol i 18 checkpointa nisu pogrešni i nisu mijenjani.** Stari manifest sačuvan je bajtno identično. [Zaseban E4 odd ispravak i retrospektivni 24D pilot](EBOSS_DR16_A03_E4_ODD_MULTIPOLE_ERRATUM_24D_2026-09-27.md) i [njegov originalno vezani manifest svih 18 slučajeva](../source_data/eboss_dr16_a03_e4_odd_multipole_summary_erratum_and_24d_drift_2026-09-27.json) točno obnavljaju oba multipola, šest s-binova i sve tri usporedbe; [independent source-only CI `36299976302` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36299976302), s negativnim kontrolama. Numerički zaključak E4 i stanje A-03/A-04 ostaju nepromijenjeni. Opaženi odd je SEALED.


## Konačni arhivski status — v1.11 (27. 9. 2026.)

**E4 lokalni mock-only 18/18 ZAVRŠEN; dva originalna upload artefakta neovisno bajtno auditirana i arhivirana.** [Izvorni E4 JSON, 444 242 bajta](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json), SHA256 `ec45931f815dad8845bae41111f11b3e11cc0a5912f14c3e2b22c2f090e1045f`, Git blob `32246a495506deac5b9dfe71bea51c3db33202b1`. [Izvorni WSL stdout, 11 200 bajtova](../source_data/eboss_dr16_a03_e4_local_stdout_2026-09-27.log), SHA256 `32d93ec0aae6214881ba7852dfe148d2e82a284de5389aeba1dd2c2ba06798a4`, Git blob `31eb941a5a0fb65a8efca2a25530099b0bcc057d`. [Zaseban E4 per-case manifest](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json) bilježi izvorne bajtove, svih 18 SHA checkpointa, fiksne random seedove, RR podršku i opisne mjere konvergencije. [Neovisni source-only original-upload CI `36299224051`](https://github.com/dvlahek/stress-energy-closure/actions/runs/36299224051) **PASS**, uključujući negativne kontrole. Nije ponovno čitao velike gzip/FITS izvore izvan korisnikova WSL-a; izvorni WSL log ima svih 72 ponovno provjerenih SHA otisaka i redoslijed, nezavisno uparenih s prethodno arhiviranim originalnim SHA/byte zapisima.

Originalni A-02/E2/E3 roditelji reproducirani su u svih 18 slučajeva: 72 originalna selected-array SHA, 18 E2 baseline `ξ` SHA i 18 E3 full-1200 `ξ` SHA, svih 144 forward/reverse originalnih E2 baseline histogram SHA/count/norm. Originalnih 1200 randoma očuvano je za oba tracera u svakoj kapi/ID-ju. Iz originalnih eligible populacija neovisno su ponovno generirani svih 36 dodatnih PCG64 seedova i 1200/2400/4800 sorted eligible-index SHA. Na svih **54 baseline mreža** RR potpora je **144/144** i svih **432** forward/reverse pair zapisa prolazi izvorne numeričke kontrole. Izvornih 1200 već je imalo 144/144 u E3: E4 ne treba prikazivati kao otklanjanje originalnog 1200 RR defekta. E3 rupe bile su samo na dijelu poduzoraka 600/900.

**Važan negativan nalaz za numeričku konvergenciju:** Unatoč punoj RR potpori, medijan po devet fiksnih realizacija relativne `L1` promjene cijele `6×24 ξ` mreže pri prijelazu 2400 → 4800 iznosi **0,438741 NGC** i **0,434716 SGC**. Izvorni E4 `relative_L1` koristi nazivnik kao zbroj `abs(ξ_2400)` na istih 144 pozitivnih RR ćelija. Medijan najviše apsolutne šest-binske promjene po slučaju za `ell=1` jest **0,135582 NGC / 0,092292 SGC**, a za `ell=3` **0,281884 NGC / 0,114848 SGC**. Ove opisne promjene znače da **konvergencija estimatora do 4800 randoma nije pokazana**. Nisu fizički izmjerena selekcijska sistematika, stvarna survey nesigurnost, wake signal, fizički odd null, p-vrijednost ili detekcija. Devet mockova nije inferencijska 18D kovarijanca.

Sljedeći zadatak nije retroaktivno birati novi random seed ili proglašavati 4800 stabilnim. Treba zasebno definirati fizički ograničen released-tracer empirical-window/selection-systematics acceptance u A-03 i stvarni 18D inferencijski dizajn u A-04, prije pristupa opaženom odd vektoru. Autore NE kontaktirati, bez novih preuzimanja ili promjene `main`. **Observed odd SEALED.**

Sljedeći odjeljci čuvaju izvorni preregistracijski dizajn i povijesnu pre-execution WSL naredbu. Njihov tekst `E4 lokalni FITS run još nije izvršen` odnosi se na stanje prije izvornog E4 testa, ne na sadašnji status; naredbu ne ponavljati nad već uspješno arhiviranim izvornim checkpointom.


**Status 27. 9. 2026.:** preregistriran novi A-03E4 prije bilo kojeg E4 FITS retka ili dodatnog izvlačenja randoma. [Fiksni E4 protokol](../source_data/eboss_dr16_a03_e4_nested_1200_2400_4800_mock_random_density_protocol_2026-09-27.json) i [E4 runner](../scripts/audit_eboss_dr16_a03_e4_nested_random_density.py) nadovezuju se na [nepromijenjeni, neovisno arhivirani E3 18/18](../source_data/eboss_dr16_a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection_report_2026-09-26.json). Originalna E3 arhiva, njegov protokol, E2, A-02 i svih 72 SHA-pinned gzip izvora nisu mijenjani. **E4 lokalni FITS run još nije izvršen.**

## Problem i opseg

E3 je provjerio dvije disjunktne polovice 600 randoma, ugniježđenih 900 i originalnih 1200, sve uvjetno na 1200. U NGC-u svaka 600 polovica ima RR rupe, a dvije 900 realizacije nisu imale svih 144 ćelija. Zbog toga E3 još ne ispituje konvergenciju prema većoj gustoći. E4 unaprijed određuje **originalnih 1200, dodatnih 1200 (=2400) i ukupno 4800** iz istih ranije certificiranih **full random FITS** izvora. Ovo je **mock-only finite-random** test. Nije novi model selekcije stvarnih ELG/LRG uzoraka, produkcijska geometrijska maska, fizička wake injekcija ili opservacijska inferencija.

Nema novih preuzimanja, dodatnih mockova, novih ID-jeva, promjene redshifta ili maske, promjene znanstvenih binova, kontakta s autorima niti pristupa opaženim galaksijama, randomima ili odd vektoru. Devet fiksnih ID-jeva ostaje 0001,0125,0250,0375,0500,0625,0750,0875,1000, NGC i SGC zasebno.

## Izbor randoma prije znanstvenog izračuna

Iz izvornog full FITS random kataloga uzimamo originalnu A-02 populaciju 0.9 ≤ z < 1.0 nakon nepromijenjene validacije proizvoda WEIGHT_SYSTOT × WEIGHT_CP × WEIGHT_NOZ × WEIGHT_FKP. Izvornih 1200 odabire se *točno* originalnim A-02 seed računom 93127+10000×capIndex+100×tracerIndex+1 i starim numpy.choice algoritmom. Originalni odabir, sva četiri numerička vektora, chunk/weight/source metapodaci, cross-LS ξ i osam originalnih forward/reverse weighted-histogram otisaka moraju imati iste hashove kao arhivirani roditelji **prije** prihvaćanja bilo kojeg novog E4 rezultata.

Dodatni seed, unaprijed fiksiran neovisno o E3 ishodima, jest 20502027+100000×mockID+1000×capIndex+100×tracerIndex. Iz eligible komplementa originalnih 1200 izvlačimo 3600 bez ponavljanja. Prvih 1200 RNG-odabranih dodatnih indeksa dodaje se originalnih 1200 za razinu 2400; svih 3600 dodatnih indeksa daje razinu 4800. U svim razinama za brojanje parova vraća se izvorni sortirani redoslijed eligible retka. LRG i ELG randomi izvlače se **zasebno** po traceru, ID-ju i kapi. Originalni 600D po traceru se nikada ne mijenja.

**Jedini E4 scenarij jest baseline.** Originalni E2 ±5 % weight-stress ostaje zaključan u E2/E3, ali ne proširujemo ga na 2400/4800 u ovom zadatku. Cilj je izolirati konačan broj randoma prije bilo kakve zasebne nove sistematske hipoteze.

## Kriteriji i zapis rezultata

Na sve tri razine računaju se originalne četiri zasebno normalizirane cross-LS komponente i stvarna obrnuta tracer orijentacija, u istim 6×24 binovima. Ne dopuštamo lažni forward/reverse paritet, promjenu izvornog DD histograma, novi science cut, nepodržanu RR nulu ni naknadno biranje ID-ja/kape. Zapisujemo punih 18 fiksnih slučajeva, broj svih pozitivnih RR ćelija, nedostajuće indekse, SHA podrške, sve parne SHA i normalizacije te opisne razlike ξ **samo na presjeku RR potpore** između 1200, 2400 i 4800. Puni ell=0,1,2,3 rezultat moguć je samo ako svih 144 RR ćelija ima pozitivnu potporu.

Ne zahtijevamo unaprijed da je 4800 savršeno popunjen ili da svaka razlika numerički pada monotono. Izvorni devet-mock rezultati nisu 18D kovarijanca i E4 neće proizvesti p-vrijednosti, sigma, fizički null ili detekciju. Ako nema potpore ili preduvjeti zakažu, čuvamo neizmijenjeni STOP i prestajemo, bez post-hoc izmjene seedova.

## Jedina lokalna WSL naredba — isključivo nakon zelenog source-only CI-a

Pokrenuti iz WSL terminala. Naredba čuva vidljiv nebufferirani stdout/stderr i ispravan Python exit status kroz pipefail:

~~~bash
bash -lc 'set -euo pipefail; cd "$HOME/stress-energy-closure"; test "$(git branch --show-current)" = "audit/eboss-elg-bit8-ra-orientation-20260925"; git pull --ff-only; source .venv/bin/activate; python -u scripts/audit_eboss_dr16_a03_e4_nested_random_density.py --self-test; mkdir -p eboss_workspace/a03_empirical_window; python -u scripts/audit_eboss_dr16_a03_e4_nested_random_density.py 2>&1 | tee eboss_workspace/a03_empirical_window/a03_e4_local_stdout_2026-09-27.log'
~~~

Ne izvršavati ponovno nad postojećim E4 INCOMPLETE_STOP s pogreškama. Umjesto toga sačuvati i poslati izvorni STOP JSON i stdout. Nakon dovršetka predati **neizmijenjene** datoteke:

- eboss_workspace/a03_empirical_window/a03_e4_nested_random_density_1200_2400_4800.json
- eboss_workspace/a03_empirical_window/a03_e4_local_stdout_2026-09-27.log

Tek nakon neovisne provjere izvornog byte/SHA, točnih 18 slučajeva, svih 72 pred-FITS gzip hashova, 36 originalnih 1200 sample hashova i 36 fiksnih 2400/4800 ugniježđenja, arhivirati zaseban E4 izvještaj i manifest. **Ne arhivirati ili proglasiti E4 PASS samo na temelju zelenog sintetičkog CI-a.**

Sljedeći fizički A-03 empirijski prozor zahtijeva zasebne tracer-specifične released-random kontrole i prihvatljiv selection-systematics budget. A-04 zahtijeva stvarnu inferencijsku 18D kovarijancu. Oba ostaju otvorena i opaženi odd vektor ostaje SEALED.
