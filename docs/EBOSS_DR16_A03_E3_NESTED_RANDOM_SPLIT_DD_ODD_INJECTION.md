# EinsteinVlasovNP A-03E3 — konačni random split i sintetička odd injekcija na mock parovima

**Status 26. 9. 2026.:** protokol zaključan **prije** prvog novog E3 mock FITS retka, izvorni E2 18/18 izvještaj bajtno arhiviran, source-only/sintetički CI [36270704908](https://github.com/dvlahek/stress-energy-closure/actions/runs/36270704908) **PASS**. **Stvarni E3 test na devet mockova još nije izvršen.** Bez kontakta s autorima, bez novih datoteka i bez pristupa opaženom odd vektoru.

**Kanonski preregistrirani protokol:** `source_data/eboss_dr16_a03_e3_nested_mock_random_split_synthetic_dd_injection_protocol_2026-09-26.json`.  
**Runner:** `scripts/audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection.py`.  
**Izlaz:** `eboss_workspace/a03_empirical_window/a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection.json`.

## Nultni uvjeti prije E3

Sve se radi na fiksnih devet realistic EZmock ID-jeva `0001,0125,0250,0375,0500,0625,0750,0875,1000`, za NGC i SGC zasebno. Sva 72 originalna gzip izvora mora proći **punu ponovnu SHA256 provjeru prije bilo kojeg E3 FITS headera ili retka**. Ne mijenjati početni originalni `600D` galaxy uzorak niti originalno odabrane `1200R` po traceru/kapi/ID-ju. Izvorni četiri A-02 sample array otiska, svi originalni E2 ELG chunkovi, baseline `6×24 ξ` i svi E2 originalni baseline/±5 % forward/reverse weighted-pair histogram SHA moraju se reproducirati. Izvorni E2 JSON, njegov manifest, protokol i runner imaju unaprijed zaključane SHA otiske. Nikakav novi download ni ponovno pinanje izvora nije dopušteno.

## E3a: ugniježđeni random i dvije komplementarne polovice

Za svaki fiksni mock/cap/tracer radimo **zasebnu** reproducibilnu permutaciju već odabranih originalnih 1200 RANDOM redaka. Zamrznuta sjemenka:

\[
{\rm seed}=20402026+100000\,{\rm mockID}
+1000\,{\rm capIndex}+100\,{\rm tracerIndex}.
\]

Nakon permutacije zadržava se originalni sortirani A-02 redoslijed prije računanja parova. Fiksne razine su `half_A_600` (prvih 600), `half_B_600` (preostalih 600), `three_quarter_A_900` (prvih 900, dakle sadrži A) i `full_1200` (točno svih originalnih 1200). LRG i ELG randomi poduzorkuju se neovisno po traceru, ali oboje na istoj zadanoj razini. Ove dvije polovice su disjunktne **uvjetno na originalnih 1200 redaka**; nisu dvije neovisne kozmološke realizacije.

Za svaku razinu računaju se ista četiri neovisno normalizirana cross-LS člana i njihova stvarno obrnuta tracer orijentacija, te svi 6×24 binovi uz pozitivnu RR potporu. Ako manja gustoća ima RR rupe, zapisujemo točne nedostajuće `(s,\mu)` indekse, SHA maske potpore te opisne razlike `ξ` **samo na presjeku podržanih ćelija**. Nikad ne popunjavamo nedostajuće ćelije nulama niti projiciramo potpuni šesterobinski `ℓ=0,1,2,3` rezultat bez svih 144 pozitivnih RR ćelija. Originalnih 1200 mora ponovno imati 144/144 potpore i exact E2 SHA.

Na obje polovice i punih 1200 ponovno reproduciramo *unaprijed zadani* E2 `±5 %` ELG RANDOM-only radijalni model; `900` je samo baseline kontrola. Za svaki slučaj prijavljujemo koliko se umjetno inducirani odd odziv mijenja s konačnim brojem randoma. Ne lovimo povoljniji znak, ID, kapu, bin ili split.

## E3b: točno poznata sintetička DD odd injekcija

Na već izračunatom, originalnom **mock-galaxy DD histogramu** injektiramo fiksni potpisani signal

\[
h_{D_LD_E}^{\alpha}(s,\mu_j)=h_{D_LD_E}(s,\mu_j)
[1+\alpha\mu_j],\qquad \alpha\in\{+0.02,-0.02\}.
\]

Ovdje je \(\mu_j\) sredina originalnog potpisanog \(\mu\)-bina. DD normalizacija, sve stvarne galaksije, DR/RD/RR, svi randomi i sve pozicije ostaju nepromijenjeni. Pri nezavisnoj zamjeni LRG↔ELG injektirani se predznak obrće. Izračunati inkrement estimatora na **podržanim** ćelijama mora odgovarati analitičkoj vrijednosti

\[
\Delta\xi(s,\mu_j)=
\alpha\mu_j\frac{DD(s,\mu_j)/N_{DD}}
{RR(s,\mu_j)/N_{RR}},
\]

s rezidualom manjim od \(10^{-11}\) i zasebnom tracer-reversal kontrolom. Kompletnu \(\ell=0,1,2,3\) projekciju izračunavamo samo kada svih 144 RR ćelija ima potporu.

Ova injekcija je **post-binning matematički test DD člana**, nije neovisna fizička galaksijska realizacija ili Einstein–Vlasov wake. Uspješan oporavak dokazuje ispravnost *ovako definiranog* estimatorskog prijenosa na promatranom konačnom paru-kernelu. Ne dokazuje fizičku osjetljivost opaženog surveyja niti može kalibrirati detekciju.

## Jedina lokalna WSL naredba nakon zelenog CI-a

`bash -lc 'set -euo pipefail; cd "$HOME/stress-energy-closure"; test "$(git branch --show-current)" = "audit/eboss-elg-bit8-ra-orientation-20260925"; git pull --ff-only; source .venv/bin/activate; python -u scripts/audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection.py --self-test; mkdir -p eboss_workspace/a03_empirical_window; python -u scripts/audit_eboss_dr16_a03_e3_nested_random_split_dd_odd_injection.py 2>&1 | tee eboss_workspace/a03_empirical_window/a03_e3_local_stdout_2026-09-26.log'`

Sačuvati neizmijenjeni završni JSON i stdout. Naredba izričito provjerava radnu granu, dopušta samo fast-forward `git pull`, aktivira postojeći `.venv`, prvo izvršava sintetički `--self-test`, zatim stvarni runner s nebufferiranim ispisom te čuva **stvarni Python exit status** kroz `pipefail`. Istodobno zapisuje puni stdout/stderr u `eboss_workspace/a03_empirical_window/a03_e3_local_stdout_2026-09-26.log`, bez skrivanja WSL ispisa. Ako se pojavi `INCOMPLETE_STOP`, poslati isti JSON i izvornu pogrešku. Nema mijenjanja seedova, broja parova, ID-jeva ili prizivanja opaženog odd signala radi „popravljanja” ishoda.


## Zasebni post-preregistracijski tehnički popravak (bez promjene protokola)

Nakon izvornog zelenog preregistracijskog CI-a `36270704908`, pri statičkom auditu pronađen je nedjelotvoran izraz `anti_delta < 0`: najveća apsolutna vrijednost ne može biti negativna. Postojeća izravna provjera forward/reverse `ξ` pariteta bila je zasebno aktivna. [Tehnički commit `354952f`](https://github.com/dvlahek/stress-energy-closure/commit/354952f76b78dfbf48b7a7a0bb07133c2f78ff6) zamjenjuje nedjelotvoran uvjet nezavisnom analitičkom provjerom reverse DD inkrementa (prag `1e-11`) i eksplicitno bilježi rezidual jednakosti forward i zrcalnog reverse inkrementa (prag `1e-8`). Sintetički self-test sada provjerava obje dijagnostike. **Originalni E3 JSON protokol, roditeljski A-02/E2 izvještaji, prethodni neuspjesi, ID-jevi, seedovi, binovi, `±0.02` i `±5 %` nisu mijenjani.** Izvorni CI verificira pre-fix commit; za novi runner treba njegov zaseban HEAD CI i lokalni `--self-test` prije stvarnog FITS izračuna.

## Predaja izvornog E3 rezultata i nezavisni acceptance audit

Predati dva neizmijenjena artefakta: originalni `eboss_workspace/a03_empirical_window/a03_e3_mock_galaxy_nested_random_split_and_synthetic_dd_odd_injection.json` i `eboss_workspace/a03_empirical_window/a03_e3_local_stdout_2026-09-26.log`. Ne formatirati niti ponovno serijalizirati JSON prije zapisa njegovih originalnih bajtova, veličine i SHA256. Prva neovisna provjera traži isti zaključani protokol i roditeljske SHA, pre-FITS svih 72 gzip rehash potvrda u izvornom stdoutu, očuvanih svih 18 fiksnih `ID/cap` slučajeva bez promjene redoslijeda, (600_A), (600_B), (900_A), (1200) za oba tracera, disjunktnih i ugniježđenih indeksa/seedova i točnih originalnih A-02/E2 full-1200 baseline/±5 % parnih i `ξ` SHA. Za sve smanjene razine zasebno provjeriti stvarnu RR podršku i eksplicitne nedostajuće indekse; nema nula ni punih multipola na nepotpunoj mreži. Za obje `±0.02` injekcije provjeriti forward i reverse analitičke inkremente, zrcalni paritet i dopuštene projekcije. `STOP` ili neuspjeh ostaju sačuvani u izvornom obliku, bez retuninga ili zamjene slučaja.

Tek nakon neovisnog prolaska audita pohraniti **bajtno identičan izvorni E3 JSON** i **zaseban manifest** u `source_data/`, ažurirati kanonski plan i draft PR #1. Do tada E3 lokalni status ostaje pending i nema fizičke A-03 certifikacije, A-04 18D inferencije ili otvaranja opaženog odd signala.

**Preostala ograničenja:** 600/900/1200 daje konačnu **poduzoračku stabilnost uvjetno na originalnih 1200**, ne pravi limit gustoće produkcijskih randoma. Veće 2400/4800 mock-only random uzorke iz *istih* SHA-pinned FITS izvora i njihovu konvergenciju trebalo bi zasebno preregistrirati ako E3 pokaže preostali finite-RR problem. Devet realizacija nije 18D inferencijska kovarijanca, E3 nije fizički null, p-vrijednost, sigma ni granica isključenja.
