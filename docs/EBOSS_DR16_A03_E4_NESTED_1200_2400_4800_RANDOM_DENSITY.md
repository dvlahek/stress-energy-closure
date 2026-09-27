# EinsteinVlasovNP A-03E4 — ugniježđena provjera gustoće mock randoma

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
