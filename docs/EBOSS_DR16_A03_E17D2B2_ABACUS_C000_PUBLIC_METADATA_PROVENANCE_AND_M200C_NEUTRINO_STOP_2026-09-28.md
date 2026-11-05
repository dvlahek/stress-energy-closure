# A-03E17D2b2 — AbacusSummit c000 kao izvor halo stabala: nominalno bliska kozmologija, stvarni podaci još nisu prihvaćeni

**28. 9. 2026. · Originalni i neovisni audit JAVNE DOKUMENTACIJE PASS. Nema preuzimanja ili čitanja stvarnoga ASDF kataloga/stabla, novih CLASS izračuna, opaženih galaksija, mockova ili odd podataka. Puni fizički LRG×ELG bispektar BLOCKED, observed odd SEALED.**

## Problem i važna nova mogućnost

Nakon E17D2a prostornoga NFW profila, E17D2b0 vanjske statističke koncentracije iz DM14 i E17D2b1 dvije **nekalibrirane** simboličke povijesti mase preostaje pribaviti fizički neovisno simulirane i kozmološki usklađene individualne halo merger trees. Nije dopušteno tvrditi da dva prijavljena Wechslerova testna parametra predstavljaju stvarne LRG/ELG formacijske povijesti.

[AbacusSummit, službena tablica kozmologija](https://abacussummit.readthedocs.io/en/latest/cosmologies.html) uključuje primarni **c000** Planck-2018 model s jednim masivnim neutrinskim stanjem od 60 meV, nominalnim h=.6736, ωb=.02237, ωcdm=.1200, ns=.9649 i javno dokumentiranim halo katalozima i merger trees. [Glavni rad o AbacusSummitu](https://academic.oup.com/mnras/article/508/3/4017/6366248) i [rad o kvaliteti merger trees](https://academic.oup.com/mnras/article/512/1/837/6541858) dokumentiraju znanstvenu provenijenciju simulatora i veza halo–progenitor. [Službena stranica podataka](https://abacussummit.readthedocs.io/en/latest/data-products.html) navodi približni sekundarni halo redshift direktorij z=.95 i primarne z=.8 i 1.1. To je korisnije od još jedne izmišljene masene povijesti, ali **nije automatski gotov uzorak** za našu fiziku ili opažene galaksije.

[Prospektivni protokol E17D2b2](../source_data/eboss_dr16_a03_e17d2b2_abacussummit_c000_metadata_bridge_prereg_2026-09-28.json), Git blob c1079b666fa8b91af42a8191dbddf56966ea3221, commit 71d4944, prije novoga metadata audita SHA-pina originalni E8 4000q CSV, originalni E16, izvorni E17A CLASS joint, E17D2b1 joint/independent/manifest/protokol te dvije izvorne Python datoteke koje stvarno određuju CLASS kozmologiju i custom ncdm PSD. Navodi službene javne URL-ove i vremenski status njihova čitanja; **ne tvrdi da je SHA-pinao živu web dokumentaciju ni konkretan stvarni dataset**.

## Reproducirani brojčani audit — šest nominalnih podudaranja nisu fizička ekvivalencija

Originalni [CLASS source](../code/class_response_optimize.py) koristi H0=67.36, ωb=.02237, ωcdm=.1200, log-amplitudu izraženu kao As=exp(3.044)×10⁻¹⁰, τreio=.054, Nur=2.0328 i Tncdm=.71611. Izvorni [CLASS argument builder](../code/wake_two_tracer_fisher.py) fiksira ns=.9649, Nncdm=1, mν=.06 eV te učitava vlastitu E8 datoteku ncdm fazno-prostorne raspodjele. E17A numerički koristi ovaj kod na F0/F+/F−, oba originalna E16 kratka kraka i tri izvorna z čvora. Ovo su *izvorne postavke našega projekta*, ne naknadno izabrana zamjena.

Službena [c000 tablica](https://abacussummit.readthedocs.io/en/latest/cosmologies.html) objavljuje h=.6736, ωb=.02237, ωcdm=.1200, ns=.9649, Nur=2.0328, Nncdm=1, As=2.0830×10⁻⁹, τ=.0544 i ων=.00064420. Kod ponovno računa izvornu As iz zaista zamrznutoga Python izraza, bez prepisivanja unaprijed zaokruženoga rezultata:

| Ulaz | Izvorni CLASS (F0/F+/F−) | Službeni Abacus c000 | Status |
|---|---:|---:|---|
| h | .6736 | .6736 | nominalno jednako |
| ωb | .02237 | .02237 | nominalno jednako |
| ωcdm | .1200 | .1200 | nominalno jednako |
| ns | .9649 | .9649 | nominalno jednako |
| Nur | 2.0328 | 2.0328 | nominalno jednako |
| Nncdm | 1 | 1 | nominalno jednako kao broj vrsta |
| As | **2.0989031673191437×10⁻⁹** | 2.0830×10⁻⁹ | izvorni As veći za **0.7634741871888506 %** |
| τreio | **.054** | .0544 | apsolutna razlika −.0004 |

Iako Abacus c000 navodi ων za svoju standardnu neutrinsku distribuciju, originalni E8 F+/F− imaju *ne-termalni oblik* uz fiksne originalne momentne uvjete. **Nije dokazana ekvivalencija originalnih custom PSD perturbacija i Abacus c000 halo rasta.** U službenoj dokumentaciji [o kozmologiji i početnim uvjetima](https://abacussummit.readthedocs.io/en/latest/cosmologies.html) opisano je korištenje CDM+baryon početnoga spektra i neutrinskoga doprinosa kao *glatke* komponente faktora rasta. Zato Abacus c000 merger tree nije nelinearni neutrinski čestični halo wake za bilo koji od naših izvornih F+/F− slučajeva. Nominalna podudarnost h, ωb i mν ne zamjenjuje tu razliku. Čak i referentni F0 zahtijeva zasebno provjerenu jedinstvenu pozadinu, As/τ i način tretiranja neutrinskih perturbacija prije zajedničkoga fizičkog izraza.

## Drugi neovisni STOP — CompaSO skupina nije M200c

[Službeni AbacusSummit data-model](https://abacussummit.readthedocs.io/en/latest/data-products.html) izričito određuje da primarna halo masa dolazi iz broja čestica **N** u CompaSO L1 skupini. Vrijednost **SODensityL1** u ASDF zaglavlju je gustoća za spherical-overdensity grupiranje u odnosu na *srednju kozmološku gustoću* i ovisi o epohi. Ona nije ista stvar kao E17D2a/E17D2b0 **M200c = (4π/3)200ρcrit r200c³**. Ne smije se pripisati ime M200c polju N ili SO_radius bez validirane konverzije. Za realan merger-tree rezultat treba strogo definirati jednu konzistentnu masenu mjeru po vremenu, cleaning verziju i halo branch. Ako se radi remeasurement na 200ρcrit, treba evidentirati stvarnu epoch iz ASDF zaglavlja, radijalne podatke, halo centar, masene jedinice i način tretiranja čišćenja subhaloa. Sparse halo subsample na pojedinom snapshotu nije unaprijed dokazano dostatan za precizan individualni M200c.

## Treći STOP — z=.95 je oznaka direktorija, ne potvrđena fizička epoha

AbacusSummit [data products](https://abacussummit.readthedocs.io/en/latest/data-products.html) navode da se sekundarni z izlazi dobivaju na približnim redshiftima i da stvarni z treba čitati iz *svakoga* zaglavlja; izlazi sekundarnih epoha nemaju potpune halo/field particle pos/vel datoteke, samo halo info i PID podsampove. Merger-tree proizvodi su dokumentirani i dostupnost ovisi o konkretnom boxu/proizvodu, pri čemu [stranica pristupa](https://abacussummit.readthedocs.io/en/latest/data-access.html) upozorava da stabla nisu nužno ponuđena kroz jednostavno web sučelje, već može trebati pregledavanje preko Globusa. Dok realan ASDF i njegova provenijencija nisu očitani, ne može se tvrditi da postoje točno traženi primjerci/epohe, da su filename redshifti jednake epohe ni da je konkretna M200c formacijska povijest dostupna za naše sidro.

## Stvarno izveden originalni/neovisni audit

[Završni GitHub Actions CI 36419675549](https://github.com/dvlahek/stress-energy-closure/actions/runs/36419675549) **SUCCESS**. [Izvorni source-only program](../scripts/audit_eboss_dr16_a03_e17d2b2_abacus_c000_metadata_bridge.py) SHA-provjerava zamrznuti CLASS Python kod i originale E8/E16/E17A/E17D2b1, čita njegove originalne numeričke konstante i izračunava As iz originalnoga log-pivota. Uspoređuje ih s unaprijed zaključanom objavljenom c000 tablicom, izrijekom bilježi svih šest nominalnih podudaranja i dvije ne-nulte razlike. Negativna kontrola odbija sintetičko zaglavlje čak i ako mu se naivno pripišu podobna imena, a zasebno odbija naivno preimenovanje halo mase CompaSO L1 u M200c.

[Neovisni pure-stdlib certifikat](../scripts/audit_eboss_dr16_a03_e17d2b2_independent_c000_metadata_replay.py), bez uvoza izvornoga runnera, NumPyja, SciPyja, CLASS-a ili providerovih datoteka, čita originalni Python **AST**, izvodi As preko **Decimal exp** i reproducira svih osam matematičkih metadata usporedbi te hard-stop polja, uz maksimalni prikazani skalirani ostatak **0**. Namjerno izmijenjen izvorni SHA odbijen je prije usporedbe. To je neovisni audit *metapodataka i provjere provenijencije*, **nije** independent halo-merger-tree simulacija ili opažena znanstvena mjera.

[Trajni SHA manifest](../source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28/archive_manifest.json), Git blob de7c5850d1a776c9cf0f301ba60611808b162169, čuva SHA izvorne [metadata usporedbe](../source_data/eboss_dr16_a03_e17d2b2_archived_CI_2026_09_28/e17d2b2_original_abacus_c000_documented_metadata_bridge.json)  b8da4d6c33235ed82e76f548f2eb052e0fcd751a42672f48d8fffff153ac6926 i neovisnoga certifikata 10c81a9c513b6206beb4cbf1bacdc9ee51946cbe8b62bb534e91d20da88f2d96. Arhiva pripada audit grani. Opaženi odd nije otvaran.

## Strogi uvjeti za sljedeći izvorno fizički korak

Budući E17D2b3 treba zasebno **preregistrirati i odobriti pristup konkretnim simulacijskim podacima**, a zatim kao prvi realni data gate validirati službeni c000 box i phase, trajni dataset ID i checksums, *stvarne* ASDF z/headere i kozmologiju, radijalnu masenu definiciju 200ρcrit, čist/raw halo katalog i stvarni progenitor/descendant branch kroz relevantne epohe. Treba provjeriti može li točno izabrani proizvod uopće dati c200c(z) i mass history istoga haloa. Nakon toga treba odvojeno riješiti nepodudaranje As/τ i izvornih custom neutrinskih distribucija; Abacus c000 glatka neutrinska komponenta **ne** generira samokonzistentan neutrinski halo wake našega Einstein–Vlasov modela. Fizička LRG/ELG HOD/selection i stvarni galaktički 3pt window/covariance ostaju dodatni zasebni problemi.

**Puni finite-K galaktički B, eBOSS ξ, S/N, detekcijska značajnost i prijenos 24D dvotočkastog odd vektora u bispektar BLOCKED; observed odd SEALED. Bez novih opaženih galaksija/kataloga/mockova/CLASS, seedova/cuts ili kontakta s autorima. A03 pair-z RR, A04, SGC reverse i 48k ostaju otvoreni. Main netaknut; draft PR #1 bez mergea.**
