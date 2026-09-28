# A-03E17D2b3b3 — javni Abacus Globus izvor: točan c000 ph000 inventar i ispravak semantike broja stavki

**28. 9. 2026. · [Izvorni CI 36451864487 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36451864487) i [korektivni R1 CI 36452295024 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36452295024). Samo javni, Git-pinnani portal manifest i njegov source code. Nijedna stvarna halo ASDF, PID/particle, merger-tree ni observed galaxies nisu otvoreni ili preuzeti. Observed odd SEALED.**

## Problem i novo provjerljivo otkriće

Abacus source-bundled simulacijska header kopija u E17D2b3b1/b2 omogućila je precizno fiksiranje c000 ph000 nominalnog sekundarnog izlaza z0.950: kopirani originalni Redshift=0.952838237036305, SODensityL1=228.52306365966797 u jedinicama mean density, ParticleMassHMsun=2109081520.453063. Međutim, prije preuzimanja stvarnoga halo kataloga moramo potvrditi da točan c000 ph000 proizvod uopće postoji u službenom javnom inventaru i znati u koju se točno mapu ulazi preko Globusa.

Službeni Abacus [opis pristupa podacima](https://abacussummit.readthedocs.io/en/latest/data-access.html) navodi NERSC Community File System te Globus za selektivan transfer, a za providerov file-level checksum zahtijeva provjeru izvornoga lokalnog checksums.crc32 prema GNU/POSIX cksum. U ovoj fazi nisu korišteni korisnički Globus token, providerov live listing, transfer, prijava ili download halo kataloga.

## Pinned javni portal source i arhivirani inventar

[Originalni prospectively registered b3 protokol](../source_data/eboss_dr16_a03_e17d2b3b3_official_Globus_portal_product_inventory_prereg_2026-09-28.json), prvo Git blob 087055f6e71a50d67727e7e29f9ea1754d7d178e, imao je jedan nedostajući *već zamrznuti* SHA roditelj u parent map. Prvi CI 36451722972 fail-closed u preflightu prije dohvaćanja upstream sourcea. Prereg je prije upstream source checkouta dopunjen u Git blob **42dab76d5e4fb7c8325ece9ea4465775b384da21**, uz zapis originalnog bloba i greške. Originalni sadržaj, odabrani c000 ph000 i z key, podatkovni scope i znanstvene zabrane nisu promijenjeni.

U [službenom javnom portalu na Git commitu f877298344a6988d55f57a7c39f5be14dce63ff3](https://github.com/abacusorg/abacussummit-data-portal/tree/f877298344a6988d55f57a7c39f5be14dce63ff3) (commit 27. 8. 2026.) dobiveni su isključivo source kontrolirani:
- [statički portal manifest](https://github.com/abacusorg/abacussummit-data-portal/blob/f877298344a6988d55f57a7c39f5be14dce63ff3/web/portal/static/data/simulations.table.json), točno 5 097 824 bajta, Git blob cb235ad72ab2516f872bd4ca79b9447a3f4b30e3, SHA256 47dc0f2a2e263c371ac509b14814c77c3bf5c868130eec27c28867fdc70207a2;
- [portal konfiguracija](https://github.com/abacusorg/abacussummit-data-portal/blob/f877298344a6988d55f57a7c39f5be14dce63ff3/web/portal/portal.conf), Git blob 8ae1eaa5a56ddc72dc5313e9fe52ff4233968286, public Globus collection ID **ffc65d7a-0bf9-11ec-90b4-41052087bc27** i endpoint base /;
- [portal views](https://github.com/abacusorg/abacussummit-data-portal/blob/f877298344a6988d55f57a7c39f5be14dce63ff3/web/portal/views.py), Git blob 8c85c76a22746b43b8db2f5e1748a8defbc35f6f, koji potvrđuju origin_id/origin_path query i izgradnju putanja raw/cleaned iz izvornoga manifesta;
- [izvorni graditelj manifesta](https://github.com/abacusorg/abacussummit-data-portal/blob/f877298344a6988d55f57a7c39f5be14dce63ff3/build_manifest.py), Git blob 88dc7c198c1332d2a0204bb0d3c698fd1dcf28a5.

[Izvorni b3 SHA-pinnani runner](../scripts/audit_eboss_dr16_a03_e17d2b3b3_pinned_Globus_portal_inventory.py) potvrdio je identitet točnoga simulacijskog retka, originalne veličine direktorija i raw/cleaned zasebne putanje; četiri sintetičke negativne kontrole odbijaju izmijenjen simname, count, product root i collection UUID. [Originalni arhivirani report](../source_data/eboss_dr16_a03_e17d2b3b3_archived_CI_2026_09_28/e17d2b3b3_original_pinned_portal_c000_ph000_z095_aggregate_inventory.json) SHA256 **39b83642737da1f29e71aade85a2ed7d854dbdb3ff02582b7f1d61d5d92cf340** i [originalni manifest](../source_data/eboss_dr16_a03_e17d2b3b3_archived_CI_2026_09_28/archive_manifest.json) Git blob 4fae0b9418e65e6dba823e22f35fc5d5a9749824 trajno su očuvani na audit grani.

## Važan R1 ispravak: 35 su stavke direktorija, ne verificirana 35 ASDF superslaba

Izvorni b3 prikaz i rani korisnički sažetak naveli su 35 halo_info datoteka. Zatim smo pregledali [stvarni izvorni kod build_manifest.py](https://github.com/abacusorg/abacussummit-data-portal/blob/f877298344a6988d55f57a7c39f5be14dce63ff3/build_manifest.py#L55-L65). On iterira kroz sve stavke direktorija bez filtriranja sufiksa .asdf, broji njihove ukupne stavke i zbraja njihove veličine. Dakle, broj 35 u manifestu **nije** broj dokazano pojedinačnih ASDF haloa. Službeni format ima checksums.crc32 u svakome datotečnom direktoriju, a moguće su i druge pomoćne stavke. Nije dopušteno iz unfiltered counta tvrditi ni 35 ni 34 stvarna ASDF superslaba bez pravoga provider directory listinga.

[R1 unaprijed zaključani erratum](../source_data/eboss_dr16_a03_e17d2b3b3r1_portal_directory_entries_not_ASDF_counts_erratum_prereg_2026-09-28.json), Git blob 534142340197c34160cb970a53cc86d8f7bcf57c, čuva staru izvornu b3 arhivu i eksplicitno bilježi pogrešan raniji opis. [Neovisni R1 AST audit izvornog source buildera](../scripts/audit_eboss_dr16_a03_e17d2b3b3r1_directory_entry_count_erratum.py) potvrđuje da nema extension filtera, a fiktivno dodani .asdf filter i promijenjeni original count se odbijaju. [CI 36452295024 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36452295024), [R1 trajni report SHA256 08255cd1... i novi SHA manifest](../source_data/eboss_dr16_a03_e17d2b3b3r1_archived_CI_2026_09_28/archive_manifest.json), Git blob ba51044967ca5ad70e4173d62361d0c38fa839f0, ispravljaju **interpretaciju**, bez mutacije originalnih b3 artefakata.

| Produkt u fiksnome c000 ph000 z0.950 snapshotu | Sve stavke mape, uključujući eventualni checksum | Ukupna veličina svih stavki direktorija |
|---|---:|---:|
| raw halo_info | **35** | **75 722 740 486 B** |
| halo_pid_A | **35** | **7 933 832 501 B** |
| cleaned_halo_info | **35** | **18 363 791 798 B** |
| cleaned_rvpid | **35** | **488 961 946 B** |

Ova tablica **ne** navodi stvarni broj .asdf datoteka. Za sve četiri mape on ostaje **NEVERIFICIRAN**, kao i broj njihovih stvarnih izdvojenih podataka, pojedinačne veličine, providerovi checksumovi i trenutačna Globus raspoloživost. Izvorni snapshot je vezan uz javni portal commit od 27. kolovoza 2026., ne uz live direktorij na 28. rujna 2026.

## Izvedivi sljedeći korisnički korak: otvoriti usku raw mapu, ne preuzeti cijeli product

Prema SOURCE-PINNANOM portalu, točna raw Globus source mapa je /AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/. [Otvori tu mapu u Globus pregledniku](https://app.globus.org/file-manager?origin_id=ffc65d7a-0bf9-11ec-90b4-41052087bc27&origin_path=%2FAbacusSummit_base_c000_ph000%2Fhalos%2Fz0.950%2Fhalo_info%2F). Odvojena [cleaned_halo_info mapa](https://app.globus.org/file-manager?origin_id=ffc65d7a-0bf9-11ec-90b4-41052087bc27&origin_path=%2Fcleaning%2FAbacusSummit_base_c000_ph000%2Fz0.950%2Fcleaned_halo_info%2F) jest pomoćni raw/cleaned proizvod i nije zamjena za već pripremljen raw ASDF verifier. Ova je poveznica *izvedena iz javno zamrznutog source koda*, nije ovjerena live Globus sesija.

Prije transfera u korisničkom autentificiranom Globusu potrebno je pregledati realni **naziv i veličinu samo jednog dostupnog halo_info_NNN.asdf** i originalni checksums.crc32 iz iste stvarne mape; ako nema checksum registra, ne proglašavati datoteku valjanom. Odabrati samo jedan superslab i checksum file, ne svih 35 stavki niti čitav terabajtni base tar. Mogućnost transfera i odredište ovise o korisnikovoj Globus/NERSC autorizaciji. Nema API tokena niti automatskih korisničkih prijava u ovom stageu.

Kada točno određena stvarna datoteka i **provjereno autentičan** checksum manifest budu lokalno dostupni, [postojeći E17D2b3b2R1 verifier](../scripts/audit_eboss_dr16_a03_e17d2b3b2r1_local_halo_ASDF_checksum_header.py) prvo provjerava cijeli lokalni SHA256 i GNU/POSIX cksum+byte count; potom uspoređuje stvarni ASDF header Redshift, ScaleFactor, OmegaNow_m, SODensityL1 i ParticleMassHMsun/Msun s ranije SHA-pinnanim source-copy vrijednostima. Checksum manifest unesen od korisnika bez neovisne providerove potvrde sam po sebi ne potvrđuje provider origin. Zasebna remeasurement M200c(a) / tree/raw-cleaned/HOD/nu wake ostaje BLOCKED čak i nakon uspješnoga header testa.

**Nema novoga halo download, observed odd SEALED, originalni E8 F0/F± 4000q, E16 576 geometrija, 72 sources i 24 contrast frozen; nema Einstein–Vlasov galaktičkog B ili eBOSS SNR/detekcije. Main netaknut, PR #1 draft bez mergea.**
