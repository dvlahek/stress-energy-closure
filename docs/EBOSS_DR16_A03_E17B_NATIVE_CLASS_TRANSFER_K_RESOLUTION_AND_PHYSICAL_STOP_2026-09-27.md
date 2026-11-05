# A-03E17B — numerička provjera native CLASS transfer-k potpore oba kratka kraka

**27. 9. 2026. · Odvojeni zaključani E17B protokol. E17A linearni CLASS rezultat ostaje izvorno neizmijenjen. Puni fizički halo/tracer bispektar BLOCKED.**

## Problem

E17A je na izvornim 576 E16 zatvorenih trokuta dobio originalni CLASS `P_{cb}(|k_1|,z)`, `P_{cb}(|k_2|,z)`, linearni direct-vTk transfer i newtonovski `delta_cb` za tri izvorna neutrinska stanja i tri redshifta. E17A nije imao neovisnu provjeru numeričke osjetljivosti linearne interpolacije direct-vTk s native CLASS transfer-k mreže od 81 čvora. Neovisni E17A replay potvrđuje istu 81-čvornu realizaciju, ne njezinu rezolucijsku konvergenciju.

[E17B prospektivni protokol](../source_data/eboss_dr16_a03_e17b_native_transfer_k_resolution_prereg_2026-09-27.json), Git blob `d0389e336de9e9964cade589fc17594a595f2e17`, zasebno je commitan prije novih numerički rafiniranih CLASS izvođenja. Fiksirani su originalni CLASS commit `e85808324f51fc694d12e3ed7439552a3c3f9540`, E8 4000q izvor, E16 12×48 geometrija, izvorna E17A tri state JSON SHA, z=(0.945,0.950,0.955), tracer redoslijed LRG→ELG te originalni R16 direct-vTk rang. Ne mijenjaju se F±, masa, redshift matching, 30% cap, izvorni E14/E15 `S` ili Eq20 konvencija.

## Prospektivno odabrana CLASS provjera

Za svako zamrznuto neutrinsko stanje zasebno se pokreće CLASS s istim fizikalnim parametrima i s dvije **isključivo numeričke** preciznosti. Prva ima `k_per_decade_for_pk=40`, `k_per_decade_for_bao=140`; druga `80` i `280`. Originalni E17A 81-node izlaz ostaje zamrznuta referenca. Za svaku rafiniranu native mrežu provjeravaju se rast broja čvorova, raspon svakoga stvarnog `|k_i|`, stroga monotonost i jednaka mreža za tri originalna z čvora. Izvornim `CLASS.pk_cb_lin` izračunavaju se oba stvarna kratka kraka izravno, bez interpolacije tri E14 sidra.

Rafinizirani `t_ncdm[0]-t_cdm` i `delta_cb` dobivaju se iz stvarno spremljenih CLASS native nizova i interpolacijom samo unutar verificiranoga razmaka. U originalnom E17A ostaje zasebna R16-filtrirana direktna relativna LOS brzina. **Jedan direktni kratkovalni transfer nije izmjerena halo brzina i nije zamjena za originalni kondicionirani R16 rang.**

Spremljeni su puni native k/transfer nizovi, svi kračni `P_cb` i interpolirani transferi za oba stupnja te dijagnostike za svih 3×576×2×3 kračnih/redshift uzoraka. Razdvojene su razlike E17A→ultra i high→ultra. Blizu nule koristi se unaprijed zaključana skalirana norma; obje sirove vrijednosti i najgori geometrijski slučaj ostaju zapisani. Preregistrirani engineering kriteriji 1% E17A→ultra, 0,3% high→ultra i 0,1% za originalna E14 kratkovalna `P_cb` sidra **nisu fizičke granice pogreške**. Ne mijenjati ih prema dobivenom rezultatu.

Zaseban neovisni standard-library program bez importiranja originalnog E17B koda iz izvornih pohranjenih native CLASS nizova nanovo računa sve interpolacije, SHA provjerava E17A/novi E17B i ponovno izračunava maksimalne metrike. Ovo je matematički neovisni replay pohranjenih novih CLASS podataka, ne drugi fizički halo izračun niti drugi neovisni CLASS izvršni backend. Originalni E17A već ima odvojeni svježi CLASS replay.

## Granica fizičke interpretacije

Promjena CLASS rezolucije ispituje ukupnu osjetljivost linearnog numeričkog rješenja i prijenosne interpolacije. Čak i dobro slaganje tri mreže nije matematički certificirana pogreška i ne stvara retardirani Vlasov halo/tracer kernel. Fizikalna sprega `long density/wind → both short halo/tracer legs`, LRG/ELG HOD i bias, evolution/magnification/GR nuisancei, fizički eBOSS trotočkasti selekcijski prozor i neovisna 3pt kovarijanca ostaju zasebno otvoreni. Svi zaključci o punom `B`, `xi`, S/N i značajnosti ostaju **BLOCKED**.

Opaženi odd je **SEALED**. Ne preuzimati nove kataloge/mockove, ne uvoditi nove seedove ili science cuts, ne kontaktirati autore, ne mijenjati `main` i ne mergeati draft PR #1. E17B neće promijeniti izvorne E8–E17A podatke. Numeričke neuspjehe čuvati kao zasebne CI događaje, bez preimenovanja u PASS.

## CI i izvorni rezultat

Prvi [CI 36344819090](https://github.com/dvlahek/stress-energy-closure/actions/runs/36344819090) prošao je E8–E17A SHA preflight, izvorni CLASS build i oba rafinirana FD CLASS koraka, ali je zaustavljen Python `AttributeError` zbog pogrešno imenovanoga direktorijskog atributa `e14.E14DIR` u novom E17B čitaču. Nije proizveo izvorni numerički izvještaj ni završni audit. Naknadni commit `f577afa4a2276f2e29a44bacb366cf915f640f6c` zamjenjuje samo čitački put s postojećim `e17.E14DIR`; izvorni protokol, preciznosti, F stanja i acceptance kriteriji nepromijenjeni su. Taj prvi CI ostaje zabilježen kao **FAIL**.

## Završni izvori i stvarni numerički nalaz

[Završni CI 36344909227](https://github.com/dvlahek/stress-energy-closure/actions/runs/36344909227) **SUCCESS** na nepromijenjenom izvornom E17B protokolu i originalnim E8–E17A roditeljima. Originalne CLASS native k mreže imaju **81 → 238 → 475** čvorova za sva tri zamrznuta stanja. E17B svih 576 geometrija, oba kratka kraka i z=(.945,.95,.955) ponovno evaluira i uspoređuje; `P_cb` se poziva izravno na stvarnom `|k_i|` u svakom tieru.

Sljedeća tablica daje **najveće** odstupanje preko FD, F+ i F− i cijeloga zaključanog kračnog/z skupa; ovo su unaprijed definirane **skalirane numeričke dijagnostike**, ne pogreške stvarnog fizičkog bispektra:

| Polje | Izvorni E17A 81 → ultra 475 | high 238 → ultra 475 |
|---|---:|---:|
| `P_cb` | `1.5070636903e-5` (0,0015071 %) | `1.8116202781e-6` (0,0001812 %) |
| direct `theta_ncdm−theta_cdm` iz vTk | `0.0010909734617` (0,10910 %) | `0.0001585246214` (0,015852 %) |
| Newtonian `delta_cb` | `0.0008990670834` (0,089907 %) | `0.0001119249893` (0,011192 %) |

Najveća razlika rafiniranog `P_cb` na **tri originalna E14 kratkovalna sidra** iznosi `1.3452757763e-5`. To je numerički kontroliran rezultat pri izmjeni precision parametara, ne izmjena originalnih E14 arhivskih bajtova. Sva tri zamrznuta stanja ispunila su unaprijed zaključane E17B *engineering candidate* pragove, bez retuninga.

[Neovisni standard-library replay](../source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_independent_native_array_scalar_replay.json) ponovio je **41.472** native-array interpolacije i **31.104** usporedbe izvornih E17A uzoraka. Najveći skalirani interpolation replay gap `2.1883002059282083e-16`, a maksimalni discrepancy svih originalnih QA maksimuma `0`. Certifikat SHA256 `6b1e99a6d695378453dd47127e642c0d286d08dd75187b1cf0e00c3005d4ea4e` nije drugi neovisni CLASS izvršni backend: rekoristira autentične originalne rafinirane native transfer nizove.

[Originalni objedinjeni izvještaj](../source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/e17b_original_joint_native_CLASS_transfer_resolution.json) SHA256 `156763dfb81e2a8c5c9e3615f5022749fa681f458a2753fbfa5f46b462e20658`. [Potpuni arhivski manifest](../source_data/eboss_dr16_a03_e17b_archived_CI_2026_09_27/archive_manifest.json), Git blob `88920dd194febb002f280d676b9df42687b2f949`, SHA-zaključava i tri puna state JSON izvještaja (FD `ff9411b6a1bdf7c3d2b989f61b7741254389faa05bf61e0c8487b77eddb4ce2d`, F+ `66b6a344906621f6655bed09ad72dd0eabbae35bd3c47be483b558cdeb593d48`, F− `1dba5ba3ff8f3790a8ef3231e5da02380163e4672b81c1a40a122e5dd92b711d`). Originalni i neovisni rezultati sada su trajno u `source_data/` iste audit grane, ne samo u privremenim GitHub Actions artefaktima.

**Fizički status nakon numeričkog PASS-a:** E17A i E17B linearni CLASS oba kraka DONE; E17 puni finite-K retarded Einstein–Vlasov halo/tracer response **BLOCKED**. Poznate realne linearne `P_cb`, `delta_cb` i direct-vTk vrijednosti ne određuju dinamičke halo/tracer spregne koeficijente. Nismo proizveli niti prijavili cijeli galaktički `B`, eBOSS 3pt prozor, `xi`, S/N ili značajnost. Opaženi odd SEALED, A03/A04 zasebno otvoreni, izvorni 24D ostaje dvotočkasta statistika, `main` netaknut i PR #1 draft.

