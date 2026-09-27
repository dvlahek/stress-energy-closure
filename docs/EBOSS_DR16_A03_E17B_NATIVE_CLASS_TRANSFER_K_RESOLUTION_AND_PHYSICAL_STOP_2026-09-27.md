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

