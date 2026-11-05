# eBOSS DR16 A-02 — devet uparenih EZmock galaxy×random cross-LS realizacija

**Status 26. 9. 2026.:** A-01 zatvoren i arhiviran; A-02 protokol preregistriran prije prvog čitanja mock-galaxy FITS redaka drugih osam unaprijed odabranih ID-jeva; sintetički/source-only CI prošao; **lokalni A-02 parovi još nisu računani**. Glavni kanonski projektni plan: `EINSTEIN_VLASOV_NP_MASTER_PLAN.md` v1.1. Sve promjene samo u draft PR #1, bez promjene `main` i bez otvaranja opaženog odd signala.

## Zašto A-01 sada ima zatvoren status

Četiri stare `0001` mock-galaxy gzip datoteke i 32 dodatne iz unaprijed izabranih `0001,0125,0250,0375,0500,0625,0750,0875,1000` već su potvrđene punim SHA256 izvorišnih bajtova. Arhiva i odvojeni manifest: `source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_report_2026-09-26.json` i `source_data/eboss_dr16_nine_ezmock_galaxy_source_sha_uploaded_manifest_2026-09-26.json`.

Korisnikov završni mock-random report ima točno **45 653 bajta**, SHA256 `229a107a44fa1e8665c39acc9050982f8a0fcc39d15cc35e5bbc9bc272116c47`, Git blob `315b6344d9ea2a2d5f6b877fd8cc9299e9b6e24c`. Bajtno identično arhiviran u `source_data/eboss_dr16_nine_ezmock_matched_random_sha_report_2026-09-26.json`; odvojeni, prije mock redaka zaključani manifest nalazi se u `source_data/eboss_dr16_nine_ezmock_matched_random_sha_uploaded_manifest_2026-09-26.json`. Svih **36/36** ranijih punih random SHA256 referenci iz uspješnog originalnog audita od 24. 9. potvrđeno je na lokalnim komprimiranim bajtovima; 4 iz cachea, 32 sa službenog ranije fiksiranog URL-a. Ukupno 1 786 953 149 komprimiranih random bajtova. U A-01 nisu otvoreni FITS redci, parovi ni opaženi odd vektor.

## Što A-02 smije napraviti

**Protokol:** `source_data/eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport_protocol_2026-09-26.json`. **Runner:** `scripts/audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport.py`.

Prije ijednog novog mock FITS headera/retka uspoređuje točne lokalne A-01 galaxy i random JSON-ove s bajtno identičnim Git arhivama, provjerava zasebne manifeste i ranije reference te *ponovno hešira svih 72 potpunih gzip datoteka*: 36 mock galaxies i 36 same-ID/same-cap/tracer-specific randoms. Ako jedan izvor nedostaje, promijeni se, otisak ne odgovara, nije u ranije odobrenoj lokalnoj mapi ili je symlink: stop bez čitanja mock FITS redaka. Runner **ne preuzima** niti zamjenjuje datoteke i ne prihvaća nova prva SHA otiske.

Zatim reproducira `0001` originalni `600D/1200R` uzorak, njegovih osam izvornih odabranih array SHA256 i sve izvorne forward/reversed weighted-pair histogram SHA256 za oba capa. Tek nakon toga radi identičan algoritamski test na ostalih osam unaprijed izabranih ID-jeva. NGC/SGC zasebno, LRG/ELG s vlastitim objavljenim randomima; bez nove zajedničke tvrde maske.

Svi ID-jevi imaju isti ranije fiksirani z-interval `[0.9,1.0)`, deterministički no-replacement `600D/1200R` po traceru i capu, korijen sjemenke `93127` i izvornu kapu/tracer/role formulu, midpoint LOS, `theta_min=0.05°`, `s=[20,40,60,80,100,120,140] h^-1 Mpc`, 24 potpisane `mu` ćelije, `WEIGHT_SYSTOT*WEIGHT_CP*WEIGHT_NOZ*WEIGHT_FKP` i originalnu konvenciju numeričke nule. Čuvaju se originalne ELG chunk oznake, bez chunk-conditioned resamplinga.

Za svaki fiksni ID/cap izvodi četiri nezavisno normalizirana cross-LS člana, fizički zasebnu obrnutu LRG↔ELG orijentaciju, weighted pair-mirror closure, signed-`mu` ξ-mirror i paritet sirovog RR odd komponente. Pohranjuje puni pregled RR podrške i popis nepodržanih ćelija bez zero-filla. Ako fiksni uzorak ili slučaj padne, bilježi se pogreška i ne mijenja se sample count, seed, bin ili popis ID-jeva. Točan završni JSON te svi checkpointi ostaju na disku.

## Što A-02 ne smije zaključiti

A-02 je **code transport**, ne inferencijska procjena raspodjele mock odd multipola i ne physical null test. Ranije `0001` nenulte `ell=1,3` vrijednosti ne služe odabiru realization ID-ja ili analitičke konfiguracije. Centrirana kovarijanca devet mockova ima rang najviše osam, pa nije punorangovna 18D kovarijanca i ne daje poštenu “veliku sigmu”. Fizički potpuni ELG/LRG produkcijski prozor ostaje zaseban A-03 gate; opaženi odd vektor ostaje SEALED.

## Sintetički test i lokalna naredba

[Source-manifest/synthetic CI run 36244602728](https://github.com/dvlahek/stress-energy-closure/actions/runs/36244602728) prošao je bez pristupa stvarnim FITS datotekama. Provjerava dvije završne A-01 Git arhive, mapiranje svih 72 izvora, originalni `0001` parent report i sintetički negativni slučaj izmijenjenog starog random SHA; ne računa mock galaksije u CI-u.

**Jedna WSL naredba, kad želiš krenuti s A-02:**

`cd ~/stress-energy-closure && git pull --ff-only && source .venv/bin/activate && python -u scripts/audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport.py --self-test && python -u scripts/audit_eboss_dr16_nine_ezmock_galaxy_cross_ls_code_transport.py`

Izlaz: `eboss_workspace/local_ninemock_galaxy_cross_ls/nine_mock_galaxy_cross_ls_code_transport.json`. Pošalji neizmijenjeni završni JSON **i** završni ispis; ako program stane, pošalji isti JSON s `INCOMPLETE_STOP` i konkretnom greškom. Svaki uspješan ID/cap ima atomarno checkpointiran status. Pri idućem pozivu najprije provjeriti taj izvještaj, ne pretpostaviti da je test uspio samo zato što je CI zelen.
