# E17D2b3b7R1 — postojeći numerički rezultat, bez ponovnoga WSL računanja

**29. 9. 2026. · Korisnički dostavljen JSON pregledan offline; Ubuntu crash cause UNDETERMINED.** Ovo je post-run izvještaj, NIJE retroaktivna izmjena E17D2b3b7 preregistracije ili originalnog numeričkog čitača.

## Već sačuvani lokalni rezultat

Korisnik je poslao mali lokalni `halo_info_000_id_N_locked_rows.json` (1593 B, SHA256 `79025a040a583d461fb11e45ee2a747e13553c2036e046f2d48ed42f752b918f`) te `e17d2b3b7_full.log` (1389 B, SHA256 `c04beda0290f4ed271b072b5f16a30f445053a189db5b250561dbce603cde097`). JSON je učitan kao valjan JSON i fiksna tri indeksa te `N × 2109081520.453063 Msun/h` aritmetika su neovisno provjereni **na prenesenom malom JSON-u**, ne na izvornoj 2.2 GB ASDF datoteci. `N` statistika samo za superslab 000: `rows=11676687`, `min=35`, `max=324518`, `sum_N=2497001263`, `zero_N=0`. Fiksna tri indeksna retka: `(0,id=22001000,N=53,M_L1=111781320584.01234 Msun/h)`, `(5838343,id=2416060793000000,N=61,M_L1=128653972747.63684 Msun/h)`, `(11676686,id=4917001651002000,N=44,M_L1=92799586899.93477 Msun/h)`. Ovaj rezultat je **USER-UPLOADED/USER-LOCAL**, ne neovisno GitHub CI izvršenje stvarnog ASDF-a. Cjeloviti korisnikov lokalni JSON i log nisu kopirani u javni GitHub, samo ovaj audit sažetak i njihove hash vrijednosti.

## Što log zapravo pokazuje

Poslani log javlja osam sintetičkih `NEGATIVE_REJECT` kontrola PASS i pri sljedećem pokretanju `ValueError: E17D2B3B7_FAIL_CLOSED: do not overwrite any earlier local output` u `input_preflight`. To znači da je novi pokušaj pravilno odbio **postojeću izlaznu datoteku** prije numeričkog čitanja; ne dokazuje pad WSL kernela, neuspjeh originalnog prvog čitanja, korupciju ASDF-a ni OOM. U dostavljenom logu NEMA originalnog uspješnog stdout markera ni povijesti Windows/WSL kernel događaja; uzrok korisnikova rušenja Ubuntua ostaje nepoznat. NIKAD ne brisati/brisati rezultatski JSON da bi prošao rerun. Originalni fail-closed čitač i B7 prereg ostaju nepromijenjeni.

## Operativna odluka

**STOP: nema novog B7 stvarnog numeričkog pokretanja.** Sačuvati postojeći lokalni JSON i log. Ako Ubuntu trenutno ne odgovara, iz Windows PowerShell terminala prvo očitati samo `wsl --list --verbose` i `wsl --status`; tek ako je WSL i dalje zaglavljen i nema drugih poslova, korisnik može izvršiti `wsl --shutdown` (ovo prekida SVE WSL distribucije/procese) te ponovno otvoriti Ubuntu. Za dokaz uzroka pada potreban je zaseban stvarni Windows/WSL/OOM trag; ne nagađati iz Python file-exists iznimke. Bez novih ASDF transfera, otvaranja `SO_radius`, drugih halo datoteka, čestičnih podataka ili eBOSS opažanja.

U korisničkom rezultatu `min_N=35`, dok službena zaključana opisna dokumentacija kaže `more than 35 particles`; tretirati ovo samo kao izvorno-dokumentacijsku granicu za budući **source-only** audit, bez retroaktivnog filtriranja ili novih science cuts. L1 assigned `N*m_p` NIJE `M200c`; nema merger tree ni fizičkog neutrinskog wakea/bispektra. Opaženi 24D odd SEALED, originalni E8/E16/48k zamrznuti, `main` se ne dira, PR ostaje draft.
