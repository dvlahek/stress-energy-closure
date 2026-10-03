# A-03E7 — što smo stvarno naučili iz dovršenog full-galaxy 4800R/48000R mocka

**Status, 27. 9. 2026.:** retrospektivna analiza već dostavljenog rezultata. **Ne** novi unaprijed prijavljen fizikalni acceptance, ne procjena slučajne pogreške 48 000 randoma, ne dokaz valjanosti eBOSS survey-window operatora i ne Einstein–Vlasov detekcija. Nema novih FITS parova, random drawova, seedova, rezova, mock preuzimanja ili opaženog odd pristupa. Izvorni `main` i draft PR #1 nisu dirani.

[Izvorni korisnički full-mock završetak i 32 pair-SHA/4 LS/16 multipole/26-slice audit](../source_data/eboss_dr16_mock0001_full_galaxy_4800_48000_completed_source_only_audit_2026-09-27.json) zaključava izvještaj od **529771 bajtova** s SHA256 `187e21e206cdb6ded1108887a2b6670113d50b752e438d28212cfc3d187572da`. [Kompaktna numerička dijagnostika](../source_data/eboss_dr16_a03_e7_mock0001_ls_cancellation_source_only_summary_2026-09-27.json) i [ponovljiv source-only Python audit](../scripts/audit_eboss_dr16_a03_e7_fullmock_pair_conditioning.py) ne ovise o opaženom katalogu. Audit ima izričitu zaštitu izvornog SHA, oba kapa, cijeloga DD uzorka, 144/144 RR podrške i fizikalnog `mu` clipa. Sintetički zatvara sve 3! = 6 redoslijeda zamjena, provjerava paritet/Legendre projekciju i odbija promjenu DD. [CI `36306590429` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36306590429) je **samo synthetic/source-only**; stvarni rezultat ovog dokumenta neovisno je izračunat iz korisničkog JSON-a, ne iz originalnih FITS parova.

## Uzrok velikog relativnog pomaka `xi` nije sama veličina promjene parova

Neka su `D,U,V,R` četiri originalna **zasebno normalizirana** weighted `DD,DR,RD,RR` histograma u svakoj `(s,mu)` ćeliji. Za fiksni puni DD uzorak vrijedi

\[
\xi=1+\frac{D-U-V}{R}.
\]

Promjena `4800R→48000R` u istom EZmock0001 ulazu mijenja `U,V,R`, dok je `D` bitwise jednak. Redoslijedno neovisan *algebarski* Shapleyjev doprinos `U,V,R` definiran je kao aritmetička sredina teleskopskih promjena `xi` preko svih šest permutacija zamjene iz rjeđe u gušću realizaciju. Njihov zbroj egzaktno rekonstruira `delta xi` do numeričke tolerancije `2e-13`. Ovo je dijeljenje međudjelovanja koje nastaje zbog nazivnika; **nije** statistički bootstrap, mjera uzročnosti ili zasebna varijanca DR/RD/RR.

Za 48k definiramo samo lokalni, deskriptivni omjer poništavanja

\[
\kappa_{\rm cancel}=
\frac{|D|+|U|+|V|+|R|}{|D-U-V+R|}.
\]

To **nije** kondicijski broj inferencijske kovarijancijske matrice i ne daje strogu gornju ogradu slučajne pogreške. Njegova velika vrijednost pokazuje zašto mali pomaci pojedinačnih normaliziranih brojanja mogu proizvesti velike *relativne* promjene konačnog `xi` blizu nule.

| Već izmjerena source-only dijagnostika | NGC | SGC |
|---|---:|---:|
| Relativni L1 promjene pune `xi` mreže, 4800→48000 prema 48k | 1.674993 | 1.460288 |
| Srednji `|xi|` na 48000R (144 ćelije) | 0.026404 | 0.024530 |
| Srednji `|delta xi|` (144 ćelije) | 0.044226 | 0.035821 |
| Relativni L1 DR histograma | 0.02434 | 0.01766 |
| Relativni L1 RD histograma | 0.04173 | 0.01205 |
| Relativni L1 RR histograma | 0.04346 | 0.02092 |
| Medijan `kappa_cancel` u 144 ćelije | 242.21 | 237.35 |
| P90 `kappa_cancel` | 980.41 | 1166.37 |
| Najveći šest-bin `|delta xi_1|` | 0.025465 | 0.051958 |
| Najveći šest-bin `|delta xi_3|` | 0.039177 | 0.047268 |

Nisu dopuštene dvije pogrešne interpretacije: veliki relativni `L1` **ne dokazuje** da ni 48000R ne konvergira, ali 4800R **ne može služiti kao stabilna referenca ovom mocku** bez fizikalno opravdanog tolerancijskog budžeta. Smjer i veličina `xi_1,xi_3` variraju po izvornim šest `s` binova i kapi, bez univerzalne pozitivne EV signature. SGC brojanje reverse velikih parova u V2.1 je algebarski mirror, a izvorni neovisni v1 SGC `R1D2` mismatch ostaje nerazjašnjen.

## Pravi fizički budžet nije moguće izvući iz ovog mocka

U [A-03E6 teorijskom bridgeu](EBOSS_DR16_A03_E6_THEORY_TO_EMPIRICAL_WINDOW_BRIDGE_2026-09-27.md) postojeći `code/build_lrg_elg_wake_template.py` računa **shape-only** dipol i potom zasebno maksimalno normalizira `wake` i `dop`. Njegov fitted amplitude nije predviđena dimenzijski normalizirana eBOSS `xi_1(s,z)`. `code/build_lrg_elg_physical_odd_basis.py` računa pred-prozorske Hankel integrale i uvjetnu standardnu relativističku kombinaciju, ali opciono zahtijeva zasebno kalibrirane `b_LRG,b_ELG,s_LRG,s_ELG,fevo_LRG,fevo_ELG`, uz DESI ideal-shell `s²ds` biniranje. Ne smijemo DESI oblik, odabrani `frac`, free amplitude, ili `z=0.95` nazvati eBOSS fizičkom prognozom.

Minimalni, **još nedostajući** redoslijed za fizički A-03 je: najprije source-preserving apsolutni EV teorijski odziv `xi_L(s,z)` s jasno dokumentiranom vezom gravitacijskoga stanja i LRG/ELG odgovora; zasebno kalibrirani standardni odd i `L=0,2,4` even ulazi; stvarno uvjetovani tracer-specific RR operator po kapi i paru u `z` uz poznate ELG chunk/depth ovisnosti **samo u granici koju objavljeni randomi stvarno kodiraju**; zatim validacija fizičkog even→odd leakaga. Originalni E6 `6×120` blokovi marginaliziraju `z` i sintetički sučeljavaju `L=0,...,4`, ali još nisu validirani realni `M(s,z)`.

Tek tada ima smisla registrirati konačni `p`-dimenzijski eBOSS observable, neku neovisnu `C_{\rm eBOSS}` te fizički bias kriterij `|\delta A|/\sigma_A`. Devet postojećih EZmock ID-jeva daju rang centrirane kovarijance **najviše osam** za zajednički 24D vektor, a DESI 120-mock 18D kovarijanca je drugi uzorak i drugi observable. Nema numeričkog praga `epsilon_A` koji bi se legitimno odabrao *nakon* pregledanog 4800→48000 odstupanja.

**Operativni rezultat:** skupo ponovno brojanje istih parova zamrznuti. Nove random seedove, fizičke rezove, prerano čitanje opaženog odd vektora i `main` promjene ne autorizirati. Sljedeća znanstvena aktivnost treba rješavati fizičku normalizaciju i pair-redshift uvjetovanje, ne optimirati post-hoc veličinu već pregledanog mock efekta.
