# A-03E17C — zamrznuti kinetički transport na oba kraka i fizički STOP

**27. 9. 2026. · Originalni i neovisni source-only račun PASS; stvarni retardirani halo/tracer finite-K bispektar BLOCKED.** E17C slijedi E17A/E17B. Ne pokreće ponovno E8–E16, ne mijenja zamrznute distribucije i ne otvara opaženi eBOSS odd.

## Problem i zaključani protokol

E17A/E17B izračunali su izvorni linearni CLASS P_cb, direct-vTk i delta_cb na oba stvarna kratka kraka svih 576 E16 zatvorenih trokuta. Time nije određena fizička funkcija odziva od dugovalne gustoće/vjetra preko halo neutrinske perturbacije do dvaju LRG/ELG kratkih krakova. Izvorni statički Eq20 E8b/E9 ne smije se nazvati dinamičkim Vlasovljevim halo odzivom. E17C zato ne dodaje proizvoljnu finite-K simetrizaciju originalnom E14/E15 reduciranom izvoru.

[Prospektivni E17C protokol](../source_data/eboss_dr16_a03_e17c_restricted_retarded_kinetic_propagator_prereg_2026-09-27.json), Git blob 33defef2020f2931756e3ef7a3959d905d539165, commit 4274d282, zaključen je prije numeričkog računa. SHA-pina originalni E8 4000q CSV, E14/E15, E16 originalni i neovisni JSON, E17A/E17B originalne joint izvještaje, svih 12×48 originalnih E16 geometrija, tri originalna F stanja i redoslijed LRG→ELG. Nema novog CLASS-a, FITS/mock preuzimanja, seedova, rezova ili reoptimizacije neutrinskih distribucija.

## Točno izračunani ograničeni objekt

Za izotropnu homogenu izvornu distribuciju, slobodan bezsudarni transport između vanjskih perturbacija u fiksnoepohnoj vodećoj nerelativističkoj aproksimaciji daje

\[
(\partial_\eta+i\mathbf{k}\cdot\mathbf v(q))\delta f=0,\qquad
U(\eta,\eta';\mathbf{k},\mathbf q)=\Theta(\eta-\eta')\exp[-i\mathbf{k}\cdot\mathbf v(q)(\eta-\eta')].
\]

Ovo je **homogeni slobodni propagator**, ne gravitacijski force/source član, samokonzistentno Einstein–Vlasov rješenje, relativistička transportna hijerarhija ni galaktički bispektar. Bez halo potencijala i njegova vremenskog razvoja ne može se provesti puna retadirana konvolucija. Numerički E17C računa samo za nenegativni slobodni lag i ne tvrdi da je posebno testirao negativne fizikalne vremenske lagove.

Na izvornom konačnom q rasponu po stanju računa se

\[
C_F(s)=\frac{\int_{q_{\min}}^{q_{\max}}\!dq\,q^2F(q)\operatorname{sinc}(sq)}
{\int_{q_{\min}}^{q_{\max}}\!dq\,q^2F(q)}.
\]

Za oba stvarna kraka vrijedi s_i=s_0 |k_i|/(0.05 h/Mpc), uz unaprijed izabrane isključivo dijagnostičke bezdimenzijske vrijednosti s_0=(0,0.25,1). One **nisu** stvarne vremenske udaljenosti, fizikalni halo-wind, novi odabir redshifta ili empirijski fit. Izvorni q raspon nije ekstrapoliran. Originalni [izvršivi kod](../scripts/audit_eboss_dr16_a03_e17c_restricted_retarded_streaming.py) daje svih 3×576×2×3 = **10.368** originalnih state/leg/phase vrijednosti, s ne-normaliziranim integralima, normaliziranim C_F i QA. Sredina originalnih F+ i F− na izvornom q gridu jednaka je F0 do skaliranog 1.0757463211989661e-16. To ne implicira jednakost zasebno normaliziranih karakterističnih funkcija niti detekciju neutrinskog wakea.

## CI, popravci i neovisni rezultat

Prvi [CI 36346000174](https://github.com/dvlahek/stress-energy-closure/actions/runs/36346000174) ostaje FAIL: nova provjera je tražila nepostojeću E16 -0.6 kratku orijentaciju. Originalna mreža je (-1,0,0.6,1). Tehnički popravak koristi egzaktni analitički reflektirani krak samo za nedostajuću kontrolu, bez dodavanja novih science orijentacija. Drugi [CI 36346041734](https://github.com/dvlahek/stress-energy-closure/actions/runs/36346041734) izračunao je sva tri izvorna reporta, ali neovisni audit je pokušao otvoriti direktorij FD umjesto izvornog JSON imena; ostaje FAIL. Odvojeni ispravak samo čitačkog puta ne dira fiziku, izvorne SHA roditelje ni pragove.

Završni [CI 36346079827](https://github.com/dvlahek/stress-energy-closure/actions/runs/36346079827) **SUCCESS**. Sva tri originalna stanja imaju 576 zatvorenih trokuta i 3456 vrijednosti po stanju. Najveći E16 closure ostatak 1.2643861424099487e-17 h/Mpc, nulta faza i Fourierov povratni smjer daju odstupanje nula, a kratko-kraku identitet 2.7755575615628914e-17. Riječ je o matematičkim kontrolama ograničenoga propagatora, ne o bispektralnoj pogrešci.

Zasebni [neovisni standard-library audit](../scripts/audit_eboss_dr16_a03_e17c_independent_scalar_replay.py), bez uvoza originalnog E17C koda, NumPyja ili CLASS-a, iz izvornog q CSV-a nanovo gradi skalarne trapezne težine i analitičku E16 geometriju te provjerava svih 10.368 vrijednosti. Najveće skalirano odstupanje 3.3306690738754696e-16; negativna kontrola izmijenjenoga SHA otiska PASS. To je nezavisna numerička rekonstrukcija *istih* izvornih podataka, ne drugi fizički halo solver.

[Trajni manifest](../source_data/eboss_dr16_a03_e17c_archived_CI_2026_09_27/archive_manifest.json), Git blob 341b2bf174920254aa867d433b383bf3df6cbd5a, arhivira tri puna izvorna state JSON-a, izvorni joint SHA256 ad655ada4f64e5033fe900d5821f58847e3c96e72e081cd1a92eadcf9a7236f4 i neovisni certifikat SHA256 3dacec420697c5989db7e043fccf09b59bf6017951c9a24dfe20f87aa5dcd841. Svi su izvorni rezultati trajno u source_data/ na audit grani, ne samo u kratkotrajnom Actions artefaktu.

## Fizička granica

Za stvarni finite-K retardirani halo odziv potrebni su unaprijed specificirani potencijal/masena povijest haloa, vremenski razvoj u konzistentnom gaugeu, gravitacijski source, dinamička Vlasovljeva reakcija te fizički izveden long–short coupling **obaju** kratkih halo/tracer krakova. E17C nije izveo te podatke iz CLASS linearnih transfera. Dodatno su potrebni nezavisno kalibrirani high-z LRG/ELG bias/HOD, evolution/magnification/relativistički i selekcijski nuisancei te stvarni eBOSS fizički 3pt prozor i neovisna 3pt kovarijanca. Originalni 24D odd dvotočkasti vektor nije bispektar. A03 pair-z RR, A04 kovarijanca, originalni SGC reverse i 48k finite-random konvergencija ostaju otvoreni.

**Stvarni status:** E17A/B linearni oba kraka DONE; E17C ograničeni source-only transport i neovisni audit DONE; puni retardirani halo/tracer response i fizički B BLOCKED. Nema eBOSS xi, S/N, značajnosti ni unblindanja. Opaženi odd SEALED, main netaknut i PR #1 ostaje draft.
