# EinsteinVlasovNP A-03E2 — blind mock-galaxy ELG RANDOM radial-weight stress

**Status (26. 9. 2026.): A-03E2 ZAVRŠEN 18/18.** Protokol je bio zaključan prije lokalnog testa, svih 72 punih originalnih mock gzip SHA provjereno je prije E2 FITS redaka, svi originalni A-02 sample/pair/ξ SHA reproducirani, svih 54 scenarija imaju 144/144 RR ćelija i tracer-reversal closure. Točan izvorni izvještaj i zasebni manifest sada su bajtno arhivirani i zaštićeni source-only CI auditom. **A-03 empirijski fizički prozor i A-04 kovarijanca nisu time certificirani.** Autore ne kontaktirati; ne otvarati opaženi odd vektor.

## Što je već potvrđeno

- A-01: 36 realističnih mock-galaxy i 36 same-realization/cap/tracer random punih gzip SHA256 otisaka.
- A-02: svih devet fiksnih mock ID-jeva \`0001,0125,0250,0375,0500,0625,0750,0875,1000\` × NGC/SGC prošlo je stvarni cross-LS code-transport, originalni \`0001\` reprodukciju, sve \`144/144\` RR ćelije i tracer reversal. Izvještaj je zamrznut pod punim SHA256 \`15f7668fd483d8e1329fbb9684bbdf07cb49d1974d264f9ceebd85e0156daa27\`.
- A-03 E0/E1: devet-mock released-RANDOM-only uvjetni finite-bin operator i unaprijed zadana sintetička even/odd injekcija prošli su, ukupno 18/18, s arhiviranim izvještajem i odvojenim manifestom. Sintetički radijalni even ramp nije fizički izmjeren wake ili stvarna sistematika.

## Što točno radi E2

**Preregistrirani protokol:** \`source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_protocol_2026-09-26.json\`.  
**Runner:** \`scripts/audit_eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress.py\`.  
**Source-only/synthetic CI:** https://github.com/dvlahek/stress-energy-closure/actions/runs/36266710895 — **success**. CI nije pristupio lokalnim pravim mock FITS datotekama i nije izračunao E2 stvarne mock slučajeve.

Runner najprije provjerava originalni A-02 report, njegov izvorni protokol, nezamijenjeni kod i arhiviranu odluku o empirijskom putu. U WSL-u zatim ponovno hešira *svih 72* kompletnih izvornih gzip datoteka **prije ikojeg FITS headera ili retka**. Nijedan originalni otisak, ID ni uzorak ne smije se zamijeniti novim. Za svaki ID i kap ponovno deterministički odabire prethodnih \`600D/1200R\` po traceru, provjerava izvorne sample+chunk dijagnostike, ponovno računa sva četiri LS člana u obje orijentacije i uspoređuje svih osam originalnih weighted-pair histogram SHA i forward \`xi\` SHA iz A-02.

Tek nakon točne reprodukcije izvorne realizacije izračunava dva **namjerno pogrešna mock-random modela**, bez novih maski i bez mijenjanja mock galaksija:

\[
u(z)=2\frac{z-0.9}{0.1}-1,\qquad
w^{(\pm)}_{{\rm ELG},R}(z)=w_{{\rm ELG},R}(z)[1\pm0.05\,u(z)].
\]

Fiksni interval ostaje \([0.9,1.0)\), faktor je ograničen na \([0.95,1.05]\), mijenjaju se *isključivo* već odabrane ELG RANDOM težine. LRG randomi, sve galaksije, sve pozicije, ELG chunkovi, redshift i binovi ostaju identični. Svaki član i svaka orijentacija ponovno dobivaju **vlastitu** normalizaciju \(\sum w_1\sum w_2\). Očekivano nepromijenjeni članovi i broj prihvaćenih parova moraju se točno podudarati; originalni i perturbarani \(\xi(s,\mu)\) moraju proći tracer reversal.

Zatim projektira \(\ell=0,1,2,3\) na šest **istih** \(s\)-binova s fizički ograničenim \(\mu\in[-1,1]\) i bilježi opisne razlike \(\Delta\xi_\ell^{\pm}\) za svih 18 slučajeva. Ne postoji minimalni prag za „dobar” efekt i ne biramo povoljan znak, bin, kapu ili mock ID. Bilo koja ćelija bez RR potpore prekida projekciju tog slučaja uz poštenu evidenciju greške; bez zero-filla ili re-selekcije.

**Izlazni lokalni JSON:** \`eboss_workspace/a03_empirical_window/a03_e2_mock_galaxy_elg_random_radial_weight_stress.json\`. Izvještaj se atomski ažurira po slučaju. Ako job stane, sačuvati izvorni \`INCOMPLETE_STOP\` JSON i završni ispis; ne pokretati novu varijantu s promijenjenim seedovima ili ID-jevima. Uspješni checkpointovi mogu se ponovno upotrijebiti tek nakon ponovne verifikacije svih 72 izvora.

## Završni originalni rezultati, odvojeno od preregistracije

**Neizmijenjeni korisnički izvještaj:** `source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_report_2026-09-26.json`, **413 189 bajta**, SHA256 `6ef86f904bb0cdc407adabe2e125146ce4d0671ea1ac932e13d11fc481cd596a`, Git blob `8596d7322f2ec7ef164e89641958c3912f687722`. Binarno/bajtno identičan originalnom lokalnom prijenosu. **Manifest:** `source_data/eboss_dr16_a03_e2_mock_galaxy_elg_random_radial_weight_stress_uploaded_manifest_2026-09-26.json`, Git blob `7ec2536429604227b2428637ae227b53b7d42cb4`. Neovisni audit izvornog A-02 roditeljskog izvještaja i E2 sva 54 slučaja, uključujući negativne provjere mogućeg izmijenjenog SHA ili otvaranja observed odd, [prošao je u CI-u](https://github.com/dvlahek/stress-energy-closure/actions/runs/36269865486). Stvarni E2 lokalan je korisnikov rad; CI samo reaudita arhivirani izvještaj, ne ponovno računa FITS parove.

**Svi fiksni ulazi i slučajevi:** 9 ID-jeva × NGC/SGC = 18/18, tri scenarija svaki, 0 prijavljenih grešaka, svih `144/144` RR-podržanih ćelija; maksimalni `|ξ_{m LRG→ELG}(s,μ)−ξ_{m ELG→LRG}(s,−μ)| = 1.7763568394002505×10⁻¹⁵`. Najveća stvarna primijenjena random-težinska promjena ostaje manje od `±0.05`, a originalni A-02 podaci i svi odabrani redci ostaju nepromijenjeni.

**Deskriptivni odgovor na umjetnu +5 %-rampu** (po 9 mock ID-jeva zasebno u svakoj kapi, statistika je medijan `max_s |Δξ_ℓ(s)|`, a u zagradi najveća zabilježena vrijednost po tim ID-jevima):

| Kapa | Dipol ℓ=1 | Oktupol ℓ=3 |
|---|---:|---:|
| NGC | `0.00993128` (`0.0209470`) | `0.0188281` (`0.0389065`) |
| SGC | `0.00568632` (`0.0101764`) | `0.00479366` (`0.0143389`) |

Kod −5 % rampe medijani su NGC `0.00953158` (ℓ=1), `0.0184347` (ℓ=3), SGC `0.00556621` (ℓ=1), `0.00482831` (ℓ=3). Brojevi su **apsolutne promjene bezdimenzijskog mock ξ**, nisu postoci, σ-vrijednosti niti izmjerena stvarna selekcijska pogreška. Razlike u uzorkovanju, konačnoj gustoći randoma i samom obliku proizvoljnog probnog `u(z)` nisu opservacijski kalibrirane.

**Sljedeći otvoreni zadatak:** zasebno zaključati A-03E3 empirijski blind mock-only injection/recovery i provjeru konvergencije na gustoću randoma, osobito odvojiti stvarnu osjetljivost na uvjetni prozor od finite-random i sparse-pair šuma. Tek nakon toga planirati A-04 produkcijski 18D opservabilni vektor i odgovarajući neovisni ansambl. Nijedan novi observed odd pristup nije odobren.

## Što to nije

E2 mjeri **deskriptivnu osjetljivost mock cross-LS estimatora na umjetno uvedenu pogrešku u radijalnim težinama randoma**. To nije kalibrirana stvarna nesigurnost released random selekcije, nije dokaz kompletne povijesne maske, nije statistički neovisan galaxy-based physical null ansambl za 18D kovarijancu, nije \(S/N\), p-vrijednost, wake detekcija ni exclusion limit. E2 ne čita nijedan redak opaženih galaksija ili opaženog odd vektora. A-03 ostaje otvoren za valjanu empirijsku physical-window/injection-recovery kalibraciju, a A-04 za zasebnu inferencijsku kovarijancu.

**Trenutni sljedeći korak:** E2 izvještaj i manifest su već zaprimljeni, neovisno verificirani i arhivirani. Ne pokretati ponovni E2 job radi lovljenja efekta; sljedeći zadatak je preregistrirati E3 random-density/split/injection-recovery pilot prije novog lokalnog rada. Ne preuzimati nove kataloge i ne koristiti observed data.
