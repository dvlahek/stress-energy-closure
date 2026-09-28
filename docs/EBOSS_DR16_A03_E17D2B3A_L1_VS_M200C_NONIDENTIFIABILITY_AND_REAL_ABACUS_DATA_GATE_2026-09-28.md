# A-03E17D2b3a — zašto L1 masa i radijus ne određuju M200c bez radijalne raspodjele

**28. 9. 2026. · Izvorni i neovisni matematički audit PASS, [CI 36421065060](https://github.com/dvlahek/stress-energy-closure/actions/runs/36421065060). Stvarni AbacusSummit ASDF/header, halo i merger-tree podaci NISU pročitani niti preuzeti. Opaženi odd SEALED i fizički LRG×ELG bispektar BLOCKED.**

## Problem

E17D2b2 je iz službene dokumentacije utvrdio da je AbacusSummit c000 *potencijalan* javni izvor merger trees s kozmologijom nominalno bliskom izvornom CLASS-u. Međutim, [službeni model proizvoda](https://abacussummit.readthedocs.io/en/latest/data-products.html) koristi CompaSO L1 halo grupu i povezane epohno definirane SO pragove gustoće u odnosu na srednju gustoću. Naš dosadašnji E17D2a/E17D2b0 NFW model zahtijeva M200c i r200c, tj. prosječnu gustoću mase unutar sfere **200 puta kritičnu** gustoću u stvarnoj epohi. Prije validiranoga spajanja tih masenih definicija treba razjasniti koje su informacije nužne za prijelaz s L1 na M200c.

E17D2b3a ispituje jači, idealizirani slučaj: čak ako bi L1 masa i radijus te referentni prag srednje gustoće bili poznati savršeno i SO radijus stvarno zadovoljavao jednadžbu praga, može li se M200c jednoznačno odrediti bez unutarnjega profila? **Ne može**, što izravno pokazuju dvije točno konstruirane pozitivne sferne raspodjele. To nije dokaz da sve ostale CompaSO halo_info varijable ili pune pozicije čestica ne bi pomogle. To također nije mjerenje stvarnih Abacus halo masa.

[Prospektivni protokol E17D2b3a](../source_data/eboss_dr16_a03_e17d2b3a_mass_reference_nonidentifiability_prereg_2026-09-28.json), Git blob 81554b0f7cdefa599f302d828ab3529812bb3a59, commit 00dec005, zaključan je **prije novoga izračuna**. SHA-zaključava originalne E8/E16/E17A/E17D2b1, E17D2b2 originalni i neovisni izvještaj te E17D2b2 manifest/protokol/runner. Nije promijenjen originalni skup 576 trokuta, F0/F+/F−, E14/E15 72 source/24 contrast, opažena odd zabrana niti originalni CLASS.

## Ograničeni analitički kontraprimjer

Radi reproduktivnosti se **samo u matematičkom testu** uzimaju javno dokumentirane nominalne c000 konstante h=.6736, ωb=.02237, ωcdm=.1200, ων=.00064420, z=.95 i ilustrativni SO prag \(\Delta_{\mathrm{L1,mean}}=200\). Vrijednost z=.95 jest *testna oznaka*, ne očitano stvarno ASDF zaglavlje. Vrijednost praga 200 jest matematički odabir, ne tvrdnja o stvarnom SODensityL1 zaglavlju. Srednja referentna gustoća definirana je kao total matter u jednostavnoj flat matter+Λ algebri i **nije validirana stvarna referentna gustoća CompaSO**. Zanemarena je radijacija i detaljna termodinamika masivnih neutrina; ovo nije CLASS pozadina.

\[
\Omega_{m0}^{\rm toy}=\frac{\omega_b+\omega_{\rm cdm}+\omega_\nu}{h^2}
=0.3151918679932973,
\qquad
\Omega_m^{\rm toy}(z{=}0.95)=0.7733861452097148.
\]

Budući da \(\rho_m(z)=\Omega_m(z)\rho_{\rm crit}(z)\), u **ovoj** sintetičkoj konvenciji ciljni 200-kritični prag prema srednjoj gustoći iznosi

\[
\Delta_{200c,\rm mean}^{\rm toy}
=\frac{200}{\Omega_m^{\rm toy}(z)}
=258.6030293389431
> \Delta_{\rm L1,mean}^{\rm toy}=200.
\]

Postavimo identične normirane \(R_{\rm L1}=1\) i \(M_{\rm L1}=1\). Za dva različita pozitivna modela s \(p\in\{1,2\}\) definirajmo

\[
M_p(<r)=M_{\rm L1}\left(\frac{r}{R_{\rm L1}}\right)^p
\quad (0\le r\le R_{\rm L1}),
\qquad
M_p(<r)=M_{\rm L1}\quad(r\ge R_{\rm L1}),
\]

\[
\rho_p(r)=\frac{p M_{\rm L1}}{4\pi R_{\rm L1}^{p}}
r^{p-3},\quad 0<r<R_{\rm L1}.
\]

Oba modela imaju nenegativnu integrabilnu gustoću, pozitivnu monotono rastuću obuhvaćenu masu i istu **točnu** L1 masu/radijus/prag. Središnji cusp nije evaluiran kao konačna \(\rho(0)\); ukupan integrirani maseni profil je konačan i \(M(0)=0\). Odrezana sferna raspodjela i skok gustoće na radijusu služe samo kontraprimjeru i nisu tvrdnja da stvarni CompaSO ili NFW halo ima taj rub.

Prosječna gustoća unutar \(x=r/R_{\rm L1}\) je \(\Delta_{\rm L1,mean}^{\rm toy}x^{p-3}\) puta sintetička srednja gustoća. Kako je ovdje \(\Delta_{200c,\rm mean}^{\rm toy}>\Delta_{\rm L1,mean}^{\rm toy}\), jedinstveni traženi korijen je *unutar* L1 radijusa:

\[
x_{200c}(p)=
\left(\frac{\Delta_{\rm L1,mean}^{\rm toy}}
{\Delta_{200c,\rm mean}^{\rm toy}}\right)^{1/(3-p)},
\qquad
\frac{M_{200c}(p)}{M_{\rm L1}}=x_{200c}(p)^p.
\]

| Bezdimenzijska veličina — isključivo matematički test | \(p=1\) | \(p=2\) |
|---|---:|---:|
| Identična početna L1 masa | 1 | 1 |
| Identični početni L1 radijus | 1 | 1 |
| \(r_{200c}/R_{\rm L1}\) | 0.8794237574740149 | 0.7733861452097149 |
| \(M_{200c}/M_{\rm L1}\) | 0.8794237574740149 | 0.5981261296023422 |

**Rezultat:** premda je idealizirana ulazna L1 masa/radijus identična, razlika 200-kritičnih masa je **0.2812976278716727 L1 mase**. Ovo je identifikabilnost pod *nedovoljno specificiranim* radijalnim profilom, a ne statistička greška realnoga Abacus kataloga. Ne može se prenijeti 28.1 % na stvarne c000 haloe ili na opaženi LRG×ELG signal.

## Neovisni certifikat i dokumentirani real-data STOP

[Završni CI 36421065060](https://github.com/dvlahek/stress-energy-closure/actions/runs/36421065060) **SUCCESS**. [Izvorni program](../scripts/audit_eboss_dr16_a03_e17d2b3a_l1_m200c_nonidentifiability.py) zatvara točne radijalne srednje gustoće i masene jednadžbe i odbija tri sintetičke negativne kontrole: izmišljenu provenijenciju ASDF, direktorij z=.95 kao zamjenu za stvarni z iz zaglavlja i CompaSO L1 masu/radijus kao automatski M200c. [Neovisni pure-stdlib program](../scripts/audit_eboss_dr16_a03_e17d2b3a_independent_decimal_bisection.py) zasebno izvodi flat matter+Λ sintetičku aritmetiku s **Decimal 55-znamenkastom preciznošću**, 200-koračnom bisekcijom nalazi oba SO korijena te neovisno reproducira pune originalne izvještaje i odbija SHA-tampering. Najveći skalirani original–nezavisni gap je **2.220446049250313e−16**; maseni 200-kritični residual originala je do 2.220446049250313e−16 i neovisne Decimal rekonstrukcije do reda 10⁻⁵⁴.

[Trajni originalni izvještaj](../source_data/eboss_dr16_a03_e17d2b3a_archived_CI_2026_09_28/e17d2b3a_original_l1_200critical_nonidentifiability.json), SHA256 0b098291122eb51cbd3e113e7edd2276ba7c9b7c38669008486f0f73f897277e; [neovisni certifikat](../source_data/eboss_dr16_a03_e17d2b3a_archived_CI_2026_09_28/e17d2b3a_independent_decimal_bisection_200critical_replay.json), SHA256 401a3a6aa8773d165029a5dfa9b2b2123edf24eff8193b64979208c8f4308e9f; [trajni SHA manifest](../source_data/eboss_dr16_a03_e17d2b3a_archived_CI_2026_09_28/archive_manifest.json), Git blob a5557fc2ae5ec49704fbd70a8efc980d279e5555. Originalni i nezavisni bajtovi i njihov manifest trajno su pohranjeni u source_data/ **audit grane**; nisu samo prolazni Actions artefakti.

## Što mora biti poznato prije stvarnog E17D2b3b

Za stvarnu primjenu treba novi zasebni data-access/provenance protokol. Iz službeno provjerenoga proizvoda moraju biti poznati trajni dataset i konkretan c000 box/phase, raw ili cleaned inačica, providerovi checksumovi i lokalni SHA256, **stvarni z iz svakog ASDF zaglavlja**, masene jedinice, ParticleMassHMsun, N, SODensityL1 i njegova stvarna gustoćna referenca. [Abacus data products](https://abacussummit.readthedocs.io/en/latest/data-products.html) posebno upozoravaju da SO_radius može imati specijalnu vrijednost ako prag nije dosegnut. Ne smije se koristiti bez odgovarajućega sentinel/status tumačenja.

Za **realni M200c kroz isti merger tree** treba validiran providerov M200c s dokazivom radijalnom definicijom ili dovoljno kompletna pozicijska čestična okolina istoga haloa u svakoj stvarnoj epohi za konzistentno ponovno mjerenje 200ρcrit radijusa. Potreban je ispravan centar, periodic wrapping, obuhvat relevantnih halo i field čestica, prava gustoćna referenca, usklađen cleaning i stvarni progenitor/descendant ID. Primarni 10% čestični RV podsampovi i sekundarni PID-only podsampovi sami po sebi nisu **točno** puno 200ρcrit mass remeasurement. Službeni release navodi providerov **POSIX/GNU cksum** u checksums.crc32; ne smije se uspoređivati kao da je zlib.crc32.

Čak i nakon pravih halo stabala, originalna CLASS \(A_s\) i \(\tau_{\rm reio}\) razlikuju se od c000, a Abacusov smooth-neutrino N-body nije nelinearni neutrinski Vlasovljev wake za E8 F+/F−. Treba neovisna masivno-neutrinska i LRG/ELG HOD/selection kalibracija te stvarni retarded Einstein–Vlasov long–short oba kraka. A03 pair-z RR, A04, originalni SGC reverse, 48k i zasebni eBOSS 3-point window/covariance i dalje su otvoreni. Originalni 24D **dvotočkasti** odd nije bispektar.

**Nije čitan stvarni Abacus halo, particle, PID, ASDF ni merger tree; nisu pokretani CLASS, novi mockovi ili galaktički katalog. Nema fizičkog B, eBOSS ξ, S/N ni detekcije. Observed odd SEALED, main netaknut, draft PR #1 bez mergea.**
