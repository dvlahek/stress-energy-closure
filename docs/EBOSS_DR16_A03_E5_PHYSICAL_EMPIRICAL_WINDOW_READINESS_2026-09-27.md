# A-03E5 — fizički ograničen empirijski LRG×ELG prozor: kriteriji spremnosti

**Stanje 27. 9. 2026.:** ovo je plan i *source-only STOP gate* nastao **nakon** E4. Nije registrirani novi FITS eksperiment, dokaz konačne fizičke maske, nova numerička tolerancija ni dopuštenje za čitanje opaženog odd vektora. [Izvorna E5/A-04 odluka](../source_data/eboss_dr16_a03_e5_a04_physics_first_readiness_protocol_2026-09-27.json) zaključava izvorne roditelje i otvorene uvjete. Korisnička odluka ostaje: bez kontakta s izvornim autorima, novih preuzimanja i `main` promjena.

## Problem koji E4 nije riješio

Arhivirani [E4 18/18 rezultat](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_report_2026-09-27.json) i [neovisni manifest](../source_data/eboss_dr16_a03_e4_mock_galaxy_nested_random_density_1200_2400_4800_uploaded_manifest_2026-09-27.json) imaju pozitivnu RR potporu u svih 144 ćelija za sve 54 realizacija×kapa×razina. To nije novi uspjeh pokrivenosti u odnosu na originalnih 1200: i originalni A-02/E3 imaju 144/144. E3 sparse RR rupe odnose se na smanjene poduzorke 600 i 900.

Prethodno zadani E4 opisni medijan relativne L1 promjene cijelog 6×24 estimatora `ξ` na 2400→4800 iznosi 0,438741305822 u NGC i 0,434715736514 u SGC. To **ne certificira konvergenciju estimatora ni kod 4800 randoma**. Nije dopušteno retroaktivno proglasiti razliku prihvatljivom novim pragom ili ponavljati nove seedove do povoljnog ishoda. Devet ID-jeva ostaje ista kodno-transportna skupina; nije uzorak za inferencijsku 18D kovarijancu.

## Što točno možemo zvati empirijskim prozorom

U odabranom putu radimo samo s neizmijenjenim, autentificiranim objavljenim LRG i ELG tracer-specifičnim randomima. Empirijski parni prozor čini parno vagana raspodjela `R_LRG R_ELG`, uz neovisne normalizacije cross-Landy–Szalay članova. Time se **ne rekonstruira** povijesna kontinuirana LRG MANGLE × ELG BRICKMASK maska. Kasniji javni skript bitova 8–11 nije dokaz izvršne 2020. produkcijske verzije; dvije objavljene loše ploče ne daju njihovu izvornu efektivnu geometriju. Ne izmišljati zajedničku tvrdu masku, radijus ploče niti ponovno primjenjivati veto na već objavljene clustering datoteke.

Postojeći E0/E1 random-only operator te E2/E3/E4 mock-galaxy testovi su prethodne **numeričke kontrole**, ne mjerenje fizičkog even-to-odd leakagea. E2 ±5 % ELG-random weight ramp je arbitrarna stres proba; E3 post-DD ±0,02 injekcija je analitička kontrola, ne populacija fizičkih wake galaksija.

## Budući acceptance mora biti fizički i unaprijed definiran

Prije novog input-only ili mock-only testa treba zaključati konačni odd observable i redoslijed komponenti, teorijski fizički wake predložak i njegovu propagaciju kroz `R_LRG R_ELG`; zasebne kape NGC/SGC, ELG `eboss21/22/23/25` chunkove i stvarno raspoloživo imaging-depth/radial `n(z)` uvjetovanje; originalne verzije, težine, geometriju i LOS. Potrebno je unaprijed izabrati neovisne kontrole finite-random šuma i selekcijskog curenja u **jedinicama konačnog observablea**, zajedno s opravdanim budžetom fizikalne nesigurnosti. Numerički prag mora slijediti iz teorijske osjetljivosti i neovisne inferencijske kovarijance, ne iz pregledanih E4 promjena. Danas **nema** takvog fizički obrazloženog broja; ostaje `NOT_DEFINED`.

U odvojenom unaprijed zaključanom protokolu definirati poznati fizički wake/null injection i recovery kroz isti released-tracer pair-window operator, kao i even-to-odd propuštanje zbog dopuštenih radijalnih/depth/chunk perturbacija. Koristiti neovisnu validaciju novih pravila koja su zamišljena nakon E4, uz sve slučajeve, ne izdvajati one s manjim odstupanjem. Ako nisu identificirane stvarne depth/selection veze u released randomima, to je ograničenje claim-a, ne razlog da se veze izmisle.

**Fail-closed uvjet:** do zaključanog konačnog observablea, fizičkog modela, numeričkog budžeta, nezavisne validacije i A-04 kovarijance, A-03 ostaje **PHYSICAL_UNCERTIFIED**. Zeleni E5 source-only CI smije potvrditi samo da se taj STOP provodi i roditelji nisu promijenjeni.

## Odvojen A-04 i opservacijski seal

A-04 definicija i kovarijancijski preduvjeti nalaze se u [zasebnom dokumentu](EBOSS_DR16_A04_EBOSS_OBSERVABLE_COVARIANCE_READINESS_2026-09-27.md). Postojeći eBOSS pilot je 6 `s` binova × 2 odd multipola × 2 kape = **24** komponente ukupno, odnosno 12 po kapi. Ne preuzimati gotov **DESI DR1** 18D dipolni vektor ili 120-mock kovarijancu kao eBOSS rezultat. Devet eBOSS mock realizacija ima rang centrirane kovarijance najviše osam.

**Ništa nije unblindano.** CI ili ovaj plan ne čitaju opažene galaksije, opažene randome ni odd vektor, ne preuzimaju nove mockove i ne odobravaju buduće preuzimanje. Originalni E4, raniji izvještaji i draft PR #1 ostaju netaknuti.
