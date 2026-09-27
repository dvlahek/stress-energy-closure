# A-04 — eBOSS DR16 observable i inferencijska kovarijanca: definicija prije mock akvizicije

**Status 27. 9. 2026.:** ugovor spremnosti, NE konačno registrirani inferencijski observable. Nema novih eBOSS mock preuzimanja, parnih izračuna, opservacijskog odd pristupa ili prijenosa DESI rezultata. Povezano s [A-03E5 fizičkim acceptance uvjetima](EBOSS_DR16_A03_E5_PHYSICAL_EMPIRICAL_WINDOW_READINESS_2026-09-27.md) i [izvornim E5/A-04 source-only STOP protokolom](../source_data/eboss_dr16_a03_e5_a04_physics_first_readiness_protocol_2026-09-27.json).

## Dimenzionalnost je otvorena znanstvena odluka

Postojeći eBOSS A-02/E2/E3/E4 code-transport pilot koristi **jedan z-interval [0,9,1,0)**, šest `s` binova s rubovima 20,40,…,140 Mpc/h, `ℓ=1,3` i dvije kape NGC/SGC. To daje **6×2=12 komponenti po kapi** ili 24 komponenti za zajednički, eksplicitno definiran dvokapni vektor. Njegova 6×24 `s-μ` mreža je dodatni interni estimator, ne već registrirani konačni 18D znanstveni observable.

U zasebnoj DESI DR1 liniji [120-mock kovarijancijski artefakt](../source_data/lrg_elg_covariance_r4_120_checkpoint_2026-09-23.json) opisuje 18D *dipol* kao **tri redshift intervala × šest separacijskih binova**, uz DESI izvore i pravila. [DESI windowed finite-mock rezultat](../source_data/lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json) pripada istoj DESI liniji. Te matrice, izračunati Hartlapov faktor, likelihood ili statistički zaključci nisu primjenjivi na eBOSS DR16 i njegov jedan high-z pilot. Nije legitimno izbaciti šest eBOSS komponenti poslije E4 da se dobije isti broj 18.

Prije eBOSS A-04 računanja treba zasebno, neovisno o opaženoj odd vrijednosti i dobivenim E4 pomacima, odabrati fizikalno motiviranu dimenziju i točan redoslijed komponenti, redshift i `s` rubove, `ℓ`, cap kombinaciju, LOS, orijentaciju LRG→ELG, weights/fiducial, jedinice, nuisance matricu i linearnu projekciju Einstein–Vlasov wake/null predloška kroz validirani empirical pair window. Ako cilj ostane **18D**, definicija mora biti nova zasebno registrirana eBOSS definicija koja ne mijenja retroaktivno stari 6×24 transportni test. Do te odluke stoji `EBOSS_18D_UNDEFINED`.

## Zašto još nema kalibrirane eBOSS kovarijance

Devet izvornih matched realistic EZmock ID-jeva predstavlja devet realizacija, ne 18 neovisnih realizacija zbog dvije kape. Centrirana sample kovarijanca njihovih 9 punih vektora ima rang najviše 8. Zbog toga ne možemo invertirati ni 12D per-cap, ni 18D, ni 24D joint sample kovarijancu bez dodatne, zasebno opravdane strukture. Jednostavno dodavanje random redaka na istu mock-galaxy realizaciju ne stvara nove nezavisne svemirske realizacije i ne popravlja ovaj rang.

Za budući eBOSS ansambl treba **prije preuzimanja** zaključati nezavisne ID-jeve, isti DR16 release, odgovarajuće same-ID LRG/ELG galaksije i randome za obje kape, source-vintage i sva pravila izvlačenja te identičan A-03 prozor. Odvojeni metadata-only inventar i realan disk/CPU/network budžet prethode bilo kakvom bulk prijenosu. Nove datoteke i novi sample trenutno **nisu odobreni**.

Budući inferencijski protokol mora registrirati linearni/nelinearni nuisance prostor, rang i eigen-spektralnu stabilnost, osjetljivost na originalne ID-jeve, procjenu stvarnog `N_eff`, zamrznuto pravilo regularizacije/shrinkage i coverage/tail kalibraciju s neovisnim null i physical-injection realizacijama. Hartlapova korekcija zahtijeva opravdane Gaussian/Wishart i nezavisne mock pretpostavke i dovoljno realizacija u odnosu na **stvarni eBOSS p**. Izvorni DESI faktor nije ulazna konstanta za eBOSS. Uz `N` običnih nezavisnih null mockova empirijski rep nema rezoluciju bolju od približno `1/(N+1)`; ni samo `N=120` ne omogućuje neposrednu empirijsku kalibraciju vrlo rijetkih repova poput nominalnih 5σ.

## Kontrolna granica

A-04 ostaje `BLOCKED` dok ne postoje: potpuno zamrznuti eBOSS observable, fizički propagiran A-03 pair window i budžet sistematika, odvojen dovoljno velik i opravdano nezavisan eBOSS ansambl, numerički stabilna kovarijanca i unaprijed fiksiran inferencijski/null/coverage postupak. **Opaženi odd vektor ostaje SEALED.** Ovaj dokument ne tvrdi da se podaci čitaju, da su odobreni novi downloadi ili da je izgledan signal.
