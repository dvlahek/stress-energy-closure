# A-03E23 — izvorni high-z eBOSS estimand, javni LRG/ELG HOD i stroga granica redshift transporta

**29. 9. 2026. · Znanstvena provjera objavljene literature i analitički redshift rezultat, bez novog lokalnog numeričkog, mock, CLASS, FITS ili opaženog odd računanja.** Ovaj dokument razrješava moguću pogrešku: objavljeni *puni* eBOSS LRG×ELG efektivni redshift ne smije se zamijeniti s originalnim projektom namjerno definiranim high-z \`[0.9,1.0)\` prozorom. Uporaba književnih HOD masa ostaje modelno i uzorački uvjetovana. Nema novih science cuts ni retuninga originalnih E8–E22.

## 1. Koji je zapravo naš opservable?

[Izvorni E6 physics/window STOP](EBOSS_DR16_A03_E6_THEORY_TO_EMPIRICAL_WINDOW_BRIDGE_2026-09-27.md), Git blob \`b9515cb95405dc271a8be767c78ab70c35c1ff71\`, eksplicitno navodi da **postojeći E0/E1 LRG×ELG empirijski RR operator** vrijedi na izvornom \`[0.9,1.0)\` high-z intervalu i da je **parni redshift marginaliziran**. Sklop E0/E1 ima dvije odvojene kape, 120 fine-s, 6 coarse-s, originalni mid-point LOS, LRG→ELG orijentaciju, ℓodd=1,3 i 24D **pilot**. Njegov izvorni NPZ nije otvoren u ovom E23. Korisnik je prethodno odabrao autentične published tracer-specific randome bez kontaktiranja autora; ne rekonstruiramo povijesnu zajedničku masku.

Zbog toga je originalni z=.95 **unutrašnji matematički E8/E10/E13/E17/E21 čvor** u namjernom high-z prozoru, NE potvrđeni pair-weighted \`z_eff\`. Naš E13 izravni CLASS transport računan je samo na z=.945,.95,.955. Ti čvorovi obuhvaćaju samo podinterval širine .01 unutar izvorne širine .1. Nije izračunana uniformna kontrola redshift evolucije fizikalnog galaxy/wake izvora kroz cijeli high-z interval.

Vanjska literatura pripada DRUGAČIJIM uzorcima i težinama:
- [Wang i sur., MNRAS 498 (2020) 3470–3483, DOI 10.1093/mnras/staa2593](https://doi.org/10.1093/mnras/staa2593), za objavljeni full LRG+ELG multitracer configuration-space joint BAO/RSD i cross parove na razmacima 25–150 h⁻¹ Mpc navode \`z_eff=.77\`; to NIJE redshift izvornog high-z E0/E1 operatora ni naš ℓodd/s-bin pair estimator.
- [Zhao i sur., MNRAS 532 (2024) 783–804, DOI 10.1093/mnras/stae1452](https://doi.org/10.1093/mnras/stae1452), u svojoj *širokoj*, CMASS+eBOSS LRG (z=.6–1.0) × ELG analizi navode Cross(NGC) \`z_eff=.763\`, Cross(SGC) \`z_eff=.774\`; njihovi odvojeni LRG/ELG, footprint i pair-k ponderi nisu naš high-z 24D. Ti brojevi ne zamjenjuju originalni z=.95 bez nove sample-matched kalibracije.
- [Alam, Peacock, Kraljic, Ross i Comparat, MNRAS 497 (2020) 581–595, DOI 10.1093/mnras/staa1956](https://doi.org/10.1093/mnras/staa1956), za *šire* eBOSS multitracer HOD razmatraju z=.7–1.1 i napominju pad gustoće LRG na z>.9. To upozorava na znanstvenu snagu high-z LRG×ELG testa, ali samo po sebi ne određuje broj/parno ponderiranje naših originalnih mock/real visokoz parova.

**Ispravka interpretacije:** iz objavljenih širih sample \`z_eff≈.77\` ne slijedi da je originalni E8/E10 \`z=.95\` „pogrešna epoha”. Slijedi samo da external broad-sample published effective-z nije dokaz efektivnog redshifta našega high-z estimanda. Svaki budući redshift pomak teorije mora biti odabran isključivo prema vlastitom nezavisno utvrđenom pair supportu, nikada prema opaženom odd rezultatu.

## 2. Matematički dopuštena redshift aproksimacija ima dokaziv *uvjet*, ne numeričku potvrdu

Neka pozitivna, normalizirana pair-weight mjera \(\mathrm d\omega(z)\) u jednoj fiksnoj cap/fine-s/μ ćeliji ima potporu \([a,b)\), \(a=.9,\ b=1\). Neka \(T(s,\mu,z)\) predstavlja jedan fizički predwindow izvor (EV + standardni nuisance kao zasebno specificirani članovi), s derivacijama po redshiftu. Ako je \(\bar T=\int T(z)\mathrm d\omega(z)\) i \(L=\sup_{z\in[a,b]}|\partial_z T|\), za originalni matematički \(z_0=.95\) vrijedi stroga elementarna ocjena

\[
|\bar T-T(z_0)|\le L \int |z-z_0|\,\mathrm d\omega(z)\le 0.05 L.
\]

Ona vrijedi za **pozitivno ponderirani izvor u jednoj ćeliji** i pod uvjetom da je dokazani L konačan. Niti jedan raniji E8/E10/E13/E17 nije dao takav globalni L za fizikalni eBOSS galaxy source, pa izraz **ne postavlja prihvatljivu vrijednost finite-z systematic errora**.

Ako se iz odgovarajuće *iste* pozitivne mjere nezavisno dobije \(z_{\rm eff}=\int z\,\mathrm d\omega\) i \(M_2=\sup|\partial_z^2 T|\), Taylorov prvi red pri tom srednjaku iščezava i slijedi

\[
|\bar T-T(z_{\rm eff})|
\le \tfrac12 M_2 \operatorname{Var}_\omega(z)
\le \tfrac{(b-a)^2}{8}M_2
=0.00125\,M_2.
\]

Posljednja nejednakost koristi maksimalnu varijancu pozitivne mjere na intervalu širine .1. To NIJE pogreška našega 24D outputa: signed Legendre projekcije, even→odd RR miješanje, selekcija i normalizacija mogu imati signed kernels i stvarni cap/bin-specific pair-z. Njihov fizički operator-norm/z budget nije certificiran. Bez z-uvjetovanog \`R_L R_E\` ili dokazane \(\partial_z T\) / \(\partial_z^2T\) kontrole ne postoji opravdan numerički fizikalni tolerancijski prag.

## 3. Objavljeni HOD podaci daju realan red veličine, ali ne naš \(M_{200c}\)

Alam i sur. (2020) prijavljuju modelne *srednje host-halo mase* za svoj eBOSS z=.7–1.1 multitracer HOD:
- LRG: \(1.9\times10^{13}\,h^{-1}M_\odot\), modelni satelitski udio približno 17 %.
- ELG pod HMQ (high-mass quenching) HOD: \(1.1\times10^{12}\,h^{-1}M_\odot\), modelni satelitski udio približno 17 %.
- ELG pod drukčijim \(\operatorname{erf}\) HOD: \(2.9\times10^{12}\,h^{-1}M_\odot\), modelni satelitski udio približno 12 %.

Njihov halo catalogue/occupation framework temelji se na N-body FoF host halos, a naši [E17D2b0 uvjetni DM14 snapshotovi](EBOSS_DR16_A03_E17D2B0_PUBLISHED_DM14_POPULATION_HALO_SNAPSHOTS_AND_ASSEMBLY_STOP_2026-09-28.md), Git blob \`a3b1e44c2264902649f9a04a7b33474ed0ca5c87\`, na izričitoj **\(M_{200c}\)** i **\(c_{200c}(M,z)\)** konvenciji te zasebnoj DM14 Planck-2013 pozadini. Nije dopušteno ubaciti FoF srednje mase u \(M_{200c}\) Poisson/NFW halo bez konverzije i njezine nesigurnosti. Štoviše, njihov sample-broad mean mass ne određuje high-z conditional \(p(M,\mathrm{central/satellite},\mathrm{env}\mid z,\mathrm{cap},\mathrm{selection})\). Jedna objavljena srednja masa \(\int M\,p(M)\mathrm dM\) ne fiksira ni \(\int \gamma_F(M,z)\,p(M)\mathrm dM\) za nepoznati nelinearni dinamički trenje/galaktički odziv. Oboje zahtijeva masu/vrstu haloa, okruženje, vremensku povijest i response.

Alam i sur. za svoju multitracer analizu prijavljuju dodatnu 1-halo conformity strukturu pri malim razmacima ispod približno 5 h⁻¹ Mpc, gdje nezavisna HOD ocupacija LRG i ELG ne opisuje u potpunosti cross. To **nije** detekcija nove EV fizike i nije izravna kalibracija naših 20–140 h⁻¹ Mpc odd binova. Ali ne opravdava ni prešutno dodjeljivanje nezavisnih halo/selection koeficijenata dvama tracerima bez modelne pretpostavke.

E17D2b0 već daje stvarnu numeričku *uvjetnu* prostornu informaciju: na svojih 18 deklariranih DM14 masenih/z/koncentracijskih QA snapshotova i originalnih E16 valnih brojeva \(\max|1-u|=1.2605532853371404\times10^{-4}\) ([zamrznuti original joint](../source_data/eboss_dr16_a03_e17d2b0_archived_CI_2026_09_28/e17d2b0_original_joint_DM14_conditional_mass_concentration_spatial_snapshots.json), Git blob \`cc278b67858722f092f8b1e21b449afb6375c9ac\`). To ograničava samo *oblikovnu razliku tog prostornog NFW izvora u toj deklariranoj familiji i k-mreži*, NE ukupni Einstein–Vlasov wake, vremenski integrirani odgovor ili stvarne LRG hoste. Ne valja zbog te male razlike dalje fino podešavati statički profil dok su masa/evolucija/selection fizički neodređeni.

## 4. Praktični znanstveni prijelaz

Ostaje originalni zamrznuti *high-z* LRG×ELG \([.9,1)\) 24D pilot, sa zasebnim NGC/SGC i originalnim LRG→ELG midpoint LOS. E23 **nije** novi science cut, novi z-bin, novi HOD fit ni odobrenje da se koristi published broad-sample \(z_{\rm eff}\). Fizički prewindow \(\xi_\ell(s,z)\) mora imati jedan **state-specific** konzistentan originalni Einstein–Vlasov + halo/tracer model kroz relevantnu potporu, ne samo originalne tri uske z točke i dvije DM14 QA mase. Za projekt prve dvije alternativne rute E22 i dalje trebaju: (1) fizički derived drag/halo velocity bias + Doppler/selection c_L,c_E za **isti high-z sample** ili (2) stvarni mixed B i LOS-odd kernel uz isti izvorni halo/history/HOD.

**Pri budućem fizičkom pragu:** nezavisni published tracer-specific randomi mogu podržati z-uvjetovani pair-window po kapama i originalnim chunk/depth podskupovima *ako* se takvo proširenje prospektivno registrira i dopušta stvarni pristup source randomima; E23 ne čita niti preuzima ikakav stvarni redak. Do tada A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED, originalni observed odd 24D SEALED. Ne prenositi broad-sample HOD/redshift na high-z bez dokaza, ne uvoditi science cut radi povoljnije E18 random-drift ili 24D significance. No new CLASS/FITS/Abacus/WSL/mock/data, originalni E8/E10/E13/E16/E17D2/E18/E19/E20/E21/E22 nepromijenjeni, main untouched, PR #1 draft.
