# A-03E17D1 — izvorna vremenska povijest kao nedostajući ulaz za retardirani Vlasovljev odziv

**27. 9. 2026. · Originalni i neovisni source-only CI PASS. Puni fizički halo/tracer finite-K bispektar BLOCKED; opaženi eBOSS odd SEALED.**

## Problem

E17A/E17B su numerički provjerili linearne CLASS oba kratka kraka. E17C je dobio homogeni faktor slobodnog leta zamrznute neutrinske distribucije. E17D0 je izračunao linearizirani Newtonov Vlasovljev gravitacijski impulsni odziv po jedinici *simboličkog vanjskog potencijala* na originalnoj konačnoj q mreži. Nijedan od tih podataka ne specificira vremensku povijest vanjskoga halo potencijala, stvarnu formaciju haloa, nelinearnu samogravitaciju ili LRG/ELG spregu.

E17D1 izolira upravo tu nedostajuću vremensku informaciju. S istom originalnom E8–E17D0 teorijskom podlogom uspoređuje dvije unaprijed definirane, pozitivne, jednako normirane vremenske povijesti vanjskoga *testnog* izvora. One su matematičke provjere linearnoga retardiranog operatora, **nisu** fizički procijenjeni halo potencijali, starosti svemira, eBOSS templates ni novi science cuts.

[Prospektivni protokol E17D1](../source_data/eboss_dr16_a03_e17d1_causal_external_history_nonidentifiability_prereg_2026-09-27.json), Git blob 791502d5385e35b358fb2d2afebf5fac77bd6ad5, commit 1578795d, zabilježen je prije prvoga novog E17D1 numeričkog računanja. SHA-zaključava originalni E8 4000q CSV, E16 originalnih 12×48 egzaktnih trokuta, E17C/E17D0 joint i E17D0 neovisni original te E17D0 runner Git blob. Izvorni F0/F+/F−, originalni E14/E15 72 source/24 contrast, LRG→ELG i oba stvarna kratka kraka nisu promijenjeni.

## Što je stvarno izračunano

Na *bezdimenzijskom testnom vremenu* u, uz originalnu E17D0 fiksnoepohno-Newtonovu normaliziranu impulsnu funkciju H_F(s), definiramo kauzalni operator

\[
\mathcal R_F(\sigma;g,u)=
\int_0^u dw\,g(w)\,H_F\!\bigl(\sigma(u-w)\bigr),
\qquad
\sigma_i=s_0\frac{|\mathbf k_i|}{0.05\,h/\mathrm{Mpc}}.
\]

Izvorne E16 vrijednosti s_0=(0,0.25,1) unaprijed su zaključane i ne znače fizički odabrana konformalna vremena. Posljednja jednadžba je **ograničeni, zamrznuto-epohni, linearni vanjsko-silni Volterrin funkcional**, a ne izračunata FRW evolucija stvarnog fizičkog gravitacijskog potencijala.

Dva pozitivna testna izvora su

\[
g_E(u)=4\,\boldsymbol{1}_{[0,\,1/4]}(u),\qquad
g_L(u)=4\,\boldsymbol{1}_{[3/4,\,1]}(u).
\]

Njihovi su vremenski integrali jednaki jedinici. Promijenjen je samo raspored djelovanja u simboličkom testnom vremenu. Pri u=1 oba izvora već su djelovala; pri u=1/2 kasni izvor još nije započeo i zbog gornje granice u u retardiranom integralu mora dati točno nulu. Ovo je kauzalni negativni kontrolni slučaj, ne opažanje.

Za izvor amplitude 4 na intervalu [a,b] i u≥a, uz b' = min(b,u), definiramo t_lo=u-b' i t_hi=u-a. Analitička vremenska integracija daje

\[
\mathcal R_F(\sigma;g,u)=
-\frac{4}{I_F\sigma}
\int_{q_{\min}}^{q_{\max}}\!dq\,q F'_F(q)
\left[
\operatorname{sinc}(\sigma q t_{\mathrm{hi}})
-
\operatorname{sinc}(\sigma q t_{\mathrm{lo}})
\right],
\qquad
I_F=\int_{q_{\min}}^{q_{\max}}\!dq\,q^2F_F(q).
\]

Za σ=0 i prije početka izvora rezultat je točno nula. Numerički izvorni [E17D1 izvršivi kod](../scripts/audit_eboss_dr16_a03_e17d1_causal_external_history.py) upotrebljava istu izvornu komadno-linearnu 4000q distribuciju i GL2 intervalnu kvadraturu kao E17D0, bez q ekstrapolacije. Pohranjuje oba povijesna odziva na u=1 i kasni uzročni kontrolni odziv na u=1/2 za svaku od 3×576×2×3 = **10.368** izvornih state/leg/phase kombinacija.

## Neovisna matematika i numerički QA

[Neovisni standard-library program](../scripts/audit_eboss_dr16_a03_e17d1_independent_finite_q_history_replay.py), bez uvoza originalnog E17D1 runnera, NumPyja, CLASS-a ili opaženih podataka, nanovo rekonstruira E16 skalarne modove iz kutnog skalarnog produkta i originalnu 4000q komadno-linearnu distribuciju. Računa *isti* vremenski integral zasebnom integracijom po q po dijelovima. Za A(q)=sinc(σq t_hi)−sinc(σq t_lo):

\[
\mathcal R_F=
\frac{4}{I_F\sigma}
\left\{
\int_{q_{\min}}^{q_{\max}}\!dq\,F_F(q)
  [A(q)+qA'(q)]
-\bigl[q F_F(q) A(q)\bigr]_{q_{\min}}^{q_{\max}}
\right\}.
\]

Originalni konačni q_max=20 rub nije zamijenjen granicom u beskonačnosti. Dodatni originalni kontrolni test numerički integrira E17D0 impulsni H_F u vremenu 48-točkovnom Gaussovom kvadraturom na dva krajnja zamrznuta E16 trokuta, oba kraka, s_0=.25 i 1, sva tri izvorna stanja i oba izvora. Testovi oddness u predznaku σ, egzaktne E16 geometrijske zamjene i uzročnog pre-source nultog odziva ostaju odvojeni od bilo kakve interpretacije galaktičkog B.

## Rezultat zamrznutoga testa

[Završni originalni i neovisni CI 36349740207](https://github.com/dvlahek/stress-energy-closure/actions/runs/36349740207) **SUCCESS**; svi unaprijed zaključani engineering/numerički kriteriji prošli su bez retuninga. Po svakom F stanju originalni izvještaj ima 576 E16 trokuta, oba kratka kraka i tri s_0 faze, odnosno 3.456 slučajeva s oba testa izvora.

| Provjera | Izvorni rezultat | Značenje |
|---|---:|---|
| FD prva originalna geometrija, prvi krak, s_0=.25: \(\lvert\mathcal R_E-\mathcal R_L\rvert\) | **0.16497643780585003** | Zaključani kontraprimjer jedinstvenosti *ograničenoga testnog modela* |
| Ista razlika za F+ | 0.164852342680776 | Originalni zamrznuti F+ bez fitanja |
| Ista razlika za F− | 0.16510053290901458 | Originalni zamrznuti F− bez fitanja |
| Kasni izvor prije početka, u=.5 | točno 0 u sva tri stanja | Retardirana kauzalna kontrola |
| Izvorna neparnost predznaka σ i odziv s_0=0 | točno 0 QA odstupanje | Matematičke simetrije istoga funkcionala |
| Maksimalni originalni time-Gauss48 vs. analitički integral gap | 1.942890293094024e-16 | Vrijednosti samo u preregistriranom spot-check skupu |
| Maksimalna nezavisna full finite-q integracija-po-dijelovima razlika | **5.244277234695005e-13** | Neovisni skalarni replay svih 10.368 uzoraka, *oba* izvora |

Oba pozitivna testna izvora imaju *isti ukupni integral* i djeluju na točno iste izvorne F/CLASS roditelje, ali daju različit retardirani odgovor. To je izravna demonstracija: **u ovom ograničenom modelu poznati zamrznuti neutrinski izvor i njegov linearni impulsni kernel nisu dovoljni za jedinstven prisilni odziv ako povijest gravitacijske sile nije specificirana**. Brojevi 0.165 i 5e-13 nisu fizikalni eBOSS signal, detekcijska značajnost, svemirska vremenska skala niti bound na stvarni bispektar. Ova konstrukcija ne dokazuje da bi *potpuno specificiran* fizički Einstein–Vlasov halo imao nejedinstveno rješenje.

Neovisni replay ima zasebnu SHA-tamper negativnu kontrolu koja je odbila izmijenjeni originalni izvještaj. [Originalni joint rezultat](../source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_original_joint_causal_history_nonidentifiability.json) SHA256 1e52880db6c60cff1fe3cc1ca0a00e553c63ca08601399334d049e8403e18781; [neovisni certifikat](../source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/e17d1_independent_full_finite_q_history_replay.json) SHA256 b082bdd916e1854a3ba886dcbb5576e704f485e874cc964e65bba089fd161d17. [Trajni arhivski manifest](../source_data/eboss_dr16_a03_e17d1_archived_CI_2026_09_27/archive_manifest.json), Git blob e5c02cd78ba85dea308f6ad5324d26aebc6627e2, pohranjuje originalna tri puna izvještaja, joint izvještaj i nezavisni izvještaj zajedno sa SHA256 bajtovima na audit grani.

## Što ostaje fizički otvoreno

Sljedeći **fizički** korak mora doći iz modela halo potencijala i njegovih početnih/vanjskih uvjeta, koji su neovisno specificirani prije usporedbe s opaženim odd podacima. Tek zatim može se računati samokonzistentan retardirani gravitacijski Vlasov odgovor, stvarna evolucija u FRW/gauge pozadini, oba kratka halo/tracer kraka i pripadni LRG/ELG high-z HOD/bias/evolution/magnification/GR nuisancei. Ograničeni E17D1 ne uklanja nijedan od tih nedostajućih fizičkih ulaza, ali izravno pokazuje zašto odabir proizvoljne fiksne halo povijesti može promijeniti odgovor čak i uz jednak ukupni izvor.

Puni finite-K halo/tracer B, stvarni eBOSS trotočkasti selekcijski prozor i neovisna 3pt kovarijanca ostaju **BLOCKED**. Originalna 24D dvotočkasta odd statistika nije bispektar. A03 pair-z RR, A04, izvorni neovisni SGC reverse i 48k random-konvergencija ostaju otvoreni. Nema novih CLASS/FITS/mockova, seedova, science cuts, kontaktiranja autora, eBOSS ξ, S/N, značajnosti ili unblindanja. Opaženi odd SEALED, main netaknut, draft PR #1 bez mergea.
