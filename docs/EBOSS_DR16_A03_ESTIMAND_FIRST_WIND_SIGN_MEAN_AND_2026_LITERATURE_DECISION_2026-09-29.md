# A-03 — prvo identificirati opservable: uvjetni wake ≠ neuvjetovani eBOSS odd

**29. 9. 2026. · Fizikalna interpretacija već zamrznutih rezultata, ne nova numerička preregistracija, novi FITS test, novi signal ni odobrenje za unblinding.** Ovaj dokument čita postojeći E10, E6 i E18 te provjerenu vanjsku literaturu. Ne mijenja E8 F0/F+/F− 4000q, E16 576 trokuta, E4 originalni mock-only rezultat, E18 opisni rezultat, E0/E1 empirijski random operator ili E17D0/D1 ograničeni vanjski odziv. Nema WSL-a ni novog velikog podatkovnog posla.

## 1. Problem: E18 kvantificira šum **drugoga** statističkog objekta

Izvorni [E18 devet-mock rezultat](../source_data/eboss_dr16_a03_e18_posthoc_nine_mock_24d_paired_random_vs_mock_scatter_descriptive_result_2026-09-29.json) pokazuje da je joint-24D, na originalnih devet fiksnih E4 realizacija s istim 600D galaksijama, 2400→4800R paired RMS promjene \`0.0947419874387\` uz opisni među-ID high-R RMS \`0.180712396171\`; njihov omjer \`0.524269443858\` NIJE sigma, fizička varijanca ni limit. Za 1200→4800R opisni omjer je \`1.11383414719\`. Originalni E18 numerički rezultat neovisno je auditiran u [CI 36541582286](https://github.com/dvlahek/stress-energy-closure/actions/runs/36541582286), pet negativnih kontrola PASS, bez novih FITS. To je stupanj numeričke pouzdanosti staroga **neponderiranoga dvotočkastoga eBOSS** pilota, ne dokaz da je njegov očekivani EV wake mean različit od nule.

Originalni [E10 paritetni fizikalni audit](EBOSS_DR16_A03_E10_PHYSICS_OF_CONDITIONAL_WIND_AND_UNCONDITIONAL_ODD_NO_GO_2026-09-27.md), [archived E10 frozen summary](../source_data/eboss_dr16_a03_e10_frozen_conditional_angular_wind_parity_summary_2026-09-27.json), Git blob \`a9e0f49a58b9fd5e0ad975d4d14ba6521a648c28\`, ima pri originalnom ilustrativnom \`k=0.05 h/Mpc\` i matematičkom \`z=.95\` uvjetne Fourierove odd koeficijente po \`(b_LRG-b_ELG)P_cb\`. Pri pozitivnom E9 +1σ LOS vjetru razlika \`F+−F−\`: \`T_1=-9.197287256397608×10^-5\` i \`T_3=-6.911159088564392×10^-5\`. Nisu fizički predwindow \`xi_l(s,z)\`, nisu binirani eBOSS 24D broj i ne smiju se numerički uspoređivati s E18 \`0.0947\` bez Hankela, apsolutnoga bias/P_cb/galaxy odziva, redshift modela i istoga pair-windowa.

## 2. Precizni simetrični dvograni fizički rezultat

Za **jednu** zamrznutu distribuciju F i fiksne \`k,z,|v|\` pišemo \`T_F(+v)=t_F\`, \`T_F(-v)=-t_F\`, kako u E10 slijedi iz promjene predznaka gravitacijskoga wake faznoga pomaka, uz iste amplitude i suprotni LOS smjer. Ako su dva smjera jednako zastupljena i galaxy/tracer selection ne ovisi o predznaku, slijedi

\[
\langle T_F \rangle = \tfrac12 t_F+\tfrac12(-t_F)=0,\qquad
\langle T_F^2\rangle=t_F^2.
\]

Za jednostavan **isključivo teorijski** model nenegativnih relativnih selekcijskih težina \`w_+,w_-\` za dva smjera,

\[
\langle T_F\rangle_{\rm selected}
=\frac{w_+t_F-w_-t_F}{w_++w_-}
=\eta\,t_F,\qquad
\eta=\frac{w_+-w_-}{w_++w_-},\quad |\eta|\le1.
\]

To je identitet uvjetnog modela, **ne izmjerena eBOSS selekcijska asimetrija, ograda punoga nelinearnog galaktičkog odziva ili nova fitted nuisance-parametrizacija**. Trenutačni E0/E1 \`R_LRG R_ELG\` izvori nisu razdvojeni po predznaku rekonstruiranoga neutrino–CDM vjetra i ne daju \`w_+/w_-\`. Bez fizički specificiranog halou/traceru vezanoga vjetra i selekcijske korelacije \`eta\` nije identificiran.

Neka \(\mathcal H\) označava linearni Fourier→konfiguracijski prijelaz s **ispravno fiksiranom** orijentacijom, a \(M\) zadani, za oba predznaka ISTI, linearni E6 conditional finite-bin RR operator. U ovom ograničenom modelu,

\[
\frac{M\mathcal H T_F(+v)+M\mathcal H T_F(-v)}{2}
=M\mathcal H\frac{T_F(+v)+T_F(-v)}2=0.
\]

Zato **sam wind-independent RR window ne pretvara simetrični, sign-odd EV intrinsic mean u nenulti EV mean**. Istodobno E6 \(\sum_{L=0,2,4}M_{oL}\xi_L\) može dati nenulti *standardni even→odd finite-bin doprinos*. Takav doprinos ne smije se automatski imenovati EV wake. U stvarnom surveyu window, halo i galaksije mogu korelirati s vjetrom, \(\mathcal H\) i \(M\) mogu imati dodatne z/selection ovisnosti, a standardni relativistički i nelinearni članovi mogu preživjeti. Gore navedena nula stoga NIJE univerzalni teorem o stvarnom eBOSS odd signalu.

## 3. Potrebna fizička veza koja trenutačno ne postoji

Postojeći [E6 mathematical/synthetic bridge](EBOSS_DR16_A03_E6_THEORY_TO_EMPIRICAL_WINDOW_BRIDGE_2026-09-27.md) zaključao je 5 prewindow \(\ell=0,\ldots,4\) → 2 output \(\ell=1,3\), 6 coarse-s binova, obje kape, 24D pilot, source-pinned E0/E1 \`6×120\` blokove i sintetičku even/odd algebru. Nije otvorio originalni 250104-B RR NPZ niti izveo fizički normaliziran EV \(\xi_L(s,z)\). Postojeći \`code/build_lrg_elg_wake_template.py\` zasebno maksimalno normira wake i doppler, s \`absolute_wake_prediction=False\`. Postojeća DESI physical odd baza vanjski očekuje LRG/ELG bias/magnification/evolution parametre i ideal-shell \`s²ds\` umjesto stvarnoga eBOSS para. **Ništa od toga nije fizički normaliziran eBOSS wake mean.**

Prije bilo kakve afirmativne tvrdnje o originalnom **unweighted 24D** EV testu potreban je istovremeno:
- teorijski definiran *neuvjetovani* LRG×ELG two-point EV odgovor (uključujući uvjete pod kojima long-wind sign ne poništi mean), s kvantitativnim halo/galaxy tracer odgovorom i standardnim GR/wide-angle/even izvorima
- bezdimenzijski \(\xi_L(s,z)\) s apsolutnom normalizacijom i exact LRG→ELG midpoint LOS orijentacijom
- empirijski dokaziv, redshift i zasebno NGC/SGC/chunk/depth uvjetovani \`R_L R_E\` operator, samo u opsegu objavljenih tracer randoma
- nakon toga eksplicitno fiksiran opservable i neovisna eBOSS covariance/selection/finite-random coverage kalibracija.

Bez prve stavke dodatna gustoća randoma, nove mock realizacije ili samo E6 \`shape-only\` predložak **ne zatvaraju fizičku hipotezu**. Devet mockova daju 24D centrirani covariance rank najviše osam. E18 omjer nije fizički kriterij za dopuštenje unblindinga.

## 4. Znanstvene putanje su različiti eksperimenti

**Postojeći neponderirani eBOSS 2pt odd:** zadržati ga kao estimator i izveden QA. EV interpretacija kao nenulti sign-sensitive wake mean zahtijeva samostalan, izvan opaženoga odd vektora deriviran i validiran wind–halo–galaxy/selection korelacijski model. U odsutnosti takva modela dopušten je zaseban standardni odd ili općenito ograničeni empirijski test, ali nije dopuštena etiketa „EV wake detection”. Nema nove selekcije na temelju pregledanih E18 pomaka.

**Velocity-conditioned/weighted estimand:** E10 matematički identificira \(\langle \operatorname{sign}(v_{\rm LOS})T_F(v)\rangle=t_F\) u jednostavnom sign modelu. No taj *novi* estimator zahtijeva nezavisan rekonstruirani relativni neutrino–CDM vektorski smjer, halo/galaxy mapping, provjerljiv wind-conditioned pair-window i zasebne mockove. Ne „pridodavati” sign=+1 svakom opaženom paru, ne proglašavati originalnu E9 +1σ brzinu promatranom individualnom brzinom.

**Parity-odd squeezed three-point:** to je zaseban opažaj i nije originalni eBOSS 24D dvotočkasti vektor. E14/E15 source-only reducirani \(i\mu_K S\) i E16/E17 oba kratka kraka nisu puni finite-K retardirani galaktički bispektar. E17D0/D1 vanjski potencijal i history nisu nezavisno fizikalno specificirani kao puni halo/galaxy coupling; stvarni 3pt prozor/kovarijanca ne postoje. Bez toga se ni ova ruta ne zove detekcijom.

Vanjska usporedba, jasno odvojena od naših podataka: Zhu i Castorina, Phys. Rev. D **101**, 023525 (2020), DOI [10.1103/PhysRevD.101.023525](https://doi.org/10.1103/PhysRevD.101.023525), formuliraju relativni neutrino–CDM signal pomoću paritetno neparnoga bispektra. Nascimento i LoVerde, JCAP **04 (2026) 019**, DOI [10.1088/1475-7516/2026/04/019](https://doi.org/10.1088/1475-7516/2026/04/019), istražuju detekciju wakea s rekonstruiranim relativnim velocity poljem i galaxy/matter cross-correlations. Njihove prognoze pripadaju njihovim survey/tracer/field pretpostavkama i **nisu** mjerenja našega eBOSS 24D ni potvrda E8 signala.

## 5. Operativna odluka i najkraći znanstveni put

**Ne otvarati opaženi eBOSS odd vektor i ne naručivati velik matched-mock ensemble dok god ne postoji nenulti, fizikalno identificiran EV estimand koji ta opažajna definicija stvarno mjeri.** Nastavak research rada usmjeriti na izvod wind-sign-averaged halo–LRG/ELG tracer responsea ili na zasebnu teorijsku identifikabilnost Einstein–Vlasov \(F_+,F_-\) koja ne pretendira na empirijsku wake detekciju. Ako izvod pokaže nulu pod relevantnim, neovisno specificiranim pretpostavkama, to je znanstvena granica postojećega 2pt testa, a ne razlog za artificial survey cuts ili tune.

Sačuvati korisnikovu raniju odluku o empirijskom objavljenom random-only putu bez kontaktiranja autora. Originalni F0/F+/F−/E8/E16/48k, E18 opisni rezultat i parent SHA netaknuti. Nema WSL naredbe, novog Abacus transfera, novih eBOSS izvora, novih seedova, 18D post-hoc truncationa, S/N, opaženog odd ni main mergea. Draft PR #1 ostaje otvoren.
