# A-03E22 — dva fizička kanala iz istoga wakea, zajednička sila i pravi eBOSS estimand

**29. 9. 2026. · Problem-first teorijska sinteza postojećih E13/E17D1/E20/E21 i izravno provjerenoga objavljenog halo modela.** Nova E22 računanja provjeravaju isključivo originalne source-only koeficijente i sintetičke granične slučajeve. Ovaj dokument nije nova FRW Einstein–Vlasov simulacija, ne identificira fizikalan halo profil, ne dobiva LRG/ELG HOD, ne otvara opaženi 24D odd i ne predstavlja eBOSS detekciju.

## 1. Zajednički fizikalni izvor: stvarna povijest halo potencijala

Zamrznuti E17D0 izračunava linearizirani gravitacijski Vlasovljev impulsni odgovor po jedinici *vanjskog testnog potencijala*. Zamrznuti E17D1 zatim dokazuje da dvije različite pozitivne simboličke povijesti vanjske sile s jednakim integralom daju različite retardirane rezultate: pri prvome originalnom FD testu apsolutna razlika je 0.16497643780585003. To je broj bezdimenzijskoga testnog funkcionala, **ne fizička sila, halo drag ili eBOSS amplituda**. Fizički potencijal haloa, njegovo okruženje i početni uvjeti zato su zajednički nedostajući ulazi za OBA dolje opisana kanala.

Formalno, za poznatu početnu fazno-prostornu raspodjelu F i uz stvarni, tek specificirani izvor \(\Phi_{h+\mathrm{env}}(\boldsymbol k,\eta')\) kauzalni linearni neutrinski odgovor može se zapisati kao

\[
\delta\rho_{\nu,F}^{\,w}(\boldsymbol k,\eta)
=\int_{\eta_i}^{\eta}d\eta'\,
\mathcal K_F(\boldsymbol k;\eta,\eta')
\Phi_{h+\mathrm{env}}(\boldsymbol k,\eta') .
\]

Kernel \(\mathcal K_F\) označava pravilno normalizirani Vlasovljev odgovor u odgovarajućoj FRW/gauge konvenciji. Gornji izraz je **definicija traženoga fizičkog closurea**, ne tvrdnja da je taj fizički proizvod već izračunat. Originalne F0/F+/F− raspodjele ostaju zamrznute. Vanjska literatura Okoli, Scrimgeour, Afshordi i Hudson, *Dynamical friction in the primordial neutrino sea*, MNRAS 468 (2017) 2164–2175, DOI https://doi.org/10.1093/mnras/stx560, izvodi FD-like halo-model neutrinsko dinamičko trenje s one-halo i two-halo doprinosima. Njihov aproksimacijski zakon masa/neutrino-raspodjela i numeričke prognoze **nisu** ista stvar kao naš F± epoch-matched Einstein–Vlasov response, pa ih nije dopušteno preslikati samo množenjem mase na četvrtu potenciju.

## 2. Linearni kanal: drag → brzina haloa → opažena brojnost

Za fizički definiran halo mase M u fizikalnim koordinatama force-per-mass iz stvarnoga wakea formalno je

\[
\boldsymbol a^{\,w}_{h,F}
=-\frac{1}{M}\int d^3R\,\rho_h(\boldsymbol R)
\nabla_{\!R}\Phi_{\nu,F}^{\,w}(\boldsymbol R).
\]

U linearnom, lokalnom i isotropnom **toy** limitu može se približiti kao \(\boldsymbol a^{\,w}_{h,F}=\gamma_F(M,t)[\boldsymbol v_\nu-\boldsymbol v_h]\), s \(\gamma_F\) nepoznatim. Za stvaran halo u općem slučaju potreban je uzročni, prostorno i vremenski nelinearan odgovor, a ne konstanta \(\gamma_F\). Ako se linearizira oko zajedničke CB brzine, rezultat za relativni halo/galaxy velocity može se samo formalno pisati

\[
\boldsymbol u_{a,F}(k,t)
=\int_{t_i}^{t}dt'\,
G_{a,F}(k;t,t')\gamma_{a,F}(k,t')\,
\boldsymbol v_{\nu c,F}(k,t')+\cdots .
\]

Halo-masa/HOD integracija i okolinski **parni** response nisu automatski faktorizirani u neovisne LRG i ELG mase. E17D1 nije specificirao stvarno \(\gamma_{a,F}\), \(G_{a,F}\), početno \(\boldsymbol u\), povijest halo-potencijala ili tracer-mapping.

Za opaženu galaktičku brojnost pojednostavljena linearna struktura je

\[
\Delta_a=b_a\delta_{cb}
-\mathcal H^{-1}\partial_\parallel
v_{g,a,\parallel}
+\mathcal D_a v_{g,a,\parallel}
+\text{drugi relativistički, gravitacijski i selekcijski članovi}.
\]

\(\mathcal D_a\) je Dopplerov koeficijent opažajne brojnosti koji ovisi o evoluciji/magnification i geometriji; nije izmjereni neutrinski coupling. Ako wake inducira \(\boldsymbol u_{a,F}=\beta_{a,F}\boldsymbol v_{\nu c,F}\), onda Kaiserov dio \(-\mathcal H^{-1}\partial_\parallel u_{a,\parallel}\) daje pri linearnom potencijalnom vjetru **realni parni \(\mu^2\)** cross uz gustoću. Ne daje automatski imaginarni odd. Nediferencirani Dopplerov/selection-like dio može ga dati s \(c_{a,F}=\mathcal D_a\beta_{a,F}+c_{a,F}^{\,selection}\), ali **ni \(\beta\), ni \(\mathcal D\), ni \(c^{selection}\) nisu izračunati za originalne eBOSS LRG/ELG**. Čista real-space lokalna brojnost u isotropnom modelu nema vanjski LOS vektor za nediferencirani \(v_\parallel\).

Za originalni E12/E13 Fourierov predznak \(P_{\delta_{cb},v_F}=+i\mu C_{v,F}\), E21 identitet ostaje

\[
\operatorname{Im}P_{LE,F}^{\,linear}
=\mu(A_L c_{E,F}-A_E c_{L,F})C_{v,F}(k).
\]

Jedan međutracerski odd constrainira samo tu jednu kombinaciju, ne oba pojedinačna \(c\). Usporedba F± u istim fizičkim jedinicama dodatno zahtijeva da se koeficijenti \(c_{a,F}\) fizički definiraju po stanju; držanje njihove numeričke vrijednosti fiksnom bio bi **dodatni modelni uvjet**, ne rezultat E13.

## 3. Kvadratni kanal: wake-matter coupling → mixed B → LOS

U drugom, matematički odvojenom ograničenom modelu galaktički wake izvor ima oblik

\[
S_{a,F}(\boldsymbol k,\eta)=\int_{\eta_i}^{\eta}d\eta'
\int\frac{d^3q}{(2\pi)^3}
\Gamma^{F}_{a,i}(\boldsymbol k,\boldsymbol q;\eta,\eta')
v_{\nu c,i}^{(1)}(\boldsymbol q,\eta')
\delta_{cb}^{(1)}(\boldsymbol k-\boldsymbol q,\eta')
+\cdots .
\]

Stvarni \(\Gamma^F_{a,i}\) mora doći iz istoga, fizikalno specificiranog halo-potencijala, njegovih jednadžbi gibanja i okolinski uvjetovanoga HOD-a, ne iz E17D1 simboličkih testnih povijesti. Kao u E20, vodeći cross s linearnom gustoćom sadrži miješani \(B_{v_i\delta\delta}\). Za **zajednički centrirana linearna Gaussova polja** \(\langle v^{(1)}\delta^{(1)}\delta^{(1)}\rangle=0\). To ne isključuje gravitacijom inducirani tree-level korelator iz prvih nelinearnih polja:

\[
B^{tree}_{v\delta\delta}
=
\langle v^{(2)}\delta^{(1)}\delta^{(1)}\rangle
+\langle v^{(1)}\delta^{(2)}\delta^{(1)}\rangle
+\langle v^{(1)}\delta^{(1)}\delta^{(2)}\rangle .
\]

Ti drugi redovi imaju svoje fizikalne \(\nu\)-CB i CB transfer/response kernele. Zamrznuti E17A obrađuje linearne kratke krakove, a E17D1 testnu vanjsku povijest, **ne izračunava** nijedan puni nelinearni halo–galaxy \(B^{tree}\) ili \(\Gamma^F_{a,i}\). U isotropnom parity-invariant real-space modelu skalara bez izdvojenoga opažačkog LOS-a čak ni nenulti mixed B ne mora dati \(\operatorname{Im}P_{S\delta}\). Treba posebno izvesti realni redshift-space/selection LOS-neparni operator. Nemojte zbrajati ovaj efekt s linearnim c-kanalom kao dva neovisna fizička signala dok isti wake nije konsistentno razložen bez dvostrukog brojanja.

## 4. Zajednički wake nije ni dva neovisna pomaka ni univerzalna nula

Okoli i suradnici u dijelu o signal-to-noise izričito ističu da susjedni haloi dijele velik dio velike neutrinske wake strukture. Stoga razlika njihovih izoliranih jednopojedinačnih halo-pomaka nije valjan opažajni parni signal. Za doslovno **istu glatku zajedničku silu** uzorkovanu na dva položaja razmaknuta za \(\boldsymbol r\), egzaktni Fourierov identitet daje

\[
\Delta\boldsymbol a_{\rm com}(\boldsymbol k)
=2i\sin(\boldsymbol k\!\cdot\!\boldsymbol r/2)
\boldsymbol a_{\rm com}(\boldsymbol k),\qquad
|\Delta\boldsymbol a_{\rm com}|
\le |\boldsymbol k\!\cdot\!\boldsymbol r|\,
|\boldsymbol a_{\rm com}|.
\]

To je geometrijska long-wavelength supresija **diferencijalne sile**, ne granica cijeloga neutrinskoga odd signala. Za četiri originalna duga E13 K i šest E4 s-centara samo formalni raspon \(Ks\) jest 0.03–0.65, a ne izmjereni common-wake leakage ili dokazan limit za svih 24 komponenti. Različiti biasi ili različiti stvarni Dopplerovi koeficijenti mogu proizvesti međutracerski odd i kada su inducirane brzine dvaju tracera zajedničke. [Izvorni E22 source-only CI](https://github.com/dvlahek/stress-energy-closure/actions/runs/36568450115) potvrđuje oba matematička granična slučaja bez halo podataka.

## 5. Posebna zamka pri prenošenju Okolijeve prognoze na naš 24D

U objavljenome FD/halo-model radu Okoli i suradnika fazni doprinos u njihovoj Eq. (33) sadrži \(\boldsymbol k\cdot\boldsymbol v_{\nu c}\). Njihov signed imaginarni cross-power (Appendix A, Eq. A8) zbog te faze mijenja predznak s koherentnim neutrinskim vjetrom. Istodobno prognoza kvadrata statističke osjetljivosti u njihovoj Eq. (34) koristi **usrednjeni kvadrat** vremenske derivacije faze, \(\langle|\dot\phi_{\boldsymbol k}|^2\rangle\). Pod eksplicitnim uvjetom simetričnih predznaka i wind-sign-independent uzorkovanja moguće je imati \(\langle\dot\phi\rangle=0\) uz \(\langle|\dot\phi|^2\rangle>0\). To je matematička razlika dvaju estimanda, ne kritika njihove conditional/per-mode statistike.

Njihova prognoza za standardni multi-tracer RSD, na njihovim FD/mass/halo/survey pretpostavkama, ne identificira očekivanu srednju vrijednost **našega originalnog neponderiranoga eBOSS 24D odd estimatora** niti može služiti kao naša detekcijska značajnost. Naš originalni E13 \(P_{\delta v}\ne0\) dopušta zasebno linearni Doppler/selection-like cross, ali sam po sebi ne pretvara *kvadratni* conditional wake phase u nenulti Gaussian-leading \(\langle v\delta\delta\rangle\). Druga istraživačka ruta, Zhu i Castorina, PRD 101, 023525 (2020), DOI https://doi.org/10.1103/PhysRevD.101.023525, namjerno koristi parity-odd bispektar/rekonstrukciju relativne brzine: to **nije** naš originalni 24D 2pt opservable.

## 6. Što smo doista numerički provjerili i što ostaje otvoreno

[E22 prije nove QA zaključani protokol](../source_data/eboss_dr16_a03_e22_dual_channel_halo_response_comparison_protocol_2026-09-29.json), Git blob 5832921806cdbf4dd5ab3fba06e3ce6903d4f025, jasno razlikuje prethodno pregledanu literaturu/E13 rezultate od novih graničnih testova. [E22 QA program](../scripts/audit_eboss_dr16_a03_e22_dual_channel_halo_response.py) na originalnih četiri E21 K čvorova ponovno računa fizički konzistentne F± brzinske omjere, numerički provjerava samo sintetički konstantni causal drag ODE, dvije vrste tracer null/reversal uvjeta i E20 Gaussian/parity gate. [CI 36568450115](https://github.com/dvlahek/stress-energy-closure/actions/runs/36568450115) SUCCESS, sedam negativnih kontrola. Numerički F+−F− kontrast \((C_{v,+}-C_{v,-})/[(C_{v,+}+C_{v,-})/2]\) na originalnim K=.001,.002,.003,.005 iznosi 0.0071215920, 0.0048086717, 0.0018832604 i 0.0018243680. To su **0.712%, 0.481%, 0.188%, 0.182%** originalnoga linearnog velocity-density transfera, ne ograda fizičkih LRG×ELG odd razlikâ ako \(c_{a,F}\) nije poznat. [Mali E22 izvorno izvedeni status](../source_data/eboss_dr16_a03_e22_original_E21_four_K_dual_channel_physical_gate_source_only_summary_2026-09-29.json) pohranjuje omjere i CI provenance, ali izvorni GitHub Actions ZIP artefakt nije ponovno byte-hashiran u ovoj analizi.

**Potrebni fizički uvjeti za linearni kanal:** stvarna F- i epošno određena halo/env povijest, galaxy velocity/selection mapping i opažajni \(\mathcal D_a\), HOD i zajednički pair-conditioned wake, full-k/z te standardni GR/wide-angle nuisance. **Za kvadratni kanal:** isti fizički izvor plus drugi red Vlasov/CB, puni \(B_{v_i\delta\delta}\), halo \(\Gamma_i\), LOS-odd projekcija i zasebna definicija opservablea ako prelazimo na 3pt. Za oba nedostaju još originalni A-03 physical pair window i opravdano nezavisni A-04 inference mocks. Ne kontaktirati autore prema korisnikovoj prethodnoj odluci; nema novih stvarnih podataka, seedova/cutova, retuninga E8/E16/F±/48k, WSL-a, opaženih galaksija/odd, fizičkoga S/N, promjene main ni mergea draft PR-a.
