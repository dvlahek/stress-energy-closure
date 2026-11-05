# A-03E20 — miješani trotočkasti izvor kao nužan uvjet ograničenog kvadratnog wake-odziva

**29. 9. 2026. · Novi matematički rezultat i provjerena sintetička algebra, NE fizički halo–galaxy bispektar ili eBOSS detekcija.** E19 je pokazao da postojeći marginalni LRG×ELG randomi ne određuju korelaciju galaktičke selekcije s predznakom neutrino–CDM vjetra. Sljedeći fizikalni problem jest identificirati prostorno promjenjiv halo/tracer odziv koji uopće može preživjeti usrednjavanje vjetra u neponderiranoj dvotočkastoj korelaciji. Ne uvodimo novi galaxy cut, novi estimator ili novi opservacijski uzorak.

## 1. Precizno ograničen odzivni model

Uvedimo centrirano CB gustoćno polje δ, relativno vektorsko polje v_i i zajednički wake-like kvadratni izvor S. Za dvije stvarne galaktičke tracerske fluktuacije definiramo *modelni ansatz*

\[
\delta_L(\boldsymbol k)=b_L\delta(\boldsymbol k)+\lambda_L S(\boldsymbol k),\qquad
\delta_E(\boldsymbol k)=b_E\delta(\boldsymbol k)+\lambda_E S(\boldsymbol k)
\]

i

\[
S(\boldsymbol k)=\int\frac{d^3q}{(2\pi)^3}
\Gamma_i(\boldsymbol k,\boldsymbol q;\hat n,z)
v_i(\boldsymbol q)\delta(\boldsymbol k-\boldsymbol q).
\]

Za realni S vrijedi \(\Gamma_i(-k,-q)=\Gamma_i(k,q)^*\). Lokalni primjer S(x)=v(x)·∇δ(x) ima \(\Gamma_i=i(k-q)_i\). To je samo kvadratni gradijentni prototip. **Fizikalni \(\Gamma_i\), stvarni visokoz LRG/ELG \(b_a,\lambda_a\), halo potencijal i retardirana povijest nisu izvedeni u E17D0/D1.** U E10 je fazni odgovor nelinearan u okupaciji ovisnoj o \(|v|\), stoga se njegovi koeficijenti ne smiju bez novog izvoda poistovjetiti s konstantnim \(\lambda_a\) ovoga modela.

S izvornom LRG→ELG Fourier konvencijom \(P_{LE}=\langle\delta_L\delta_E^*\rangle\), realnim koeficijentima i \(Q=P_{S\delta}\), algebarski do prvoga reda u wake-odzivu dobivamo

\[
P_{LE}=b_Lb_EP_{\delta\delta}
+\lambda_Lb_EQ+b_L\lambda_EQ^*+O(\lambda^2),
\]

\[
\operatorname{Im}P_{LE}^{(1)}
=(\lambda_Lb_E-b_L\lambda_E)\operatorname{Im}Q,\qquad
Q(k)=\int\frac{d^3q}{(2\pi)^3}
\Gamma_i(k,q)\,B_{v_i\delta\delta}(q,k-q,-k).
\]

U zadnjoj formuli delta-funkcija očuvanja impulsa iz definicije miješanog bispektra izostavljena je radi preglednosti. Zamjena L↔E kompleksno konjugira cross-power i mijenja predznak njegovoga imaginarnog dijela. Samo u posebnom slučaju **jednakoga wake-response koeficijenta** \(\lambda_L=\lambda_E\) vodeći antisimetrični koeficijent proporcionalan je \(b_E-b_L\), kao u posebnoj izvornoj E10 shemi uz njezinu konvenciju predznaka. Opći ansatz s različitim \(\lambda_a\) ne podržava univerzalnu tvrdnju da je signal nula uvijek kada su \(b_L=b_E\).

## 2. Gaussova nula nije pretpostavka neovisnosti brzine i gustoće

Ako su v_i i δ zajednički centrirana Gaussova *linearno evoluirana* polja, svaki miješani treći moment je nula: \(B_{v_i\delta\delta}=0\). To ostaje točno čak i ako je njihova **dvotočkasta korelacija nenulta**. U gornjem kvadratnom ansatzu zato \(Q=0\) i linearni wake dio \(P_{LE}^{(1)}\) ne može dati neponderirani EV odd mean. Ovaj rezultat NIJE tvrdnja da cijeli galaktički odd power ili svi neutrinski učinci nestaju; doprinosi \(O(\lambda^2)\), ne-Gaussova gravitacijska evolucija, drugačiji linearni velocity/selection operatori i standardni relativistički Doppler mogu ostati.

Strogo sintetička provjera za **zajedničku predznačnu simetriju**, koja sama po sebi nije Gaussova simulacija: X,Y su neovisni jednako vjerojatni ±1, V=X,D1=X,D2=X+Y. Sva tri srednjaka su nula, tri različita parna cross-momenta su 1, ali \(\langle VD1D2\rangle=0\). Obrnuti *moment-only* svjedok V=X,D1=Y,D2=XY ima također nulte srednjake i sve tri parne cross-momente jednake nuli, ali \(\langle VD1D2\rangle=1\). To demonstrira zašto sam power-spectrum ne određuje miješani treći moment; drugi svjedok **nije** prostorno konzistentan Einstein–Vlasov bispektar.

## 3. Još jedan nužan uvjet: neparna LOS projekcija

Nenulti \(B_{v\delta\delta}\) sam po sebi još uvijek ne jamči mjerljivi \(l=1,3\). Ako su S=v·∇δ i δ stvarni skalari u homogenu, izotropnu i paritetno simetričnu real-space ansamblu *bez izdvojenoga LOS vektora*, njihova \(\xi_{S\delta}(\boldsymbol r)\) je parna pod \(\boldsymbol r\to-\boldsymbol r\), pa je \(Q=P_{S\delta}\) realan. Stoga \(\operatorname{Im}Q=0\), čak i kada neka miješana trotočkasta statistika nije nula. Za neponderirani **redshift-space** odd potrebna je fizički specificirana projekcija koja razlikuje \(+\hat n\) i \(-\hat n\), odnosno odgovarajuća tracer/velocity/selection asimetrija. Postojeći E6 finite-bin even→odd leakage može zasebno biti nenulti, ali nije dokaz našega EV mehanizma.

## 4. Izvorni dokaz, njegova granica i sljedeći fizički rad

[Prije novog E20 QA zaključani izvorni protokol](../source_data/eboss_dr16_a03_e20_quadratic_wake_mixed_bispectrum_structural_protocol_2026-09-29.json), originalni Git blob 300f73be14306ccf45342592e842e30bda63b8c6, prethodi [E20 source-only standard-library programu](../scripts/audit_eboss_dr16_a03_e20_quadratic_wake_mixed_bispectrum.py). Program provjerava pet SHA-pinned roditelja, egzaktne racionalne četverostanjske momente, imaginarni tracer contrast i LRG/ELG reversal, reality primjera \(i(k-q)\) te pet negativnih kontrola. [Završni E20 CI 36563347448 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36563347448). Manji raniji CI 36563041019 nije prošao jer je novi čitač zatražio nepostojeći naziv E19 seal-polja; ispravljen je samo naziv polja u čitaču. Druga završna promjena bila je izričito označiti konačni sign-symmetric diskretni svjedok kao *ne-Gaussov*, bez promjene matematike. Originalni E8/E10/E19/E17D1 i prereg nisu promijenjeni. [Trajni E20 kratki izvorni CI sažetak](../source_data/eboss_dr16_a03_e20_archived_exact_moment_and_cross_power_structural_QA_2026-09-29.json) zaseban je zapis iz reproducibilnog koda i verificiranih CI logova, ne tvrdnja da je ponovno otvoren i byte-hashiran zip artefakt.

Relevantna vanjska literatura daje kontekst, a ne numerički rezultat našega projekta: [Zhu i Castorina, PRD 101, 023525 (2020)](https://doi.org/10.1103/PhysRevD.101.023525) obrađuju paritetno neparni relativno-brzinski bispektar; [Nascimento i LoVerde, *Neutrino winds on the sky* (2023)](https://arxiv.org/abs/2307.00049) izvode formalizam u kojem lokalna anisotropija neutrino–CDM korelacije ulazi u trotočkaste korelacije. Nijedan od tih rezultata nije naš dobiveni \(B_{v\delta\delta}\) ili fizički eBOSS 24D prediction.

**Operativni zaključak:** za ovaj kvadratni model potreban je stvarno izveden, epošno i kozmološki ujednačen \(B_{v_i\delta\delta}\) s \(\Gamma_i\) koji uključuje relevantan halo/galaxy response i LOS, s validiranom LRG/ELG razlikom, puni k/z razvoj, Hankel i A-03 prozor. E17A linearni CLASS kratki krakovi i E17D1 dva jednako integrirana testna halo-potencijala to ne određuju. Ne zamijeniti nepoznati nelinearni kernel proizvoljnim EdS \(F_2\), ne postavljati opaženi vjetar na +1σ E9, ne koristiti E18 opisni random/scatter omjer kao physical leakage i ne naručivati dodatne mockove prije fizikalno identificiranoga estimanda. Opaženi eBOSS 24D odd SEALED, A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED, bez WSL, novih galaksija/mockova/Abacus podataka, cut/seed retuninga ili main mergea.
