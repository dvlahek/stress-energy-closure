# A-03E10 — fizička granica: uvjetni wake dipol nije srednji eBOSS dipol

**27. 9. 2026. — konkretan rezultat:** originalni E8 fizički distribuirani `F_0,F_+,F_-` i E9 odvojeno evoluirani CLASS rezultati daju pre-window **uvjetni imaginarni Fourierov LRG×ELG odd multipol pri zadanome smjeru LOS vjetra**, u jedinicama `(b_LRG-b_ELG)P_cb`. Za model dvaju jednako zastupljenih suprotnih vjetrova, **srednji intrinzični odd jest točno nula**, iako su uvjetni dipol, oktupol i njihov RMS nenulti. To je nužna provjera pariteta koja sprječava lažan prelazak E9 pozitivnoga `+1σ` vjetra u neuvjetovanu eBOSS 24D korelacijsku funkciju. *Nije dokaz da je puni stvarni eBOSS odd nula:* stvarni tracer selection, nelinearno coupling, svjetlosna geometrija i even→odd prozor nisu modelirani.

## Fizikalno značenje: što dopuštaju izvorne jednadžbe

[Okoli et al., MNRAS 468 (2017), DOI 10.1093/mnras/stx560, Appendix A Eq. (A7)–(A12)](https://academic.oup.com/mnras/article/468/2/2164/3063204) izveli su za lokalno gotovo uniforman neutrino–CDM vjetar dodatni imaginarni član dvotracerskoga Fourierova cross-powera. Za našu orijentaciju `C_LE=<δ_LRG^s δ_ELG^{s*}>` i `Δb=b_LRG-b_ELG`, njihov Eq. (A8) daje, u originalnoj linijskoj RSD aproksimaciji,

\[
\operatorname{Im}P_{LE}(k,\mu\mid\mathbf v)
 = \Delta b\,P_{cb}(k,z)\,\mu^2\,
 \frac{\dot\phi_{\mathbf k}(z,\mathbf v)}{H(z)} .
\]

`μ=hat{k}·hat{n}` jest LOS Fourier kut; `φ_{\mathbf k}` dodatno ovisi o **`mathbf v·hat{k}`** i o rezonantnoj okupaciji `f_s(q_{\rm res})`. Za pozitivan uvjetni koherentan model iz E9 neka je `mathbf v=v_LOS hat n`. Tada kutni nastavak E9 **ne smije** posvuda postaviti `μ=1`: u točnoj, ranije zaključanoj E8 occupancy konvenciji,

\[
q_{\rm res,s}(z,\mu)=
\frac{m_\nu|v_{{\rm LOS},s}(z)\mu|}
     {c\,T_{\nu0}(1+z)},\qquad
\phi_s(k,z,\mu)=
\phi_s(k,z,1)\,\mu\,
\frac{f_s(q_{\rm res,s}(z,\mu))}
     {f_s(q_{\rm res,s}(z,1))}.
\]

Za promjenu predznaka vjetra `mathbf v→-mathbf v`, rezonantni `q_{\rm res}` ostaje isti, a cijela `φ`, njezina derivacija i imaginarni `P_{LE}` mijenjaju predznak. To je fizikalni izvor paritetne nule u simetričnom modelu.

Normalizirani **Fourierov** multipol izvora definiramo samo kao

\[
T_\ell^s(k,z\mid\mathbf v)
=\frac{2\ell+1}{2}\int_{-1}^{1}d\mu\,P_\ell(\mu)\,
\mu^2\,\frac{d\phi_s(k,z,\mu\mid\mathbf v)}{d\ln a},
\qquad \ell=1,3.
\]

**Ovo nije `ξ_ℓ(s)`.** Tek uz stvarni `P_{cb}(k,z)`, `Δb`, puni `k` oblik, Hankelovu transformaciju s odgovarajućim fazama i empirijski validirani `R_LRG R_ELG` operator može se fizički izračunati korelacijski predložak. Pri reverziji tracera kompleksno-konjugirani cross-power ima suprotan imaginarni predznak. Za jednak bias `Δb=0`, Eq.(A8) odd član iz ove specifične fizike nula je bez obzira na jakost wake izvora. Sve te činjenice provjerene su kao negativne kontrole.

## Izvršni rezultati na **izvornim**, ne regeneriranim E8/E9 bajtovima

[Prije ovog E10 testa zaključan source-only protokol](../source_data/eboss_dr16_a03_e10_wind_parity_angular_odd_bridge_protocol_2026-09-27.json) definira redshift `z=.95` (samo matematički midpoint, ne eBOSS measured `z_eff`), komovirajući `k=.05h/Mpc`, izvornu E9 `+1σ` koherentnu LOS brzinu **za svako F stanje zasebno**, redshift stencilu `±.005`, velocity-source z step `0.004(1+z)`, `256/512` Gauss–Legendre čvorova, exact E8 `4000q` interpolaciju, sve `ℓ=1,3` i nulte/symmetry kontrole.

Izvorni E8 CSV SHA256 `bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0`, izvorni [E9 3-state CI report](https://github.com/dvlahek/stress-energy-closure/actions/runs/36325438404) SHA256 `450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890`; **ni E8 ni E9 izvorni rezultati nisu promijenjeni**. E10 [reproducibilni kod](../scripts/audit_eboss_dr16_a03_e10_angular_wind_parity_bridge.py), [audit sažetak](../source_data/eboss_dr16_a03_e10_frozen_conditional_angular_wind_parity_summary_2026-09-27.json), [GitHub CI `36328824782` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36328824782), puni JSON SHA256 `cc55a914a8be0174eb7c1e49e0e1968537cedfcf4a429d52ad66ab86cb44a488`. `512` čvorova:

| Conditioned `+1σ` LOS wind, `T_ℓ` per `Δb P_cb` | `ℓ=1` | `ℓ=3` |
|---|---:|---:|
| FD | 0.00026888572322902307 | 0.00017614445289763344 |
| `F_+` | 0.00022278876746842203 | 0.00014150296987877026 |
| `F_-` | 0.00031476164003239810 | 0.00021061456076441418 |
| `F_+−F_-` | −0.00009197287256397608 | −0.00006911159088564392 |

U istome uvjetnom modelu negativni LOS vjetar daje `T_ell(-v)=-T_ell(+v)` za svaki originalni `F`, **ne drugi neutrinski distribucijski model**. S `Pr(+v)=Pr(-v)=1/2`, i bez odabira haloa/tracera ovisnog o tom predznaku,

\[
\mathbb E_{\pm v}[T_\ell]=0,\qquad
\mathbb E_{\pm v}[T_\ell^2]=T_\ell(+v)^2,\quad\ell=1,3.
\]

Dakle, srednja vrijednost i varijanca nisu ista statistika. U objavljenom Appendixu A11–A12 upravo je potrebna statistika promjenjivoga koherentnog vjetra `<|dot(phi)|²>`; ne smije se zamijeniti kvadratom srednjega od unaprijed pozitivnoga smjera bez znanja selekcijske korelacije.

## Koje smo fizičke pretpostavke isključili

Ovo je **strogo simetrični model** znakova pri istoj magnitudi i pri identičnom, predznakom neuvjetovanom galaxy tracer selection. U stvarnom eBOSS-u LOS prozor i selekcija, korelacije wind–halo–galaxy, nelinearno velocity–density mode coupling, posebni relativistički odd članovi i finite-s/mu even→odd leakage mogu dati drukčiju distribuciju ili odvojeno nenulti observable. Njihov doprinos nije ovim E10 auditom izračunan ili ograničen; zato ne tvrdimo opći teorem `ξ_odd=0` za galaksije.

E10 kod **ne otvara** originalne opažene galaksije, opažene randome, observed odd niti puni eBOSS E0/E1 RR NPZ. Zeleni CI ne certificira NGC/SGC fizički `R_LRG R_ELG` operator. Izvorni SGC neovisni reverse problem i finite-random 48k konvergencija ostaju odvojeni.

## Jedna fizički valjana izlazna ruta i što nedostaje

U istoj simetričnoj dvogranoj raspodjeli valjan **uvjetni**/velocity-weighted statistički objekt jest

\[
S_{\ell}^{(v)}=
\mathbb E_{\pm v}\!\left[
\frac{v_\mathrm{LOS}}{|v_\mathrm{LOS}|}T_\ell
\right]=T_\ell(+v).
\]

Zato je moguće imati `S_ell^(v)≠0` i istodobno `E[T_ell]=0`. Međutim, eBOSS 24D originalna dvotočkasta `ξ_ℓ(s)` **ne sadrži** tu velocity-wind težinu. Praktična realizacija `S_ell^(v)` zahtijeva *neovisan, validiran i nezamagljen* rekonstruktor relativnoga neutrino–CDM vektorskog polja i halo/tracer conditional selection model, novi unaprijed zaključani estimator i njegov selekcijski prozor te dovoljno nezavisnih mockova za novu kovarijancu. Nije dopušteno umjetno pripisati `v>0` svim galaksijama ili koristiti originalni `σ_R16` kao stvarnu brzinu haloa.

Druga ruta je iz originalnog unconditional estimatorovog cilja fizički izvesti stvarnu non-Gaussian/selection-induced **wind–density–tracer correlator** s predznakom i validirati ga u neovisnom ensembleu, ne iz originalnog E9 scalar `+1σ` ranga. Ako takva korelacija ostane nula, originalni 24D estimator nema numerički predefiniran **nenulti** EV odd mean u našem modelu i opaženi odd se ne smije unblindati radi `>3σ` testa.

**Reprodukcijski erratum bez science retuninga:** prvi E10 CI `36328749137` ispravno je stao na SHA gateu jer je *ponovno izračunan* frozen E8 CSV u drugom NumPy okruženju bio fizički gotovo jednak, ali imao drugi bajtni SHA256 `935e4ab1...`; original `bf8f48de...` je ostao neizmijenjen. E10 CI nije relaksirao toleranciju. Ispravljen workflow koristi **originalni E8 Actions artefakt** `36307484339` i originalni E9 artefakt `36325438404`, validira oba stroga SHA256 i tek tada izvršava E10. `256→512` najveći relativni kvadraturni gap `1.615576e−6` nasuprot unaprijed registriranom QA `1e−4`, nula upozorenja.

**Odluka:** E10 source-only paritetni audit **PASS**. Korak prema stvarnom eBOSS EV predlošku **BLOCKED BY ESTIMAND/SELECTION**, čak i prije fizičkoga RR windowa i A04 covariance. Opaženi odd **SEALED**, ne preuzimati nove kataloge ili mockove bez odluke, ne promijeniti originalne cuts/seedove i ne mergeati draft PR u `main`.
