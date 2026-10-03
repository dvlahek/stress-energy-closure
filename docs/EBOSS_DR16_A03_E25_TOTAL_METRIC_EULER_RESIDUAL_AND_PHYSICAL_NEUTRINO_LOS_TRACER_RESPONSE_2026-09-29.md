# A-03E25 — potpuni gravitacijski potencijal, Eulerov ostatak i neutrinski LOS odziv galaksija

**29. 9. 2026. · Dva izvoda i novi izvorno zaključani matematički QA, a NE kalibrirani neutrinski halo–galaxy model.** E24 je pokazao da standardni parni auto/cross clustering, unutar izričito ograničenog dvopoljnog modela, ne određuje zajednički predznak fizički još nepoznatih line-of-sight odziva LRG i ELG na relativnu brzinu neutrino–CDM. Sada postavljamo uži fizikalni problem: koji član u opaženoj brojnosti smije predstavljati **stvarno** neutrinsko halo ubrzanje, a koji je već uključen u isti gravitacijski potencijal?

## 1. Izvorna opažajna jednadžba prije uporabe Eulerove relacije

[Bonvin, Lepori i sur. (2023), *A case study for measuring the relativistic dipole of a galaxy cross-correlation with the Dark Energy Spectroscopic Instrument*, MNRAS 525, 4611–4627, odjeljak 2.1, jednadžba (2)](https://academic.oup.com/mnras/article/525/3/4611/7251503) zapisuju puni gauge-invariant linearni kontrast brojnosti uključujući standardni density/RSD, magnification/lensing, Doppler, ubrzanje i gravitacijski redshift. Za **lokalni** dio u njihovoj konvenciji, uz bezdimenzijsku peculiar velocity \(V=v_{\rm phys}/c_{\rm light}\), pozitivni \(\mathcal H=aH\), konformalni vremenski izvod \(^{\prime}\) i **ukupni** metric potential \(\Psi_{\rm tot}\), smijemo izdvojiti

\[
\Delta_a^{\rm loc}\supset
 B_a(z)\,V_{a,\parallel}
 +\frac{V_{a,\parallel}^{\prime}+\partial_\parallel\Psi_{\rm tot}}{\mathcal H},
\]

\[
B_a=1-\frac{\mathcal H'}{\mathcal H^2}
       +\frac{5s_a-2}{r\mathcal H}
       -5s_a+f_a^{\rm evo}.
\]

Ovdje \(s_a\) znači uobičajeni magnification bias, a \(f_a^{\rm evo}\) evolution bias *u konvenciji citiranoga izvora*. Ovo nisu naša izmjerena LRG/ELG high-z svojstva. Standardni Kaiserov član \(-\mathcal H^{-1}\partial_\parallel V_{a,\parallel}\) i ostali lokalni, integrirani i opažački potencijalni članovi postoje odvojeno. Ne preuzimamo samo Dopplerov član, zanemarujući pritom ubrzanje ili gradijent potencijala. [Bonvin i sur., *Dipolar modulation in the size of galaxies: the effect of Doppler magnification* (2017), jednadžba (48)](https://academic.oup.com/mnras/article/472/4/3936/4082100) neovisno ilustrira istu strukturu u odgovarajućoj dipolnoj kombinaciji.

Definirajmo **Eulerov ostatak u odnosu na isti ukupni potencijal**

\[
\mathcal E_a\equiv
V_{a,\parallel}^{\prime}+\mathcal H V_{a,\parallel}
+\partial_\parallel\Psi_{\rm tot}.
\]

Egzaktno unutar izdvojenih linearnih lokalnih članova slijedi

\[
\boxed{\Delta_a^{\rm loc}\supset
       (B_a-1)V_{a,\parallel}
       +\frac{\mathcal E_a}{\mathcal H}.}
\]

Ovo je jednadžbeni identitet, ne pretpostavka da je \(\mathcal E_a\ne0\). Definiramo \(D_a=B_a-1\). To je stroga razlika od ranije E22 shematske oznake \(c_a\), koja nije eksplicitno razdvajala potpuni potencijal, trenutnu brzinu i ubrzanje. Koeficijent \(D_a\) vrijedi samo uz gore navedene definicije \(s_a,f_a^{\rm evo}\) i opažajnu konvenciju.

## 2. Dva različita fizikalna slučaja i zabrana dvostrukog brojanja

**Geodezijska gibanja u punom metric wakeu.** Ako se u linearnom point-tracer limitu \(\Psi_{\rm tot}=\Psi_{\rm bg}+\Psi_{\nu}^{w}\) i galaksija giba po Eulerovoj jednadžbi za **taj ukupni** potencijal, tada je \(\mathcal E_a=0\), premda \(\Psi_{\nu}^{w}\ne0\) i neutrini mogu promijeniti \(V_{a,\parallel}\). U tom slučaju nema zasebnog \(\mathcal E/\mathcal H\) doprinosa. **To NE poništava** promijenjeni \(D_a V_a\), standardne potencijalne redshift/light-cone članove ni cjelokupni observable. Gravitacijski wake već jest u metric potentialu.

Ako umjesto toga u Eulerovu ostatku prešutno koristimo samo \(\Psi_{\rm bg}\), dobivamo

\[
\mathcal E_{a,\rm bg}=
V_a'+\mathcal H V_a+\partial_\parallel\Psi_{\rm bg}
=\mathcal E_{a,\rm tot}-\partial_\parallel\Psi_{\nu}^{w}.
\]

Za geodezijski \(\mathcal E_{a,\rm tot}=0\), prividni \(\mathcal E_{a,\rm bg}=-\partial_\parallel\Psi_{\nu}^{w}\) **nije nova sila koju možemo dodatno pribrojiti** uz \(\partial_\parallel\Psi_{\rm tot}\). Izvorna opažajna jednadžba eksplicitno vraća \(+\partial_\parallel\Psi_{\nu}^{w}/\mathcal H\), pa se ta dva člana egzaktno poništavaju. To je E25 fizički stop protiv lažne pojačane amplitude iz dvostrukoga brojanja.

**Fizički zaseban efektivni halo/galaxy Eulerov ostatak.** Halo ima konačnu veličinu, vlastitu povijest i okolinski uvjetovanu tracersku populaciju. Za grubo usrednjeno središte haloa ili galaktički velocity-bias field može se pojaviti nenulti efektivni \(\mathcal E_{a,\rm tot}\) u odnosu na odabranu large-scale metric-fluid jednadžbu. Ali on se mora neovisno izvesti iz prostornoga prosjeka gravitacijske sile, tidalnih/halo-satellitnih i selection uvjeta i **ne smije se pretpostaviti samo zato što postoji neutrinski wake**. Ne implicira automatski kršenje principa ekvivalencije pojedinačnih čestica. Ova granica konzistentna je s kontekstom [Okoli, Scrimgeour, Afshordi i Hudson (2017)](https://doi.org/10.1093/mnras/stx560), koji izvode gravitacijsko trenje *haloa* i ističu zajedničku wake okolinu, ali njihove FD halo-model aproksimacije nisu originalno F± kalibrirano \(\mathcal E_{a,\rm tot}\).

## 3. Minimalni, dimenzijski ispravan fizički prijenos prema E24 \(\chi_F\)

Neka je \(U_{\nu c,F}\equiv v_{\nu c,F}^{\rm km/s}/c_{\rm light}\) bezdimenzijska relativna LOS brzina. Samo ako neovisni fizikalni izvod u istom k/z/LOS režimu i za istu high-z selekciju opravdava linearne koeficijente

\[
\delta V_{a,\parallel,F}=\beta_{a,F}U_{\nu c,F},\quad
\frac{\delta\mathcal E_{a,F}}{\mathcal H}
=\epsilon_{a,F}U_{\nu c,F},\quad
\delta\Delta_{a,F}^{\rm selection}
=s^{w}_{a,F}U_{\nu c,F},
\]

tada je ukupni minimalni signed response za ovaj **izdvojeni lokalni** kanal

\[
d_{a,F}=D_a\beta_{a,F}+\epsilon_{a,F}+s^w_{a,F},
\qquad
c^{\rm E21}_{a,F}
=\frac{d_{a,F}}{c_{\rm light}}\quad[(\mathrm{km/s})^{-1}].
\]

\(s^w\) ovdje znači *neutrinski wind-conditioned selection response*, a nije standardni magnification \(s_a\) iz \(B_a\). Koeficijenti \(\beta,\epsilon,s^w,D\) **nisu** izračunati ili izmjereni za originalne high-z LRG/ELG i potencijalno mogu biti nula. Kod geodezijskoga *punog* potencijala \(\epsilon=0\), ali \(\beta\) i drugi metrčki terms ne moraju biti nula. Vremenski retarded \(\beta\) mora biti izveden iz izvora poput E17D1 s stvarnom halo+env poviješću, ne iz jednoga statičkog NFW snapshota.

Originalni E12/E13/E21 spektar je u \((\mathrm{km/s})\,\mathrm{Mpc}^3\):

\[
P_{\delta_{cb},v_{\nu c,\parallel}^{\rm km/s}}(K,\mu)
=i\mu\,C_{v,F}(K).
\]

Za bezdimenzijsku brzinu \(\,C_{U,F}=C_{v,F}/c_{\rm light}\,\), pa se E24 kontrast sad može zapisati bez skrivene pretvorbe jedinica:

\[
\boxed{
\operatorname{Im}P_{LE,F}^{\rm local,w}
=\mu\,[A_Ld_{E,F}-A_Ed_{L,F}]
\frac{C_{v,F}}{c_{\rm light}}.}
\]

\(A_L,A_E\) su realni μ-even *baseline* responzi, koji u stvarnom full-number-count modelu još sadrže baseline RSD i druge dopuštene članove. Ovo je samo lokalni source-to-odd prijenos, a **ne** cijela opažena relativistička cross-korelacija. Treba dodatno uključiti promjenu metric potencijala \(\Psi_\nu^w\) u drugim source i integrated terms, standardne GR Doppler, gravitation/evolution/magnification, wide-angle i even→odd survey leakage u originalnoj LRG→ELG konvenciji.

Kaiserova samo-derivacijska brzinska korekcija \(-\partial_\parallel V_\parallel/\mathcal H\) daje drugi imaginarni Fourierov faktor i zato, za originalni \(P_{\delta v}=i\mu C_v\), realni **parni** \(\mu^2\) cross. Nediferencirani \(D_a V_\parallel\) ili neovisni opažajni \(\mathcal E/\mathcal H\) mogu davati imaginarni odd uz specificirane fizičke koeficijente. Niti nenulti \(\beta\) sam po sebi jamči nenulti \(\chi_F=A_Lc_E-A_Ec_L\): dva različita odziva mogu se u međutracerskom kontrastu poništiti.

## 4. Egzaktne provjere i fizikalno ograničenje rezultata

[E25 protokol](../source_data/eboss_dr16_a03_e25_total_metric_euler_residual_and_linear_los_response_protocol_2026-09-29.json), Git blob **01c8449bf7388b2391b53dd5c7ba069a23c2186c**, zaključen je *nakon* što su E24 teorija i relativistička literatura pregledani, ali *prije* nove E25 egzaktne aritmetike. [E25 standard-library source-only kod](../scripts/audit_eboss_dr16_a03_e25_total_metric_euler_residual.py) provjerava šest originalnih immutable E17D1/E21/E22/E23/E24 izvora, originalna 3 F × 4 long-K E21 brzinska polja i egzaktne racionalne toy slučajeve za oba Eulerova granična režima, Kaiser/LOS parity, LRG/ELG reversal i \(c_{\rm light}=299792.458\,\mathrm{km/s}\).

[CI 36594651403](https://github.com/dvlahek/stress-energy-closure/actions/runs/36594651403) **SUCCESS**, deset negativnih kontrola. Na geodezijskom racionalnom primjeru nenulti \(\partial\Psi_w\) daje \(\mathcal E_{\rm tot}=0\), ali \(\mathcal E_{\rm bg}=-\partial\Psi_w\), te se potonji višak poništava nakon vraćanja punoga potencijala. Drugi, namjerno **sintetički**, slučaj s neovisnim nenultim ostatkom pokazuje njegovo pravilno pojavljivanje u brojnosti. Nije izveden stvarni halo non-geodesic force, nisu obrađeni full Einstein–Vlasov metric transfer ili kompletnost gauge/observer terms na nekoj realizaciji. [Trajni E25 izvorno izvedeni QA sažetak](../source_data/eboss_dr16_a03_e25_source_only_total_metric_euler_and_unit_response_summary_2026-09-29.json) zapisuje sintetički dokaz i CI provenance, ne tvrdi da su bajtovi originalnog Actions ZIP-a ponovno preuzeti i hashirani.

**Fizički STOP i sljedeći pravi izvod:** U zasebnom forward modelu za originalni \([.9,1)\) sample neovisno derivirati punu, state-specific \(\Psi_{\nu,F}^{w}(x,\eta)\) iz stvarnog zajedničkog halo+environment izvora, usklađeni \(\delta V_{a,F}\) i \(\mathcal E_{a,F}\) u *istom* potencijalu i gaugeu, te magnification/evolution/wind-selection \(s_a,f_a^{evo},s^w_a\). Tek tada iz njih slijede \(\beta,\epsilon,d,\chi\), puni z/k i fizički \(\xi_{\ell=1,3}\). E17D1 dokazuje zašto bez halo history to još nije identificirano; E23 zašto published broader-sample HOD i z_eff ne mogu odigrati ulogu originalne high-z selekcije. Ako ne možemo zatvoriti taj model bez opaženoga odd, zadržati strogo teorijski rezultat i ne proizvesti fake eBOSS template. **Observed odd SEALED, A03 PHYSICAL_UNCERTIFIED, A04 BLOCKED, bez WSL/CLASS/FITS/ASDF/Abacus/mocks, cut/seed retuninga, kontakta autora ili main/PR mergea.**
