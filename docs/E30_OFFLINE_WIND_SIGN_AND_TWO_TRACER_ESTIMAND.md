# E30 — uvjetna ravnoteža vjetra i dvije različite odd strukture

**29. 9. 2026. · Offline matematička sinteza i nova precizacija uvjeta, bez novoga fizikalnog mjerenja.** Polazište su raniji [E10](https://github.com/dvlahek/stress-energy-closure/blob/audit/eboss-elg-bit8-ra-orientation-20260925/docs/EBOSS_DR16_A03_E10_PHYSICS_OF_CONDITIONAL_WIND_AND_UNCONDITIONAL_ODD_NO_GO_2026-09-27.md), [E19](https://github.com/dvlahek/stress-energy-closure/blob/audit/eboss-elg-bit8-ra-orientation-20260925/docs/EBOSS_DR16_A03_E19_CONDITIONAL_FOURIER_SHAPE_VS_MARGINAL_RR_IDENTIFIABILITY_2026-09-29.md), [E20](https://github.com/dvlahek/stress-energy-closure/blob/audit/eboss-elg-bit8-ra-orientation-20260925/docs/EBOSS_DR16_A03_E20_QUADRATIC_WAKE_MIXED_BISPECTRUM_AND_LOS_PARITY_GATE_2026-09-29.md), [E21](https://github.com/dvlahek/stress-energy-closure/blob/audit/eboss-elg-bit8-ra-orientation-20260925/docs/EBOSS_DR16_A03_E21_ORIGINAL_E13_LINEAR_LOS_CHANNEL_AND_NORMALIZATION_2026-09-29.md), [E22](https://github.com/dvlahek/stress-energy-closure/blob/audit/eboss-elg-bit8-ra-orientation-20260925/docs/EBOSS_DR16_A03_E22_TWO_PHYSICAL_CHANNELS_OKOLI_PHASE_SHARED_WAKE_AND_SIGNED_ESTIMAND_2026-09-29.md), E24–E26 i E29. Raniji dokazi ostaju raniji dokazi; ovdje je NOVO izričito uvjetovanje sign ravnoteže po halo/galaxy/pair stanju i njegova operacionalizacija kao fizikalnoga gatea. Nema novog CLASS, Abacus, FITS, mock ni observed odd rezultata.

## 1. Problem: jednako mnogo `+v` i `-v` ne određuje sredinu općeg halo-pair odziva

Neka je `S=sign(v_nucdm,LOS)` i neka `X` sadrži sve varijable koje ulaze u uvjetnu amplitudu: masu/profil i povijest haloa, lokalnu CB gustoću, okoliš, tracer occupation, položaj/separaciju, redshift, eventualno amplitudu `|v|` te selekciju. Za **samo S-odd dio** fizičkoga izvora pišemo `T_F(S,X)=S t_F(X)`. Ako je `p_F(S,X)` stvarna zajednička raspodjela u ciljanoj galaktičkoj populaciji, a `W_F(S,X)>=0` stvarna nepoznata selekcijska/pair težina, tada je točno

\[
m_F^{\rm intrinsic}={\mathbb E_{p_F}[W_F S t_F(X)]\over\mathbb E_{p_F}[W_F]}
=\mathbb E_{q_F}\big[t_F(X)\,\eta_F(X)\big],\qquad
\eta_F(X)=\mathbb E_{q_F}[S\mid X],\quad q_F\propto p_F W_F.
\]

Stoga je **dovoljan** uvjet nule `eta_F(X)=0` za gotovo svaki `X`. Za **konstantni** `t_F(X)=t_F`, dovoljno je i `E_q[S]=0`. U općem slučaju samo `Pr_q(S=+1)=Pr_q(S=-1)=1/2` nije dovoljan uvjet: `S` može korelirati s `t_F(X)` preko halo/tracer okoliša čak i kada sam izbor `W` ne ovisi eksplicitno o `S`. To je precizno ograničenje, ne univerzalna tvrdnja da u eBOSS-u korelacija postoji.

Egzaktan svjedok: `S=+/-1`, `X=S`, po jedna jednako vjerojatna grana, `W=1` i `t(X)=X`. Marginalni `S` je savršeno uravnotežen, ali `E[S t(X)]=1`. Izvorni E10 razmatrao je **istu fiksnu amplitudu** za dva suprotna smjera, pa njegova nula ostaje valjana u vlastitom eksplicitnom modelu. E19 je već dokazao da wind-untagged objavljeni randomi ne identificiraju `W_+/W_-`; sada je vidljivo da je dodatno potreban i stvarni *joint* halo/tracer-pair zakon `p(S,X)`, ne samo dva globalna broja randoma.

U fiksnome dvogranu scalar slučaju `m_F=eta_F t_F`, `eta_F=(w_+-w_-)/(w_++w_-)`. Razlika dvaju stanja općenito je `eta_- t_- - eta_+ t_+`, **ne** `eta (t_--t_+)` bez zasebno opravdane zajedničke `eta`. Ovaj identitet ne kalibrira `eta_F` ili `t_F` za eBOSS; E10 `+1σ` vjetar također je zadavan *po F stanju*, pa njegove uvjetne Fourierove amplitude nisu gotov state-matched fizički galaxy-template.

## 2. Zašto E10 nula ne isključuje linearni E21 dvotočkasti kanal

Jednakost `E[v]=0` ne implicira `E[delta_cb v]=0`. U izvornom E26 jednom adijabatskom linearnom modu gustoća i relativna brzina korelirani su preko istog primordijalnog moda; izvorni E21 cross `P_{delta,v}=+i mu C_{v,F}` je nenulti. Za izričito ograničen E21 model s realnim mu-even baseline koeficijentima vrijedi

\[
\operatorname{Im}P_{LE,F}^{\rm linear}(k,\mu,z)
=\mu\,[A_Lc_{E,F}-A_Ec_{L,F}]\,C_{v,F}(k,z)
\equiv \mu\chi_F C_{v,F}.
\]

Ovo je korelacija DVA linearna polja (`delta`, `v`), ne E10 srednja uvjetne faze preko dvije unaprijed izabrane brzine. Egzaktni centrirani toy svjedok `delta=d`, `v=d+e/2`, `d,e=+/-1` neovisni, daje `E[delta]=E[v]=0`, `E[delta v]=1`, a `E[v delta delta]=0`. Posljednja nula je primjer centrirane zajedničke centralne simetrije i slaže se s Gaussovim vodećim ograničenjem E20 za kvadratni source `S~v delta`. Nije fizikalna simulacija niti dokaz potpunog nelinearnog Gaussova/halo odd null-a. Nenulti nelinearni `B_{v delta delta}` i dalje mora proći vlastiti LOS-odd/paritetni gate.

Po E25 fizički nediferencirani LOS koeficijent može se *tek nakon* potpunog halo/galaxy izvoda zapisati `c_{a,F}=[(B_a-1)beta_{a,F}+epsilon_{a,F}+s^w_{a,F}]/c_light`. `epsilon=0` za točno geodezijsko gibanje u **ukupnom** potencijalu; isti neutrinski wake ne smije se ponovno dodati kao Eulerova sila. Svi `beta,epsilon,s^w,B_a` za naš high-z uzorak zasad su neodređeni. Kod istog frakcijskog `c_L/A_L=c_E/A_E` imamo `chi_F=0` čak i za nenulti pojedinačni odziv. Standardni samo-derivacijski Kaiserov dio je mu-even i sam ne proizvodi taj signed odd.

## 3. Strukturna granica razlikovanja F+ i F- prije A03/A04

Ako se `chi_F(k,z)` ostavi kao **potpuno slobodna state-specific funkcija**, tada za svaki čvor gdje `mu A_L C_{v,F}` nije nula možemo konstruirati, primjerice uz `c_L=0`, vrijednost `c_E=ImP_target/(mu A_L C_{v,F})` koja daje isti zadani linearni odd za oba F stanja. To je egzaktna *neidentifikabilnost neopravdano slobodnoga modela*, ne dokaz da fizikalni `chi_F` može imati proizvoljan oblik: stvaran kauzalni halo/tracer model upravo mora ograničiti tu funkciju. Ako je dopušten samo zajednički nepoznati scalar-amplitude odgovor, oblik može nositi informaciju jedino ako dva F templatea nakon istoga fizičkog prozora nisu kolinearna modulo nuisance. Stvarni test `rank([N,t_F])>rank(N)` zahtijeva eBOSS-specific predwindow `xi`, originalni 24D window, standardni Doppler nuisance i valjanu kovarijancu. Četiri izvorna duga E21 `K` čvora i E19 jedan uvjetni Fourierov čvor nisu taj test.

Ovdje postoje **tri odvojena ranga**: E26 rank jedan zajedničke adijabatske primordijalne faze na fiksnom modu; E24 odd-only jedan kontrast `chi=A_L c_E-A_E c_L` na fiksnom modu; A04 najviše rank osam centrirane 24D sample kovarijance iz originalnih devet nezavisnih mock ID-jeva. Nijedan rang ne zamjenjuje drugi. Dvije kape nisu novih devet nezavisnih svemira; originalni DESI 18D/120-mock proizvod nije eBOSS A04.

Konačno, fizički odd model mora nositi standardne GR/Doppler, wide-angle i čak→neparno survey-window doprinose **odvojeno** od EV dijela. Sintetička leakage amplituda nije izmjerena EV amplituda; numerički even→odd curenje ne uklanja E10 sign-balanced nulu zajedničkoga linearnog operatora.

## 4. Offline kontrole i stvarni status

[E30 egzaktna QA skripta](e30_offline_estimand_qa.py) s Python `Fraction` provjerava 21 stavku: uravnoteženi marginalni znak nasuprot nenultom density–wind crossu, nulti centrirani mixed third, kondicionalnu amplitudu, različite state-specific `eta`, dvotracerski null i reverziju, common-c sign duality, 24D/9-mock rang osam te even-window negativnu kontrolu. [Strojni rezultat](e30_offline_estimand_result.json) je isključivo matematički test, ne novi CLASS ili cosmological mock. `E30_OFFLINE_QA_PASS` odnosi se na aritmetiku i logiku teorema.

**Znanstveni zaključak:** originalni neponderirani eBOSS 24D može dobiti fizički nenulti EV dio samo uz nezavisno specificiran i kalibriran **(a)** stvarni uvjetni halo-pair `eta_F(X)` i source ili **(b)** stvarni linearni signed `chi_F(k,z)` te standardnu relativističku projekciju, ili **(c)** odvojeni nelinearni mixed bispectrum/LOS mehanizam. Niti jedan još nije izračunat za originalni sample. Stoga `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`, observed galaxy rows i originalni 24D odd `SEALED`. Bez novih science cuts, F± retuninga, ASDF/FITS/CLASS/WSL i main/PR mergea.
