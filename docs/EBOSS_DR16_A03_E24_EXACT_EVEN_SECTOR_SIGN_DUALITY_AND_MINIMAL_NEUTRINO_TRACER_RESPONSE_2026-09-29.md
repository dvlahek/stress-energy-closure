# A-03E24 — parni clustering ne određuje predznak linearnoga neutrinskog LOS odziva

**29. 9. 2026. · Novi izričito ograničeni, egzaktni teorijski rezultat.** E21 je iz originalnog E13 CLASS direct-vTk dao nenulti linearni izvor \(P_{\delta_{cb}v_{\nu c,\parallel}}\) na četiri zamrznuta duga K, ali nije izračunao kako LRG i ELG galaksije fizički reagiraju na relativni vjetar. E22 razdvaja hipotetski linearni opažajni Doppler/selection odziv od kvadratnog \(\Gamma v\delta\) wakea. E23 je ponovno potvrdio izvorni high-z \([0.9,1.0)\) eBOSS RR prozor i odbio preslikavanje HOD/redshift rezultata širih eBOSS uzoraka na taj high-z 24D pilot. E24 odgovara na uži problem: **može li standardni parni auto/cross clustering, čak i kada bi bio savršeno modeliran, odrediti predznak nedostajućega neutrinskog responsea?**

## 1. Egzaktni dvopoljni model i konvencija

Definiramo samo transparentnu *restriktivnu* determinističku komponentu stvarnih opaženih polja za svaki F, k, z i LOS \(\mu\):

\[
\delta_a^s=A_a(k,\mu,z)\,\delta_{cb}+c_{a,F}(k,z)\,v_{\nu c,\parallel,R16,F},
\qquad a=L,E.
\]

\(A_a\) su realni i parni u \(\mu\); \(c_{a,F}\) su realni i parni u \(\mu\), ali fizički **NEPOZNATI**. Ovo nije standardno Kaiserovo automatsko izdvajanje: samo \(-\mathcal H^{-1}\partial_\parallel v_\parallel\) daje parni \(\mu^2\). Nefazni, nediferencirani \(c_a v_\parallel\) bi trebao fizički izveden relativistički Doppler/selection ili drugi dopušteni opažajni mehanizam iz stvarne halo–galaxy evolucije. E24 ne tvrdi da on u originalnim eBOSS galaksijama postoji.

Izvorni E12/E13 daje \(\langle\delta_{cb} v_{\nu c,\parallel}^{*}\rangle=i\mu C_{v,F}\) i \(\langle v_{\nu c,\parallel}\delta_{cb}^{*}\rangle=-i\mu C_{v,F}\). \(P_{\delta\delta}\) i \(P_{vv}\) su realni i parni u \(\mu\), a \(P_{vv}\ge0\). S originalnom LRG→ELG orijentacijom, direktna, BEZ aproksimacije po veličini \(c\), algebra unutar ovoga dvopoljnog modela daje

\[
\begin{aligned}
P_{LL}&=A_L^2 P_{\delta\delta}+c_L^2P_{vv},\\
P_{EE}&=A_E^2 P_{\delta\delta}+c_E^2P_{vv},\\
\operatorname{Re}P_{LE}
 &=A_LA_EP_{\delta\delta}+c_Lc_EP_{vv},\\
\operatorname{Im}P_{LE}
 &=\mu(A_Lc_E-A_Ec_L)C_{v,F}.
\end{aligned}
\]

To nisu puna relativistička opažena brojnost ni potpuni eBOSS odd: standardni gravitacijski redshift, magnification/evolution, wide-angle, stochastic, selekcijska parnost i drugi velocity–density korelatori ovdje nisu uključeni. Dodavanje dodatnih kompleksnih baseline \(A_a\) ili izravno koreliranih polja može promijeniti zaključak.

## 2. Nova egzaktna neidentifikabilnost: globalni predznak

Unutar gore specificiranoga modela **simultano** \(c_L\to-c_L,\ c_E\to-c_E\) čuva sva tri parna spectra, čak i njihove članke reda \(c^2P_{vv}\), jer sadrže isključivo \(c_L^2,c_E^2,c_Lc_E\). Istodobno \(\operatorname{Im}P_{LE}\to-\operatorname{Im}P_{LE}\). To je točna dvostruka sign-degeneracija fizički nepoznatoga neutrinskog odziva iz **samo parnoga** 2pt skupa unutar ograničenog modela. Nije riječ samo o Gaussovu limitu, numeričkoj toleranciji ili premalom broju mockova. Parni auto/cross može pri redu \(c^2\) pridonijeti informaciji o apsolutnim vrijednostima ako su svi drugi inputi poznati, ali ne razbija simultanu sign-degeneraciju.

Pri prvom redu u \(c\), \(\partial_{c_a}P_{LL}\), \(\partial_{c_a}P_{EE}\) i \(\partial_{c_a}\operatorname{Re}P_{LE}\) svi su nula u \(c=0\). Neparni \(\operatorname{Im}P_{LE}\) zadržava linearni kontrast

\[
\chi_F(k,z)=A_L c_{E,F}-A_E c_{L,F}.
\]

Jedan LE odd po fiksnom k,z,\(\mu\) može identificirati samo \(\chi_F\) (uz poznati \(C_{v,F}\), ostale nuisance i validnu statistiku), ne oba koeficijenta pojedinačno. Pomak \((c_L,c_E)\mapsto(c_L,c_E)+t(A_L,A_E)\) točno čuva **odd kontrast**, ali općenito mijenja parne članove \(O(c^2)\); ne tvrditi da je to neidentifikabilnost kombiniranoga *punog* čak i visokopreciznog parnog+odd skupa.

## 3. Zašto objavljeni HOD/velocity bias nije taj kontrast

[Alam et al., MNRAS 497, 581 (2020)](https://doi.org/10.1093/mnras/staa1956) uklapaju multitracer HOD preko širih \(0.7\le z\le1.1\) eBOSS uzoraka. [Ávila et al., MNRAS 499, 5486 (2020)](https://doi.org/10.1093/mnras/staa2951) za DR16 ELG razdvajaju satelitsku disperziju, velocity bias i mogući halo-radijalni infall. [Yuan et al., MNRAS 510, 3301 (2022)](https://doi.org/10.1093/mnras/stab3355) eksplicitno razlikuju centralni/satelitski *halo-relative* velocity-scatter parametar od occupation/secondary halo biases te koriste eBOSS sample \(0.7\le z\le1.1\). Nijedan takav HOD/velocity-scatter parametar nije iz definicije \(\partial\delta_a^{s}/\partial v_{\nu c,\parallel}\) uz fiksirane gravitacijske, selekcijske i gauge uvjete. Parni RSD/satelitska disperzija i koordinatni Kaiserov velocity bias ne daju automatski \(\chi_F\).

Za high-z \([.9,1)\) fizički \(\chi_F\) formalno zahtijeva izvorno specificiran \(p_a(M,c,\mathrm{assembly},\mathrm{env},\mathrm{cen/sat}\mid z,\mathrm{cap},\mathrm{selection})\), vlastiti neutrinski halo velocity response i **opažajni** velocity/selection koeficijent. Iako se pri determinističkom long-mode limitu može agregirati u jedan \(\chi_F\), konačni međutracerski bin ima parno/okolinski uvjetovan zajednički wake; zasebni marginalni HOD-ovi sami ne definiraju pun \(p_{LE}\) parova ni wind–selection korelaciju. Izvorni E19 dokaz marginal-RR neidentifikabilnosti nije razriješen HOD srednjim masama.

Kao zaseban egzaktan **čisto matematički** HOD-joint svjedok, dva 2×2 nenegativna modela s marginals \(p_L=p_E=(1/2,1/2)\) daju \(\mathbb E[XY]=+1\) za zajedničke razrede ili \(-1\) za suprotne razrede \(X,Y=\pm1\). Dakle ni identični *jednotracerski* halo marginals ne specificiraju parni okoliš. Ova četiri diskretna stanja nisu stvarne mase eBOSS haloa, nisu validirani high-z 2pt modeli i ne dokazuju da izmjereni standardni cross-clustering nema nikakvu informaciju o tom paru.

## 4. Izvorna QA i operativna odluka

[E24 protokol zaključen nakon originalnoga E21/E23 znanja, ali prije NOVIH E24 sintetičkih testova](../source_data/eboss_dr16_a03_e24_even_sector_linear_wind_response_sign_identifiability_protocol_2026-09-29.json), originalni Git blob \`3fd79b8cfa116b05c5ab55a248f0f95f2ec7a45b\`, razdvaja posthoc znanstveni izvod od unaprijed testiranih egzaktnih aritmetičkih kontrola. [E24 standard-library program](../scripts/audit_eboss_dr16_a03_e24_even_sector_wind_sign_identifiability.py) SHA-pina originalni E21/E22/E23, ponovno provjerava sva tri originalna F i četiri duga E21 K bez novog CLASS-a, te sa strogo racionalnim parametrima provjerava egzaktne LE/EL i LOS predznake, \(O(c)\) even-null, \(O(c^2)\) čak i uz globalnu promjenu predznaka, rank-null i toy halo-pair witness.

Prvi CI 36583036306 **FAIL CLOSED iz čisto tehničkoga razloga**: novi E24 čitač tražio je ASCII \`M200c\` u originalnom E23 tekstu koji uredno piše LaTeX \`M_{200c}\`; originalni E23 Git SHA, E21/E22, znanstveni izvod i protokol nisu mijenjani. Samo source-reader keyword promijenjen je na stvarni originalni literal. Završni [E24 CI 36583153519](https://github.com/dvlahek/stress-energy-closure/actions/runs/36583153519) **SUCCESS**, devet negativnih kontrola i svih originalnih 3×4 izvora; originalni mali CI zip s egzaktnom QA nije ponovno byte-hashiran za ovaj dokument. [Trajni E24 rezultat](../source_data/eboss_dr16_a03_e24_source_only_even_odd_linear_response_identifiability_summary_2026-09-29.json) izričito je izvorno izvedeni sažetak i provenijencija, a ne reprodukcija bajtova toga zipa.

Znanstveni idući cilj sada je **neovisno identificirati \(\chi_F(k,z)\) ili fizički izvedene pojedinačne \(c_{a,F}\)** za originalni high-z LRG×ELG selection, primjerice iz neovisno specificiranoga relativno-brzinskoga halo+galaxy responsea uz učvršćenu relativističku number-count projekciju. Ako je \(\chi_F=0\), ovaj ograničeni linearni odd mehanizam je nula premda \(P_{\delta v}\ne0\). Ako je \(\chi_F\ne0\), treba pun originalni k/z/LOS/window/nuisance i A04 inferencijsku kovarijancu prije ijedne empirijske tvrdnje. E20 kvadratni kanal još zahtijeva zaseban puni \(B_{v\delta\delta}\), \(\Gamma_i\) i LOS. **Ne otvarati originalni observed 24D odd, ne prenositi broad-sample HOD, ne uvoditi cut/seeds, nove ASDF/FITS/mocks/CLASS/WSL, kontakt autora ili main merge.**
