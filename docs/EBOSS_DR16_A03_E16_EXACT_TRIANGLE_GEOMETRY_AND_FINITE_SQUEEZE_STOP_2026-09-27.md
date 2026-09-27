# A-03E16 — zatvoren Fourierov trokut i granica reducirane squeezed aproksimacije

**27. 9. 2026. · Originalni E16 i neovisni matematički audit PASS. Puni fizički eBOSS bispektar nije identificiran.** Ovaj korak slijedi E15, koji je provjerio samo reducirani `B_source/(i mu_long Delta b)=S(k_short,K_long,ell_short)`. Originalni E8–E15 izvori, uključujući svih 72 E14 slučaja i svih 24 kontrasta `F_+−F_-`, ostaju izvorno SHA-pinned i trajno arhivirani na audit grani. Opaženi odd je SEALED.

## Problem koji smo zatvorili

E14/E15 pohranjuju `k_short` i `K_long`, ali *ne* određuju dva stvarna kratka vektora unutar trokuta. Zbog toga jedan imaginarni faktor `i mu_long S` nije jedinstveni bispektar `B(k1,k2,K)`. Za fizikalni forward model trebaju oba kraka, kut među njima, azimut i kutove prema liniji gledanja, linearna polja `P_cb(k1),P_cb(k2)` i fizikalno izvedena sprega dugog s kratkim modovima.

E16 [protokol zaključan prije računanja](../source_data/eboss_dr16_a03_e16_exact_triangle_geometry_finite_squeeze_prereg_2026-09-27.json) ne bira novi izvor `F_\pm`, nove redshift/galaxy rezove ni numeričke acceptance pragove prema rezultatu. Polazi iz originalnog E14 full SHA256 `0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8` i E15 full SHA256 `6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983`. Predeclared originalna mreža ima `k=(.05,.075,.1) h/Mpc`, `K=(.001,.002,.003,.005) h/Mpc`, `ell_short=(1,3)`, tri stanja i ilustrativni `z=.95`. Nema novoga CLASS računanja, FITS podataka ili izmjene izvornog eBOSS estimatora.

## Egzaktna 3D geometrija

Za `n=(0,0,1)` i srednji kratki vektor `k` definiramo stvarne krakove

\[
\mathbf k_1=\mathbf k-\frac{\mathbf K}{2},\qquad
\mathbf k_2=-\mathbf k-\frac{\mathbf K}{2},\qquad
\mathbf k_1+\mathbf k_2+\mathbf K=0.
\]

Neka su `r=K/k≤0.1`, `mu_s=khat·n`, `mu_L=Khat·n`, `phi` relativni azimut, te

\[
c=\widehat{\mathbf k}\cdot\widehat{\mathbf K}
=\mu_s\mu_L+\sqrt{1-\mu_s^2}\sqrt{1-\mu_L^2}\cos\varphi.
\]

Tada, bez squeezed truncation,

\[
\frac{k_{1,2}^2}{k^2}=1+\frac{r^2}{4}\mp rc,\quad
\mu_1=\frac{\mu_s-r\mu_L/2}{\sqrt{1+r^2/4-rc}},\quad
\mu_2=\frac{-\mu_s-r\mu_L/2}{\sqrt{1+r^2/4+rc}}.
\]

Razlikovati `mu_s`, `mu_L`, `mu_1`, `mu_2` i ranije projicirani `ell_short`; ni jedan od njih nije automatski originalni E15 dugi LOS multipol `L`.

E16 [izvršivi originalni kod](../scripts/audit_eboss_dr16_a03_e16_exact_closed_triangle_geometry.py) provjerio je 48 **prije zaključanih** trodimenzijskih konfiguracija na svakoj od 12 kombinacija `(k,K)`, ukupno 576. Provjera sadrži egzaktno `k1+k2+K=0`, oba originalno orijentirana kratka LOS kuta i zamjenu kratkih krakova. Originalni [E16 full source JSON](../source_data/eboss_dr16_a03_e16_original_72_source_exact_triangle_geometry_2026_09_27.json), SHA256 `5d81f95a5571175f06eb032b865285ba34804701693f76af392be7fc31fb06db`, trajno je spremljen na GitHub [CI 36339038115 PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36339038115).

## Što konačni `K/k` sigurno mijenja, a što još ne znamo

Samo za *poznati* kvazistatički E14 faktor `1/k^2`, zamjena središnjega kratkoga vektora stvarnim krakom daje čisto geometrijski omjer

\[
G_i=\frac{k^2}{k_i^2}
=\frac{1}{1+r^2/4\mp rc},
\qquad
\frac{1}{(1+r/2)^2}\leq G_i\leq\frac{1}{(1-r/2)^2}.
\]

Pri najvećem originalnom `r=.1` pojedinačni omjer leži između `0.907029478...` i `1.108033241...`, pa maksimalno *pojedinačno* odstupanje iznosi **10.8033241%**. To **nije** error budget niti granica pogreške stvarnog bispektra: mogući fizički doprinosi dvaju krakova, njihova relativna težina, dinamički Vlasov kernel, linearni `P_cb(k_i)`, tracer coupling i gauge/LOS termini još nisu izračunani.

E14 pohranjuje samo tri izvorne CLASS `P_cb(k)` vrijednosti po stanju, na `k=.05,.075,.1 h/Mpc`. E16 dokazuje da u **8 od 12** originalnih `(k,K)` parova neki zatvoreni kratki krak izlazi iz toga pohranjenoga raspona. Na rubnim `k=.05` i `k=.1` to se događa za svako originalno `K>0`. Interpolacija ili izmišljena ekstrapolacija iz tri točke ne zatvara ovaj fizikalni ulaz.

Čak i kada je svaki geometrijski krak poznat, reduciran `S` ne određuje njegovu fizičku amplitudu. Dva jednostavna **isključivo matematička** dovršetka imaju isti originalni E15 squeezed limit i isti Hermitian reality pod reverzijom svih vektora:

\[
\frac{B_A}{i\mu_L\Delta b}=S,\qquad
\frac{B_B}{i\mu_L\Delta b}
=S\,G_{\rm sym}(c,r),\qquad
G_{\rm sym}=\frac12\left(\frac{k^2}{k_1^2}+\frac{k^2}{k_2^2}\right).
\]

Za `r>0` funkcije se razlikuju, a za `r→0` obje daju `S`. Druga je *toy simetrizacija samog kvazistatičkog geometrijskog faktora*, nije objavljena ili kalibrirana Einstein–Vlasov dinamika. Najveća razlika tog toy simetričnoga omjera na predeclared 576 konfiguracija jest **0.753136%**, što također nije bound na stvarnu fizikalnu pogrešku. Isti ulazi, isti squeezed limit i reality dopuštaju različite konačne oblike, pa iz originalnog E14/E15 **nije moguće jedinstveno rekonstruirati** puni bispektar.

## Neovisno verificirano i trajno spremljeno

[E16 neovisni čisti Python audit](../scripts/audit_eboss_dr16_a03_e16_independent_exact_triangle_replay.py) **ne uvozi** originalnu E16 implementaciju niti koristi NumPy, CLASS ili FITS. Iz originalnih E14/E15 i SHA-zaključanog E16 ponovno izvodi svih 12×48 skalarnih geometrija, točne `G_i` granice, 72 trostanjskih izvora, 24 izvorna `F_+−F_-` kontrasta i osam slučajeva izvan originalnog `P_cb` supporta. [Neovisni CI 36339155299 PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36339155299) daje maksimalni skalirani mismatch `2.220446049250313e−16`; njegov originalni [JSON SHA256 `8a6ee7610a0412de7bcaf5333950f41add0973a000bd245d49b41b77efc8e10d`](../source_data/eboss_dr16_a03_e16_independent_original_scalar_geometry_replay_2026_09_27.json) također je trajno na audit grani. Ovo je neovisna **matematička** certifikacija istih izvornih podataka, ne nova neovisna fizikalna simulacija.

## Zaključana granica idućeg fizičkog koraka

E16 geometrija jest dovršena. **Puna fizika trokuta nije**. Za sljedeći korak treba pri originalnim kratkim i dugim osima dobiti po stanju `P_cb(k1,z)`, `P_cb(k2,z)` i njihove relevantne transfer/LOS response, izvesti kako oba kratka kraka nose originalni direct-vTk rank-conditioned wake i zasebno testirati konačni `K/k` kernel. Držanje samo pozitivnoga `+1sigma` LOS vjetra ili množenje staroga `S` geometrijskim `G_sym` nije zamjena za taj izvod.

Tek zatim smijemo modelirati LRG/ELG high-z HOD/bias, magnifikaciju/evolution/GR nuisance, empirijski NGC/SGC chunk/depth/z **triple** selection prozor i novu nezavisnu eBOSS bispektralnu kovarijancu. Originalni unconditional 24D `xi_odd` nije nova trotočkasta statistika. Originalni SGC v1 independent reverse kvar i 48k finite-random nekonvergencija ostaju odvojeni. Observed odd SEALED, nema novih mock/katalog preuzimanja, novih science seed/cut, promjene `main` ili mergeanja PR-a.
