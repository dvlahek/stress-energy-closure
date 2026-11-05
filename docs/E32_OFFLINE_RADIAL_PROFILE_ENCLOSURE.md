# E32 — isti halo percentili ne određuju u(k), ali dopuštaju konzervativnu ogradu

**29. 9. 2026. · Novi konstruktivni matematički rezultat; samo ilustrativni numerički primjeri.** [Službena AbacusSummit specifikacija](https://github.com/abacusorg/AbacusSummit/blob/4b1959c710cb0c49aa305c6213a228aa2a4587ff/docs/data-products.rst) navodi L2 masene percentile 10,25,33,50,67,75,90,95,98 i 100% te sirovo pohranjene omjere prema `r100_L2com`. Ovo NIJE `M200c` profil. Naših 20 dvovremenskih lokalnih PID parova nema pročitane, dekodirane vrijednosti tih percentila. Nema novih ASDF čitanja; Ubuntu se ranije srušio.

## 1. Problem fizičkoga E29 izraza

U E29 vanjskom Bornovu odgovoru nije dovoljno znati ukupni `N` i centar. Za evoluirajući halo treba `u_test(k,z_obs) integral dz' [M(z')/M_obs] u_src(k,z') K_F(k,z')`. Čak i uz savršene iste-object L2 percentile u dva snimka, raspodjela mase između spremljenih radijusa nije zadana. Na kratkim skalama ta je raspodjela važna. Osim toga, L2 percentili nisu puni L1 `M200c`/environment, pa ni njihov eventualno idealan u(k) ne zatvara fizički E28 drag.

## 2. Konstruktivna nejednoznačnost čak za SVIH deset percentila

Za nenegativnu sferno usrednjenu distribuciju mase ukupne mase `M`, neka je `Q(p)` radijus koji obuhvaća udio `p` mase i neka je `p_i=(.10,.25,.33,.50,.67,.75,.90,.95,.98,1)`. Pozitivna masa u percentilnom binu je `w_i=p_i-p_(i-1)` i radijus joj leži između `r_(i-1)` i `r_i=Q(p_i)` (uz `r_0=0`). Sferni Fourierov profil glasi

\[
u(k)={1\over M}\int_0^\infty \operatorname{sinc}(kr)\,dM(r),\quad
\operatorname{sinc}(x)={\sin x\over x}.
\]

Stvorimo model A s cijelim `w_i` na desnom rubu `r_i`. U modelu B polovica svakoga `w_i` stoji na sredini binova `(r_(i-1)+r_i)/2`, a druga polovica na `r_i`. U oba je **svih deset** `Q(p_i)` egzaktno jednako `r_i`, uključujući `r100` i ukupnu masu, no `u_A(k)` i `u_B(k)` općenito se razlikuju. Ovo je matematički svjedok s pozitivnim sfernim ljuskama, ne model za stvaran Abacus halo.

Za ilustrativne bezdimenzijske `r_i=i/10` i `r100=1`, QA daje:

| `q=k r100` | `u_A(q)` | `u_B(q)` | `|u_A-u_B|` |
|---:|---:|---:|---:|
| 1 | 0.957256 | 0.960730 | 0.003475 |
| 3 | 0.679086 | 0.702689 | 0.023602 |
| 8 | 0.143809 | 0.176547 | 0.032738 |
| 20 | 0.030933 | 0.064369 | 0.033436 |

Ove vrijednosti nemaju jedinice halo-ubrzanja, nisu `k=8 Mpc^-1` E28 slučaj i ne pokazuju fizičku UV divergenciju. Dokazuju samo nepotpunost percentilne informacije o Fourierovu profilu.

## 3. Što ipak možemo strogo ograničiti iz percentila

Identitet `sinc(x)=integral_0^1 cos(tx)dt` daje `|sinc'(x)|<=1/2` za svaki realni `x`. Neka su `m_i=(r_(i-1)+r_i)/2` i `U_mid(k)=sum_i w_i sinc(k m_i)`. U svakom binu pravi radijus odstupa od sredine najviše `(r_i-r_(i-1))/2`. Zato je za **sferno usrednjeni pozitivni** profil dokaziva ograda

\[
\boxed{\left|u(k)-U_{\rm mid}(k)\right|
\le {\lvert k\rvert\over 4}\sum_i w_i (r_i-r_{i-1})}. 
\]

Interval se dodatno može presjeći s `[-1,1]`. To je sigurna, ali ne nužno oštra ograda; ne zahtijeva pretpostavku NFW-a ili interpolaciju nepoznate mase. Numerička širina raste barem u ovoj konzervativnoj formuli s `|k|`, pa isti percentili mogu slabo ograničavati kratke valove. Za ilustrativni `q=8` kod nalazi interval `[0.009286,0.409286]`, koji obuhvaća oba gore konstruirana profila; pri `q=20` dopušten je i predznak obaju smjerova (interval približno `[-0.4022,0.5978]`).

**Primjenjivost:** ako kasnije sigurno očitamo stvarne L2 percentile, možemo bez profilnog fitanja dati interval za *sferno usrednjeni L2* `u(k)`. Ograda ne daje puni anisotropni `u(vector k)`, stvarni ukupni L1 host `M200c`, same-object fizikalno dinamičku povijest ili nelinearni Born error. Pogotovo ne pretvara 94.22589% UV-osjetljivi E28 finite-band račun u fizički total force.

## 4. Kontrole i odluka

[E32 QA skripta](e32_offline_radial_quantile_qa.py) prošla je 41 uvjet, uključujući svih deset egzaktno istih percentila, različit transform, pozitivnost, intervalne ograde i namjerno pokvareni prvi percentil. Izvorni [E32 JSON](e32_offline_radial_quantile_result.json) dokumentira samo ilustrativnu konstrukciju. Nema razloga ponovno otvarati golema ASDF polja samo zato da dobijemo `r10/r50/r90`: ona ne određuju stvarni E29 profil. `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`, observed eBOSS odd `SEALED`.
