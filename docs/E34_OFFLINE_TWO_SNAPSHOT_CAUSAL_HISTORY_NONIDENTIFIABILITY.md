# E34 — dva identična halo snimka ne određuju retardirani izvor

**29. 9. 2026. · Matematička konstrukcija i nezavisna offline QA; nije novi fizički E28 rezultat.** Nadovezuje se na E29 strukturni vanjski Bornov izraz `u_test(k,t_obs) ∫dt K_F(k,t) [M(t)/M_obs] u_src(k,t)`. Postojećih 20 lokalno PID-uparenih Abacus haloa ne daje kontinuirane masene profile, stvarni neutrinov vjetar niti početni neutrinski wake.

## 1. Problem i konstrukcija

Čak i kada bismo za **isti objekt** točno znali profil, masu i radijalnu brzinu u dvjema epohama, srednji dio njegove povijesti nije jednoznačno određen. To nije samo poteškoća mjerenja Abacusovih L2 percentila, nego matematička granica dvovremenskog forward modela.

Za dokaz uzmimo pozitivnu sfernu ljusku stalne mase `M`, fiksnog središta i radijusa

\[
r_{\pm}(t)=r_0\pm\varepsilon\sin^2(\pi t/T),\qquad
0\le t\le T,\quad 0<\varepsilon<r_0.
\]

Obje povijesti imaju iste radijuse `r_0` i nultu radijalnu brzinu u oba snimka. Radijalne brzine iznose `v_r=±(eps*pi/T) sin(2*pi*t/T)`. Za test-funkciju `phi(r)` vrijedi `d phi[r(t)]/dt=phi'[r(t)]v_r(t)`, pa transport pozitivne ljuske konzervira masu i zadovoljava slabu jednadžbu kontinuiteta. Ta kinematička konstrukcija **nije** samodosljedno riješen Poisson–Vlasov ili Einstein–Vlasov halo: potrebna sila koja ostvaruje trajektoriju nije izvedena iz vlastite gravitacije.

Sferni, masom normirani Fourierov profil ljuske je `u_±(k,t)=sinc(k*r_±(t))`. Fiksiramo istu krajnju testnu ljusku `u_test=sinc(k*r_0)`, istu masu i isti ilustrativni kauzalni kernel `K(t)=exp[-(T-t)/T]`. Integrirani rezultati su `J_±=u_test ∫_0^T K(t)u_±(k,t)dt`.

Pri **isključivo bezdimenzijskom** `T=r_0=k=1`, `eps=0.1`, vrijedi

| Ilustrativna dijagnostika | Vrijednost |
|---|---:|
| `J_+` | `0.439549571374331` |
| `J_-` | `0.455163833342087` |
| `J_- − J_+` | `0.0156142619677557` |

Predznak ne ovisi o finoj numerici: na intervalnom rasponu `r∈[.9,1.1]` funkcija `sinc(r)` strogo je padajuća, `K(t)>0`, a unutarnji pomak `g(t)>0` za `0<t<T`. Stoga `J_- > J_+` egzaktno. Simpsonov jaz između 512 i 1024 podintervala jest `1.43e−13`. Ovo nisu E28 jedinice, F± razlike, halo-drag ni opažajni predložak.

## 2. Što bi bilo potrebno za *ograničiti* pogrešku povijesti

Ako neovisni fizički model daje vremensku Lipschitzovu konstantu `L = sup_t |du_src/dt|`, za linearnu interpolaciju `u_lin` između poznatih profila vrijedi

\[
|u_{\rm src}(t)-u_{\rm lin}(t)|
\leq \frac{2L}{T}t(T-t).
\]

Za isti fiksni vanjski linearni kernel i relativnu masu `m(t)≥0`, slijedi

\[
\left|u_{\rm test}\int_0^T K(t)m(t)
[u_{\rm src}(t)-u_{\rm lin}(t)]dt\right|
\le |u_{\rm test}|\frac{2L}{T}
\int_0^T |K(t)|m(t)t(T-t)dt.
\]

Zadnja ograda pretpostavlja da `L`, `K` i `m` dolaze iz **istog** fizičkog modela. Dva snimka sama ne daju `L`. Ni naš Abacusov centarski pomak ne ograničava radijalnu promjenu profila ili raniji neutrinski wake. Za ilustrativni slučaj QA koristi siguran `|sinc'(x)|≤1/2` i `|dr/dt|≤eps*pi/T`, pa `L≤|k|eps*pi/(2T)`.

Za generički neidentički nulti **potpisani** kernel može se napraviti unutarnja perturbacija lokalizirana ondje gdje kernel ima fiksni predznak. Ipak, prikazana numerička vrijednost vrijedi samo za pozitivni toy `K`, ne za stvarni oscilatorni E28 `K_F`.

## 3. Znanstvena posljedica i status

Stvarni PID match i čak hipotetski točni endpoint profili **ne identificiraju** E29 `∫u_src K_F`. Nužan je fizikalno specificiran i dinamički dopušten kontinuirani halo+environment model ili neovisno opravdan bound na propuštenu povijest. Samo interpoliranje dvaju `M200c`/L2 percentile snimaka nije validacija. Ovaj dokaz pojačava E29 STOP, ali ne dokazuje da stvarna E28 F± razlika ima točno ovaj iznos ili predznak.

Izvorni `e34_offline_causal_history_result.json` i `e34_offline_causal_history_qa.py` sadrže 51/51 PASS, zasebne parent SHA256, negativnu kontrolu predznaka i provjeru slabe kontinuitetne relacije. `A03 PHYSICAL_UNCERTIFIED`, `A04 BLOCKED`, observed eBOSS odd **SEALED**, bez ASDF/FITS/WSL, novih simulacija, GitHub izmjena i main mergea.
