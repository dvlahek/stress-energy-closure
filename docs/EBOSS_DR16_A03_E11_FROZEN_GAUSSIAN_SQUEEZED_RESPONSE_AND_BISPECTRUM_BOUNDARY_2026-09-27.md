# A-03E11 — nenulti uvjetovani squeezed response, nulti unconditional dvotočkasti mean

**Datum i znanstveni status: 27. 9. 2026.** Originalni E8 source `F_0,F_+,F_-`, E9 odvojeni CLASS density-derived R16 neutrino–CDM proxy-vjetrovi i E10 uvjetna Fourierova Eq.(A8) angular odd amplituda daju izvorno reproducibilan *Gaussian/Wick long–short response coefficient* `Γ_ell`. Taj je koeficijent **nenult** za svaki originalni state, ali strogo simetrični wind model i dalje daje nulti neuvjetovani srednji intrinzični dvotočkasti odd. E11 nije opaživi eBOSS bispektar, galaktički `xi_1/xi_3` ni detekcija.

## Zašto je potrebno prije sljedećeg eBOSS testa

U [Okoli et al. (2017), Eq. (A7)–(A12)](https://doi.org/10.1093/mnras/stx560) lokalno uniforman neutrino vjetar stvara imaginarni dvotracerski Fourierov cross-power. Njegov predznak mijenja se pod `v_nu c → −v_nu c`. Prethodni [E10 izvorno zatvoreni test](EBOSS_DR16_A03_E10_PHYSICS_OF_CONDITIONAL_WIND_AND_UNCONDITIONAL_ODD_NO_GO_2026-09-27.md) izračunao je uvjetni `T_ell(r=1)`, ali dokazao i `E[T_ell]=0` za simetrični vjetar kada selekcija ne preferira njegov predznak. Zato ne smijemo preslikati E9 `+1sigma` predznak u originalnu unconditional 24D `xi_odd`.

[Zhu i Castorina, PRD 101 (2020) 023525](https://doi.org/10.1103/PhysRevD.101.023525) pokazali su da neutrino–CDM relativna brzina može ostaviti antisimetričan bispektar koji povezuje različita kozmička polja. [Nascimento i Loverde, JCAP 04 (2026) 019](https://doi.org/10.1088/1475-7516/2026/04/019) također naglašavaju zahtjev za realističnim modelom i odgovarajućim tracerima. **Naš E11 nije rekonstrukcija njihovih estimatora**: ovdje računamo samo vlastiti source-only response koeficijent predanoga modela, koji je nužan, ali ne i dovoljan uvjet za ostvariv bispektar.

## Zamrznuta fizika i jednadžba

[Prospektivni E11 protokol](../source_data/eboss_dr16_a03_e11_gaussian_wind_squeezed_response_prereg_2026-09-27.json) zaključan je **nakon E10 rezultata, prije E11 Gaussian quadrature**. SHA originalnog E8 4000q CSV-a `bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0`, originalnog E9 tri-CLASS JSON-a `450b35f68a329f545f21bbe7ec02e36941e4c303037ccb7c011830c315d15890` i originalnog E10 angular JSON-a `cc55a914a8be0174eb7c1e49e0e1968537cedfcf4a429d52ad66ab86cb44a488` provjereni su **bajt-po-bajt prije pokretanja**. Nije bilo ponovnog CLASS-a, preuzimanja kataloga/mokova, promjene izvornoga `F_\pm` ni nove optimizacije.

Neka je `r~N(0,1)` **isti koherentni Gaussian rank** na svim trima originalnim stanjima i redshift čvorovima. Originalni E9 density-transfer-derived proxy `sigma_LOS,s(z)` pomičemo u `v_parallel,s(z,r)=r sigma_LOS,s(z)`. To je linearno-Gaussian **model vjetra**, a ne izmjereno halo polje ili dokaz da njegov density-derived proxy točno odgovara CLASS vTk neutrino brzini. Originalni E10 E8-occupancy/k-scaling ostaje točno isti:

\[
\phi_s(k_{\rm com},z,\mu,r)=
\phi_s(k_{\rm com},z,\mu=1,r=1)
\;r\mu\;
\frac{f_s\!\left(q_{\ast,s}(z)\,|r\mu|\right)}
     {f_s\!\left(q_{\ast,s}(z)\right)},\qquad
q_{\ast,s}(z)=\frac{m_\nu\sigma_{\rm LOS,s}(z)}
{cT_{\nu0}(1+z)}.
\]

`f_{\rm FD}(q)=1/(e^q+1)`, `f_\pm(q)=f_{\rm FD}(q)\times F_\pm(q)/F_0(q)`. Potpuno smo zadržali E8 4000q interpolaciju i originalni točni Eq.(20) Fermi–Dirac coefficient, **ne** nerešen `mu=1` approximation iz tiskane literature.

Za originalni `k_com=.05h/Mpc`, matematički midpoint `z=.95` (ne measured eBOSS `z_eff`), `z` derivacija `±.005` i zamrznuti velocity-source `0.004(1+z)`:

\[
T_{\ell,s}(r)=\frac{2\ell+1}{2}\int_{-1}^{1}
d\mu\;P_\ell(\mu)\,\mu^2\,
\frac{d\phi_s(k_{\rm com},z,\mu,r)}{d\ln a},
\quad \ell=1,3,
\]

\[
\boxed{\Gamma_{\ell,s}=\mathbb E_{r}[rT_{\ell,s}(r)]
=\mathbb E_{r}\!\left[\frac{dT_{\ell,s}(r)}{dr}\right]},
\quad r\sim{\cal N}(0,1).
\]

Druga jednakost jest Gaussian integration-by-parts/Stein identity, numerički provjerena pri zaključanim `Δr=.002,.001`. Pod **dodatnom** pretpostavkom da stvarni dugi mod `D_L(\mathbf K)` i naš lokalni `r` imaju zajedničku Gaussian linearnu kovarijancu i dobro definiranu zajedničku geometriju, dobivamo *schematic* squeezed source response `P_{D_L,r}(K,z)\,\Gamma_{\ell,s}`. Faktor `P_{D_L,r}` može biti imaginarnog/longitudinalnog predznaka `i\hat{\mathbf K}\cdot\mathbf n`, a u E11 **nije** izračunan niti opažen. Taj proizvod dodatno zahtijeva odgovarajuću realizaciju short-mode cross-powera i odgovore haloa/galaksija, pa nije trenutačni fizički bispektar.

## Reproducirani izvorni rezultati

[CI `36331206915` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36331206915) učitao je originalne E8, E9 i E10 Actions artefakte. [E11 kod](../scripts/audit_eboss_dr16_a03_e11_gaussian_squeezed_wind_response.py) računao je četiri unaprijed zaključane kombinacije `256/512` Gauss–Legendre i `48/96` Gauss–Hermite, oba Stein derivative koraka i originalni E10 `r=1` closure. Arhivirani originalni rezultat SHA256 `cade0959c112ca57e9cd5fd4433a49cbec3be0b74eb1f0a1d6eab9350a756ebd` ima 14 823 bajta. [Independently audited E11 sažetak](../source_data/eboss_dr16_a03_e11_frozen_Gaussian_squeezed_source_only_audited_summary_2026-09-27.json) dodatno čuva točne source SHAs, originalnu GitHub CI ID i vlastitu, bez E11 implementacije ponovno izračunanu `96x512` kvadraturu.

| `Γ_ell=E[rT_ell]`, originalni `512x96` | `ell=1` | `ell=3` |
|---|---:|---:|
| FD | 0.0002600278965058659 | 0.0001684107348834053 |
| `F_+` | 0.0002037604735313876 | 0.00012544186994843285 |
| `F_-` | 0.0003160261251109056 | 0.00021117602154044362 |
| `F_+−F_-` | −0.0001122656515795180 | −0.00008573415159201077 |

**Svi intrinzični unconditional odd meanovi iz simetričnoga modela ostaju nula** unutar strojne preciznosti, uključujući diferenciju `F_+−F_-`. To nije međusobno proturječje s nenultim `Γ`, jer `rT(r)` ima paran wind paritet. Neovisna rekonstrukcija iz originalnih E8/E9/E10 bajtova vraća svih šest `Γ` uz max apsolutno odstupanje `2.87314e−18`, E10 `r=1` uz max `3.79471e−19`, neuvjetovani srednji predznak `≤3.38813e−21`; nije pokrenut CLASS niti je otvoren observed odd.

**Numerička QA:** originalni CI max relative 48/96 GH, 256/512 angular i .002/.001 Stein korak `0.00137017`, nula upozorenja prema **prije izračuna zaključanim** tehničkim pragovima. Time je samo numeric source-only reproducibility zadovoljena; fizički eBOSS acceptance ili Einstein–Vlasov numerička truncation/hierarchy convergence nije ovime definiran ni verificiran.

## Što još ne možemo zaključiti i iduća kontrola

`Γ_s` **ne može** zamijeniti originalnu unconditional eBOSS `xi_1` ili njegovu 24D kovarijancu. Nismo računali fizički `P_{D_L,r}(K,z)`, njegov sign/LOS tensor, `K<<k` triangle, full `k_s,z` short response, realan LRG/ELG mass/bias/HOD/magnification/evolution selection, redshift-conditioned finite-bin pair/triple window ni neovisnu bispectrum/mock covariance. Čak i za izvorni model treba provjeriti je li E9 density-derived neutrino–CDM velocity proxy konzistentan s izravnim CLASS `vTk` transferom za **iste zamrznute F±**. Bez te provjere `Γ` je uvjetni *proxy-wind* odziv, ne empirijski kalibrirani ili točan vektorski neutrino vjetar.

**Sljedeća prospektivna E12 odluka:** na unaprijed fiksiranoj dugovalnoj `K` mreži i izvornom CLASS commitu, source-only provjeriti originalni E9 density-derivative wind nasuprot direktnim `vTk` i metric source terms, zatim izvesti gauge/convention-consistent `P_{D_L,r}(K,z)` za definirani *long-mode cb density* te zasebno zapisati koje dodatne galaktičke/hod pretpostavke trebaju za eBOSS. Ako direct/proxy mismatch nije razriješen, nema numeričkoga fizičkog squeezed bispektra u eBOSS jedinicama. [CLASS službena uputa za dTk i vTk](https://github.com/lesgourg/class_public/blob/master/explanatory.ini) razlikuje density od velocity transfera. Opaženi odd SEALED; originalni SGC independent-reverse kvar i 48k random konvergencija ostaju otvoreni; `main` nije mijenjan.
