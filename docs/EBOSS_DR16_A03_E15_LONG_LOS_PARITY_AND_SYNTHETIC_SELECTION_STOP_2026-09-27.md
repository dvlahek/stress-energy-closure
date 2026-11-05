# A-03E15 — long-LOS paritet izvornog E14 source-only odziva i granica survey-averaginga

**27. 9. 2026. | Status: source-only E15 i neovisna matematička provjera PASS. Stvarni eBOSS bispektar, trotočkasti prozor i kovarijanca i dalje nisu izračunani.**

E12 je na izvorna tri zamrznuta `F_0,F_+,F_-` provjerio CLASS direct `vTk` prema E9 density-derived proxyju. E13 je odvojeno kalibrirao originalni `R=16 h^-1 Mpc` direct-vTk Gaussian source rank i dugo `P_{\delta_{cb},r}(K)`. E14 je iz originalnoga trostanjskog CLASS `P_{cb}(k)` izveo **72 fiksna source-only produkta**, `3` stanja × `3` kratka `k` × `4` duga `K` × `2` kratka odd multipola; neovisna matematička reprodukcija provjerila je svih `72` i `24` izvornu razliku `F_+−F_-`. Sva izvorna E8–E14 numerika i SHA-manifesti trajno su na audit grani, ne samo u privremenim GitHub Actions artefaktima. Opaženi odd ostaje SEALED.

## Problem koji E14 nije riješio

Arhivirani [E14 originalni izvor](../source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json) ne pohranjuje konačni opaživi bispektar, nego **realni reducirani koeficijent**

\[
S_{\ell_s}^{(s)}(k_s,K_L)
=\frac{B_{{\rm source},\ell_s}^{(s)}}
{i\mu_K\,(b_{\rm LRG}-b_{\rm ELG})},
\qquad \mu_K=\hat{\boldsymbol K}\cdot\hat{\boldsymbol n},
\]

u jedinicama `Mpc^6`. Indeks `ell_s=1,3` već je ranije projektiran **kratkom** LOS kutnom integracijom. E14 nije modelirao zatvaranje punog trokuta `k_3=-k_s-K_L`, azimut ili zajedničku selekciju dugog i dvaju kratkih modova. Ne smije se iz `S` zaključiti da je izračunan skalarni, kutno usrednjeni bispektar ili stvarna eBOSS `xi_odd`.

[E15 protokol](../source_data/eboss_dr16_a03_e15_long_los_parity_and_selection_leakage_prereg_2026-09-27.json) zaključan je **nakon** E14 outputa i prije E15 kutnoga računanja. Točno čuva originalni E14 SHA256 `0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8`, predeklarirane `72+24` komponente, originalne `k_s=(.05,.075,.1)h/Mpc`, `K_L=(.001,.002,.003,.005)h/Mpc`, fiksni `z=.95` samo kao matematički čvor i dvije kutne kvadrature `32/64`.

## Provjerena nužna long-LOS geometrija

E15 vraća **samo faktor koji je E14 eksplicitno izdvojio**:

\[
\frac{B_{{\rm source},\ell_s}^{(s)}
(k_s,K_L;\mu_K)}{b_{\rm LRG}-b_{\rm ELG}}
=i\mu_K S_{\ell_s}^{(s)}(k_s,K_L).
\]

Dugi LOS multipol `L` ovog ograničenog modela definiran je zasebno od kratkog `ell_s`:

\[
B^{(L)}_{\ell_s}/(i\,\Delta b)
=\frac{2L+1}{2}
\int_{-1}^{1} d\mu_K\,P_L(\mu_K)\mu_K S_{\ell_s}.
\]

Po Legendreovoj ortogonalnosti, **`L=1` daje `S`**, a `L=0,2,3,4` su **nula** u ovom posebno pretpostavljenom separabilnom izvornom modelu. Pri `K_L→-K_L` imaginarni član mijenja predznak, kako zahtijeva kompleksna konjugacija Fourierova bispektra realnih polja. Zamjena LRG/ELG redoslijeda također mijenja predznak antisimetriziranog imaginarnog člana. To nije dokaz da stvarna nelinearna bispektralna kutna funkcija nema druge `L` multipole, jer E14 ne sadrži punu geometriju trokuta.

### Strogo sintetički selection negative control

Radi prepoznavanja mogućega pogrešnog usrednjavanja, a **ne** kalibracije eBOSS prozora, unaprijed je odabrano `epsilon=0.1`. Za ilustrativne pozitivne kutne težine `W_even=1+epsilon P2(mu_K)` i `W_odd=1+epsilon P1(mu_K)`, normirana kutna sredina `<B>_W` daje:

\[
\langle B/(i\Delta b)\rangle_{W_{\rm even}}=0,\qquad
\langle B/(i\Delta b)\rangle_{W_{\rm odd}}=
\frac{\epsilon}{3} S_{\ell_s}.
\]

Za te iste težine `L=1` koeficijent jest `S(1+2\epsilon/5)` za `W_even`, odnosno `S` za `W_odd`. **`W_odd` je samo toy geometrija.** Niti njegova `epsilon`, niti rezultat `S/30` nisu procjena stvarnoga eBOSS LRG/ELG chunk/depth/redshift leakagea ili dozvoljena acceptance tolerancija.

## Dovršeni brojčani audit na nepromijenjenom E14

[Izvorna E15 CI `36338418964` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36338418964) izvodi [E15 originalni kod](../scripts/audit_eboss_dr16_a03_e15_long_los_parity_and_synthetic_selection.py) nad [E14 byte-exact source](../source_data/eboss_dr16_a03_e14_archived_CI_2026_09_27/e14_frozen_multik_direct_rank_unitbias_source_only.json). [Puni izvorni E15 JSON](../source_data/eboss_dr16_a03_e15_original_E14_long_LOS_parity_and_toy_selection_2026_09_27.json) od `164466` bajtova, SHA256 `6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983`, **trajno je committan** na audit granu. Uz `32/64` GL max `6.405466435044362e-16` skaliranoga numeričkog gapa, `0` tehničkih upozorenja.

Za originalni **ilustrativni**, neizmjereni `k_s=.05h/Mpc,K_L=.005h/Mpc,ell_s=1`, razlika `F_+−F_-`:
- reducirani izvorni E14 `S=-429268.1132787536 Mpc^6`;
- dugovalni `L=1` imaginarni koeficijent po jediničnom `Δb` `-429268.11327875365 Mpc^6`;
- neponderirani `L=0` `0` unutar numeričke preciznosti;
- **isključivo toy** `W_odd=1+0.1P1` kutna sredina `-14308.937109291786 Mpc^6`.

[Neovisna E15 CI `36338533938` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36338533938) ne učitava originalnu E15 implementaciju niti NumPy. [Čisti Python standard-library audit](../scripts/audit_eboss_dr16_a03_e15_independent_analytic_parity_replay.py) iz originalnoga E14 i arhiviranog E15 provjerava svih `72` slučajeva i `24` kontrasta izravnim analitičkim Legendre identitetima, max skalirani gap `1.1960918366860087e−15`; izvorni audit SHA256 `29973819c235e7c5dbbcaefe5e68b8ebd33677669b4dbf8bd6ec30f1768007a9` također je [trajno spremljen](../source_data/eboss_dr16_a03_e15_independent_original_72_angular_replay_2026_09_27.json). To je nezavisna **algebarska** kontrola istog modela, ne novi neovisni fizikalni halo simulator.

## Nužan idući fizički korak i STOP

Prije nego se `i\mu_K S` pretvori u pravi bispektar, treba fiksirati punu zatvorenu geometriju `(\boldsymbol k,\boldsymbol K,-\boldsymbol k-\boldsymbol K)`, odgovore **obaju** kratkih krakova, azimut i relativne orijentacije u Fourierovoj trotočki. Trenutačni E14/E15 `k_s^{-2}` nastao je iz kvazistatičke, *ne* vremenski retardirane Eq20 aproksimacije, a `K_L/k_s` doseže `0.1`; same izvorne 72 algebraičke vrijednosti nisu fizičko ograničenje izostavljenih konačnih `K/k` članova. Tek potom dolaze stvarni LRG/ELG high-z bias/HOD, standardne relativističke i wide-angle nuisance komponente, stvarni NGC/SGC chunk/depth/z triple selection window i dovoljno neovisnih mockova za novu bispektralnu kovarijancu.

Originalni unconditional eBOSS 24D odd vektor nije bispektar, opaženi odd ostaje SEALED, originalni SGC reverse audit i finite-random 48k konvergencija ostaju otvoreni. Nema novih FITS/mock preuzimanja, originalnih `F_\pm` promjena, dodatnih survey cutova, seed tuninga, pristupa opaženom odd vektoru ili promjene `main`; draft PR ostaje otvoren.
