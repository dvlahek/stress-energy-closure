# A-03E8 — apsolutno normaliziran neutrinski izvor i fizički dimenzioniran uvjetni gravitacijski kernel

**27. 9. 2026. — ostvareno:** zamrznute stress–energy-matched distribucije `F_+`, `F_-` rekonstruirane su **u originalnoj CLASS normalizaciji**, a njihov **neutrinski statički gravitacijski odziv po jedinici halo-potencijala** određen je u SI jedinicama, uz izričite fizikalne pretpostavke. [Izvršivi E8 rekonstruktor](../scripts/build_eboss_dr16_a03_e8_frozen_resonant_occupation.py), [E8b SI kernel](../scripts/build_eboss_dr16_a03_e8b_conditional_static_halo_green_kernel.py), [SHA-zaključani E8a/E8b protokoli](../source_data/eboss_dr16_a03_e8_frozen_resonant_occupancy_protocol_2026-09-27.json) i [stvarni sažetak](../source_data/eboss_dr16_a03_e8_frozen_kinetic_resonance_and_static_SI_green_result_2026-09-27.json) služe za audit. Source-only CI [E8 `36307484339` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36307484339) i [E8b `36307751800` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36307751800) čuvaju generirane 4000-q CSV i numeričke JSON artefakte. **Nije** dobivena fizikalno normalizirana eBOSS galaktička `xi_1,xi_3` niti je otvoren opaženi odd.

## Problem koji raniji shape-only model nije riješio

Originalni `code/build_lrg_elg_wake_template.py` oblik `wake(s)` dijeli njegovim najvećim apsolutnim članom. Time čuva relativni profil, ali gubi amplitudu. `code/wake_phase7_template.py` dodatno nasljeđuje DESI-BGS literaturno kalibriranu Fisherovu skalu; ona nije izravno eBOSS LRG×ELG galaktičko predviđanje. Nema osnove taj free-amplitude oblik proglasiti apsolutnim `xi_1` na `z\in[0.9,1.0)`.

Naš fizički prvi korak stoga je konstrukcija istoga originalnog sedmerodimenzijskoga null smjera iz `code/class_response_optimize.py`, istih sedam prethodno izabranih `COEFF`, `q=linspace(0,20,4000)`, `m_\nu=0.06` eV, `z_{\rm match}=1100`, `30 %` kapa. Nema novog optimizatora, novoga random seeda, novoga uzorka ili post-hoc smjera.

\[
F_0(q)=\frac{2}{(2\pi)^3}\frac1{e^q+1},
\qquad F_\pm(q)=F_0(q)\pm0.30\,\Delta F_{\mathrm{fixed}}(q),
\quad \max_q|\Delta F_{\mathrm{fixed}}/F_0|=1.
\]

Obje distribucije ostaju pozitivne, a `n`, `rho`, `P` kod `z_match` se podudaraju do maksimalnoga relativnog odstupanja `3.579686991212226e-16`. Izvršivi izlaz 4000-q CSV ima SHA256 `bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0`. To je **apsolutna source-convention normalizacija kinetičkog stanja**, ne kalibracija galaktičkog bias modela.

## Fizički neutrinski response bez proizvoljne maksimalne normalizacije

[Okoli et al., MNRAS 468 (2017), §3.2, Eq. (14)–(20)](https://academic.oup.com/mnras/article/468/2/2164/3063204) izvode imaginarni dio linearnoga statičkoga halo-odziva iz rezonantnog boundary-terma. U njihovoj Fermi–Dirac verziji taj član sadrži `[exp(q_*)+1]^{-1}`. Za **istu** halo-geometriju, isti gravitacijski potencijal i `v_\parallel`, uz izotropne pozitivne `F_\pm` s istom degeneracijskom konvencijom, FD occupancy zamjenjuje se `F_\pm(q_*)/F_0(q_*)` puta FD source. Ova generalizacija koristi linearnost boundary-terma u izvornoj distribuciji. Ne tvrdi da javni rad izračunava naš `F_\pm` ili da se halo-velocitet u dvama fizičkim kozmologijama nužno ne mijenja.

\[
q_*=\frac{m_\nu |v_\parallel|}{c\,T_{\nu0}(1+z)},\qquad
{\cal R}_\pm(q_*)=\frac{F_\pm(q_*)}{F_0(q_*)}.
\]

Pod statičkim **fizičkim**, ne komovirajućim, valnim brojem `k_phys` i `N_\nu=1` za izvorni jedan `ncdm` izračun, koeficijent **gravitacijskoga polja po jedinici Fourierova potencijala** duž `k` jest

\[
\left.\frac{g_{\nu,k}}{\Phi_k}\right|_{\rm FD} =
\frac{2N_\nu G m_\nu^{4} v_\parallel}
     {\hbar^{3}k_{\rm phys}}\;
\frac{1}{e^{q_*}+1},
\qquad
\left.\frac{g_{\nu,k}}{\Phi_k}\right|_\pm =
{\cal R}_\pm(q_*)\left.\frac{g_{\nu,k}}{\Phi_k}\right|_{\rm FD}.
\]

Ovaj scalar duž `k` uključuje predznak `v_\parallel`, a jedinica mu je `m^{-1}`. Fizički valni broj budućega eBOSS računanja mora se povezati s komovirajućim kao `k_phys=(1+z)k_com`, uz potpunu provjeru konvencija potencijala i redshift-time evolucije. Ova statička formula po svojoj definiciji nije eBOSS survey response matrix.

**Uvjetni ilustrativni izračun, ne izmjerena eBOSS brzina niti `z_eff`:** `z=0.95`, `|v_parallel|=200 km/s`, `k_phys=0.05 h/Mpc`, `h=0.6736` daju `q_*=0.12204700265` i

| Statički per-Phi gravitacijski odziv (m⁻¹) | Vrijednost |
|---|---:|
| Izvorni FD | 1.2816113834564112e-27 |
| Matched `F_+` | 1.0442538734468939e-27 |
| Matched `F_-` | 1.5189688934659290e-27 |
| `F_+−F_-` | −4.747150200190350e-28 |

Točno vrijedi `(K_+ + K_-)/2=K_FD` i `(K_+−K_-)/K_FD=-0.3704048092478419` za ovaj uvjetni primjer. Nema posebne procjene signifikantnosti u tim brojevima. [E8b actual source-only CI JSON](https://github.com/dvlahek/stress-energy-closure/actions/runs/36307751800) ima SHA256 `6b8f2b7688a746331ba2ae2776f5ac70ebdbc2b1a4e3dcc6472593a7ef745020`, uz negativne kontrole za zamjenu `F_\pm`, predznak `v_\parallel`, `v=0` i zabranjeni `k_phys≤0`.

## Zašto to još nije apsolutni eBOSS dipol

Statički `g_{\nu,k}/\Phi_k` ne određuje stvarni random-field `v_{\nu c}(k,z)`, njegovu vremensku evoluciju, kozmološki `P_{cb}`, halo/cross-tracer odgovor ili opažanu LRG→ELG projekciju. Publikacija eBOSS LRG×ELG multitracer analize [Wang et al., MNRAS 498 (2020) 3470](https://academic.oup.com/mnras/article/498/3/3470/5897383) prikazuje kombinirani rezultat pri `z_eff=0.77`, što nije kalibracija za našu odvojenu high-z `[0.9,1.0)` selekciju. Nema fizikalne osnove posuditi DESI-BGS `b`, `f_evo`, `s`, HOD ili literaturni 2017 Fisherov normalizacijski faktor te ih proglasiti našim released eBOSS uzorkom.

Za eBOSS `xi_{1,3}` trebamo zasebno opravdan **evoluirani** ne-termalni `F_\pm` CLASS/Einstein–Vlasov relativni transfer/halo potential i wake-to-galaxy map, eBOSS LRG/ELG bias/evolution/magnification i selection response **po kapi i paru**, te [fizički z-uvjetovan empirical RR prozor](EBOSS_DR16_A03_E6_THEORY_TO_EMPIRICAL_WINDOW_BRIDGE_2026-09-27.md). Tek uz zasebno validiranu [A-04 eBOSS 24D ili unaprijed novu pD kovarijancu](EBOSS_DR16_A04_EBOSS_OBSERVABLE_COVARIANCE_READINESS_2026-09-27.md) moguće je predregistrirati dopuštenu numeričku pogrešku i fizički `A_{\rm EV}/\sigma_A`. Devet starih mock ID-jeva i DESI 120-mock kovarijanca nisu taj `C_{\rm eBOSS}`.

**Izvršni status:** prvi dvije source-only fizikalne etaže (kinetički izvor, uvjetni dimenzionalni halo Green coefficient) završene. Novi CLASS/FITS izračun još nije pokrenut i nije impliciran iz zelnog E8 CI-a. Originalni opaženi odd **SEALED**, survey-window A-03 **PHYSICAL_UNCERTIFIED**, inferencijski A-04 **BLOCKED**; izvorni v1 SGC neovisni reverse kvar nije riješen ovim teorijskim računanjem. Ne otvarati opažene odd podatke, ne preuzimati nove mockove niti miješati ovaj `k_{\rm phys}` s komovirajućim eBOSS `k`.
