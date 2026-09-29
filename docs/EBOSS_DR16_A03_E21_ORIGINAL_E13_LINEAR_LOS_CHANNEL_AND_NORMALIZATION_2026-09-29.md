# A-03E21 — izvorni CLASS linearni LOS kanal postoji; galaktički odziv nije identificiran

**29. 9. 2026. · Posthoc znanstvena usporedba već izvršenoga E13; nova E21 aritmetika/QA preregistrirana je nakon pregleda E13 izvora i prije E21 CI.** [Originalni E21 protokol](../source_data/eboss_dr16_a03_e21_linear_LOS_velocity_two_tracer_response_normalization_protocol_2026-09-29.json), Git blob 7fc84392a915bf3b242bbc12a27e78267dd643f1. [E21 mali izvorno izvedeni rezultat](../source_data/eboss_dr16_a03_e21_archived_F0_Fpm_linear_LOS_velocity_cross_unit_response_2026-09-29.json), Git blob cc3312c5313cce60fad4a88c3a56ef3d4e675e38. Ni originalni E8/E10/E12/E13/E19/E20 ni njihove kozmologije ili podaci nisu promijenjeni. Nema novoga CLASS, FITS, ASDF, WSL, opaženih galaksija ili opaženoga eBOSS odd.

## 1. E20 je no-go samo za svoj kvadratni model

E20 je pokazao da za odabrani kvadratni izvor S∼vδ prvi cross-power član ovisi o miješanom B_vδδ, koji je nula za zajednički centrirana Gaussova linearna polja. Međutim, drukčiji, **linearni LOS-projicirani** fenomenološki odgovor ne ovisi o tom trećem momentu. To ne znači da je svaki linearni v_nu−v_cdm galaktički član fizički dopušten ili nenulti: stvarni redshift-space odabir, halo ubrzanje/velocity bias, gauge i tracer evolucija moraju ga izvesti i odrediti.

Definiramo samo kao transparentan model

\[
\delta_L^s=A_L(k,\mu)\delta_{cb}+c_L\,v_{r,\mathrm{LOS},R16},\qquad
\delta_E^s=A_E(k,\mu)\delta_{cb}+c_E\,v_{r,\mathrm{LOS},R16}.
\]

\(A_L,A_E\) su stvarni, parni u μ i zasad neodređeni standardni tracer odzivi. \(c_L,c_E\) su stvarni i nepoznati koeficijenti \((\mathrm{km/s})^{-1}\); **E13 ne mjeri, ne izračunava i ne dokazuje da su nenulti**. \(v_{r,\mathrm{LOS},R16}\) jest zamrznuta izvorna R16-filtrirana linearna CLASS neutrino–CDM relativna LOS brzina, ne opažena brzina pojedinačnog LRG/ELG haloa.

Izvorni [E12 direct-vTk kod](../scripts/audit_eboss_dr16_a03_e12_direct_class_vTk_long_cross.py) i [E13 izvješća](../source_data/eboss_dr16_a03_e13_archived_CI_2026_09_27/e13_rank_matched_direct_vTk_gaussian_fixed_mode_source_only.json) u točnoj konvenciji E12 daju

\[
P_{\delta_{cb},r_F}(K,\mu)=i\mu\,C_{r,F}(K),\quad
r_F=v_{r,\mathrm{LOS},R16,F}/\sigma_{R16,F},\quad
P_{\delta_{cb},v_F}=i\mu\,C_{v,F},\quad
C_{v,F}=\sigma_{R16,F}C_{r,F}.
\]

Za \(\langle\delta_L^s\delta_E^{s*}\rangle\), uz \(P_{v\delta}=P_{\delta v}^*\), slijedi

\[
\operatorname{Im}P_{LE}
=\mu(A_Lc_E-c_LA_E)\,C_{v,F}(K).
\]

To je formula *po jedinici nepoznatoga fenomenološkog odziva*, ne fizički izračunat neutrinski wake u eBOSS-u. Zamjena L↔E i μ→−μ mijenjaju predznak; jednaki A i jednaki c daju nulu. Samo jednaki A nisu dovoljni za nulu ako se c razlikuju. Postojeći u literaturi standardni relativistički Doppler također može stvarati imaginarni različito-tracerski cross-power; njegov *galaktički* LOS velocity field ne smije se poistovjetiti s originalnim neutrino–CDM relativnim v poljem. Usporedi [relativistički odd cross-power](https://academic.oup.com/mnras/article/501/2/2547/6041033) i [baryon-CDM velocity-bias analogiju](https://doi.org/10.1103/PhysRevD.94.063508), koje nisu naša kalibracija neutrinskih koeficijenata.

## 2. Stvarna aritmetika iz četiri originalna E13 K čvora

Izvorni direct-vTk E13 SHA-pinned FD/F+/F− source JSON-ovi za matematički z=0.95 (NIJE release eBOSS z_eff), K=(.001,.002,.003,.005) h/Mpc i isti top-hat R16 daju \(\sigma_{R16,\rm FD}=135.6779814380\), \(\sigma_{R16,+}=135.8620668110\), \(\sigma_{R16,-}=135.4853716298\) km/s. Ti σ su state-specific iz originalnoga direct CLASS inputa. Zamrznuti F± shape na k=.05 iz E10 NIJE long-K C ovog E21.

| K (h/Mpc) | F+−F− za C_r (Mpc³) | F+−F− za C_v ((km/s) Mpc³) |
|---:|---:|---:|
| 0.001 | +138.1313173730 | +30715.8842886183 |
| 0.002 | +173.4247690434 | +55675.7625894658 |
| 0.003 | −130.6403258667 | +37370.1576352790 |
| 0.005 | −233.6816249870 | +60749.7454510555 |

Dakle, na dva originalna K čvora sama usporedba po različitim vlastitim state-normalizacijama mijenja **predznak međustanjskog kontrasta**. To nije potvrda fizičkoga predznaka galaktičkoga wakea: držati \(\sigma_F\) i interpretaciju koeficijenata dosljednom. Ako je fizički koeficijent \(c_{a,v}\) isti pri čistoj promjeni varijable, u originalnoj normiranoj reprezentaciji odgovara mu \(c_{a,r}=c_{a,v}\sigma_F\). Fiksirati c_r kroz sva F stanja značilo bi drugi fizički model od fiksiranja c_v.

[E21 standard-library provjera](../scripts/audit_eboss_dr16_a03_e21_original_E13_linear_LOS_velocity_units.py) računala je originalnih 3 stanja × 4 K × 2 jedinice, unakrsno provjerila originalni E13 aggregate, svih pet negativnih kontrola i tracer/μ parity. [Originalni E21 CI 36564575568 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36564575568), bez novih CLASS/observed podataka i bez ponavljanja izvornog E13 testa. Nije provedena statistička hipoteza, novi science cut ili survey forecast.

## 3. Zašto to nije u proturječju s E10/E19

E10 u strogoj fiksoj **koherentnoj** dvogranoj \(\pm v\) konstrukciji ima \(\langle T(+v)\rangle_{\pm}=0\) uz od predznaka nezavisnu selekciju. E21 računa drukčiji objekt: prostorni Fourierov **dvotočkasti** \(P_{\delta v}=\langle\delta(k)v(k)^*\rangle\), koji može biti nenulti čak i za \(\langle v\rangle=\langle\delta\rangle=0\). Zato samo E10 jednopointno predznačno usrednjavanje ne dokazuje da je svaki linearni relativno-brzinski galaxy cross-power nula. Ali E21 također ne uklanja E19 problem: empirijski marginalni RR ne određuje fizikalne, tracer- i velocity-sign-uvjetovane \(c_L,c_E\) niti puni halo selection.

## 4. Fizički STOP

Dobiveni C_v(K) čvorovi su stvarni **linearni izvorni CLASS velocity–density cross** izražen u dosljednim jedinicama. Nisu fizički halo–galaxy coupling, predwindow \(\xi_\ell(s,z)\), realni galaxy Doppler, neutrinski wake 24D oblik ili njegov unblindani S/N. K čvorovi su dugi, a E10 k=.05 je zaseban uvjetni izvor; interpolation/extrapolation 4 K na punu eBOSS 24D mrežu nije autorizirana. Fizički isti-likelihood model mora zasebno odrediti c_a ili pokazati zašto bi bili nula, standardne relativističke/wide-angle/evolution nuisance, puni k/z, stvarni released random prozor i valjanu neovisnu kovarijancu. Alternativni E20 kvadratni kanal ostaje zaseban i zahtijeva B_vδδ plus LOS-odd projiciranje.

**Odluka:** dvije teorijske rute su otvorene, a ne spojene ad hoc: (a) fizički specificiran linearni neutrino–CDM LOS tracer response, za koji E21 već daje originalnu linearnu cross-power podlogu; (b) stvarni nelinearni mixed-bispectrum/halo kernel iz E20. Nijedna ruta trenutno nema validirani LRG/ELG c_a ili Γ; ne otvarati observed 24D odd, ne nabavljati nove mockove niti spajati PR.
