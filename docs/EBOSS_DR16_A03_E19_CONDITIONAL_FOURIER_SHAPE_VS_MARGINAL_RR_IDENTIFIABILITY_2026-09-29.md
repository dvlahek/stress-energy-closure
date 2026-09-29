# A-03E19 — uvjetni Fourierov oblik i neidentificiran predznak vjetra

**29. 9. 2026. · Eksplicitno retrospektivni matematički rezultat iz već viđenog originalnog E10, bez novih fizikalnih/opažajnih podataka.** [Izvorni E10](../source_data/eboss_dr16_a03_e10_frozen_conditional_angular_wind_parity_summary_2026-09-27.json) Git blob a9e0f49a58b9fd5e0ad975d4d14ba6521a648c28 i [novi source-only E19 rezultat](../source_data/eboss_dr16_a03_e19_posthoc_conditional_Fourier_shape_and_marginal_RR_nonidentifiability_2026-09-29.json) Git blob fc5044684b84ff3de86c236988617ead79d8e7e6 ostaju jasno razdvojeni od stvarnog survey testa.

## 1. Stvarna matematička razlika između originalnih uvjetnih F stanja

Pri jedinom zamrznutom E10 čvoru k=0.05 h/Mpc i matematičkom z=0.95, uz originalni +1σ koherentan LOS vjetar, Fourierovi koeficijenti imaju isključivo relativnu normalizaciju po (b_LRG−b_ELG)P_cb.

| F stanje | Uvjetni T1 | Uvjetni T3 | T3/T1 |
|---|---:|---:|---:|
| FD | 0.00026888572322902307 | 0.00017614445289763344 | 0.6550903885201916 |
| F+ | 0.00022278876746842203 | 0.00014150296987877026 | 0.6351440940523486 |
| F− | 0.00031476164003239810 | 0.00021061456076441418 | 0.6691239782037476 |

Originalni uvjetni dvokomponentni vektori F+/F− nisu sasvim proporcionalni: determinant je 2.382851535110157×10⁻⁹ u kvadratu izvornih koeficijentnih jedinica, a kut među normaliziranim pozitivno usmjerenim oblicima 0.023841179533334908 rad. Razlika omjera T3/T1 iznosi 0.033979884151398965 (5.35 % F+ omjera). Originalni E10 256/512 kutni numerički QA gap 1.6155759955262566×10⁻⁶ nije fizikalna nesigurnost. Ovo je oblikovni rezultat na jednom (k,z), ne mjerenje galaktičke razlučivosti nakon k/z usrednjavanja, Hankela, tracer physics, relativističkih nuisancea i eBOSS operatora. Nema procjene S/N.

## 2. Konstruktivna neidentifikabilnost iz istoga marginalnog RR

Razmotrimo samo ograničeni model dvaju jednako vjerojatnih predznaka lokalnog vjetra, T(+v)=t i T(−v)=−t za fiksan F, jednake |v| i iste ostale uvjete. Neka pozitivni R(x) predstavlja samo marginalnu težinu para pri koordinati x, bez podatka o latentnom vjetru. Za svaki η(x)∈[−1,1] definiramo nenegativne, nepoznate uvjetne težine

\[
W_\pm(x)=R(x)[1\pm\eta(x)].
\]

Za sve η marginalni parni prozor jednak je [W_+(x)+W_−(x)]/2=R(x). Ali ponderirani uvjetni EV mean je

\[
\bar T_{\rm sel}(x)=
\frac{W_+(x)t(x)+W_-(x)[-t(x)]}{W_+(x)+W_-(x)}
=\eta(x)t(x).
\]

Dakle isti objavljeni marginalni RR u ovoj teorijskoj konstrukciji dopušta η=0 i η≠0, uz različit intrinzični srednji izvor. Objavljeni eBOSS LRG/ELG randomi ne nose nezavisni neutrino–CDM wind-sign tag. Ne smije se iz njih izmišljati η niti iz same granice |η|≤1 izvoditi granica ukupnog galaktičkog EV signala. Ako oba smjera koriste isti linearni operator, simetrični intrinzični mean i nakon njega ostaje nula. Standardni relativistički i finite-bin even→odd doprinosi nisu time isključeni niti se smiju proglasiti EV wakeom.

## 3. Što je dokazano, a što je idući fizički zadatak

E19 source-only standard-library kod računa originalne uvjetne omjere/angle/determinant iz SHA-pinned E10, odbija izmijenjen E10, lažni unconditional mean, krivi rezultat i negativne latentne težine; dodatno provjerava jednaki sintetički marginal R za više različitih η. Sintetički R=7.25 nema veze s opaženim randomima. Kod ne otvara E0/E1 RR NPZ, FITS, ASDF, opažene galaksije ni opaženi odd; nije novi E8/E10 re-run.

Sama originalna 24D neponderirana galaktička korelacija ne postaje nenulta EV predikcija samo zbog uvjetno nekolinearnih F± koeficijenata. Potreban je fizički izveden ili nezavisno validiran vjetrom uvjetovani halo–LRG/ELG odziv koji određuje stvarno selekcijsko usrednjavanje, zatim apsolutni prewindow xi_l(s,z), odgovarajući tracer-specific z/selection pair operator i neovisna A-04 kovarijanca. Ili zasebno registrirati i modelirati velocity-conditioned/3pt opservable. Originalni E4/E8/E10/E16/E18, korisnikov izbor objavljenog empirijskog prozora bez kontakta autora, observed odd SEALED, audit grana i draft PR ostaju nepromijenjeni.
