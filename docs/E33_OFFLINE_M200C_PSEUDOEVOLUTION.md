# E33 — dvije epohe ne znače automatski fizikalni rast M200c

**29. 9. 2026. · Offline izvedeni približni kozmološki odnos i eksplicitno hipotetski statični profili.** Stvarni raniji i kasniji Abacus header redshift iz postojećega E29 malog JSON-a jesu `z_e=1.027997062494471` i `z_l=0.952838237036305`. Za ovu dijagnostiku koristimo točno označenu aproksimaciju ravnoga matter+Lambda pozadinskog modela `Omega_m=0.315192`, `Omega_Lambda=0.684808`, **ne** novu točnu massive-neutrino CLASS evoluciju. Ulazni JSON ostaje immutable; nema novih podataka.

Kritična gustoća razmjerna je `H(z)^2`. U navedenoj aproksimaciji

\[
{\rho_{c,l}\over\rho_{c,e}}=
{\Omega_m(1+z_l)^3+\Omega_\Lambda\over
 \Omega_m(1+z_e)^3+\Omega_\Lambda}=0.915023317906.
\]

Kritična gustoća između snimaka pada oko **8.50%**. Sferni overdensity rub `200 rho_c` zato se pomiče i kada je stvarni fizikalni profil 100% statičan. Za isključivo matematički halo `rho(r) proportional r^(-gamma)` s `0<gamma<3`, pri fiksnoj fizikalnoj raspodjeli vrijedi `mean(rho,<r) proportional r^(-gamma)` pa

\[
R_{200c}\propto\rho_c^{-1/\gamma},\qquad
M_{200c}\propto\rho_c^{-(3-\gamma)/\gamma}.
\]

Za hipotetski `gamma=2` isti statični profil daje **+4.5403% prividnoga M200c**, a za `gamma=2.5` **+1.7920%**, uz nula stvarnog akretiranog materijala po konstrukciji. Ti brojevi nisu mase niti profili naših 20 haloa. Stvarna pseudo-evolucija NFW profila i haloska granica ovise o pravoj raspodjeli, povijesti i kozmologiji.

**Posljedica za E29:** čak i da sutra imamo dva prava `M200c`, njihova razlika bez zajedničke fizičke radijalne/lagrangijske definicije ne bi automatski mjerila akreciju koju E29 source-history kernel treba. Danas imamo samo `N` dodijeljenih L1 čestica, koji uopće nije `M200c`. Dosadašnji medijan relativne promjene `N=+3.09%` zato se ne smije uspoređivati s gore navedenim hipotetskim pseudo-rastom kao procjena istinske akrecije ili E28 Wechslerova `alpha`.

[E33 QA skripta](e33_offline_mass_definition_qa.py) prošla je 12 kontrola uključujući granični slučaj `Omega_m=0`, istu epohu i statični mass-radius identitet. [E33 JSON](e33_offline_mass_definition_result.json) uključuje SHA256 prethodnog malog E29 ulaza. Nema novih Abacus ASDF/FITS redaka, stvarnog M200c mjerenja, novog E28 pojasnog rezultata, observed odd ili A04 inferencije.
