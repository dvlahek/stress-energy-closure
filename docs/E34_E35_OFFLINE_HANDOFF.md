# EinsteinVlasovNP — E34/E35 offline nastavak, 29. 9. 2026.

## Nova dva rezultata (matematička, ne mjerenje)

**E34:** Pozitivna sferna ljuska stalne mase s dvije različite, kontinuitetno konzistentne međuvremenske putanje ima iste radijuse, mase, profile, centre i nulte radijalne brzine na oba krajnja snimka, ali daje različite retadirane vanjski Born-profile integrale. Pri isključivo ilustrativnom bezdimenzijskom `K=exp(t-1)`, `r_±=1±0.1sin²(pi t)`, `k=1`, nalaz `J_+=0.439549571374331`, `J_-=0.455163833342087`, razlika `0.0156142619677557`. Egzaktan predznak slijedi iz monotonosti sinc i pozitivnog toy kernela. 51 kontrola PASS. To NIJE stvarni E28 `K_F`, samogravitirajući halo, fizička F± razlika ni numerički total drag. Ako fizička teorija *neovisno* ograniči `|du/dt|≤L`, derivirana je konzervativna vremenska interpolacijska ograda.

**E35:** Za pozitivne sferne masene profile sa poznatim drugim momentima postoji stroga `|1-u_src u_test|≤k²(<r²>src+<r²>test)/6`, uz bezuvjetni plafon 2. U fiksnom vanjskom linearnom retardiranom integralu, apsolutna error ograda zahtijeva **momente source profila za svako doprinosno vrijeme**, a ne samo u dvije krajnje epohe. Ilustrativne pozitivne smjese daju pri bezdimenzijskom k=1 stvarni faktor-error `0.153739548347` i momentnu ogradu `0.1655`; endpoint-only negativna kontrola daje stvarni error `0.337494275778` nasuprot pogrešnoj granici `0.0075`. 85 kontrola PASS. Ni jedna od ovih brojki nije E28 fizički k, R, sila ili fizički kmax.

## Znanstvena integracija s E29 i prioritet

E29 20 od 20 lokalno PID-uparenih haloa i njihova prostorna/metapodatkovna QA ostaju valjani u svom lokalnom opsegu, no ne sadrže same-object, dinamički konzistentan `M200c(t)`, kontinuirani `u_src(k,t)` niti početni neutrinski wake. E32 dokazuje da i svi L2 percentili ostavljaju neodređen profil između percentila; E33 da promjena `M200c` može biti pseudo-evolucija. E34 dodatno dokazuje da i hipotetski dva **potpuna** krajnja profila ne identificiraju retardiranu povijest. E35 pokazuje koji konkretan dodatni input (uniformno opravdan vremenski radijalni moment, ili još jači fizički profil/transport) može omogućiti kontrolirani **low-k** test bez izmišljanja NFW fitanja. E28R1 ostaje uvjetni finite-band rezultat s 94.22589% FD low-mass alpha.4 iz k>1/Mpc. Ne proglašavati UV certifikat niti inferirati fizičku F± galaktičku razliku.

**Prvi sljedeći fizički input nije novo otvaranje Abacus ASDF percentile stupaca.** Potrebna je vanjski opravdana dinamička halo+environment povijest, početni neutrinski odziv ili njegova ograda i fizikalni kratkovalni Born-validity test. Zasebno E30: eBOSS 24D linearni LOS kanal treba neovisno izračunan high-z `chi_F=A_L c_E-A_E c_L` i točan number-count/window transfer. Postojećih devet eBOSS mockova ne daje invertibilnu 24D inferencijsku kovarijancu.

## Pravila i reprodukcija

GitHub `dvlahek/stress-energy-closure` ostaje na radnoj grani `audit/eboss-elg-bit8-ra-orientation-20260925`, draft PR #1; **ništa nije poslano na GitHub i main nije mijenjan**. Originalna F0/F+/F−, E8 4000q, E16, E27/E28, E4, E7 48k i LOS/cuts/seedovi zamrznuti. Observed eBOSS galaxy rows i odd vektor **SEALED**. A-03 PHYSICAL_UNCERTIFIED, A-04 BLOCKED.

Za lokalnu reprodukciju samo standardni Python 3: `python e34_offline_causal_history_qa.py` i `python e35_offline_lowk_profile_bound_qa.py` u direktoriju otpakiranog E34/E35 paketa. Skripte odbijaju prepisati postojeći rezultat ako se bajtovi razlikuju. U paketu su navedeni prethodni mali parent dokumenti i E29 JSON za potpuni parent-SHA audit; nisu uključene velike ASDF/FITS datoteke niti opaženi odd.
