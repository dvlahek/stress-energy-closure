# E27 — prvi uvjetni fizičko-vremenski neutrinski wake, ne galaktički signal

**29. 9. 2026.** Ovo je nova numerička integracija izvornog F0/F+/F− Vlasovljeva kernela kroz fizičko konformalno vrijeme na unaprijed zadanoj uvjetnoj halo-povijesti. Nije puni samokonzistentni Einstein–Vlasov, kalibrirani LRG/ELG halo drag, ukupna sila, opaženi eBOSS odd ili detekcijska značajnost.

## Točno računsko pitanje

Originalni E17D0 dao je samo impulsnu težinu po jedinici vanjskoga potencijala. Originalni E17D1 imao je simboličko bezdimenzijsko vrijeme. E17D2b1 je sačuvao dvije nekalibrirane Wechsler-like M200c putanje α=.4 i .8 za dvije deklarirane fizičke masene reference, ali nije izračunao fizički retardirani wake. E27 kombinira iste zamrznute podatke, **ne mijenjajući originalne F raspodjele**.

Konačan, potpuno specificiran Bornov benchmark koristi izvorni CLASS H_F(z)/c u Mpc⁻¹, samo stare originalne z čvorove za linearnu interpolaciju na intervalu z=.95..1, izvorni E8 konačni q∈[0,20], mν=.06 eV i Tν0=.71611×2.7255 K×8.617333262e−5 eV/K. Kratkovalni k=.05 h_CLASS/Mpc, h_CLASS=.6736, vanjski point-halo profil U(k)=1, fiksirana relativna brzina haloa +200 km/s uz µ=+1, početni wake na z=1 postavljen na nulu. Zbog vremenskoga početnog uvjeta rezultat uključuje samo izvor unutar z=1→.95, ne prethodnu halo memoriju.

U originalnoj E17D0 konvenciji impulsni kernel je

\[
\mathcal K_F(s)=-
\frac{\int_0^{20} dq\,q^2 F'_F(q)\operatorname{sinc}'(sq)}
{\int_0^{20}dq\,q^2 F_F(q)}.
\]

Fizička faza slobodnog leta je s=k(Tν0/mν)∫(1+z') dz'/H_CLASS,F(z'), a pomak putujućega haloa Δx=(v_h/c)∫dz'/H_CLASS,F(z'). U komovirajućoj Fourierovoj konvenciji volumena vanjski Poissonov potencijal točkastoga haloa je Ψ_h(k,z)=−4πG Mα(z)/(a c² k²), s izvornom uvjetnom Mα(z)=M0 exp[−α(z−.95)/1.95]. Zato je računani Fourierov izvor

\[
\frac{\delta n_{\nu,F}^{w}(\boldsymbol k,.95)}{\bar n_{\nu,F}}
=-\frac{4\pi GM_0}{c^2 k}\frac{m_\nu}{T_{\nu0}}
\int_{.95}^{1}\frac{dz}{H_{{\rm CLASS},F}(z)}
\frac{M_\alpha(z)}{M_0}\,
\mathcal K_F[s_F(k;z,.95)]\,
\exp[i k\mu\Delta x_F(z)] .
\]

Ishod je u Mpc³ po odabranom Fourierovu modu, a **ne** prostorna gustoća u pojedinom halo radijusu ni potpuno gravitacijsko ubrzanje. Fizikalna jednadžba je ograničena na vanjski Newtonov linearni Born, zanemaruje samogravitaciju wakea, povratnu reakciju haloa, okoliš i metric perturbation closure. Originalne E17D2b1 M200c povijesti ostaju nekalibrirane, jer Wechslerov objavljeni oblik potječe iz druge masene definicije.

## Novi numerički rezultati

[Prvi E27 CI 36599144647](https://github.com/dvlahek/stress-energy-closure/actions/runs/36599144647) SUCCESS nakon dvaju usko tehničkih, fail-closed pokušaja prije samoga fizikalnog izračuna. Prvi je bio skraćeni Git blob originalnoga E17D2b1 slučaja 03 u novome čitaču, drugi NumPy 1.26 API naziv trapezoid umjesto trapz. Nisu mijenjani izvori, parametri ili jednadžba.

Za izvorno niže maseno sidro M0=1.4903129657228018e12 fizičkih Msun, pri z=.95, µ=+1 i v_h=200 km/s, imaginarni konačno-vremenski odgovor [Mpc³] jest:

| F stanje | α=.4 | α=.8 |
|---|---:|---:|
| FD | 2.62774557197819e−5 | 2.60767355606882e−5 |
| F+ | 2.62773211932041e−5 | 2.60766021764594e−5 |
| F− | 2.62775903238812e−5 | 2.60768690218459e−5 |

Za FD razlika dviju **uvjetnih α** povijesti iznosi 0.763849% referentnoga α=.4. Druga izvorna fizička masena referenca daje deset puta veći rezultat: to je točno linearno M0 skaliranje vanjskoga point-source Bornova modela, ne univerzalni zakon stvarnoga neutrinskog haloa. Originalni F± izvori nisu retunirani.

Prva pojedinačna vremenska 129-vs-257 kvadratura imala je relativni imaginarni jaz približno 2.2e−5, veći od maloga F± kontrasta. Zato prvi E27 run sam NIJE opravdao objavu preciznoga međustanjskog kontrasta. [Nakon prvih rezultata eksplicitno registriran E27R1 preciznosni protokol](../source_data/eboss_dr16_a03_e27r1_posthoc_Fpm_contrast_resolution_protocol_2026-09-29.json), Git blob a1bdd5264ebe914f0896cc7b60ec8919bb37434e, prije novih 513/1025 računanja zamrznuo je istu izvornu masu, oba α i zahtjev relativnog jaza manje od 1% **same potpisane razlike**.

[Precizni E27R1 CI 36599391672](https://github.com/dvlahek/stress-energy-closure/actions/runs/36599391672) SUCCESS:

| α | Im(F−−F+) [Mpc³], 1025 čvorova | (Im_F−−Im_F+)/Im_FD | 513–1025 relativni jaz kontrasta |
|---|---:|---:|---:|
| .4 | 2.69124460564041e−10 | 1.02417184096e−5 | 4.6196774e−6 |
| .8 | 2.66839249508282e−10 | 1.02329166299e−5 | 4.5996901e−6 |

Kontrast je time **numerički razlučen unutar deklariranoga originalnog 4000q i konačno-vremenskoga Bornova modela**. Međustanjskih približno 0.0010% nije fizikalni limit ukupnoga neutrinskog signala, ni greška opservablea ni detekcijski S/N. Statički E8b rezonantni kernel računa drugi fizički opservable i drugo k značenje, stoga njegov raniji F kontrast ne treba izravno numerički uspoređivati s ovim.

## Što nas doista zaustavlja

Za silu i ubrzanje na halo treba cijeli 3D k integral s fizički konzistentnim konačnim halo profilom, okruženjem i dužom poviješću izvora; jedan k nije ukupni drag. Za high-z LRG×ELG χ treba stvarna, same-sample halo history/HOD, observer Doppler/selection i prethodno definirani signed pair window. Originalni Abacus z=.953 L1 N ne daje istodobno kalibrirani M200c(a) ni punu particle history; F± halos nisu neovisno simulirani. Opaženi 24D odd ostaje SEALED, A03 fizička galaktička prognoza i A04 inferencija BLOCKED. Nema novog CLASS, ASDF, FITS, mockova, WSL, rezova/seedova, autora ni mergea; main netaknut, PR draft.
