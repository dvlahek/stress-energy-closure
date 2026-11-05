# A-03E6 — teorijski wake predložak u eBOSS finite-bin RR prozoru

**Status, 27. 9. 2026.:** matematičko/sintetičko sučelje je implementirano i provjereno. **Fizički predviđena eBOSS amplituda, konačni inferencijski pD vektor, fizički budget sistematika i eBOSS kovarijanca nisu izračunani.** [Prije novog E6 testa zaključani protokol](../source_data/eboss_dr16_a03_e6_physical_template_to_empirical_window_bridge_protocol_2026-09-27.json) nastao je *nakon* E4 i eksplicitno ne dopušta naknadni E4 prihvatni prag. Nema novih FITS podataka, dodatnih randoma, pristupa opaženom odd vektoru, kontakta s autorima ili promjene `main`.

## Što stvarno imamo u teorijskom kodu

- [`code/build_lrg_elg_wake_template.py`](../code/build_lrg_elg_wake_template.py) računa predznak i oblik Hankelova dipola iz zamrznutog Einstein–Vlasov hidden-state smjera, ali zatim `wake` i `dop` zasebno dijeli njihovim najvećim apsolutnim vrijednostima. Izlaz je `shape-only`; ukupnu LRG−ELG amplitudu prepušta fitu. To nije apsolutna prognoza `ξ_1` u eBOSS jedinicama.
- [`code/build_lrg_elg_physical_odd_basis.py`](../code/build_lrg_elg_physical_odd_basis.py) računa DESI pred-prozorske `wake`, relativističku `ν_1` i wide-angle bazu. Postupak usrednjavanja `s² ds` pretpostavlja idealne sferne ljuske, a kod izričito zahtijeva vanjsku kalibraciju LRG/ELG bias/evolution/magnification koeficijenata i DESI prozor. To nije eBOSS operativni operator. [DESI z-razlučena 18D baza](../code/build_lrg_elg_zresolved_physical_basis.py) s tri z-bina također nije postojeći eBOSS high-z 24D pilot.
- [Postojeći eBOSS RR operator](../scripts/build_eboss_dr16_conditional_window.py) već daje *uvjetni finite-bin* linearni odgovor `M_{cap,ell_out,ell_in}` sa šest izlaznih s-binova i 120 fine-s stupaca. Originalni [E0/E1 izvorni operator/protokol](../source_data/eboss_dr16_a03_empirical_rr_window_injection_protocol_2026-09-26.json) ima 203 polja, uz odvojene NGC/SGC i devet mock ID-jeva. Njegove originalne 250104 NPZ bajtove, SHA256 `4f63e04ce4c5cfdb9404e5a56dee44907cb9a28cc015c6ed27ac718195feae6b`, **nismo ponovno otvorili** u E6.

## Definirana matematička veza

Za poznati pred-prozorski model na 120 fine-s binova, pri odgovarajućoj fizički opravdanoj redshift aproksimaciji,

\[
\xi_{\mathrm{pre}}^c(s_i,\mu)=\sum_{L=0}^{4}\xi_L^c(s_i)P_L(\mu).
\]

Očekivani coarse-s odd output postojećeg conditional RR operatora jest

\[
y_{o,J}^c=\sum_{L=0}^{4}\sum_{i=0}^{119}
M_{oL}^{c}[J,i]\,\xi_{L}^c(s_i),
\quad o\in\{1,3\},\quad c\in\{\mathrm{NGC},\mathrm{SGC}\}.
\]

`M` se određuje samo postojećim objavljenim tracer-specifičnim `R_\mathrm{LRG}R_\mathrm{ELG}` parovima, uz neovisnu normalizaciju za svaki mu-bin i točne Legendre integrale istog eBOSS signed-mu grida. Sve `ξ` komponente i output moraju biti **bezdimenzijski korelacijski koeficijenti**, uz zamrznutu orijentaciju **LRG→ELG** i midpoint LOS. Pri reverziji tracera i signed-mu definicije odd doprinos mijenja predznak; nije dopuštena prešutna zamjena orijentacije.

Razdvajanje nužno zadržava **svih pet** ulaznih multipola:

\[
y_o^c=\underbrace{\sum_{L=1,3}M_{oL}^c\xi_L^c}_{\text{intrinzični odd}}
+\underbrace{\sum_{L=0,2,4}M_{oL}^c\xi_L^c}_{\text{even-to-odd finite-bin odgovor}}.
\]

Za konstantnu even `ξ_L(s)` i točnu signed-mu ortogonalnost drugi doprinos nestaje uz kontrolirane edge-korekcije. Za `ξ_L` koja varira unutar šireg s-bina i mu-ovisnu RR raspodjelu on općenito ne nestaje. Sintetička matematička provjera **nije** izmjereni eBOSS fizički leakage niti dokaz da je zanemarenje ulaznih `L≥5` multipola fizički opravdano.

### Pilot dimenzija i nedostajući z

Zadržani tehnički redoslijed je NGC `ℓ1` šest s-binova, NGC `ℓ3` šest binova, zatim SGC `ℓ1` i SGC `ℓ3`, ukupno **24D**. To je samo fiksno dijagnostičko preslikavanje, **ne nova konačna inferencijska selekcija**. Konačni eBOSS observable može imati drugu, prije mjerenja fizikalno obrazloženu definiciju, ali iz postojećeg 24D pilota se ne smije nakon E4 proizvoljno uklanjati šest komponenti radi oznake „18D”.

**Važna fizička granica:** eBOSS E0/E1 RR operator odnosi se na `[0.9,1.0)` i marginalizira parni redshift. Realni Einstein–Vlasov `ξ_L(s,z)` i standardni relativistički/selection doprinos ovise o `z`. Prije fizičke prognoze nužno je opravdati effective-z aproksimaciju, **uz zasebnu cap i pair-selection provjeru**, ili u novom, zasebno unaprijed registriranom poslu validirati z-uvjetovani RR operator. `z_\mathrm{eff}=0.95` ne smije se jednostavno pretpostaviti kao izmjereni effective-z.

## Kako se iz `δξ` dobiva fizički budžet

Tek nakon što postoje fizički normaliziran model, definiran konačni vektor, kalibrirana nezavisna **eBOSS** kovarijanca `C`, nuisance matrica `B` i stvarno opravdana sistematska promjena `δy` istoga observablea, vrijedi simbolički

\[
W=C^{-1},\qquad
t_\perp=t-B(B^{T}WB)^{-1}B^{T}Wt,
\]
\[
\delta A=\frac{t_\perp^{T}W\,\delta y}{t_\perp^{T}W t_\perp},
\qquad
\sigma_A=(t_\perp^{T}Wt_\perp)^{-1/2}.
\]

`δA/σ_A` je odgovarajuća *algebarska* mjera pomaka procjene amplitude. Njeno numeričko ograničenje tek treba fizički i statistički definirati **prije** nove neovisne validacije, ne odabrati prema E4 originalnom odstupanju. Ako je `t_\perp` degeneriran ili je kovarijanca nevaljana, STOP. Devet eBOSS ID-jeva s rank≤8 i DESI DR1 120-mock kovarijanca nisu dopuštena zamjena za validan eBOSS `C`. E4 `2400→4800` razlike na istim mock galaksijama nisu neovisne realizacije kozmičke kovarijance.

## Što sada konkretno prolazi

[Source-only E6 izvršivi audit](../scripts/audit_eboss_dr16_a03_e6_theory_window_bridge.py) prije bilo kakvog testa SHA-verificira sve originalne roditeljske Git blabove, uključujući nepromijenjeni eBOSS conditional RR kod, originalni E4 i transparentni odd erratum, i detektira da je raspoloživi DESI wake predložak samo oblik. Zatim *isključivo na determinističkom sintetičkom pozitivnom 120×240 RR* verificira po 10 `6×120` blokova po kapi, svih 24 izlaza, even/odd razlaganje, linearitet, paritet pri refleksiji, nulti/konstantni even te nenulti radial-gradient finite-bin odgovor. Odvojeno provjerava samo **sintetičku algebru** nuisance-projected matched-filter `δA/σ_A`, uključujući singularne i degenerirane negativne kontrole. Odbija DESI `shape-only` kao apsolutni eBOSS model, proizvoljni effective-z, izostavljeni even član, pogrešni 18D rez i opaženi pristup.

[CI `36300514824` PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36300514824). On **ne izračunava fizički wake vektor** i ne čita originalni veliki FITS ili originalni eBOSS RR NPZ. Stoga se sada ne smije tiskati numerički `A_\mathrm{EV}`, `σ_A`, fizička tolerancija finite-random pogreške, eBOSS p-vrijednost ili exclusion/detection tvrdnja. A-03 fizički prozor, A-04 nezavisni eBOSS ansambl i opaženi odd seal ostaju otvoreni prema planu.
