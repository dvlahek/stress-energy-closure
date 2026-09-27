# A-03E8c — bezdimenzijski odziv jednog moda i provjera normalizacije objavljene formule

**Datum i status, 27. 9. 2026.** Ovo je retrospektivna, original-SHA vezana *source-only* fizikalna provjera E8/E8b izvora, **ne** novi CLASS/FITS izračun, ne konačni eBOSS signal i ne novi acceptance. [Originalni E8/E8b sažetak](../source_data/eboss_dr16_a03_e8_frozen_kinetic_resonance_and_static_SI_green_result_2026-09-27.json), [točni erratum i izvorni Git blabove](../source_data/eboss_dr16_a03_e8c_exact_eq20_phase_and_bibliography_erratum_2026-09-27.json), [ponovljiva E8c skripta](../scripts/audit_eboss_dr16_a03_e8c_conditional_phase_and_eq20_conventions.py) i [manifest svih devet ranijih ilustrativnih slučajeva](../source_data/eboss_dr16_a03_e8c_original_nine_case_conditional_phase_source_only_manifest_2026-09-27.json) čine jednu zaključanu dijagnostiku. [E8c CI 36315647049 PASS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36315647049) provjerava SHA, fizički/komovirajući `k`, sign-reversal i odbijanje neoznačenog `k`; zasebni pregled 9/9 korisničkih E8+E8b JSON slučajeva ima SHA256 rezultata `58cbb49e1d910d1dd102644dbcec2845c04987ac2fb87cb10ef19a4d10db0f79`.

## Što je izračunato i zašto je to korisno

Za izvorni jedan masivni ncdm, isti izotropni `F_\pm` i *fiksiran* statički halo-potencijal, komponentu relativne brzine i **fizički** valni broj, [Okoli et al. 2017 Eq. (20), prvi točan izraz](https://doi.org/10.1093/mnras/stx560) daje gravitacijsko polje wakea `g_{\nu,k}=K_\nu\Phi_k` duž `k`. Originalno halo gravitacijsko polje iste ravninske mode ima amplitudu `k_phys |Phi_k|`. Stoga definiramo

\[
\alpha_s(k_{\rm phys},z,v_\parallel)\equiv
\frac{K_s(k_{\rm phys},z,v_\parallel)}{k_{\rm phys}}
=
\frac{2NGm_\nu^4\,v_\parallel}{\hbar^3 k_{\rm phys}^2}
\,f_s(q_*),\quad
q_*=\frac{m_\nu |v_\parallel|}{c\,T_{\nu0}(1+z)}.
\]

`s=FD,+,-` označuje originalnu **okupaciju bez zajedničkog CLASS `2/(2pi)^3` prefaktora**; taj prefaktor u originalnim omjerima `F_\pm/F_0` egzaktno otpada. Za jedan Fourierov mod `g_\nu/g_h` je u linearnoj perturbaciji kvadraturno (imaginarno) pomaknut za `i\alpha`, uz predznak ovisan o definiciji Fourierove transformacije. **Ovo je lokalni, trenutačni per-unit-potential statički koeficijent, ne galaktički RSD fazni parametar nakon evolucije.**

Kada je `a=1/(1+z)`, `k_phys=k_com/a`, iz iste jednadžbe slijedi stroga, **isti-mod** algebarska jednakost

\[
\alpha_s=
\frac{2NGm_\nu^4 a^2v_\parallel}{\hbar^3 k_{\rm com}^2}
\,f_s(q_*).
\]

Ista je to `a²/k_com²` struktura [objavljene Eq. (32)–(33)](https://academic.oup.com/mnras/article/468/2/2164/3063204). Ali ne smije se bez provjere zamijeniti fizička oznaka `0.05\,h/Mpc` komovirajućom oznakom iste numeričke vrijednosti. Različiti su modovi.

| Ilustrativni izvorni `z=.95,v=200 km/s`, `k_phys=.05 h/Mpc` | FD | F+ | F− | F+ − F− |
|---|---:|---:|---:|---:|
| `alpha_s` | 0.0011741803783 | 0.0009567193488 | 0.0013916414078 | −0.0004349220590 |

Ako **umjesto toga** `0.05 h/Mpc` označuje `k_com` pri `z=.95`, tada je `k_phys=.0975 h/Mpc` i FD `alpha=0.0003087916840`. To je matematička `(1+.95)^{-2}=0.2629848784` promjena pri jednakoj fizičkoj `v` i `F(q_*)`, **ne mjerenje z-ovisnosti wakea**. Izvorni `z=.95`, `v=200` i taj `k` ostaju samo dijagnostički primjeri, ne izmjereni parametri odabranih eBOSS galaksija.

## Bibliografski i fizički erratum — ne miješati `f_{\rm FD}` i autorov `mu`

Točan DOI objavljenoga rada Okoli, Scrimgeour, Afshordi, Hudson, MNRAS 468 (2017) 2164–2175, jest **[10.1093/mnras/stx560](https://doi.org/10.1093/mnras/stx560)**. Originalni SHA-zaključani E8 protokol krivo navodi `stx539`; original **ne mijenjamo**, nego ovu ispravku dodajemo kao naknadni audit.

Ozbiljnija je granica normalizacije. **Doslovno prvi** oblik njihove Eq. (20) ima `f_FD(q_*)=1/(exp(q_*)+1)`. Njihova sljedeća *približna* Eq. (20) piše faktor `mu` i navodi `0.7 ≲ mu < 1`, a kasniji forecast uzima `mu=1`. Za fizikalni `q_*>=0`, `f_FD<=1/2`. Zato iz prikazanih jednadžbi **ne proizlazi numerička jednakost** tih koeficijenata, osim uz dodatnu nenavedenu definiciju ili konverziju faktora. Ne pretpostavljati da je autorima sigurno promaknuo faktor dva, niti proizvoljno namještati naš `N` ili `F_\pm` kako bismo dobili `mu=1`. Pri izvornom `q_*=.122047` je `f_FD=.4695260669`, pa bi stroga formalna zamjena `f_FD → mu=1` skalirala instantni kernel za **2.1298072045**, bez dokazivanja da je to ispravna fizička kalibracija. Naši E8b SI brojevi ostaju jasno označeni kao **točna occupancy grana prvoga Eq. (20)**, ne kao potvrda autorove približne kalibracije.

## Što je još nužno da nastane LRG×ELG dipol

Radovljeva Eq. (32)–(34) povezuje **vremenski ovisnu** fazu s imaginarnim RSD cross-powerom. U našem izračunu nedostaju posebno `v_{\nu c}(k,z)` i `d\alpha/dt` za evoluirane originalne `F_\pm`, zadana fizička LRG/ELG bias/evolution/magnification veza i isti eBOSS tracer-specific, pair-z/depth/chunk uvjetovan empirical RR operator. `alpha` i izvorna `F_+/F_-` razlika po jednom modu **nisu** `xi_1`, `xi_3`, `A_{\rm EV}` niti `S/N`. Dodatno A-04 traži stvarnu, nezavisnu eBOSS kovarijancu; devet izvornih mock ID-jeva ima rang centrirane kovarijance najviše osam. Originalni SGC independent reverse `R1D2` problem ostaje numerički otvoren. 48k random konvergencija također nije zasebno certificirana.

**Odluka:** nema novih full-pair WSL recountova, novih post-pilot seedova ili mock preuzimanja, nema otvaranja observed odd, niti prepisivanja stare E8/E8b fizike. Ovo je nova audit razina na originalnom draft PR-u i grani, `main` ostaje netaknut.
