# E26 — originalni adijabatski density–wind rang i granica Dopplerove identifikabilnosti

**29. 9. 2026.** Izvorno utemeljen linearni rezultat. Nije fizička galaktička prognoza ni eBOSS detekcija.

## 1. Što daju originalni E12/E13/E21

Za originalni pojedinačni adijabatski primordijalni mod R(k) na fiksnom F,K,z=.95 imamo delta_cb=D_F R i v_nucdm,LOS=i mu V_F R. Izvorni direct-vTk filter definira V_F=-299792.458 theta_rel,F W_R16/k_Mpc u km/s. Originalni E12 daje P_R=2 pi² Delta_R²/k_Mpc³, E13 zasebno daje svježi D,theta_rel,W, a E21 originalni C_v. U originalnoj Fourierovoj konvenciji P_delta,v=+i mu C_v uz C_v=-D V P_R.

Zasebno rekonstruirani P_delta,delta=D² P_R i P_v,v=mu² V² P_R daju egzaktno det P= P_delta,delta P_v,v - mu² C_v²=0. To je rang jedan *izvornog zajedničkog linearnog adijabatskog moda*, ne tvrdnja o stvarnim nelinearnim halovima, šumu, nezavisnim primordijalnim modovima ili svim opažajnim poljima.

## 2. Provjera na originalnim izvorima

[E26 protokol](../source_data/eboss_dr16_a03_e26_original_linear_source_rank_protocol_2026-09-29.json) Git blob 4b09681d188bd324b9d9f3ba459bd5608ba41195 zaključan je prije nove numeričke QA, ali nakon znanstvenog izvoda. [Izvorni E26 kod](../scripts/audit_eboss_dr16_a03_e26_original_adiabatic_density_wind_rank.py) ne definira P_vv iz E21 crossa: primordijalnu amplitudu uzima iz originalnoga E12, a zasebni D,theta,W iz originalnoga E13, pa cross provjerava naspram originalnoga E21. [CI 36596911178](https://github.com/dvlahek/stress-energy-closure/actions/runs/36596911178) SUCCESS na originalnim FD/plus/minus i K=.001,.002,.003,.005 h/Mpc za mu=1,.5,0. Maksimalni relativni det reziduali: FD 6.988414963062396e-16, plus 2.3199683601096716e-16, minus 8.198357900668836e-16. Svih osam negativnih kontrola prošlo je. Rezultat je source-only, bez novih CLASS ili opaženih podataka.

## 3. Standardni Doppler i uvjet zasebnoga fizičkog signala

Rekonstrukcija v_nucdm=i mu(V/D)delta_cb iz ISTE linearne density realizacije ne dodaje novu nezavisnu početnu fazu na istom K. Ipak, dva fizički opravdana *predloška* kroz više K/z mogu biti linearno nezavisna i kada dijele primordijalnu fazu. Uz nuisance stupce N, neutrinski template t_nu nakon stvarnoga prozora može se odvojiti samo ako rank([N,t_nu])>rank(N) te postoji valjana kovarijanca istoga estimatora. Fizikalni koeficijenti galaxy responsea iz E25 moraju biti nezavisno kalibrirani.

E12 originalni RUNNER čita t_cdm, ali arhivirani četiri-K E12 JSON čuva theta_nu-theta_cdm, ne zaseban numerički t_cdm. Zato standardni CB-Doppler predložak i njegov numerički preklop s neutrinskim kroz originalni K/z i 24D prozor ovdje NISU izračunani. Ne uvoziti proizvoljni 1/K, ne interpolirati četiri duga K na puni eBOSS odd i ne proglasiti nul-determinantu dvaju izvornih polja opservacijskom nerazlučivošću.

## Fizički STOP

Prije eBOSS inference trebaju high-z halo+environment povijest, pravi neutrinski beta/epsilon/selection i baseline galaxy Doppler, puna k/z i signed parna selekcija te neovisna A04 kovarijanca. E20 kvadratni mixed-bispectrum kanal je zaseban. Originalni observed 24D odd SEALED; A03 physical template i A04 inferencija BLOCKED. Nema novih WSL, CLASS, FITS, Abacus, mockova, cutova/seedova ni mergea draft PR #1. Main se ne mijenja.
