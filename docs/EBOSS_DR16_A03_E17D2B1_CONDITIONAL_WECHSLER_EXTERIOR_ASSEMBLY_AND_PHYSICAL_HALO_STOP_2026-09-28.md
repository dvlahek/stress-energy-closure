# A-03E17D2b1 — uvjetni rast mase haloa i uzročno nedostajuća vanjska gravitacijska povijest

**28. 9. 2026. · Izvorni i neovisni source-only audit PASS. Puni fizički LRG×ELG halo/tracer bispektar BLOCKED, opaženi odd SEALED.** Ni jedna opažena galaksija, odd multipol, CLASS izračun, mock, novi science cut ili seed nije korišten; originalni E8–E17D2b0 roditelji i E14/E15 72/24 ostaju nepromijenjeni.

## Problem

E17D2a je numerički provjerio statičku radijalnu NFW gustoću, vanjski Newtonov potencijal i oba stvarna izvorna kratka kraka. E17D2b0 je iz [Dutton–Macciò 2014](https://doi.org/10.1093/mnras/stu742) objavljene *populacijske* c200c(M200c,z) relacije izračunao 18 uvjetnih snapshotova na dva deklarirana masena numerička referentna sidra, tri izvorna samo matematička z čvora i ilustrativnim koncentracijskim pomacima. Ti snapshotovi nisu praćeni jedan te isti halo i ne daju individualni M(a) ni derivaciju potencijala.

E17D2b1 na tu neidentifikabilnost odgovara ograničenim testom: dvije unaprijed definirane **različite** monotone funkcije prirasta mase imaju **isto** izvorno DM14 populacijsko sidro M200c, c200c i r200c pri z0=.95. Ispitujemo njihovu razliku izvan sidra i analitički dlnM/dlna i vanjski ∂Φ/∂z. To su uvjetne modelne putanje s nepoznatim parametrima, **ne rekonstrukcije stvarnih haloa**.

## Objavljeni funkcionalni oblik i stroga razlika konvencija

[Wechsler i sur., *Concentrations of Dark Halos from Their Assembly Histories*, Astrophys. J. 568 (2002) 52, DOI:10.1086/338765](https://arxiv.org/abs/astro-ph/0108151) u svojim numeričkim simulacijama razmatraju eksponencijalni oblik prirasta mase s karakterističnim parametrom vremena nastanka; [dokumentacija fizičkog modela Galacticus](https://galacticus.readthedocs.io/en/latest/physics/darkMatterHaloMassAccretionHistory.html) eksplicitno prikazuje isti oblik s referentnom epohom a0. E17D2b1 koristi samo **objavljeni funkcionalni oblik**, s unaprijed postavljenim slobodnim parametrom:

\[
M_{\alpha}(a)=M_0\exp\!\left[-\alpha\left(\frac{a_0}{a}-1\right)\right],
\quad a_0=\frac{1}{1+0.95},\quad\alpha=2a_c\in\{0.4,0.8\}.
\]

Originalni Wechslerovi prirasti mase i objavljena koncentracijsko-formacijska veza koriste **Mvir i cvir**. Naše fiksno DM14 sidro koristi **M200c i c200c**. Primjena istoga *oblika* M(a) na M200c stoga je **dodatna nekalibrirana matematička pretpostavka** i ne tvrdi se da je povijest empirijski validirana za tu masenu konvenciju. Niti se koristi ili konvertira Wechslerova cvir–a_c relacija u DM14 c200c. Vrijednosti α=.4 i .8 nisu posterior distribucija, mjerene epohe formacije niti statistička ograda LRG/ELG populacije; zaključane su kao dva različita izvora u QA testu.

[Prospektivni protokol E17D2b1](../source_data/eboss_dr16_a03_e17d2b1_conditional_wechsler_exterior_halo_assembly_prereg_2026-09-28.json), Git blob 28e33fad09db24b8b19739f50322169180442988, commit 25b6240, zaključan je prije novoga E17D2b1 računanja. SHA-pina originalni E8 4000q, originalni E16, E17D1, E17D2a, E17D2b0 joint i neovisni certifikat, E17D2b0 trajni manifest i originalni runner te točne izvorne E17D2b0 medijalne z=.95 case_04 i case_13 JSON bajtove.

## Samo izvana polje i tri unaprijed zaključana z čvora

Uz a=(1+z)⁻¹ ista uvjetna funkcija je

\[
M_\alpha(z)=M_0\exp\!\left[-\alpha\frac{z-z_0}{1+z_0}\right],
\quad \frac{d\ln M_\alpha}{dz}=-\frac{\alpha}{1+z_0},
\quad \frac{d\ln M_\alpha}{d\ln a}=\alpha\frac{1+z}{1+z_0}.
\]

Prijavljeni uzorci su **isključivo** originalni z=(.945,.95,.955), bez proširenja survey geometrije ili identifikacije tih čvorova kao novih eBOSS z_eff. Svaka od dvije izvorne masene reference 10¹² i 10¹³ hDM14⁻¹ M☉ dobiva dvije odvojene α putanje, četiri tročvorna testna modela.

Uvjetna DM14 ravna matter+Λ pozadina ima hDM14=.671, Ωm=.3175 i ΩΛ=.6825; ne pokreće se novi CLASS, a originalni hCLASS=.6736 i mν=.06 eV **nisu proglašeni kozmološki ekvivalentnima**. U svakom modelnom čvoru primjenjuje se definicija

\[
\rho_{\rm crit,DM14}(z)=\frac{3H_{\rm DM14}(z)^2}{8\pi G},
\qquad r_{200c,\rm phys}(z)=
\left[\frac{3M_\alpha(z)}{4\pi\,200\rho_{\rm crit,DM14}(z)}\right]^{1/3}.
\]

Za svaku masu testni fizički radijus R=3r200c(z0) je **fiksan** kroz tri čvora i provjeren da je strogo izvan halo r200c(z) za obje putanje. Za bilo koju sfernu raspodjelu čija je masa sadržana unutar r200c(z) tada Newtonov teorem o sfernim ljuskama na tom vanjskom radijusu daje trenutni prostorni potencijal i radijalnu gravitacijsku silu

\[
\Phi_{\rm ext}(R,z)=-\frac{GM_\alpha(z)}R,
\qquad g_r(R,z)=-\frac{GM_\alpha(z)}{R^2}.
\]

Izvodi i njihovi predznaci pripadaju **uvjetnom trenutnom Newtonovu vanjskom problemu pri fiksnom fizičkom R**. Primjerice dln|Φext|/dln a=dlnM/dln a, ali to nije puni FRW-gauge potencijal, samokonzistentna Einstein–Vlasov evolucija, valni Fourierov izvor ili retardirana sila na originalne neutrinske tokove. Fizička koncentracija c200c(z) **izvan sidra nije izračunana**, pa model ne izmišlja unutarnji NFW profil na drugim vremenima. Originalni E16 576 kutova samo su SHA-gated; novi kod ih ne kopira 576 puta za sferni vanjski problem i ne izmišlja two-leg bispektar.

## Numerika, neovisnost i rezultat

[Završni GitHub Actions CI 36387669876](https://github.com/dvlahek/stress-energy-closure/actions/runs/36387669876) **SUCCESS** nakon dvaju zabilježenih tehničkih failova. [Izvorni program](../scripts/audit_eboss_dr16_a03_e17d2b1_wechsler_conditional_exterior_assembly.py) provjerava SHA točnih originalnih sidara, četiri pozitivne monotone putanje, egzaktni povrat početne mase/radijusa, 200ρcrit mass–radius closure, exteriornost R i Newtonovu silu. Analitičke derivacije uspoređuje s centralnom konačnom razlikom na isključivo zaključanim z=.945,.95,.955. [Neovisni pure-stdlib replay](../scripts/audit_eboss_dr16_a03_e17d2b1_independent_wechsler_exterior_replay.py) zasebno parsira izvorna DM14 sidra i izvodi masu, H(z), ρcrit, r200c, Φext, gradijent i derivacije te odbija namjerno promijenjeni SHA.

Prvi [CI 36387598676](https://github.com/dvlahek/stress-energy-closure/actions/runs/36387598676) FAIL posljedica je stroge decimalne usporedbe 67.1/100 s .671 u neovisnoj preliminarnoj provjeri; tehnički fix 93419d4 prihvatio je samo binary-floating-point razliku 1e−15. Drugi [CI 36387635924](https://github.com/dvlahek/stress-energy-closure/actions/runs/36387635924) FAIL posljedica je razlike velikoga/maloga slova u provjeri teksta namjerno odbijenoga SHA; fix c8ebb39 promijenio je samo testni niz, ne provjeru sadržaja ili SHA. Izvorni parametrizirani modeli, izvori, fizikalni STOP i preregistrirani numerički pragovi nisu mijenjani.

Za **obje** početne mase rezultat je

| Matematički slučaj | α=.4 | α=.8 |
|---|---:|---:|
| M(z0) relativno prema istom originalnom sidru | 1 | 1 |
| c200c(z0) iz originalnoga DM14 median | isto sidro | isto sidro |
| Φext(R,z0) relativno prema istom originalnom sidru | 1 | 1 |
| dlnM/dln a pri z0 | **0.4** | **0.8** |
| M(z=.955)/M(z0) | 0.998974884764343 | 0.9979508203899325 |
| M(z=.945)/M(z0) | 1.0010261671752625 | 1.0020533873695967 |

Maksimalni prijavljeni izvorni mass/r200/exterior-reclosure relativni ostatak po četiri putanje je reda 1e−15; centralna razlika mase ima maksimalni očekivani O(Δz²) relativni ostatak **7.012932e−7** pri α=.8. [Neovisni certifikat](../source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28/e17d2b1_independent_stdlib_exterior_assembly_replay.json) potvrdio je sva četiri modela na 12 originalnih čvorova, najviši skalirani preklop **2.1410651029896144e−13**, a namjerno promijenjeni SHA odbijen je prije usporedbe fizike. Ovo je neovisna algebra i numerički replay istoga uvjetnog modela, ne nezavisna halo simulacija.

[Originalni joint izvještaj](../source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28/e17d2b1_original_joint_conditional_exterior_assembly_rates.json) SHA256 19cfd3e3a9382d07b4a978522e437780b79467dcbf8f8adcf45d3804b46ad3b9. [Trajni SHA manifest](../source_data/eboss_dr16_a03_e17d2b1_archived_CI_2026_09_28/archive_manifest.json), Git blob 2b6779fdcc39f9c6b5516b64044e6f942f93ebba, arhivira sva četiri originalna puna JSON-a, joint original i neovisni izvještaj na audit grani.

## Zaključak i precizni fizički STOP

Isto M(z0), c200c(z0), r200c(z0) i *vanjsko* Φext(R,z0) ne određuju dlnM/dln a čak ni u zadanoj jednoparametarskoj conditional Wechsler-formi bez dodatne kalibracije parametra rasta. Dokaz je ograničen na **dvije odabrane matematičke putanje**, ne implicira da obje odgovaraju statistički validiranim halo histories s istim c200c niti da je puni Einstein–Vlasov Cauchyjev problem nejedinstven. Posebno originalni Wechsler Mvir/cvir i naši DM14 M200c/c200c ne smiju se spojiti bez zasebne validacije.

Za stvarni fizički E17D2b potreban je vanjski halo formation/merger-tree ili validiran statistički mass-history prior kompatibilan s M200c/c200c i originalnom masivno-neutrinskom CLASS kozmologijom, high-z LRG/ELG halo mass/HOD/selection, evoluirajući unutarnji halo potencijal i konzistentan relativistički Einstein–Vlasov/FRW forced response na stvarnim vremenima. Tek nakon toga mogu se računati fizički long–short oba kratka tracer kraka i zasebni eBOSS 3pt prozor/kovarijanca. E17D1 bezdimenzijski u **nije** pretvoren u fizički cosmic time. Originalni 24D odd dvotočkasti vektor nije bispektar.

**Puni galaktički finite-K B, eBOSS ξ, S/N i detekcijska značajnost BLOCKED; opaženi odd SEALED; A03 pair-z RR, A04 independent 3pt covariance, SGC independent reverse i 48k finite-random još otvoreni. Main netaknut, PR #1 ostaje draft bez mergea.**
