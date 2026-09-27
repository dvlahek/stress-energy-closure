# A-03E17D2a — fizički NFW prostorni halo izvor, Poisson QA i granica retardiranog modela

**28. 9. 2026. · Izvorni i neovisni numerički prostorni model PASS; stvarni halo vremenski odziv i puni LRG×ELG finite-K bispektar BLOCKED.** Opaženi odd SEALED. Originalni E8–E17D1, E16 12×48 geometrija i E14/E15 72/24 nisu mijenjani.

## Znanstveni problem

E17D1 je numerički pokazao da dvije vanjske simboličke sile jednakoga ukupnog impulsa mogu dati različit retardirani odziv, čak i kada je originalna neutrinska distribucija i njezin ograničeni Vlasovljev kernel identičan. Ta demonstracija ne specificira stvarni halo. Prvi fizički uvjet koji se može zasebno ispuniti bez opaženog odd signala jest **pozitivna raspodjela halo mase, konačna ukupna masa i gravitacijski potencijal koji rješava Poissonovu jednadžbu**. Ona sama još ne specificira povijest formacije ili stvarnu populaciju galaksija.

E17D2a uzima literaturno definiran NFW prostorni profil [Navarro, Frenk & White 1997, DOI:10.1086/304888](https://doi.org/10.1086/304888), a za standardni normalizirani Fourierov transform radijalno odrezanog profila slijedi [Cooray & Sheth 2002, arXiv:astro-ph/0206508](https://arxiv.org/abs/astro-ph/0206508). Radijalni *truncation at r200* ovdje je unaprijed eksplicitna modelna pretpostavka; to nije potvrda da svaki stvarni eBOSS halo ima oštri rub.

[Prospektivni E17D2a protokol](../source_data/eboss_dr16_a03_e17d2a_truncated_nfw_spatial_poisson_halo_family_prereg_2026-09-28.json), Git blob e7e0281f0b6b205214db5faf03a7c364521b3f2c, commit c3f20ae, zaključan je prije prvog E17D2a numeričkog testa. Čuva originalni E8 4000q i E16 SHA, E17D0 originalni joint, E17D1 originalni/nezavisni SHA, E17D1 Git blob manifesta/protokola, izvorne FD/F+/F− i originalni LRG→ELG redoslijed. Nikakvi neutrinski izvori nisu ponovno fitani ili računani.

## Prostorna masa i Newtonov potencijal: potpuno uvjetan halo model

Neka su M200 ukupna pozitivna masa unutar r200, r_s>0 karakteristični radijus, c=r200/r_s>0 i y=r/r_s. Za 0<y<c:

\[
A(c)=\ln(1+c)-\frac{c}{1+c}, \qquad
\rho_{\rm NFW}(r)=\frac{M_{200}}{4\pi r_s^3 A(c)}
                   \frac{1}{y(1+y)^2},
\]
\[
\frac{M(<r)}{M_{200}}=\frac{A(y)}{A(c)},\qquad
A(y)=\ln(1+y)-\frac{y}{1+y}.
\]

Za y≥c ovdje je halo gustoća modelno odrezana i kumulativna masa jednaka M200. Oštri rez gustoće ne stvara površinski maseni delta-sloj, jer su potencijal i njegova prva radijalna derivacija kontinuirani na r200. Pozitivna radijalna gustoća ima integrabilni NFW centralni cusp, ali model ne evaluira divergentnu gustoću u r=0.

Uz uvjet da fizički Newtonov potencijal iščezava u beskonačnosti, za y≤c vrijedi

\[
\frac{\Phi(r)}{G M_{200}/r_s}
=-\frac{\ln(1+y)/y-1/(1+c)}{A(c)}.
\]

Za y≥c potencijal je -1/y u istoj normalizaciji. U središtu je potencijal konačan, -c/[(1+c)A(c)], a radijalna derivacija zadovoljava y² d[Φ/(GM200/r_s)]/dy=A(y)/A(c). Time je statička Poissonova jednadžba u ovom pojedinačnom sfernom profilu zadovoljena. **To nije vremenski ovisno rješenje Einstein–Vlasova ni izračunani halo neutrinski wake.**

## Prostorni Fourierov kernel: tri stvarna originalna vektora, bez bispektra

Za pozitivnu masu definiramo x=k_phys r_s,phys=k_com r_s,com i normalizirani *prostorni halo density profile*

\[
u(x,c)=\frac{1}{A(c)}
\int_0^c dy\frac{y}{(1+y)^2}\,
\operatorname{sinc}(xy),\qquad u(0,c)=1.
\]

Analitička formula s Si/Ci dana je u prospektivnom protokolu i [izvornom runneru](../scripts/audit_eboss_dr16_a03_e17d2a_nfw_spatial_poisson_halo.py). Za svaki *nenulti* fizički Fourierov mod statička Poissonova jednadžba daje

\[
\widetilde\Phi(k_{\rm phys})=
-\frac{4\pi G M_{200}}{k_{\rm phys}^{2}}\,
u(k_{\rm phys} r_{s,\rm phys},c).
\]

Kod sprema isključivo bezdimenzijski oblik \(-u(x,c)/x^2\), te originalne E16 |k1|, |k2| i K zasebno. **Tri spremljene prostorne težine nisu umnožene u B i ne opisuju stvarni long–short gravitacijski response.** Nulta Fourierova frekvencija potencijala za izoliranu masu ima standardnu 1/k² singularnost, pa se pri x=0 testira samo dobro definirana u(0)=1.

Šest unaprijed odabranih numeričkih QA slučajeva c=(2,4,8) i α=k0 r_s,com=(0.1,1), uz k0=0.05 h/Mpc, definira x_i=α|k_i|/k0 i x_L=α K/k0. **Ovo nisu** mjerenja koncentracije, r_s, M200, redshift ili fiducijalna HOD/halo kalibracija. Svaki originalni slučaj koristi svih 576 E16 geometrija, oba egzaktna kratka kraka i originalni dugi K: 6×576×3=10.368 zasebnih prostornih vrijednosti u(x,c).

## Preregistrirani QA i neovisni audit

Originalni program koristi analitičku Si/Ci NFW formulu, zatvorene E16 vektore, monotoni kumulativni halo M(<r), unutarnji i vanjski potencijal, izostanak površinskog singularnog sloja i kontinuitet radijalne sile. U izvornom frozen E16 skupu nema svih zrcaljenih kutova (npr. -0.6), pa kratko-kraku simetriju za te nedostajuće orijentacije račun provjerava **analitičkim odrazom**, bez dodavanja novih science orijentacija. Pozitivna radijalna masa jamči |u|≤1, ali ne zahtijeva u≥0 za sve x zbog mogućeg Fourierova osciliranja.

[Neovisni stdlib program](../scripts/audit_eboss_dr16_a03_e17d2a_independent_nfw_radial_replay.py), bez uvoza izvornoga programa, NumPyja ili SciPyja, ponovno izrađuje originalne E16 analitičke modove i zasebno izravno integrira pozitivan NFW radijalni Fourierov integrand 4096-panelnom Simpsonovom kvadraturom. Za svaki jedinstven radijalni x zasebno provjerava i 2048/4096 podjelu. Provjerava i originalne potpune SHA izvještaje, svih 10.368 short/short/long vrijednosti, radijalnu masu/potencijal te neovisnu negativnu kontrolu namjerno krivoga SHA otiska. Ovo je numerički neovisan *prostorni Poisson/NFW* audit, ne druga cosmological N-body ili Einstein–Vlasov simulacija.

Prvi [CI 36354065327](https://github.com/dvlahek/stress-energy-closure/actions/runs/36354065327) ostaje **FAIL**. Tehnički uzrok: provjera kumulativne mase u originu pokušala je izračunati A(0) preko helpera koji ispravno zahtijeva strogo pozitivan argument koncentracije c, iako je egzaktno M(<0)=0. Jedini ispravak, commit 7e099bd, dodao je egzaktni M(<0)=0 ogranak. Nije promijenio matematički model, koncentracije, α, izvorne E16 modove, ranije izvore ni prospektivne QA pragove.

Završni [CI 36354103086](https://github.com/dvlahek/stress-energy-closure/actions/runs/36354103086) **SUCCESS**. Šest slučajeva i neovisni full-original 10.368-mode replay prošli su preregistrirane kriterije. Neovisna maksimalna skalirana razlika originalne analitičke Si/Ci formule i izravne 4096-panelne radijalne integracije je **1.2443794644456727e−12**. Kontinuitet radijalnoga potencijala i sile na r200 u originalnim analitičkim QA uvjetima ima odstupanje 0, originalni najviši E16 trokutni closure 1.2643861424099487e−17 h/Mpc, a kratko-kraku geometrijski/refleksijski ostatak po numeričkim slučajevima ne prelazi 1.887380e−15. Maksimalno odstupanje originalne male-x druge-momentne dijagnostike i analitičkog transformera zabilježeno je zasebno kao numerički test, ne kao fizikalni halo error budget.

[Trajni manifest svih osam rezultata](../source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28/archive_manifest.json), Git blob ccb17f155211d88b943306bf5d31a35ffb4811f3, trajno SHA-arhivira šest originalnih punih c/α JSON izvještaja, izvorni [joint](../source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28/e17d2a_original_joint_nfw_static_spatial_poisson_family.json) SHA256 040c9711bfdfee6682acda4b48bc578b64e88c11d143d35a9a9ae46ff4c74a31 i neovisni [certifikat](../source_data/eboss_dr16_a03_e17d2a_archived_CI_2026_09_28/e17d2a_independent_stdlib_direct_radial_fourier_replay.json) SHA256 436d689762582e74fda83e3df0f1752fe2bbdf595b2c189f13aabce018ef2287. To su izvorni izračunati bajtovi i neovisni replay u source_data/ audit grane, ne privremeni Actions artefakti.

## Fizički STOP

Dovršeni E17D2a specificira samo **uvjetnu familiju statičkih prostorno konzistentnih pojedinačnih halo profila**. Iz NFW oblika sami po sebi ne proizlaze M200, c(M,z), stvarni r_s, populacija/evolucija LRG/ELG, halo formation history i vremenska funkcija Ψ(k,η). E17D2a ne tvrdi da je NFW dokazan jedinstven opis svakog eBOSS haloa, ne promiče E17D1 simboličke funkcije g_E/g_L u fizičku povijest i ne kombinira prostorni u(x,c) s E17D0 source-om da stvori navodni fizički B.

Prije stvarnoga E17D2b fizičkog retarded halo modela trebaju neovisno fizički definirani ili kalibrirani M200 i c(M,z), halo initial/formation conditions, evoluirajući gravitacijski potencijal uz usklađenu Einstein–Vlasov/FRW gauge dinamiku, stvarni second-order long–short coupling za **oba** kratka kraka i high-z LRG/ELG HOD/bias/evolution/magnification/relativistički nuisancei. Zatim zasebno trebaju stvarni eBOSS physical triple-selection/window i neovisna 3pt kovarijanca istoga estimatora. Originalni 24D odd dvotočkasti vektor nije bispektar. A03 pair-z RR/A04, originalni SGC independent reverse i 48k random konvergencija ostaju otvoreni.

**Puni fizikalni galaktički B, eBOSS ξ, S/N, detekcijska značajnost i unblinding BLOCKED; opaženi odd SEALED. Bez novih CLASS, FITS, kataloga, mockova, seedova, science cuts ili kontakta s autorima. Main netaknut; draft PR #1 bez mergea.**
