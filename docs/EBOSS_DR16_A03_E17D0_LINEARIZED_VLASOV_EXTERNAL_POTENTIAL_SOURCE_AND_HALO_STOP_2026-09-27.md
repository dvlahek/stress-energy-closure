# A-03E17D0 — izvorni Vlasovljev odziv po jedinici vanjskog potencijala, ne puni halo bispektar

**27. 9. 2026. · Završni originalni i neovisni CI PASS; puni finite-K retardirani halo/tracer bispektar BLOCKED.** E17D0 nastavlja E17C bez promjene E8–E17C, svih 576 E16 trokuta i triju originalnih F distribucija. Opaženi eBOSS odd ostaje SEALED.

## Problem i fizička granica

E17C je dobio homogeni bezsudarni propagator, ali slobodan let ne uključuje gravitacijsku silu. E17A/E17B poznaju originalni CLASS linearni P_cb i direct-vTk oba stvarna kratka kraka, ali ne određuju vanjski potencijal haloa i njegovu vremensku povijest niti LRG/ELG spregu. E17D0 izdvaja sljedeći izvodiv korak: **linearnu relativnu neutrinsku gustoću kao funkcional vanjskog Newtonova potencijala**, bez pretpostavljanja njegove vrijednosti.

[Prospektivni E17D0 protokol](../source_data/eboss_dr16_a03_e17d0_linearized_vlasov_external_potential_source_prereg_2026-09-27.json), Git blob 66df071f0750201e792ed2b3ffc81e28b13dfce6, commit e58663e4, zaključan je prije novog E17D0 računa. Originalni 4000q F0/F+/F− SHA, E16, E17A, E17B i E17C originalni/nezavisni roditelji, 12×48 geometrija i LRG→ELG redoslijed ostaju nepromijenjeni. Nema novih CLASS, FITS/mock podataka, seedova, rezova ili tuninga.

## Izvedena jednadžba i točno računani objekt

U vodećem nerelativističkom Newtonovu opisu na fiksnoj pozadinskoj epohi, s originalnom bezdimenzijskom komoving q varijablom i vanjskim **simboličkim** Newtonovim potencijalom ψ_ext, linearizirana Vlasovljeva jednadžba ima oblik

\[
\partial_\eta\delta f_F
+i k\mu\frac{qT_{\nu0}}{a m_\nu}\delta f_F
=i\,\frac{a m_\nu}{T_{\nu0}}k\mu\,F'(q)\,\psi_{\rm ext}(k,\eta).
\]

Ovo je ograničena Newtonova vanjsko-silna specijalizacija; puna relativistička masivno-neutrinska hijerarhija uključuje vremenski razvoj metričkih potencijala i konzistentnu Einsteinovu spregu. Referenca za širi formalizam: Ma i Bertschinger, *Astrophysical Journal* 455 (1995), doi:10.1086/176550, [astro-ph/9506072](https://arxiv.org/abs/astro-ph/9506072). Ne tvrdi se da E17D0 implementira njihovo puno relativističko rješenje.

Za impuls zadanoga ψ_ext homogeni propagator iz E17C daje kutno integriranu bezdimenzijsku težinu po jedinici simboličkog impulsa potencijala:

\[
H_F(s)=-
\frac{\displaystyle\int_{q_{\min}}^{q_{\max}}
 dq\,q^2F'(q)\,\operatorname{sinc}'(sq)}
{\displaystyle\int_{q_{\min}}^{q_{\max}}dq\,q^2 F(q)}.
\]

sinc(x)=sin(x)/x. Impuls bi nosio zaseban faktor \(k\,a m_\nu/T_{\nu0}\int d\eta\,\psi_{\rm ext}\), ali mu **nije dodijeljena fizička amplituda ni halo povijest**. E17D0 ne uzima proizvoljni potencijal, ne spaja H_F s E14/E15 S i ne proglašava H_F fizičkim bispektrom.

Sva F stanja potječu iz originalnog E8 CSV-a s q∈[0,20] i 4000 točaka. Račun koristi originalni komadno-linearni F(q) unutar izvornog raspona, njegove egzaktne intervalne nagibe i unaprijed zaključanu GL2 kvadraturu, bez ekstrapolacije, glatkanja ili reoptimiranja. Za svaku od 576 geometrija izravno koristi oba stvarna |k1|, |k2| i izvorne isključivo dijagnostičke bezdimenzijske s0=(0,0.25,1), s_i=s0 |k_i|/(0.05 h/Mpc). Ti s_i nisu stvarni vremenski lagovi, a z=.945/.95/.955 iz E17A nije pretvoren u novu E17D0 vremensku evoluciju.

## Neovisna kontrola konačnoga q raspona

Neovisni audit računa matematički istu težinu **bez diferenciranja F** pomoću integracije po dijelovima:

\[
H_F(s)=
\frac{-\bigl[q^2F(q)\,\operatorname{sinc}'(sq)\bigr]_{q_{\min}}^{q_{\max}}
+\displaystyle\int_{q_{\min}}^{q_{\max}}\!dq\,F(q)
\left[2q\,\operatorname{sinc}'(sq)+s q^2\,\operatorname{sinc}''(sq)\right]}
{\displaystyle\int_{q_{\min}}^{q_{\max}}\!dq\,q^2F(q)}.
\]

**Rub pri q_max=20 eksplicitno je zadržan.** Originalna konačna mreža nije proizvoljno proglašena beskonačnom ni njezin rub izbačen bez dokaza. Na istom rasponu analitička kratko-vremenska kontrola jest \(H_F(s)/s\to -1 + [q^3F(q)]_{\min}^{\max}/(3I_F)\). Nulta faza, H_F(-s)=-H_F(s), kratko-kraku zamjena i matematički K/k→0 test zaključani su prije računanja. Neovisni [stdlib audit](../scripts/audit_eboss_dr16_a03_e17d0_independent_finite_boundary_replay.py) izvorni CSV parsira zasebno i rekonstruira analitičke E16 krakove bez uvoza originalnog [E17D0 koda](../scripts/audit_eboss_dr16_a03_e17d0_external_potential_vlasov_source.py), NumPyja ili CLASS-a. Ovo nije drugi Einstein–Vlasov solver.

## Originalni i neovisni rezultati

[CI 36348300285](https://github.com/dvlahek/stress-energy-closure/actions/runs/36348300285) **SUCCESS**. Sve tri originalne distribucije prošle su 3×576×2×3 = **10.368** fixed original state/leg/phase vrijednosti. Maksimalni QA preko triju stanja:

| Kontrola | Najveći rezultat | Opseg |
|---|---:|---|
| Nulta faza H_F(0) | 0 | Originalni ograničeni funkcional |
| H_F(s)+H_F(−s) | 0 | Matematička odd-s kontrola |
| Kratko-kraka zamjena | 2.775558×10⁻¹⁷ | Originalni E16 sampled/analitički odraz |
| Trokutni closure | 1.264386×10⁻¹⁷ h/Mpc | E16 originalna geometrija |
| Squeezed izvorna matematička dijagnostika | 1.669617×10⁻⁶ | k·(1±ε/2), ε=10⁻⁵ |
| Kratko-s identitet s konačnim rubom | 2.171836×10⁻⁸ | Skalirani originalni QA |
| F+/F− sredina ne-normaliziranog izvora prema FD | 8.673617×10⁻¹⁹ | Linearnost originalnoga F para |
| Neovisna finite-q integracija po dijelovima | 5.163541×10⁻¹³ | Svih 10.368 izvornih uzoraka |

Neovisni audit prošao je i SHA-tamper negativnu kontrolu. Originalni I_F nazivnici za FD, F+ i F− nisu matematički identični do svih znamenaka, pa se **linearnost ne-normaliziranoga izvora ne smije automatski prepisati kao identitet zasebno normaliziranih H_F**. To je razlog za odvojenu evidenciju J_F i H_F.

[Izvorni zajednički izvještaj](../source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27/e17d0_original_joint_external_potential_source_only.json) SHA256 b9c8f21973755c64d8708a4ccb4d86ee56bbc12d448ffaa098b7495eb85661ae; [neovisni izvještaj](../source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27/e17d0_independent_finite_boundary_scalar_replay.json) SHA256 2258c53b35d62ec9f00475894af7fcffdfb1cd98ecc7badef5bdfc131ca7e1af. [Trajni SHA manifest](../source_data/eboss_dr16_a03_e17d0_archived_CI_2026_09_27/archive_manifest.json), Git blob 705948d755b16ce4f8b5bd9b32d444221f064bee, čuva originalna tri puna state izvještaja, zajednički izvještaj i neovisni certifikat na audit grani.

## Fizički STOP i sljedeći uvjeti

E17D0 dokazuje i računa samo ograničeni linearni **odgovor po jedinici vanjskog Newtonova potencijala** na originalnoj zamrznutoj distribuciji. Nema specificiranog stvarnog ψ_ext(k,η), halo formacijske/masene povijesti, fizičkih vremenskih lagova, FRW pune evolucije, samokonzistentne Einstein–Vlasov metrike, nelinearne long–short spregne funkcije ni LRG/ELG high-z HOD/bias/evolution/magnification/relativističke kalibracije. Stoga dva izvorno matematički dopuštena finite-K nastavka iz E16 nisu fizički razlučena ovim računom; E14/E15 S i 72/24 izvora ostaju neizmijenjeni.

Za puni E17D potreban je nezavisno fizički specificiran vanjski/halo potencijal s dinamičkim izvorom, zatim relativistički/retardirani Vlasov i self-gravity, oba kratkovalna halo/tracer kraka, redshift/selection i zasebna survey 3pt kovarijanca. A03 fizički pair-z RR, A04 neovisna eBOSS kovarijanca, SGC independent-reverse i 48k finite-random konvergencija i dalje su otvoreni. Originalni eBOSS 24D dvotočkasti odd **nije** trotočkasti estimator.

**Opaženi odd SEALED; nikakav fizički galaktički B, eBOSS ξ, S/N, significance ili unblinding nije izveden.** Main je netaknut i PR #1 ostaje draft.
