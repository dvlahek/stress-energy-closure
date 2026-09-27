# A-03E17A — oba stvarna kratka kraka: originalni CLASS i fizička granica

**27. 9. 2026. · E17A linearni CLASS i neovisni numerički audit PASS. E17 puni finite-(K) fizički bispektar BLOCKED.** Ovo je nastavak [E16 točno zatvorenoga trokuta](EBOSS_DR16_A03_E16_EXACT_TRIANGLE_GEOMETRY_AND_FINITE_SQUEEZE_STOP_2026-09-27.md), a ne novi izvor neutrinskih distribucija, novi eBOSS estimator ili opservacijski rezultat.

## Problem i prethodno zaključana odluka

E14/E15 daju reducirani source-only (S(k,K,\ell)) uz izvorni (i\mu_{\rm long}\Delta b). E16 je pokazao da taj rezultat ne određuje dva stvarna kratka kraka niti jednoznačnu finite-(K) trotočkastu amplitudu. U osam od dvanaest izvornijih ((k,K)) kombinacija barem jedan krak izlazi iz skupa samo triju pohranjenih E14 (P_{cb}(k)) vrijednosti. Ekstrapolacija iz tih triju brojeva nije fizikalni izračun.

[E17 protokol](../source_data/eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk_prereg_2026-09-27.json) zaključen je **prije novoga CLASS izvršavanja** u [commitu a0911e3](https://github.com/dvlahek/stress-energy-closure/commit/a0911e3303b4620f51d3338d9be2db4147c45216). Neizmijenjeni su originalni E8 CSV SHA256 `bf8f48deb9514c5101d5ce314bfbd397485102a5487da1039f141e3cad6143d0`, E14 puni SHA256 `0f1efb45593cd51d9ffa1818157bdd52baba525a8fc4bb9cc30c12d98ba91fd8`, E15 puni SHA256 `6062aed4b581f2db89192bd1c1cf2639a74b6a313db9e10dca9e53041a95d983`, E16 originalni SHA256 `5d81f95a5571175f06eb032b865285ba34804701693f76af392be7fc31fb06db` i E16 neovisni SHA256 `8a6ee7610a0412de7bcaf5333950f41add0973a000bd245d49b41b77efc8e10d`. Originalni CLASS commit ostaje `e85808324f51fc694d12e3ed7439552a3c3f9540`. Smjer F+ i F−, masa 0,06 eV, deformacija 30 % i originalni LRG→ELG redoslijed nisu promijenjeni.

## Što je stvarno izračunano

Na istoj E16 geometriji vrijedi

\[
\mathbf k_1=\mathbf k-\mathbf K/2,\qquad
\mathbf k_2=-\mathbf k-\mathbf K/2,\qquad
\mathbf k_1+\mathbf k_2+\mathbf K=0.
\]

Sačuvano je svih **12 izvornih `(k,K)` parova × 48 orijentacija = 576 zatvorenih trokuta**. Za FD, F+ i F− zasebno, pri `z = 0.945, 0.950, 0.955`, iz izvornoga CLASS-a izravno je pozvan `pk_cb_lin(h |k_i|,z)` za **svaki od dvaju krakova**. Stoga je nedostajuća linearna (P_{cb}) podrška iz E16 dobivena novim izvornim CLASS pozivima; nije nadomještena interpolacijom ili ekstrapolacijom triju originalnih E14 vrijednosti.

Za oba stvarna `|k_i|` dodatno su izravno uzeti CLASS linearni transferi `t_ncdm[0]-t_cdm` i newtonovski `delta_cb`. Evaluirani su **linearnom interpolacijom isključivo unutar stvarno provjerenoga CLASS transfer-k grida od 81 čvora** pri svakom od tri z. Pohranjeni su parovi indeksa između kojih je evaluacija provedena, originalni LOS kutovi obaju krakova te E12 konvencija jednog linearnog brzinskog moda. *Single-mode transfer na kratkom kraku nije originalna R16-filtrirana koherentna dugovalna LOS brzina*. Ne tvrdimo da je 81-točkasti transfer grid zasebno numerički konvergirao ili da je nova interpolacija nelinearni halo odgovor.

Odvojeno, s izvornim state-specific E13 direct-vTk **R16 LOS rankom** i izvornom pozitivnom zamrznutom Eq20 okupacijom izračunata je uvjetna *kvazistatička, trenutna* faza na svakom stvarnom kraku, uključujući centriranu tri-z razliku `dphase/dln(a)`. Izvještaji to eksplicitno označavaju `Eq20_*_MODEL_ONLY`. Nije nametnut prosjek, polovica, simetrizacijski faktor ili linearna kombinacija krakova kao fizički (B).

## Numerički zaključak, neovisna provjera i trajni izvori

[Izvorni CLASS program](../scripts/audit_eboss_dr16_a03_e17_two_leg_CLASS_direct_vTk.py) i [odvojeni neovisni program](../scripts/audit_eboss_dr16_a03_e17_independent_fresh_CLASS_scalar_replay.py) ne dijele originalnu E17 implementaciju geometrije i interpolacije. Neovisni program sam konstruira svaki trokut, ponovno pokreće ista tri zamrznuta CLASS stanja i uspoređuje svih **1728 state×geometry redaka** odnosno **10368 `P_cb` i 10368 direct-vTk uzoraka**. Ovo je neovisno izvođenje numeričkih poziva i algebre, **ne neovisna fizikalna simulacija nelinearnog Vlasov halo sustava**.

[CI run 36340288916](https://github.com/dvlahek/stress-energy-closure/actions/runs/36340288916) završio je **SUCCESS** nakon tehničkoga ispravka naziva JSON polja u neovisnom čitaču. Raniji [CI 36340175889](https://github.com/dvlahek/stress-energy-closure/actions/runs/36340175889) nije proglašen PASS: izvorni CLASS dijelovi prošli su, ali neovisni čitač zatražio je `leg1_h_Mpc` umjesto zapisanoga `k1_h_Mpc`. Popravljen je samo čitač; fizikalni izvor, originalni E16 trokut i preregistrirani pragovi nisu mijenjani.

[Neovisni puni certifikat](../source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_independent_fresh_CLASS_scalar_replay.json) SHA256 `816570daa6acf3f81fffed75f4b3f5b3ba7e446630cfb93a28a1b7d1fdd974f6` potvrđuje:

- `max_CLASS_Pcb_relative_gap = 1.8839397268061077e-15`, `max_DIRECT_theta_relative_gap = 6.4358730989037295e-16` i `max_cb_density_relative_gap = 4.970422703779554e-16`;
- `max_triangle_closure_h_Mpc = 1.2643861424099487e-17` te `max_geometric_relative_gap = 2.920478102429548e-16`;
- originalna E14 tri-`k` (P_{cb}) sidra po stanju imaju najveći relativni gap **0**;
- diagnostic `K→10^{-5}K` za **linearni (P_{cb}) i oba kraka**, a ne za bispektar, daje najveći relativni gap `4.302313543387438e-7`.

[Izvorni objedinjeni E17 izvještaj](../source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/e17_original_joint_three_state_CLASS_both_short_legs_source_only.json) SHA256 `094c49cb89d6dd231a5cbc16ecdfe455bcb3938ddec375fd28c81c8098ab2ad9`. [Potpuni arhivski manifest](../source_data/eboss_dr16_a03_e17_archived_CI_2026_09_27/archive_manifest.json) trajno čuva izvornu FD datoteku SHA256 `cd89f798899e12d267fedd1b01d3ef1ddf8b1d8dbe5d14845b5e25d6c9c9af34`, F+ `4c95de74fd7ea03388cecb94bbe9fd06180865b045c5d68b925b69a3f2b3ba10`, F− `9ed7a5f219367570c37d01ba4b19031152b1930e551a2f77033fa9fa50b9ba52`, izvorni zajednički i neovisni izvještaj. Pet originalnih/nezavisnih JSON-ova nije ostalo samo u prolaznim CI artefaktima: [commit 85b543c](https://github.com/dvlahek/stress-energy-closure/commit/85b543c746ab667fd9e6d505356d38245e4b3cd6) trajno ih je pohranio u `source_data/` na audit grani.

## Simetrija, predznaci i granica fizičkoga zaključka

Zatvoreni trokut i podjela krakova izračunani su bez squeezed skraćenja. Zamjena `k1↔k2` permutira njihove stvarne valne brojeve i pripadne linearne CLASS vrijednosti. Realni, izotropni linearni CLASS (P_{cb}(-\mathbf k_i)=P_{cb}(\mathbf k_i)) zadovoljava potrebni reality uvjet. Iz E15 je zadržana algebraička činjenica da se formalni reducirani antisimetrični faktor `Delta_b(LRG,ELG)` mijenja u minus pri zamjeni LRG/ELG. Međutim, to **nije dokaz** predznaka ili Hermitian strukture *cijelog* galaktičkog bispektra bez fizikalne spregnutosti obaju krakova, selekcije i ostalih doprinosa. Niti provjera (P_{cb}) pri (K/k\to0) ne certificira squeezed limit nepoznatoga (B).

**Fizička prepreka ostaje točno lokalizirana:** E17A je zatvorio linearnu (P_{cb}) i direct-vTk podršku stvarnih krakova te uvjetni Eq20 model. Nije izveden nelinearni, vremenski retardirani Vlasov halo/galaxy odgovor na dugovalnu gustoću i vjetar, niti high-z LRG/ELG bias/HOD, evolution, magnification i relativistički nuisancei. Bez njih source-only `S` i dva linearna kraka ne određuju jedinstveni fizički finite-(K) (B(k_1,k_2,K)). Stoga su **fizički bispektar, eBOSS trotočkasti estimator/prozor, `xi`, signal-to-noise i značajnost BLOCKED**.

A03 fizički pair-redshift RR prozor, A04 neovisna eBOSS kovarijanca, originalni SGC neovisni reverse i finite-random konvergencija ostaju otvoreni. Devet originalnih matched mockova ne može dati punorangovnu kovarijancu izvornom 24D pilot-vektoru. Opaženi eBOSS odd ostaje **SEALED**, bez novih kataloga/mockova/seedova/cutova, bez kontakta s autorima i bez prijenosa DESI kovarijance. PR #1 ostaje draft; `main` nije mijenjan.
