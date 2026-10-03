# E17D2b3b0 — uski pristup AbacusSummit podacima: zaključani OFFLINE protokol, bez preuzimanja

**28. 9. 2026. · [CI 36443992886 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36443992886) isključivo za SHA-pinnani offline data-access gate. Stvarni podaci: BLOCKED_DOCUMENTATION_ONLY_NO_REAL_ABACUS_DATA. Observed odd SEALED. Main netaknut, PR #1 draft.**

## Znanstveni razlog i praktična odluka

E17D2b2 je pokazao da dokumentirani AbacusSummit c000 ima šest nominalno podudarnih CLASS kozmoloških parametara, ali drukčiji As i tau, glatku neutrinsku N-body aproksimaciju bez F+/F− wakea i CompaSO L1 mass konvenciju nejednaku M200c. E17D2b3a pokazao je na dvije pozitivne **idealizirane** gustoće da čak savršena L1 masa i SO radijus ne određuju M200c bez informacija o unutarnjem radijalnom profilu. Sljedeći korak mora precizirati koji se stvarni podaci mogu zakonito i znanstveno smisleno upotrijebiti.

[Službeni AbacusSummit data-access opis](https://abacussummit.readthedocs.io/en/latest/data-access.html) navodi raw/off-machine originalni release **DOI 10.13139/OLCF/1811689**, očišćene halo kataloge s pripadajućim auxiliary datotekama **DOI 10.13139/OLCF/1828535** i halo light-cone **DOI 10.13139/OLCF/1825069**, koji **nije** zamjena za snapshot merger tree. Na OLCF tape arhivi puni halo-paket jedne base simulacije iznosi približno **6,6 TB**. [Javna Abacus procjena prostora](https://abacussummit.readthedocs.io/en/latest/disk-space.html) za jednu base simulaciju navodi približno **75 755 MB** halo_info za cjeloviti direktorij z0.950 i **101 105 MB** za sve sekundarne halo proizvode te epohe. To su indikativni ukupni iznosi iz službene dokumentacije, **ne** stvarna veličina određene odabrane superslab datoteke, stvarna raspoloživost konkretnog proizvoda ili lokalni slobodni prostor.

OLCF/Constellation može zahtijevati izvoz velike arhive. NERSC/Globus uz odgovarajući pristup može omogućiti uže datoteke, ali točan item, box/phase, tree file i dostupnost tek moraju biti stvarno provjereni. [abacusutils CompaSO upute](https://abacusutils.readthedocs.io/en/latest/compaso.html) dokumentiraju učitavanje **jednog superslaba** iz stvarno pronađenoga halo_info_XXX.asdf te samo potrebnih stupaca argumentom fields. Očišćeni katalog mora imati usklađene cleaning datoteke, a jedno učitano halo_info polje nije automatski dovoljno za radijalnu 200-kritičnu remeasurement. [Službena specifikacija sekundarnih izlaza](https://abacussummit.readthedocs.io/en/latest/data-products.html) objašnjava da približni z0.950 sadrži halo_info i PID subsample, bez RV/field particle datoteka; stvarni redshift dolazi **iz konkretnog ASDF zaglavlja**, a SODensityL1 je epošno definiran prema srednjoj gustoći, ne automatski 200 puta kritičnoj.

## Izvršen prospektivni dry-run

[Protokol](../source_data/eboss_dr16_a03_e17d2b3b0_abacus_real_data_access_dryrun_prereg_2026-09-28.json), Git blob **49df0cd9c110b3ef9db02a3ad0cb9a4ff20554e1**, commit **dd94d4c**, zaključan prije numeričkoga/CI izvođenja. SHA-pina originalne E8 4000q, E16 576 kutnih geometrija, E17A CLASS oba kratka kraka, E17D2b2 originalnu metadata provjeru, E17D2b3a originalni+independent dokaz i njegov manifest/protokol. Preregistrirani [BLOCKED predložak kandidata](../source_data/eboss_dr16_a03_e17d2b3b0_abacus_access_candidate_BLOCKED_template_2026-09-28.json), SHA256 **53e4d85503192ee082e124c721a9eb04fc3ea82a573600991849e29f6da49597**, navodi AbacusSummit_base_c000_ph000 **samo kao kandidata za inventar**, ne kao validiranu znanstvenu fazu. Oznaka z0.950 nije stvarni header redshift. Raw ili cleaned release, Globus path, točna superslab/tree datoteka, byte size, checksumovi i lokalni disk nisu izmišljeni ili upisani kao stvarni.

[Izvršivi Python protokol](../scripts/audit_eboss_dr16_a03_e17d2b3b0_access_contract.py) ne sadrži downloader, HTTP, ASDF parser, halo loader ili CLASS. Protokol provjerava zaključane roditeljske SHA-ove i točan predložak. Negativne kontrole odbijaju pretvaranje direktorija z0.950 u dokaz stvarnoga redshifta, zlib CRC32 kao zamjenu za **GNU/POSIX cksum**, L1 kao automatski M200c, fake real ASDF dostupnost, pokušaj tvrdnje da je izračunat fizički B ili otvaranje observed odd. Čak ako netko samostalno postavi sve JSON zastavice na true, ovaj kod **ne može** ispisati REAL_DATA_READY, jer nema neovisnoga čitanja i provjere stvarnih providerovih file byteova/headera. Za to je potreban nov zaseban preregistrirani real-data stage.

[CI 36443992886](https://github.com/dvlahek/stress-energy-closure/actions/runs/36443992886) **SUCCESS** potvrđuje samo ispravan fail-closed offline gate. [Izvorni statusni JSON](../source_data/eboss_dr16_a03_e17d2b3b0_archived_CI_2026_09_28/e17d2b3b0_offline_real_abacus_access_contract_BLOCKED.json), SHA256 **3917ce7da42b9064bdaece35e7ed1ef63dc880bf8fb48a5143a38f3882d4c4fe**, navodi **20 konkretnih otvorenih uvjeta**. [Permanentni SHA manifest](../source_data/eboss_dr16_a03_e17d2b3b0_archived_CI_2026_09_28/archive_manifest.json), Git blob **b5b6827b0687b75ed4c0b2db5642bf1e150683a9**, arhiviran je na audit grani. Nije preuzeta ili otvorena nijedna stvarna Abacus datoteka.

## WSL naredba — samo kodni self-test bez znanstvenih podataka

WSL **nije potreban za završetak E17D2b3b0**, jer je GitHub CI već PASS. Za neobveznu lokalnu provjeru na postojećem računalu upotrijebiti sljedeće naredbe. Svaka Python provjera ispisuje PASS i BLOCKED, ne ostaje nijema.

    cd ~/stress-energy-closure
    git status --short
    git fetch origin audit/eboss-elg-bit8-ra-orientation-20260925
    git switch audit/eboss-elg-bit8-ra-orientation-20260925
    git pull --ff-only origin audit/eboss-elg-bit8-ra-orientation-20260925
    python -u scripts/audit_eboss_dr16_a03_e17d2b3b0_access_contract.py --self-test

Ako su prisutne lokalne izmjene, ne koristiti prisilni reset. Naredba dohvaća samo Git kod i zaključanu malu metadata arhivu iz audit grane; ne pokreće Globus, ne čita ASDF, ne pokreće CLASS i ne pristupa opaženim galaksijama.

## Sljedeći E17D2b3b1 — zaseban stvarni data gate

Prije ijednog pravog preuzimanja zasebno potvrditi providera, odabrani raw/cleaned DOI, točan c000 box/phase, verziju i file path, veličinu i dostupan slobodan prostor, dostupnost točnog halo_info superslaba i odgovarajućega cleaned/raw tree produkta. Nakon zasebnog odobrenja pristupa i transfera provjeriti na **stvarno dobivenim byteovima** GNU/POSIX cksum i lokalni SHA256 te ASDF header redshift, kozmologiju, ParticleMassHMsun, SODensityL1 i stvarnu mass density/species referencu. Povezati stvarni progenitor/descendant i čišćenje across epochs. **L1 N puta particle mass i SO_radius nisu M200c**; potrebna je nezavisno verificirana providerova definicija M200c ili adekvatan full-particle 200rho_crit radial reconstruction s halo+field okolnim česticama, centrima i periodic slab padding.

Čak i uz valjan halo mass history, ne-termalni E8 F+/F− neutrinski odziv nije generiran Abacusovim smooth-neutrino N-bodyjem; originalni As/tau mismatch, fizička LRG/ELG HOD i stvarna survey 3-point estimator/window/covariance ostaju odvojeni uvjeti. Originalni 24D 2-point odd nije bispektar. A03 pair-z RR, A04, SGC independent reverse i 48k ostaju otvoreni.

**Nema stvarnog halo M200c(a), neutrinova wakea, punoga Einstein–Vlasov finite-K galaktičkog B, eBOSS ξ/SNR ni značajnosti. Observed odd SEALED, originalni F i 72/24 frozen, main netaknut, draft PR #1 bez mergea.**
