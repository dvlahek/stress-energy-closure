# E17D2b3b2 — stvarni L1 prag i masa čestice iz službene KOPIJE Abacus zaglavlja, bez stvarne halo datoteke

**28. 9. 2026. · CI 36447312957 PASS. Podrijetlo: SHA-pinnana službena source-bundled kopija originalnoga simulacijskog zaglavlja. Per-file halo_info ASDF, provider checksum, pojedinačni M200c, merger tree, neutrinski wake i galaktički B i dalje BLOCKED. Observed odd SEALED.**

## Problem

E17D2b3a je pokazao da jednaka idealizirana L1 masa/radijus ne određuje jednoznačno M200c bez unutarnjega radijalnog profila, ali je rabio **ilustrativan** L1 prag 200 puta mean density u jednostavnom flat matter+Λ toy modelu. Taj se ilustrativni broj nikad nije smio prenositi u realni Abacus. E17D2b3b1 zatim je iz SHA-pinnane službene header-copy arhive utvrdio da AbacusSummit_base_c000_ph000 ima actual *source-copy state* redshift **0.952838237036305**, ne nominalnu oznaku z0.950. Za sljedeći pristup pravim halo podacima trebamo istu vremensku točku i **stvarni L1 SO threshold i particle mass iz toga kopiranoga statea**.

U [službenom opisu data products](https://abacussummit.readthedocs.io/en/latest/data-products.html) stoji da je CompaSO L1 definiran epohno promjenjivim SODensityL1 **u odnosu na srednju kozmičku gustoću**, da je primarna masa L1 haloa broj članova N puta ParticleMassHMsun, te da SO_radius može biti specijalna konstantna vrijednost ako crossing nije dostignut. M200c zahtijeva 200 puta **kritičnu** gustoću i zato se L1 masa/radijus ne smiju preimenovati u M200c.

## Prospektivno zaključan ulaz i izvršenje

[Protokol E17D2b3b2](../source_data/eboss_dr16_a03_e17d2b3b2_official_header_copy_mass_reference_prereg_2026-09-28.json), Git blob **324e722200997b21a0fa97d965b8ef096e1587e6**, commit **80dcdf9**, zaključan je prije ponovnoga čitanja upstream arhive. SHA-pina zamrznute originalne E8 CSV, E16 geometriju, E17A original CLASS joint, prethodni E17D2b3b1 originalni JSON/manifest/protokol/runner i prethodni E17D2b3b0 BLOCKED report. Unaprijed zahtijeva šest konkretnih numeričkih polja, provjeru z-a preko faktora skale i dvije negativne kontrole. Nije pretpostavio nepoznate vrijednosti L1 praga ili mase čestice.

Izvorno [abacusorg/abacusutils commit 24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6](https://github.com/abacusorg/abacusutils/tree/24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6) daje samo ugrađenu komprimiranu KOPIJU simulacijskih zaglavlja, 9 486 279 bajtova, Git blob **0065029cdd515ffa76c61a27b81e9f16b2656f38**, SHA256 **1c8b690acedd2a61a6a033be3173f79bb7e8cadbcb4a9e97e5db10ee4f10308e** i izvorni službeni Blosc dekoder iste verzije. Nije preuzet niti pročitan stvarni halo_info superslab ili merger tree. [CI 36447312957](https://github.com/dvlahek/stress-energy-closure/actions/runs/36447312957) je **SUCCESS**.

## Izmjereni source-copy metapodaci

| Polje za AbacusSummit_base_c000_ph000, metadata key z0.950 | Vrijednost |
|---|---:|
| Redshift | **0.952838237036305** |
| ScaleFactor | 0.5120751842291016 |
| OmegaNow_m | **0.774150059509772** |
| SODensityL1 u jedinicama srednje kozmičke gustoće | **228.52306365966797** |
| ParticleMassHMsun | **2109081520.453063 Msol/h** |
| ParticleMassMsun | 3131059264.330557 Msol |

U ovoj kopiranoj pozadini

\[
\frac{200\rho_{\mathrm{crit}}}{\rho_m} =
\frac{200}{\Omega_m(z)} =
258.3478455412757.
\]

Stoga je **omjer L1 i 200-kritičnog praga gustoće**, ne omjer halo masa:

\[
\boxed{
\frac{\Delta_{\mathrm{L1,mean}}}
{\Delta_{200c,\mathrm{mean}}}
=\frac{228.52306365966797\,\Omega_m(z)}{200}
=0.8845557166574369.
}
\]

To je izravna usporedba *dvaju različitih definicijskih pragova* pri istoj source-copy epohi. Iz omjera nije moguće izvesti realni M200c/M_L1 bez individualnoga profila gustoće. Napose, broj 0.88456 nije očekivani halo-mass correction factor, galaxy bias, LRG/ELG signal ili predikcija eBOSS statističke značajnosti.

[SHA-pinnani izvorni runner](../scripts/audit_eboss_dr16_a03_e17d2b3b2_header_copy_l1_mass_reference.py) dvaput dekodira **iste** službene source-copy byteove, ponovno zatvara redshift/scale relation, izračunava prag neovisnim Decimal 40-digit putem i odbija sintetičko krivotvorenje SODensityL1=0 i redshifta prepisanoga iz nazivnog z0.950. Maksimalni prikazani relativni gap neovisne decimalne provjere omjera pragova je 1.1102230246251565e−16. Drugo dekodiranje istih byteova nije zasebna N-body realizacija niti independent sample of physical halos.

[Trajni izvorni JSON](../source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28/e17d2b3b2_pinned_header_copy_l1_threshold_vs_200crit.json) SHA256 **efa1878250ca395ab65207f485f214854b817188613acaae1ce668f58af78a38** i [trajni SHA manifest](../source_data/eboss_dr16_a03_e17d2b3b2_archived_CI_2026_09_28/archive_manifest.json) Git blob **0563c53fdb966bb9382a3b910e94e237fef90bf8** pohranjeni su u audit granu nakon završnog PASS-a.

## Stvarna halo datoteka i dalje je zaseban acceptance gate

Javni providerov [NERSC/Globus pristup](https://abacussummit.readthedocs.io/en/latest/data-access.html) omogućava uske skupove, ali konkretan c000 ph000 z0.950 halo_info superslab, pripadajući provider checksums.crc32 i tree datoteka nisu potvrđeni ili preuzeti u ovom koraku. Javni web portal nije dao provjerljiv file-item URL. Nema osnove izmisliti halo_info_000.asdf kao *stvarno dostupnu datoteku*. Za realni file gate treba potvrditi originalni provider item, sirove byteove i stvarnu veličinu, provjeriti GNU/POSIX cksum i lokalni SHA256 **na istoj datoteci** i pročitati njen ASDF header. Tada usporediti Redshift, ScaleFactor, OmegaNow_m, SODensityL1 i ParticleMassHMsun s gore zaključanim source-copy vrijednostima. Ako se ne poklapaju, FAIL-CLOSED. Ako se poklapaju, dokazana je samo kompatibilnost headera te datoteke s official-copy referencom, ne M200c ili halo merger history.

Za stvarni M200c kroz epohe treba nezavisno verificirani provider M200c ili puna relevantna halo+field pozicijska čestična okolina, periodičnost, centar i 200critical proračun. Sekundarni halo PID-only subsample nije dovoljan za točan full-particle M200c remeasurement. I dalje su otvoreni originalni CLASS As/tau nasuprot c000, originalni custom F+/F− neutrinski Vlasov wake nasuprot Abacus smooth-neutrino N-body, high-z LRG/ELG HOD/selection, stvarni eBOSS 3pt window/covariance, A03 pair-z RR, A04, SGC independent reverse i 48k. Originalni 24D dvotočkasti odd nije bispektar.

**Nema stvarne halo ASDF datoteke, individualnih masa, physical finite-K galaxy B, observed ξ/SNR ili detekcije. Originalni F i 72/24 frozen, observed odd SEALED, main netaknut, PR #1 draft bez mergea.**


## E17D2b3b2R1 — pripremljen stvarni-file provjerivač, testiran samo na sintetičkim bajtovima

[Prospektivni R1 protokol](../source_data/eboss_dr16_a03_e17d2b3b2r1_real_ASDF_byte_header_verifier_prereg_2026-09-28.json), Git blob **b5e4bed15c9841bdfa8933a4d6960c9eda701e36**, commit f3dfa3b, SHA-pina ovaj izvorni E17D2b3b2 report i manifest, originalni E17D2b3b1, E8/E16 i službeni upstream Blosc Python izvor. Prije implementacije definira obvezni path/polje/checksum gate te šest negativnih kontrola. [Pripremljeni lokalni verifier](../scripts/audit_eboss_dr16_a03_e17d2b3b2r1_local_halo_ASDF_checksum_header.py) **ne sadrži downloader** i ne treba pokretati dok stvarna providerova datoteka, pravi checksum manifest i SHApinnani dekoder ne postoje. GitHub [synthetic-only CI 36447886533](https://github.com/dvlahek/stress-energy-closure/actions/runs/36447886533) **SUCCESS** testirao je samo malu privremenu datoteku koju je sam stvorio. Nije kopirao stvarne Abacus byteove u Actions ili repozitorij.

Budući korisnički odabir stvarnoga superslaba mora biti unutar točne putanje AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/halo_info_NNN.asdf i uz checksums.crc32 iz te konkretne mape. Verifier odbija simboličke poveznice, pogrešnu simulaciju i putanju, nedostajuće ili duplicirane checksum zapise, nedosljedan GNU/POSIX cksum ili broj bajtova, a uvijek izračunava lokalni puni SHA256. Uz to koristi provjereni službeni Blosc extension i čita **samo ASDF header**, ne masene/galaktičke halo stupce. Za Redshift, ScaleFactor, OmegaNow_m, SODensityL1, ParticleMassHMsun i ParticleMassMsun zahtijeva podudaranje s prethodnim SHA-pinnanim source-copy brojevima. Samo ako se sve provjere slažu, može izdati **LOCAL_BYTES_MATCH_SUPPLIED_CKSUM_AND_OFFICIAL_COPIED_HEADER_PROVIDER_MANIFEST_ORIGIN_UNATTESTED**. To nije potvrda da je datoteka stvarno došla od providera: checksum manifest unesen lokalno bez neovisno provjerenoga providerova podrijetla može biti krivotvoren. Rezultat se čuva isključivo u **lokalnom** outputu, nikad automatski na GitHubu.

Dopušteni budući terminalni poziv, samo nakon zasebno odobrenog pribavljanja točne stvarne datoteke i njezina autentificiranoga checksum manifesta, ima oblik:

    python -u scripts/audit_eboss_dr16_a03_e17d2b3b2r1_local_halo_ASDF_checksum_header.py --file "<stvarna putanja>/AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/halo_info_NNN.asdf" --checksum-manifest "<ista halo_info mapa>/checksums.crc32" --official-extension "<upstream Git-pinnani abacusutils>/abacusnbody/data/asdf.py" --output "<lokalni izlaz>/halo_header_verification.json"

Gornji NNN i putanje **namjerno nisu izmišljeni** kao stvarno dostupne provider datoteke; CI ne koristi način rada --file. Prije prvoga stvarnog poziva treba provjeriti što provider zaista ima i zasebno odobriti taj pristup. Uspjeh ove provjere još uvijek ne dopušta M200c iz L1 N/SO_radius ili full galaxy B. Observed odd ostaje SEALED.
