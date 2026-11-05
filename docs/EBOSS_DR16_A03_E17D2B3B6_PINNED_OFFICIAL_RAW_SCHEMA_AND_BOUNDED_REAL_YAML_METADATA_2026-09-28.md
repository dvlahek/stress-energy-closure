# E17D2b3b6 — službena raw halo shema i uski YAML metadata gate

**28. 9. 2026. · Source-only provjera PASS; korisnikov stvarni lokalni metadata probe čeka izvršenje.**

## Problem i rezultat za E17D2b3b5

Korisnikov stvarni `AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/halo_info_000.asdf` prošao je lokalnu provjeru byte-counta, GNU/POSIX `cksum`, punog SHA256 i dvaju čitača stvarnog zaglavlja. Dostavljeni SHA256 realne datoteke je `7f3755a469b3c25c869a0394b06192780d4a03717fd918d01c7459b0405724c3`. [Arhiviran je korisnikov izvještaj](../source_data/eboss_dr16_a03_e17d2b3b5_user_reported_real_halo_000_byte_and_dual_header_local_acceptance_2026-09-28.json), Git blob `c0ac072263f3d8b19d01caae091cb743c22a4037`. Ovo je **USER-REPORTED LOCAL** rezultat, nije provjera stvarne ASDF datoteke na GitHubu ili neovisna providerska atestacija. Nijedan stvarni halo niz nije pročitan.

## Prospektivna odluka prije prvog stvarnog halo polja

[Zaključani b6 protokol](../source_data/eboss_dr16_a03_e17d2b3b6_official_schema_real_ASDF_YAML_descriptors_only_prereg_2026-09-28.json), Git blob `e1f43da1a73aad3ef85c76a114fbffadc76d4761`, dopušta **samo tekstualni YAML prefiks do prvog samostalnoga `...`**, maksimalno **262 144 bajta**. Točna datoteka ostaje 000 i istodirektorijski checksum manifest sa SHA256 `b1db66f78cf17d9c2620a0be0917547eb5dc17a9e8f10bd02b0b211691fcb6a5`. Dopušteni su samo nazivi stupaca i opisni metapodaci oblika, tipa i bloka. Zabranjen je pristup bilo kojoj vrijednosti realnog halo polja, `np.asarray`, rezanje polja, Blosc dekompresija ili novo preuzimanje.

[Source-only rezultat i CI reference](../source_data/eboss_dr16_a03_e17d2b3b6_source_only_official_schema_and_synthetic_QA_result_2026-09-28.json), Git blob `2684cef727a9cf03752707adf89d22b7e1c5005c`.

## Što zapravo kaže službeni izvor

Na upstream commitovima [AbacusSummit `4b1959c...`](https://github.com/abacusorg/AbacusSummit/tree/4b1959c710cb0c49aa305c6213a228aa2a4587ff) i [abacusutils `24ab0dda...`](https://github.com/abacusorg/abacusutils/tree/24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6), Git blobovi data-products `f7c847b174365744b23cb6d1ebaeb010a5bd6ca7`, CompaSO `a5170574467501e0b9adb4a71e8e432a32843dbf`, loader `adb16aee1cbac863db5301c6e380938ee7f76547` i official Blosc extension `8c5e5135736409bb0fab77bff06ad1259951189d` izvorno potvrđuju:

- `z0.950` sekundarni output: katalog haloa i poduzorkovani particle PID podaci; nije službeno obećan odgovarajući halo RV/field položajni proizvod za tu epohu. Iz drugog primarnog `z` ne smijemo imputirati istovremeni realni profil.
- Raw `halo_info` ima **zaseban ASDF komprimirani binarni blok po stupcu**, stoga uski pregled tekstualne sheme ne mora dekomprimirati binarni sadržaj.
- `N` je broj čestica dodijeljenih CompaSO L1 halou. Kataloška masa `N * ParticleMassHMsun` jest L1 assigned mass, ne puna `M_{200c}`. CompaSO L1/SO nalazi ovise o česticama L0 grupe, mogu imati kompetitivnu podjelu čestica.
- `SO_radius` označava radijus prema CompaSO SO postupku, no kod nedosegnutog SO crossing postoji **posebna konstantna granična vrijednost**. L1 SO prag `SODensityL1` izražen je u srednjoj kozmičkoj gustoći, ne fiksno `200 rho_crit`.
- Raw nenormalizirani položaji/radijusi izraženi su u jedinicama kutije. Pomoćni abacusutils loader množi ih s `BoxSize` da dobije komovne Mpc/h. Centar `SO_central_particle` ne smije se bez dodatne odluke zamijeniti `x_L2com`.
- `abacusutils` ima `cleaned=True` kao zadano; time može zahtijevati cleaned proizvode koji nisu lokalno preneseni. b6 ne poziva taj loader.

**Posljedica:** jedna raw sekundarna halo_info datoteka sama ne daje potpun individualni sferni `200rho_crit` profil. Budući eventualni NFW, percentilni ili L1-only proxy mora biti jasno deklariran kao modelna/uvjetna procjena, nikad izmjereni pravi `M_{200c}`. Za puni profil treba zasebno odobreni source/data gate s odgovarajućim česticama i eksplicitno određenim centrom, granicom i epoškom gustoćom.

## CI i naredba za lokalnu metapodatkovnu provjeru

[Source-only GitHub CI `36475043792` SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36475043792) provjerio je sve četiri upstream Git blob reference, 12 izvornih semantičkih ograničenja, sintetički YAML fixture i šest negativnih kontrola. **Nijedna stvarna ASDF datoteka nije ni preuzeta ni otvorena u CI-u.**

Opća repository integrity provjera pri prvom b6 commitu pogrešno je tretirala dva upstream official docs path stringa kao lokalne datoteke. [Usko pravilo na audit grani](../scripts/check_repository_integrity.py) sada izuzima samo ta dva **točno Git-pinnana vanjska** putna zapisa uz obveznu provjeru upstream repo/commita; prvotni b6 prereg ostaje neizmijenjen. [CI `36475172756` SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36475172756).

Za korisnikov već postojeći lokalni ASDF predviđen je [b6 bounded YAML reader](../scripts/audit_eboss_dr16_a03_e17d2b3b6_bounded_real_ASDF_YAML_schema_metadata.py), Git blob `718d2dde75f29dfbf2e32ab54cb624fc2a6688f8`. Naredba se izdaje korisniku nakon git fast-forward audita; rezultat ostaje samo mali lokalni JSON. `id`, `N`, `SO_radius` i `SO_central_particle` obvezni su nazivi polja, dodatni nazivi mogu biti samo metadata-only; ostali realni halo redci ili binarni blokovi ne otvaraju se.

**STATUS:** b6 SOURCE-ONLY/SYNTHETIC PASS, b6 ACTUAL LOCAL YAML FIELD-DESCRIPTOR PROBE PENDING. Ne tvrditi da je već poznata stvarna lokalna shema dok je korisnik ne izvrši. E8 4000q, originalni F0/F+/F−, E16 576/72/24, 48k zamrznuti, observed odd 24D SEALED. Bez S/N, promjene main ili merge PR-a.
