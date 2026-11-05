# A-03E17D2b3b1 — službena kopija zaglavlja identificira z0.950 kao z=0.952838237036305

**28. 9. 2026. · [CI 36445799032 SUCCESS](https://github.com/dvlahek/stress-energy-closure/actions/runs/36445799032) za točno SHA-pinnanu službenu KOPIJU SIMULACIJSKIH ZAGLAVLJA. Nijedan stvarni halo_info ASDF, katalog, halo particle ili tree nije otvoren. Observed odd SEALED. Puni fizički bispektar BLOCKED.**

## Problem i izvor

E17D2b3b0 nije mogao prihvatiti naziv sekundarnog direktorija z0.950 kao stvarni redshift: [Abacus dokumentacija](https://abacussummit.readthedocs.io/en/latest/data-products.html) objašnjava da sekundarni izlazi ne moraju pogoditi nominalnu epohu. [Službeni abacusutils metadata modul](https://abacusutils.readthedocs.io/en/latest/metadata.html) sadrži komprimiranu arhivu kopija izvornih simulacijskih zaglavlja. Ta javna, unaprijed proizvedena arhiva omogućuje prvi source-only test bez preuzimanja ogromnih halo kataloga.

[Prospektivni protokol](../source_data/eboss_dr16_a03_e17d2b3b1_official_embedded_header_metadata_prereg_2026-09-28.json) zaključava originalne E8/E16/E17A/E17D2b2/E17D2b3a/E17D2b3b0 roditelje i upstream [abacusorg/abacusutils commit 24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6](https://github.com/abacusorg/abacusutils/tree/24ab0dda5fea9ae406b1afdacaf4bbf989de9bc6). Točan upstream Git blob arhive je 0065029cdd515ffa76c61a27b81e9f16b2656f38, veličina 9 486 279 B i izmjereni SHA256 1c8b690acedd2a61a6a033be3173f79bb7e8cadbcb4a9e97e5db10ee4f10308e. Službeni kompresijski dekoder koristi se s istoga zaključanog upstream commita.

Prva verzija prerega imala je blob a282b157a80f88a1e0bf95c01942453d9a59efe4. CI preflight 36445390661 i integrity provjera 36445325404 otkrili su tehnički nedostatak: link na upstream skriptu bio je naveden kao da je lokalni fajl i SHA već zamrznutoga E17D2b2 roditelja nije bio upisan u parent map. **Prije ikakvog čitanja upstream metapodataka** dopunjen je prereg u bloba 36c2bcda723b95565e604723e651beef08aad992, commit 4ec8c2a, s točnim zapisom o amandmanu. Nisu mijenjani simulacija, z label, parametri ni kriteriji. CI 36445592859 i 36445617020 zatim su prošli SHA-pinnani download kopije, ali su odbili čitanje jer osnovni ASDF nije registrirao službeni Abacus Blosc tip kompresije. Završni workflow je to tehnički popravio registracijom SHA-pinnanoga službenoga extensiona, bez promjene izvornoga podatka.

## Stvarni rezultat u službenoj KOPIJI zaglavlja

| Izvorno kopirano polje | Vrijednost |
|---|---:|
| Simulacija | AbacusSummit_base_c000_ph000 |
| Točan metadata state key | z0.950 |
| Redshift iz kopije zaglavlja | **0.952838237036305** |
| ScaleFactor iz istog statea | **0.5120751842291016** |
| Redshift izveden iz 1/a − 1 | 0.9528382370363053 |
| OmegaNow_m iz istoga statea | 0.774150059509772 |

Izvorni header-copy redshift i vrijednost iz faktora skale slažu se unutar ranije zaključane kontrole. Naziv direktorija z0.950 razlikuje se od stvarne zapisane epohe za 0.002838237036305. Nije dopušteno zbog toga mijenjati izvorne matematičke E17A z čvorove (.945, .95, .955) ili nove science cuts, ili interpretirati z kao već potvrđen redshift konkretnoga halo ASDF-a. Izvorni E17D2b3a kontraprimjer i njegove ilustrativne kozmološke vrijednosti ostaju zaseban, zamrznuti matematički test, ne stvarna konverzija Abacus L1 u M200c.

[Izvršivi SHA-pinnani runner](../scripts/audit_eboss_dr16_a03_e17d2b3b1_official_header_copy.py) dekodira točno jedan source-copy state, provjerava Redshift protiv ScaleFactor i dodatno ponavlja dekodiranje **istih byteova**; nije neovisna simulacija niti čitanje stvarnog haloa. [Originalni rezultat](../source_data/eboss_dr16_a03_e17d2b3b1_archived_CI_2026_09_28/e17d2b3b1_pinned_official_embedded_header_copy_state_z095.json), SHA256 aacd9b4f6e3dd63fb353c72488669e7eb01e365d48c220664a0dfda8dc9e0673, i [trajan SHA manifest](../source_data/eboss_dr16_a03_e17d2b3b1_archived_CI_2026_09_28/archive_manifest.json), Git blob 4eb928e90a98dd35add7d4b8ebe7aedbcf14e0dc, spremljeni su na audit grani.

## Što je i dalje BLOCKED

Službena kopija sadržava nazive polja ParticleMassHMsun i SODensityL1, ali **ne potvrđuje sadržaj, stvarni z i cksum odabrane konkretne halo_info datoteke**. Za fizički E17D2b3b2 potreban je novi zasebni real-data pristupni protokol, stvarno verificiran c000 box/phase/raw-or-cleaned/product path, službeni GNU/POSIX cksum plus lokalni SHA256 odabranih stvarnih file byteova, stvarni ASDF header z, ParticleMassHMsun, gustoćnu konvenciju i dosljedni merger-tree branch. KompaSO L1 N/SO_radius nije M200c; potreban je providerov validirani 200-kritični maseni proizvod ili dostatan halo+field pozicijski obuhvat za remeasurement. Izvorni CLASS As/tau, originalni F+/F− custom ncdm i Abacus smooth-neutrino N-body ostaju fizikalno nejednaki za originalni neutrinski wake, uz nepoznati stvarni high-z LRG/ELG HOD i realan eBOSS 3pt window/covariance.

**Nije otvoren observed odd, nije preuzet stvarni ASDF halo/tree, nema M200c(a), fizičkoga Einstein–Vlasov long–short galaktičkog B, eBOSS ξ/SNR ili detekcijske značajnosti. Originalni 72/24 i F± frozen, A03 pair-z RR/A04/SGC reverse/48k otvoreni, main netaknut i PR #1 draft bez mergea.**
