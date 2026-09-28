# A-03 E17D2b3b4 — USER-SUPPLIED Globus raw halo_info listing and exact checksum manifest

**28. 9. 2026. · Source-only metadata gate.** Observed eBOSS 24D odd **SEALED**. No real ASDF bytes, observed galaxies, CLASS, particles, merger trees or new mocks were obtained or read. Only the existing public Git-pinned Abacus portal manifest plus the user's pasted Globus listing and uploaded 1,386-byte `checksums.crc32` were used.

## Problem and prior source

The prior E17D2b3b3/R1 source-pinned portal table gave 35 **directory entries**, aggregate 75,722,740,486 bytes, in `AbacusSummit_base_c000_ph000/halos/z0.950/halo_info/`. It did **not** establish 35 ASDFs or any individual filenames/checksums. Its public source is `abacusorg/abacussummit-data-portal` commit `f877298344a6988d55f57a7c39f5be14dce63ff3`, table blob `cb235ad72ab2516f872bd4ca79b9447a3f4b30e3`; earlier correction and original reports stay immutable.

## User intake, provenance and prospective gate

The user manually pasted an authenticated-looking Globus File Manager listing under **AbacusSummit Public Data Access** showing `checksums.crc32` and `halo_info_000.asdf` through `halo_info_033.asdf`. Its visible sizes were rounded. The user independently uploaded `checksums.crc32`, 1,386 exact bytes and 34 GNU/POSIX `cksum` records. This is **USER-SUPPLIED evidence**, not a connector-authenticated provider session, signed provider manifest or independent origin attestation.

Before the new aggregate audit, [the E17D2b3b4 prospective prereg](../source_data/eboss_dr16_a03_e17d2b3b4_user_globus_checksum_intake_prereg_2026-09-28.json) was committed as Git blob `a6ad151ff116d8b20208bebe76d4423db589b76f`. It froze the original portal aggregate, the user's input hash, exact selected `halo_info_000.asdf` CRC/size, dual-parser audit and fail-closed controls. The byte-identical archived [user-provided checksum manifest](../source_data/eboss_dr16_a03_e17d2b3b4_user_supplied_checksums.crc32) has Git blob `0bca2fd0a28c9e163ad6198f67616da4355f4c81` and SHA256 `b1db66f78cf17d9c2620a0be0917547eb5dc17a9e8f10bd02b0b211691fcb6a5`. Its source is not upgraded to provider-attested merely because it is committed.

## Exact audit result

Independent regex and non-regex tokenization of all 34 records agree, including every filename, CRC and byte count. Filenames are exactly `halo_info_000.asdf` through `halo_info_033.asdf` with no gaps or duplicates. Summed exact record sizes give **75,722,739,100 bytes**. Adding the exact **1,386 bytes** of `checksums.crc32` gives **75,722,740,486 bytes**, identical to the prior source-pinned public portal snapshot directory aggregate. The 35 prior directory entries are explained by the user's 34 ASDF records plus one checksum manifest; this does not independently authenticate the user-uploaded bytes.

**Preregistered one-file choice:** `halo_info_000.asdf`, manifest POSIX CRC `3502514872`, exact recorded size `2,223,483,833` bytes. The file itself has **NOT** been downloaded or read. No per-file SHA256 of real ASDF exists yet. No physical individual halo mass, `M200c(a)`, halo tree, neutrino F+/F− wake or physical galaxy bispectrum follows from this source-only audit.

[Original report](../source_data/eboss_dr16_a03_e17d2b3b4_archived_source_only_2026_09_28/e17d2b3b4_original_user_checksum_ledger.json) SHA256 `086da6da6958af67bfbb02ff59fe6b9dc4a81b71497b76671644f4bf4d883de8`, [independent report](../source_data/eboss_dr16_a03_e17d2b3b4_archived_source_only_2026_09_28/e17d2b3b4_independent_user_checksum_ledger.json) SHA256 `0aac5757798a4aaf3c6351f72ef41e68c093ea297dfd74b79ccfec296ff21d67`, and [archive SHA manifest](../source_data/eboss_dr16_a03_e17d2b3b4_archived_source_only_2026_09_28/archive_manifest.json) are all retained on the audit branch. [Pure-stdlib replay and six synthetic negative controls](../scripts/audit_eboss_dr16_a03_e17d2b3b4_user_Globus_checksum_manifest.py) do not read ASDF or connect to Globus. The exact GitHub Actions workflow is `.github/workflows/eboss_a03_e17d2b3b4_user_globus_checksum_source_only.yml`; its run conclusion must be checked separately.

## Next real-data gate (NOT executed)

Register a separate exact one-file transfer/access protocol only after checking that the authenticated Globus collection ID and directory are the previously pinned official source. The only planned real candidate is `halo_info_000.asdf` and its same-directory original checksum manifest, **not** the 34-file folder, particles, a 6.6 TB base tar or the cleaned/tree products. Preserve the original local checksum manifest. On any later separately authorized local transfer, compute the full file SHA256 and GNU/POSIX `cksum`, compare exact byte count and CRC, then pass the file through existing E17D2b3b2R1 fail-closed local ASDF header verifier with the original pinned Abacus decoder and reference source-copy header. A user-provided checksum manifest alone never proves independent provider origin. Even a passing header check does not establish `M200c` or `M200c(a)`. Those require separately identified real per-halo profile and same-object time/merger-tree data.

**STOP:** observed 24D odd SEALED; original E8 F0/F+/F− 4000q, E16 576 triangles, 72 source /24 contrast, 48k unchanged. No physical nonlinear wake/B, new science cut/seed/download or significance. `main` unchanged and PR #1 draft.
