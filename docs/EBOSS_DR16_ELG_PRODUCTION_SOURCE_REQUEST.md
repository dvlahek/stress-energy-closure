# Request for exact DR16 ELG supplemental-mask production sources

**Status (26 September 2026):** Prepared for the user's approval. **Not sent**, no recipient contacted, no file or catalogue modified. Any received material must be archived at its raw original byte identity and source attribution **before** inspecting or applying it. This draft does not authorize unblinding the observed odd-sector statistic.

**Suggested recipients:** Dr Anand Raichoor and Dr Arnaud de Mattia, as named coauthors of the DR16 ELG LSS catalogue study. The [2021 paper](https://academic.oup.com/mnras/article/500/3/3254/5942664) explicitly offers the supplemental Python on request; the [public `brickmask` README](https://github.com/cheng-zhao/brickmask/blob/b9eb684a579b56ec3dbdb46549224be7e3fa2830/README.md) credits Dr de Mattia for sharing the original extra-maskbit script. An acknowledged source contribution is **not** independent proof of exact historical execution.

**Subject:** eBOSS DR16 ELG LSS: request for original extra-mask Python and excluded-plate rule

Dear Dr Raichoor and Dr de Mattia,

I am checking the selection provenance of the published eBOSS DR16 ELG LSS data and random catalogues for a cross-tracer analysis. Your catalogue paper points to `brickmask` v1.0 for bits 1–7 and notes that the additional Python for bits 8–11 and two excluded `eboss22` plates is available on request.

Could you please share the original supplemental Python and configuration used for the DR16 release, or a versioned reference to that exact production implementation? In particular, I am looking for the coordinate convention and 37-pixel HEALPix rule for bit 8; the ordering, boundaries and overlap treatment for polygon bits 9–11; and the exact region/overlap exclusion for PLATE-MJD `9430-58112` and `9395-58113`.

It would also help to have the release-specific invocation or acceptance rules for both the ELG data and random catalogues, including the stage at which the excluded plates were applied. If the original scripts are no longer available, a versioned production mask reference or corresponding pre-veto output would be useful.

The goal is to reproduce the published selection without introducing new vetoes or changing the released catalogues. I will preserve the original source files and their checksums.

Thank you for your help.

Best regards,

Dino Vlahek

## Evidence to attach if requested

- Paper-cited `brickmask` v1.0 is the annotated tag `3a94b72f4ee39c215cab713cf464d54106899838`, commit `4c0f940934ff8c4e6b0f0c709b3093383abff8d0`, with tag message “Version for eBOSS DR16”. It does not contain the later extra-mask Python helper.
- A later public `brickmask` checkout `b9eb684a579b56ec3dbdb46549224be7e3fa2830` contains `scripts/eBOSS_ELG_extra.py`. Its README (Git blob `98514607195e01f035e4898f3a647af402077e53`) explicitly credits Arnaud de Mattia for sharing the original extra-maskbit script. This makes a direct source request well targeted but does **not** authenticate that later file as the actual 2020 DR16 production executable.
- The released full ELG FITS catalogue has already been SHA-checked locally, with a separately frozen `RA`/`DEC`/`mskbit`-only aggregate bit-8 correspondence audit. Do not send galaxy positions or other confidential research material.
- The two bad plate identifiers are published in §3.2 item (xii) of Raichoor et al.; their continuous production veto geometry is not established by identifiers alone.

**Do not send automatically.** Any subsequently received files should be saved as exact raw bytes and SHA-pinned before parsing or applying them, with clear separation between source-only inference and observed-odd analysis.

## Acceptance rules for any response

A reply, public helper link or reconstructed notebook is **not automatically a certified historical production source**. Record the sender/URL, release and version, original raw bytes, full-file SHA256 and exact invocation/configuration first. A-03 retains two explicitly different conclusions: authenticated historical mask production (only with original versioned execution/acceptance evidence) or a strictly limited empirical pair-window from existing published tracer-specific randoms (requiring independent blind validations). Neither outcome is an automatic A-05 observed-odd unblinding gate.

The corresponding frozen decision document is `source_data/eboss_dr16_a03_selection_window_provenance_decision_2026-09-26.json`.
