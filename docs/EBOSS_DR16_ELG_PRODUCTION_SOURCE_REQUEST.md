# Request for the original DR16 ELG supplemental-mask production sources

**Status:** Draft request only. Not sent. Do not treat a reply or an offered file as authenticated until its release/version, raw bytes and source are recorded. This source request changes no selection, published random, veto or observed odd statistic.

**Suggested recipient:** Authors of Raichoor et al., MNRAS 500 (2021), 3254–3274, especially the DR16 ELG LSS catalogue and mask-production authors. The paper's Data Availability states that custom Python lines for bits 8–11 and the two `eboss22` low-quality plates are available on request.

**Subject:** Reproducibility request: original eBOSS DR16 ELG angular-mask Python and bad-plate geometry

Dear Dr Raichoor and colleagues,

I am checking the angular-selection provenance of the public eBOSS DR16 ELG clustering data and random catalogues for a cross-tracer analysis. Your MNRAS 500 (2021) paper cites `brickmask` v1.0 for mask bits 1–7 and notes that custom Python lines for bits 8–11 and the two low-quality `eboss22` plates are available on request.

Could you please share, or point me to, the version of those additional Python lines and their configuration that was actually used for the DR16 release? In particular, it would help to have the coordinate convention used for the 37-pixel HEALPix bit-8 correction; the ordering and boundary treatment of the three extra polygon masks (bits 9–11); and the original region/plate-overlap exclusion for PLATE-MJD 9430-58112 and 9395-58113.

If available, I would also appreciate the corresponding source/version or invocation notes for applying these operations to the released ELG data and randoms, or a versioned pre-veto mask reference against which the implementation could be checked. I am trying to reproduce the published selection without introducing any additional veto or changing the released catalogues.

Thank you for your help.

Best regards,

Dino Vlahek

## Evidence to attach if requested

- Paper-cited `brickmask` v1.0 is the annotated tag `3a94b72f4ee39c215cab713cf464d54106899838`, commit `4c0f940934ff8c4e6b0f0c709b3093383abff8d0`, with tag message “Version for eBOSS DR16”. It does not contain the later extra-mask Python helper.
- A later public `brickmask` checkout `b9eb684a579b56ec3dbdb46549224be7e3fa2830` contains `scripts/eBOSS_ELG_extra.py`. That source is useful for comparison, but its role in the actual 2020 DR16 release is not proven.
- The released full ELG FITS catalogue has already been SHA-checked locally, with a separately frozen `RA`/`DEC`/`mskbit`-only aggregate bit-8 correspondence audit. Do not send galaxy positions or other confidential research material.
- The two bad plate identifiers are published in §3.2 item (xii) of Raichoor et al.; their continuous production veto geometry is not established by identifiers alone.

**Do not send automatically.** Any subsequently received files should be saved as exact raw bytes and SHA-pinned before parsing or applying them, with clear separation between source-only inference and observed-odd analysis.
