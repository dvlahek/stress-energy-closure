# DR16 ELG production-lineage evidence and remaining selection gates

**Audit status, 2026-09-26:** public helper provenance and earlier released bit-8 labels are verified within their stated scopes. Historical release **execution**, the non-bitcoded `eboss22` bad-plate exclusion, continuous angular selection and the LRG×ELG physical pair window are **not certified**. This file does not authorize unblinding or any new selection.

## What is directly established

The byte-identical public `cheng-zhao/brickmask` extra helper at commit `b9eb684a579b56ec3dbdb46549224be7e3fa2830` has Git blob `d66c5e6ed6a8260b5e8dcc840795485345392ced`. The published helper itself, not a rewritten interpretation, assigns a bit 8 from its list of 37 HEALPix pixels. It adds bits 9, 10, 11 to `ELG_centerpost.ply`, `ELG_TDSSFES_62arcsec.pix.snap.balk.ply`, and `ebosselg_badphot.26Aug2019.ply`, respectively. All three polygon **source bytes** are separately SHA256-pinned in `source_data/eboss_dr16_eleven_official_polygon_sha_2026-09-25.json`. The publicly pinned C source under `-DEBOSS` uses bit 0 for rounded-pixel validity and bit 2 for the rounded-versus-truncated XYBUG correction. The four fixed example FITS images passed source-byte, WCS and CRC checks and have independently archived **aggregate** bit histograms. None of these facts proves how every released random position was processed.

The separately preregistered, SHA-gated audit of the official `eBOSS_ELG_full_ALLdata-vDR16.fits` decoded only `RA`, `DEC`, `mskbit` and found 15 bit-8-positive labels among 269,178 rows. At those same positions, the fixed native-RA candidates reproduced all 15 positives and all 269,163 negatives; the original literal call did not reproduce the positives. Exact report SHA256 `950fc1edb1c84d5d57545bb9b7617e632859f31953f9a9102547d5832d188751`, archived in `source_data/eboss_dr16_elg_full_bit8_label_geometry_audit_2026-09-25.json`. This authenticates a **released catalogue label correspondence at observed positions**, not the actual historical executable or complete inter-object angular mask. The original upstream helper visibly passes radians with `lonlat=True` and also reflects RA. Do not edit or substitute the upstream call based on the observed 15 labels.

A separate fixed-stratum pilot queried three extra polygons on 80,000 **already released** observed/mock0001 clustering random rows in both caps, with zero sampled memberships. This cannot establish the full-catalogue membership, and did not include full BRICKMASK pixels or bit 8. Its evidence is `source_data/eboss_dr16_elg_extra_paired_random_pilots_2026-09-25.json`.

## Independent published release contract

Raichoor et al., *MNRAS* **500** (2021), 3254–3274, §3.1–3.2, Table 5, and Data Availability:
https://academic.oup.com/mnras/article/500/3/3254/5942664 .
The paper distinguishes early target-selection bits 1–5 from later LSS angular vetoes, states that vetoes are stored in catalogue `mskbit` **except** the separately removed pair of low-quality `eboss22` plates, and specifies bits 1–7 via public `brickmask`. Its availability statement says bits 8–11 and two low-quality plates require additional Python lines available on request. That publication establishes an **algorithm family and exceptions**, not that the version of public helper audited here was the one executed by the original DR16 release. Paper Table 5 target-removal counts and survey areas must not be identified with per-brick pixel histograms or unconditionally summed when veto regions overlap.

The published ELG randoms sample a tracer-specific angular selection and receive shuffled data redshifts conditionally by chunk and imaging depth. See Raichoor et al. §3.1 and the public eBOSS ELG anisotropic clustering description:
https://academic.oup.com/mnras/article/499/4/5527/5917997 .
LRG randoms have a distinct MANGLE/completeness selection. The LRG×ELG random cross-pairs `R_L R_E`, with independently normalized cross-count terms, sample the published pair-window. Constructing a shared thresholded HEALPix hard mask is **not** equivalent. No newly reconstructed veto should be blindly applied to a previously selected published clustering catalogue.

## The actual provenance blockers

| Gate | Current evidence | Missing item / acceptance criterion |
|---|---|---|
| ELG BRICKMASK bits 1–7 | Pinned C `brickmask` code and four audited image samples; full four-family filename list of 19,381 entries | Authenticated version/configuration, source SHA inventory as needed and original published *data and random* acceptance invocation. Four images do not certify the remaining 19,377. |
| Extra bit 8 | Pinned 37-pixel public script and 15/15 released catalogue label match for native RA at exact catalogue positions | Authenticated release production Python/version or independently verified pre-veto output containing positive **and negative** discriminating bit-8 positions and an explicit RA/DEC convention. Do not silently 'repair' public helper source. |
| Extra polygon bits 9–11 | Exact public source mapping and three SHA-pinned polygon files; fixed random-only 80k pilot | Original release acceptance logic, boundary/overlap handling and exhaustive or separately registered released-random checks, without retroactive veto. |
| `eboss22` low-quality plates | Published **non-bitcoded** exclusion with exact PLATE-MJD identifiers `9430-58112` and `9395-58113`; published 13.9 deg² / 3,123 targets | Original production *region geometry and overlap/priority rule*, with release-versioned application to data/randoms. Identifiers alone are not a plate-footprint polygon. |
| ELG radial/selection randoms | Released tracer-specific catalogues and already passed random-only fine weighted `n(z)` checks | Historical angular/radial construction and imaging-depth-bin assignment provenance, or a clearly limited released-random empirical window validation with explicit scope. |
| Joint cross-estimator | Published random-only and nine-mock window pilots already archived | Mock-galaxy cross-estimator closure and sufficient joint 18D covariance under the **released** tracer-specific randoms; no observed odd vector. |

The exact known source/aggregate relationships and deliberately open statuses are machine-readable in `source_data/eboss_dr16_elg_production_lineage_gate_2026-09-26.json`. An offline AST/SHA gate, `scripts/audit_eboss_dr16_elg_public_source_lineage.py`, verifies the public helper and archived early evidence without opening FITS, pixels, galaxies, randoms or odd data. Synthetic tampering of bit/RA source fails closed. [Its CI check passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36223996369).

**Immediate source request, not a new selection:** obtain the historical production script/configuration or a versioned official mask reference for bits 8–11 and both low-quality `eboss22` plates, including exact release, original byte fingerprints, RA/DEC convention and its application to data and randoms. If unavailable, document the missing provenance explicitly and keep any future released-random-defined cross-window validation conceptually distinct from a continuous physical-mask reconstruction. Preserve PR #1 as draft and the observed odd vector unopened.

## Paper-cited 2020 code tag compared with later public source

The paper's Data Availability explicitly points to
`cheng-zhao/brickmask` tag **`v1.0`** for bits 1–7. Its annotated
Git tag object `3a94b72f4ee39c215cab713cf464d54106899838`
has message **“Version for eBOSS DR16”**, date 2020-05-09 and
resolves to commit `4c0f940934ff8c4e6b0f0c709b3093383abff8d0`.
This is a stronger *published source-version anchor* than our
previously pinned, later public checkout
`b9eb684a579b56ec3dbdb46549224be7e3fa2830`
(2023-11-08). It is NOT independent proof that the historical
release production job actually ran the tag with any particular
configuration. The v1.0 Git tree does **not** contain the later
`scripts/eBOSS_ELG_extra.py` helper; therefore the 2023 vendored
copy cannot be treated as a byte-authenticated 2020 production
script for supplemental bits8–11.

Original MIT-licensed source files and each relevant version's
license have been archived in
`source_data/upstream/brickmask_v1_0/` and
`source_data/upstream/brickmask_b9eb684/`, with each copied
Git blob independently verified as byte-identical to its source.
The frozen, machine-readable comparison is
`source_data/eboss_dr16_elg_brickmask_v1_0_source_comparison_protocol_2026-09-26.json`.
In both **reviewed 8-bit, same-rounded-pixel and same-truncated-pixel
code paths**, validity is bit0 and the xybug correction uses
bit2 = 4. The tag hardcodes exactly the four published
`eboss21,22,23,25` source filenames. The v1.0 program writes
`VETOMASK` and `MCHUNK` (the latter is a `1<<i` family
bitset); the later public `-DEBOSS` branch has a configurable
maskbit column and optional single `SUBID`, overwritten by the
latest valid checked subsample. These are NOT interchangeable
raw output fields in overlapping-family cases. The released
full catalogue's `mskbit` field has not yet been linked
byte-by-byte to the v1.0 `VETOMASK` production output.

A source-only auditor,
`scripts/audit_eboss_dr16_elg_brickmask_v1_0_source.py`,
checks all **eight** exact upstream source Git blobs, selected
unmodified C assignment/WCS anchors and a fully synthetic
exhaustive `256×256 = 65,536` rounded/truncated byte-pair
bit-code algebra test. It also exercises a two-valid-family
counterexample for `MCHUNK` versus `SUBID`, and a tampered
historical bit2 macro must fail. [Source-only CI passed](https://github.com/dvlahek/stress-energy-closure/actions/runs/36224705917).
It does not test equivalence of the two actual WCS
world-to-pixel functions or a production FITS run; it opens
no observed or random rows and does not touch the odd vector.

## Newly recovered paper-specified bad-plate identifiers

Raichoor et al. §3.2 item (xii) specify the exact two
**`eboss22` PLATE-MJD pairs `9430-58112` and
`9395-58113`**. Table 5 attributes 13.9 deg² and 3,123
removed targets to the pair's non-bitcoded plate exclusion.
Primary paper: https://academic.oup.com/mnras/article/500/3/3254/5942664 .
The discovery is separately recorded in
`source_data/eboss_dr16_elg_eboss22_badplate_published_ids_2026-09-26.json`.
Knowing the two IDs closes the **identifier** provenance
item only. It does not determine an exact continuous veto
polygon or the original handling of plate overlap and
the published data/random selections. In particular,
do not impose a guessed full circular plate radius, reapply
these plate vetoes to released clustering products, or
compare odd data to adjudicate geometry.

The remaining external source request can now be narrow:
the actual 2020 DR16 supplemental Python/configuration
for bits8–11, and the exact mask geometry/application
to **PLATE-MJD 9430-58112 and 9395-58113** for both
released ELG data and randoms. Accept only an attributable
versioned source/reference and record its original byte
fingerprints before any new selection-level check.
