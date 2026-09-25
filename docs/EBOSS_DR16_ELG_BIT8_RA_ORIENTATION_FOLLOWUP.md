# eBOSS DR16 ELG bit 8: independent RA-orientation follow-up

Status (2026-09-25): **source-only geometry checked; official production bit-8 convention NOT identified.** This is separate from the existing source/API unit discrepancy and does not certify a physical mask, reapply any veto, open observed odd data or modify the DESI results.

## Problem and limitations of the previous roundtrip

The byte-pinned upstream `cheng-zhao/brickmask` script at commit `b9eb684a579b56ec3dbdb46549224be7e3fa2830` forms
`theta=radians(90-dec)` and `phi=radians(360-ra)`, then calls `healpy.ang2pix(1024,theta,phi,nest=False,lonlat=True)`.
That call passes radians to an option that expects degrees. The previous frozen source-only test recovered all 37 listed pixel indices using `lonlat=False` after constructing test RA from `360-phi`, but that demonstrates the internal roundtrip of the *reflected-RA* transformation, **not** its compatibility with celestial RA in the public survey footprint.

A separate source-only follow-up is frozen in
`source_data/eboss_dr16_elg_bit8_ra_orientation_protocol_2026-09-25.json`.
It records the exploratory arithmetic that motivated the follow-up **before** the CI replay. The pinned upstream script is read and AST-parsed to extract the entire list of 37 source pixels, and its existing Git blob and integer-array SHA256 are verified. The follow-up never opens a random or galaxy catalogue.

## Independent native-RA versus reflected-RA check

For each published HEALPix pixel, derive its centre with
`hp.pix2ang(1024, pixel, nest=False, lonlat=True)`. The returned
longitude is used as *native HEALPix longitude*; treating that as celestial RA is a separate, explicitly testable hypothesis. At these same 37 native-longitude centres, evaluate four alternatives: the upstream call as written, the units-only repair with `360-RA`, native-RA spherical radians, and native-RA longitude/latitude degrees. Native-RA radians and native-RA degrees agree and recover all source indices. The reflected-RA units-only mapping returns no pixel in the 37-pixel source list at these native-longitude centres.

As a **diagnostic only**, the published approximate ELG NGC rectangle
`126<RA<166`, `13.8<DEC<32.5` and SGC rectangle
`-43<RA<45`, `-5<DEC<5` provide an external sky-location check.
With native HEALPix longitude as RA, 10 listed pixel centres are in the approximate NGC rectangle, 26 in SGC and one is just outside the rounded NGC declination limit (DEC approximately 32.531 degrees). With reflected RA, zero centres are in the approximate NGC rectangle, 26 in SGC and 11 outside. These rectangles are **not** the exact DR16 footprint; their counts must never define a physical cut, accepted area or pair window.

The 37 equal-area nside-1024 pixels have a combined geometric area of approximately `0.121304 deg^2`, compatible with the published one-decimal `0.1 deg^2`. Raichoor et al. (MNRAS 500, 3254–3274, 2021, doi:10.1093/mnras/staa3336), Sec. 3.2 and Table 5, describe the rejection of 37 pixels for bit 8 and list 15 removed targets. We have not independently checked those 15 catalogue targets.

Code: `scripts/audit_eboss_dr16_elg_bit8_ra_orientation.py`.
The frozen synthetic self-test and all-pixel source-only CI replay passed at commit
`5f3d88e7a9b9ba46ad0b78e627d92435477fbf4e`:
https://github.com/dvlahek/stress-energy-closure/actions/runs/36151465886

## What this establishes and what remains open

The existing script has a documented radian/`lonlat=True` API discrepancy. Its `360-RA` transformation also requires an independent physical-coordinate check: the former 37/37 reflected roundtrip is insufficient to validate the actual sky orientation of bit 8. These are **source-code and source-list** checks. They do not show how the official released DR16 clustering catalogues or published production mask were generated, so do not claim that those catalogues were incorrectly masked.

The independent reference gate requires an exact-version official production bitmap, or an independently verified **pre-veto bit-coded** output with coordinates, bit-8 flags, source provenance and SHA256. Test both positive and negative positions that discriminate the literal, reflected-RA and native-RA conventions; include released-mask and chunk handling only after the production reference is authenticated. An aggregate area, the number 15, or post-veto clustering catalogue survivors alone cannot distinguish candidate algorithms. If the official production reference is unavailable, record that unresolved limitation; do **not** manufacture it or silently replace the upstream call.

The LRG MANGLE veto/sector rules and four ELG brickmask FITS families can be audited in parallel with this reference search. The published LRG and ELG clustering randoms retain their own selection functions, and the cross-LS pair window remains the corresponding LRG×ELG `R_L R_E`; do not invent a shared mask from coarse random occupancy.

The observed eBOSS odd vector remains unopened. Random-only LS closure, 11/11 polygon SHA provenance and the 36/36 fine weighted n(z) check are already complete and do not need repetition.

## Next local WSL gate: official *reference filename* discovery only

The separately frozen protocol is
`source_data/eboss_dr16_elg_bit8_reference_discovery_protocol_2026-09-25.json`.
The runner `scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py` reads
only the official **DR16 root HTML index** using a 16-MiB bound. It checks
the previously recorded exact SHA256 of that index before accepting any
filenames. It does not repeat any polygon downloads or FITS parsing.

Run in the user's existing WSL checkout of this audit branch:

```bash
cd ~/stress-energy-closure
git pull --ff-only
source .venv/bin/activate
python scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py --self-test
python scripts/audit_eboss_dr16_elg_bit8_reference_discovery.py \
  --out eboss_workspace/official_mask_inventory/bit8_reference_index.json
```

A successful index match means **filenames only**, not that any independent
bit-8 reference was found. A changed official HTML listing fails closed and
must be reviewed separately; do not silently substitute its SHA. Inspect the
reported candidates for a genuinely bit-coded *pre-veto* official production
output. Post-veto `eBOSS_ELG_clustering_*` files and the ordinary SDSS
`EBOSS_TARGET1` ELG selection flags cannot certify the published ELG
BRICKMASK extra bit 8.

The networking-free synthetic self-test passed in
[workflow 36152677996](https://github.com/dvlahek/stress-energy-closure/actions/runs/36152677996).
The official root listing has **not** been retrieved by that CI test.
No reference bitmap was authenticated, and no physical mask/odd-data gate
was opened.
