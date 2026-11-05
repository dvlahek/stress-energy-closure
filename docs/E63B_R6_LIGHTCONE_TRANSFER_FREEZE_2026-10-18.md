# E63B-R6 to E63D frozen lightcone-transfer record

**Science decision date:** 18 October 2026  
**Freeze record created:** 5 November 2026  
**Repository:** `dvlahek/stress-energy-closure`  
**Branch:** `audit/eboss-elg-bit8-ra-orientation-20260925`

## Final decision

E62 remains a PASS for the periodic Quijote N-body truth calibration.

E63B-R6 is a valid frozen FAIL for transfer of the unchanged reconstruction
method into the OuterRim/CosmoDC2 observer lightcone.

E64 is not authorized.

Observed eBOSS galaxy rows and the sealed observed odd vector remain SEALED.

A future lightcone-aware estimator must be declared as a separate
methodological branch and cannot replace, relabel or tune this frozen result.

## R6 geometry and retrieval

The geometry-only repair locked the center at:

- RA = 62.5 deg
- Dec = -35.5 deg

before truth evaluation.

The hybrid R6 retrieval produced:

- raw rows = 5,339,896
- exact-cone unique centrals = 200,493
- removed duplicates = 44,213
- minimum frozen outer-sector coverage = 1,460

The large local cache is not committed.

Its SHA-256 is:

`74d0402bfb2dc87b7ba5f174fb43adb141db4dc21c7f7b1b8acc9a20eefa1f82`

## Frozen science configuration

- selected sources = 54,893
- probes = 512
- probe seed = 202610630001
- smoothing scale R = 16 Mpc/h
- primary cutoff = 256 Mpc/h
- check cutoff = 192 Mpc/h
- truth Pearson/sign gates = 0.70 / 0.70
- cutoff Pearson/sign gates = 0.90 / 0.90

## E63B frozen result

Status:

`FAIL_BUFFERED_LIGHTCONE_TRUTH_TRANSFER_QUANTIFIED`

LOS:

- truth Pearson = 0.344988532283740
- truth sign agreement = 0.570312500
- cutoff Pearson = 0.993446128600096
- cutoff sign agreement = 0.888671875

Cartesian truth Pearson:

- x = 0.773132587327396
- y = 0.713213483349367
- z = -0.675174854395608

The lightcone transfer gate therefore fails.

## Independent velocity-convention audit

The catalog-redshift diagnostic independently confirms the original
`(+vx,+vy,+vz)` convention:

- Pearson = 0.999946973784688
- robust sign agreement = 1.000000000
- slope = 0.999769754322148

No velocity-axis or sign correction is licensed.

## E63C local-basis diagnostic

The parent E63B LOS result was exactly reproduced before decomposition.

Truth Pearson:

- radial `e_r` = 0.344988532283740
- tangential `e_RA` = 0.731278918539692
- tangential `e_Dec` = -0.152256467512530

The transfer failure is anisotropic.

## E63D fixed-split diagnostic

Fixed diagnostic groups:

- z_true < 0.95: 249
- z_true >= 0.95: 263
- Dec < -35.5 deg: 231
- Dec >= -35.5 deg: 281

Radial truth Pearson:

- lower-z =
  0.350374139904348
- higher-z =
  0.237780060983288

Higher-redshift radial degradation is present, but neither the fixed
redshift split nor the fixed north/south split supplies a single
explanation for the full failure.

No new science gate or p-value was introduced.

## Frozen interpretation

The E62 periodic-box PASS is not sufficient evidence of observational
transferability.

Under the frozen E63B conditions the estimator is numerically stable with
respect to the 192/256 Mpc/h cutoff comparison, but does not recover the
truth velocity field sufficiently well in the observer lightcone.

The E63B-D diagnostic branch is therefore closed.

**Final branch decision:** preserve E63B FAIL, do not proceed to E64,
and keep observed eBOSS data and the observed odd vector SEALED.
