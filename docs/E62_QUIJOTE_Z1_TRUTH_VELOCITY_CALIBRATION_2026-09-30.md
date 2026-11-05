# E62 — Quijote z=1 N-body truth velocity calibration

## Purpose

E61 locked the complete eligible eBOSS ELG random pool as the deterministic survey-selection operator. E62 tests a different question: does the exact E59/E61 density-to-velocity reconstruction recover the direction of a **known true N-body halo velocity field**?

This stage deliberately uses an independent complete periodic simulation box. It is a reconstruction-method truth calibration, not an eBOSS cut-sky transfer calibration and not an Einstein–Vlasov signal-amplitude calibration.

## Frozen public dataset

Use Quijote:

- cosmology: `fiducial`
- realization: `0`
- halo catalogue: `FoF`
- snapshot: `2` (nominal z=1)
- source directory in the public Quijote data tree: `/Halos/FoF/fiducial/0/`

The Quijote documentation states that FoF positions are `GroupPos/1e3` in Mpc/h and peculiar velocities are `GroupVel*(1+z)` in km/s. The simulation box is 1 Gpc/h. The public halo-statistics table fixes `Ntot=165107` for the fiducial z=1 sample; E62 therefore uses exactly the 165107 most massive FoF halos, without post-result threshold tuning.

Public documentation:

- https://quijote-simulations.readthedocs.io/en/latest/halos.html
- https://quijote-simulations.readthedocs.io/en/latest/halo_statistics.html
- https://quijote-simulations.readthedocs.io/en/latest/access.html

## Frozen reconstruction

The primary reconstruction is unchanged from E59:

- R = 16 Mpc/h top-hat smoothing
- primary cutoff = 256 Mpc/h
- cutoff check = 192 Mpc/h
- unit tracer-number weights
- periodic minimum-image geometry
- 4096 deterministic probe halos, seed `202609620000`

Because the box is complete and periodic, there is no survey selection mask and no random-catalogue subtraction.

The primary truth gate is applied separately to x, y, and z:

- Pearson(reconstruction, true velocity) >= 0.70
- sign agreement >= 0.70
- Pearson(192,256) >= 0.90
- sign agreement(192,256) >= 0.90

All three axes must pass. Thresholds are frozen before any E62 result.

## Local data preparation

Do not download a Quijote snapshot. Only the FoF halo catalogue directory is needed.

Use the official Quijote Globus collection and transfer only:

`/Halos/FoF/fiducial/0/`

to a local directory, for example:

`~/stress-energy-closure/eboss_workspace/quijote/Halos/FoF/fiducial/0/`

If the project environment does not provide `readfof`, install the stable PyPI package with `python -m pip install Pylians` in the project venv before running E62.

## Run

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1

python3 -u scripts/e62_quijote_z1_truth_velocity_calibration.py --self-test
```

Expected:

```text
E62_SYNTHETIC_PERIODIC_TRUTH_CALIBRATION_SELF_TEST_PASS
```

Then:

```bash
python3 -u scripts/e62_quijote_z1_truth_velocity_calibration.py \
  --run \
  --catalog-dir "$HOME/stress-energy-closure/eboss_workspace/quijote/Halos/FoF/fiducial/0" \
  2>&1 | tee e62_quijote_z1_truth_velocity_calibration.log
```

Output:

`source_data/e62_quijote_z1_truth_velocity_calibration_result.json`

## Interpretation

A PASS authorizes only a later end-to-end survey/lightcone truth-transfer test. It does not show that the eBOSS ELG cut-sky reconstruction has the same truth correlation, does not validate an ELG HOD, does not simulate the custom F+/F- neutrino distributions, and does not calibrate an absolute Einstein–Vlasov amplitude.

Observed eBOSS galaxy rows and the observed odd vector remain SEALED.
