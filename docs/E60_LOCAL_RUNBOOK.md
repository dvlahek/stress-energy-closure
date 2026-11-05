# E60 local run

Copy `e60_elg_velocity_reconstruction_mask_robustness.py` into
`~/stress-energy-closure/scripts/`.

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1

python3 -u scripts/e60_elg_velocity_reconstruction_mask_robustness.py --self-test
```

Expected:

```text
E60_SYNTHETIC_MASK_ROBUSTNESS_SELF_TEST_PASS
```

Then run:

```bash
python3 -u scripts/e60_elg_velocity_reconstruction_mask_robustness.py --run \
  2>&1 | tee e60_elg_velocity_reconstruction_mask_robustness.log
```

The script checkpoints after every random replica. It uses only the already-local nine EZmocks and keeps observations sealed.

Output:
`source_data/e60_elg_velocity_reconstruction_mask_robustness.json`.
