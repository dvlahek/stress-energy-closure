# E59 local run (v2)

Copy `e59_mock_velocity_tag_repeatability.py` to `~/stress-energy-closure/scripts/`.

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1

python3 -u scripts/e59_mock_velocity_tag_repeatability.py --self-test
```

Expected:

```text
E59_SYNTHETIC_VELOCITY_RECONSTRUCTION_SELF_TEST_PASS
```

Then:

```bash
python3 -u scripts/e59_mock_velocity_tag_repeatability.py --run \
  2>&1 | tee e59_mock_velocity_tag_repeatability.log
```

The run uses only already-local EZmocks. It re-hashes all 72 compressed inputs before FITS rows. Observed galaxy rows and the observed odd vector remain sealed.

Output: `source_data/e59_mock_velocity_tag_repeatability.json`.
