# E54 local run

Requires local E53 Stage-2 JSON with PASS_FULLD_4800_AND_48000_CONDITIONED_SAMPLE_SCALING.

Smoke test:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e54_local_conditioned_mock0001_48k_random_replica_floor.py --self-test
```

Production:

```bash
python3 -u scripts/e54_local_conditioned_mock0001_48k_random_replica_floor.py --run \
  2>&1 | tee e54_48k_random_replica_floor.log
```

Checkpoint: source_data/e54_conditioned_mock0001_48k_random_replica_floor.json
