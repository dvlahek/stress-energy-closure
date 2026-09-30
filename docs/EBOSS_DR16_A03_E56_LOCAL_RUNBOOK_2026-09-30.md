# E56 local run

Requirements already present locally:

- `source_data/e54_conditioned_mock0001_48k_random_replica_floor.json`
- `source_data/e55_conditioned_full_eligible_ninemock_two48k_result.json`
- the fixed E55 script/module in `scripts/`

Copy:
`e56_local_conditioned_full_eligible_ninemock_four48k.py`
to `~/stress-energy-closure/scripts/`.

Smoke test:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e56_local_conditioned_full_eligible_ninemock_four48k.py --self-test
```

Expected:

`E56_SYNTHETIC_FOUR_REPLICA_SELF_TEST_PASS`

Production:

```bash
python3 -u scripts/e56_local_conditioned_full_eligible_ninemock_four48k.py --run \
  2>&1 | tee e56_ninemock_four48k.log
```

E56 reuses all completed R1/R2 work. For mock0001 it reuses E54 R1-R4, so new
pair counting is required only for R3/R4 of the other eight mock IDs.

Checkpoint/result:

`source_data/e56_conditioned_full_eligible_ninemock_four48k_result.json`

The run checkpoints after every random-dependent pair term and can be resumed
with the same command.
