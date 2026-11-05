# E55 local run

Copy `e55_local_conditioned_full_eligible_ninemock_two48k.py` into `~/stress-energy-closure/scripts/`.

Keep the local completed E54 result at:
`source_data/e54_conditioned_mock0001_48k_random_replica_floor.json`.
Mock0001 R1/R2 will then be reused exactly.

Smoke test:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e55_local_conditioned_full_eligible_ninemock_two48k.py --self-test
```

Production:

```bash
python3 -u scripts/e55_local_conditioned_full_eligible_ninemock_two48k.py --run \
  2>&1 | tee e55_ninemock_fullD_two48k.log
```

The JSON checkpoints after DD and every random-dependent pair term:
`source_data/e55_conditioned_full_eligible_ninemock_two48k_result.json`.

The run can be resumed with the same command after interruption.
