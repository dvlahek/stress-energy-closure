# E58 local run

Required parent:
`source_data/e57_conditioned_full_eligible_ninemock_six48k_result.json`

Smoke test:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e58_local_conditioned_full_eligible_ninemock_seven48k.py --self-test
```

Expected:
`E58_SYNTHETIC_SEVEN_REPLICA_SELF_TEST_PASS`

Production:

```bash
python3 -u scripts/e58_local_conditioned_full_eligible_ninemock_seven48k.py --run \
  2>&1 | tee e58_ninemock_seven48k.log
```

Only R7 is newly pair-counted. R1-R6 and D1D2 are reused from E57.

Checkpoint/result:
`source_data/e58_conditioned_full_eligible_ninemock_seven48k_result.json`
