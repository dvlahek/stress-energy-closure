# E53 local run

Copy `e53_local_conditioned_mock0001_sample_scaling.py` into
`~/stress-energy-closure/scripts/`.

First run only the safe code self-test:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e53_local_conditioned_mock0001_sample_scaling.py --self-test
```

Then Stage 1 only:

```bash
python3 -u scripts/e53_local_conditioned_mock0001_sample_scaling.py --run \
  2>&1 | tee e53_stage1_fullD_4800.log
```

Do NOT add `--include-48000` on the first run. Stage 1 checkpoints every term.
After inspecting Stage 1, Stage 2 can resume the same JSON and add 48000R:

```bash
python3 -u scripts/e53_local_conditioned_mock0001_sample_scaling.py --run --include-48000 \
  2>&1 | tee e53_stage2_fullD_48000.log
```

Checkpoint/result:
`source_data/e53_conditioned_mock0001_sample_scaling.json`
