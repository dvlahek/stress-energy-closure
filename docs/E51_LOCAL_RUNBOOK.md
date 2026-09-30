# E51 local run

Copy `e51_local_conditioned_marked_ls_ninemock.py` to
`~/stress-energy-closure/scripts/`.

The script reuses the already authenticated local A02 nine-mock sources. It
rehashes the frozen sources before FITS row access and then rereads only the
same deterministic 600D/1200R samples. It never touches observed catalogues.

Run:

```bash
cd ~/stress-energy-closure
source .venv/bin/activate 2>/dev/null || true
export PYTHONUNBUFFERED=1
python3 -u scripts/e51_local_conditioned_marked_ls_ninemock.py \
  2>&1 | tee e51_local_conditioned_marked_ls.log
```

A checkpoint is written after every completed cap/mock case to
`source_data/e51_conditioned_marked_ls_ninemock_result.json`.

For a code-only smoke test:

```bash
python3 -u scripts/e51_local_conditioned_marked_ls_ninemock.py --self-test
```
