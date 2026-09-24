# Local eBOSS DR16 fine-RR runbook (data blind)

This is a reproducible local route for checking one cap and one candidate redshift
bin without rerunning the full eight-case GitHub workflow. It reads **only public
eBOSS random FITS**, verifies their SHA256 against the committed input audits,
and uses the same pinned Corrfunc stack, multiplicative eBOSS weights, fixed
0.05-degree angular exclusion, midpoint LOS, signed-mu grid, separation edges
and deterministic half-sample seeds as the GitHub workflow. Do **not** use
galaxy data or odd-sector outputs during this input-only gate.

## WSL shell

From the existing repository clone, while on `main`:

```bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate
python -m pip install -r requirements-desi.txt
python scripts/audit_eboss_dr16_fine_rr_allbins.py --self-test
```

The first local shard below is **SGC, index 3**, the already validated
diagnostic interval `0.9 <= z < 1.0`. This is a relatively small,
independent environment check. To run one of the previously missing lower
candidate intervals, change `--z-indices` to `0`, `1` or `2`.
`--caps` accepts `NGC` or `SGC`; `--threads` sets Corrfunc threads.

```bash
mkdir -p eboss_workspace/local_rr
set -o pipefail
python scripts/audit_eboss_dr16_fine_rr_allbins.py \
  --caps SGC --z-indices 3 --threads 4 \
  --reuse-verified-cache \
  --cache-dir eboss_workspace/local_rr/fits \
  --out eboss_workspace/local_rr/SGC_z3.json \
  2>&1 | tee eboss_workspace/local_rr/SGC_z3.log
```

The first run may download the two public SGC random FITS files. With
`--reuse-verified-cache`, their SHA256-verified copies remain on disk and
are hashed again before reuse in subsequent local shards. Without the flag,
the original CI behavior downloads fresh catalogues and removes them after
inspection. Make sure the WSL disk and RAM have sufficient free space.
Run only one expensive shard at a time unless available resources are known.

For a completed shard, the JSON status is `local_rr_shard_complete`.
`SGC_z3.npz` contains `fine_rr_SGC_z3_full`,
`fine_rr_SGC_z3_half_A`, `fine_rr_SGC_z3_half_B` and their individual
weighted normalizations, with the same naming as the full CI artifact.
The script checks that every fine `(s,mu)` cell has support through the
counter result, computes an independent direct coarse RR histogram, and
requires a fine-to-coarse normalized `L1` residual below `1e-8`.
(The observed number of positive cells is also recorded, not silently
assumed.) All separations and redshift cuts remain **candidate settings**.

After each completed cap/bin the script atomically writes a matching
`SGC_z3.checkpoint.npz` and `SGC_z3.checkpoint.json`. If a later case
or final report fails, earlier completed counts can still be recovered.
Do not treat a `partial_checkpoint` as a certified full survey window.
The CI workflow also uploads these checkpoints when a job exits early.

## Comparing and sharing results

Use the final JSON and NPZ together, plus the relevant workflow run URL and
repository commit. Compare the input file hashes, cap, bin index, fixed
normalization, `fine_to_independent_coarse` and the two disjoint-half
differences. The local run may differ at low floating-point precision from
the CI run because thread count and compiled environment can change the
order of weighted summation. Compare with the numerical closure tolerance,
not byte-identical NPZ files.

Running a local shard does **not** validate an eBOSS physical survey window,
the joint mock covariance or the observed odd sector. The real galaxy odd
data vector must remain unopened until the prospective protocol is frozen.
