#!/usr/bin/env bash
set -euo pipefail

# DESI DR1 full-sample five-tracer robustness rerun with 256 mark permutations.
# This keeps the primary analysis settings fixed and changes only the number
# of permutations from 32 to 256. Results are written to a separate directory.

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

DATA_NGC="BGS_BRIGHT-21.5_NGC_clustering.dat.fits"
DATA_SGC="BGS_BRIGHT-21.5_SGC_clustering.dat.fits"
RAND_NGC="BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits"
RAND_SGC="BGS_BRIGHT-21.5_SGC_0_clustering.ran.fits"
TEMPLATES="phase7_template/phase7_templates.csv"
OUTDIR="phase7_multitracer_fullsample_perm256"

for f in "$DATA_NGC" "$DATA_SGC" "$RAND_NGC" "$RAND_SGC" "$TEMPLATES"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing required input: $f" >&2
    exit 1
  fi
done

python code/desi_dr1_phase7_multitracer_fullsample.py \
  --data "$DATA_NGC" "$DATA_SGC" \
  --random "$RAND_NGC" "$RAND_SGC" \
  --templates "$TEMPLATES" \
  --outdir "$OUTDIR" \
  --zmin 0.10 --zmax 0.40 --dz-proxy 0.02 --ntracer 5 \
  --analysis-z-edges 0.10,0.20,0.30,0.40 \
  --sep-edges 20,40,60,80,100,120,140 \
  --random-factor 2.0 --neighbors-per-anchor 48 \
  --theta-min-deg 0.05 --jackknife 30 --permutations 256 \
  --seed 20260913

printf '\nCompleted 256-permutation robustness run.\nOutput: %s/%s\n' "$ROOT" "$OUTDIR"
