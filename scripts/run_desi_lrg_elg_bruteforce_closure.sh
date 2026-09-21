#!/usr/bin/env bash
set -euo pipefail

# Small deterministic NGC z=0.8-0.9 brute-force closure.
# Uses one random realization because this is an estimator algebra/geometry
# validation, not a covariance or random-density test.
#
# Usage:
#   bash scripts/run_desi_lrg_elg_bruteforce_closure.sh [data_root] [nthreads]

ROOT="${1:-data/desi_dr1_lrg_elg}"
NTHREADS="${2:-8}"
OUT="${OUT:-lrg_elg_bruteforce_pair_closure}"

python code/audit_lrg_elg_bruteforce_pair_closure.py \
  --lrg-data "${ROOT}/LRG_NGC_clustering.dat.fits" \
  --elg-data "${ROOT}/ELG_LOPnotqso_NGC_clustering.dat.fits" \
  --lrg-random "${ROOT}/LRG_NGC_0_clustering.ran.fits" \
  --elg-random "${ROOT}/ELG_LOPnotqso_NGC_0_clustering.ran.fits" \
  --outdir "${OUT}" \
  --zlo 0.80 \
  --zhi 0.90 \
  --ndata 1200 \
  --nrandom 6000 \
  --seed 1729 \
  --sep-edges 20,40,60,80,100,120,140 \
  --mu-bins 24 \
  --theta-min-deg 0.05 \
  --distance-cosmology desi \
  --nthreads "${NTHREADS}" \
  --block 256
