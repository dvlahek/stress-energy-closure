#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-data/desi_dr1_lrg_elg}"
OUT="${2:-lrg_elg_rr_window_r4}"
NTHREADS="${3:-10}"
NRANDOM="${4:-4}"

lr=()
er=()
for ((r=0; r<NRANDOM; r++)); do
  lr+=("${ROOT}/LRG_NGC_${r}_clustering.ran.fits" "${ROOT}/LRG_SGC_${r}_clustering.ran.fits")
  er+=("${ROOT}/ELG_LOPnotqso_NGC_${r}_clustering.ran.fits" "${ROOT}/ELG_LOPnotqso_SGC_${r}_clustering.ran.fits")
done

python code/build_lrg_elg_zresolved_rr_window.py \
  --lrg-random "${lr[@]}" \
  --elg-random "${er[@]}" \
  --outdir "${OUT}" \
  --z-edges 0.80,0.90,1.00,1.10 \
  --output-sep-edges 20,40,60,80,100,120,140 \
  --fine-step 1.0 \
  --mu-bins 240 \
  --theta-min-deg 0.05 \
  --distance-cosmology desi \
  --nthreads "${NTHREADS}"
