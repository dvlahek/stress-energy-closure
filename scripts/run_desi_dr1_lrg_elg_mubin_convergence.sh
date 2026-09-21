#!/usr/bin/env bash
set -euo pipefail

# DESI DR1 full-sky r4 mu-bin convergence test.
# Baseline production uses 240 mu bins; this script computes 120 and 480
# with all other estimator choices frozen.
#
# Usage:
#   bash scripts/run_desi_dr1_lrg_elg_mubin_convergence.sh [data_root] [nthreads]
#
# Defaults:
#   data_root=data/desi_dr1_lrg_elg
#   nthreads=16

ROOT="${1:-data/desi_dr1_lrg_elg}"
NTHREADS="${2:-16}"
NRANDOM=4

ld=("${ROOT}/LRG_NGC_clustering.dat.fits" "${ROOT}/LRG_SGC_clustering.dat.fits")
ed=("${ROOT}/ELG_LOPnotqso_NGC_clustering.dat.fits" "${ROOT}/ELG_LOPnotqso_SGC_clustering.dat.fits")
lr=()
er=()
for ((r=0; r<NRANDOM; r++)); do
  lr+=("${ROOT}/LRG_NGC_${r}_clustering.ran.fits" "${ROOT}/LRG_SGC_${r}_clustering.ran.fits")
  er+=("${ROOT}/ELG_LOPnotqso_NGC_${r}_clustering.ran.fits" "${ROOT}/ELG_LOPnotqso_SGC_${r}_clustering.ran.fits")
done

for MU in 120 480; do
  OUT="lrg_elg_exact_zresolved_r4_desi_mu${MU}"
  echo "===== DESI r4 full-sky mu-bins=${MU} ====="
  python code/desi_dr1_lrg_elg_exact_zresolved.py \
    --lrg-data "${ld[@]}" \
    --elg-data "${ed[@]}" \
    --lrg-random "${lr[@]}" \
    --elg-random "${er[@]}" \
    --outdir "${OUT}" \
    --z-edges 0.80,0.90,1.00,1.10 \
    --mu-bins "${MU}" \
    --theta-min-deg 0.05 \
    --distance-cosmology desi \
    --nthreads "${NTHREADS}" \
    --skip-reverse
done

echo "DONE mu-bin convergence: 120 and 480; baseline is frozen 240."
