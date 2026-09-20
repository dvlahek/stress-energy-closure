#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-data/desi_dr1_lrg_elg}"
OUT="${2:-lrg_elg_exact_zresolved_r4_desi}"
NTHREADS="${3:-16}"
NRANDOM="${4:-4}"
REGION="${5:-both}"

lr=(); er=(); ld=(); ed=()
case "$REGION" in
  both)
    ld=("${ROOT}/LRG_NGC_clustering.dat.fits" "${ROOT}/LRG_SGC_clustering.dat.fits")
    ed=("${ROOT}/ELG_LOPnotqso_NGC_clustering.dat.fits" "${ROOT}/ELG_LOPnotqso_SGC_clustering.dat.fits")
    ;;
  NGC)
    ld=("${ROOT}/LRG_NGC_clustering.dat.fits")
    ed=("${ROOT}/ELG_LOPnotqso_NGC_clustering.dat.fits")
    ;;
  SGC)
    ld=("${ROOT}/LRG_SGC_clustering.dat.fits")
    ed=("${ROOT}/ELG_LOPnotqso_SGC_clustering.dat.fits")
    ;;
  *) echo "REGION must be both, NGC, or SGC" >&2; exit 2 ;;
esac

for ((r=0; r<NRANDOM; r++)); do
  case "$REGION" in
    both)
      lr+=("${ROOT}/LRG_NGC_${r}_clustering.ran.fits" "${ROOT}/LRG_SGC_${r}_clustering.ran.fits")
      er+=("${ROOT}/ELG_LOPnotqso_NGC_${r}_clustering.ran.fits" "${ROOT}/ELG_LOPnotqso_SGC_${r}_clustering.ran.fits")
      ;;
    NGC)
      lr+=("${ROOT}/LRG_NGC_${r}_clustering.ran.fits")
      er+=("${ROOT}/ELG_LOPnotqso_NGC_${r}_clustering.ran.fits")
      ;;
    SGC)
      lr+=("${ROOT}/LRG_SGC_${r}_clustering.ran.fits")
      er+=("${ROOT}/ELG_LOPnotqso_SGC_${r}_clustering.ran.fits")
      ;;
  esac
done

python code/desi_dr1_lrg_elg_exact_zresolved.py \
  --lrg-data "${ld[@]}" \
  --elg-data "${ed[@]}" \
  --lrg-random "${lr[@]}" \
  --elg-random "${er[@]}" \
  --outdir "$OUT" \
  --z-edges 0.80,0.90,1.00,1.10 \
  --mu-bins 240 \
  --theta-min-deg 0.05 \
  --distance-cosmology desi \
  --nthreads "$NTHREADS"
