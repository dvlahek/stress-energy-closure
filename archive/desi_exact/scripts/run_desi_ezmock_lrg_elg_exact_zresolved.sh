#!/usr/bin/env bash
set -euo pipefail

START=1
COUNT=20
NRANDOM=1
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTROOT="mocks_zresolved_desi"
NTHREADS=16
REGION="both"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nrandom) NRANDOM="$2"; shift 2 ;;
    --inroot) INROOT="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --region) REGION="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

END=$((START+COUNT-1))
for ((m=START; m<=END; m++)); do
  indir="${INROOT}/mock${m}"
  out=$(printf "%s/exact_z_%04d" "$OUTROOT" "$m")
  result="${out}/lrg_elg_exact_zresolved_odd_multipoles.csv"
  if [[ -s "$result" ]]; then
    echo "SKIP mock${m}: $result exists"
    continue
  fi

  lr=(); er=(); ld=(); ed=()
  case "$REGION" in
    both)
      ld=("${indir}/LRG_ffa_NGC_clustering.dat.fits" "${indir}/LRG_ffa_SGC_clustering.dat.fits")
      ed=("${indir}/ELG_LOP_ffa_NGC_clustering.dat.fits" "${indir}/ELG_LOP_ffa_SGC_clustering.dat.fits")
      ;;
    NGC)
      ld=("${indir}/LRG_ffa_NGC_clustering.dat.fits")
      ed=("${indir}/ELG_LOP_ffa_NGC_clustering.dat.fits")
      ;;
    SGC)
      ld=("${indir}/LRG_ffa_SGC_clustering.dat.fits")
      ed=("${indir}/ELG_LOP_ffa_SGC_clustering.dat.fits")
      ;;
    *) echo "REGION must be both, NGC, or SGC" >&2; exit 2 ;;
  esac
  for ((r=0; r<NRANDOM; r++)); do
    case "$REGION" in
      both)
        lr+=("${indir}/LRG_ffa_NGC_${r}_clustering.ran.fits" "${indir}/LRG_ffa_SGC_${r}_clustering.ran.fits")
        er+=("${indir}/ELG_LOP_ffa_NGC_${r}_clustering.ran.fits" "${indir}/ELG_LOP_ffa_SGC_${r}_clustering.ran.fits")
        ;;
      NGC)
        lr+=("${indir}/LRG_ffa_NGC_${r}_clustering.ran.fits")
        er+=("${indir}/ELG_LOP_ffa_NGC_${r}_clustering.ran.fits")
        ;;
      SGC)
        lr+=("${indir}/LRG_ffa_SGC_${r}_clustering.ran.fits")
        er+=("${indir}/ELG_LOP_ffa_SGC_${r}_clustering.ran.fits")
        ;;
    esac
  done

  python code/desi_dr1_lrg_elg_exact_zresolved.py \
    --lrg-data "${ld[@]}" \
    --elg-data "${ed[@]}" \
    --lrg-random "${lr[@]}" \
    --elg-random "${er[@]}" \
    --outdir "$out" \
    --z-edges 0.80,0.90,1.00,1.10 \
    --mu-bins 240 \
    --theta-min-deg 0.05 \
    --distance-cosmology desi \
    --nthreads "$NTHREADS" \
    --skip-reverse
done
