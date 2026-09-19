#!/usr/bin/env bash
set -euo pipefail

# Run the exact LRG x ELG odd estimator on downloaded DESI DR1 EZmocks.
# Usage:
#   bash scripts/run_desi_ezmock_lrg_elg_exact.sh --start 1 --count 20 --nrandom 1
#   bash scripts/run_desi_ezmock_lrg_elg_exact.sh --start 1 --count 1000 --nrandom 4

START=1
COUNT=20
NRANDOM=1
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTROOT="mocks"
NTHREADS=""
ZMIN="0.80"
ZMAX="1.10"
MUBINS="240"
THETA="0.05"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nrandom) NRANDOM="$2"; shift 2 ;;
    --inroot) INROOT="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    -h|--help)
      cat <<EOF
Usage: $0 [--start N] [--count N] [--nrandom N] [--inroot DIR] [--outroot DIR] [--nthreads N]

Defaults:
  --start 1
  --count 20
  --nrandom 1
  --inroot data/desi_dr1_ezmock_dark_v1
  --outroot mocks

Final covariance must use the same NRANDOM as the frozen data vector.
EOF
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if (( NRANDOM < 1 || NRANDOM > 18 )); then
  echo "--nrandom must be between 1 and 18" >&2
  exit 2
fi
END=$((START + COUNT - 1))
if (( START < 1 || END > 1000 )); then
  echo "Requested mock range ${START}..${END} is outside 1..1000" >&2
  exit 2
fi

for ((m=START; m<=END; m++)); do
  indir="${INROOT}/mock${m}"
  out=$(printf "%s/exact_%04d" "$OUTROOT" "$m")
  result="${out}/lrg_elg_exact_odd_multipoles.csv"

  if [[ -s "$result" ]]; then
    echo "SKIP mock${m}: $result exists"
    continue
  fi

  lrg_random=()
  elg_random=()
  for ((r=0; r<NRANDOM; r++)); do
    lrg_random+=(
      "${indir}/LRG_ffa_NGC_${r}_clustering.ran.fits"
      "${indir}/LRG_ffa_SGC_${r}_clustering.ran.fits"
    )
    elg_random+=(
      "${indir}/ELG_LOP_ffa_NGC_${r}_clustering.ran.fits"
      "${indir}/ELG_LOP_ffa_SGC_${r}_clustering.ran.fits"
    )
  done

  required=(
    "${indir}/LRG_ffa_NGC_clustering.dat.fits"
    "${indir}/LRG_ffa_SGC_clustering.dat.fits"
    "${indir}/ELG_LOP_ffa_NGC_clustering.dat.fits"
    "${indir}/ELG_LOP_ffa_SGC_clustering.dat.fits"
    "${lrg_random[@]}"
    "${elg_random[@]}"
  )
  for f in "${required[@]}"; do
    if [[ ! -s "$f" ]]; then
      echo "Missing input for mock${m}: $f" >&2
      exit 3
    fi
  done

  echo "=== exact mock${m} -> $out ==="
  cmd=(
    python code/desi_dr1_lrg_elg_exact.py
    --lrg-data
      "${indir}/LRG_ffa_NGC_clustering.dat.fits"
      "${indir}/LRG_ffa_SGC_clustering.dat.fits"
    --elg-data
      "${indir}/ELG_LOP_ffa_NGC_clustering.dat.fits"
      "${indir}/ELG_LOP_ffa_SGC_clustering.dat.fits"
    --lrg-random "${lrg_random[@]}"
    --elg-random "${elg_random[@]}"
    --outdir "$out"
    --zmin "$ZMIN" --zmax "$ZMAX"
    --mu-bins "$MUBINS"
    --theta-min-deg "$THETA"
    --skip-reverse
  )
  if [[ -n "$NTHREADS" ]]; then
    cmd+=(--nthreads "$NTHREADS")
  fi
  "${cmd[@]}"
done

echo "DONE exact mocks ${START}..${END}, nrandom=${NRANDOM}"
