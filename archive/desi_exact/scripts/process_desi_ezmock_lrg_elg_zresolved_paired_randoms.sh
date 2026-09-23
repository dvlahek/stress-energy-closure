#!/usr/bin/env bash
set -euo pipefail

START=1
COUNT=10
NTHREADS=16
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTR1="mocks_zresolved_desi_pair_r1"
OUTR4="mocks_zresolved_desi_pair_r4"
KEEP_INPUTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --inroot) INROOT="$2"; shift 2 ;;
    --out-r1) OUTR1="$2"; shift 2 ;;
    --out-r4) OUTR4="$2"; shift 2 ;;
    --keep-inputs) KEEP_INPUTS=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

END=$((START+COUNT-1))
for ((m=START; m<=END; m++)); do
  r1=$(printf "%s/exact_z_%04d/lrg_elg_exact_zresolved_odd_multipoles.csv" "$OUTR1" "$m")
  r4=$(printf "%s/exact_z_%04d/lrg_elg_exact_zresolved_odd_multipoles.csv" "$OUTR4" "$m")
  indir="${INROOT}/mock${m}"

  if [[ -s "$r1" && -s "$r4" ]]; then
    echo "SKIP mock${m}: paired r1/r4 outputs exist"
    if (( KEEP_INPUTS == 0 )) && [[ -d "$indir" ]]; then rm -rf -- "$indir"; fi
    continue
  fi

  echo "===== PAIRED DESI-COSMOLOGY mock${m} ====="
  bash scripts/download_desi_ezmock_lrg_elg.sh \
    --start "$m" --count 1 --nrandom 4 --outroot "$INROOT"

  if [[ ! -s "$r1" ]]; then
    bash scripts/run_desi_ezmock_lrg_elg_exact_zresolved.sh \
      --start "$m" --count 1 --nrandom 1 \
      --inroot "$INROOT" --outroot "$OUTR1" --nthreads "$NTHREADS"
  fi

  if [[ ! -s "$r4" ]]; then
    bash scripts/run_desi_ezmock_lrg_elg_exact_zresolved.sh \
      --start "$m" --count 1 --nrandom 4 \
      --inroot "$INROOT" --outroot "$OUTR4" --nthreads "$NTHREADS"
  fi

  if [[ ! -s "$r1" || ! -s "$r4" ]]; then
    echo "ERROR mock${m}: paired output missing; preserving $indir" >&2
    exit 4
  fi

  if (( KEEP_INPUTS == 0 )); then rm -rf -- "$indir"; fi
  echo "DONE paired mock${m}"
done

echo "DONE paired DESI-cosmology r1/r4 mocks ${START}..${END}"
