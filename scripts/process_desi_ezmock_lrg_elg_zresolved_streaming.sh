#!/usr/bin/env bash
set -euo pipefail

START=1
COUNT=40
NRANDOM=1
NTHREADS=16
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTROOT="mocks_zresolved"
KEEP_INPUTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nrandom) NRANDOM="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --inroot) INROOT="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --keep-inputs) KEEP_INPUTS=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

END=$((START+COUNT-1))
for ((m=START; m<=END; m++)); do
  out=$(printf "%s/exact_z_%04d" "$OUTROOT" "$m")
  result="${out}/lrg_elg_exact_zresolved_odd_multipoles.csv"
  indir="${INROOT}/mock${m}"

  if [[ -s "$result" ]]; then
    echo "SKIP mock${m}: z-resolved output exists"
    if (( KEEP_INPUTS == 0 )) && [[ -d "$indir" ]]; then rm -rf -- "$indir"; fi
    continue
  fi

  echo "===== STREAM ZRES mock${m} ====="
  bash scripts/download_desi_ezmock_lrg_elg.sh \
    --start "$m" --count 1 --nrandom "$NRANDOM" --outroot "$INROOT"

  bash scripts/run_desi_ezmock_lrg_elg_exact_zresolved.sh \
    --start "$m" --count 1 --nrandom "$NRANDOM" \
    --inroot "$INROOT" --outroot "$OUTROOT" --nthreads "$NTHREADS"

  if [[ ! -s "$result" ]]; then
    echo "ERROR mock${m}: expected z-resolved result missing; preserving $indir" >&2
    exit 4
  fi
  if (( KEEP_INPUTS == 0 )); then
    echo "CLEAN mock${m} inputs: $indir"
    rm -rf -- "$indir"
  fi
  echo "DONE ZRES mock${m}: $result"
done

echo "DONE z-resolved mocks ${START}..${END}, nrandom=${NRANDOM}"
