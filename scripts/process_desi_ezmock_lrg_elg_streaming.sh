#!/usr/bin/env bash
set -euo pipefail

# Stream DESI DR1 EZmock LRG x ELG realizations through the exact estimator.
# Downloads one mock at a time, runs the estimator, verifies the output, then
# deletes the large FITS inputs. This keeps disk usage roughly constant.
#
# Example:
#   bash scripts/process_desi_ezmock_lrg_elg_streaming.sh --start 1 --count 20 --nrandom 1 --nthreads 16

START=1
COUNT=20
NRANDOM=1
NTHREADS=16
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTROOT="mocks"
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
    -h|--help)
      cat <<EOF
Usage: $0 [--start N] [--count N] [--nrandom N] [--nthreads N] [--inroot DIR] [--outroot DIR] [--keep-inputs]

Defaults:
  --start 1
  --count 20
  --nrandom 1
  --nthreads 16
  --inroot data/desi_dr1_ezmock_dark_v1
  --outroot mocks

By default each mock input directory is deleted only after a successful exact output is verified.
EOF
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

END=$((START + COUNT - 1))
if (( START < 1 || END > 1000 )); then
  echo "Requested mock range ${START}..${END} is outside 1..1000" >&2
  exit 2
fi
if (( NRANDOM < 1 || NRANDOM > 18 )); then
  echo "--nrandom must be between 1 and 18" >&2
  exit 2
fi

mkdir -p "$INROOT" "$OUTROOT"

for ((m=START; m<=END; m++)); do
  out=$(printf "%s/exact_%04d" "$OUTROOT" "$m")
  result="${out}/lrg_elg_exact_odd_multipoles.csv"
  indir="${INROOT}/mock${m}"

  if [[ -s "$result" ]]; then
    echo "SKIP mock${m}: exact output already exists"
    if (( KEEP_INPUTS == 0 )) && [[ -d "$indir" ]]; then
      echo "CLEAN existing inputs: $indir"
      rm -rf -- "$indir"
    fi
    continue
  fi

  echo "===== STREAM mock${m} ====="
  bash scripts/download_desi_ezmock_lrg_elg.sh \
    --start "$m" --count 1 --nrandom "$NRANDOM" --outroot "$INROOT"

  bash scripts/run_desi_ezmock_lrg_elg_exact.sh \
    --start "$m" --count 1 --nrandom "$NRANDOM" \
    --inroot "$INROOT" --outroot "$OUTROOT" --nthreads "$NTHREADS"

  if [[ ! -s "$result" ]]; then
    echo "ERROR mock${m}: expected result missing; preserving inputs for debugging: $indir" >&2
    exit 4
  fi

  if (( KEEP_INPUTS == 0 )); then
    echo "CLEAN mock${m} inputs: $indir"
    rm -rf -- "$indir"
  fi

  echo "DONE mock${m}: $result"
done

echo "DONE streamed mocks ${START}..${END}, nrandom=${NRANDOM}"
