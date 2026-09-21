#!/usr/bin/env bash
set -euo pipefail

# Run one of five fixed 40-mock production chunks for the DESI-fiducial r4
# full-sky z-resolved estimator.
#
# Usage:
#   bash scripts/run_desi_r4_200_chunk40.sh 1 [nthreads]
#   ...
#   bash scripts/run_desi_r4_200_chunk40.sh 5 [nthreads]

CHUNK="${1:-}"
NTHREADS="${2:-10}"
OUTROOT="${OUTROOT:-mocks_zresolved_desi_r4_200}"
INROOT="${INROOT:-data/desi_dr1_ezmock_dark_v1}"

if [[ ! "$CHUNK" =~ ^[1-5]$ ]]; then
  echo "Usage: $0 CHUNK(1..5) [nthreads]" >&2
  exit 2
fi

START=$(( (CHUNK - 1) * 40 + 1 ))
COUNT=40
END=$(( START + COUNT - 1 ))

echo "DESI r4 production chunk $CHUNK/5: mocks $START..$END"
echo "OUTROOT=$OUTROOT"
echo "NTHREADS=$NTHREADS"

bash scripts/process_desi_ezmock_lrg_elg_zresolved_streaming.sh \
  --start "$START" \
  --count "$COUNT" \
  --nrandom 4 \
  --nthreads "$NTHREADS" \
  --inroot "$INROOT" \
  --outroot "$OUTROOT"

echo "DONE chunk $CHUNK/5: mocks $START..$END"
