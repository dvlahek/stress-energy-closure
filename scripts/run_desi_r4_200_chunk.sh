#!/usr/bin/env bash
set -euo pipefail

# Run one of eight fixed 25-mock production chunks for the DESI-fiducial r4
# full-sky z-resolved estimator.
#
# Usage:
#   bash scripts/run_desi_r4_200_chunk.sh 1 [nthreads]
#   ...
#   bash scripts/run_desi_r4_200_chunk.sh 8 [nthreads]

CHUNK="${1:-}"
NTHREADS="${2:-16}"
OUTROOT="${OUTROOT:-mocks_zresolved_desi_r4_200}"
INROOT="${INROOT:-data/desi_dr1_ezmock_dark_v1}"

if [[ ! "$CHUNK" =~ ^[1-8]$ ]]; then
  echo "Usage: $0 CHUNK(1..8) [nthreads]" >&2
  exit 2
fi

START=$(( (CHUNK - 1) * 25 + 1 ))
COUNT=25
END=$(( START + COUNT - 1 ))

echo "DESI r4 production chunk $CHUNK/8: mocks $START..$END"
echo "OUTROOT=$OUTROOT"
echo "NTHREADS=$NTHREADS"

bash scripts/process_desi_ezmock_lrg_elg_zresolved_streaming.sh \
  --start "$START" \
  --count "$COUNT" \
  --nrandom 4 \
  --nthreads "$NTHREADS" \
  --inroot "$INROOT" \
  --outroot "$OUTROOT"

echo "DONE chunk $CHUNK/8: mocks $START..$END"
