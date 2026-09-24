#!/usr/bin/env bash
# Process one of four nonoverlapping sets of 30 EZmock realizations (mock1–mock120).
set -euo pipefail

CHUNK="${1:-}"
NTHREADS="${2:-16}"
INROOT="${INROOT:-data/desi_dr1_ezmock_dark_v1}"
OUTROOT="${OUTROOT:-mocks_zresolved_desi_r4_120}"

if [[ ! "$CHUNK" =~ ^[1-4]$ ]]; then
  echo "Usage: $0 CHUNK(1..4) [nthreads]" >&2
  exit 2
fi

START=$(( (CHUNK - 1) * 30 + 1 ))
bash scripts/process_desi_ezmock_lrg_elg_zresolved_streaming.sh \
  --start "$START" \
  --count 30 \
  --nrandom 4 \
  --nthreads "$NTHREADS" \
  --inroot "$INROOT" \
  --outroot "$OUTROOT"
