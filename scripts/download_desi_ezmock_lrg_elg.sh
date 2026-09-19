#!/usr/bin/env bash
set -euo pipefail

# Download DESI DR1 EZmock FFA LRG x ELG clustering catalogs.
# Usage:
#   bash scripts/download_desi_ezmock_lrg_elg.sh --start 1 --count 20 --nrandom 1
#   bash scripts/download_desi_ezmock_lrg_elg.sh --start 1 --count 1000 --nrandom 4

START=1
COUNT=20
NRANDOM=1
OUTROOT="data/desi_dr1_ezmock_dark_v1"
BASE="https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/EZmock/dark/v1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nrandom) NRANDOM="$2"; shift 2 ;;
    --outroot) OUTROOT="$2"; shift 2 ;;
    --base) BASE="$2"; shift 2 ;;
    -h|--help)
      cat <<EOF
Usage: $0 [--start N] [--count N] [--nrandom N] [--outroot DIR]

Defaults:
  --start 1
  --count 20
  --nrandom 1
  --outroot data/desi_dr1_ezmock_dark_v1

NRANDOM must be between 1 and 18.
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
if (( START < 1 || START > 1000 )); then
  echo "--start must be between 1 and 1000" >&2
  exit 2
fi
END=$((START + COUNT - 1))
if (( END > 1000 )); then
  echo "Requested end mock $END exceeds mock1000" >&2
  exit 2
fi

download_one() {
  local url="$1"
  local dst="$2"
  mkdir -p "$(dirname "$dst")"
  if [[ -s "$dst" ]]; then
    echo "HAVE $dst"
    return
  fi
  echo "GET  $url"
  curl --http1.1 -fL --retry 8 --retry-delay 3 --retry-all-errors --continue-at - -o "$dst" "$url"
}

for ((m=START; m<=END; m++)); do
  src="${BASE}/mock${m}"
  dst="${OUTROOT}/mock${m}"
  mkdir -p "$dst"
  echo "=== mock${m} ==="

  for tracer in LRG_ffa ELG_LOP_ffa; do
    for cap in NGC SGC; do
      fn="${tracer}_${cap}_clustering.dat.fits"
      download_one "${src}/${fn}" "${dst}/${fn}"
      for ((r=0; r<NRANDOM; r++)); do
        fn="${tracer}_${cap}_${r}_clustering.ran.fits"
        download_one "${src}/${fn}" "${dst}/${fn}"
      done
    done
  done
done

echo "DONE mocks ${START}..${END}, nrandom=${NRANDOM}"
