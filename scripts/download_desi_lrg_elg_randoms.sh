#!/usr/bin/env bash
set -euo pipefail

N="${1:-4}"
OUT="${2:-data/desi_dr1_lrg_elg}"
BASE="https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5"

mkdir -p "$OUT"

fetch_one () {
  local rel="$1"
  local out="$OUT/$rel"
  if [ -s "$out" ]; then
    echo "exists: $out"
    return
  fi
  echo "downloading: $rel"
  curl -fL --retry 4 --retry-all-errors --retry-delay 5 \
    --connect-timeout 20 -o "$out" "$BASE/$rel"
}

for ((i=0; i<N; i++)); do
  fetch_one "LRG_NGC_${i}_clustering.ran.fits"
  fetch_one "LRG_SGC_${i}_clustering.ran.fits"
  fetch_one "ELG_LOPnotqso_NGC_${i}_clustering.ran.fits"
  fetch_one "ELG_LOPnotqso_SGC_${i}_clustering.ran.fits"
done

echo "Downloaded/verified random realizations 0..$((N-1)) in $OUT"
