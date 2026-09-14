#!/usr/bin/env bash
set -euo pipefail

# Local Linux runner for the Phase-7 EZmock five-tracer placebo covariance control.
# Usage:
#   bash scripts/run_ezmock_placebo_local.sh 1
#   bash scripts/run_ezmock_placebo_local.sh 1 8
#   KEEP_FITS=1 bash scripts/run_ezmock_placebo_local.sh 1 2
#
# Large realization-specific random FITS are not required once a validated
# NGC/SGC random pair has been saved under shared_random_mock2/.  The launcher
# builds a compact angular-only cache once; each mock then redraws redshifts
# from its own data in the same narrow-z/cap cells.  This preserves the survey
# footprint without freezing mock2 radial selection and removes ~GB/mock of
# repeated random-catalog downloads.

FIRST="${1:-1}"
LAST="${2:-$FIRST}"
KEEP_FITS="${KEEP_FITS:-0}"
ROOT_BASE='https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/EZmock/bright/v1'
OUTROOT="${OUTROOT:-phase7_ezmock_local}"
TEMPLATE_DIR="${TEMPLATE_DIR:-phase7_template}"
TEMPLATE_CSV="$TEMPLATE_DIR/phase7_templates.csv"
SEED="${SEED:-20260913}"
DOWNLOAD_ATTEMPTS="${DOWNLOAD_ATTEMPTS:-40}"
ARIA_CONNECTIONS="${ARIA_CONNECTIONS:-8}"
RANDOM_FACTOR="${RANDOM_FACTOR:-2.0}"
SHARED_RANDOM_DIR="${SHARED_RANDOM_DIR:-$OUTROOT/shared_random_mock2}"
SHARED_RANDOM_CACHE="${SHARED_RANDOM_CACHE:-$SHARED_RANDOM_DIR/shared_random_angular.npz}"
SHARED_RANDOM_NGC_COUNT="${SHARED_RANDOM_NGC_COUNT:-1000000}"
SHARED_RANDOM_SGC_COUNT="${SHARED_RANDOM_SGC_COUNT:-500000}"

mkdir -p "$OUTROOT" "$TEMPLATE_DIR"

if command -v python >/dev/null 2>&1; then
  PYTHON_BIN="${PYTHON_BIN:-python}"
elif command -v python3 >/dev/null 2>&1; then
  PYTHON_BIN="${PYTHON_BIN:-python3}"
else
  echo 'ERROR: neither python nor python3 was found' >&2
  exit 2
fi
if ! command -v curl >/dev/null 2>&1; then
  echo 'ERROR: curl not found' >&2; exit 2
fi

echo "Using Python: $PYTHON_BIN ($($PYTHON_BIN --version 2>&1))"
if command -v aria2c >/dev/null 2>&1; then
  echo "Using aria2c segmented downloads with $ARIA_CONNECTIONS connections"
else
  echo 'aria2c not found; using resumable curl fallback'
  echo 'For faster downloads: sudo apt install -y aria2'
fi

"$PYTHON_BIN" - <<'PY'
import importlib.util
need=['numpy','scipy','astropy']
miss=[m for m in need if importlib.util.find_spec(m) is None]
if miss:
    raise SystemExit('Missing Python packages: '+', '.join(miss)+'. Create/activate a venv and install them there.')
PY

if [[ ! -f "$TEMPLATE_CSV" ]]; then
  echo "Template not found at $TEMPLATE_CSV"
  echo 'ERROR: local EZmock runs now expect the fixed validated Phase-7 template.' >&2
  echo 'Place phase7_templates.csv under phase7_template/; do not rebuild CLASS for this transport run.' >&2
  exit 3
fi

# Build the compact angular cache once if the user has preserved the validated
# mock2 random pair.  The source FITS remain untouched.
if [[ ! -s "$SHARED_RANDOM_CACHE" && \
      -s "$SHARED_RANDOM_DIR/mock_NGC.ran.fits" && \
      -s "$SHARED_RANDOM_DIR/mock_SGC.ran.fits" ]]; then
  echo "Building compact shared angular random cache: $SHARED_RANDOM_CACHE"
  "$PYTHON_BIN" code/desi_phase7_build_shared_random_cache.py \
    --ngc "$SHARED_RANDOM_DIR/mock_NGC.ran.fits" \
    --sgc "$SHARED_RANDOM_DIR/mock_SGC.ran.fits" \
    --out "$SHARED_RANDOM_CACHE" \
    --ngc-count "$SHARED_RANDOM_NGC_COUNT" \
    --sgc-count "$SHARED_RANDOM_SGC_COUNT" \
    --seed 20260914
fi

if [[ -s "$SHARED_RANDOM_CACHE" ]]; then
  echo "Using shared angular random cache: $SHARED_RANDOM_CACHE"
  echo 'Random redshifts will be redrawn separately for each mock within narrow-z/NGC-SGC cells.'
else
  echo 'Shared random cache not available; realization-specific random FITS will be downloaded.'
fi

fetch_one () {
  local url="$1" out="$2" part="${2}.part"
  local attempt rc have dir name

  dir=$(dirname "$out")
  name=$(basename "$out")
  mkdir -p "$dir"

  # Files produced by the current launcher are complete once renamed from .part.
  if [[ -s "$out" && ! -e "$part" ]]; then
    echo "REUSE existing download: $out ($(stat -c '%s' "$out") bytes)"
    return 0
  fi

  if [[ -e "$out" && ! -s "$out" ]]; then
    rm -f "$out"
  fi

  if command -v aria2c >/dev/null 2>&1; then
    echo "DOWNLOAD segmented: $url"
    aria2c \
      --continue=true \
      --max-connection-per-server="$ARIA_CONNECTIONS" \
      --split="$ARIA_CONNECTIONS" \
      --min-split-size=1M \
      --file-allocation=none \
      --max-tries=0 \
      --retry-wait=5 \
      --connect-timeout=60 \
      --timeout=60 \
      --summary-interval=10 \
      --console-log-level=notice \
      --dir="$dir" \
      --out="${name}.part" \
      "$url"
    test -s "$part"
    mv "$part" "$out"
    rm -f "${part}.aria2"
    echo "DOWNLOAD COMPLETE: $out ($(stat -c '%s' "$out") bytes)"
    return 0
  fi

  for attempt in $(seq 1 "$DOWNLOAD_ATTEMPTS"); do
    have=0
    if [[ -f "$part" ]]; then
      have=$(stat -c '%s' "$part" 2>/dev/null || echo 0)
    fi
    echo "DOWNLOAD attempt $attempt/$DOWNLOAD_ATTEMPTS: $url"
    echo "  existing bytes: $have"

    set +e
    if [[ "$have" -gt 0 ]]; then
      curl --http1.1 -fL --connect-timeout 60 --max-time 900 --continue-at - -o "$part" "$url"
      rc=$?
    else
      curl --http1.1 -fL --connect-timeout 60 --max-time 900 -o "$part" "$url"
      rc=$?
    fi
    set -e

    if [[ "$rc" -eq 0 && -s "$part" ]]; then
      mv "$part" "$out"
      echo "DOWNLOAD COMPLETE: $out ($(stat -c '%s' "$out") bytes)"
      return 0
    fi

    have=0
    if [[ -f "$part" ]]; then
      have=$(stat -c '%s' "$part" 2>/dev/null || echo 0)
    fi
    echo "  interrupted (curl rc=$rc), retained $have bytes; reconnecting..."
    if [[ "$rc" -eq 33 ]]; then
      echo 'ERROR: server refused HTTP range resume (curl rc=33).' >&2
      echo "Partial file retained at: $part" >&2
      return 33
    fi
    sleep 5
  done

  echo "ERROR: download did not complete after $DOWNLOAD_ATTEMPTS resumable attempts: $url" >&2
  echo "Partial file retained at: $part" >&2
  return 18
}

for M in $(seq "$FIRST" "$LAST"); do
  ROOT="$ROOT_BASE/mock${M}"
  WORK="$OUTROOT/work_mock${M}"
  OUT="$OUTROOT/mock_${M}"
  mkdir -p "$WORK" "$OUT"

  if [[ "${FORCE:-0}" != '1' && -s "$OUT/mock_$(printf '%02d' "$M")_summary.json" && -s "$OUT/mock_$(printf '%02d' "$M")_vector.csv" ]]; then
    echo "SKIP completed mock $M -> $OUT"
    continue
  fi

  echo "==== EZmock $M ===="
  fetch_one "$ROOT/BGS_ffa_NGC_clustering.dat.fits" "$WORK/mock_NGC.dat.fits"
  fetch_one "$ROOT/BGS_ffa_SGC_clustering.dat.fits" "$WORK/mock_SGC.dat.fits"

  RANDOM_ARGS=()
  VALIDATE_RANDOMS=0
  if [[ -s "$SHARED_RANDOM_CACHE" ]]; then
    RANDOM_ARGS=(--shared-random-cache "$SHARED_RANDOM_CACHE")
  else
    fetch_one "$ROOT/BGS_ffa_NGC_0_clustering.ran.fits" "$WORK/mock_NGC.ran.fits"
    fetch_one "$ROOT/BGS_ffa_SGC_0_clustering.ran.fits" "$WORK/mock_SGC.ran.fits"
    RANDOM_ARGS=(--random "$WORK/mock_NGC.ran.fits" "$WORK/mock_SGC.ran.fits")
    VALIDATE_RANDOMS=1
  fi

  "$PYTHON_BIN" - "$WORK" "$VALIDATE_RANDOMS" <<'PY'
from astropy.io import fits
from pathlib import Path
import sys
w=Path(sys.argv[1]); vr=int(sys.argv[2])
files=[w/'mock_NGC.dat.fits',w/'mock_SGC.dat.fits']
if vr:
    files += [w/'mock_NGC.ran.fits',w/'mock_SGC.ran.fits']
for p in files:
    with fits.open(p,memmap=True) as h:
        h.verify('exception')
        if len(h)<2 or h[1].data is None or len(h[1].data)==0:
            raise RuntimeError(f'invalid FITS: {p}')
        print('FITS_OK',p,len(h[1].data))
PY

  "$PYTHON_BIN" code/desi_phase7_ezmock_placebo_realization.py \
    --data "$WORK/mock_NGC.dat.fits" "$WORK/mock_SGC.dat.fits" \
    "${RANDOM_ARGS[@]}" \
    --templates "$TEMPLATE_CSV" --outdir "$OUT" --mock-id "$M" \
    --zmin 0.10 --zmax 0.40 --dz-proxy 0.02 --ntracer 5 \
    --analysis-z-edges 0.10,0.20,0.30,0.40 \
    --sep-edges 20,40,60,80,100,120,140 \
    --random-factor "$RANDOM_FACTOR" --neighbors-per-anchor 48 --theta-min-deg 0.05 \
    --seed "$SEED"

  if [[ "$KEEP_FITS" != '1' ]]; then
    rm -rf "$WORK"
  fi
  echo "DONE mock $M -> $OUT"
done

N=$(find "$OUTROOT" -type f -name 'mock_*_vector.csv' | wc -l | tr -d ' ')
echo "Completed vectors currently available: $N"

if [[ "$N" -ge 40 ]]; then
  echo 'At least 40 realizations available; running aggregate.'
  "$PYTHON_BIN" code/desi_phase7_ezmock_aggregate.py \
    --mock-root "$OUTROOT" \
    --real-vector source_data/wake_phase7_multitracer_real_vector.csv \
    --outdir "$OUTROOT/aggregate"
  echo "Aggregate: $OUTROOT/aggregate/summary_ezmock_placebo_covariance.json"
else
  echo 'Aggregate requires >=40 completed realizations. Run additional mock ranges when ready.'
fi
