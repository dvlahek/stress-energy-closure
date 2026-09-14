#!/usr/bin/env bash
set -euo pipefail

# Local Linux runner for the Phase-7 EZmock five-tracer placebo covariance control.
# Usage:
#   bash scripts/run_ezmock_placebo_local.sh 1
#   bash scripts/run_ezmock_placebo_local.sh 1 8
#   KEEP_FITS=1 bash scripts/run_ezmock_placebo_local.sh 1 2
#
# The script processes one realization at a time and removes the large FITS
# inputs after each successful realization unless KEEP_FITS=1. Small output
# vectors/windows remain under phase7_ezmock_local/.
#
# DESI's public HTTP endpoint can close long transfers before completion.
# Downloads are therefore resumable. If aria2c is available it is preferred
# because segmented HTTP-range downloads are much faster on this endpoint;
# curl -C - remains the fallback.

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
  echo 'Building fixed Phase-7 template requires the pinned CLASS Python module.'
  if ! "$PYTHON_BIN" - <<'PY'
import importlib.util,sys
sys.exit(0 if importlib.util.find_spec('classy') else 1)
PY
  then
    cat >&2 <<'EOF'
ERROR: CLASS/classy is not installed in the active Python environment.
Install the pinned CLASS version once inside your virtual environment:
  rm -rf class_public
  git clone https://github.com/lesgourg/class_public.git class_public
  git -C class_public fetch origin e85808324f51fc694d12e3ed7439552a3c3f9540
  git -C class_public checkout e85808324f51fc694d12e3ed7439552a3c3f9540
  python -m pip install ./class_public
Then rerun this script.
EOF
    exit 3
  fi
  "$PYTHON_BIN" code/wake_phase7_template.py \
    --outdir "$TEMPLATE_DIR" --mass 0.06 --z-match 1100 --frac 0.30 \
    --s-min 20 --s-max 140 --s-step 2
fi

fetch_one () {
  local url="$1" out="$2" part="${2}.part"
  local attempt rc have dir name

  dir=$(dirname "$out")
  name=$(basename "$out")
  mkdir -p "$dir"

  # Preserve a partial file left by an older launcher.
  if [[ -s "$out" && ! -e "$part" ]]; then
    echo "Preserving existing partial download: $out -> $part"
    mv "$out" "$part"
  elif [[ -e "$out" && ! -s "$out" ]]; then
    rm -f "$out"
  fi

  if command -v aria2c >/dev/null 2>&1; then
    echo "DOWNLOAD segmented: $url"
    # aria2c continues an existing .part file and keeps its .aria2 control file
    # across interrupted runs.  Multiple HTTP ranges avoid the very slow
    # single-stream transfer observed from the DESI LBL endpoint.
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
      curl --http1.1 -fL \
        --connect-timeout 60 --max-time 900 \
        --continue-at - -o "$part" "$url"
      rc=$?
    else
      curl --http1.1 -fL \
        --connect-timeout 60 --max-time 900 \
        -o "$part" "$url"
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

  # A completed realization is never downloaded again unless FORCE=1.
  if [[ "${FORCE:-0}" != '1' && -s "$OUT/mock_$(printf '%02d' "$M")_summary.json" && -s "$OUT/mock_$(printf '%02d' "$M")_vector.csv" ]]; then
    echo "SKIP completed mock $M -> $OUT"
    continue
  fi

  echo "==== EZmock $M ===="
  fetch_one "$ROOT/BGS_ffa_NGC_clustering.dat.fits" "$WORK/mock_NGC.dat.fits"
  fetch_one "$ROOT/BGS_ffa_SGC_clustering.dat.fits" "$WORK/mock_SGC.dat.fits"
  fetch_one "$ROOT/BGS_ffa_NGC_0_clustering.ran.fits" "$WORK/mock_NGC.ran.fits"
  fetch_one "$ROOT/BGS_ffa_SGC_0_clustering.ran.fits" "$WORK/mock_SGC.ran.fits"

  "$PYTHON_BIN" - "$WORK" <<'PY'
from astropy.io import fits
from pathlib import Path
import sys
w=Path(sys.argv[1])
for p in [w/'mock_NGC.dat.fits',w/'mock_SGC.dat.fits',w/'mock_NGC.ran.fits',w/'mock_SGC.ran.fits']:
    with fits.open(p,memmap=True) as h:
        h.verify('exception')
        if len(h)<2 or h[1].data is None or len(h[1].data)==0:
            raise RuntimeError(f'invalid FITS: {p}')
        print('FITS_OK',p,len(h[1].data))
PY

  "$PYTHON_BIN" code/desi_phase7_ezmock_placebo_realization.py \
    --data "$WORK/mock_NGC.dat.fits" "$WORK/mock_SGC.dat.fits" \
    --random "$WORK/mock_NGC.ran.fits" "$WORK/mock_SGC.ran.fits" \
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
