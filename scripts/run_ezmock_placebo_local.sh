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

FIRST="${1:-1}"
LAST="${2:-$FIRST}"
KEEP_FITS="${KEEP_FITS:-0}"
ROOT_BASE='https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/EZmock/bright/v1'
OUTROOT="${OUTROOT:-phase7_ezmock_local}"
TEMPLATE_DIR="${TEMPLATE_DIR:-phase7_template}"
TEMPLATE_CSV="$TEMPLATE_DIR/phase7_templates.csv"
SEED="${SEED:-20260913}"

mkdir -p "$OUTROOT" "$TEMPLATE_DIR"

if ! command -v python >/dev/null 2>&1; then
  echo 'ERROR: python not found' >&2; exit 2
fi
if ! command -v curl >/dev/null 2>&1; then
  echo 'ERROR: curl not found' >&2; exit 2
fi

python - <<'PY'
import importlib.util
need=['numpy','scipy','astropy']
miss=[m for m in need if importlib.util.find_spec(m) is None]
if miss:
    raise SystemExit('Missing Python packages: '+', '.join(miss)+'. Install with: pip install numpy scipy astropy')
PY

if [[ ! -f "$TEMPLATE_CSV" ]]; then
  echo "Template not found at $TEMPLATE_CSV"
  echo 'Building fixed Phase-7 template requires the pinned CLASS Python module.'
  if ! python - <<'PY'
import importlib.util,sys
sys.exit(0 if importlib.util.find_spec('classy') else 1)
PY
  then
    cat >&2 <<'EOF'
ERROR: CLASS/classy is not installed.
Install the pinned CLASS version once:
  rm -rf class_public
  git clone https://github.com/lesgourg/class_public.git class_public
  git -C class_public checkout e85808324f51fc694d12e3ed7439552a3c3f9540
  pip install ./class_public
Then rerun this script.
EOF
    exit 3
  fi
  python code/wake_phase7_template.py \
    --outdir "$TEMPLATE_DIR" --mass 0.06 --z-match 1100 --frac 0.30 \
    --s-min 20 --s-max 140 --s-step 2
fi

fetch_one () {
  local url="$1" out="$2"
  rm -f "$out"
  echo "DOWNLOAD $url"
  curl --http1.1 -fL --retry 8 --retry-all-errors --retry-delay 5 \
    --connect-timeout 60 --max-time 1800 -o "$out" "$url"
  test -s "$out"
}

for M in $(seq "$FIRST" "$LAST"); do
  ROOT="$ROOT_BASE/mock${M}"
  WORK="$OUTROOT/work_mock${M}"
  OUT="$OUTROOT/mock_${M}"
  mkdir -p "$WORK" "$OUT"

  echo "==== EZmock $M ===="
  fetch_one "$ROOT/BGS_ffa_NGC_clustering.dat.fits" "$WORK/mock_NGC.dat.fits"
  fetch_one "$ROOT/BGS_ffa_SGC_clustering.dat.fits" "$WORK/mock_SGC.dat.fits"
  fetch_one "$ROOT/BGS_ffa_NGC_0_clustering.ran.fits" "$WORK/mock_NGC.ran.fits"
  fetch_one "$ROOT/BGS_ffa_SGC_0_clustering.ran.fits" "$WORK/mock_SGC.ran.fits"

  python - "$WORK" <<'PY'
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

  python code/desi_phase7_ezmock_placebo_realization.py \
    --data "$WORK/mock_NGC.dat.fits" "$WORK/mock_SGC.dat.fits" \
    --random "$WORK/mock_NGC.ran.fits" "$WORK/mock_SGC.ran.fits" \
    --templates "$TEMPLATE_CSV" --outdir "$OUT" --mock-id "$M" \
    --zmin 0.10 --zmax 0.40 --dz-proxy 0.02 --ntracer 5 \
    --analysis-z-edges 0.10,0.20,0.30,0.40 \
    --sep-edges 20,40,60,80,100,120,140 \
    --random-factor 1.0 --neighbors-per-anchor 48 --theta-min-deg 0.05 \
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
  python code/desi_phase7_ezmock_aggregate.py \
    --mock-root "$OUTROOT" \
    --real-vector source_data/wake_phase7_multitracer_real_vector.csv \
    --outdir "$OUTROOT/aggregate"
  echo "Aggregate: $OUTROOT/aggregate/summary_ezmock_placebo_covariance.json"
else
  echo 'Aggregate requires >=40 completed realizations. Run additional mock ranges when ready.'
fi
