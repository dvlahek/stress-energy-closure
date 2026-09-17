#!/usr/bin/env bash
set -euo pipefail

# End-to-end Linux runner for the DESI DR1 full-sample five-tracer
# robustness check with 256 mark permutations.
#
# It reproduces the frozen primary configuration and changes only
# --permutations 32 -> 256. Outputs go to a separate directory.
#
# Usage from the repository root:
#   bash scripts/run_desi_perm256_full_linux.sh
#
# Optional environment variables:
#   PYTHON=python3
#   OUTDIR=phase7_multitracer_fullsample_perm256
#   WORKDIR=.desi_perm256_work

PYTHON="${PYTHON:-python3}"
OUTDIR="${OUTDIR:-phase7_multitracer_fullsample_perm256}"
WORKDIR="${WORKDIR:-.desi_perm256_work}"
CLASS_COMMIT="e85808324f51fc694d12e3ed7439552a3c3f9540"

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
mkdir -p "$WORKDIR"

if ! command -v git >/dev/null 2>&1; then
  echo "git is required" >&2
  exit 1
fi
if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required" >&2
  exit 1
fi
if ! command -v "$PYTHON" >/dev/null 2>&1; then
  echo "$PYTHON is required" >&2
  exit 1
fi

# Install Python dependencies into the active environment.
"$PYTHON" -m pip install --upgrade pip setuptools wheel
"$PYTHON" -m pip install numpy scipy astropy cython

# Fetch and install the exact CLASS revision if it is not already staged.
CLASS_DIR="$WORKDIR/class_public"
if [[ ! -d "$CLASS_DIR/.git" ]]; then
  git init "$CLASS_DIR"
  git -C "$CLASS_DIR" remote add origin https://github.com/lesgourg/class_public.git
  git -C "$CLASS_DIR" fetch --depth 1 origin "$CLASS_COMMIT"
  git -C "$CLASS_DIR" checkout --detach FETCH_HEAD
fi
ACTUAL_CLASS_COMMIT="$(git -C "$CLASS_DIR" rev-parse HEAD)"
if [[ "$ACTUAL_CLASS_COMMIT" != "$CLASS_COMMIT" ]]; then
  echo "Unexpected CLASS commit: $ACTUAL_CLASS_COMMIT" >&2
  exit 1
fi
"$PYTHON" -m pip install "$CLASS_DIR"

# Build the fixed hidden-state and Doppler templates used by the analysis.
TEMPLATE_DIR="$WORKDIR/phase7_template"
mkdir -p "$TEMPLATE_DIR"
"$PYTHON" code/wake_phase7_template.py \
  --outdir "$TEMPLATE_DIR" --mass 0.06 --z-match 1100 --frac 0.30 \
  --s-min 20 --s-max 140 --s-step 2

# Download the exact public DESI DR1 BGS files used by the frozen workflow.
PRIMARY='https://data.desi.lbl.gov/public/dr1'
AWS='https://desidata.s3.amazonaws.com/dr1'
DATA_ROOT='survey/catalogs/dr1/LSS/iron/LSScats/v1.5'

fetch_one() {
  local rel="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "Using existing $out"
    return
  fi
  rm -f "$out"
  if curl -fL --retry 1 --retry-delay 3 --connect-timeout 15 --max-time 300 \
      -o "$out" "$PRIMARY/$rel"; then
    echo "Downloaded $rel from DESI primary"
  else
    echo "DESI primary unavailable for $rel; using official DESI AWS mirror"
    rm -f "$out"
    curl -fL --retry 4 --retry-all-errors --retry-delay 5 \
      --connect-timeout 20 --max-time 1800 \
      -o "$out" "$AWS/$rel"
  fi
  test -s "$out"
}

DATA_NGC="$WORKDIR/BGS_BRIGHT-21.5_NGC_clustering.dat.fits"
DATA_SGC="$WORKDIR/BGS_BRIGHT-21.5_SGC_clustering.dat.fits"
RAND_NGC="$WORKDIR/BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits"
RAND_SGC="$WORKDIR/BGS_BRIGHT-21.5_SGC_0_clustering.ran.fits"

fetch_one "$DATA_ROOT/BGS_BRIGHT-21.5_NGC_clustering.dat.fits" "$DATA_NGC"
fetch_one "$DATA_ROOT/BGS_BRIGHT-21.5_SGC_clustering.dat.fits" "$DATA_SGC"
fetch_one "$DATA_ROOT/BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits" "$RAND_NGC"
fetch_one "$DATA_ROOT/BGS_BRIGHT-21.5_SGC_0_clustering.ran.fits" "$RAND_SGC"

"$PYTHON" - <<PY
from astropy.io import fits
from pathlib import Path
files = [
    Path(r"$DATA_NGC"), Path(r"$DATA_SGC"),
    Path(r"$RAND_NGC"), Path(r"$RAND_SGC")
]
for p in files:
    with fits.open(p, memmap=True) as h:
        h.verify('exception')
        assert len(h) > 1 and h[1].data is not None
        print('FITS_OK', p, len(h[1].data))
PY

# Run the frozen full-sample analysis with only the permutation count changed.
"$PYTHON" code/desi_dr1_phase7_multitracer_fullsample.py \
  --data "$DATA_NGC" "$DATA_SGC" \
  --random "$RAND_NGC" "$RAND_SGC" \
  --templates "$TEMPLATE_DIR/phase7_templates.csv" \
  --outdir "$OUTDIR" \
  --zmin 0.10 --zmax 0.40 --dz-proxy 0.02 --ntracer 5 \
  --analysis-z-edges 0.10,0.20,0.30,0.40 \
  --sep-edges 20,40,60,80,100,120,140 \
  --random-factor 2.0 --neighbors-per-anchor 48 \
  --theta-min-deg 0.05 --jackknife 30 --permutations 256 \
  --seed 20260913

# Record exact provenance next to the robustness output.
git rev-parse HEAD > "$OUTDIR/REPOSITORY_COMMIT.txt"
printf '%s\n' "$CLASS_COMMIT" > "$OUTDIR/CLASS_COMMIT.txt"
cp "$TEMPLATE_DIR/template_summary.json" "$OUTDIR/template_summary.json"

printf '\nDONE\nRepository: %s\nOutput: %s/%s\n' "$ROOT" "$ROOT" "$OUTDIR"
