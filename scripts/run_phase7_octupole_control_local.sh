#!/usr/bin/env bash
set -euo pipefail

# Reproduce the DESI DR1 luminosity-rank octupole control.
# The reference dipole is the original 32-permutation realization used by this control.
#
# Usage:
#   bash scripts/run_phase7_octupole_control_local.sh
#
# Optional:
#   OUTDIR=phase7_octupole_control SEED=20260913 bash scripts/run_phase7_octupole_control_local.sh

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

OUTDIR="${OUTDIR:-phase7_octupole_control}"
SEED="${SEED:-20260913}"
DATA_DIR="${DATA_DIR:-phase7_desi_catalogs}"
mkdir -p "$OUTDIR" "$DATA_DIR"

PRIMARY='https://data.desi.lbl.gov/public/dr1'
AWS='https://desidata.s3.amazonaws.com/dr1'
RELROOT='survey/catalogs/dr1/LSS/iron/LSScats/v1.5'

fetch_one () {
  local rel="$1" out="$2"
  if [[ -s "$out" ]]; then
    echo "REUSE $out"
    return 0
  fi
  rm -f "${out}.part"
  echo "DOWNLOAD $rel"
  if ! curl -fL --retry 3 --retry-all-errors --retry-delay 5 \
      --connect-timeout 30 --max-time 3600 -o "${out}.part" "$PRIMARY/$rel"; then
    echo "Primary DESI endpoint unavailable; trying AWS mirror"
    rm -f "${out}.part"
    curl -fL --retry 6 --retry-all-errors --retry-delay 5 \
      --connect-timeout 30 --max-time 3600 -o "${out}.part" "$AWS/$rel"
  fi
  test -s "${out}.part"
  mv "${out}.part" "$out"
}

DNGC="$DATA_DIR/BGS_BRIGHT-21.5_NGC_clustering.dat.fits"
DSGC="$DATA_DIR/BGS_BRIGHT-21.5_SGC_clustering.dat.fits"
RNGC="$DATA_DIR/BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits"
RSGC="$DATA_DIR/BGS_BRIGHT-21.5_SGC_0_clustering.ran.fits"

fetch_one "$RELROOT/BGS_BRIGHT-21.5_NGC_clustering.dat.fits" "$DNGC"
fetch_one "$RELROOT/BGS_BRIGHT-21.5_SGC_clustering.dat.fits" "$DSGC"
fetch_one "$RELROOT/BGS_BRIGHT-21.5_NGC_0_clustering.ran.fits" "$RNGC"
fetch_one "$RELROOT/BGS_BRIGHT-21.5_SGC_0_clustering.ran.fits" "$RSGC"

python - "$DNGC" "$DSGC" "$RNGC" "$RSGC" <<'PY'
from astropy.io import fits
from pathlib import Path
import sys
for name in sys.argv[1:]:
    p=Path(name)
    with fits.open(p, memmap=True) as h:
        h.verify('exception')
        if len(h)<2 or h[1].data is None or len(h[1].data)==0:
            raise RuntimeError(f'invalid FITS: {p}')
        print('FITS_OK', p, len(h[1].data))
PY

python code/desi_dr1_phase7_octupole_control.py \
  --data "$DNGC" "$DSGC" \
  --random "$RNGC" "$RSGC" \
  --outdir "$OUTDIR" \
  --reference-dipole source_data/phase7_octupole_control/reference_dipole_32perm.csv \
  --zmin 0.10 --zmax 0.40 --dz-proxy 0.02 --ntracer 5 \
  --analysis-z-edges 0.10,0.20,0.30,0.40 \
  --sep-edges 20,40,60,80,100,120,140 \
  --random-factor 2.0 --neighbors-per-anchor 48 \
  --theta-min-deg 0.05 --jackknife 30 --permutations 32 \
  --seed "$SEED"

printf '%s\n' "$SEED" > "$OUTDIR/STOCHASTIC_SEED.txt"
git rev-parse HEAD > "$OUTDIR/REPOSITORY_COMMIT.txt"

echo
echo "DONE: $OUTDIR"
echo "Primary result:"
python - "$OUTDIR/summary_octupole_control.json" <<'PY'
import json,sys
s=json.load(open(sys.argv[1]))
o=s['octupole_null_test']; d=s['dipole_reproduction']['reference_check']
print('octupole_status =', o['status'])
print('octupole_empirical_p =', o['global_empirical_pvalue'])
print('octupole_mahalanobis =', o['global_mahalanobis'])
print('octupole_max_abs_single_bin_z =', o['max_abs_single_bin_z_diagnostic'])
print('dipole_reference_available =', d.get('available'))
if d.get('available'):
    print('dipole_reference_max_abs_difference =', d['max_abs_difference'])
    print('dipole_reference_correlation =', d['correlation'])
PY
