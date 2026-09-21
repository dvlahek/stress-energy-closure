#!/usr/bin/env bash
set -euo pipefail

# Bootstrap the Python environment for DESI DR1 LRGxELG production on Ubuntu/WSL.
# Run from the repository root:
#   bash scripts/bootstrap_ubuntu_desi.sh

if ! command -v sudo >/dev/null 2>&1; then
  echo "sudo is required for Ubuntu system packages." >&2
  exit 2
fi

sudo apt-get update
sudo apt-get install -y \
  ca-certificates curl git \
  build-essential gcc g++ make pkg-config \
  libgsl-dev libgomp1 \
  python3 python3-dev python3-venv \
  tmux

python3 -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip setuptools wheel cython
python -m pip install -r requirements-desi.txt

python - <<'PY'
import importlib.metadata as md
import numpy as np
import scipy
import astropy
import Corrfunc
import pycorr
import cosmoprimo
from cosmoprimo.fiducial import TabulatedDESI

cosmo = TabulatedDESI()
r1 = float(np.asarray(cosmo.comoving_radial_distance(1.0)))
print("DESI environment smoke test: PASS")
for package in ["numpy", "scipy", "astropy", "pycorr", "Corrfunc", "cosmoprimo"]:
    try:
        print(f"{package}: {md.version(package)}")
    except Exception:
        print(f"{package}: imported")
print(f"TabulatedDESI chi(z=1) = {r1:.6f} Mpc/h")
PY

python code/desi_dr1_lrg_elg_exact_zresolved.py --help >/dev/null
echo "Estimator CLI smoke test: PASS"
echo "Activate later with: source .venv/bin/activate"
