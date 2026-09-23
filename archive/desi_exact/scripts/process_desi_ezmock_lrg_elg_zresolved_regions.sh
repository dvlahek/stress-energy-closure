#!/usr/bin/env bash
set -euo pipefail

START=1
COUNT=40
NRANDOM=4
NTHREADS=16
INROOT="data/desi_dr1_ezmock_dark_v1"
OUTNGC="mocks_zresolved_desi_ngc"
OUTSGC="mocks_zresolved_desi_sgc"
KEEP_INPUTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --start) START="$2"; shift 2 ;;
    --count) COUNT="$2"; shift 2 ;;
    --nrandom) NRANDOM="$2"; shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --inroot) INROOT="$2"; shift 2 ;;
    --out-ngc) OUTNGC="$2"; shift 2 ;;
    --out-sgc) OUTSGC="$2"; shift 2 ;;
    --keep-inputs) KEEP_INPUTS=1; shift ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

END=$((START+COUNT-1))
for ((m=START; m<=END; m++)); do
  ngc=$(printf "%s/exact_z_%04d/lrg_elg_exact_zresolved_odd_multipoles.csv" "$OUTNGC" "$m")
  sgc=$(printf "%s/exact_z_%04d/lrg_elg_exact_zresolved_odd_multipoles.csv" "$OUTSGC" "$m")
  indir="${INROOT}/mock${m}"

  if [[ -s "$ngc" && -s "$sgc" ]]; then
    echo "SKIP mock${m}: NGC and SGC outputs exist"
    if (( KEEP_INPUTS == 0 )) && [[ -d "$indir" ]]; then rm -rf -- "$indir"; fi
    continue
  fi

  echo "===== REGIONAL DESI mock${m} ====="
  bash scripts/download_desi_ezmock_lrg_elg.sh \
    --start "$m" --count 1 --nrandom "$NRANDOM" --outroot "$INROOT"

  if [[ ! -s "$ngc" ]]; then
    bash scripts/run_desi_ezmock_lrg_elg_exact_zresolved.sh \
      --start "$m" --count 1 --nrandom "$NRANDOM" --region NGC \
      --inroot "$INROOT" --outroot "$OUTNGC" --nthreads "$NTHREADS"
  fi
  if [[ ! -s "$sgc" ]]; then
    bash scripts/run_desi_ezmock_lrg_elg_exact_zresolved.sh \
      --start "$m" --count 1 --nrandom "$NRANDOM" --region SGC \
      --inroot "$INROOT" --outroot "$OUTSGC" --nthreads "$NTHREADS"
  fi

  if [[ ! -s "$ngc" || ! -s "$sgc" ]]; then
    echo "ERROR mock${m}: regional output missing; preserving $indir" >&2
    exit 4
  fi
  if (( KEEP_INPUTS == 0 )); then rm -rf -- "$indir"; fi
  echo "DONE regional mock${m}"
done

echo "DONE regional NGC/SGC mocks ${START}..${END}, nrandom=${NRANDOM}"
