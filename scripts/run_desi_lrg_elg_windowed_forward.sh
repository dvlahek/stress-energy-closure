#!/usr/bin/env bash
set -euo pipefail

RR="${1:-lrg_elg_rr_window_r4/lrg_elg_zresolved_rr_window_counts.npz}"
OUT="${2:-lrg_elg_windowed_forward_r4}"

python code/build_lrg_elg_zresolved_windowed_forward.py \
  --rr-counts "${RR}" \
  --prewindow-template source_data/lrg_elg_r4_inference_inputs/linked_standard_and_wake.csv \
  --linked-summary source_data/lrg_elg_r4_inference_inputs/linked_standard_summary.json \
  --outdir "${OUT}" \
  --mass 0.06 \
  --z-match 1100 \
  --frac 0.30
