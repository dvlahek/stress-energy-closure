# Blinded eBOSS LRG×ELG cross-LS random-only pilot

The published DR16 LRG and ELG clustering randoms independently sample
their tracer-specific angular/radial selection functions. A cross-LS
estimator uses four **oriented** pair terms and the tracer-specific
LRG×ELG random-pair window. This is not a new forced common angular cut.

This local pilot tests estimator algebra only: it creates disjoint
pseudo-D and pseudo-R subsamples from each of the *same parent random*
catalogues. Do not interpret pseudo-xi as measured galaxy odd multipoles,
physical mock-galaxy closure, a null p-value or a detection.

## Fixed inputs and supported test

The prospective protocol is
`source_data/eboss_dr16_random_only_cross_ls_protocol_2026-09-25.json`.
The first run uses SGC and mock 0001. For both released observed randoms
and the corresponding realistic mock randoms it selects z in [0.9,1.0),
applies the same validated four-factor weights, and uses one fixed
permutation to split 600 pseudo-D and 600 pseudo-R positions per tracer.
The script SHA-checks all four raw FITS against pinned prior provenance
and requires the original per-catalogue fine weighted n(z) JSON from
the successful local 36-comparison control. No observed or mock
galaxy FITS, odd measurements, wake templates or nuisance fits are read.

## Reproduce in the already validated WSL checkout

```bash
cd ~/stress-energy-closure
git pull --ff-only origin main
source .venv/bin/activate
set -o pipefail

python scripts/audit_eboss_dr16_random_only_cross_ls_pilot.py --self-test

mkdir -p eboss_workspace/local_pair_window
python scripts/audit_eboss_dr16_random_only_cross_ls_pilot.py \
  --ensemble-json eboss_workspace/local_nz/source/ninemock_window_ensemble.json \
  --cap SGC --mock-id 1 \
  --fine-checkpoint-dir eboss_workspace/local_nz \
  --observed-cache-dir eboss_workspace/local_rr/fits \
  --mock-cache-dir eboss_workspace/local_pair_window/mock_fits \
  --out-dir eboss_workspace/local_pair_window \
  2>&1 | tee eboss_workspace/local_pair_window/sgc_mock0001.log
```

Success is indicated by `EBOSS_RANDOM_ONLY_CROSS_LS_PILOT_OK SGC 1`,
with the machine-readable output
`eboss_workspace/local_pair_window/random_only_cross_ls_pilot.json`.
The diagnostics include four normalized cross pair terms, independent
reverse-order parity checks, and RR-supported xi mirror closure.
Unsupported RR cells are NaN, never silently filled for inference.

## Mandatory remaining gates

Separately pin and interpret the official LRG MANGLE sector/veto
semantics and ELG brickmask/image/extra-veto selection. Confirm each
released clustering random sample against its own published selection.
Do not reapply published vetoes blindly. Validate the full-density
cross-RR window, the pair-level selection and *mock-galaxy* cross
estimator closure with adequate mocks before calculating an observed
odd-sector data vector or a statistical significance.
