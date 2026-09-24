# Numerical source data

This directory contains numerical results and validation products used in the manuscript. The result-to-code mapping is given in [REPRODUCIBILITY.md](../REPRODUCIBILITY.md). Earlier implementations and intermediate outputs are retained on the [research archive branch](https://github.com/dvlahek/stress-energy-closure/tree/archive/research-development-2026-09-24).

## Einstein–Vlasov response and forecasts

The `Fig1_`–`Fig5_` tables contain the numerical results for source-matched evolution, inverse reconstruction, mass dependence and kinetic hierarchy. `ExtendedData_direct_memory_convergence.csv` records the direct/memory comparison. The `CLASS_`, `hidden_channel_retention_direct.json`, `ksz_wake_screening_summary.json` and `wake_` products contain the corresponding transfer-level, forecast and consistency calculations.

## DESI DR1 luminosity-rank analysis

The reported five-tracer analysis is represented by:

- `phase7_desi_fullsample_perm256_summary.json` and `wake_phase7_multitracer_real_vector_perm256.csv`: DESI measurement and permutation calibration.
- `phase7_validation_manifest.json`: analysis definition and validation references.
- `phase7_gfinder_massproxy_summary.json`, `phase7_seed_control_summary.json`, `phase7_abacus_final_summary.json` and `phase7_abacus_injection_amplitude_sweep.json`: robustness and mock validation.
- `phase7_ezmock_local/` and `phase7_octupole_control/`: random-rank mock and higher-multipole controls.

The conservative fitted coefficient is $A_{\rm wake}=-0.0768409\pm0.0851660$, with empirical two-sided permutation $p=0.37354$.

## DESI DR1 LRG–ELG cross-correlation

The exact-pair analysis uses the products in `lrg_elg_r4_inference_inputs/` and `lrg_elg_windowed_forward_r4/`. The latter includes the window matrix, convolved templates and individual leave-one-out mock statistics.

- `lrg_elg_windowed_finite_mock_r4_120_final_2026-09-23.json`: final numerical inference.
- `lrg_elg_exact_final_manifest_2026-09-23.json`: analysis definition and result references.
- `lrg_elg_covariance_r4_120_checkpoint_2026-09-23.json`: covariance diagnostics.
- `lrg_elg_bruteforce_pair_closure_2026-09-21.json`, `lrg_elg_mubin_convergence_r4.json` and `lrg_elg_zresolved_random_density_audit_10pair.json`: estimator and numerical validation.

The window-convolved fit gives $A_{\rm wake}=0.00327466\pm0.00148162$ and nominal $Z=2.21019$. Empirical calibration with 120 mocks gives two-sided equivalents of $1.84\sigma$ and $1.96\sigma$ for the two specified statistics. The inference is interpreted as an approximately two-sigma indication, not as a detection.
