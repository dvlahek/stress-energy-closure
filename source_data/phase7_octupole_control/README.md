# DESI DR1 octupole control

This directory contains the measured luminosity-rank octupole vector, jackknife covariance and numerical summary for the 32-permutation higher-multipole consistency test.

The octupole has empirical permutation $p=0.81818$ and maximum single-bin $|z|=1.3221$. No physical wake template is fitted to the octupole.

The dipole output from the same run agrees exactly with `reference_dipole_32perm.csv`. We retain this reference because the primary dipole measurement was subsequently recalibrated using 256 permutations. The equality check applies to the recorded 32-permutation realization, not to the later null-corrected vector.

The retained numerical and provenance files have a separate, verifiable `SHA256SUMS.txt`. Run `sha256sum -c SHA256SUMS.txt` from this directory to check all six entries. `SHA256SUMS_ORIGINAL_RUN.txt` preserves the original eight-entry checksum list, including four raw run products that are not in this compact release. The current runner regenerates those raw products; the published `summary_octupole_control_compact.json` is a separate compact summary and does not have the same output schema as `summary_octupole_control.json`.

The original-run checksums document the historical run. A rerun using current descriptive metadata need not reproduce the raw summary byte for byte. See `docs/OBSERVABLE_OCTUPOLE_CONTROL.md` for the estimator and covariance conventions.
