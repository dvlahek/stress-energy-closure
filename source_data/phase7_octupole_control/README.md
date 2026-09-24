# DESI DR1 octupole control

This directory contains the measured luminosity-rank octupole vector, jackknife covariance and numerical summary for the 32-permutation higher-multipole consistency test.

The octupole has empirical permutation $p=0.81818$ and maximum single-bin $|z|=1.3221$. No physical wake template is fitted to the octupole.

The dipole output from the same run agrees exactly with `reference_dipole_32perm.csv`. We retain this reference because the primary dipole measurement was subsequently recalibrated using 256 permutations. The equality check applies to the recorded 32-permutation realization, not to the later null-corrected vector.

The directory also includes the measured octupole data vector, octupole covariance, run commit, seed and file checksums. See `docs/OBSERVABLE_OCTUPOLE_CONTROL.md` for the estimator definition and covariance convention.
