# CLASS response-optimization diagnostic

This diagnostic tests whether the small Planck-anchored tensor B-mode separation in the manuscript can be substantially increased by changing the stress-energy-matched radial relic distribution or the relic mass while preserving the same 30% pointwise deformation cap.

## Setup

- Planck 2018 reference cosmology, as in the manuscript forecast.
- Pinned CLASS commit `e85808324f51fc694d12e3ed7439552a3c3f9540`.
- One non-cold species with masses 0.03, 0.06, 0.10, 0.18, 0.30 and 0.60 eV.
- Smooth ten-function deformation basis around the CLASS Fermi-Dirac spectrum.
- Exact numerical null-space matching of particle density, energy density and pressure at z=1100.
- Final positive pairs satisfy a 30% pointwise fractional-deformation cap.
- Optimization target: ideal full-sky cosmic-variance-limited tensor B-mode S/N over 20 <= ell <= 300. Instrument noise, foregrounds, lensing residuals and cosmological-parameter degeneracies are intentionally omitted, so this is an optimistic upper-bound diagnostic.

## Full CLASS mass sweep

The largest actual signal among the tested masses occurs at 0.10 eV. The full CLASS final-pair calculation gives

- ideal full-sky cosmic-variance S/N = 0.0389;
- maximum symmetric B-mode separation = 0.0274% near ell = 232--233;
- separation at ell = 80 = about -0.0021%;
- maximum n/rho/P matching error = 4.8e-16.

All six tested masses remain below S/N = 0.04. The mass-sweep values are stored in `source_data/CLASS_response_optimization_mass_sweep.csv`.

The small-probe Jacobian strongly over-predicts the achievable final response (for the 0.10 eV case it predicts S/N about 6.5 while the confirmed final-pair value is about 0.039). The Jacobian-guided construction should therefore be interpreted only as a search heuristic, not as a linear-response upper bound. The mass sweep is also not a proof that no other positive stress-energy-matched distribution can produce a larger signal.

## Resolution convergence of the best tested case

For the 0.10 eV optimized pair, independent CLASS resolution levels give:

| calculation | ideal CV S/N | max separation |
|---|---:|---:|
| default | 0.038974 | 0.0274245% |
| resolution 1 | 0.038935 | 0.0274125% |
| resolution 2 | 0.038943 | 0.0274129% |

The maximum pointwise change between resolution levels 1 and 2 over 20 <= ell <= 300 is 1.77e-5 percentage points. The optimized signal is therefore numerically stable under the tested hierarchy and momentum-resolution refinements.

## Manuscript consequence

This diagnostic does not support reframing the CLASS calculation as a detectability forecast. Direct response optimization and a semi-relativistic mass sweep increase the best stable separation only modestly, from the order-0.02% benchmark to order 0.027%. The manuscript should retain the cosmological calculation as a physical-realization control demonstrating that the source-degenerate response distinction survives in a standard Boltzmann observable, not as an observationally accessible prediction.
