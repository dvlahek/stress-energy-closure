# Observational forecast provenance

This calculation is an observationally anchored forecast, not a fit to observational data.

## External inputs

### Planck 2018 reference cosmology

The background cosmological parameters are anchored to the Planck 2018 base-LambdaCDM TT,TE,EE+lowE+lensing central values:

- Planck Collaboration VI, *Planck 2018 results. VI. Cosmological parameters*, Astronomy & Astrophysics **641**, A6 (2020), DOI: 10.1051/0004-6361/201833910.
- Values used here are rounded to the published central values: H0=67.36 km s^-1 Mpc^-1, omega_b=0.02237, omega_cdm=0.1200, n_s=0.9649, ln(10^10 A_s)=3.044 and tau=0.054.
- The reference massive-neutrino convention follows the standard CLASS Planck-2018 setup: one 0.06 eV ncdm species, N_ur=2.0328 and T_ncdm/T_gamma=0.71611.

No Planck maps, likelihood files or measured CMB bandpowers are ingested by this workflow.

### Boltzmann solver

The observable calculation uses the public CLASS code:

- https://github.com/lesgourg/class_public
- pinned commit: `e85808324f51fc694d12e3ed7439552a3c3f9540`

CLASS documents that `use_ncdm_psd_files = 1` reads an arbitrary non-cold-dark-matter phase-space distribution from a file and that `tensor_method = exact` includes ncdm species in the tensor calculation. The workflow therefore changes only the homogeneous relic phase-space distribution; it does not modify CLASS perturbation equations.

### Motivation for nonthermal relic spectra

The test pair is synthetic and is not claimed to be a measured neutrino distribution. Nonthermal relic-neutrino spectra are nevertheless a legitimate cosmological possibility and have been studied explicitly, for example:

- G. Barenboim, J. Froustey, C. Pitrou and H. Sanchis, *Primordial neutrinos fade to gray: Constraints from cosmological observables*, Physical Review D **111**, 123549 (2025).
- G. Barenboim, H. Sanchis, W. H. Kinney and D. Rios, *Bound on thermal y distortion of the cosmic neutrino background*, Physical Review D **110**, 123535 (2024).
- K. Ala-Mattinen, M. Heikinheimo, K. Kainulainen and K. Tuominen, *Momentum distributions of cosmic relics: Improved analysis*, Physical Review D **105**, 123005 (2022), DOI: 10.1103/PhysRevD.105.123005.

## Constructed stress-energy-matched pair

The code starts from the standard Fermi-Dirac shape f_FD(q)=1/(exp(q)+1), with q=p/T_ncdm. Two smooth positive deformations F_+ and F_- are constructed so that the particle-density, energy-density and pressure moments are matched at z=1100 for m=0.06 eV. The deformation is restricted to a 30% pointwise departure from the Fermi-Dirac baseline in the reference run.

The pair is generated from a finite smooth basis. A null-space construction imposes the three moment constraints, and a response proxy is used only to select a well-separated direction inside that null space. CLASS then independently computes the full cosmological tensor response for each distribution.

The generated `pair_summary.json` records the numerical moment mismatch and all reference parameters. The generated PSD files and CLASS input files are included in the workflow artifact.

## Observable and interpretation

The workflow computes the primordial tensor B-mode spectrum with `modes=t`, `output=tCl,pCl`, `tensor_method=exact` and `r=0.01`. The value of r fixes the overall normalization; the reported relative difference between the two stress-energy-matched spectra is the relevant comparison.

The result should be described as a Planck-anchored CLASS forecast for a controlled nonthermal relic pair. It is not a detectability claim and is not a likelihood analysis of current CMB data.
