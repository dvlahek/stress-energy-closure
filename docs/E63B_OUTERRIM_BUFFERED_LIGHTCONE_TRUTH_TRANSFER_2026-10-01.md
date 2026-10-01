# E63B — buffered OuterRim lightcone truth transfer

E63 source preflight passed: the public IRSA CosmoDC2 table exposes true sky/redshift, 3D position and velocity truth. E63B now performs the first preregistered lightcone truth test.

The test deliberately isolates **periodic snapshot -> lightcone geometry/evolution** from later survey-mask/HOD effects.

- source: CosmoDC2 / Outer Rim
- center: RA=55 deg, Dec=-41 deg
- source cone: radius 9.5 deg, 0.75<=z_true<1.16
- candidate retrieval floor: halo_mass >= 5e12 Msun
- science tracer: central galaxies only, top-mass-ranked to the exact E62 number density 165107/(1000 Mpc/h)^3
- probes: 512 deterministic selected tracers in radius 2.5 deg, 0.9<=z_true<1.0
- source/probe angular gap: 7 deg
- reconstruction: unchanged R=16 Mpc/h top-hat kernel, 256 Mpc/h primary cutoff, 192 Mpc/h check

The radial and angular buffers are chosen prospectively so every probe's full 256 Mpc/h support lies inside the source shell/cone. Therefore no survey random subtraction is used in E63B.

Primary gate is the **line-of-sight** velocity because this is the eventual conditioned tag:
Pearson(rec_LOS,true_LOS)>=0.70, sign agreement>=0.70, and 192/256 cutoff Pearson/sign>=0.90.

A PASS advances only to E64, where cut-sky/radial selection and ELG-like sampling must be introduced with velocity truth still retained. Observed eBOSS galaxies and observed odd remain sealed.
