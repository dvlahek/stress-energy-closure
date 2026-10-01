# E63 source preflight — public OuterRim lightcone truth source

E62 passed direct N-body velocity truth in a periodic Quijote z=1 box. Before fixing a survey/lightcone transfer test, E63 first verifies that a public lightcone source exposes all quantities needed for a truth-labelled cut-sky reconstruction.

The chosen source is CosmoDC2 Mock V1 through IRSA TAP. CosmoDC2 is based on the Outer Rim N-body simulation, covers 440 deg2, and the public table exposes true sky coordinates/redshift together with 3D positions and velocity components.

This stage is deliberately bounded. It queries only a 2-degree-radius region centered at RA=55 deg, Dec=-41 deg and 0.9<=z_true<1.0, with TOP 4097 rows. The patch is not the final E63 scientific geometry. The query only certifies schema and access.

A PASS authorizes a second preregistration containing the actual tracer definition, survey/lightcone geometry, radial selection and truth gates. No observed eBOSS catalogue or observed odd vector is accessed.
