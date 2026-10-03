# E38 — one explicit certificate combining time-history uncertainty and finite halo size

E35 gave a low-k finite-size bound and E37 gave the missing time-history bound once `L_k` is known.

For

`A(k)=u_test(k) integral W(k,t) u_src(k,t) dt`

compare with `A_point=integral W dt`. Let `u_lin` interpolate the two endpoint source transforms. E37 bounds `|u_src-u_lin|` by `B_time(k,t)`. E35 gives endpoint/test bounds

`delta_i(k)=min(2,k^2<r^2>_i/6)`.

A sufficient absolute certificate is

`|A-A_point| <= integral |W| [B_time + delta_test + (1-s)delta0 + s delta1] dt`, with `s=t/T`.

The use of `|W|` is deliberate: accidental cancellations cannot manufacture a small error bar.

Toy QA at dimensionless `k=0.6` passed 205 checks. Temporal, spatial and combined bound components are 0.0906491, 0.0119008 and 0.10255.

For EinsteinVlasovNP this states exactly what E29 still lacks: not another endpoint metadata field, but a physically justified `L_k` (or stronger dynamical envelope) plus physically meaningful radial second-moment bounds.
