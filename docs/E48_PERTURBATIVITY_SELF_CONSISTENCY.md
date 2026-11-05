# E48 — self-consistency test of the E47 unconditioned linear channel

E47R1 is numerically clean, but its physical response is tiny. E48 asks a stricter question: **can the response coefficient required by E45 remain inside the linear response regime used to derive E47?**

The E25/E47 local term is

`delta Delta_a = d_a U`, with `U=v_rel,LOS/c`.

Define `D=max(|d_L|,|d_E|)`. The frozen E21 R16 LOS one-sigma velocity is about `135.862 km/s`, hence

`U_1sigma = 4.531871e-04`.

Even the extremely permissive condition `D U_1sigma <= 1` requires

`D <= 2206.6`,

and a 10% perturbative condition gives `D <= 220.7`.

## Optimized order-unity tracer geometry

For a diagnostic reference geometry `b_L=2`, `b_E=1`, `f=1`, E48 optimized over **all** directions `(d_L,d_E)` with `max(|d_L|,|d_E|)=D`. This is not a measured high-z eBOSS bias model.

For ~83% symmetric sign recovery:

- NGC requires `D >= 1.33e+05`, giving `D U_1sigma = 60.50`;
- SGC requires `D >= 1.58e+05`, giving `D U_1sigma = 71.44`.

The optimal direction in both caps is the maximally favorable anti-aligned response `d_L=-d_E`, so this is already the best case within that reference geometry.

For ~94% recovery the one-sigma local response becomes about `243.8` in NGC and `106.3` in SGC.

These values are not merely larger than a cautious linear regime; they are tens to hundreds.

## How extreme would the tracer geometry need to be to rescue perturbativity?

Take the most favorable anti-aligned response, set `f=1`, and allow the one-sigma local response to become as large as **unity**. To reach ~83% recovery still requires

- `b_L+b_E >= 250.1` in NGC,
- `b_L+b_E >= 299.8` in SGC.

If the linear term is required to stay below 10%, those sums rise to about `2511` and `3009`.

E48 therefore does **not** claim a mathematical no-go for arbitrary tracer functions. It establishes something narrower and more useful: the original unconditioned E21/E47 linear two-tracer route cannot reach the observed E45 practical-sensitivity scale with order-unity tracer geometry while remaining self-consistent as a linear local response.

That makes the E19 velocity-conditioned route the natural next experiment. E47 should be retained as a negative/limitation result rather than further optimized against the sealed observation.
