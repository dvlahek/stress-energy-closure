# K01H — Physically regulated phase-space unlocking window

**Date:** 2026-11-06  
**Branch:** `experiment/k01-matched-f-detectability-20261105`  
**Scope:** physics synthesis of K01F/K01G with the existing frozen DM14 NFW halo profiles. No observations, no new self-gravitating simulation.

## 1. The problem

K01G shows that the hidden-state response is moment-protected at small retarded phase

[
s sim k,v_{
u,m th},Delta t_{m memory}ll1,
]

and becomes progressively unlocked as (s) approaches order unity.

The earlier E28 finite-band recoil was dominated by (k>1,{m Mpc}^{-1}), but its total amplitude could not be certified because 94% of the force came from the broad (1<kle8,{m Mpc}^{-1}) band.

The key question is whether hidden-state sensitivity starts before the halo profile enters a genuinely UV-sensitive regime.

## 2. K01G unlocking scales for the original short history

For the same (z=1	o0.95) memory interval used in E27,

[
s_{max}=0.02351
]

at

[
k=0.05,h/{m Mpc}.
]

The exact frozen source kernel gives the first thresholds

[
s_{0.1%}simeq0.189,
qquad
s_{1%}simeq0.636,
qquad
s_{10%}simeq0.989.
]

Using (spropto k), this corresponds to

[
k_{0.1%}simeq0.40,h/{m Mpc},
]

[
oxed{k_{1%}simeq1.35,h/{m Mpc}},
]

[
oxed{k_{10%}simeq2.10,h/{m Mpc}}.
]

With (h=0.6736), the last two are approximately

[
k_{1%}simeq0.91,{m Mpc}^{-1},
qquad
k_{10%}simeq1.42,{m Mpc}^{-1}.
]

## 3. Existing frozen halo scales

The original DM14 conditional halo anchors at (z=0.95) are:

### Low-mass anchor

[
M_{200c}=1.4903	imes10^{12}M_odot,
]

[
c_{200c}=5.4399,
qquad
r_{200c}^{m phys}=0.1674,{m Mpc},
]

[
r_s^{m com}=0.05999,{m Mpc}.
]

### High-mass anchor

[
M_{200c}=1.4903	imes10^{13}M_odot,
]

[
c_{200c}=4.5634,
qquad
r_{200c}^{m phys}=0.3606,{m Mpc},
]

[
r_s^{m com}=0.15408,{m Mpc}.
]

## 4. Unlocking occurs before halo-profile suppression

At the 1% unlocking scale (ksimeq0.91,{m Mpc}^{-1}),

[
k r_ssimeq0.055
]

for the low-mass halo and

[
k r_ssimeq0.140
]

for the high-mass halo.

At the 10% unlocking scale (ksimeq1.42,{m Mpc}^{-1}),

[
k r_ssimeq0.085
]

for the low-mass halo and

[
k r_ssimeq0.219
]

for the high-mass halo.

Both remain well below unity.

The corresponding exact truncated-NFW Fourier amplitudes are:

| (k) [Mpc(^{-1})] | (u_{m low}(k)) | (u_{m high}(k)) |
|---:|---:|---:|
| 0.91 | 0.9958 | 0.9797 |
| 1.42 | 0.9898 | 0.9514 |
| 2.00 | 0.9798 | 0.9062 |
| 3.00 | 0.9553 | 0.8029 |
| 5.00 | 0.8815 | 0.5572 |
| 8.00 | 0.7289 | 0.2881 |

Thus the hidden-state response starts to become percent-level while the halo forcing profile is still nearly unsuppressed.

## 5. Physically preferred window

This suggests a new target window

[
oxed{
0.8lesssim klesssim2 {m Mpc}^{-1}
}
]

for the original short-memory benchmark.

This range has three desirable properties:

1. the hidden-state source contrast is beginning to unlock;
2. the NFW profile remains close to unity for the low-mass halo and only moderately suppressed for the high-mass halo;
3. the calculation avoids the very high-(k) tail that dominated the old E28 total-force integral.

This is a **mesoscopic wake window**, not the asymptotic UV.

## 6. Reinterpretation of E28

E28 integrated all the way to (8,{m Mpc}^{-1}), so its total force became dominated by scales whose physical reliability was difficult to certify.

K01H suggests that the scientifically cleaner target is not the all-(k) total drag.

The target should instead be a **band-limited, physically resolved wake observable** centered near the memory-unlocking transition.

This changes the question from

> What is the total neutrino drag after integrating to an arbitrary high-(k) cutoff?

to

> Does the hidden-state difference appear in a finite, halo-resolved (k)-window before UV sensitivity becomes dominant?

That question can be answered without assigning physical meaning to the full E28 force.

## 7. Main physical prediction

The project now predicts a scale-dependent transition:

### Large scales

[
kll k_{m unlock}:
qquad
F_+simeq F_-
	ext{ dynamically because the response is moment-protected.}
]

### Intermediate scales

[
ksim k_{m unlock}:
qquad
	ext{phase mixing begins to resolve hidden structure while the halo profile remains physical.}
]

### Deep small scales

[
k r_sgtrsim1:
qquad
	ext{halo structure and nonlinear dynamics require a genuinely resolved treatment.}
]

This creates a physically meaningful transition scale between moment closure and phase-space resolution.

## 8. Candidate core concept

A useful name for the mechanism is

[
oxed{	ext{gravitational phase-space unlocking}}
]

or more conservatively

[
oxed{	ext{retarded memory unlocking of collisionless hidden states}.}
]

The control parameter is

[
s(k,Delta t)
sim
k,v_{
u,m th},Delta t_{m memory}.
]

Finite moment closure controls the low-(s) expansion, while phase mixing exposes higher information in (F) as (s) grows.

## 9. Next decisive physics calculation

The next calculation should be a **finite-window wake contrast**, not a total-force integral.

Define

[
mathcal A_F[k_1,k_2]
=
int_{k_1}^{k_2}
dln k,
W(k),
mathcal R_F(k),
]

where (W(k)) contains the physical halo profile and the retarded Vlasov kernel.

The primary hidden observable should be

[
oxed{
R_{m window}
=
rac{
mathcal A_-[k_1,k_2]
-
mathcal A_+[k_1,k_2]
}{
|mathcal A_{FD}[k_1,k_2]|
}.
}
]

The recommended first window is

[
oxed{
kin[0.8,2.0],{m Mpc}^{-1}.
}
]

This calculation directly tests the K01G unlocking mechanism in a physically regulated halo-scale band.

It should replace further attempts to infer a total halo drag from the old (0.01)–(8,{m Mpc}^{-1}) E28 integral.
