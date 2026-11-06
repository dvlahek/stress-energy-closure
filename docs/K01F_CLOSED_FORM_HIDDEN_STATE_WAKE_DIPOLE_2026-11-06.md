# K01F — Closed-form hidden-state wake dipole

**Date:** 2026-11-06  
**Branch:** `experiment/k01-matched-f-detectability-20261105`  
**Scope:** analytic/static single-halo BKT wake physics; no observations, survey covariance, or self-gravitating simulation.

## 1. Starting point

For a halo moving with constant velocity in a static background, the generalized late-time halo-frame wake for an arbitrary isotropic collisionless distribution (F(q)) can be written, up to the common overall normalization,

[
deltaho_{
u,F}^{H}(r,mu)
=
rac{mathcal C}{r}
W_F(mu;alpha),
]

where

[
W_F(mu;alpha)
=
int_0^infty ds,
F!left(sqrt{s^2+alpha^2+2alpha smu}ight),
]

[
alpha = rac{m_
u a v_H}{c,T_{
u,0}},
qquad
mu=hat{mathbf r}cdothat{mathbf v}_H.
]

This follows directly from the static constant-velocity limit of the moving-halo Vlasov/BKT solution.

## 2. Exact dipole reduction

Define the Legendre dipole coefficient

[
W_{1,F}(alpha)
=
rac{3}{2}
int_{-1}^{1}dmu,
mu,W_F(mu;alpha).
]

For fixed (s), use

[
q^2=s^2+alpha^2+2alpha smu,
]

so that

[
mu=rac{q^2-s^2-alpha^2}{2alpha s},
qquad
dmu=rac{q,dq}{alpha s}.
]

Swapping the (s) and (q) integrations yields a remarkable simplification. For (q>alpha), the inner (s)-integral cancels exactly. For (q<alpha), it equals (-2q^2/alpha^2) before the Legendre prefactor.

Therefore

[
oxed{
W_{1,F}(alpha)
=
-rac{3}{alpha^2}
int_0^alpha dq,q^2F(q).
}
]

Hence

[
oxed{
deltaho_{
u,1}^{H,F}(r)
=
-rac{3mathcal C}{ralpha^2}
int_0^alpha dq,q^2F(q).
}
]

The wake dipole is therefore controlled only by particles with comoving momentum below the halo-speed threshold (q<alpha).

This is the collisionless analogue of the familiar result that the dynamical-friction contribution is controlled by the population slower than the perturber.

## 3. Hidden-state theorem for the frozen pair

For the frozen source pair, define

[
Delta F(q)=F_-(q)-F_+(q).
]

The hidden dipole difference is exactly

[
oxed{
Delta W_1(alpha)
=
-rac{3}{alpha^2}
int_0^alpha dq,q^2Delta F(q).
}
]

Define the positive low-momentum hidden-occupation ratio

[
oxed{
mathcal H(alpha)
=
rac{
int_0^alpha dq,q^2[F_-(q)-F_+(q)]
}{
int_0^alpha dq,q^2F_{FD}(q)
}.
}
]

Since all three dipoles have the same overall negative sign in the adopted (mu=hat rcdothat v_H) convention,

[
mathcal H(alpha)
=
rac{|W_{1,-}|-|W_{1,+}|}{|W_{1,FD}|}.
]

Thus the static wake dipole distinguishes the hidden states if and only if their cumulative low-(q) occupations differ below the halo-speed threshold.

Finite matching of total number, energy and pressure moments does not constrain this truncated cumulative functional.

## 4. Frozen E8 prediction at z=0.95

For

[
m_
u=0.06 {m eV},
qquad z=0.95,
]

the threshold is

[
alpha(v_H)
=
rac{m_
u v_H}
{cT_{
u,0}(1+z)}.
]

The frozen E8 distributions give:

| (v_H) [km/s] | (alpha) | (mathcal H(alpha)) |
|---:|---:|---:|
| 50  | 0.03051 | 0.22964 |
| 100 | 0.06102 | 0.25988 |
| 200 | 0.12205 | 0.32493 |
| 300 | 0.18307 | 0.39245 |
| 500 | 0.30512 | 0.51015 |
| 800 | 0.48819 | 0.53795 |
| 1000 | 0.61024 | 0.41173 |

For the existing 200 km/s benchmark,

[
oxed{
mathcal H(200,{m km/s})=0.32493.
}
]

Therefore, in the static point-halo BKT wake,

[
oxed{
|W_{1,-}|-|W_{1,+}|
=
0.32493,|W_{1,FD}|.
}
]

This is a direct wake-dipole prediction, not a resonant source proxy and not a scattering proxy.

## 5. Relation to earlier K01 anchors

At the same (z=0.95), (v=200) km/s benchmark:

- E8 pointwise resonant occupation contrast:
  [
  (F_+-F_-)/F_{FD}=-0.37040;
  ]

- K01F exact static wake-dipole magnitude contrast:
  [
  (|W_{1,-}|-|W_{1,+}|)/|W_{1,FD}|=0.32493;
  ]

- K01C deterministic local scattering proxy:
  [
  (D_--D_+)/|D_{FD}|=0.43399.
  ]

These are three different observables and must not be identified numerically.

Their common physical ordering is

[
F_- 	ext{ has greater relevant low-}q	ext{ occupation than }F_+,
]

and correspondingly the (F_-) state produces the larger wake/drag response in all three constructions.

This sign/order consistency is the physically important cross-check.

## 6. Velocity-selective tomography

The threshold

[
q<alpha(v_H)
]

means that halo velocity acts as a cumulative phase-space filter.

Changing (v_H) changes the upper integration limit in

[
int_0^{alpha(v_H)}q^2F(q),dq.
]

Therefore a population of halos with different peculiar velocities probes different cumulative sectors of the hidden distribution.

For the frozen E8 deformation the hidden wake response rises from roughly 23% at 50 km/s to above 50% around 500–800 km/s before declining again.

A simple scan of the frozen pair places the maximum of (mathcal H) near

[
alphasimeq0.42,
qquad
v_Hsimeq 6.9	imes10^2 {m km/s},
]

with

[
mathcal H_{max}simeq0.56.
]

This suggests a new interpretation: velocity-conditioned wakes provide a form of **phase-space tomography** of collisionless hidden structure.

## 7. Main physical statement

The strongest K01 result is now:

> Matching a finite set of stress-energy moments does not fix the velocity-oriented gravitational wake of a collisionless medium. In the static moving-halo limit, the wake dipole is exactly proportional to the cumulative phase-space occupation below the halo-speed threshold.

Mathematically,

[
oxed{
mathcal M_i[F_+]=mathcal M_i[F_-]
;
otRightarrow;
int_0^alpha q^2F_+(q)dq
=
int_0^alpha q^2F_-(q)dq.
}
]

The moving gravitational perturbation therefore resolves information that is invisible to finite moment closure.

## 8. Interpretation boundary

This is not a claim that the complete stress-energy tensor is identical for the two states at all epochs and scales.

The frozen construction matches selected source moments at the source epoch.

K01F is also not yet a survey observable. It is an analytic wake-level observable in the external-potential/BKT approximation.

## 9. Next physics stage

**K01G:** generalize the exact static dipole result to an expanding universe and a realistic halo-velocity history.

The key question is no longer whether the hidden wake exists. K01F answers that analytically.

The next question is how much of the static hidden dipole survives after the retarded time kernel, cosmological expansion, and the distribution of halo velocities are included.

Only after K01G should a self-gravitating simulation be used as a validation layer.
