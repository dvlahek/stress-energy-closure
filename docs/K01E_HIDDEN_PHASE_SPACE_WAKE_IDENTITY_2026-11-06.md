# K01E — Hidden phase-space structure in the neutrino wake

**Date:** 2026-11-06  
**Branch:** `experiment/k01-matched-f-detectability-20261105`  
**Scope:** physics derivation only; no observed data, no survey covariance, no detection significance.

## 1. Physical question

Can two collisionless neutrino distributions that are matched in a finite set of macroscopic moments produce different gravitational wakes around the same moving halo?

The frozen E8 states satisfy, to machine precision,

[
mathcal M_i[F_+]=mathcal M_i[F_-]
]

for the selected number/energy/pressure source moments, while (F_+(q)
eq F_-(q)) pointwise.

The question is therefore not whether the background moments differ. They do not, by construction. The question is whether the retarded collisionless response is uniquely determined by those moments.

## 2. Generalization of the single-halo Vlasov solution

For a non-relativistic collisionless species in an external CDM potential, the linearized Boltzmann/Vlasov equation in the time variable (d t/deta=a^2) is

[
rac{partial f_1}{partialeta}
+rac{mathbf q}{m_
u}cdot
abla_x f_1
=
m_
u a^2(eta),
abla_xphicdot
abla_q F(q),
]

where (F(q)) is any isotropic background distribution.

The standard moving-point-mass solution therefore does not require a Fermi–Dirac shape. Replacing the thermal (f_{FD}) by an arbitrary isotropic (F) gives, up to the common phase-space normalization convention,

[
deltaho_{
u,F}(eta,mathbf x)
=
C(eta)
int_0^eta deta',
rac{a(eta')}{(eta-eta')^2}
F!left(q_star(eta',mathbf x)ight),
]

with

[
q_star(eta',mathbf x)
=
rac{m_
u}{T_{
u,0}}
rac{|mathbf x-mathbf x_H(eta')|}
{eta-eta'}.
]

The constant (C) contains the common factors (GM m_
u^4), powers of (a), and the chosen normalization of (F).

Hence for the frozen pair,

[
oxed{
deltaho_{
u,-}-deltaho_{
u,+}
=
C(eta)
int_0^eta deta',
rac{a(eta')}{(eta-eta')^2}
left[F_-(q_star)-F_+(q_star)ight].
}
]

This is the central K01E identity.

## 3. Why finite moment closure is insufficient

The source matching imposes a finite set of linear constraints

[
int dq,w_i(q),[F_-(q)-F_+(q)]=0,
qquad i=1,ldots,N_m.
]

The wake is a different linear functional:

[
mathcal W_{eta,mathbf x}[Delta F]
=
int_0^eta deta',
K(eta,eta',mathbf x),
Delta F(q_star(eta',mathbf x)).
]

Finite moment matching would force (mathcal W[Delta F]=0) only if every wake functional of interest were contained in the finite span of the matched moment functionals. The moving-halo kernel instead samples/evaluates (Delta F) along a continuum of kinematically selected momenta (q_star).

Therefore, generically,

[
oxed{
{mathcal M_i[F_+]=mathcal M_i[F_-]}_{i=1}^{N_m}
;
otRightarrow;
deltaho_{
u,F_+}
=
deltaho_{
u,F_-}.
}
]

This is a statement about finite moment closure, not a claim that the complete stress-energy tensor is identical at all times or all scales.

## 4. Static constant-velocity form

For a halo moving with constant velocity in a static background, the late-time halo-frame wake becomes

[
deltaho_{
u,F}^{H}(mathbf x)
propto
int_0^infty du,
F!left(|umathbf x+oldsymbolalpha|ight),
]

where

[
oldsymbolalpha=
rac{m_
u a}{T_{
u,0}}mathbf v_H.
]

Thus the hidden-state contrast is directly

[
oxed{
Deltadeltaho_
u^{H}(mathbf x)
propto
int_0^infty du,
Delta F!left(|umathbf x+oldsymbolalpha|ight).
}
]

This is a line/ray transform of the hidden phase-space difference. Matching a finite number of radial moments does not determine this transform.

## 5. Resonant anchor for the frozen E8 pair

Using the frozen parameters

- (m_
u=0.06) eV,
- (z=0.95),
- (|v_parallel|=200) km/s,
- (T_{
u,0}=0.71611,T_{CMB}=1.9518) K,

the non-relativistic resonant momentum is

[
q_{m res}
=
rac{m_
u |v_parallel|}
{c,T_{
u,0}(1+z)}
=
0.1220470.
]

At that momentum the frozen E8 construction gives

[
rac{F_+}{F_{FD}}=0.8147976,
qquad
rac{F_-}{F_{FD}}=1.1852024,
]

and therefore

[
oxed{
rac{F_+-F_-}{F_{FD}}
=
-0.3704048.
}
]

So the matched pair differs by about 37% of the FD occupancy exactly in the low-(q) region selected by a 200 km/s wake resonance.

This is a source-level kinetic contrast, not an absolute galaxy signal.

## 6. Independent nonlinear/local anchor

The deterministic K01C truncated-NFW scattering calculation gave

[
rac{D_- - D_+}{|D_{FD}|}
=
0.43399
]

for the preregistered primary (b_{max}=r_{200}) local scattering proxy.

The 37.0% E8 resonant occupation contrast and 43.4% K01C local scattering contrast are **not the same observable** and must not be equated numerically. Their common ordering,

[
F_- > F_+ quad	ext{in the relevant low-}q	ext{ region},
]

together with

[
D_- > D_+,
]

is nevertheless physically consistent with the wake being sensitive to hidden low-momentum occupation.

## 7. Primary K01E observable

Define the halo-frame Legendre dipole of the neutrino wake,

[
delta_{
u,1}^{F}(r)
=
rac{3}{2}
int_{-1}^{1}dmu,
mu,
delta_{
u,F}^{H}(r,mu),
]

where (mu=hat{mathbf r}cdothat{mathbf v}_H).

The hidden-state observable is

[
oxed{
Deltadelta_{
u,1}(r)
=
delta_{
u,1}^{F_-}(r)
-
delta_{
u,1}^{F_+}(r).
}
]

A useful dimensionless response ratio is

[
oxed{
R_{m hid}(r)
=
rac{Deltadelta_{
u,1}(r)}
{|delta_{
u,1}^{FD}(r)|}.
}
]

This is the primary physics target for K01F.

## 8. Falsifiable prediction

For halo velocities whose kinematic support substantially samples the E8 low-(q) region around (qsim0.1)–0.3, the frozen construction predicts a non-zero hidden-state wake dipole.

The sign and radial structure of (Deltadelta_{
u,1}(r)) must be computed from the full line integral above; they must not be inferred globally from a single resonant point.

A null result after evaluating the full wake functional would falsify the proposed observable consequence for the chosen matched pair even though the distributions remain microscopically different.

## 9. Scientific interpretation

The target statement is:

> A finite set of matched stress-energy moments does not uniquely determine the retarded gravitational response of a collisionless medium; a moving gravitational perturbation can resolve hidden phase-space structure through its velocity-selective wake kernel.

This should not be strengthened to “identical full stress-energy tensors have different Einstein evolution.” The present construction matches selected source moments at the frozen source epoch.

## 10. Next physics stage

**K01F:** evaluate the generalized moving-halo wake integral for (F_{FD},F_+,F_-) and compute the radial dipole profiles (delta_{
u,1}^{F}(r)) and (R_{m hid}(r)).

The first K01F calculation should remain external-potential / BKT linear-response physics. SWIFT or another self-gravitating engine is reserved for later validation, not for defining the effect.
