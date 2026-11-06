# K01G — Retarded memory unlocking of hidden phase-space structure

**Date:** 2026-11-06  
**Branch:** `experiment/k01-matched-f-detectability-20261105`  
**Scope:** physics reinterpretation of the existing E17D0/E27 retarded Vlasov source using the frozen E8 matched states. No new observations, no survey covariance, no self-gravity.

## 1. Why E27 looked tiny while K01F is large

K01F gives a large late-time/static hidden wake dipole because the stationary moving-halo response probes the cumulative low-momentum population.

The existing E27 calculation instead imposed:

- zero incoming wake at (z=1),
- observation at (z=0.95),
- one large-scale mode (k=0.05,h/{m Mpc}),
- a short retarded memory interval,
- an external point-halo Born source.

E27 therefore probes an **early-time / small-phase transient**, not the stationary wake.

Its precise hidden contrast was

[
rac{operatorname{Im}R_- - operatorname{Im}R_+}
{operatorname{Im}R_{FD}}
=
1.02417	imes10^{-5}
]

for the alpha=0.4 halo-history anchor, and (1.02329	imes10^{-5}) for alpha=0.8.

The smallness is physically explained by finite-moment matching.

## 2. Retarded source kernel

The frozen E17D0 source is

[
H_F(s)
=
-rac{1}{I_2[F]}
int dq,q^2F'(q),j_0'(sq),
]

with

[
I_2[F]=int dq,q^2F(q),
]

and

[
j_0'(x)
=
-rac{x}{3}
+rac{x^3}{30}
-rac{x^5}{840}
+O(x^7).
]

Integrating term by term by parts, with negligible finite-support boundary terms,

[
H_F(s)
=
-s
+rac{s^3}{6}rac{I_4[F]}{I_2[F]}
-rac{s^5}{120}rac{I_6[F]}{I_2[F]}
+O(s^7),
]

where

[
I_n[F]=int dq,q^nF(q).
]

## 3. Moment-protected transient

The frozen E8 states satisfy the same number-density moment:

[
I_2[F_+]=I_2[F_-]=I_2[F_{FD}]
]

to machine precision.

Therefore the leading term

[
H_F(s)=-s+O(s^3)
]

is exactly state-independent.

The hidden-state difference begins only at cubic order:

[
H_-(s)-H_+(s)
=
rac{s^3}{6}
left[
rac{I_4[F_-]}{I_2[F_-]}
-
rac{I_4[F_+]}{I_2[F_+]}
ight]
+O(s^5).
]

Dividing by the leading FD response (H_{FD}simeq-s),

[
oxed{
rac{H_-(s)-H_+(s)}{H_{FD}(s)}
=
-rac{s^2}{6}
Delta!left(rac{I_4}{I_2}ight)
+O(s^4).
}
]

For the frozen E8 pair,

[
rac{I_4}{I_2}igg|_{FD}=12.93920,
]

[
rac{I_4}{I_2}igg|_{+}=13.03174,
qquad
rac{I_4}{I_2}igg|_{-}=12.84666,
]

so

[
oxed{
rac{H_- - H_+}{H_{FD}}
simeq
0.030847,s^2
qquad(sll1).
}
]

This is the K01G memory-unlocking law.

## 4. E27 phase budget

E27 used

[
s(z)
=
k_{m com}
rac{T_{
u,0}}{m_
u}
int_{z_{m obs}}^{z}
rac{1+z'}{H(z')}dz'.
]

For the frozen FD CLASS background,

[
z_{m obs}=0.95,quad z_{m init}=1.0,
quad k=0.05,h/{m Mpc},
]

the endpoint phase is only

[
oxed{s_{max}=0.02351.}
]

The retarded integral is weighted over (0<s<s_{max}), so its effective phase is smaller than the endpoint.

Using the exact frozen source kernel, the E27 contrast of approximately (1.02	imes10^{-5}) corresponds to an effective source phase near

[
s_{m eff}simeq0.0183.
]

The small-phase law predicts

[
0.030847,s_{m eff}^2
simeq1.0	imes10^{-5},
]

which quantitatively explains the archived E27 result.

## 5. Phase-space unlocking

The exact frozen source kernel shows that the hidden contrast grows rapidly once the phase is no longer small.

Representative exact source-level thresholds are:

| hidden source contrast | first (s) |
|---:|---:|
| (10^{-4}) | (ssimeq0.057) |
| (10^{-3}) | (ssimeq0.189) |
| (1%) | (ssimeq0.636) |
| (5%) | (ssimeq0.833) |
| (10%) | (ssimeq0.989) |

The precise response becomes oscillatory at larger (s), so this table is a first-threshold diagnostic rather than a monotonic asymptotic law.

## 6. Scale dependence at fixed short history

For the same E27 history (z=1	o0.95), (spropto k).

Since

[
s=0.02351
]

at (k=0.05,h/{m Mpc}), the corresponding first-threshold scales are approximately

[
oxed{
k_{0.1%}simeq0.40,h/{m Mpc},
}
]

[
oxed{
k_{1%}simeq1.35,h/{m Mpc},
}
]

[
oxed{
k_{10%}simeq2.10,h/{m Mpc}.
}
]

Thus the matched hidden-state signal is parametrically suppressed on large scales and is progressively unlocked on smaller scales or over longer retarded histories.

## 7. Connection to E28

The earlier E28 finite-band halo-recoil calculation was strongly dominated by short scales: about 94% of the FD finite-band acceleration in the low-mass alpha=0.4 benchmark came from (1<kle8,{m Mpc}^{-1}).

K01G now gives a physical explanation for why hidden-state sensitivity naturally migrates toward those scales:

- large-scale / short-memory response is moment-protected,
- higher (k) accumulates phase faster,
- finite moment matching ceases to shield the detailed distribution once (s=O(1)).

The E28 UV sensitivity still prevents using its total amplitude as a certified physical halo drag. K01G does **not** rehabilitate the E28 amplitude. It explains why the hidden-state dependence was expected to become stronger in the same small-scale regime.

## 8. Relation to K01F

K01F and K01G describe opposite temporal regimes of the same collisionless response.

### Early-time / small-phase regime

[
sll1:
qquad
rac{Delta H}{H}
propto s^2.
]

Finite moment matching protects the response.

### Late-time / stationary regime

The moving-halo dipole approaches the cumulative low-momentum functional

[
W_{1,F}(alpha)
=
-rac{3}{alpha^2}
int_0^alpha q^2F(q)dq,
]

which gives a 32.5% hidden dipole contrast for the frozen pair at 200 km/s.

Therefore the large static K01F result and tiny short-memory E27 result are not contradictory.

They are the two ends of a **memory-unlocking transition**.

## 9. Main physical statement

A more complete statement of the project is now:

> Finite moment matching produces a moment-protected transient in collisionless gravitational response. At short retarded times or large spatial scales the response is nearly indistinguishable, because the leading Vlasov kernel is fixed by the matched moments. Phase mixing at longer memory or higher wavenumber progressively unlocks sensitivity to the hidden phase-space distribution.

Symbolically,

[
oxed{
	ext{moment closure}
;Rightarrow;
	ext{temporary large-scale degeneracy},
}
]

but

[
oxed{
	ext{phase mixing}
;Rightarrow;
	ext{hidden-state observability}.
}
]

## 10. Scientific significance

This is stronger than simply stating that two distributions give different responses.

The new mechanism predicts **where and when** the difference should appear:

[
oxed{
ssim
k,v_{
u,m th},Delta t_{m memory}
}
]

acts as the control parameter.

The hidden state is effectively invisible for (sll1) and becomes dynamically resolvable as (s) approaches order unity.

This provides a direct bridge between:

- finite-moment non-uniqueness,
- collisionless memory,
- phase mixing,
- scale dependence,
- and velocity-oriented neutrino wakes.

## 11. Next physics stage

K01H should no longer ask whether a hidden response exists.

It should ask:

[
oxed{
	ext{What realistic halo histories and spatial scales maximize }
s
	ext{ while remaining physically modelable?}
}
]

The immediate target is a physically regulated halo wake across (ksim0.1)–a few (h/{m Mpc}), with finite halo profile and long causal history.

This is the scale/time window where K01G predicts the transition from moment-protected behavior to observable hidden-state response.
