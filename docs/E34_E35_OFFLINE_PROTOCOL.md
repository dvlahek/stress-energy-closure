# E34–E35 — offline causal-history and low-k profile-control protocol

Date: 2026-09-29. Scope: mathematical follow-up to already completed E29–E33 and the frozen E27/E28 **conditional** externally forced Born model. Both studies below are illustrative analytic constructions and independent standard-library arithmetic QA. This is neither new halo/survey data nor physical calibration of E28.

## Immutable inputs and operational seal

Read only the already supplied small E29/E30–E33 artefacts for provenance. No ASDF, FITS, CLASS, halo arrays, new mocks, GitHub mutation, WSL execution or observed eBOSS rows/odd multipoles. Original E8 F0/F+/F−, E16, E27/E28 sources, E4 mocks, cuts, random seeds, and LOS orientation remain untouched. Draft PR #1 and `main` unchanged.

## E34 — two endpoint snapshots do not identify a retarded profile history

For a positive spherical thin shell of constant mass M, set `r_±(t)=r0 ± eps sin²(pi t/T)`, `t ∈ [0,T]`, `r0>eps>0`, same centre, same mass, same endpoint shell radii and zero radial velocity at both endpoints. The shell transport is a weak solution of radial mass continuity, because each mass element obeys `dr/dt=v_r(t)` and total mass is conserved. Compute `u_±(k,t)=sinc(k r_±(t))` and a common endpoint test factor `u_test=sinc(k r0)`. Integrate illustrative, fixed positive causal weight `K(t)=exp(-(T-t)/T)`; `J_± = u_test ∫ K(t) u_±(k,t)dt` (T=r0=k=1 and eps=.1, dimensionless). Show `J_+ != J_-` although both snapshots and endpoint radial velocities match; this toy kernel is **not** the E28 physical `K_F`. Provide a rigorous sufficient sign argument using monotonic sinc on [0.9,1.1]. For an arbitrary nontrivial signed continuous physical K, the endpoint-to-history nonuniqueness is generic under localized interior perturbations; do not assert a specific E28 numerical gap without its kernel.

Derive a conditional stability statement: if a normalized profile obeys a separately established temporal Lipschitz constant `L= sup |du/dt|` between endpoints, then the linear endpoint interpolant `u_lin` satisfies `|u(t)-u_lin(t)| <= 2 L t(T-t)/T` and `|∫K(u-u_lin)dt| <= 2L/T ∫|K|t(T-t)dt`. For equal endpoint u this is conservative. Do **not** infer L from two samples or use a toy L as a bound for Abacus objects. Check 64,128,256,512,1024-panel Simpson consistency, continuity identity for test functions, endpoint equality and deliberately wrong sign.

## E35 — a profile-independent low-k bound with explicit inputs

For a positive, normalized spherical radial mass distribution with finite `<r²>`, `u(k)=E[sinc(kr)]`, and `sinc(x)=E_{μ∈[-1,1]} cos(x μ)`, derive `0<=1-u(k)<=k²<r²>/6`, `|u(k)|<=1`. Thus for independent positive source/test profiles `|1-u_src u_test| <= k²(<r²>src+<r²>test)/6`. Within a frozen linear externally forced response `A(k)=u_test ∫K(k,t)m(t)u_src(k,t)dt`, positive `m(t)` and potentially signed K, the absolute point-mass replacement error is bounded by `k²/6 ∫|K|m(t)(<r²>test+<r²>src(t))dt`. The theorem assumes independently justified radii/moments for **all source times**, not merely endpoints, and does not bound higher-order Vlasov corrections, environmental mass, anisotropy, or E28 UV force.

Test finite positive toy radial mixtures and both E34 shell histories across a fixed dimensionless k grid. Include an adversarial **same-endpoints** shell with a large transient radius to prove endpoint-only moments cannot be inserted into a uniform-in-time bound. All numerical wave numbers/radii here are dimensionless and must never be interpreted as an E28 k cutoff.

## Reporting

Save E34 and E35 original scripts and separate JSON reports with parent SHA256 fingerprints. Fail closed on a single failed assertion, preserve successful originals on replay, distinguish exact symbolic proof from toy quadrature. Package reproducible script/JSON/notes/handoff into a ZIP. Do not infer a physical halo drag, a calibrated F+ vs F− template, eBOSS significance, or permission to unblind. A-03 PHYSICAL_UNCERTIFIED; A-04 BLOCKED; observed odd SEALED.
