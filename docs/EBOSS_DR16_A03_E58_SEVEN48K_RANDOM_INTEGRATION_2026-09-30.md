# E58 — seven-replica 48k random-integration closure

E57 passes SGC but NGC misses the frozen median gate by only 0.001413: 0.251413 > 0.25. All per-field maxima already pass comfortably.

E58 adds exactly one R7 replica to every mock/cap and keeps all physics, geometry, tag fields, injection basis, cuts and the E56/E57 gate unchanged. R=7 is the minimum integer implied by the E57 NGC result under standard 1/sqrt(R) Monte-Carlo scaling,

```
ceil[6*(0.2514132726578298/0.25)^2] = 7.
```

The final-mean gate remains:
- median effective-random / between-mock SD <= 0.25,
- every field <= 0.50,
- both NGC and SGC must pass.

This is numerical integration control only. Observed galaxy rows and the observed odd vector remain sealed; no physical velocity reconstruction, covariance inversion, p-values or detection/exclusion inference are authorized.

**Anti-chasing rule:** if E58 still fails, stop replica-by-replica extension and move to deterministic or higher-density random integration.
