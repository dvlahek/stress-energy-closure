# E57 — six-replica random-integration closure

E56 materially improved the numerical integration but missed the already frozen median precision target by a small amount in both caps. NGC had median/max effective-random-to-between-mock ratios 0.289/0.443 and SGC 0.289/0.424. The per-field maximum criterion 0.50 already passed in both caps; only the median <=0.25 criterion remained open.

E57 does not change that target. Using standard Monte-Carlo 1/sqrt(R) scaling, the minimum integer number of 48k replicas predicted to bring both medians below 0.25 is R=6. Five replicas project to about 0.259 in each cap, while six project to about 0.236. This is a numerical resource decision made while observed galaxy rows and the observed odd vector remain sealed.

R1-R4 are reused exactly from E56. E57 adds only R5/R6. The final statistic is the six-replica mean within each mock. The frozen E56 final-mean gate remains: median effective-random/between-mock SD <=0.25, every field <=0.50, and both caps must pass.

Passing E57 would close only the conditional random-subset integration gate. It would not certify a physical velocity tag, absolute conditioned amplitude, finite-parent-pool common bias, covariance inversion, significance, detection or exclusion.
