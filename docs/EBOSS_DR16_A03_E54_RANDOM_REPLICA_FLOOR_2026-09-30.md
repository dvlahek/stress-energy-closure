# E54 — finite-random floor at fixed full mock galaxies

E53 Stage 2 established that the conditioned F+/F- fingerprint is stable but the marked-LS background still decreases strongly from 4800R to 48000R. E54 therefore keeps mock0001 full eligible galaxies, the four E51 technical sign fields, the E19 F+/F- injection basis, the pair geometry and projection fixed, and varies only the Monte-Carlo random catalogue.

Four fresh 48k random catalogues per cap/tracer are drawn with predeclared deterministic PCG64 seeds from the same frozen source-eligible random pools. Replicas can overlap because the source pool is finite; E54 reports those overlaps and does not treat sd/sqrt(4) as a rigorous independent-sample error.

The primary amplitude-calibrated coordinate is m_d=(d.b)/(d.d), d=q_minus-q_plus. Adding lambda*d shifts m_d by exactly lambda. The preregistered engineering gate is median replica SD <=0.01 across tag fields in each cap, max field SD <=0.02, post-window angle range <=0.002 rad, and F-state difference-norm-ratio range <=0.01. This is only a technical random-integration gate relative to unit lambda, not a physical eBOSS sensitivity threshold.

Local script SHA256: 658e91ad7d3434385844f149762d88f336566fd0928dfe248236b7d9b33b7f87
Local package SHA256: fa091a20df15989fc6ac7efb10f2029f58faaeeb1843802a539358a88429df73
