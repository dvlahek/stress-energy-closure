# E60 — ELG reconstruction mask robustness

E59 passed its preregistered internal repeatability screen for ELG and failed it for LRG. The strongest remaining E59-specific caveat is that both split halves shared a single frozen R1 48k random selection field. A common random-field error can artificially increase half-half agreement.

E60 changes no physical model and does not revisit tracer selection. It fixes ELG from E59, replays the same 512 probe midpoints and R=16 Mpc/h velocity kernel, and recomputes the full ELG reconstruction with the exact E58 R1--R7 48k random subsets. All 21 replica pairs are compared. Random-replica RMS is also compared with the seven-replica mean field.

For a direct common-bias check, mock0001 in each cap is additionally evaluated with the complete source-eligible ELG random pool. This is intentionally only an anchor because the full pool is much more expensive.

Passing E60 closes the numerical random-mask caveat only. It does not validate the reconstructed field against true velocity. That remains the next gate.
