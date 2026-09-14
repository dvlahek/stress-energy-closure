# Phase-7 local EZmock placebo ensemble

This directory archives compact outputs from the local DESI DR1 BGS EZmock placebo-covariance pipeline.

The public DR1 EZmock BGS products used here do not expose luminosity columns. The tracer split is therefore an equal-count rank constructed from the released `RAN_NUM_0_1` field within the same narrow-redshift and NGC/SGC cells used by the real-data luminosity proxy. These files are a survey-geometry/covariance/systematics control, not a luminosity-matched mock likelihood and not a halo-mass wake constraint. The Abacus layer remains the physical luminosity-ranked mock validation.

## Archived realizations

- `mock_02_*` — validated local smoke realization with 2x random density and realization-specific released EZmock random catalogs. It is retained only as a transport/estimator validation and is **excluded from the final homogeneous covariance ensemble**.
- `mock_03_*` — first accepted realization of the homogeneous production ensemble. It uses one validated full NGC/SGC survey-random pair as the fixed selection-function pool, while `make_random_proxy` preserves each current mock's counts in narrow-z/cap cells. This is the production random-geometry mode used for mock 3 onward.

The earlier experimental shared-angular-cache/redshift-remapping route was rejected because it produced order-unity estimator artifacts and distorted the forward window. That implementation has been removed from the repository.

## Local-file SHA256 provenance

- `mock_02_summary.json`: `d453bcad27d4385cefcd5f7db3bde01ce82ff106373fa694c328aee9d11bf68c`
- `mock_02_vector.csv`: `68ddd23f03a82c77c8e9ffd2316c30dfac68819f1ea718f7986d8d719d2bf75f`
- `mock_02_window.csv`: `9ff6e47cb8351bc45fcc89b4f6a35137ef5341ace8cd9afc38b747fa646f0fc3`
- `mock_03_summary.json`: `f049d40aeeacd2a942a4dcf0d680fbf173615d0052dd070835463b8852c61c0d`
- `mock_03_vector.csv`: `25418a1722cf53c4e7ae056f50a8e7687d43655c221a03e3ebfe58d50e47e95e`
- `mock_03_window.csv`: `aa696db174eec5a6be47a40940e0950cfc63c58c3a7420bc941e5c83e0a0b015`

Mock 3 passed the explicit sanity gate with `max|xi0| = 0.1048502089` and `max|xi1| = 0.0160603752`.
