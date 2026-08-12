# PhysioClinStats 0.2.1

- Relaxed the two LAPACK-derived medsens sensitivity-curve comparisons
  (acme_control/d0, ade_control/z0) in the pinned mediation fixture test from
  1e-12 to 1e-8. These curves are RNG-invariant analytic functions computed
  through matrix/QR operations and drift by up to ~2.1e-10 relative across
  BLAS/LAPACK implementations (measured OpenBLAS vs reference LAPACK); the
  fixture was pinned under reference LAPACK. RNG-reproduced quantities keep
  1e-12. Fixes the r-universe binary-check failure.

# PhysioClinStats 0.2.0

## Cross-tool validation

- Added strict offline parity fixtures for uncertainty, SCED, mixed/MMRM,
  survival, recovery, imputation, and causal inference. The gate rejects
  missing or mutated fixtures, stale surface inventory, direction drift,
  unexpected warnings, and unexecuted required rows.
- Added deterministic publication CSV, JSON, and Markdown tables plus an
  independent 100-row and 100-mutation audit.

## Causal-inference workflows

- Added `causalMediation()`, a guarded `mediation::mediate()` workflow that
  verifies retained model rows, preserves the caller's RNG state, reports
  control, treated, and average ACME/ADE estimates, and records supported
  `medsens()` diagnostics.
- Added `targetTrialProtocol()` for explicit seven-component target-trial
  declarations with a separate time-zero rule and provenance representations
  for static and dynamic strategies.
- Added `targetTrialEmulate()` for intention-to-treat and per-protocol
  clone-censor-weight analyses of long person-period data. Raw and truncated
  weights, censor reasons, positivity and balance diagnostics, and
  clone-clustered robust uncertainty are retained in the result.
- Added guarded WeightIt and ipw parity paths for binary baseline IPTW.
- Added pinned offline backend fixtures and a 100-seed independent numerical
  validation with deliberate future-information, cumulative-product-grouping,
  and treatment-coding mutations.

These workflows estimate contrasts only under their stated consistency,
exchangeability, positivity, time-zero, and model assumptions. They do not
turn observational data into a randomized trial or prove a biological
mechanism.

# PhysioClinStats 0.1.0

Initial release of PhysioClinStats as a standalone package in the x-biosignal
ecosystem. This package provides the clinical inference layer for
rehabilitation and physiological signal studies, with heavy modelling backends
kept as optional, guarded dependencies.

## New Features

- Optional-backend detection and guarding for the clinical inference engine.
  Modelling packages are loaded only when present, so the package installs
  and loads without pulling in heavy dependencies.
  - `clinStatsBackends()` reports availability of the optional modelling
    backends (`mmrm`, `emmeans`, `lme4`, and `SingleCaseES`) as a named
    logical vector, letting callers check support before running an analysis.
  - `requireBackend()` asserts that a named optional backend is installed and
    otherwise stops early with an actionable install message.

## Documentation

- Package-level documentation describing the intended scope: mixed-effects /
  MMRM longitudinal models, single-case (N-of-1) designs, estimands with
  multiple imputation, and per-estimate uncertainty.
