# Emulate a declared target trial

Applies eligibility at the declared time zero and estimates either an
intention-to-treat contrast or a per-protocol clone-censor-weight
contrast. The result remains conditional on consistency,
exchangeability, positivity, correct treatment/censoring/outcome models,
and a correctly aligned time zero. A Cox hazard ratio is non-collapsible
and is not a marginal risk ratio.

## Usage

``` r
targetTrialEmulate(
  data,
  protocol,
  id,
  time,
  treatment,
  outcome,
  event = NULL,
  baseline_covariates = NULL,
  time_varying_covariates = NULL,
  estimand = c("per_protocol", "intention_to_treat"),
  weight_backend = c("internal", "WeightIt", "ipw"),
  stabilized = TRUE,
  numerator_covariates = NULL,
  weight_truncation = c(0.01, 0.99),
  max_weight = Inf,
  conf.level = 0.95
)
```

## Arguments

- data:

  Long person-period data.

- protocol:

  A
  [`targetTrialProtocol()`](https://x-biosignal.github.io/PhysioClinStats/reference/targetTrialProtocol.md)
  object.

- id, time, treatment, outcome:

  Column names. For survival outcomes, `outcome` is the subject-level
  outcome/censoring time repeated on each long row and `event` is its
  repeated status. Without `event`, `outcome` is the repeated
  fixed-horizon binary outcome.

- event:

  Optional binary event-status column.

- baseline_covariates, time_varying_covariates:

  Declared predictors for treatment/adherence models.

- estimand:

  `"per_protocol"` or `"intention_to_treat"`.

- weight_backend:

  Propensity backend. Optional backends are validation paths and never
  silent fallbacks.

- stabilized:

  Use time-based numerator probabilities.

- numerator_covariates:

  Additional numerator-model covariates.

- weight_truncation:

  Lower and upper quantiles applied after raw cumulative weights are
  retained.

- max_weight:

  Additional positive upper cap applied after quantile truncation.

- conf.level:

  Confidence level.

## Value

An `AnalysisResult` containing effects, models, clone rows, both raw and
analysis weights, diagnostics, censoring counts, and provenance.
