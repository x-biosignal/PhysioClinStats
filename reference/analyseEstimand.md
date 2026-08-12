# Analyse an estimand under multiple imputation for dropout

End-to-end analysis of a declared
[`defineEstimand`](https://x-biosignal.github.io/PhysioClinStats/reference/defineEstimand.md)
estimand: multiply-impute the dropout-missing longitudinal response (in
wide form, so the within-subject correlation is respected), fit an MMRM
([`fitMMRM`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md))
to each completed data set, and pool the fixed effects with Rubin's
rules
([`poolEstimates()`](https://x-biosignal.github.io/PhysioClinStats/reference/poolEstimates.md)).
Treatment-policy and hypothetical strategies use the standard MAR
imputation; other strategies warn that they need a bespoke imputation
model.

## Usage

``` r
analyseEstimand(
  data,
  estimand,
  response,
  treatment,
  time,
  subject,
  covariates = NULL,
  m = 20,
  method = "norm",
  seed = NULL,
  covariance = "unstructured"
)
```

## Arguments

- data:

  Long-format data (subject x time rows) with the response, treatment,
  time, subject, and any (subject-level) covariates.

- estimand:

  An `"estimand"` from
  [`defineEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/defineEstimand.md).

- response, treatment, time, subject:

  Column names.

- covariates:

  Optional subject-level covariate column names.

- m:

  Number of imputations (default 20).

- method:

  mice method (default `"norm"`, appropriate for a continuous
  longitudinal response).

- seed:

  Optional integer seed.

- covariance:

  MMRM covariance structure (default `"unstructured"`).

## Value

An `AnalysisResult` (type `"estimand_analysis"`) carrying the estimand
attributes and the pooled fixed-effect table.

## References

ICH E9(R1); Rubin 1987; Mallinckrodt 2008 (MMRM).

## See also

[`defineEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/defineEstimand.md),
[`multipleImputation()`](https://x-biosignal.github.io/PhysioClinStats/reference/multipleImputation.md),
[`poolEstimates()`](https://x-biosignal.github.io/PhysioClinStats/reference/poolEstimates.md)
