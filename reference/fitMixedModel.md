# Fit a linear mixed-effects model

Fits a linear mixed model with
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) and returns tidy
fixed-effect and random-effect tables in an `AnalysisResult`.
Denominator degrees of freedom (and the associated t tests) use
`lmerTest` when it is installed - either the Satterthwaite or the
Kenward-Roger approximation.

## Usage

``` r
fitMixedModel(
  data,
  formula,
  random = NULL,
  method = c("REML", "ML"),
  df = c("satterthwaite", "kenward-roger")
)
```

## Arguments

- data:

  A data.frame in long format.

- formula:

  A model formula. It may carry the random-effects terms directly (e.g.
  `y ~ x + (x | subject)`), or supply them via `random`.

- random:

  Optional one-sided random-effects formula (e.g. `~ (1 | subject)`)
  appended to `formula` when the latter has none.

- method:

  `"REML"` (default) or `"ML"`.

- df:

  Denominator-df method for the fixed-effect tests: `"satterthwaite"`
  (default) or `"kenward-roger"` (both need `lmerTest`).

## Value

An `AnalysisResult` (`type = "mixed_model"`) whose `estimate` is the
fixed-effect coefficient vector, with `result$fixed` (the tidy
fixed-effect table), `result$random` (variance components),
`result$sigma`, `result$fit` (the fitted model, for
[`estimatedMarginalMeans()`](https://x-biosignal.github.io/PhysioClinStats/reference/estimatedMarginalMeans.md))
and `result$formula`.

## References

Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear
mixed-effects models using lme4. Journal of Statistical Software, 67(1).
Kuznetsova A, Brockhoff PB, Christensen RHB (2017). lmerTest package.
Journal of Statistical Software, 82(13).

## See also

[`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md),
[`estimatedMarginalMeans()`](https://x-biosignal.github.io/PhysioClinStats/reference/estimatedMarginalMeans.md)

## Examples

``` r
if (requireNamespace("lme4", quietly = TRUE)) {
  data(sleepstudy, package = "lme4")
  fitMixedModel(sleepstudy, Reaction ~ Days + (Days | Subject))
}
#> <AnalysisResult> mixed_model 
#>   estimate: 251.40510,  10.46729 
#>   method: REML 
#>   fields: fixed, random, sigma, fit, formula, df_method 
#>   provenance: 1 entr(ies)
```
