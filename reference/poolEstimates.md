# Pool estimates across imputations (Rubin's rules)

Combines the coefficient estimates from a list of models fitted to
multiply imputed data sets using Rubin's rules with the Barnard-Rubin
degrees-of- freedom adjustment. The pooled estimate, total variance, df,
confidence interval, and fraction of missing information reproduce
[`mice::pool`](https://amices.org/mice/reference/pool.html).

## Usage

``` r
poolEstimates(fits, conf_level = 0.95)
```

## Arguments

- fits:

  A list of fitted models (each with
  [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html)), or a list of tidy data
  frames with columns `term`, `estimate`, and `std.error` (and
  optionally `df.residual`).

- conf_level:

  Confidence level for the pooled interval (default 0.95).

## Value

An `AnalysisResult` (type `"pooled_estimates"`) whose `result$estimates`
is a data frame of pooled `estimate`, `std.error`, `df`, `statistic`,
`p.value`, `conf.low`, `conf.high`, and `fmi` per term.

## References

Rubin 1987; Barnard & Rubin 1999.
[`mice::pool`](https://amices.org/mice/reference/pool.html).

## See also

[`multipleImputation()`](https://x-biosignal.github.io/PhysioClinStats/reference/multipleImputation.md),
[`analyseEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/analyseEstimand.md)

## Examples

``` r
set.seed(1)
fits <- lapply(1:5, function(i) lm(mpg ~ hp + wt,
  data = mtcars[sample(nrow(mtcars), replace = TRUE), ]))
poolEstimates(fits)
#> <AnalysisResult> pooled_estimates 
#>   estimate: 37.90362826, -0.02973744, -4.27704284 
#>   method: rubin_pool 
#>   fields: estimates 
```
