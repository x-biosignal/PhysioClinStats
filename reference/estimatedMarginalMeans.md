# Estimated marginal means

Computes estimated marginal (least-squares) means from a fitted mixed or
MMRM model with `emmeans`, optionally with a follow-up contrast,
returning the marginal-mean table (and contrasts) in an
`AnalysisResult`.

## Usage

``` r
estimatedMarginalMeans(model, specs, contrasts = NULL, level = 0.95, ...)
```

## Arguments

- model:

  An `AnalysisResult` from
  [`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md)
  /
  [`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md),
  or a raw model object `emmeans` understands.

- specs:

  A formula or character spec for the marginal means, e.g.
  `~ treatment | time`.

- contrasts:

  Optional contrast method passed to
  [`emmeans::contrast`](https://rvlenth.github.io/emmeans/reference/contrast.html)
  (e.g. `"pairwise"`, `"trt.vs.ctrl"`); `NULL` (default) returns means
  only.

- level:

  Confidence level for the intervals (default 0.95).

- ...:

  Further arguments forwarded to
  [`emmeans::emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.html).

## Value

An `AnalysisResult` (`type = "emmeans"`) with `result$emmeans` (the
marginal-mean table) and, when requested, `result$contrasts`.

## References

Lenth RV (2024). emmeans: Estimated Marginal Means. R package.

## See also

[`pairwiseContrasts()`](https://x-biosignal.github.io/PhysioClinStats/reference/pairwiseContrasts.md),
[`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md)

## Examples

``` r
if (requireNamespace("mmrm", quietly = TRUE) &&
    requireNamespace("emmeans", quietly = TRUE)) {
  data(fev_data, package = "mmrm")
  fit <- fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
                 covariates = c("RACE", "SEX"))
  estimatedMarginalMeans(fit, ~ ARMCD | AVISIT)
}
#> <AnalysisResult> emmeans 
#>   estimate: c("PBO", "TRT", "PBO", "TRT", "PBO", "TRT", "PBO", "TRT"), c("VIS1", "VIS1", "VIS2", "VIS2", "VIS3", "VIS3", "VIS4", "VIS4"), c("33.33186", "37.10609", "38.17145", "41.90375", "43.67397", "46.75452", "48.38576", "52.78422"), c("0.7577720", "0.7641827", "0.6111675", "0.6006462", "0.4583865", "0.5058662", "1.1740752", "1.1710810"), c("148.1457", "143.1765", "147.0330", "143.4844", "129.8027", "130.1341", "134.0814", "132.6244"), c("31.83442", "35.59555", "36.96364", "40.71650", "42.76710", "45.75373", "46.06366", "50.46780"), c("34.82930", "38.61663", "39.37926", "43.09101", "44.58085", "47.75531", "50.70786", "55.10063") 
#>   method: estimated marginal means 
#>   fields: emmeans, contrasts, emm_object 
#>   provenance: 1 entr(ies)
```
