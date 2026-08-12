# Conformal prediction band for a single new patient

Convenience wrapper around split-conformal
[`conformalInterval()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalInterval.md)
for one new observation, returning the guaranteed-coverage predicted
band.

## Usage

``` r
conformalPredict(model, calib, newx, alpha = 0.1, seed = NULL)
```

## Arguments

- model:

  A fitted model (see
  [`conformalInterval()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalInterval.md)).

- calib:

  The calibration data frame (response + predictors).

- newx:

  A one-row data frame of predictors for the new patient.

- alpha:

  Miscoverage level (default 0.1).

- seed:

  Optional integer random seed.

## Value

An `AnalysisResult` with `type = "conformal"` for the single patient.

## See also

[`conformalInterval()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalInterval.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(80)); d$y <- 1.5 * d$x + rnorm(80)
fit <- lm(y ~ x, data = d[1:40, ])
conformalPredict(fit, d[41:80, ], data.frame(x = 0.5))
#> <AnalysisResult> conformal 
#>   fields: point, lower, upper, alpha, method, quantile 
#>   provenance: 1 entr(ies)
```
