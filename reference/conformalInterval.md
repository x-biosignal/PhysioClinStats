# Conformal prediction interval for a fitted model

Distribution-free predictive intervals with finite-sample marginal
coverage at least `1 - alpha`, requiring only exchangeability. `"split"`
conformal (Lei et al. 2018) scores an independent calibration set with
the fitted model and takes the conformal quantile of the absolute
residuals; `"jackknife_plus"` (Barber et al. 2021) refits the model
leave-one-out over the calibration data and combines the leave-one-out
predictions and residuals.

## Usage

``` r
conformalInterval(
  model,
  calib,
  newx,
  alpha = 0.1,
  type = c("split", "jackknife_plus"),
  seed = NULL
)
```

## Arguments

- model:

  A fitted model supporting `predict(model, newdata=)` (and, for
  `"jackknife_plus"`, `update(model, data=)`), e.g. an `lm`.

- calib:

  A data frame with the response and predictors: an independent
  calibration set for `"split"`, or the full data to refit over for
  `"jackknife_plus"`.

- newx:

  A data frame of predictor values to predict for.

- alpha:

  Miscoverage level (default 0.1 for 90% intervals).

- type:

  `"split"` (default) or `"jackknife_plus"`.

- seed:

  Optional integer random seed recorded in the provenance.

## Value

An `AnalysisResult` (from PhysioCore) with `type = "conformal"` whose
`result` holds `point`, `lower`, `upper` (one per `newx` row), the
`alpha`/`method`/`quantile`, and a provenance log capturing the seed.

## Details

The coverage guarantee requires the calibration data to be exchangeable
with, and (for `"split"`) independent of, the data the model was trained
on. When the model uses a transformed response (e.g. `log(y) ~ x`) the
interval is on that transformed scale.

## References

Lei, J. et al. (2018). JASA 113(523):1094-1111. Barber, R.F. et al.
(2021). Ann Statist 49(1):486-507.

## See also

[`conformalPredict()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalPredict.md),
[`bayesEstimate()`](https://x-biosignal.github.io/PhysioClinStats/reference/bayesEstimate.md)

## Examples

``` r
set.seed(1)
train <- data.frame(x = rnorm(60)); train$y <- 2 * train$x + rnorm(60)
calib <- data.frame(x = rnorm(60)); calib$y <- 2 * calib$x + rnorm(60)
fit <- lm(y ~ x, data = train)
res <- conformalInterval(fit, calib, data.frame(x = c(-1, 0, 1)))
PhysioCore::resultValue(res)$lower
#> [1] -3.6257266 -1.6733436  0.2790395
```
