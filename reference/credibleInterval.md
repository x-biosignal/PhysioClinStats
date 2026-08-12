# Credible interval from posterior samples

The highest-density interval (HDI, the narrowest interval carrying
`level` probability) or the equal-tailed quantile interval, from a
vector of posterior draws.

## Usage

``` r
credibleInterval(x, level = 0.95, method = c("hdi", "quantile"))
```

## Arguments

- x:

  Numeric vector of posterior samples.

- level:

  Credible level (default 0.95).

- method:

  `"hdi"` (default, narrowest) or `"quantile"` (equal-tailed).

## Value

Numeric `c(lower, upper)`.

## See also

[`bayesEstimate()`](https://x-biosignal.github.io/PhysioClinStats/reference/bayesEstimate.md)

## Examples

``` r
set.seed(1)
credibleInterval(rnorm(2000), method = "hdi")
#> [1] -2.001003  2.075245
```
