# Celeration (trend) line for a single-case phase

Fits a trend line to one phase, either by the White & Haring
split-middle method (the standard SCED hand method) or by ordinary least
squares.

## Usage

``` r
scedCelerationLine(
  value,
  time = NULL,
  method = c("split_middle", "ols"),
  ratio_period = 7
)
```

## Arguments

- value:

  Numeric phase observations, or pass `time`/`value` explicitly.

- time:

  Optional numeric session index (default `seq_along(value)`).

- method:

  `"split_middle"` (default) or `"ols"`.

- ratio_period:

  Sessions spanned by the reported celeration ratio (default 7, a weekly
  ratio for daily data).

## Value

An `AnalysisResult` (`type = "sced_celeration"`) whose `estimate` is the
slope (celeration per session), with `result$intercept`,
`result$fitted`, `result$bounce` (the ratio of the largest positive to
the largest negative residual magnitude) and `result$celeration_ratio`
(the multiplicative change over `ratio_period` sessions).

## References

White OR, Haring NG (1980). Exceptional Teaching, 2nd ed. Columbus, OH:
Merrill.

## See also

[`scedTwoSDBand()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTwoSDBand.md),
[`scedABAB()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedABAB.md)

## Examples

``` r
scedCelerationLine(c(4, 6, 5, 8, 7, 10, 9))
#> <AnalysisResult> sced_celeration 
#>   estimate: 1 
#>   method: celeration (split_middle) 
#>   fields: estimate, ci_lower, ci_upper, intercept, fitted, residuals, bounce, celeration_ratio, method, ratio_period 
#>   provenance: 1 entr(ies)
```
