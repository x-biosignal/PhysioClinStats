# Kaplan-Meier survival fit with attainment-time summary

Fits a Kaplan-Meier estimator with
[`survival::survfit`](https://rdrr.io/pkg/survival/man/survfit.html),
optionally stratified by a grouping variable, and returns the median and
quartile attainment times with their confidence intervals.

## Usage

``` r
survivalFit(
  data,
  time = "time",
  event = "event",
  group = NULL,
  conf_level = 0.95
)
```

## Arguments

- data:

  A data.frame with the survival columns.

- time, event:

  Column names of the follow-up time and the event indicator (1 = event,
  0 = censored).

- group:

  Optional grouping-variable column name for stratified curves.

- conf_level:

  Confidence level (default 0.95).

## Value

An object of class `physio_km` (a list) with `fit` (the `survfit`),
`quantiles` (median / Q1 / Q3 attainment times with CIs, per group) and
`summary` table, plus `print` and `plot` methods.

## References

Kaplan EL, Meier P (1958). Nonparametric estimation from incomplete
observations. JASA, 53(282), 457-481.

## See also

[`coxModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/coxModel.md),
[`timeToMilestone()`](https://x-biosignal.github.io/PhysioClinStats/reference/timeToMilestone.md)

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  survivalFit(survival::lung, "time", "status", group = "sex")
}
#> Kaplan-Meier fit
#> Call: survfit(formula = form, data = data, conf.int = conf_level)
#> 
#>         n events median 0.95LCL 0.95UCL
#> sex=1 138    112    270     212     310
#> sex=2  90     53    426     348     550
#> 
#> Attainment times (median / quartiles):
#>  group  Q1 median  Q3 median_lower median_upper
#>  sex=1 144    270 457          212          310
#>  sex=2 226    426 687          348          550
```
