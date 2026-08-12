# Cox proportional-hazards model with a PH check

Fits a Cox proportional-hazards model with
[`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html), tests
the proportional-hazards assumption with
[`survival::cox.zph`](https://rdrr.io/pkg/survival/man/cox.zph.html),
and warns when the global test is significant. Returns hazard ratios
with confidence intervals.

## Usage

``` r
coxModel(
  data,
  formula = NULL,
  time = "time",
  event = "event",
  covariates = NULL,
  conf_level = 0.95,
  ph_alpha = 0.05
)
```

## Arguments

- data:

  A data.frame.

- formula:

  A model formula whose response is `survival::Surv(time, event)`, or
  supply `time`/`event`/`covariates` instead.

- time, event, covariates:

  Alternative to `formula`: the time and event column names and a
  character vector of covariate columns.

- conf_level:

  Confidence level (default 0.95).

- ph_alpha:

  Significance level for the global PH warning (default 0.05).

## Value

An `AnalysisResult` (`type = "cox_model"`) whose `estimate` is the
hazard-ratio vector, with `result$coefficients` (coef, HR, SE, z, p,
CI), `result$ph_test` (the `cox.zph` table), `result$ph_violated` and
`result$fit`.

## References

Cox DR (1972). Regression models and life-tables. JRSS B, 34(2).
Grambsch PM, Therneau TM (1994). Proportional hazards tests and
diagnostics based on weighted residuals. Biometrika, 81(3).

## See also

[`survivalFit()`](https://x-biosignal.github.io/PhysioClinStats/reference/survivalFit.md),
[`milestoneHazard()`](https://x-biosignal.github.io/PhysioClinStats/reference/milestoneHazard.md)

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  coxModel(survival::lung, time = "time", event = "status",
           covariates = c("age", "sex", "ph.ecog"))
}
#> <AnalysisResult> cox_model 
#>   estimate: 1.0111282, 0.5754446, 1.5899912 
#>   method: Cox proportional hazards 
#>   fields: coefficients, ph_test, ph_violated, fit 
#>   provenance: 1 entr(ies)
```
