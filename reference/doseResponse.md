# Therapy dose-response model

Fits an outcome (or change from baseline) as a function of a continuous
/ ordinal therapy dose, with optional confounder adjustment, subject
random intercepts and a non-linear dose form. Returns the fitted dose
effect, a predicted dose-response curve and the coefficient table.

## Usage

``` r
doseResponse(
  data,
  outcome,
  dose,
  covariates = NULL,
  subject = NULL,
  baseline = NULL,
  form = c("linear", "log", "spline", "emax"),
  df = 3,
  n_grid = 50
)
```

## Arguments

- data:

  A `data.frame` with the columns named below.

- outcome:

  Name of the outcome column.

- dose:

  Name of the therapy-dose column (hours / repetitions / sessions).

- covariates:

  Optional character vector of confounder columns to adjust for (added
  linearly).

- subject:

  Optional subject-id column; when given, a random intercept per subject
  is fit via
  [`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md)
  (needs `lme4`).

- baseline:

  Optional baseline column; when given it is added as an adjustment
  covariate (ANCOVA-style change analysis).

- form:

  Dose form: `"linear"`, `"log"` (needs positive dose), `"spline"`
  (natural cubic spline, `df` terms) or `"emax"`
  (`E0 + Emax*dose/(ED50+dose)`, fit by
  [`stats::nls()`](https://rdrr.io/r/stats/nls.html); no
  subject/covariates).

- df:

  Spline degrees of freedom for `form = "spline"` (default 3).

- n_grid:

  Number of dose points in the predicted curve (default 50).

## Value

A
[PhysioCore::AnalysisResult](https://x-biosignal.github.io/PhysioCore//reference/AnalysisResult.html)
of `type = "dose_response"`: `estimate` is the dose effect (slope,
spline coefficients, or `E0`/`Emax`/`ED50`); `result` holds
`coefficients`, the `curve` (`dose`, `predicted`), the `fit` and the
`form`.

## References

Ruberg SJ (1995) Dose response studies. *J Biopharm Stat* 5:1-14.

## See also

[`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md),
[`proportionalRecoveryRule()`](https://x-biosignal.github.io/PhysioClinStats/reference/proportionalRecoveryRule.md)

## Examples

``` r
set.seed(1)
d <- data.frame(dose = rep(c(0, 5, 10, 20), each = 15))
d$change <- 0.4 * d$dose + rnorm(nrow(d), 0, 2)
dr <- doseResponse(d, outcome = "change", dose = "dose", form = "linear")
PhysioCore::resultValue(dr)$coefficients
#>          term  estimate         se  statistic            p
#> 1 (Intercept) 0.1382792 0.34475004  0.4010999 6.898205e-01
#> 2        dose 0.4087947 0.03009225 13.5847166 1.137395e-19
```
