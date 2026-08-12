# Fit a mixed model for repeated measures (MMRM)

Fits the Mallinckrodt MMRM for a longitudinal randomised trial: fixed
effects for treatment, categorical time and their interaction (plus any
covariates), with a within-subject covariance over the repeated time
points. Uses the `mmrm` package when available - with Kenward-Roger or
Satterthwaite adjusted degrees of freedom - and falls back to
[`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html) otherwise. The
fallback reproduces the `mmrm` fixed effects for the unstructured
covariance (its primary use) and approximates the homogeneous
`cs`/`ar1`/`toeplitz` structures; it reports between-within
(containment) degrees of freedom and emits a message, since
Kenward-Roger / Satterthwaite df need `mmrm`.

## Usage

``` r
fitMMRM(
  data,
  response,
  treatment,
  time,
  subject,
  covariates = NULL,
  covariance = c("unstructured", "ar1", "compound-symmetry", "toeplitz"),
  df = c("kenward-roger", "satterthwaite")
)
```

## Arguments

- data:

  A long-format data.frame with one row per subject-time.

- response, treatment, time, subject:

  Column names of the outcome, the treatment factor, the categorical
  time factor and the subject id.

- covariates:

  Optional character vector of extra fixed-effect covariate columns
  (e.g. baseline value, stratifiers).

- covariance:

  Within-subject covariance: `"unstructured"` (default), `"ar1"`,
  `"compound-symmetry"` or `"toeplitz"`.

- df:

  Denominator-df method: `"kenward-roger"` (default) or
  `"satterthwaite"` (used by the `mmrm` path only).

## Value

An `AnalysisResult` (`type = "mmrm"`) whose `estimate` is the
fixed-effect coefficient vector, with `result$coefficients` (estimate,
SE, df, t, p), `result$backend` (`"mmrm"` or `"gls"`), `result$fit` and
`result$formula`.

## References

Mallinckrodt CH et al. (2008). Recommendations for the primary analysis
of continuous endpoints in longitudinal clinical trials. Drug
Information Journal, 42. Sabanes Bove D et al. (2023). mmrm: Mixed
Models for Repeated Measures. R package.

## See also

[`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md),
[`estimatedMarginalMeans()`](https://x-biosignal.github.io/PhysioClinStats/reference/estimatedMarginalMeans.md)

## Examples

``` r
if (requireNamespace("mmrm", quietly = TRUE)) {
  data(fev_data, package = "mmrm")
  fitMMRM(fev_data, "FEV1", "ARMCD", "AVISIT", "USUBJID",
          covariates = c("RACE", "SEX"))
}
#> <AnalysisResult> mmrm 
#>   estimate: 30.77747548,  1.53049977,  5.64356535,  0.32606192,  3.77423004,  4.83958845, 10.34211288, 15.05389826, -0.04192625, -0.69368537,  0.62422703 
#>   method: MMRM (unstructured, kenward-roger df, mmrm) 
#>   fields: coefficients, backend, fit, formula, covariance 
#>   provenance: 1 entr(ies)
```
