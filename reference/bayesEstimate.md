# Bayesian estimate of a mean (conjugate Normal-Normal or Stan)

For a numeric sample, returns the analytic Normal-Normal conjugate
posterior for the mean: a normal prior `N(prior_mean, prior_sd^2)`
combined with a normal likelihood of known scale `sigma` gives a normal
posterior in closed form (a flat prior, `prior_sd = Inf`, reproduces the
sampling posterior `N(mean(y), sigma^2/n)`). For a model `formula` the
estimation is delegated to rstanarm/brms when installed, and errors
informatively otherwise – the analytic numeric-vector path always runs.

## Usage

``` r
bayesEstimate(
  y,
  data = NULL,
  prior_mean = 0,
  prior_sd = Inf,
  sigma = NULL,
  level = 0.95,
  seed = NULL
)
```

## Arguments

- y:

  A numeric vector (analytic path) or a model `formula` (Stan path).

- data:

  Data frame for the `formula` path.

- prior_mean, prior_sd:

  Normal prior mean and SD (default 0 and `Inf`, a flat prior).

- sigma:

  Known likelihood SD; `NULL` uses the sample SD of `y`.

- level:

  Credible level (default 0.95).

- seed:

  Optional integer random seed recorded in the provenance.

## Value

An `AnalysisResult` (from PhysioCore) with `type = "bayes"` whose
`result` holds `posterior_mean`, `posterior_sd`, `ci_lower`, `ci_upper`,
the `level`/`method`, plus a provenance log with the seed.

## References

Gelman, A. et al. (2013). Bayesian Data Analysis, 3rd ed.

## See also

[`credibleInterval()`](https://x-biosignal.github.io/PhysioClinStats/reference/credibleInterval.md),
[`conformalInterval()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalInterval.md)

## Examples

``` r
set.seed(1)
bayesEstimate(rnorm(30, mean = 5), prior_mean = 0, prior_sd = 10)
#> <AnalysisResult> bayes 
#>   fields: posterior_mean, posterior_sd, ci_lower, ci_upper, level, method 
#>   provenance: 1 entr(ies)
```
