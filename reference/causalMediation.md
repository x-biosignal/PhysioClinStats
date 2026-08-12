# Model-based causal mediation analysis

Calls
[`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html)
after checking that the mediator and outcome models describe the same
ordered observations. Estimates rely on sequential ignorability,
consistency, positivity, absence of a treatment-induced mediator-outcome
confounder, and correct mediator and outcome models. An indirect
association, even when statistically significant, is not evidence that a
biological mechanism has been established.

## Usage

``` r
causalMediation(
  model_m,
  model_y,
  treat,
  mediator,
  control.value = 0,
  treat.value = 1,
  covariates = NULL,
  sims = 1000L,
  boot = FALSE,
  boot.ci.type = "perc",
  conf.level = 0.95,
  robustSE = FALSE,
  sensitivity = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- model_m, model_y:

  Supported fitted mediator and outcome models, or `AnalysisResult`
  objects whose `result$fit` contains those models.

- treat, mediator:

  Treatment and mediator variable names.

- control.value, treat.value:

  Treatment contrast supplied to the backend.

- covariates:

  Optional named list fixing covariates for a conditional mediation
  estimand.

- sims:

  Positive integer number of simulations or bootstrap replicates.

- boot:

  Use the nonparametric bootstrap.

- boot.ci.type:

  Bootstrap confidence-interval type.

- conf.level:

  Confidence level in `(0, 1)`.

- robustSE:

  Request heteroskedasticity-consistent uncertainty for supported
  `lm`/`glm` models.

- sensitivity:

  Run
  [`mediation::medsens()`](https://rdrr.io/pkg/mediation/man/medsens.html)
  for an all-linear model pair.

- seed:

  Optional integer seed. The caller's global RNG state is restored.

- ...:

  Further arguments passed to
  [`mediation::mediate()`](https://rdrr.io/pkg/mediation/man/mediate.html).

## Value

An `AnalysisResult` with a complete ACME/ADE/total/proportion table, the
backend object, settings, formulas, assumptions, and sensitivity
diagnostics.
