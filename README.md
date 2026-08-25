# PhysioClinStats

<!-- badges: start -->
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioClinStats)](https://x-biosignal.r-universe.dev/PhysioClinStats)
<!-- badges: end -->

Clinical inference engine for the [x-biosignal](https://github.com/x-biosignal)
ecosystem: mixed-effects/MMRM longitudinal models, single-case (N-of-1) designs,
estimands with multiple imputation, causal mediation, declared target-trial
emulation, and per-estimate uncertainty. Heavy modelling backends are optional
(guarded) Suggests.

## Installation

```r
install.packages("PhysioClinStats",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

## Causal mediation

```r
if (requireNamespace("mediation", quietly = TRUE)) {
  data(jobs, package = "mediation")
  model_m <- lm(job_seek ~ treat + econ_hard + sex + age, data = jobs)
  model_y <- lm(
    depress2 ~ treat + job_seek + econ_hard + sex + age,
    data = jobs
  )
  result <- causalMediation(
    model_m, model_y,
    treat = "treat", mediator = "job_seek",
    sims = 1000, seed = 481
  )
  PhysioCore::resultValue(result)$effects
}
```

The estimates require sequential ignorability, consistency, positivity, no
treatment-induced mediator-outcome confounder, and correctly specified models.
An indirect effect is not proof of a biological mechanism.

## Target-trial emulation

```r
protocol <- targetTrialProtocol(
  eligibility = function(data) data$eligible,
  treatment_strategies = list(never = 0, always = 1),
  assignment = "clone eligible participants",
  time_zero = 0,
  follow_up = 12,
  outcome = "recovery by week 12",
  causal_contrast = "always versus never",
  analysis_plan = "clone-censor-weight"
)

result <- targetTrialEmulate(
  longitudinal_data, protocol,
  id = "participant", time = "week", treatment = "treated",
  outcome = "recovered",
  baseline_covariates = c("age", "baseline_score"),
  time_varying_covariates = "current_score",
  estimand = "per_protocol"
)
PhysioCore::resultValue(result)$effects
```

Time zero and every protocol component are declared rather than inferred.
Inspect positivity, balance, censoring, and raw-versus-analysis weight
diagnostics before interpreting a contrast. Cox hazard ratios are
non-collapsible and are not marginal risk ratios.

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev). Community and
policy documents live in the umbrella repository:

- [Code of Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
