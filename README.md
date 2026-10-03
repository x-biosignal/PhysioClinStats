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
# the containers build on Bioconductor, so its repositories are needed too
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages("PhysioClinStats",
  repos = c("https://x-biosignal.r-universe.dev", BiocManager::repositories()))
```

## Causal mediation

```r
library(PhysioClinStats)

if (requireNamespace("mediation", quietly = TRUE) &&
    requireNamespace("PhysioCore", quietly = TRUE)) {
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

# A small synthetic longitudinal cohort (no external data needed): one eventual
# outcome per participant, observed over a shared weekly follow-up grid.
set.seed(1)
weeks <- 0:12
make_participant <- function(i) {
  baseline <- rnorm(1, 50, 10)
  treated <- rbinom(1, 1, 0.5)
  recovered <- rbinom(1, 1, plogis((baseline - 50) / 12 + 1.2 * treated))
  data.frame(
    participant = sprintf("p%03d", i),
    week = weeks,
    eligible = TRUE,
    age = round(rnorm(1, 60, 8), 1),
    baseline_score = round(baseline, 1),
    current_score = round(baseline + cumsum(rnorm(length(weeks), 1.5 * treated, 2)), 1),
    treated = treated,
    recovered = recovered
  )
}
longitudinal_data <- do.call(rbind, lapply(seq_len(150), make_participant))

result <- targetTrialEmulate(
  longitudinal_data, protocol,
  id = "participant", time = "week", treatment = "treated",
  outcome = "recovered",
  baseline_covariates = c("age", "baseline_score"),
  time_varying_covariates = "current_score",
  estimand = "per_protocol"
)
if (requireNamespace("PhysioCore", quietly = TRUE)) {
  PhysioCore::resultValue(result)$effects
}
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
