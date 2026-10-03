# Single-case (N-of-1) inference with PhysioClinStats

PhysioClinStats is the clinical inference engine of the ecosystem:
mixed-effects and MMRM longitudinal models, single-case (N-of-1)
designs, estimands with multiple imputation, and per-estimate
uncertainty. The heavy modelling backends (`lme4`, `mmrm`, `emmeans`,
`mice`, …) are optional, guarded `Suggests`.

This vignette uses only the **distribution-free single-case tools**,
which need no modelling backend, so it runs anywhere. Each estimator
reproduces its reference implementation (`SingleCaseES`, `scan`)
exactly.

``` r

library(PhysioClinStats)
```

## 1. A single case (AB design)

We have one participant measured repeatedly: a baseline phase (A) and an
intervention phase (B). Higher is better (e.g. comfortable gait speed,
scaled).

``` r

A <- c(20, 20, 26, 25, 22)          # baseline
B <- c(28, 25, 30, 29, 31, 30)      # intervention
```

## 2. Non-overlap effect sizes

The non-overlap family summarises how cleanly the intervention points
separate from baseline. `resultValue()` pulls the numbers out of the
returned `AnalysisResult`.

``` r

pnd <- scedPND(A, B)   # percentage of non-overlapping data
pem <- scedPEM(A, B)   # percentage exceeding the median
nap <- scedNAP(A, B)   # non-overlap of all pairs (= AUC), with a CI

PhysioExperiment::resultValue(pnd)$estimate
#> [1] 0.8333333
PhysioExperiment::resultValue(pem)$estimate
#> [1] 1
PhysioExperiment::resultValue(nap)[c("estimate", "ci_lower", "ci_upper", "p_value")]
#> $estimate
#> [1] 0.95
#> 
#> $ci_lower
#> [1] 0.6006868
#> 
#> $ci_upper
#> [1] 0.995108
#> 
#> $p_value
#> [1] 0.01371083
```

Tau rescales NAP to \[-1, 1\]; Tau-U additionally corrects for a
baseline trend (and exposes both the `SingleCaseES` and `scan`
denominator conventions).

``` r

PhysioExperiment::resultValue(scedTau(A, B))$estimate
#> [1] 0.9
PhysioExperiment::resultValue(scedTauU(A, B, method = "parker"))$estimate
#> [1] 0.8
```

## 3. Two-SD band

A Shewhart-style band flags runs of intervention points beyond 2 SD of
the baseline mean – a transparent visual-analysis aid.

``` r

band <- scedTwoSDBand(A, B)
res <- PhysioExperiment::resultValue(band)
res[c("upper", "lower", "first_run_at")]
#> $upper
#> [1] 28.1857
#> 
#> $lower
#> [1] 17.0143
#> 
#> $first_run_at
#> [1] 3
```

## 4. Frame the question as an estimand

Even for a single case it helps to state the estimand explicitly (ICH
E9(R1)): the treatment, population, endpoint, intercurrent-event
strategy and summary measure.

``` r

defineEstimand(
  treatment = "task-oriented training",
  population = "one post-stroke participant",
  endpoint = "comfortable gait speed",
  intercurrent = list(event = "missed sessions",
                      strategy = "treatment-policy"),
  summary_measure = "non-overlap of all pairs"
)
#> <estimand> ICH E9(R1)
#>   treatment:   task-oriented training 
#>   population:  one post-stroke participant 
#>   endpoint:    comfortable gait speed 
#>   intercurrent event: missed sessions 
#>   strategy:    treatment-policy 
#>   summary:     non-overlap of all pairs
```

## Where to go next

- Longitudinal group models –
  [`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md),
  [`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md),
  [`estimatedMarginalMeans()`](https://x-biosignal.github.io/PhysioClinStats/reference/estimatedMarginalMeans.md)
  – require their modelling backends and guard for them at run time.
- Estimands with missing data –
  [`multipleImputation()`](https://x-biosignal.github.io/PhysioClinStats/reference/multipleImputation.md),
  [`poolEstimates()`](https://x-biosignal.github.io/PhysioClinStats/reference/poolEstimates.md),
  [`analyseEstimand()`](https://x-biosignal.github.io/PhysioClinStats/reference/analyseEstimand.md).
- Causal tools –
  [`causalMediation()`](https://x-biosignal.github.io/PhysioClinStats/reference/causalMediation.md),
  [`targetTrialEmulate()`](https://x-biosignal.github.io/PhysioClinStats/reference/targetTrialEmulate.md).
- [`conformalPredict()`](https://x-biosignal.github.io/PhysioClinStats/reference/conformalPredict.md)
  gives distribution-free prediction intervals.

See
[`?PhysioClinStats`](https://x-biosignal.github.io/PhysioClinStats/reference/PhysioClinStats-package.md)
and each function’s help page for details.
