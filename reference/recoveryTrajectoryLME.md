# Fit a population recovery-trajectory NLME

Fits a nonlinear mixed-effects recovery curve to longitudinal panel data
with subject-level random effects (partial pooling), and returns the
population fixed effects together with per-subject predicted asymptote,
rate, and time-to-90%-recovery.

## Usage

``` r
recoveryTrajectoryLME(
  data,
  subject,
  time,
  outcome,
  model = c("exponential", "asymptotic", "logistic")
)
```

## Arguments

- data:

  A long-format data frame.

- subject, time, outcome:

  Column names (character) for the grouping factor, the time variable,
  and the response.

- model:

  `"exponential"` (default), `"asymptotic"` (an alias), or `"logistic"`.

## Value

An `AnalysisResult` (type `"recovery_trajectory"`) whose `result` holds
the `fixed_effects`, `random_effects`, a per-subject data frame
(`asymptote`, `rate`, `time_to_90`), and the fitted `nlme` object.

## Details

Models: `"exponential"`/`"asymptotic"` use the self-starting asymptotic
curve \\y = A + (R_0 - A)\\e^{-e^{lrc} t}\\ (`SSasymp`); `"logistic"`
uses `SSlogis`.

## References

Pinheiro & Bates 2000 (nlme); Lindstrom & Bates 1990.

## See also

[`proportionalRecoveryRule()`](https://x-biosignal.github.io/PhysioClinStats/reference/proportionalRecoveryRule.md),
[`latentClassGrowth()`](https://x-biosignal.github.io/PhysioClinStats/reference/latentClassGrowth.md)

## Examples

``` r
set.seed(1)
df <- do.call(rbind, lapply(1:12, function(s) {
  A <- 60 + rnorm(1, 0, 5); rate <- 0.3 * exp(rnorm(1, 0, 0.2))
  t <- 0:8; data.frame(subject = s, time = t,
    y = A * (1 - exp(-rate * t)) + rnorm(9, 0, 2))
}))
fit <- recoveryTrajectoryLME(df, "subject", "time", "y")
#> Warning: Iteration 1, LME step: nlminb() did not converge (code = 1). Do increase 'msMaxIter'!
```
