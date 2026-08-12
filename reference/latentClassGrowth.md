# Latent-class growth analysis of longitudinal trajectories

Fits a finite mixture of polynomial growth trajectories to panel data
(each subject belongs to one latent class), selecting the number of
classes by BIC. Returns per-subject class assignments, posterior class
probabilities, and the normalised entropy of the classification.

## Usage

``` r
latentClassGrowth(
  data,
  subject,
  time,
  outcome,
  n_classes = 1:3,
  degree = 1,
  seed = NULL
)
```

## Arguments

- data:

  A long-format data frame.

- subject, time, outcome:

  Column names for the grouping factor, the time variable, and the
  response.

- n_classes:

  Integer vector of class counts to compare (default `1:3`); the
  BIC-minimising count is selected.

- degree:

  Polynomial degree of the per-class growth curve (default 1, linear).

- seed:

  Optional RNG seed (flexmix uses random starts).

## Value

An `AnalysisResult` (type `"latent_class_growth"`) whose `result` holds
the selected `n_classes`, the `bic` table, the per-subject `assignment`,
the `posterior` probabilities, the `entropy`, and the fitted flexmix
`model`.

## References

Nagin 2005 (group-based trajectory modelling); Gruen & Leisch 2008
(flexmix).

## See also

[`recoveryTrajectoryLME()`](https://x-biosignal.github.io/PhysioClinStats/reference/recoveryTrajectoryLME.md),
[`proportionalRecoveryRule()`](https://x-biosignal.github.io/PhysioClinStats/reference/proportionalRecoveryRule.md)

## Examples

``` r
set.seed(1)
df <- do.call(rbind, lapply(1:40, function(s) {
  fast <- s <= 20; t <- 0:6
  y <- (if (fast) 5 * t else 0.5 * t) + rnorm(7, 0, 1.5)
  data.frame(subject = s, time = t, y = y)
}))
latentClassGrowth(df, "subject", "time", "y", n_classes = 1:3)
#> <AnalysisResult> latent_class_growth 
#>   estimate: 2 
#>   method: flexmix_lcga 
#>   fields: n_classes, bic, assignment, posterior, entropy, class_sizes, model 
```
