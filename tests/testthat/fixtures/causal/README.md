# WS8-11 causal-inference fixtures

These fixtures are offline numerical references. Regenerate the RDS files with:

```r
source("tests/testthat/fixtures/causal/generate-fixtures.R")
```

- `mediation-jobs.rds` uses the `jobs` framing data shipped by `mediation`,
  linear model formulas recorded in the object, seed 481, and 500
  quasi-Bayesian simulations. The direct `mediate()` and `medsens()` outputs
  are stored with R and package versions.
- `baseline-iptw.rds` uses a locally generated, seeded binary-treatment
  dataset and stores direct base-R logistic, WeightIt, and ipw ATE weights.
- `target-trial-*.csv` is an author-reconstructed, hand-auditable static
  sustained-treatment example. Every expected clone row states its adherence,
  artificial-censor decision, denominator/numerator interval probability, and
  cumulative raw weight. No package cloning or weighting function generated
  these CSV tables.

The target-trial example uses denominator probability 0.5 and numerator
probability 0.4 on every adherent interval, hence interval weight 0.8. A clone
is censored immediately before its first treatment deviation and carries its
pre-deviation cumulative weight on the audit row.
