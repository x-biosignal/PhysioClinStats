# Pairwise contrasts of estimated marginal means

Convenience wrapper computing all pairwise differences of the estimated
marginal means for `specs`.

## Usage

``` r
pairwiseContrasts(model, specs, adjust = "tukey", level = 0.95, ...)
```

## Arguments

- model:

  An `AnalysisResult` from
  [`fitMixedModel()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMixedModel.md)
  /
  [`fitMMRM()`](https://x-biosignal.github.io/PhysioClinStats/reference/fitMMRM.md),
  or a raw model object `emmeans` understands.

- specs:

  A formula or character spec for the marginal means, e.g.
  `~ treatment | time`.

- adjust:

  Multiplicity adjustment passed to `emmeans` (default `"tukey"`).

- level:

  Confidence level for the intervals (default 0.95).

- ...:

  Further arguments forwarded to
  [`emmeans::emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.html).

## Value

An `AnalysisResult` (`type = "emmeans_contrasts"`) with
`result$contrasts`.

## See also

[`estimatedMarginalMeans()`](https://x-biosignal.github.io/PhysioClinStats/reference/estimatedMarginalMeans.md)
