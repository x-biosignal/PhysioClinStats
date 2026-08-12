# Non-overlap of all pairs (NAP)

The probability that a randomly chosen intervention observation exceeds
a randomly chosen baseline observation, with ties counted as one half
(Parker & Vannest 2009) - equivalent to the area under the ROC curve.
The confidence interval is the score interval of Newcombe (2006) and the
standard error the unbiased estimator, both as in SingleCaseES.

## Usage

``` r
scedNAP(
  A_data,
  B_data = NULL,
  improvement = c("increase", "decrease"),
  confidence = 0.95
)
```

## Arguments

- A_data, B_data:

  Numeric baseline (A) and intervention (B) observations, or pass a
  data.frame to the phase orchestrator
  [`scedABAB()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedABAB.md).

- improvement:

  `"increase"` (default) if higher scores are better, or `"decrease"` if
  lower scores are the therapeutic goal.

- confidence:

  Confidence level for the interval (default 0.95).

## Value

An `AnalysisResult` (`type = "sced_nap"`) whose `estimate` is the NAP,
with an analytic `uncertainty` interval and `result$se`,
`result$p_value`.

## References

Parker RI, Vannest KJ (2009). An improved effect size for single-case
research: non-overlap of all pairs. Behavior Therapy, 40(4), 357-367.
Newcombe RG (2006). Confidence intervals for an effect size measure
based on the Mann-Whitney statistic. Statistics in Medicine, 25.

## See also

[`scedTau()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTau.md),
[`scedTauU()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTauU.md)

## Examples

``` r
scedNAP(c(20, 20, 26, 25), c(28, 25, 30, 29))
#> <AnalysisResult> sced_nap 
#>   estimate: 0.90625 
#>   method: NAP 
#>   uncertainty: analytic 95% [0.4821987, 0.9886549]
#>   fields: estimate, ci_lower, ci_upper, se, z, p_value, n_A, n_B, improvement 
#>   provenance: 1 entr(ies)
```
