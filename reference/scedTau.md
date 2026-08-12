# Tau non-overlap (Tau)

The rank-correlation non-overlap effect size, \\\tau = 2\\\mathrm{NAP} -
1\\ (Parker et al. 2011), with the SE and interval derived from
[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md).

## Usage

``` r
scedTau(
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

An `AnalysisResult` (`type = "sced_tau"`) whose `estimate` is Tau in
\[-1, 1\].

## See also

[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md),
[`scedTauU()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTauU.md)

## Examples

``` r
scedTau(c(20, 20, 26, 25), c(28, 25, 30, 29))
#> <AnalysisResult> sced_tau 
#>   estimate: 0.8125 
#>   method: Tau 
#>   uncertainty: analytic 95% [-0.03560253, 0.9773097]
#>   fields: estimate, ci_lower, ci_upper, se, z, p_value, n_A, n_B, improvement 
#>   provenance: 1 entr(ies)
```
