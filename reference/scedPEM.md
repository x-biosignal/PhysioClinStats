# Percentage exceeding the median (PEM)

The proportion of intervention-phase observations beyond the baseline
median, counting ties as one half (Ma 2006).

## Usage

``` r
scedPEM(A_data, B_data = NULL, improvement = c("increase", "decrease"))
```

## Arguments

- A_data, B_data:

  Numeric baseline (A) and intervention (B) observations, or pass a
  data.frame to the phase orchestrator
  [`scedABAB()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedABAB.md).

- improvement:

  `"increase"` (default) if higher scores are better, or `"decrease"` if
  lower scores are the therapeutic goal.

## Value

An `AnalysisResult` (`type = "sced_pem"`) whose `estimate` is the PEM.

## References

Ma H-H (2006). An alternative method for quantitative synthesis of
single-subject research: percentage of data points exceeding the median.
Behavior Modification, 30(5), 598-617.

## See also

[`scedPND()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedPND.md),
[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md)

## Examples

``` r
scedPEM(c(20, 20, 26, 25), c(28, 25, 30, 29))
#> <AnalysisResult> sced_pem 
#>   estimate: 1 
#>   method: PEM 
#>   fields: estimate, ci_lower, ci_upper, baseline_median, n_A, n_B, improvement 
#>   provenance: 1 entr(ies)
```
