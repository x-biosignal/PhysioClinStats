# Percentage of non-overlapping data (PND)

The proportion of intervention-phase observations that exceed the most
extreme baseline observation (Scruggs, Mastropieri & Casto 1987).

## Usage

``` r
scedPND(A_data, B_data = NULL, improvement = c("increase", "decrease"))
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

An `AnalysisResult` (`type = "sced_pnd"`) whose `estimate` is the PND in
\[0, 1\].

## References

Scruggs TE, Mastropieri MA, Casto G (1987). The quantitative synthesis
of single-subject research. Remedial and Special Education, 8(2).

## See also

[`scedPEM()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedPEM.md),
[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md),
[`scedTauU()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTauU.md)

## Examples

``` r
scedPND(c(20, 20, 26, 25), c(28, 25, 30, 29))
#> <AnalysisResult> sced_pnd 
#>   estimate: 0.75 
#>   method: PND 
#>   fields: estimate, ci_lower, ci_upper, n_A, n_B, improvement 
#>   provenance: 1 entr(ies)
```
