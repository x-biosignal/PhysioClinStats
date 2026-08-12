# Two-standard-deviation band decision rule

Draws the Nelson (1984) 2-SD band from the baseline (A) phase - the
baseline mean plus and minus two baseline standard deviations - and
flags the intervention (B) phase for a systematic shift: the rule fires
when a run of `consecutive` or more successive B points falls on the
same side beyond the band.

## Usage

``` r
scedTwoSDBand(
  A_data,
  B_data = NULL,
  improvement = c("increase", "decrease"),
  k = 2,
  consecutive = 2L
)
```

## Arguments

- A_data, B_data:

  Numeric baseline and intervention observations.

- improvement:

  `"increase"` (default) or `"decrease"`; sets which side of the band
  counts as improvement.

- k:

  Band half-width in baseline SDs (default 2).

- consecutive:

  Number of successive out-of-band points that trigger the flag (default
  2, Nelson's rule).

## Value

An `AnalysisResult` (`type = "sced_2sd"`) whose `estimate` is the
logical decision, with `result$mean`, `result$sd`, `result$upper`,
`result$lower`, `result$outside` (per-B-point side: -1/0/+1) and
`result$first_run_at` (index of the first triggering point, or `NA`).

## References

Nelson LS (1984). The Shewhart control chart - tests for special causes.
Journal of Quality Technology, 16(4), 237-239. Applied to SCED in Gast &
Ledford (2014).

## See also

[`scedCelerationLine()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedCelerationLine.md),
[`scedABAB()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedABAB.md)

## Examples

``` r
scedTwoSDBand(c(10, 12, 11, 9, 10), c(15, 16, 17, 16, 18))
#> <AnalysisResult> sced_2sd 
#>   estimate: TRUE 
#>   method: 2-SD band 
#>   fields: estimate, ci_lower, ci_upper, mean, sd, upper, lower, outside, first_run_at, k, consecutive, improvement 
#>   provenance: 1 entr(ies)
```
