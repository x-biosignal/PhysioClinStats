# Tau-U with baseline-trend correction

Parker, Vannest, Davis & Sauber's (2011) Tau-U combines the A-vs-B
non-overlap with a correction for baseline (phase A) trend. Two
reference conventions differ on the denominator and this function
exposes both:

## Usage

``` r
scedTauU(
  A_data,
  B_data = NULL,
  improvement = c("increase", "decrease"),
  method = c("parker", "scan"),
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

- method:

  Denominator convention, `"parker"` (default) or `"scan"`.

- confidence:

  Confidence level for the interval (default 0.95).

## Value

An `AnalysisResult` (`type = "sced_tau_u"`) whose `estimate` is Tau-U,
with `result$S`, `result$se_S`, `result$z`, `result$p_value` and the
component S values.

## Details

- `method = "parker"` (default):

  matches
  [`SingleCaseES::Tau_U`](https://jepusto.github.io/SingleCaseES/reference/Tau_U.html):
  \\(S\_{AB} - S\_{trend A}) / (m n)\\, where \\S\_{AB}\\ is the A-vs-B
  Kendall S and \\S\_{trend A}\\ the within-baseline Kendall S. Can
  exceed 1 when a baseline trend runs counter to improvement.

- `method = "scan"`:

  matches
  [`scan::tau_u()`](https://jazznbass.github.io/scan/reference/tau_u.html)'s
  "A vs. B - Trend A" row: the same S numerator over a tie-corrected
  pair count \\D\\, giving a Kendall tau-b in \[-1, 1\].

The Z statistic and p-value come from the tie-corrected variance of the
Kendall S numerator and are identical under both methods (they depend on
S, not the denominator), matching `scan`'s Z and p.

## References

Parker RI, Vannest KJ, Davis JL, Sauber SB (2011). Combining nonoverlap
and trend for single-case research: Tau-U. Behavior Therapy, 42(2),
284-299.

## See also

[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md),
[`scedTau()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTau.md)

## Examples

``` r
scedTauU(c(20, 20, 26, 25, 22, 23), c(28, 25, 24, 27, 30, 30, 29, 28))
#> <AnalysisResult> sced_tau_u 
#>   estimate: 0.7708333 
#>   method: Tau-U (parker) 
#>   fields: estimate, ci_lower, ci_upper, S, S_AB, S_trendA, D, z, p_value, method, n_A, n_B, improvement 
#>   provenance: 1 entr(ies)
```
