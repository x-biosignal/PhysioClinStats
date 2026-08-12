# Analyse an ABAB (reversal) single-case design

Orchestrates the SCED analysis of a reversal design: it separates the
phases (by default the four phases `A1, B1, A2, B2`), and for each
intervention contrast (`B1` vs `A1` and `B2` vs `A2`) computes the
non-overlap effect sizes, the 2-SD-band decision and per-phase
celeration lines.

## Usage

``` r
scedABAB(
  x,
  phases = c("A1", "B1", "A2", "B2"),
  contrasts = list(c("A1", "B1"), c("A2", "B2")),
  improvement = c("increase", "decrease"),
  value = NULL,
  phase = NULL,
  time = NULL,
  assay = NULL,
  channel = 1L
)
```

## Arguments

- x:

  A `PhysioExperiment` single-channel series, or a
  `data.frame(time, value, phase)`. Observations are ordered by `time`
  before analysis, so a `time` column need not be pre-sorted; without
  one, row order is used. Ordering matters because the Tau-U
  baseline-trend correction and the 2-SD consecutive-run rule are
  sequence-dependent.

- phases:

  Ordered character vector of the phase labels (default
  `c("A1","B1","A2","B2")`).

- contrasts:

  A list of length-2 `c(baseline, intervention)` phase-label pairs
  (default the two reversal contrasts).

- improvement:

  `"increase"` (default) or `"decrease"`.

- value, phase, time, assay, channel:

  Column/assay selectors passed to the extractor when `x` is a
  data.frame or PhysioExperiment.

## Value

An object of class `sced_abab` (a list) with `phase_data`, `contrasts`
(per contrast: NAP/Tau/Tau-U/PND/PEM results, the 2-SD decision and both
phase celeration lines) and a `summary` data.frame, with `print` and
`plot` methods.

## See also

[`scedNAP()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedNAP.md),
[`scedTwoSDBand()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedTwoSDBand.md),
[`scedCelerationLine()`](https://x-biosignal.github.io/PhysioClinStats/reference/scedCelerationLine.md)

## Examples

``` r
df <- data.frame(
  value = c(10, 11, 9, 18, 19, 20, 10, 12, 11, 21, 22, 20),
  phase = rep(c("A1", "B1", "A2", "B2"), each = 3))
res <- scedABAB(df)
res
#> Single-case ABAB analysis
#>   phases: A1, B1, A2, B2 (12 observations)
#>   improvement: increase
#> 
#>  contrast NAP Tau  Tau_U PND PEM band_flag p_value
#>  B1 vs A1   1   1 1.1111   1   1      TRUE  0.0441
#>  B2 vs A2   1   1 0.8889   1   1      TRUE  0.1072
```
