# Evaluate the proportional recovery rule with an artefact guard

Tests whether recovery follows the proportional recovery rule \\\Delta
\approx p \\(\mathrm{max} - \mathrm{initial})\\ (default \\p = 0.7\\),
classifies fitters vs non-fitters, handles ceiling effects, and runs two
diagnostics against the known statistical artefact: a
mixture-vs-single-line comparison and a regression-to-the-mean
(label-shuffle) control.

## Usage

``` r
proportionalRecoveryRule(
  initial,
  followup,
  max_score,
  expected_prop = 0.7,
  fitter_tol = 0.2,
  ceiling_frac = 0.9,
  n_shuffle = 999,
  seed = NULL
)
```

## Arguments

- initial:

  Numeric baseline scores.

- followup:

  Numeric follow-up scores (same length).

- max_score:

  Scalar maximum attainable score (the recovery ceiling).

- expected_prop:

  Expected recovery proportion (default 0.7).

- fitter_tol:

  Half-width around `expected_prop` for a "fitter" (default 0.2).

- ceiling_frac:

  Subjects with `initial >= ceiling_frac * max_score` are treated as
  ceiling cases (default 0.9) and excluded from the slope fit.

- n_shuffle:

  Label-shuffles for the regression-to-the-mean control (default 999).

- seed:

  Optional RNG seed.

## Value

An `AnalysisResult` (type `"proportional_recovery"`) whose `result`
holds the fitted `slope` (+ CI and test vs `expected_prop`), per-subject
`classification`, the `mixture_preferred` flag, and
`artefact_suspected`.

## Details

A genuine proportional rule implies the follow-up regresses on the
initial score with slope \\1 - p\\; a pure regression-to-the-mean /
coupling artefact instead reproduces the \\\Delta\\-vs-potential slope
when the follow-up scores are shuffled across subjects. When the
observed slope is not distinguishable from that shuffled null,
`artefact_suspected` is set.

## References

Prabhakaran 2008; Winters 2015; Hawe, Scott & Dukelow 2019
(proportional-recovery artefact).

## See also

[`recoveryTrajectoryLME()`](https://x-biosignal.github.io/PhysioClinStats/reference/recoveryTrajectoryLME.md),
[`latentClassGrowth()`](https://x-biosignal.github.io/PhysioClinStats/reference/latentClassGrowth.md)

## Examples

``` r
set.seed(1)
init <- runif(60, 5, 55)
fu <- init + 0.7 * (66 - init) + rnorm(60, 0, 3)   # genuine 0.7 rule
proportionalRecoveryRule(init, fu, max_score = 66)
#> <AnalysisResult> proportional_recovery 
#>   estimate: 0.7098355 
#>   method: proportional_recovery_rule 
#>   fields: slope, slope_ci, expected_prop, p_vs_expected, followup_on_initial_slope, n, n_ceiling, n_fitter, classification, observed_prop, mixture_preferred, mixture_bic, shuffle_slope_ci, artefact_suspected 
```
